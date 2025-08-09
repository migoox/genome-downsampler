#include "qmcp-solver/cpu_crawling_solver.hpp"
#if defined(__clang__) && defined(XRAY_INSTRUMENTATION_ENABLED)
#include <xray/xray_interface.h>
#endif

#if defined(__clang__) && defined(XRAY_INSTRUMENTATION_ENABLED)
#define XRAY_INSTRUMENT [[clang::xray_always_instrument]] __attribute__((noinline))
#else
#define XRAY_INSTRUMENT
#endif

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <future>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <queue>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "bam-api/bam_api.hpp"
#include "bam-api/read.hpp"
#include "qmcp-solver/solver.hpp"

template <typename T>
void printVector(std::vector<T>& v)
{
    for (auto e : v) std::cout << e << " ";
    std::cout << "\n";
}

// TODO find better name
struct HelperRead
{
    bam_api::Index endInd;
    bam_api::ReadIndex readInd;
    HelperRead(uint32_t readInd, bam_api::Index endInd) : endInd(endInd), readInd(readInd) {}
    bool operator<(const HelperRead& other) const { return endInd < other.endInd; }

    bool operator>(const HelperRead& other) const { return endInd > other.endInd; }

    bool operator<=(const HelperRead& other) const { return endInd <= other.endInd; }

    bool operator>=(const HelperRead& other) const { return endInd >= other.endInd; }

    bool operator==(const HelperRead& other) const { return endInd == other.endInd; }

    bool operator!=(const HelperRead& other) const { return endInd != other.endInd; }

    bool operator<(const bam_api::Index other) const { return endInd < other; }

    bool operator>(const bam_api::Index other) const { return endInd > other; }

    bool operator<=(const bam_api::Index other) const { return endInd <= other; }

    bool operator>=(const bam_api::Index other) const { return endInd >= other; }

    bool operator==(const bam_api::Index other) const { return endInd == other; }

    bool operator!=(const bam_api::Index other) const { return endInd != other; }

    friend std::ostream& operator<<(std::ostream& os, const HelperRead& hr)
    {
        os << "(ind: " << hr.readInd << " end: " << hr.endInd << ")";
        return os;
    }
};

template <typename T>
void shiftRightFromIndex(std::vector<T>& vec, size_t index)
{
    if (vec.empty() || index >= vec.size() - 1) return;

    for (size_t i = vec.size() - 1; i > index; --i)
    {
        vec[i] = vec[i - 1];
    }
}

template <class T>
// bigger elements on smaller index's
struct SizedPriorityQueue
{
    std::vector<T> buffer;
    uint32_t M;
    explicit SizedPriorityQueue(uint32_t M) : M(M) { buffer.reserve(M); }
    bool tryAdd(T val)
    {
        int m = buffer.size();

        if (m < M)
        {
            buffer.push_back(val);  // TMP just to resize vec, val will be removed from vec
            m++;
        }
        else if (buffer[m - 1] >= val)
            return false;
        for (int x = m - 2; x >= 0; x--)
        {
            if (buffer[x] >= val)
            {
                shiftRightFromIndex(buffer, x + 1);
                buffer[x + 1] = val;
                return true;
            }
        }
        shiftRightFromIndex(buffer, 0);
        buffer[0] = val;
        return true;
    }

    void clear() { buffer.clear(); }
    uint32_t size() const { return buffer.size(); }

    void eraseLess(int num)
    {
        if (buffer[buffer.size() - 1] >= num) return;
        for (int x = buffer.size() - 1; x >= 0; x--)
            if (buffer[x] >= num)
            {
                buffer.erase(buffer.begin() + x + 1, buffer.end());
                return;
            }
    }

    void eraseFront(int num)
    {
        if (num >= buffer.size())
        {
            buffer.clear();
            return;
        }
        buffer.erase(buffer.begin(), buffer.begin() + num);
    }

    T& operator[](size_t index) { return buffer[index]; }
};

// test time 58s
void findBestSamples_tryFindShorter(uint32_t max_coverage,
                                    const bam_api::SOAPairedReads& input_sequence,
                                    uint32_t start_ind, uint32_t end_ind,
                                    std::vector<uint32_t>& start_index_count_pref_sum,
                                    std::vector<bam_api::ReadIndex>& read_ids_by_ascending_start_id,
                                    SizedPriorityQueue<HelperRead>& best_M,
                                    SizedPriorityQueue<HelperRead>& best_M_prev)
{
    for (uint32_t x = start_ind; x <= end_ind; x++)
    {
        int subseq_last_sorted_ind = start_index_count_pref_sum[x + 1] - 1;
        int offset = static_cast<int>(std::min(
            max_coverage, start_index_count_pref_sum[x + 1] - start_index_count_pref_sum[x]));

        for (int sorted_ind = subseq_last_sorted_ind;
             sorted_ind >= 0 && sorted_ind > subseq_last_sorted_ind - offset; sorted_ind--)
        {
            bam_api::ReadIndex ind = read_ids_by_ascending_start_id[sorted_ind];
            best_M.tryAdd(HelperRead(ind, input_sequence.end_inds[ind]));
        }
    }

    for (int x = 0; x < best_M_prev.size(); x++) best_M.tryAdd(best_M_prev[x]);
}

// test time 38s
void findBestSamples(uint32_t max_coverage, const bam_api::SOAPairedReads& input_sequence,
                     uint32_t start_ind, uint32_t end_ind,
                     std::vector<uint32_t>& start_index_count_pref_sum,
                     std::vector<bam_api::ReadIndex>& read_ids_by_ascending_start_id,
                     SizedPriorityQueue<HelperRead>& best_M,
                     SizedPriorityQueue<HelperRead>& best_M_prev)
{
    bool found = false;
    uint32_t ind_offset = 0;
    do
    {
        found = false;
        for (uint32_t x = start_ind; x <= end_ind; x++)
        {
            int64_t subseq_sorted_ind =
                static_cast<int64_t>(start_index_count_pref_sum[x + 1]) - 1 - ind_offset;
            if (subseq_sorted_ind < static_cast<int64_t>(start_index_count_pref_sum[x])) continue;
            bam_api::ReadIndex ind = read_ids_by_ascending_start_id[subseq_sorted_ind];
            found |= best_M.tryAdd(HelperRead(ind, input_sequence.end_inds[ind]));
        }
        ind_offset++;
    } while (found);

    for (int x = 0; x < best_M_prev.size(); x++) best_M.tryAdd(best_M_prev[x]);
}

XRAY_INSTRUMENT
std::unique_ptr<qmcp::Solution> qmcp::CpuCrawlingSolver::solve(uint32_t max_coverage,
                                                               bam_api::BamApi& bam_api)
{
    const bam_api::SOAPairedReads& input_reads = bam_api.get_paired_reads_soa();
    const uint32_t reads_count = input_reads.get_reads_count();

    uint32_t start_ind = 0;
    uint32_t end_ind = 0;
    std::vector<uint32_t> index_samples_more_needed;
    std::vector<uint32_t> start_index_count_pref_sum;
    std::vector<bam_api::ReadIndex> read_ids_by_ascending_start_id;
    SizedPriorityQueue<HelperRead> best_M1(max_coverage);  // M stands for max_coverage
    SizedPriorityQueue<HelperRead> best_M2(max_coverage);
    SizedPriorityQueue<HelperRead>& best_M = best_M1;
    SizedPriorityQueue<HelperRead>& best_M_prev = best_M2;
    std::vector<int64_t> ret_vec(reads_count, 0);

    index_samples_more_needed.resize(input_reads.ref_genome_length, max_coverage);

    // this can be removed by implementing sort by count
    start_index_count_pref_sum.resize(input_reads.ref_genome_length + 1, 0);
    read_ids_by_ascending_start_id.resize(reads_count);

    for (int x = 0; x < reads_count; x++)
    {
        start_index_count_pref_sum[input_reads.start_inds[x] + 1]++;
        read_ids_by_ascending_start_id[x] = x;
    }

    for (int x = 1; x < start_index_count_pref_sum.size(); x++)
        start_index_count_pref_sum[x] += start_index_count_pref_sum[x - 1];
    start_index_count_pref_sum.push_back(
        start_index_count_pref_sum[start_index_count_pref_sum.size() - 1]);

    std::sort(read_ids_by_ascending_start_id.begin(), read_ids_by_ascending_start_id.end(),
              [&](int ind1, int ind2)
              {
                  if (input_reads.start_inds[ind1] == input_reads.start_inds[ind2])
                      return input_reads.end_inds[ind1] < input_reads.end_inds[ind2];
                  return input_reads.start_inds[ind1] < input_reads.start_inds[ind2];
              });

    while (start_ind < input_reads.ref_genome_length)
    {
        // find best samples
        best_M.clear();
        // findBestSamples_tryFindShorter(max_coverage, input_sequence, start_ind, end_ind,
        //                                start_inex_count_pref_sum, read_ids_by_ascending_start_id
        //                                , best_M, best_M_prev);
        findBestSamples(max_coverage, input_reads, start_ind, end_ind, start_index_count_pref_sum,
                        read_ids_by_ascending_start_id, best_M, best_M_prev);

        if (best_M.size() == 0)
        {
            start_ind = end_ind + 1;
            end_ind++;
            continue;
        }
        uint32_t N =
            std::min(index_samples_more_needed[end_ind], best_M.size());  // N samples to take
        // update coverage array and ret vec
        ret_vec[best_M[0].readInd] = 1;
        int covered_by = 1;
        for (int x = 1; x < N; x++)
        {
            ret_vec[best_M[x].readInd] = 1;
            for (int y = best_M[x - 1].endInd; y > best_M[x].endInd; y--)
                index_samples_more_needed[y] -= covered_by;
            covered_by++;
        }
        // update start, end
        start_ind = end_ind + 1;
        int it;
        for (it = best_M[N - 1].endInd; N < index_samples_more_needed[it] &&
                                        it > input_reads.start_inds[best_M[N - 1].readInd];
             it--)
            index_samples_more_needed[it] -= N;
        end_ind = it + 1;
        for (; it >= static_cast<int>(input_reads.start_inds[best_M[N - 1].readInd]); it--)
            index_samples_more_needed[it] = 0;
        best_M.eraseLess(end_ind);
        best_M.eraseFront(N);

        if (end_ind < start_ind) end_ind = start_ind;

        SizedPriorityQueue<HelperRead>& tmp = best_M_prev;
        best_M_prev = best_M;
        best_M = tmp;
    }

    return obtain_sequence(input_reads, ret_vec);
}

std::unique_ptr<qmcp::Solution> qmcp::CpuCrawlingSolver::obtain_sequence(
    const bam_api::SOAPairedReads& sequence, std::vector<int64_t>& solution)
{
    auto reduced_reads = std::make_unique<Solution>();
    for (bam_api::ReadIndex read_id = 0; read_id < sequence.get_reads_count(); read_id++)
        if (solution[static_cast<int>(read_id)] > 0) reduced_reads->push_back(read_id);

    return reduced_reads;
}