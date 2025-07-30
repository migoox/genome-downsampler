#include "qmcp-solver/cpu_crawling_solver.hpp"

#include <xray/xray_interface.h>

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

void printReads1(const bam_api::SOAPairedReads& input_sequence)
{
    for (int x = 0; x < input_sequence.ref_genome_length; x++)
    {
        for (int y = 0; y < input_sequence.get_reads_count(); y++)
        {
            if (x >= input_sequence.start_inds[y] && x <= input_sequence.end_inds[y])
                std::cout << "1 ";
            else
                std::cout << "0 ";
        }
        std::cout << "\n";
    }
}

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
    uint32_t readInd;
    HelperRead(uint32_t readInd, bam_api::Index endInd) : endInd(endInd), readInd(readInd) {}
    bool operator<(const HelperRead& other) const { return endInd < other.endInd; }

    bool operator>(const HelperRead& other) const { return endInd > other.endInd; }

    bool operator<=(const HelperRead& other) const { return endInd <= other.endInd; }

    bool operator>=(const HelperRead& other) const { return endInd >= other.endInd; }

    bool operator==(const HelperRead& other) const { return endInd == other.endInd; }

    bool operator!=(const HelperRead& other) const { return endInd != other.endInd; }

    bool operator<(const int other) const { return endInd < other; }

    bool operator>(const int other) const { return endInd > other; }

    bool operator<=(const int other) const { return endInd <= other; }

    bool operator>=(const int other) const { return endInd >= other; }

    bool operator==(const int other) const { return endInd == other; }

    bool operator!=(const int other) const { return endInd != other; }

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
struct SizedPriorityQueue_vectorShiftBased
{
    std::vector<T> buffer;
    uint32_t M;
    explicit SizedPriorityQueue_vectorShiftBased(uint32_t M) : M(M) { buffer.reserve(M); }
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
    uint32_t size() const
    {
        return buffer.size();
    }

    void eraseLess(int num)
    {
        if (buffer[buffer.size() - 1] >= num)
            return;
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

    T& operator[](size_t index)
    {
        return buffer[index];
    }
};

#define TIME_SOLUTION
#ifdef TIME_SOLUTION
#include <chrono>
using namespace std::chrono;
#endif

// test time 58s
void findBestSamples_tryFindShorter(uint32_t max_coverage, const bam_api::SOAPairedReads& input_sequence, uint32_t start_ind, uint32_t end_ind,
                                    std::vector<uint32_t>& start_inex_count_pref_sum, std::vector<uint32_t>& reads_nr_sorted_by_start_ind,
                                    SizedPriorityQueue_vectorShiftBased<HelperRead>& best_M, SizedPriorityQueue_vectorShiftBased<HelperRead>& best_M_prev)
{
    for (uint32_t x = start_ind; x <= end_ind; x++)
    {
        int subseq_last_sorted_ind = start_inex_count_pref_sum[x + 1] - 1;
        int offset = static_cast<int>(std::min(max_coverage, start_inex_count_pref_sum[x + 1] - start_inex_count_pref_sum[x]));

        for (int sorted_ind = subseq_last_sorted_ind; sorted_ind >= 0 &&
                                                      sorted_ind > subseq_last_sorted_ind - offset;
             sorted_ind--)
        {
            uint32_t ind = reads_nr_sorted_by_start_ind[sorted_ind];
            best_M.tryAdd(HelperRead(ind, input_sequence.end_inds[ind]));
        }
    }

    for (int x = 0; x < best_M_prev.size(); x++)
        best_M.tryAdd(best_M_prev[x]);
}

// test time 38s
void findBestSamples(uint32_t max_coverage, const bam_api::SOAPairedReads& input_sequence, uint32_t start_ind, uint32_t end_ind,
                     std::vector<uint32_t>& start_inex_count_pref_sum, std::vector<uint32_t>& reads_nr_sorted_by_start_ind,
                     SizedPriorityQueue_vectorShiftBased<HelperRead>& best_M, SizedPriorityQueue_vectorShiftBased<HelperRead>& best_M_prev)
{
    bool found = false;
    uint32_t ind_offset = 0;
    do
    {
        found = false;
        for (uint32_t x = start_ind; x <= end_ind; x++)
        {
            int64_t subseq_sorted_ind = static_cast<int64_t>(start_inex_count_pref_sum[x + 1]) - 1 - ind_offset;
            if (subseq_sorted_ind < static_cast<int64_t>(start_inex_count_pref_sum[x]))
                continue;
            uint32_t ind = reads_nr_sorted_by_start_ind[subseq_sorted_ind];
            found |= best_M.tryAdd(HelperRead(ind, input_sequence.end_inds[ind]));
        }
        ind_offset++;
    } while (found);

    for (int x = 0; x < best_M_prev.size(); x++)
        best_M.tryAdd(best_M_prev[x]);
}

template <typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
struct CountSortHleper
{
    std::vector<T> counts, prefSum;
    std::map<T, T> prefSumIndexOutIndMapping;
    T M,
        minVal, maxVal, lowLimitter, total, prefSumStart;
    explicit CountSortHleper()
    {
        minVal = std::numeric_limits<T>::max();
        maxVal = 0;
        total = 0;
    }

    void setM(T M)
    {
        counts.resize(M * 10);  // TODO test
        prefSum.resize(M * 10);
        this->M = M;
    }

    void setLowLimit(T lowLimit)
    {
        lowLimitter = lowLimit;
    }

    bool tryAddVal(T val)
    {
        if (total >= M && val <= minVal)
            return false;
        minVal = std::min(val, minVal);
        maxVal = std::max(val, maxVal);
        total++;
        counts[val - lowLimitter]++;
    }

    void calculatePrefSums()
    {
        T ind;
        T newMappingVal = 0;
        prefSum[maxVal - lowLimitter] = counts[maxVal - lowLimitter];  // TODO should not be 0, but deal with this case
        prefSumIndexOutIndMapping[maxVal - lowLimitter] = newMappingVal++;
        for (ind = maxVal - lowLimitter - 1; ind >= 0; ind--)
        {
            if (prefSum[ind + 1] >= M)
            {
                prefSum[ind + 1] = M;
                break;
            }
            prefSum[ind] = counts[ind] + prefSum[ind + 1];
            if (counts[ind])
                prefSumIndexOutIndMapping[ind] = newMappingVal++;
        }
        prefSumStart = ind + 1;
        for (; ind >= 0; ind--)
            counts[ind] = 0;
    }
};

void findBestSamples_countSort(uint32_t max_coverage, const bam_api::SOAPairedReads& input_sequence, uint32_t start_ind, uint32_t end_ind,
                               std::vector<uint32_t>& start_inex_count_pref_sum, std::vector<uint32_t>& reads_nr_sorted_by_start_ind,
                               SizedPriorityQueue_vectorShiftBased<HelperRead>& best_M, SizedPriorityQueue_vectorShiftBased<HelperRead>& best_M_prev)
{
    CountSortHleper<uint32_t> countSortHelper;
    uint32_t M = best_M.M;
    if (countSortHelper.M != M)
        countSortHelper.setM(M);

    countSortHelper.setLowLimit(start_ind);
    bool found = false;
    uint32_t ind_offset = 0;

    // fill count array
    do
    {
        found = false;
        for (uint32_t x = start_ind; x <= end_ind; x++)
        {
            int64_t subseq_sorted_ind = static_cast<int64_t>(start_inex_count_pref_sum[x + 1]) - 1 - ind_offset;
            if (subseq_sorted_ind < static_cast<int64_t>(start_inex_count_pref_sum[x]))
                continue;
            uint32_t ind = reads_nr_sorted_by_start_ind[subseq_sorted_ind];
            found |= countSortHelper.tryAddVal(input_sequence.end_inds[ind]);
        }
        ind_offset++;
    } while (found);

    for (int x = 0; x < best_M_prev.size(); x++)
        countSortHelper.tryAddVal(best_M_prev[x].endInd);

    countSortHelper.calculatePrefSums();

    // collect solution
    std::vector<HelperRead> sol;
    sol.reserve(M);
    int count = 0;

    found = false;
    ind_offset = 0;
    do
    {
        found = false;
        for (uint32_t x = start_ind; x <= end_ind; x++)
        {
            int64_t subseq_sorted_ind = static_cast<int64_t>(start_inex_count_pref_sum[x + 1]) - 1 - ind_offset;
            if (subseq_sorted_ind < static_cast<int64_t>(start_inex_count_pref_sum[x]))
                continue;
            uint32_t ind = reads_nr_sorted_by_start_ind[subseq_sorted_ind];
            uint32_t found_end_ind = input_sequence.end_inds[ind];
            if (countSortHelper.prefSum[found_end_ind] > 0 && countSortHelper.counts[found_end_ind] > 0)
            {
                found = true;
                sol[countSortHelper.prefSumIndexOutIndMapping[found_end_ind]] = HelperRead(ind, found_end_ind);
                count++;
            }
        }
        ind_offset++;
    } while (found);

    for (int x = 0; x < best_M_prev.size(); x++)
        if (countSortHelper.prefSum[best_M_prev[x].endInd] > 0 && countSortHelper.counts[best_M_prev[x].endInd] > 0)
        {
            found = true;
            sol[countSortHelper.prefSumIndexOutIndMapping[best_M_prev[x].endInd]] = best_M_prev[x];
            count++;
        }

    for (int x = 0; x < count; x++)
        best_M.tryAdd(sol[x]);
}

// #define DEBUG_PRINT
// ALGORYTM PEŁZAJĄCY
[[clang::xray_always_instrument]] __attribute__((noinline))
std::unique_ptr<qmcp::Solution>
qmcp::CpuCrawlingSolver::solve(uint32_t max_coverage,
                               bam_api::BamApi& bam_api)
{
    const bam_api::SOAPairedReads& input_sequence = bam_api.get_paired_reads_soa();
    const uint32_t seq_size = input_sequence.ref_genome_length;
    const uint32_t reads_count = input_sequence.get_reads_count();

    uint32_t start_ind = 0;
    uint32_t end_ind = 0;
    std::vector<uint32_t> index_samples_more_needed;
    std::vector<uint32_t> start_inex_count_pref_sum;
    std::vector<uint32_t> reads_nr_sorted_by_start_ind;
    SizedPriorityQueue_vectorShiftBased<HelperRead> best_M1(max_coverage), best_M2(max_coverage);
    SizedPriorityQueue_vectorShiftBased<HelperRead>&best_M = best_M1, &best_M_prev = best_M2;
    std::vector<int64_t> ret_vec(reads_count, 0);

    index_samples_more_needed.resize(seq_size, max_coverage);

    // this can be removed by implementing sort by count
    start_inex_count_pref_sum.resize(seq_size + 1, 0);
    reads_nr_sorted_by_start_ind.resize(reads_count);

    for (int x = 0; x < reads_count; x++)
    {
        start_inex_count_pref_sum[input_sequence.start_inds[x] + 1]++;
        reads_nr_sorted_by_start_ind[x] = x;
    }

#ifdef DEBUG_PRINT
    std::cout << "reads:\n";
    printReads1(input_sequence);

    std::cout << "startIndCount:  ";
    printVector(start_inex_count_pref_sum);
#endif

    for (int x = 1; x < start_inex_count_pref_sum.size(); x++)
        start_inex_count_pref_sum[x] += start_inex_count_pref_sum[x - 1];
    start_inex_count_pref_sum.push_back(start_inex_count_pref_sum[start_inex_count_pref_sum.size() - 1]);

#ifdef DEBUG_PRINT
    std::cout << "startIndCountPrefSum:  ";
    printVector(start_inex_count_pref_sum);
#endif

    std::sort(reads_nr_sorted_by_start_ind.begin(), reads_nr_sorted_by_start_ind.end(),
              [&](int ind1, int ind2)
              {
                  if (input_sequence.start_inds[ind1] == input_sequence.start_inds[ind2])
                      return input_sequence.end_inds[ind1] < input_sequence.end_inds[ind2];
                  return input_sequence.start_inds[ind1] < input_sequence.start_inds[ind2];
              });

#ifdef DEBUG_PRINT
    std::cout << "sortedReadsInd:  ";
    printVector(reads_nr_sorted_by_start_ind);
#endif

    while (start_ind < seq_size)
    {
#ifdef DEBUG_PRINT
        std::cout << "startInd: " << start_ind << " endInd: " << end_ind << "\n";
        std::cout << "samplesMoreNeededForIndex:  ";
        printVector(index_samples_more_needed);
#endif
        // find best samples
        best_M.clear();
        // findBestSamples_tryFindShorter(max_coverage, input_sequence, start_ind, end_ind,
        //                                start_inex_count_pref_sum, reads_nr_sorted_by_start_ind,
        //                                best_M, best_M_prev);
        findBestSamples(max_coverage, input_sequence, start_ind, end_ind,
                        start_inex_count_pref_sum, reads_nr_sorted_by_start_ind,
                        best_M, best_M_prev);

        if (best_M.size() == 0)
        {
            start_ind = end_ind + 1;
            end_ind++;
            continue;
        }
        // N samples to take
        uint32_t N = std::min(index_samples_more_needed[end_ind], best_M.size());
        // update coverage array and ret vec
#ifdef DEBUG_PRINT
        std::cout << "bestM:   ";
        printVector(best_M.buffer);
#endif
        ret_vec[best_M[0].readInd] = 1;
        int covered_by = 1;
        for (int x = 1; x < N; x++)
        {
            ret_vec[best_M[x].readInd] = 1;
            for (int y = best_M[x - 1].endInd; y > best_M[x].endInd; y--)
                index_samples_more_needed[y] -= covered_by;
            covered_by++;
        }

#ifdef DEBUG_PRINT
        std::cout << "samplesMoreNeededForIndex:  ";
        printVector(index_samples_more_needed);

#endif

#ifdef DEBUG_PRINT
        std::cout << "retVec:   ";
        printVector(ret_vec);

#endif
        // update start, end
        start_ind = end_ind + 1;
        int it;
        for (it = best_M[N - 1].endInd; N < index_samples_more_needed[it] &&
                                        it > input_sequence.start_inds[best_M[N - 1].readInd];
             it--)
            index_samples_more_needed[it] -= N;
        end_ind = it + 1;
        for (; it >= static_cast<int>(input_sequence.start_inds[best_M[N - 1].readInd]); it--)
            index_samples_more_needed[it] = 0;
        best_M.eraseLess(end_ind);
        best_M.eraseFront(N);

#ifdef DEBUG_PRINT
        std::cout << "samplesMoreNeededForIndex:  ";
        printVector(index_samples_more_needed);

        std::cout << "bestM:   ";
        printVector(best_M.buffer);
#endif
        if (end_ind < start_ind)
            end_ind = start_ind;

        SizedPriorityQueue_vectorShiftBased<HelperRead>& tmp = best_M_prev;
        best_M_prev = best_M;
        best_M = tmp;
    }

    // printVector(ret_vec);
    // std::reverse(ret_vec.begin(), ret_vec.end());
    return obtain_sequence(input_sequence, ret_vec);
}

std::unique_ptr<qmcp::Solution> qmcp::CpuCrawlingSolver::obtain_sequence(
    const bam_api::SOAPairedReads& sequence, std::vector<int64_t>& solution)
{
    auto reduced_reads = std::make_unique<Solution>();
    for (bam_api::ReadIndex read_id = 0; read_id < sequence.get_reads_count(); read_id++)
        if (solution[static_cast<int>(read_id)] > 0) reduced_reads->push_back(read_id);

    return reduced_reads;
}