#ifndef REGION_API_HPP
#define REGION_API_HPP

#include <unordered_set>
#include <vector>

#include "bam-api/aos_paired_reads.hpp"
#include "bam-api/paired_reads.hpp"
#include "bam-api/soa_paired_reads.hpp"

namespace bam_api {
class RegionApi {
   public:
    explicit RegionApi(const AOSPairedReads& paired_reads, bool is_paired = true);
    explicit RegionApi(const SOAPairedReads& paired_read, bool is_paired = true);

    const AOSPairedReads& get_paired_reads_aos();
    const SOAPairedReads& get_paired_reads_soa();
    const PairedReads& get_paired_reads() const;
    std::unordered_set<ReadIndex> add_pseudo_reads();

    std::vector<ReadIndex> find_pairs(const std::vector<ReadIndex>& ids) const;

    // testing purposes
    std::vector<uint32_t> find_input_cover() const;
    std::vector<uint32_t> find_filtered_cover(const std::vector<ReadIndex>& active_ids) const;

   private:
    bool is_paired_ = true;
    SOAPairedReads soa_paired_reads_;
    bool is_soa_loaded_ = false;
    AOSPairedReads aos_paired_reads_;
    bool is_aos_loaded_ = false;
};

}  // namespace bam_api

#endif
