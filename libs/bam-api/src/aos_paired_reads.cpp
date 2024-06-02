#include "bam-api/aos_paired_reads.hpp"
#include "bam-api/soa_paired_reads.hpp"

// AOSPairedReads methods
void bam_api::AOSPairedReads::push_back(Read&& read) {
    reads.push_back(std::move(read));
}

void bam_api::AOSPairedReads::push_back(const Read& read) {
    reads.push_back(read);
}

bam_api::Read bam_api::AOSPairedReads::get_read_by_index(ReadIndex index) const {
    return reads[index];
}

bam_api::ReadIndex bam_api::AOSPairedReads::get_reads_count() const {
    return reads.size();
}

void bam_api::AOSPairedReads::reserve(size_t size) {
    reads.reserve(size);
}

void bam_api::AOSPairedReads::clear() {
    reads.clear();
}

bam_api::AOSPairedReads& bam_api::AOSPairedReads::from(const SOAPairedReads& soa) {
    this->reserve(soa.get_reads_count());
    this->clear();

    for (ReadIndex i = 0; i < soa.get_reads_count(); ++i) {
        this->push_back(soa.get_read_by_index(i));
    }

    this->ref_genome_length = soa.ref_genome_length;
    this->read_pair_map = soa.read_pair_map;
    this->bam_id_to_read_index = soa.bam_id_to_read_index;

    return *this;
}

