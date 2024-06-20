#include "bam-api/soa_paired_reads.hpp"

#include "bam-api/aos_paired_reads.hpp"

// SOAPairedReads methods
bam_api::ReadQuality bam_api::SOAPairedReads::get_quality(ReadIndex index) const {
    return qualities[index];
}

void bam_api::SOAPairedReads::set_quality(ReadIndex index, ReadQuality quality) {
    qualities[index] = quality;
}

bam_api::Read bam_api::SOAPairedReads::get_read_by_index(ReadIndex index) const {
    return Read(ids[index], start_inds[index], end_inds[index], qualities[index],
                seq_lengths[index], is_first_reads[index]);
}

bam_api::ReadIndex bam_api::SOAPairedReads::get_reads_count() const { return ids.size(); }

void bam_api::SOAPairedReads::clear() {
    ids.clear();
    start_inds.clear();
    end_inds.clear();
    qualities.clear();
    seq_lengths.clear();
    is_first_reads.clear();
}

bam_api::SOAPairedReads& bam_api::SOAPairedReads::from(const AOSPairedReads& aos) {
    this->reserve(aos.get_reads_count());
    this->clear();

    for (const auto& read : aos.reads) {
        this->push_back(read);
    }

    this->ref_genome_length = aos.ref_genome_length;

    return *this;
}
