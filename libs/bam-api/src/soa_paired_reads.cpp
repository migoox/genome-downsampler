#include "bam-api/soa_paired_reads.hpp"
#include "bam-api/aos_paired_reads.hpp"

// SOAPairedReads methods
void bam_api::SOAPairedReads::push_back(Read&& read) {
    ids.emplace_back(read.bam_id);
    start_inds.emplace_back(read.start_ind);
    end_inds.emplace_back(read.end_ind);
    qualities.emplace_back(read.quality);
    seq_lengths.emplace_back(read.seq_length);
    is_first_reads.emplace_back(read.is_first_read);
}

void bam_api::SOAPairedReads::push_back(const Read& read) {
    ids.push_back(read.bam_id);
    start_inds.push_back(read.start_ind);
    end_inds.push_back(read.end_ind);
    qualities.push_back(read.quality);
    seq_lengths.push_back(read.seq_length);
    is_first_reads.push_back(read.is_first_read);
}

bam_api::Read bam_api::SOAPairedReads::get_read_by_index(ReadIndex index) const {
    return Read(ids[index], start_inds[index], end_inds[index], qualities[index],
                seq_lengths[index], is_first_reads[index]);
}

bam_api::ReadIndex bam_api::SOAPairedReads::get_reads_count() const { return ids.size(); }

void bam_api::SOAPairedReads::reserve(size_t size) {
    ids.reserve(size);
    start_inds.reserve(size);
    end_inds.reserve(size);
    qualities.reserve(size);
    seq_lengths.reserve(size);
    is_first_reads.reserve(size);
}

void bam_api::SOAPairedReads::clear() {
    ids.clear();
    start_inds.clear();
    end_inds.clear();
    qualities.clear();
    seq_lengths.clear();
    is_first_reads.clear();
}

bam_api::SOAPairedReads& bam_api::SOAPairedReads::operator=(const AOSPairedReads& aos) {
    this->reserve(aos.get_reads_count());
    this->clear();

    for (const auto& read : aos.reads) {
        this->push_back(read);
    }

    this->ref_genome_length = aos.ref_genome_length;
    this->read_pair_map = aos.read_pair_map;
    this->bam_id_to_read_index = aos.bam_id_to_read_index;

    return *this;
}
