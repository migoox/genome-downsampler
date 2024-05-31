#ifndef READS_GEN_HPP
#define READS_GEN_HPP
#include <cstdint>
#include <functional>
#include <random>

#include "bam-api/paired_reads.hpp"

namespace reads_gen {

// Creates discrete distribution basing on the provided `dist_func`(for x in [0, 1] interval) which
// respresents histogram function. Then basing on the created distribution paired reads are
// generated.
bam_api::AOSPairedReads rand_reads(std::mt19937& generator, bam_api::ReadIndex pairs_count,
                                   bam_api::Index genome_length, uint32_t read_length,
                                   const std::function<double(double)>& dist_func);

// Generates random paired reads basing on the uniform distribution
bam_api::AOSPairedReads rand_reads_uniform(std::mt19937& generator, bam_api::ReadIndex pairs_count,
                                           bam_api::Index genome_length, uint32_t read_length);
}  // namespace reads_gen

#endif
