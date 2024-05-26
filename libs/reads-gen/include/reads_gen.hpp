#ifndef READS_GEN_HPP
#define READS_GEN_HPP
#include <cstdint>
#include <functional>
#include <random>

#include "bam-api/bam_paired_reads.hpp"

namespace reads_gen {

// Basing on the values of the `dist_func` for x in [0, 1] the distribution is created. The
// distribution is used during random indices generations.
bam_api::AOSPairedReads rand_reads(std::mt19937& generator, bam_api::ReadIndex pairs_count,
                                   bam_api::Index genome_length, uint32_t read_length,
                                   const std::function<double(double)>& dist_func);

bam_api::AOSPairedReads rand_reads_uniform(std::mt19937& generator, bam_api::ReadIndex pairs_count,
                                           bam_api::Index genome_length, uint32_t read_length);
}  // namespace reads_gen

#endif