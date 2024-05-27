#include "../include/reads_gen.hpp"

#include <iostream>
#include <random>

#include "bam-api/bam_paired_reads.hpp"

bam_api::AOSPairedReads reads_gen::rand_reads(std::mt19937& generator,
                                              bam_api::ReadIndex pairs_count,
                                              bam_api::Index genome_length, uint32_t read_length,
                                              const std::function<double(double)>& dist_func) {
    uint32_t starts_count = genome_length - read_length + 1;
    std::vector<double> dist_data(starts_count, 0);
    double sum = 0;
    for (uint32_t i = 0; i < starts_count; ++i) {
        dist_data[i] = dist_func(static_cast<double>(i) / static_cast<double>(starts_count - 1));
        if (dist_data[i] < 0.0) {
            dist_data[i] = 0.0;
        }

        sum += dist_data[i];
    }

    for (uint32_t i = 0; i < starts_count; ++i) {
        dist_data[i] /= sum;
    }

    std::discrete_distribution<> dist(dist_data.begin(), dist_data.end());

    bam_api::AOSPairedReads result;
    result.ref_genome_length = genome_length;
    for (bam_api::ReadIndex i = 0; i < 2 * pairs_count; i += 2) {
        bam_api::Index first = dist(generator);
        bam_api::Index second = dist(generator);

        if (first > second) {
            std::swap(first, second);
        }

        if (first > genome_length - read_length * 2 && second > genome_length - read_length * 2) {
            first = genome_length - read_length * 2;
            second = genome_length - read_length;
        } else if (first + read_length > second) {
            second = first + read_length;
        }

        result.reads.push_back(bam_api::Read(i, first, first + read_length - 1, 0, true));
        result.reads.push_back(bam_api::Read(i + 1, second, second + read_length - 1, 0, false));
        result.bam_id_to_read_index.push_back(i);
        result.bam_id_to_read_index.push_back(i + 1);
    }

    return result;
}

bam_api::AOSPairedReads reads_gen::rand_reads_uniform(std::mt19937& generator,
                                                      bam_api::ReadIndex pairs_count,
                                                      bam_api::Index genome_length,
                                                      uint32_t read_length) {
    std::uniform_int_distribution<> dist_first(
        0, static_cast<int32_t>(genome_length - 2 * read_length));
    std::uniform_int_distribution<> dist_second(0,
                                                static_cast<int32_t>(genome_length - read_length));

    bam_api::AOSPairedReads result;
    result.ref_genome_length = genome_length;
    for (bam_api::BAMReadId i = 0; i < 2 * pairs_count; i += 2) {
        bam_api::Index first = dist_first(generator);
        bam_api::Index second = dist_second(generator);

        if (first > second) {
            std::swap(first, second);
        }

        if (first + read_length > second) {
            second = first + read_length;
        }

        result.reads.push_back(bam_api::Read(i, first, first + read_length - 1, 0, true));
        result.reads.push_back(bam_api::Read(i + 1, second, second + read_length - 1, 0, false));
        result.bam_id_to_read_index.push_back(i);
        result.bam_id_to_read_index.push_back(i + 1);
    }

    return result;
}
