#include "../include/bam-api/bam_api.hpp"

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <map>
#include <optional>
#include <ostream>
#include <vector>

#include "../include/bam-api/bam_paired_reads.hpp"

bam_api::SOAPairedReads bam_api::BamApi::read_bam_soa(
    const std::filesystem::path& filepath) {
    SOAPairedReads paired_reads;
    read_bam(filepath, paired_reads);
    return paired_reads;
}

bam_api::AOSPairedReads bam_api::BamApi::read_bam_aos(
    const std::filesystem::path& filepath) {
    AOSPairedReads paired_reads;
    read_bam(filepath, paired_reads);
    return paired_reads;
}

void bam_api::BamApi::read_bam(const std::filesystem::path& filepath,
                               PairedReads& paired_reads) {
    sam_hdr_t* in_samhdr = NULL;
    samFile* infile = NULL;
    int ret_r = 0;
    bam1_t* bamdata = NULL;

    bamdata = bam_init1();
    if (!bamdata) {
        std::cerr << "Failed to allocate data memory!" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // open input file
    infile = sam_open(filepath.c_str(), "r");
    if (!infile) {
        std::cerr << "Could not open " << filepath << std::endl;
        sam_close(infile);
        std::exit(EXIT_FAILURE);
    }

    // read header
    in_samhdr = sam_hdr_read(infile);
    if (!in_samhdr) {
        std::cerr << "Failed to read header from file!" << std::endl;
        sam_close(infile);
        std::exit(EXIT_FAILURE);
    }

    // Load the index
    hts_idx_t* idx = sam_index_load(infile, filepath.c_str());
    if (idx) {
        uint64_t mapped_seq_c, unmapped_seq_c;
        hts_idx_get_stat(idx, 0, &mapped_seq_c, &unmapped_seq_c);
        paired_reads.reserve(mapped_seq_c);
        hts_idx_destroy(idx);
    }

    paired_reads.ref_genome_length = in_samhdr->target_len[0];

    ReadIndex id = 0;
    std::map<std::string, Read> read_map;
    std::string current_qname;
    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
        Read current_read{
            .id = id,
            .start_ind = static_cast<ReadIndex>(bamdata->core.pos),
            .end_ind = static_cast<ReadIndex>(bamdata->core.pos +
                                              bamdata->core.l_qseq),
            .quality = bamdata->core.qual,
            .is_first_read =
                static_cast<bool>(bamdata->core.flag & BAM_FREAD1)};
        current_qname = bam_get_qname(bamdata);  // NOLINT

        auto read_one_iterator = read_map.find(current_qname);
        if (read_one_iterator != read_map.end()) {
            // if first read of pair is under index i, second is under i+1
            if (current_read.is_first_read) {
                paired_reads.push_back(current_read);
                paired_reads.push_back(read_one_iterator->second);
            } else {
                paired_reads.push_back(read_one_iterator->second);
                paired_reads.push_back(current_read);
            }

            paired_reads.read_pair_map[read_one_iterator->second.id] =
                current_read.id;
            paired_reads.read_pair_map.push_back(read_one_iterator->second.id);
        } else {
            read_map.insert({current_qname, current_read});
            paired_reads.read_pair_map.push_back(std::nullopt);
        }

        id++;
    }

    if (ret_r >= 0)
        std::cerr << "Failed to read bam file (sam_read1 error code:" << ret_r
                  << ")" << std::endl;

    // cleanup
    if (in_samhdr) {
        sam_hdr_destroy(in_samhdr);
    }
    if (infile) {
        sam_close(infile);
    }
    if (bamdata) {
        bam_destroy1(bamdata);
    }
}

void bam_api::BamApi::write_sam(const std::filesystem::path& input_filepath,
                                const std::filesystem::path& output_filepath,
                                std::vector<ReadIndex>& read_ids,
                                bool use_bam) {
    sam_hdr_t* in_samhdr = NULL;
    samFile* infile = NULL;
    samFile* outfile = NULL;
    int ret_r = 0;
    bam1_t* bamdata = NULL;

    bamdata = bam_init1();
    if (!bamdata) {
        std::cerr << "Failed to allocate data memory!" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // open input file
    infile = sam_open(input_filepath.c_str(), "r");
    if (!infile) {
        std::cerr << "Could not open " << input_filepath << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // open output file
    outfile = sam_open(output_filepath.c_str(), use_bam ? "wb" : "w");
    if (!outfile) {
        std::cerr << "Could not open " << output_filepath << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // read header
    in_samhdr = sam_hdr_read(infile);
    if (!in_samhdr) {
        std::cerr << "Failed to read header from file!" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if (sam_hdr_write(outfile, in_samhdr) < 0) {
        std::cerr << "Can't write header to bam file: " << output_filepath
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }

    ReadIndex id = 0;
    std::sort(read_ids.begin(), read_ids.end());
    auto current_read_i = read_ids.begin();

    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0 &&
           current_read_i != read_ids.end()) {
        if (id == *current_read_i) {
            if (sam_write1(outfile, in_samhdr, bamdata) < 0) {
                std::cerr << "Can't write line to bam file: " << output_filepath
                          << std::endl;
                std::exit(EXIT_FAILURE);
            }

            current_read_i++;
        }

        id++;
    }

    if (current_read_i != read_ids.end() && ret_r >= 0)
        std::cerr << "Failed to read bam file (sam_read1 error code:" << ret_r
                  << ")" << std::endl;

    // cleanup
    if (in_samhdr) {
        sam_hdr_destroy(in_samhdr);
    }
    if (infile) {
        sam_close(infile);
    }
    if (outfile) {
        sam_close(outfile);
    }
    if (bamdata) {
        bam_destroy1(bamdata);
    }
}
