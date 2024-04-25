#include "../include/bam-api/bam_api.hpp"

#include <htslib/sam.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <ostream>
#include <vector>

#include "../include/bam-api/bam_paired_reads.hpp"

bam_api::SOAPairedReads bam_api::BamApi::read_bam_soa(
    const std::filesystem::path& filepath) {
    AOSPairedReads paired_reads = read_bam_aos(filepath);
    return paired_reads.to_soa();
}

bam_api::AOSPairedReads bam_api::BamApi::read_bam_aos(
    const std::filesystem::path& filepath) {
    AOSPairedReads paired_reads;

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
        std::exit(EXIT_FAILURE);
    }

    // read header
    in_samhdr = sam_hdr_read(infile);
    if (!in_samhdr) {
        std::cerr << "Failed to read header from file!" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    paired_reads.ref_genome_length = in_samhdr->target_len[0];

    int64_t id = 0;
    std::map<std::string, Read> read_map;
    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
        Read current_read{.id = id,
                          .start_ind = bamdata->core.pos,
                          .end_ind = bamdata->core.pos + bamdata->core.l_qseq,
                          .quality = bamdata->core.qual,
                          .qname = bam_get_qname(bamdata),  // NOLINT
                          .is_first_read = static_cast<bool>(
                              bamdata->core.flag & BAM_FREAD1)};

        auto read_one_iterator = read_map.find(current_read.qname);
        if (read_one_iterator != read_map.end()) {
            paired_reads.read_pair_map[read_one_iterator->second.id] = id;
            paired_reads.read_pair_map[id] = read_one_iterator->second.id;
            paired_reads.reads.push_back(read_one_iterator->second);
            paired_reads.reads.push_back(current_read);
        } else {
            read_map.insert({current_read.qname, current_read});
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
    return paired_reads;
}

int32_t bam_api::BamApi::write_bam_aos(
    const std::filesystem::path& input_filepath,
    const std::filesystem::path& output_filepath, AOSPairedReads& paired_reads) {
    SOAPairedReads soa = paired_reads.to_soa();
    return write_bam_soa(input_filepath, output_filepath, soa);
}

int32_t bam_api::BamApi::write_bam_soa(
    const std::filesystem::path& input_filepath,
    const std::filesystem::path& output_filepath, SOAPairedReads& paired_reads) {
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
    outfile = sam_open(output_filepath.c_str(), "wb");
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

    int64_t id = 0;
    int32_t current_read_i = 0;
    uint64_t read_count = paired_reads.ids.size();
    std::vector<int64_t> ids(read_count);
    std::copy(paired_reads.ids.begin(), paired_reads.ids.end(), ids.begin());
    std::sort(ids.begin(), ids.end());

    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0 &&
           read_count > current_read_i) {
        if (id == ids[current_read_i]) {
            if (sam_write1(outfile, in_samhdr, bamdata) < 0) {
                std::cerr << "Can't write line to bam file: " << output_filepath
                          << std::endl;
                std::exit(EXIT_FAILURE);
            }
            current_read_i++;
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
    if (outfile) {
        sam_close(outfile);
    }
    if (bamdata) {
        bam_destroy1(bamdata);
    }

    return current_read_i;
}
