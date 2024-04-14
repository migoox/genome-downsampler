#include "../include/bam-api/bam_api.hpp"

#include <htslib/sam.h>

#include <cstdlib>
#include <iostream>
#include <map>
#include <ostream>

#include "../include/bam-api/bam_paired_reads.hpp"

void bam_api::BamApi::read_bam(bam_api::PairedReads& paired_reads,
                               const std::filesystem::path& filepath) {
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
    int temp = 0;
    std::map<std::filesystem::path, int64_t> read_map;
    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
        Read current_read{.id = id,
                          .start_ind = bamdata->core.pos,
                          .end_ind = bamdata->core.pos + bamdata->core.l_qseq,
                          .quality = bamdata->core.qual,
                          .qname = bam_get_qname(bamdata),
                          .is_first_read = static_cast<bool>(
                              bamdata->core.flag & BAM_FREAD1)};

        auto read_one_iterator = read_map.find(current_read.qname);
        if (read_one_iterator != read_map.end()) {
            paired_reads.read_pair_map[read_one_iterator->second] = id;
            paired_reads.read_pair_map.push_back(read_one_iterator->second);
            temp++;
        } else {
            paired_reads.read_pair_map.push_back(-1);
            read_map.insert({current_read.qname, current_read.id});
        }

        paired_reads.push_back(current_read);
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

bam_api::SOAPairedReads bam_api::BamApi::read_bam_soa(
    const std::filesystem::path& filepath) {
    SOAPairedReads ret_paired_reads;
    read_bam(ret_paired_reads, filepath);
    return ret_paired_reads;
}

bam_api::AOSPairedReads bam_api::BamApi::read_bam_aos(
    const std::filesystem::path& filepath) {
    AOSPairedReads ret_paired_reads;
    read_bam(ret_paired_reads, filepath);
    return ret_paired_reads;
}
