#include "../include/bam-api/bam_api.hpp"

#include <htslib/sam.h>

#include <iostream>
#include <ostream>

#include "bam-api/bam_data.hpp"

void bam_api::BamApi::test_func() { std::cout << "BamApi TEST\n"; }

bam_api::AOSPairedReads bam_api::BamApi::read_bam_aos(std::string filepath) {
    sam_hdr_t* in_samhdr = NULL;
    samFile* infile = NULL;
    int ret_r = 0;
    bam1_t* bamdata = NULL;

    AOSPairedReads ret_paired_reads;

    bamdata = bam_init1();
    if (!bamdata) {
        std::cerr << "Failed to allocate data memory!" << std::endl;
        // mytkom/TODO: log fatal
    }

    // open input file
    infile = sam_open(filepath.data(), "r");
    if (!infile) {
        std::cerr << "Could not open " << filepath << std::endl;
        // mytkom/TODO: log fatal
    }

    // read header
    in_samhdr = sam_hdr_read(infile);
    if (!in_samhdr) {
        std::cerr << "Failed to read header from file!" << std::endl;
        // mytkom/TODO: log fatal
    }

    int64_t id = 0;
    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
        ret_paired_reads.reads.push_back(
            Read{.id = id,
                 .start_ind = bamdata->core.pos,
                 .end_ind = bamdata->core.pos + bamdata->core.l_qseq,
                 .quality = bamdata->core.qual,
                 .qname = bam_get_qname(bamdata),
                 .is_first_read =
                     static_cast<bool>(bamdata->core.flag & BAM_FREAD1)});

        // mytkom/TODO: find paired sequence
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

    return ret_paired_reads;
}

bam_api::SOAPairedReads bam_api::BamApi::read_bam_soa(std::string filepath) {
    sam_hdr_t* in_samhdr = NULL;
    samFile* infile = NULL;
    int ret_r = 0;
    bam1_t* bamdata = NULL;

    SOAPairedReads ret_paired_reads;

    bamdata = bam_init1();
    if (!bamdata) {
        std::cerr << "Failed to allocate data memory!" << std::endl;
        // mytkom/TODO: log fatal
    }

    // open input file
    infile = sam_open(filepath.data(), "r");
    if (!infile) {
        std::cerr << "Could not open " << filepath << std::endl;
        // mytkom/TODO: log fatal
    }

    // read header
    in_samhdr = sam_hdr_read(infile);
    if (!in_samhdr) {
        std::cerr << "Failed to read header from file!" << std::endl;
        // mytkom/TODO: log fatal
    }

    int64_t id = 0;
    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
        ret_paired_reads.ids.push_back(id);
        ret_paired_reads.start_inds.push_back(bamdata->core.pos);
        ret_paired_reads.end_inds.push_back(bamdata->core.pos +
                                            bamdata->core.l_qseq);
        ret_paired_reads.qualities.push_back(bamdata->core.qual);
        ret_paired_reads.qnames.push_back(bam_get_qname(bamdata));
        ret_paired_reads.is_first_read.push_back(static_cast<bool>(bamdata->core.flag & BAM_FREAD1));
        // mytkom/TODO: handle paired seq
        // mytkom/TODO: ref genome length
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

    return ret_paired_reads;
}
