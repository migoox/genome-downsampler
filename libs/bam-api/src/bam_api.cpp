#include "../include/bam-api/bam_api.hpp"

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <optional>
#include <vector>

#include "../include/bam-api/bam_paired_reads.hpp"
#include "logging/log.hpp"

bam_api::SOAPairedReads bam_api::BamApi::read_bam_soa(
    const std::filesystem::path& filepath, uint32_t min_seq_length,
    uint32_t min_mapq) {
    SOAPairedReads paired_reads;
    LOG(logging::kDebug) << "bam-api reading " << filepath << " to soa";
    read_bam(filepath, paired_reads, min_seq_length, min_mapq);
    return paired_reads;
}

bam_api::AOSPairedReads bam_api::BamApi::read_bam_aos(
    const std::filesystem::path& filepath, uint32_t min_seq_length,
    uint32_t min_mapq) {
    AOSPairedReads paired_reads;
    LOG(logging::kDebug) << "bam-api reading " << filepath << " to aos";
    read_bam(filepath, paired_reads, min_seq_length, min_mapq);
    return paired_reads;
}

void bam_api::BamApi::read_bam(const std::filesystem::path& filepath,
                               PairedReads& paired_reads,
                               uint32_t min_seq_length, uint32_t min_mapq) {
    sam_hdr_t* in_samhdr = NULL;
    samFile* infile = NULL;
    samFile* filtered_out_file = NULL;
    int ret_r = 0;
    bam1_t* bamdata = NULL;

    bamdata = bam_init1();
    if (!bamdata) {
        LOG(logging::kError) << "Failed to allocate data memory!";
        std::exit(EXIT_FAILURE);
    }

    // open input file
    infile = sam_open(filepath.c_str(), "r");
    if (!infile) {
        LOG(logging::kError) << "Could not open " << filepath;
        sam_close(infile);
        std::exit(EXIT_FAILURE);
    }

    auto filtered_out_path = filepath;
    filtered_out_path.replace_filename("filtered_out_sequences.bam");
    filtered_out_file = sam_open(filtered_out_path.c_str(), "wb");
    if (!filtered_out_file) {
        LOG(logging::kError) << "Could not open " << filtered_out_path;
        sam_close(filtered_out_file);
        std::exit(EXIT_FAILURE);
    }

    // read header
    in_samhdr = sam_hdr_read(infile);
    if (!in_samhdr) {
        LOG(logging::kError) << "Failed to read header from file!";
        sam_close(infile);
        std::exit(EXIT_FAILURE);
    }

    if (sam_hdr_write(filtered_out_file, in_samhdr) < 0) {
        LOG(logging::kError)
            << "Can't write header to bam file: " << filtered_out_path;
        std::exit(EXIT_FAILURE);
    }

    // Reserve memory for reads if index is present
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
    hts_pos_t rlen = 0;
    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
        rlen = bam_cigar2rlen(static_cast<int32_t>(bamdata->core.n_cigar),
                              bam_get_cigar(bamdata));
        Read current_read{
            .id = id,
            .start_ind = static_cast<GenomeIndex>(bamdata->core.pos),
            .end_ind = static_cast<GenomeIndex>(bamdata->core.pos + rlen - 1),
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
            if (bamdata->core.l_qseq >= min_seq_length &&
                bamdata->core.qual >= min_mapq)
                read_map.insert({current_qname, current_read});
            paired_reads.read_pair_map.push_back(std::nullopt);
        }

        id++;
    }

    if (ret_r >= 0)
        LOG(logging::kError)
            << "Failed to read bam file (sam_read1 error code:" << ret_r << ")";

    sam_hdr_destroy(in_samhdr);
    // Reopen infile to read iterate through it second time
    sam_close(infile);
    infile = sam_open(filepath.c_str(), "r");
    if (!infile) {
        LOG(logging::kError) << "Could not open " << filepath;
        sam_close(infile);
        std::exit(EXIT_FAILURE);
    }
    // read header
    in_samhdr = sam_hdr_read(infile);
    if (!in_samhdr) {
        LOG(logging::kError) << "Failed to read header from file!";
        sam_close(infile);
        std::exit(EXIT_FAILURE);
    }

    id = 0;
    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
        if (!paired_reads.read_pair_map[id]) {
            if (sam_write1(filtered_out_file, in_samhdr, bamdata) < 0) {
                LOG(logging::kError)
                    << "Can't write line to bam file: " << filtered_out_path;
                std::exit(EXIT_FAILURE);
            }
        }

        id++;
    }

    if (ret_r >= 0) {
        LOG(logging::kError)
            << "Failed to read bam file (sam_read1 error code:" << ret_r << ")";
    } else {
        LOG(logging::kDebug) << "Read bam file have been read correctly";
    }

    // cleanup
    sam_hdr_destroy(in_samhdr);
    sam_close(infile);
    sam_close(filtered_out_file);
    bam_destroy1(bamdata);
}

uint32_t bam_api::BamApi::write_bam(
    const std::filesystem::path& input_filepath,
    const std::filesystem::path& output_filepath,
    std::vector<ReadIndex>& read_ids) {
    LOG(logging::kDebug) << "bam-api writing" << read_ids.size() << " reads to "
                         << output_filepath << " on the basis of "
                         << input_filepath;

    sam_hdr_t* in_samhdr = NULL;
    samFile* infile = NULL;
    samFile* outfile = NULL;
    int ret_r = 0;
    bam1_t* bamdata = NULL;

    bamdata = bam_init1();
    if (!bamdata) {
        LOG(logging::kError) << "Failed to allocate data memory!";
        std::exit(EXIT_FAILURE);
    }

    // open input file
    infile = sam_open(input_filepath.c_str(), "r");
    if (!infile) {
        LOG(logging::kError) << "Could not open " << input_filepath;
        std::exit(EXIT_FAILURE);
    }

    // open output file
    std::string open_mode = output_filepath.extension() == ".bam" ? "wb" : "w";
    outfile = sam_open(output_filepath.c_str(), open_mode.c_str());
    if (!outfile) {
        LOG(logging::kError) << "Could not open " << output_filepath;
        std::exit(EXIT_FAILURE);
    }

    // read header
    in_samhdr = sam_hdr_read(infile);
    if (!in_samhdr) {
        LOG(logging::kError) << "Failed to read header from file!";
        std::exit(EXIT_FAILURE);
    }

    if (sam_hdr_write(outfile, in_samhdr) < 0) {
        LOG(logging::kError)
            << "Can't write header to bam file: " << output_filepath;
        std::exit(EXIT_FAILURE);
    }

    ReadIndex id = 0;
    std::sort(read_ids.begin(), read_ids.end());
    auto current_read_i = read_ids.begin();
    uint32_t reads_written = 0;

    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0 &&
           current_read_i != read_ids.end()) {
        if (id == *current_read_i) {
            if (sam_write1(outfile, in_samhdr, bamdata) < 0) {
                LOG(logging::kError)
                    << "Can't write line to bam file: " << output_filepath;
                std::exit(EXIT_FAILURE);
            }

            reads_written++;
            current_read_i++;
        }

        id++;
    }

    if (current_read_i != read_ids.end() && ret_r >= 0)
        LOG(logging::kError)
            << "Failed to read bam file (sam_read1 error code:" << ret_r << ")";

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

    LOG(logging::kDebug) << "bam-api " << reads_written
                         << " reads have been written to " << output_filepath
                         << " on the basis of " << input_filepath;

    return reads_written;
}
