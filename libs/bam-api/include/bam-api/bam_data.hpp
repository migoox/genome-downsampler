#ifndef SEQUENCE_T_HPP
#define SEQUENCE_T_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace bam_api {

struct Read {
  int64_t id;
  int64_t start_ind;
  int64_t end_ind;
  uint32_t quality;
  std::string qname;
  bool is_first_read;
};

struct AOSPairedReads {
  std::vector<Read> reads;
  std::vector<int64_t> read_pair_map;
  int64_t ref_genome_length;
};

struct SOAPairedReads {
  std::vector<int64_t> ids;
  std::vector<int64_t> start_inds;
  std::vector<int64_t> end_inds;
  std::vector<uint32_t> qualities;
  std::vector<std::string> qnames;
  std::vector<bool> is_first_read;
  std::vector<int64_t> read_pair_map;
  int64_t ref_genome_length;
};

} // namespace bam_api

#endif
