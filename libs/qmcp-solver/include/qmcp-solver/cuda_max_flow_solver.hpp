#ifndef QMCP_CUDA_MAX_FLOW_SOLVER_HPP
#define QMCP_CUDA_MAX_FLOW_SOLVER_HPP

#include <cstdint>
#include <memory>
#include <vector>

#include "bam-api/bam_api.hpp"
#include "bam-api/paired_reads.hpp"
#include "solver.hpp"

namespace qmcp {

class CudaMaxFlowSolver : public Solver {
   public:
    typedef uint32_t Node;
    typedef uint32_t Capacity;
    typedef uint32_t Label;
    typedef uint32_t NeighborInfoIndex;
    typedef int32_t Excess;

    enum class EdgeDirection : uint8_t { Forward, Backward };

    std::unique_ptr<Solution> solve(uint32_t required_cover, bam_api::BamApi& bam_api) override;
    bool uses_quality_of_reads() override { return false; }

    void set_block_size(uint32_t block_size);
    void set_kernel_cycles(uint32_t kernel_cycles);

    static constexpr uint32_t kDefaultBlockSize = 512;
    static constexpr uint32_t kDefaultKernelCycles = 5000;
    static constexpr uint32_t kUseGlobalRelabelMin = 300;

   private:
    void clear_graph();

    // Calling this function is valid only after the import_data is called
    void create_graph(const bam_api::SOAPairedReads& sequence, uint32_t required_cover);

    static void add_edge(std::vector<std::vector<Node>>& neighbors_dict,
                         std::vector<std::vector<EdgeDirection>>& edge_dir_dict,
                         std::vector<std::vector<Capacity>>& residual_capacity_dict,
                         std::vector<std::vector<uint32_t>>& inversed_edge_offset_dict, Node start,
                         Node end, Capacity capacity);

    void global_relabel();

    // This function is responsible for first step of push-relabel algorithm
    Excess create_preflow();

    uint32_t block_size_ = kDefaultBlockSize;
    uint32_t kernel_cycles_ = kDefaultKernelCycles;

    std::filesystem::path input_filepath_;
    bam_api::SOAPairedReads input_sequence_;
    std::vector<uint32_t> max_coverage_;

    // Global-relabel BFS data
    std::vector<bool> is_visited_;

    // === Graph data ===
    std::vector<Excess> excess_func_;
    std::vector<Label> label_func_;
    std::vector<bool> is_marked_;

    // Maps node to start/end index in neighbors info arrays
    std::vector<NeighborInfoIndex> neighbors_start_ind_;
    std::vector<NeighborInfoIndex> neighbors_end_ind_;

    // Maps read index to neighbor index (not neighbors info index)
    // For example for node 1 where N(1) = [5, 3, 6]:
    // - read (1,2) -> 1
    // - read (1,4) -> 0
    // - read (1,5) -> 2
    // It's required for creating output from residual network
    std::vector<uint32_t> read_ind_to_neighbor_offset_;

    // Neighbors info array is an array that stores packed information about the
    // neighbors of all vertices:

    // Maps NeighborInfoIndex into neighbor Node
    std::vector<Node> neighbors_;

    // Maps NeighborInfoIndex into neighbor's neighbor NeighborInfoIndex
    std::vector<NeighborInfoIndex> inversed_edge_offset_;

    // Maps NeighborInfoIndex into residual capacity
    std::vector<Capacity> residual_capacity_;

    // Maps NeighborInfoIndex into EdgeDirection
    std::vector<EdgeDirection> edge_dir_;
};

}  // namespace qmcp

#endif
