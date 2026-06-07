#ifndef QMCP_ORTOOLS_FLOW_TYPES_HPP
#define QMCP_ORTOOLS_FLOW_TYPES_HPP

#include <cassert>
#include <cstdint>
#include <limits>

#include "bam-api/read.hpp"
#include "ortools/graph/ebert_graph.h"

namespace qmcp {

using OrNode = operations_research::NodeIndex;
using OrArc = operations_research::ArcIndex;
using OrFlow = operations_research::FlowQuantity;
using OrCost = operations_research::CostValue;

inline OrNode to_node(bam_api::Index index) {
    assert(index <= static_cast<bam_api::Index>(std::numeric_limits<OrNode>::max()));
    return static_cast<OrNode>(index);
}

inline OrArc to_arc(bam_api::ReadIndex read_index) {
    assert(read_index <= static_cast<bam_api::ReadIndex>(std::numeric_limits<OrArc>::max()));
    return static_cast<OrArc>(read_index);
}

inline OrFlow to_flow(int64_t quantity) {
    assert(quantity >= 0);
    assert(quantity <= std::numeric_limits<OrFlow>::max());
    return static_cast<OrFlow>(quantity);
}

inline OrFlow to_signed_flow(int64_t quantity) {
    assert(quantity >= std::numeric_limits<OrFlow>::min());
    assert(quantity <= std::numeric_limits<OrFlow>::max());
    return static_cast<OrFlow>(quantity);
}

inline OrCost to_cost(int64_t cost) {
    assert(cost >= std::numeric_limits<OrCost>::min());
    assert(cost <= std::numeric_limits<OrCost>::max());
    return static_cast<OrCost>(cost);
}

}  // namespace qmcp

#endif
