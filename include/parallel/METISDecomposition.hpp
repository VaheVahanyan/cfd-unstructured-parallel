#ifndef METISDECOMPOSITION_HPP
#define METISDECOMPOSITION_HPP

#include <vector>

#include "parallel/DomainDecomposition.hpp"

class Mesh;

/**
 * @class METISDecomposition
 * @brief Graph-based domain decomposition using the METIS library.
 *
 * Constructs a dual graph from the mesh connectivity (cells as vertices,
 * internal faces as edges) and utilizes METIS to compute a k-way partition
 * that minimizes edge-cut while balancing the computational load.
 */
class METISDecomposition final : public DomainDecomposition {
protected:
    [[nodiscard]] std::vector<int> ComputePartition(const Mesh& global_mesh,
                                                    int nproc) const override;
};

#endif  // METISDECOMPOSITION_HPP
