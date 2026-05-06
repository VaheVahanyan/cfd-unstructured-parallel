#ifndef RCBDECOMPOSITION_HPP
#define RCBDECOMPOSITION_HPP

#include <vector>

#include "parallel/DomainDecomposition.hpp"

class Mesh;

/**
 * @class RCBDecomposition
 * @brief Recursive Coordinate Bisection (RCB) domain decomposition.
 *
 * Divides the computational domain by recursively splitting cells along
 * the longest spatial axis.
 */
class RCBDecomposition final : public DomainDecomposition {
protected:
    [[nodiscard]] std::vector<int> ComputePartition(const Mesh& global_mesh,
                                                    int nproc) const override;
};

#endif  // RCBDECOMPOSITION_HPP
