#ifndef P1RECONSTRUCTION_HPP
#define P1RECONSTRUCTION_HPP

#include <cstddef>
#include <vector>

#include "reconstruction/Reconstruction.hpp"

class Mesh;
class Workspace;
class Face;
struct Cell;

/**
 * @class P1Reconstruction
 * @brief Second-order 2D spatial reconstruction with gradient caching.
 *
 * Computes unlimited cell-centered gradients via least-squares, applies
 * a minmod/Barth-Jespersen limiter to enforce monotonicity, and caches
 * the results in the Workspace. Face reconstructions use the cached
 * gradients to perform linear extrapolation.
 */
class P1Reconstruction final : public Reconstruction {
public:
    P1Reconstruction() = default;
    ~P1Reconstruction() override = default;

    void ComputeGradients(const Mesh& mesh, Workspace& workspace) const override;

    void ReconstructInteriorFace(const Mesh& mesh,
                                 const Workspace& workspace,
                                 const Face& face,
                                 PrimitiveCell& owner_state,
                                 PrimitiveCell& neighbor_state) const override;

    void ReconstructBoundaryFaceInterior(const Mesh& mesh,
                                         const Workspace& workspace,
                                         const Face& face,
                                         PrimitiveCell& interior_state) const override;

private:
    struct PrimitiveGradient final {
        PrimitiveCell dx;
        PrimitiveCell dy;
    };

    [[nodiscard]] PrimitiveCell LoadCellPrimitive(const Workspace& workspace, std::size_t cell_id) const;

    void CollectNeighborCellIds(const Mesh& mesh, const Cell& cell, std::vector<std::size_t>& neighbor_ids) const;

    [[nodiscard]] PrimitiveGradient ComputeUnlimitedGradient(const Mesh& mesh, const Workspace& workspace,
                                                             const Cell& cell,
                                                             const std::vector<std::size_t>& neighbor_ids) const;

    [[nodiscard]] double ComputeLimiter(const Mesh& mesh, const Workspace& workspace, const Cell& cell,
                                        const std::vector<std::size_t>& neighbor_ids,
                                        const PrimitiveGradient& grad) const;

    [[nodiscard]] double ComputeBarthJespersenPhi(double w_cell, double w_min, double w_max,
                                                  double w_face_candidate) const;
};

#endif  // P1RECONSTRUCTION_HPP
