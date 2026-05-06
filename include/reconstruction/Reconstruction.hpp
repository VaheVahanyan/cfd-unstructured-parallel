#ifndef RECONSTRUCTION_HPP
#define RECONSTRUCTION_HPP

#include "data/Variables.hpp"

class Mesh;
class Workspace;
class Face;

/**
 * @class Reconstruction
 * @brief Interface for face-based reconstruction on generic cell-centered meshes.
 */
class Reconstruction {
public:
    virtual ~Reconstruction() = default;

    /**
     * @brief Precompute and cache cell-centered gradients.
     *
     * Must be called once per Runge-Kutta stage before the face loop.
     *
     * @param mesh Mesh with geometry and connectivity.
     * @param workspace Workspace for storing computed gradients.
     */
    virtual void ComputeGradients(const Mesh& mesh, Workspace& workspace) const = 0;

    /**
     * @brief Reconstruct owner-side and neighbor-side states on one internal face.
     *
     * @param mesh Mesh with geometry and connectivity.
     * @param workspace Workspace containing primitive cache and cached gradients.
     * @param face Internal face.
     * @param owner_state Reconstructed state on owner side of the face.
     * @param neighbor_state Reconstructed state on neighbor side of the face.
     */
    virtual void ReconstructInteriorFace(const Mesh& mesh,
                                         const Workspace& workspace,
                                         const Face& face,
                                         PrimitiveCell& owner_state,
                                         PrimitiveCell& neighbor_state) const = 0;

    /**
     * @brief Reconstruct owner-side interior state on one boundary face.
     *
     * @param mesh Mesh with geometry and connectivity.
     * @param workspace Workspace containing primitive cache and cached gradients.
     * @param face Boundary face.
     * @param interior_state Reconstructed owner-side state adjacent to the face.
     */
    virtual void ReconstructBoundaryFaceInterior(const Mesh& mesh,
                                                 const Workspace& workspace,
                                                 const Face& face,
                                                 PrimitiveCell& interior_state) const = 0;
};

#endif  // RECONSTRUCTION_HPP