#ifndef MESHDISTRIBUTION_HPP
#define MESHDISTRIBUTION_HPP

#include <vector>

#include "geometry/Mesh.hpp"
#include "parallel/DomainDecomposition.hpp"

class MPIContext;

/**
 * @class MeshDistribution
 * @brief Distributes a global mesh from the root MPI rank to all other ranks.
 */
class MeshDistribution final {
public:
    struct LocalPartition final {
        Mesh mesh;
        std::vector<DomainDecomposition::NeighborHalo> halos;
    };

    /**
     * @brief Performs central partitioning on the root rank and distributes data via MPI.
     *
     * @param global_mesh The global mesh (must be valid on root, null elsewhere).
     * @param mpi The MPI context.
     * @param decomposer The decomposition strategy (e.g., RCB or METIS).
     * @return Local partition containing the local mesh and halo metadata.
     */
    [[nodiscard]] static LocalPartition DistributeFromRoot(
        const Mesh* global_mesh,
        const MPIContext& mpi,
        const DomainDecomposition& decomposer);
};

#endif  // MESHDISTRIBUTION_HPP
