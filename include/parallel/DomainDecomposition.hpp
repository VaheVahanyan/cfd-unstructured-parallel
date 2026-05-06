#ifndef DOMAINDECOMPOSITION_HPP
#define DOMAINDECOMPOSITION_HPP

#include <cstddef>
#include <unordered_map>
#include <vector>

#include "geometry/Mesh.hpp"

class MPIContext;

/**
 * @class DomainDecomposition
 * @brief Abstract base class for domain decomposition strategies.
 *
 * Defines the interface and common routines for extracting an MPI-local
 * mesh, ghost cells, and halo exchange metadata from a global mesh and
 * a computed partition array.
 */
class DomainDecomposition {
public:
    struct NeighborHalo final {
        int remote_rank = -1;
        std::vector<std::size_t> send_local_ids;
        std::vector<std::size_t> recv_local_ids;
    };

    struct Result final {
        Mesh local_mesh;

        std::size_t n_owned_cells = 0;
        std::size_t n_ghost_cells = 0;

        std::vector<int> global_part;
        std::vector<std::size_t> local_to_global_cell;
        std::unordered_map<std::size_t, std::size_t> global_to_local_cell;
        std::vector<std::size_t> owned_global_ids;
        std::vector<std::size_t> ghost_global_ids;

        std::vector<NeighborHalo> halos;
    };

    virtual ~DomainDecomposition() = default;

    /**
     * @brief Decomposes the global mesh for the current MPI rank.
     *
     * @param global_mesh The global mesh shared by all ranks.
     * @param mpi The MPI context providing rank and size.
     * @return The local mesh and communication metadata.
     */
    [[nodiscard]] Result Decompose(const Mesh& global_mesh, const MPIContext& mpi) const;

    /**
     * @brief Computes the partition array for the global mesh.
     *
     * @param global_mesh The global mesh to partition.
     * @param nproc The number of MPI processes.
     * @return A vector mapping global cell index to MPI rank.
     */
    [[nodiscard]] virtual std::vector<int> ComputePartition(const Mesh& global_mesh,
                                                            int nproc) const = 0;

    [[nodiscard]] Result BuildLocalResultForRank(const Mesh& global_mesh,
                                                 const std::vector<int>& global_part,
                                                 int rank) const;
};

#endif  // DOMAINDECOMPOSITION_HPP
