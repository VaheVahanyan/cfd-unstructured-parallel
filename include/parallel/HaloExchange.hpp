#ifndef HALOEXCHANGE_HPP
#define HALOEXCHANGE_HPP

#include <cstddef>
#include <vector>

#include <mpi.h>

#include "parallel/DomainDecomposition.hpp"
#include "parallel/StateSynchronizer.hpp"

class DataLayer;
class MPIContext;

/**
 * @class HaloExchange
 * @brief Reusable non-blocking MPI halo exchange for conservative ghost-cell states.
 *
 * Data packet sent per cell:
 * - rho
 * - rhoU
 * - rhoV
 * - E
 *
 * Communication pattern:
 * - reuse preallocated send, receive, and request buffers
 * - post all receives
 * - pack owned send cells
 * - post all sends
 * - wait for all requests
 * - unpack received states into local ghost cells
 *
 * The class assumes MPI calls are made by the main thread only. This matches
 * MPI_THREAD_FUNNELED and is compatible with OpenMP regions outside MPI calls.
 */
class HaloExchange final : public StateSynchronizer {
public:
    struct CellStatePacket final {
        double U[4];
    };

    HaloExchange(const MPIContext& mpi,
                 std::vector<DomainDecomposition::NeighborHalo> halos);

    ~HaloExchange();

    HaloExchange(const HaloExchange&) = delete;
    HaloExchange& operator=(const HaloExchange&) = delete;

    HaloExchange(HaloExchange&&) = delete;
    HaloExchange& operator=(HaloExchange&&) = delete;

    /**
     * @brief Exchange owned cell states with neighboring ranks and update ghost cells.
     */
    void Synchronize(DataLayer& layer) const override;

private:
    struct ExchangeBuffer final {
        std::vector<CellStatePacket> send;
        std::vector<CellStatePacket> recv;
    };

    const MPIContext* mpi_ = nullptr;
    std::vector<DomainDecomposition::NeighborHalo> halos_;
    mutable std::vector<ExchangeBuffer> buffers_;
    mutable std::vector<MPI_Request> requests_;

    MPI_Datatype packet_type_ = MPI_DATATYPE_NULL;

    void CreatePacketType();
    void DestroyPacketType();
    void InitializeBuffers();

    [[nodiscard]] static int CheckedCount(std::size_t count);

    void PackSendBuffer(const DataLayer& layer,
                        const DomainDecomposition::NeighborHalo& halo,
                        std::vector<CellStatePacket>& buffer) const;

    void UnpackRecvBuffer(DataLayer& layer,
                          const DomainDecomposition::NeighborHalo& halo,
                          const std::vector<CellStatePacket>& buffer) const;
};

#endif
