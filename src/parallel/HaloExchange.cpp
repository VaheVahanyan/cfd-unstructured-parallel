#include "parallel/HaloExchange.hpp"

#include <limits>
#include <stdexcept>
#include <utility>

#include "data/DataLayer.hpp"
#include "parallel/MPIContext.hpp"

HaloExchange::HaloExchange(const MPIContext& mpi,
                           std::vector<DomainDecomposition::NeighborHalo> halos)
    : mpi_(&mpi),
      halos_(std::move(halos)) {
    if (!mpi_) {
        throw std::runtime_error("HaloExchange: mpi context is null");
    }

    CreatePacketType();
    InitializeBuffers();
}

HaloExchange::~HaloExchange() {
    DestroyPacketType();
}

void HaloExchange::CreatePacketType() {
    if (packet_type_ != MPI_DATATYPE_NULL) {
        return;
    }

    MPI_Type_contiguous(4, MPI_DOUBLE, &packet_type_);
    MPI_Type_commit(&packet_type_);
}

void HaloExchange::DestroyPacketType() {
    if (packet_type_ != MPI_DATATYPE_NULL) {
        MPI_Type_free(&packet_type_);
        packet_type_ = MPI_DATATYPE_NULL;
    }
}

void HaloExchange::InitializeBuffers() {
    buffers_.clear();
    buffers_.resize(halos_.size());

    for (std::size_t i = 0; i < halos_.size(); ++i) {
        buffers_[i].send.resize(halos_[i].send_local_ids.size());
        buffers_[i].recv.resize(halos_[i].recv_local_ids.size());
    }

    requests_.clear();
    requests_.resize(2 * halos_.size(), MPI_REQUEST_NULL);
}

int HaloExchange::CheckedCount(const std::size_t count) {
    if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("HaloExchange: MPI message count exceeds int range");
    }

    return static_cast<int>(count);
}

void HaloExchange::PackSendBuffer(const DataLayer& layer,
                                  const DomainDecomposition::NeighborHalo& halo,
                                  std::vector<CellStatePacket>& buffer) const {
    if (buffer.size() != halo.send_local_ids.size()) {
        throw std::runtime_error("HaloExchange::PackSendBuffer: buffer size mismatch");
    }

    const auto& U = layer.U();

    for (std::size_t i = 0; i < halo.send_local_ids.size(); ++i) {
        const std::size_t cell_id = halo.send_local_ids[i];

        buffer[i].U[0] = U(cell_id, DataLayer::k_rho);
        buffer[i].U[1] = U(cell_id, DataLayer::k_rhoU);
        buffer[i].U[2] = U(cell_id, DataLayer::k_rhoV);
        buffer[i].U[3] = U(cell_id, DataLayer::k_E);
    }
}

void HaloExchange::UnpackRecvBuffer(DataLayer& layer,
                                    const DomainDecomposition::NeighborHalo& halo,
                                    const std::vector<CellStatePacket>& buffer) const {
    if (buffer.size() != halo.recv_local_ids.size()) {
        throw std::runtime_error("HaloExchange::UnpackRecvBuffer: buffer size mismatch");
    }

    auto& U = layer.U();

    for (std::size_t i = 0; i < halo.recv_local_ids.size(); ++i) {
        const std::size_t cell_id = halo.recv_local_ids[i];

        U(cell_id, DataLayer::k_rho) = buffer[i].U[0];
        U(cell_id, DataLayer::k_rhoU) = buffer[i].U[1];
        U(cell_id, DataLayer::k_rhoV) = buffer[i].U[2];
        U(cell_id, DataLayer::k_E) = buffer[i].U[3];
    }
}

void HaloExchange::Synchronize(DataLayer& layer) const {
    if (!mpi_) {
        throw std::runtime_error("HaloExchange::Synchronize: mpi context is null");
    }

    if (packet_type_ == MPI_DATATYPE_NULL) {
        throw std::runtime_error("HaloExchange::Synchronize: packet MPI datatype is not initialized");
    }

    if (halos_.empty()) {
        return;
    }

    constexpr int k_halo_tag = 1001;

    if (buffers_.size() != halos_.size() || requests_.size() != 2 * halos_.size()) {
        throw std::runtime_error("HaloExchange::Synchronize: internal buffers are inconsistent");
    }

    std::size_t request_id = 0;

    for (std::size_t i = 0; i < halos_.size(); ++i) {
        const DomainDecomposition::NeighborHalo& halo = halos_[i];
        std::vector<CellStatePacket>& recv_buffer = buffers_[i].recv;

        MPI_Irecv(
                  recv_buffer.empty() ? nullptr : recv_buffer.data(),
                  CheckedCount(recv_buffer.size()),
                  packet_type_,
                  halo.remote_rank,
                  k_halo_tag,
                  mpi_->Comm(),
                  &requests_[request_id]
                 );

        ++request_id;
    }

    for (std::size_t i = 0; i < halos_.size(); ++i) {
        const DomainDecomposition::NeighborHalo& halo = halos_[i];
        std::vector<CellStatePacket>& send_buffer = buffers_[i].send;

        PackSendBuffer(layer, halo, send_buffer);

        MPI_Isend(
                  send_buffer.empty() ? nullptr : send_buffer.data(),
                  CheckedCount(send_buffer.size()),
                  packet_type_,
                  halo.remote_rank,
                  k_halo_tag,
                  mpi_->Comm(),
                  &requests_[request_id]
                 );

        ++request_id;
    }

    MPI_Waitall(
                CheckedCount(requests_.size()),
                requests_.data(),
                MPI_STATUSES_IGNORE
               );

    for (std::size_t i = 0; i < halos_.size(); ++i) {
        UnpackRecvBuffer(layer, halos_[i], buffers_[i].recv);
    }
}
