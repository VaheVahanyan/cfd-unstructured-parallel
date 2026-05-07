#include "output/RunStatistics.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <mpi.h>

#include "config/Settings.hpp"
#include "geometry/Face.hpp"
#include "geometry/Mesh.hpp"
#include "parallel/HaloExchange.hpp"
#include "parallel/MPIContext.hpp"

RunStatistics::LocalMeshStats RunStatistics::CollectLocalMeshStats(
    const Mesh& mesh,
    const HaloExchange* halo_exchange
) {
    LocalMeshStats stats;

    stats.owned_cells = mesh.GetOwnedCellCount();
    stats.ghost_cells = mesh.GetGhostCellCount();
    stats.total_cells = mesh.GetCellCount();
    stats.nodes = mesh.GetNodeCount();
    stats.faces = mesh.GetFaceCount();

    for (const Face& face : mesh.Faces()) {
        if (face.IsInternal()) {
            ++stats.internal_faces;
        } else if (face.IsMPIBoundary()) {
            ++stats.mpi_faces;
        } else if (face.IsPhysicalBoundary()) {
            ++stats.physical_faces;
        }
    }

    if (halo_exchange) {
        stats.halo_neighbors = halo_exchange->GetNeighborCount();
        stats.halo_send_cells = halo_exchange->GetTotalSendCellCount();
        stats.halo_recv_cells = halo_exchange->GetTotalRecvCellCount();
    }

    return stats;
}

RunStatistics::ReducedStats RunStatistics::ReduceDouble(const double value,
                                                        const MPIContext& mpi) {
    ReducedStats stats;

    MPI_Allreduce(&value, &stats.sum, 1, MPI_DOUBLE, MPI_SUM, mpi.Comm());
    MPI_Allreduce(&value, &stats.min, 1, MPI_DOUBLE, MPI_MIN, mpi.Comm());
    MPI_Allreduce(&value, &stats.max, 1, MPI_DOUBLE, MPI_MAX, mpi.Comm());

    return stats;
}

std::string RunStatistics::StatisticsDirectory(const Settings& settings) {
    return settings.output_dir + "/statistics";
}

void RunStatistics::EnsureStatisticsDirectory(const Settings& settings,
                                              const MPIContext& mpi) {
    const std::string directory = StatisticsDirectory(settings);

    if (mpi.IsRoot()) {
        std::filesystem::create_directories(directory);
    }

    mpi.Barrier();
}

std::string RunStatistics::RankCsvFilename(const Settings& settings,
                                           const MPIContext& mpi) {
    std::ostringstream oss;

    oss << StatisticsDirectory(settings)
        << "/"
        << settings.simulation_case
        << "__np_" << std::setw(3) << std::setfill('0') << mpi.Size()
        << "__rank_" << std::setw(3) << std::setfill('0') << mpi.Rank()
        << ".csv";

    return oss.str();
}

void RunStatistics::PrintParallelSummary(const Settings& settings,
                                         const Mesh& mesh,
                                         const MPIContext& mpi,
                                         const HaloExchange* halo_exchange) {
    const LocalMeshStats local = CollectLocalMeshStats(mesh, halo_exchange);

    const ReducedStats owned = ReduceDouble(static_cast<double>(local.owned_cells), mpi);
    const ReducedStats ghosts = ReduceDouble(static_cast<double>(local.ghost_cells), mpi);
    const ReducedStats mpi_faces = ReduceDouble(static_cast<double>(local.mpi_faces), mpi);
    const ReducedStats halo_send = ReduceDouble(static_cast<double>(local.halo_send_cells), mpi);
    const ReducedStats halo_recv = ReduceDouble(static_cast<double>(local.halo_recv_cells), mpi);

    if (!mpi.IsRoot()) {
        return;
    }

    const double inv_size = 1.0 / static_cast<double>(mpi.Size());

    std::cout << "Parallel mesh statistics:\n";
    std::cout << "  case: " << settings.simulation_case << '\n';
    std::cout << "  mpi size: " << mpi.Size() << '\n';
    std::cout << "  omp threads: " << settings.omp_threads << '\n';
    std::cout << "  decomposition: " << settings.domain_decomposition_method << '\n';
    std::cout << "  use_morton: " << (settings.use_morton ? "true" : "false") << '\n';
    std::cout << "  owned cells total: " << static_cast<std::size_t>(owned.sum) << '\n';
    std::cout << "  owned cells min/avg/max: "
        << owned.min << " / " << owned.sum * inv_size << " / " << owned.max << '\n';
    std::cout << "  ghost cells total: " << static_cast<std::size_t>(ghosts.sum) << '\n';
    std::cout << "  ghost cells min/avg/max: "
        << ghosts.min << " / " << ghosts.sum * inv_size << " / " << ghosts.max << '\n';
    std::cout << "  MPI boundary faces total: " << static_cast<std::size_t>(mpi_faces.sum) << '\n';
    std::cout << "  MPI boundary faces min/avg/max: "
        << mpi_faces.min << " / " << mpi_faces.sum * inv_size << " / " << mpi_faces.max << '\n';
    std::cout << "  halo send cells total: " << static_cast<std::size_t>(halo_send.sum) << '\n';
    std::cout << "  halo recv cells total: " << static_cast<std::size_t>(halo_recv.sum) << '\n';
    std::cout << '\n';
}

void RunStatistics::WriteRankTimingCsv(const Settings& settings,
                                       const Mesh& mesh,
                                       const MPIContext& mpi,
                                       const HaloExchange* halo_exchange,
                                       const std::size_t steps,
                                       const double computation_time,
                                       const double wall_time) {
    EnsureStatisticsDirectory(settings, mpi);

    const LocalMeshStats stats = CollectLocalMeshStats(mesh, halo_exchange);
    const std::string filename = RankCsvFilename(settings, mpi);

    const bool file_exists = std::filesystem::exists(filename);

    std::ofstream out(filename, std::ios::app);

    if (!out) {
        throw std::runtime_error("RunStatistics: failed to open statistics CSV file");
    }

    if (!file_exists) {
        out << "case_name,rank,size,omp_threads,decomposition,use_morton,"
            << "total_cells,owned_cells,ghost_cells,nodes,faces,"
            << "internal_faces,mpi_faces,physical_faces,"
            << "halo_neighbors,halo_send_cells,halo_recv_cells,"
            << "steps,computation_time,wall_time\n";
    }

    out << settings.simulation_case << ','
        << mpi.Rank() << ','
        << mpi.Size() << ','
        << settings.omp_threads << ','
        << settings.domain_decomposition_method << ','
        << (settings.use_morton ? 1 : 0) << ','
        << stats.total_cells << ','
        << stats.owned_cells << ','
        << stats.ghost_cells << ','
        << stats.nodes << ','
        << stats.faces << ','
        << stats.internal_faces << ','
        << stats.mpi_faces << ','
        << stats.physical_faces << ','
        << stats.halo_neighbors << ','
        << stats.halo_send_cells << ','
        << stats.halo_recv_cells << ','
        << steps << ','
        << std::setprecision(17) << computation_time << ','
        << std::setprecision(17) << wall_time << '\n';

    mpi.Barrier();
}
