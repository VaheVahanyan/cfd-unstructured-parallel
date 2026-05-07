#ifndef RUNSTATISTICS_HPP
#define RUNSTATISTICS_HPP

#include <cstddef>

class HaloExchange;
class Mesh;
class MPIContext;
struct Settings;

/**
 * @class RunStatistics
 * @brief Collects and writes runtime statistics for scaling studies.
 *
 * The class writes one CSV file per MPI rank in order to avoid concurrent
 * writes to the same file. Aggregation and plotting are expected to be done
 * later by Python scripts.
 */
class RunStatistics final {
public:
    /**
     * @brief Print aggregated mesh and halo statistics to stdout on root rank.
     */
    static void PrintParallelSummary(const Settings& settings,
                                     const Mesh& mesh,
                                     const MPIContext& mpi,
                                     const HaloExchange* halo_exchange);

    /**
     * @brief Write per-rank timing and mesh statistics to a CSV file.
     */
    static void WriteRankTimingCsv(const Settings& settings,
                                   const Mesh& mesh,
                                   const MPIContext& mpi,
                                   const HaloExchange* halo_exchange,
                                   std::size_t steps,
                                   double computation_time,
                                   double wall_time);

private:
    struct LocalMeshStats final {
        std::size_t owned_cells = 0;
        std::size_t ghost_cells = 0;
        std::size_t total_cells = 0;
        std::size_t nodes = 0;
        std::size_t faces = 0;
        std::size_t internal_faces = 0;
        std::size_t mpi_faces = 0;
        std::size_t physical_faces = 0;
        std::size_t halo_neighbors = 0;
        std::size_t halo_send_cells = 0;
        std::size_t halo_recv_cells = 0;
    };

    struct ReducedStats final {
        double sum = 0.0;
        double min = 0.0;
        double max = 0.0;
    };

    [[nodiscard]] static LocalMeshStats CollectLocalMeshStats(const Mesh& mesh,
                                                              const HaloExchange* halo_exchange);

    [[nodiscard]] static ReducedStats ReduceDouble(double value,
                                                   const MPIContext& mpi);

    static void EnsureStatisticsDirectory(const Settings& settings,
                                          const MPIContext& mpi);

    [[nodiscard]] static std::string StatisticsDirectory(const Settings& settings);

    [[nodiscard]] static std::string RankCsvFilename(const Settings& settings,
                                                     const MPIContext& mpi);
};

#endif

