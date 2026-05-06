#include "parallel/RCBDecomposition.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "geometry/Cell.hpp"
#include "geometry/Mesh.hpp"

namespace
{
    [[nodiscard]] std::vector<std::size_t> MakeAllCellIds(const Mesh& mesh) {
        std::vector<std::size_t> ids(mesh.GetCellCount());
        for (std::size_t i = 0; i < ids.size(); ++i) {
            ids[i] = i;
        }
        return ids;
    }

    [[nodiscard]] double CellCoordByAxis(const Cell& cell, const int axis) {
        if (axis == 0) {
            return cell.center_x;
        }
        return cell.center_y;
    }

    [[nodiscard]] int ChooseSplitAxis(const Mesh& mesh, const std::vector<std::size_t>& cell_ids) {
        double min_x = 0.0, max_x = 0.0;
        double min_y = 0.0, max_y = 0.0;

        bool first = true;
        for (const std::size_t cell_id : cell_ids) {
            const Cell& cell = mesh.GetCell(cell_id);

            if (first) {
                min_x = max_x = cell.center_x;
                min_y = max_y = cell.center_y;
                first = false;
                continue;
            }

            min_x = std::min(min_x, cell.center_x);
            max_x = std::max(max_x, cell.center_x);
            min_y = std::min(min_y, cell.center_y);
            max_y = std::max(max_y, cell.center_y);
        }

        const double span_x = max_x - min_x;
        const double span_y = max_y - min_y;

        if (span_x >= span_y) {
            return 0;
        }
        return 1;
    }

    void BuildRCBRecursive(const Mesh& mesh,
                           const std::vector<std::size_t>& cell_ids,
                           const int proc_begin,
                           const int proc_end,
                           std::vector<int>& part) {
        const int proc_count = proc_end - proc_begin;
        if (proc_count <= 0) {
            throw std::runtime_error("RCBDecomposition: invalid processor interval");
        }

        if (proc_count == 1) {
            for (const std::size_t cell_id : cell_ids) {
                part[cell_id] = proc_begin;
            }
            return;
        }

        if (cell_ids.empty()) {
            return;
        }

        const int axis = ChooseSplitAxis(mesh, cell_ids);

        std::vector<std::size_t> sorted = cell_ids;
        std::stable_sort(sorted.begin(), sorted.end(),
                         [&](const std::size_t lhs, const std::size_t rhs) {
                             const Cell& a = mesh.GetCell(lhs);
                             const Cell& b = mesh.GetCell(rhs);
                             const double ca = CellCoordByAxis(a, axis);
                             const double cb = CellCoordByAxis(b, axis);
                             if (ca != cb) {
                                 return ca < cb;
                             }
                             return lhs < rhs;
                         });

        const int left_proc_count = proc_count / 2;
        const int right_proc_count = proc_count - left_proc_count;

        const std::size_t n = sorted.size();
        const std::size_t cut = static_cast<std::size_t>(std::llround(
                                                                      static_cast<double>(n) * static_cast<double>(
                                                                          left_proc_count) / static_cast<double>(
                                                                          proc_count)
                                                                     ));

        const std::vector<std::size_t> left(
                                            sorted.begin(),
                                            sorted.begin() + static_cast<std::ptrdiff_t>(cut)
                                           );
        const std::vector<std::size_t> right(
                                             sorted.begin() + static_cast<std::ptrdiff_t>(cut),
                                             sorted.end()
                                            );

        BuildRCBRecursive(mesh, left, proc_begin, proc_begin + left_proc_count, part);
        BuildRCBRecursive(mesh, right, proc_begin + left_proc_count, proc_end, part);
    }
} // namespace

std::vector<int> RCBDecomposition::ComputePartition(const Mesh& global_mesh, const int nproc) const {
    if (nproc <= 0) {
        throw std::runtime_error("RCBDecomposition: MPI size must be positive");
    }

    std::vector<int> part(global_mesh.GetCellCount(), -1);
    const std::vector<std::size_t> all_ids = MakeAllCellIds(global_mesh);

    BuildRCBRecursive(global_mesh, all_ids, 0, nproc, part);

    for (std::size_t i = 0; i < part.size(); ++i) {
        if (part[i] < 0 || part[i] >= nproc) {
            throw std::runtime_error("RCBDecomposition: invalid partition result");
        }
    }

    return part;
}
