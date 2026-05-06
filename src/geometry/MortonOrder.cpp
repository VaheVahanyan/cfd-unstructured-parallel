#include "geometry/MortonOrder.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "geometry/Cell.hpp"
#include "geometry/Face.hpp"
#include "geometry/Mesh.hpp"

unsigned long long MortonOrder::Part1By1(unsigned int value) {
    unsigned long long x = value;
    x &= 0x00000000ffffffffULL;
    x = (x | x << 16U) & 0x0000ffff0000ffffULL;
    x = (x | x << 8U) & 0x00ff00ff00ff00ffULL;
    x = (x | x << 4U) & 0x0f0f0f0f0f0f0f0fULL;
    x = (x | x << 2U) & 0x3333333333333333ULL;
    x = (x | x << 1U) & 0x5555555555555555ULL;
    return x;
}

unsigned long long MortonOrder::Encode2D(const unsigned int x, const unsigned int y) {
    return Part1By1(x) | (Part1By1(y) << 1U);
}

std::vector<MortonOrder::Entry> MortonOrder::BuildEntries(const Mesh& mesh,
                                                          const std::size_t begin,
                                                          const std::size_t end) {
    if (begin > end || end > mesh.GetCellCount()) {
        throw std::runtime_error("MortonOrder::BuildEntries: invalid range");
    }

    std::vector<Entry> entries;
    entries.reserve(end - begin);

    if (begin == end) {
        return entries;
    }

    double min_x = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();

    for (std::size_t i = begin; i < end; ++i) {
        const Cell& cell = mesh.GetCell(i);
        min_x = std::min(min_x, cell.center_x);
        max_x = std::max(max_x, cell.center_x);
        min_y = std::min(min_y, cell.center_y);
        max_y = std::max(max_y, cell.center_y);
    }

    const double dx = max_x - min_x;
    const double dy = max_y - min_y;
    constexpr double max_coord = 4294967295.0;

    for (std::size_t i = begin; i < end; ++i) {
        const Cell& cell = mesh.GetCell(i);

        double nx = 0.0;
        double ny = 0.0;

        if (dx > 0.0) {
            nx = (cell.center_x - min_x) / dx;
        }

        if (dy > 0.0) {
            ny = (cell.center_y - min_y) / dy;
        }

        nx = std::clamp(nx, 0.0, 1.0);
        ny = std::clamp(ny, 0.0, 1.0);

        const auto ix = static_cast<unsigned int>(nx * max_coord);
        const auto iy = static_cast<unsigned int>(ny * max_coord);

        Entry entry;
        entry.old_id = i;
        entry.code = Encode2D(ix, iy);
        entries.push_back(entry);
    }

    std::stable_sort(entries.begin(), entries.end(),
                     [](const Entry& lhs, const Entry& rhs) {
                         if (lhs.code != rhs.code) {
                             return lhs.code < rhs.code;
                         }
                         return lhs.old_id < rhs.old_id;
                     });

    for (std::size_t i = 0; i < entries.size(); ++i) {
        entries[i].new_id = begin + i;
    }

    return entries;
}

void MortonOrder::ReorderRange(Mesh& mesh,
                               const std::size_t begin,
                               const std::size_t end,
                               std::vector<std::size_t>& old_to_new) {
    if (begin == end) {
        return;
    }

    const std::vector<Entry> entries = BuildEntries(mesh, begin, end);
    std::vector<Cell> reordered(end - begin);

    for (const Entry& entry : entries) {
        old_to_new[entry.old_id] = entry.new_id;
        reordered[entry.new_id - begin] = mesh.Cells()[entry.old_id];
    }

    for (std::size_t i = 0; i < reordered.size(); ++i) {
        Cell& cell = reordered[i];
        cell.local_id = begin + i;
        mesh.Cells()[begin + i] = std::move(cell);
    }
}

void MortonOrder::ApplyToCells(Mesh& mesh) {
    const std::size_t n_cells = mesh.GetCellCount();

    if (n_cells == 0) {
        return;
    }

    std::size_t n_owned = mesh.GetOwnedCellCount();

    if (n_owned == 0 || n_owned > n_cells) {
        n_owned = n_cells;
    }

    std::vector<std::size_t> old_to_new(n_cells);
    for (std::size_t i = 0; i < n_cells; ++i) {
        old_to_new[i] = i;
    }

    ReorderRange(mesh, 0, n_owned, old_to_new);
    ReorderRange(mesh, n_owned, n_cells, old_to_new);

    for (Face& face : mesh.Faces()) {
        if (face.owner_cell_id != Face::k_invalid_cell_id) {
            face.owner_cell_id = old_to_new[face.owner_cell_id];
        }

        if (face.neighbor_cell_id != Face::k_invalid_cell_id) {
            face.neighbor_cell_id = old_to_new[face.neighbor_cell_id];
        }
    }

    for (Cell& cell : mesh.Cells()) {
        cell.face_ids.clear();
    }

    for (const Face& face : mesh.Faces()) {
        mesh.Cells()[face.owner_cell_id].face_ids.push_back(face.id);

        if (face.IsInternal() || face.IsMPIBoundary()) {
            mesh.Cells()[face.neighbor_cell_id].face_ids.push_back(face.id);
        }
    }

    mesh.Validate();
}
