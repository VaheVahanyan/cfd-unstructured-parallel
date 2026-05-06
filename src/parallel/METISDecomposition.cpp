#include "parallel/METISDecomposition.hpp"

#include <stdexcept>
#include <vector>

#include <metis.h>

#include "geometry/Face.hpp"
#include "geometry/Mesh.hpp"

std::vector<int> METISDecomposition::ComputePartition(const Mesh& global_mesh, const int nproc) const {
    if (nproc <= 0) {
        throw std::runtime_error("METISDecomposition: MPI size must be positive");
    }

    const std::size_t n_cells = global_mesh.GetCellCount();
    std::vector part(n_cells, 0);

    if (nproc == 1) {
        return part;
    }

    std::vector<std::vector<idx_t>> adjacency_list(n_cells);
    for (const Face& face : global_mesh.Faces()) {
        if (face.IsInternal()) {
            adjacency_list[face.owner_cell_id].push_back(static_cast<idx_t>(face.neighbor_cell_id));
            adjacency_list[face.neighbor_cell_id].push_back(static_cast<idx_t>(face.owner_cell_id));
        }
    }

    std::vector<idx_t> xadj;
    std::vector<idx_t> adjncy;
    xadj.reserve(n_cells + 1);

    idx_t edge_count = 0;
    for (std::size_t i = 0; i < n_cells; ++i) {
        xadj.push_back(edge_count);
        for (const idx_t neighbor : adjacency_list[i]) {
            adjncy.push_back(neighbor);
            edge_count++;
        }
    }
    xadj.push_back(edge_count);

    idx_t nvtxs = static_cast<idx_t>(n_cells);
    idx_t ncon = 1;
    idx_t nparts = nproc;
    idx_t objval = 0;

    std::vector<idx_t> metis_part(n_cells, 0);

    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 0;

    const int status = METIS_PartGraphKway(
                                           &nvtxs, &ncon, xadj.data(), adjncy.data(),
                                           nullptr, nullptr, nullptr, &nparts, nullptr, nullptr, options, &objval,
                                           metis_part.data()
                                          );

    if (status != METIS_OK) {
        throw std::runtime_error("METISDecomposition: METIS_PartGraphKway failed");
    }

    for (std::size_t i = 0; i < n_cells; ++i) {
        part[i] = metis_part[i];
    }

    return part;
}
