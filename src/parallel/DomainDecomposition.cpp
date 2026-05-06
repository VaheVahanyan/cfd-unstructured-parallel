#include "parallel/DomainDecomposition.hpp"

#include <algorithm>
#include <map>
#include <set>
#include <unordered_set>
#include <utility>

#include "geometry/Cell.hpp"
#include "geometry/Face.hpp"
#include "geometry/Node.hpp"
#include "parallel/MPIContext.hpp"

namespace
{
    struct HaloAccumulator final {
        int remote_rank = -1;
        std::vector<std::size_t> send_local_ids;
        std::vector<std::size_t> recv_local_ids;
        std::unordered_set<std::size_t> send_seen;
        std::unordered_set<std::size_t> recv_seen;
    };

    [[nodiscard]] std::vector<std::size_t> CollectOwnedGlobalIds(const Mesh& global_mesh,
                                                                 const std::vector<int>& part,
                                                                 const int rank) {
        std::vector<std::size_t> owned;
        owned.reserve(global_mesh.GetCellCount());

        for (std::size_t global_cell_id = 0; global_cell_id < global_mesh.GetCellCount(); ++global_cell_id) {
            if (part[global_cell_id] == rank) {
                owned.push_back(global_cell_id);
            }
        }
        return owned;
    }

    [[nodiscard]] std::set<std::size_t> CollectGhostGlobalIds(const Mesh& global_mesh,
                                                              const std::vector<int>& part,
                                                              const int rank) {
        std::set<std::size_t> ghosts;

        for (const Face& face : global_mesh.Faces()) {
            if (!face.IsInternal()) {
                continue;
            }

            const int owner_rank = part[face.owner_cell_id];
            const int neighbor_rank = part[face.neighbor_cell_id];

            if (owner_rank == rank && neighbor_rank != rank) {
                ghosts.insert(face.neighbor_cell_id);
            } else if (neighbor_rank == rank && owner_rank != rank) {
                ghosts.insert(face.owner_cell_id);
            }
        }
        return ghosts;
    }

    void CopyLocalCells(const Mesh& global_mesh,
                        const std::vector<std::size_t>& owned_global_ids,
                        const std::vector<std::size_t>& ghost_global_ids,
                        DomainDecomposition::Result& result) {
        auto& local_cells = result.local_mesh.Cells();

        result.local_to_global_cell.clear();
        result.global_to_local_cell.clear();

        local_cells.reserve(owned_global_ids.size() + ghost_global_ids.size());
        result.local_to_global_cell.reserve(owned_global_ids.size() + ghost_global_ids.size());

        auto copy_one = [&](const std::size_t global_cell_id) {
            Cell cell = global_mesh.GetCell(global_cell_id);
            cell.local_id = local_cells.size();
            cell.face_ids.clear();

            result.global_to_local_cell[global_cell_id] = cell.local_id;
            result.local_to_global_cell.push_back(global_cell_id);

            local_cells.push_back(std::move(cell));
        };

        for (const std::size_t global_cell_id : owned_global_ids) {
            copy_one(global_cell_id);
        }
        for (const std::size_t global_cell_id : ghost_global_ids) {
            copy_one(global_cell_id);
        }

        result.n_owned_cells = owned_global_ids.size();
        result.n_ghost_cells = ghost_global_ids.size();

        result.local_mesh.SetOwnedCellCount(result.n_owned_cells);
        result.local_mesh.SetGhostCellCount(result.n_ghost_cells);
    }

    void CopyLocalNodesAndRemapConnectivity(DomainDecomposition::Result& result,
                                            const Mesh& global_mesh) {
        std::set<std::size_t> used_global_node_ids;

        for (const Cell& local_cell : result.local_mesh.Cells()) {
            const std::size_t global_cell_id = result.local_to_global_cell[local_cell.local_id];
            const Cell& global_cell = global_mesh.GetCell(global_cell_id);

            for (const std::size_t global_node_id : global_cell.node_ids) {
                used_global_node_ids.insert(global_node_id);
            }
        }

        std::unordered_map<std::size_t, std::size_t> global_to_local_node;
        auto& local_nodes = result.local_mesh.Nodes();
        local_nodes.clear();
        local_nodes.reserve(used_global_node_ids.size());

        for (const std::size_t global_node_id : used_global_node_ids) {
            Node node = global_mesh.GetNode(global_node_id);
            node.id = local_nodes.size();
            global_to_local_node[global_node_id] = node.id;
            local_nodes.push_back(std::move(node));
        }

        for (Cell& local_cell : result.local_mesh.Cells()) {
            const std::size_t global_cell_id = result.local_to_global_cell[local_cell.local_id];
            const Cell& global_cell = global_mesh.GetCell(global_cell_id);

            local_cell.node_ids.clear();
            local_cell.node_ids.reserve(global_cell.node_ids.size());

            for (const std::size_t global_node_id : global_cell.node_ids) {
                local_cell.node_ids.push_back(global_to_local_node.at(global_node_id));
            }
        }
    }

    void AppendHaloId(std::vector<std::size_t>& ids,
                      std::unordered_set<std::size_t>& seen,
                      const std::size_t value) {
        if (seen.insert(value).second) {
            ids.push_back(value);
        }
    }

    void BuildLocalFacesAndHalos(const Mesh& global_mesh,
                                 const std::vector<int>& part,
                                 const int rank,
                                 DomainDecomposition::Result& result) {
        auto& local_faces = result.local_mesh.Faces();
        auto& local_cells = result.local_mesh.Cells();

        std::map<int, HaloAccumulator> halos_by_rank;
        std::unordered_map<std::size_t, std::size_t> global_to_local_node;

        for (const Cell& local_cell : result.local_mesh.Cells()) {
            const std::size_t global_cell_id = result.local_to_global_cell[local_cell.local_id];
            const Cell& global_cell = global_mesh.GetCell(global_cell_id);

            for (std::size_t k = 0; k < global_cell.node_ids.size(); ++k) {
                global_to_local_node[global_cell.node_ids[k]] = local_cell.node_ids[k];
            }
        }

        auto remap_face_nodes_to_local = [&](Face& face) {
            for (std::size_t& node_id : face.node_ids) {
                node_id = global_to_local_node.at(node_id);
            }
        };

        auto add_local_face = [&](Face face) {
            face.id = local_faces.size();
            local_faces.push_back(face);

            local_cells[face.owner_cell_id].face_ids.push_back(face.id);
            if (face.IsInternal() || face.IsMPIBoundary()) {
                local_cells[face.neighbor_cell_id].face_ids.push_back(face.id);
            }
        };

        for (const Face& global_face : global_mesh.Faces()) {
            if (global_face.IsPhysicalBoundary()) {
                const int owner_rank = part[global_face.owner_cell_id];
                if (owner_rank != rank) {
                    continue;
                }

                Face local_face = global_face;
                local_face.owner_cell_id = result.global_to_local_cell.at(global_face.owner_cell_id);
                local_face.neighbor_cell_id = Face::k_invalid_cell_id;
                local_face.kind = FaceKind::PhysicalBoundary;
                local_face.remote_rank = -1;
                local_face.remote_cell_id = Face::k_invalid_cell_id;

                remap_face_nodes_to_local(local_face);
                add_local_face(std::move(local_face));
                continue;
            }

            if (!global_face.IsInternal()) {
                continue;
            }

            const std::size_t global_owner = global_face.owner_cell_id;
            const std::size_t global_neighbor = global_face.neighbor_cell_id;
            const int owner_rank = part[global_owner];
            const int neighbor_rank = part[global_neighbor];

            if (owner_rank == rank && neighbor_rank == rank) {
                Face local_face = global_face;
                local_face.owner_cell_id = result.global_to_local_cell.at(global_owner);
                local_face.neighbor_cell_id = result.global_to_local_cell.at(global_neighbor);
                local_face.kind = FaceKind::Interior;
                local_face.remote_rank = -1;
                local_face.remote_cell_id = Face::k_invalid_cell_id;

                remap_face_nodes_to_local(local_face);
                add_local_face(std::move(local_face));
                continue;
            }

            if (owner_rank == rank && neighbor_rank != rank) {
                Face local_face = global_face;
                local_face.owner_cell_id = result.global_to_local_cell.at(global_owner);
                local_face.neighbor_cell_id = result.global_to_local_cell.at(global_neighbor);
                local_face.kind = FaceKind::MPIBoundary;
                local_face.boundary_tag = -1;
                local_face.remote_rank = neighbor_rank;
                local_face.remote_cell_id = global_neighbor;

                remap_face_nodes_to_local(local_face);
                add_local_face(local_face);

                HaloAccumulator& halo = halos_by_rank[neighbor_rank];
                halo.remote_rank = neighbor_rank;
                AppendHaloId(halo.send_local_ids, halo.send_seen, local_face.owner_cell_id);
                AppendHaloId(halo.recv_local_ids, halo.recv_seen, local_face.neighbor_cell_id);
                continue;
            }

            if (neighbor_rank == rank && owner_rank != rank) {
                Face local_face = global_face;
                local_face.owner_cell_id = result.global_to_local_cell.at(global_neighbor);
                local_face.neighbor_cell_id = result.global_to_local_cell.at(global_owner);
                local_face.kind = FaceKind::MPIBoundary;
                local_face.boundary_tag = -1;
                local_face.remote_rank = owner_rank;
                local_face.remote_cell_id = global_owner;

                local_face.normal_x *= -1.0;
                local_face.normal_y *= -1.0;

                remap_face_nodes_to_local(local_face);
                add_local_face(local_face);

                HaloAccumulator& halo = halos_by_rank[owner_rank];
                halo.remote_rank = owner_rank;
                AppendHaloId(halo.send_local_ids, halo.send_seen, local_face.owner_cell_id);
                AppendHaloId(halo.recv_local_ids, halo.recv_seen, local_face.neighbor_cell_id);
                continue;
            }
        }

        result.halos.clear();
        result.halos.reserve(halos_by_rank.size());

        for (auto& [remote_rank, halo_acc] : halos_by_rank) {
            (void)remote_rank;
            DomainDecomposition::NeighborHalo halo;
            halo.remote_rank = halo_acc.remote_rank;
            halo.send_local_ids = std::move(halo_acc.send_local_ids);
            halo.recv_local_ids = std::move(halo_acc.recv_local_ids);
            result.halos.push_back(std::move(halo));
        }

        std::sort(result.halos.begin(), result.halos.end(),
                  [](const DomainDecomposition::NeighborHalo& a,
                     const DomainDecomposition::NeighborHalo& b) {
                      return a.remote_rank < b.remote_rank;
                  });
    }
} // namespace

DomainDecomposition::Result DomainDecomposition::Decompose(const Mesh& global_mesh,
                                                           const MPIContext& mpi) const {
    const std::vector<int> global_part = ComputePartition(global_mesh, mpi.Size());
    return BuildLocalResultForRank(global_mesh, global_part, mpi.Rank());
}

DomainDecomposition::Result DomainDecomposition::BuildLocalResultForRank(
    const Mesh& global_mesh,
    const std::vector<int>& global_part,
    const int rank
) const {
    Result result;

    result.local_mesh = Mesh(global_mesh.GetDim());

    result.global_part = global_part;
    result.owned_global_ids = CollectOwnedGlobalIds(global_mesh, result.global_part, rank);

    {
        const std::set<std::size_t> ghost_set =
            CollectGhostGlobalIds(global_mesh, result.global_part, rank);
        result.ghost_global_ids.assign(ghost_set.begin(), ghost_set.end());
    }

    CopyLocalCells(global_mesh, result.owned_global_ids, result.ghost_global_ids, result);
    CopyLocalNodesAndRemapConnectivity(result, global_mesh);
    BuildLocalFacesAndHalos(global_mesh, result.global_part, rank, result);

    result.local_mesh.Validate();
    return result;
}
