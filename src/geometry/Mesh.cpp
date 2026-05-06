#include "geometry/Mesh.hpp"

#include <cmath>
#include <string>
#include <algorithm>

Mesh::Mesh(const int dim) : dim_(dim) {
    ValidateDimension();
}

int Mesh::GetDim() const {
    return dim_;
}

std::size_t Mesh::GetNodeCount() const {
    return nodes_.size();
}

std::size_t Mesh::GetFaceCount() const {
    return faces_.size();
}

std::size_t Mesh::GetCellCount() const {
    return cells_.size();
}

std::size_t Mesh::GetOwnedCellCount() const {
    return owned_cell_count_;
}

std::size_t Mesh::GetGhostCellCount() const {
    return ghost_cell_count_;
}

void Mesh::SetOwnedCellCount(const std::size_t count) {
    owned_cell_count_ = count;
}

void Mesh::SetGhostCellCount(const std::size_t count) {
    ghost_cell_count_ = count;
}

std::vector<Node>& Mesh::Nodes() {
    return nodes_;
}

const std::vector<Node>& Mesh::Nodes() const {
    return nodes_;
}

std::vector<Face>& Mesh::Faces() {
    return faces_;
}

const std::vector<Face>& Mesh::Faces() const {
    return faces_;
}

std::vector<Cell>& Mesh::Cells() {
    return cells_;
}

const std::vector<Cell>& Mesh::Cells() const {
    return cells_;
}

Node& Mesh::GetNode(const std::size_t node_id) {
    if (node_id >= nodes_.size()) {
        throw std::out_of_range("Mesh::GetNode: node_id is out of range");
    }
    return nodes_[node_id];
}

const Node& Mesh::GetNode(const std::size_t node_id) const {
    if (node_id >= nodes_.size()) {
        throw std::out_of_range("Mesh::GetNode: node_id is out of range");
    }
    return nodes_[node_id];
}

Face& Mesh::GetFace(const std::size_t face_id) {
    if (face_id >= faces_.size()) {
        throw std::out_of_range("Mesh::GetFace: face_id is out of range");
    }
    return faces_[face_id];
}

const Face& Mesh::GetFace(const std::size_t face_id) const {
    if (face_id >= faces_.size()) {
        throw std::out_of_range("Mesh::GetFace: face_id is out of range");
    }
    return faces_[face_id];
}

Cell& Mesh::GetCell(const std::size_t cell_id) {
    if (cell_id >= cells_.size()) {
        throw std::out_of_range("Mesh::GetCell: cell_id is out of range");
    }
    return cells_[cell_id];
}

const Cell& Mesh::GetCell(const std::size_t cell_id) const {
    if (cell_id >= cells_.size()) {
        throw std::out_of_range("Mesh::GetCell: cell_id is out of range");
    }
    return cells_[cell_id];
}

bool Mesh::IsPhysicalBoundaryFace(const std::size_t face_id) const {
    return GetFace(face_id).IsPhysicalBoundary();
}

bool Mesh::IsInternalFace(const std::size_t face_id) const {
    return GetFace(face_id).IsInternal();
}

bool Mesh::IsMPIBoundaryFace(const std::size_t face_id) const {
    return GetFace(face_id).IsMPIBoundary();
}

void Mesh::Clear() {
    nodes_.clear();
    faces_.clear();
    cells_.clear();
    owned_cell_count_ = 0;
    ghost_cell_count_ = 0;
    color_to_face_ids_.clear();
    num_colors_ = 0;
    cell_neighbor_offsets_.clear();
    cell_neighbor_ids_.clear();
}

void Mesh::Validate() {
    ValidateDimension();
    ValidateNodeIds();
    ValidateFaceIds();
    ValidateCellLocalIds();
    ValidateFaceConnectivity();
    ValidateCellConnectivity();
    ValidateGeometry();

    if (owned_cell_count_ + ghost_cell_count_ > cells_.size()) {
        throw std::runtime_error("Mesh::Validate: owned+ghost cell counts exceed total cell count");
    }

    ComputeFaceColors();
    BuildCellNeighbors();
    ComputeFaceColors();
}

void Mesh::ValidateDimension() const {
    if (dim_ < 1 || dim_ > 2) {
        throw std::runtime_error("Mesh::Validate: dim must be 1..2");
    }
}

void Mesh::ValidateNodeIds() const {
    for (std::size_t node_id = 0; node_id < nodes_.size(); ++node_id) {
        if (nodes_[node_id].id != node_id) {
            throw std::runtime_error("Mesh::Validate: node id does not match storage index");
        }
    }
}

void Mesh::ValidateFaceIds() const {
    for (std::size_t face_id = 0; face_id < faces_.size(); ++face_id) {
        if (faces_[face_id].id != face_id) {
            throw std::runtime_error("Mesh::Validate: face id does not match storage index");
        }
    }
}

void Mesh::ValidateCellLocalIds() const {
    for (std::size_t cell_id = 0; cell_id < cells_.size(); ++cell_id) {
        if (cells_[cell_id].local_id != cell_id) {
            throw std::runtime_error("Mesh::Validate: cell local_id does not match storage index");
        }
    }
}

void Mesh::ValidateFaceConnectivity() const {
    for (const Face& face : faces_) {
        if (face.node_ids.empty()) {
            throw std::runtime_error("Mesh::Validate: face has empty node_ids");
        }

        if (face.owner_cell_id == Face::k_invalid_cell_id) {
            throw std::runtime_error("Mesh::Validate: face has invalid owner_cell_id");
        }

        if (face.owner_cell_id >= cells_.size()) {
            throw std::runtime_error("Mesh::Validate: face owner_cell_id is out of range");
        }

        if ((face.IsInternal() || face.IsMPIBoundary()) &&
            face.neighbor_cell_id >= cells_.size()) {
            throw std::runtime_error("Mesh::Validate: face neighbor_cell_id is out of range");
        }

        if ((face.IsInternal() || face.IsMPIBoundary()) &&
            face.neighbor_cell_id == face.owner_cell_id) {
            throw std::runtime_error("Mesh::Validate: face owner and neighbor are identical");
        }

        if (face.IsPhysicalBoundary()) {
            if (face.neighbor_cell_id != Face::k_invalid_cell_id) {
                throw std::runtime_error("Mesh::Validate: physical boundary face must have invalid neighbor");
            }
            if (face.boundary_tag < 0) {
                throw std::runtime_error("Mesh::Validate: physical boundary face must have non-negative boundary_tag");
            }
            if (face.remote_rank >= 0) {
                throw std::runtime_error("Mesh::Validate: physical boundary face must not have remote_rank");
            }
        }

        if (face.IsInternal()) {
            if (face.boundary_tag >= 0) {
                throw std::runtime_error("Mesh::Validate: interior face must not have boundary_tag");
            }
            if (face.remote_rank >= 0) {
                throw std::runtime_error("Mesh::Validate: interior face must not have remote_rank");
            }
        }

        if (face.IsMPIBoundary()) {
            if (face.boundary_tag >= 0) {
                throw std::runtime_error("Mesh::Validate: MPI boundary face must not have boundary_tag");
            }
            if (face.remote_rank < 0) {
                throw std::runtime_error("Mesh::Validate: MPI boundary face must have remote_rank");
            }
            if (face.neighbor_cell_id == Face::k_invalid_cell_id) {
                throw std::runtime_error("Mesh::Validate: MPI boundary face must have ghost neighbor");
            }
        }

        for (const std::size_t node_id : face.node_ids) {
            if (node_id >= nodes_.size()) {
                throw std::runtime_error("Mesh::Validate: face node_id is out of range");
            }
        }
    }
}

void Mesh::ValidateCellConnectivity() const {
    for (const Cell& cell : cells_) {
        if (cell.node_ids.empty()) {
            throw std::runtime_error("Mesh::Validate: cell has empty node_ids");
        }

        if (cell.face_ids.empty()) {
            throw std::runtime_error("Mesh::Validate: cell has empty face_ids");
        }

        for (const std::size_t node_id : cell.node_ids) {
            if (node_id >= nodes_.size()) {
                throw std::runtime_error("Mesh::Validate: cell node_id is out of range");
            }
        }

        for (const std::size_t face_id : cell.face_ids) {
            if (face_id >= faces_.size()) {
                throw std::runtime_error("Mesh::Validate: cell face_id is out of range");
            }

            const Face& face = faces_[face_id];
            if (face.owner_cell_id != cell.local_id && face.neighbor_cell_id != cell.local_id) {
                throw std::runtime_error("Mesh::Validate: cell does not belong to one of its listed faces");
            }
        }
    }
}

void Mesh::ValidateGeometry() const {
    constexpr double eps = 1e-14;

    for (const Face& face : faces_) {
        if (!(face.measure > 0.0)) {
            throw std::runtime_error("Mesh::Validate: face measure must be > 0");
        }

        const double normal_norm =
            std::sqrt(face.normal_x * face.normal_x +
                      face.normal_y * face.normal_y);

        if (std::abs(normal_norm - 1.0) > 1e-10) {
            throw std::runtime_error("Mesh::Validate: face normal must be unit-length");
        }
    }

    for (const Cell& cell : cells_) {
        if (!(cell.volume > 0.0)) {
            throw std::runtime_error("Mesh::Validate: cell volume must be > 0");
        }

        const bool center_is_finite =
            std::isfinite(cell.center_x) &&
            std::isfinite(cell.center_y);

        if (!center_is_finite) {
            throw std::runtime_error("Mesh::Validate: cell center contains non-finite value");
        }

        const bool is_owned = cell.local_id < owned_cell_count_;

        if (is_owned) {
            if (cell.face_ids.size() < 2) {
                throw std::runtime_error("Mesh::Validate: owned cell must have at least 2 faces");
            }

            if (dim_ == 2 && cell.face_ids.size() < 3) {
                throw std::runtime_error("Mesh::Validate: 2D owned cell must have at least 3 faces");
            }
        } else {
            if (cell.face_ids.empty()) {
                throw std::runtime_error("Mesh::Validate: ghost cell must have at least 1 local face");
            }
        }
    }

    for (const Node& node : nodes_) {
        const bool coords_are_finite =
            std::isfinite(node.x) &&
            std::isfinite(node.y);

        if (!coords_are_finite) {
            throw std::runtime_error("Mesh::Validate: node coordinates contain non-finite value");
        }
    }

    for (const Face& face : faces_) {
        const bool center_is_finite =
            std::isfinite(face.center_x) &&
            std::isfinite(face.center_y);

        if (!center_is_finite) {
            throw std::runtime_error("Mesh::Validate: face center contains non-finite value");
        }

        const bool normal_is_finite =
            std::isfinite(face.normal_x) &&
            std::isfinite(face.normal_y);

        if (!normal_is_finite) {
            throw std::runtime_error("Mesh::Validate: face normal contains non-finite value");
        }

        if (std::abs(face.normal_x) < eps &&
            std::abs(face.normal_y) < eps) {
            throw std::runtime_error("Mesh::Validate: face normal must be non-zero");
        }
    }
}


void Mesh::ComputeFaceColors() {
    std::vector face_color(faces_.size(), -1);
    int max_color = -1;

    for (std::size_t i = 0; i < faces_.size(); ++i) {
        const Face& face = faces_[i];

        std::vector<bool> used_colors;

        auto check_cell = [&](const std::size_t cell_id) {
            if (cell_id == Face::k_invalid_cell_id) {
                return;
            }
            const Cell& cell = cells_[cell_id];
            for (const std::size_t cell_face_id : cell.face_ids) {
                if (cell_face_id == i) {
                    continue;
                }
                const int c = face_color[cell_face_id];
                if (c != -1) {
                    if (static_cast<std::size_t>(c) >= used_colors.size()) {
                        used_colors.resize(static_cast<std::size_t>(c + 1), false);
                    }
                    used_colors[static_cast<std::size_t>(c)] = true;
                }
            }
        };

        check_cell(face.owner_cell_id);
        check_cell(face.neighbor_cell_id);

        int chosen_color = 0;
        while (static_cast<std::size_t>(chosen_color) < used_colors.size() &&
            used_colors[static_cast<std::size_t>(chosen_color)]) {
            chosen_color++;
        }

        face_color[i] = chosen_color;
        max_color = std::max(max_color, chosen_color);
    }

    num_colors_ = max_color + 1;

    color_to_face_ids_.clear();
    color_to_face_ids_.resize(static_cast<std::size_t>(num_colors_));

    for (std::size_t i = 0; i < faces_.size(); ++i) {
        color_to_face_ids_[static_cast<std::size_t>(face_color[i])].push_back(i);
    }
}

std::size_t Mesh::GetCellNeighborBegin(const std::size_t cell_id) const {
    if (cell_id >= cells_.size()) {
        throw std::out_of_range("Mesh::GetCellNeighborBegin: cell_id is out of range");
    }

    return cell_neighbor_offsets_[cell_id];
}

std::size_t Mesh::GetCellNeighborEnd(const std::size_t cell_id) const {
    if (cell_id >= cells_.size()) {
        throw std::out_of_range("Mesh::GetCellNeighborEnd: cell_id is out of range");
    }

    return cell_neighbor_offsets_[cell_id + 1];
}

std::size_t Mesh::GetCellNeighborId(const std::size_t neighbor_offset) const {
    if (neighbor_offset >= cell_neighbor_ids_.size()) {
        throw std::out_of_range("Mesh::GetCellNeighborId: neighbor_offset is out of range");
    }

    return cell_neighbor_ids_[neighbor_offset];
}

std::size_t Mesh::GetCellNeighborCount(const std::size_t cell_id) const {
    if (cell_id >= cells_.size()) {
        throw std::out_of_range("Mesh::GetCellNeighborCount: cell_id is out of range");
    }

    return cell_neighbor_offsets_[cell_id + 1] - cell_neighbor_offsets_[cell_id];
}

void Mesh::BuildCellNeighbors() {
    cell_neighbor_offsets_.clear();
    cell_neighbor_ids_.clear();

    cell_neighbor_offsets_.resize(cells_.size() + 1, 0);

    std::vector<std::size_t> neighbors;
    neighbors.reserve(8);

    for (std::size_t cell_id = 0; cell_id < cells_.size(); ++cell_id) {
        neighbors.clear();

        const Cell& cell = cells_[cell_id];

        for (const std::size_t face_id : cell.face_ids) {
            const Face& face = faces_[face_id];

            if (face.IsPhysicalBoundary()) {
                continue;
            }

            std::size_t neighbor_id = Face::k_invalid_cell_id;

            if (face.owner_cell_id == cell.local_id) {
                neighbor_id = face.neighbor_cell_id;
            } else if (face.neighbor_cell_id == cell.local_id) {
                neighbor_id = face.owner_cell_id;
            }

            if (neighbor_id == Face::k_invalid_cell_id || neighbor_id == cell.local_id) {
                continue;
            }

            const bool already_added =
                std::find(neighbors.begin(), neighbors.end(), neighbor_id) != neighbors.end();

            if (!already_added) {
                neighbors.push_back(neighbor_id);
            }
        }

        cell_neighbor_offsets_[cell_id] = cell_neighbor_ids_.size();
        cell_neighbor_ids_.insert(cell_neighbor_ids_.end(), neighbors.begin(), neighbors.end());
    }

    cell_neighbor_offsets_[cells_.size()] = cell_neighbor_ids_.size();
}
