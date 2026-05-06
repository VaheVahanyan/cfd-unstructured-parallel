#include "geometry/GmshMeshBuilder.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include <gmsh.h>

#include "geometry/Cell.hpp"
#include "geometry/Face.hpp"
#include "geometry/Node.hpp"

bool GmshMeshBuilder::FaceKey::operator<(const FaceKey& other) const {
    return node_ids_sorted < other.node_ids_sorted;
}

GmshMeshBuilder::GmshSessionGuard::GmshSessionGuard(const bool finalize_on_destroy)
    : finalize_on_destroy_(finalize_on_destroy) {}

GmshMeshBuilder::GmshSessionGuard::~GmshSessionGuard() {
    if (finalize_on_destroy_) {
        gmsh::finalize();
    }
}

void GmshMeshBuilder::ValidateInputDimension(const int dim) {
    if (dim != 2) {
        throw std::invalid_argument("GmshMeshBuilder: Only 2D meshes (dim == 2) are supported in this version.");
    }
}

void GmshMeshBuilder::EnsureGmshInitialized() {
    if (!gmsh::isInitialized()) {
        throw std::runtime_error("GmshMeshBuilder: Gmsh is not initialized");
    }
}

bool GmshMeshBuilder::IsSupportedCellElementType(const int dim, const int element_type) {
    (void)dim;
    // 2: Triangle, 3: Quadrangle
    return element_type == 2 || element_type == 3;
}

bool GmshMeshBuilder::IsSupportedBoundaryElementType(const int dim, const int element_type) {
    (void)dim;
    // 1: 2-node line
    return element_type == 1;
}

int GmshMeshBuilder::GetExpectedNodeCountForElementType(const int element_type) {
    switch (element_type) {
    case 1: return 2; // Line
    case 2: return 3; // Triangle
    case 3: return 4; // Quadrangle
    default:
        throw std::runtime_error("GmshMeshBuilder: unsupported element type");
    }
}

std::vector<std::vector<int>> GmshMeshBuilder::GetLocalFacesForElementType(const int element_type) {
    switch (element_type) {
    case 2: // Triangle
        return { {0, 1}, {1, 2}, {2, 0} };
    case 3: // Quadrangle
        return { {0, 1}, {1, 2}, {2, 3}, {3, 0} };
    default:
        throw std::runtime_error("GmshMeshBuilder: unsupported cell element type");
    }
}

std::vector<std::size_t> GmshMeshBuilder::ConvertElementConnectivity(
    const std::vector<std::size_t>& element_nodes,
    const std::vector<int>& local_face_pattern
) {
    std::vector<std::size_t> face_nodes;
    face_nodes.reserve(local_face_pattern.size());

    for (const int local_id : local_face_pattern) {
        face_nodes.push_back(element_nodes[static_cast<std::size_t>(local_id)]);
    }
    return face_nodes;
}

GmshMeshBuilder::FaceKey GmshMeshBuilder::MakeFaceKey(const std::vector<std::size_t>& node_ids) {
    FaceKey key;
    key.node_ids_sorted = node_ids;
    std::sort(key.node_ids_sorted.begin(), key.node_ids_sorted.end());
    return key;
}

// --- 2D Vector Math ---

GmshMeshBuilder::Vec2 GmshMeshBuilder::Add(const Vec2& a, const Vec2& b) {
    return Vec2{a.x + b.x, a.y + b.y};
}

GmshMeshBuilder::Vec2 GmshMeshBuilder::Subtract(const Vec2& a, const Vec2& b) {
    return Vec2{a.x - b.x, a.y - b.y};
}

GmshMeshBuilder::Vec2 GmshMeshBuilder::Multiply(const double scalar, const Vec2& v) {
    return Vec2{scalar * v.x, scalar * v.y};
}

GmshMeshBuilder::Vec2 GmshMeshBuilder::Divide(const Vec2& v, const double scalar) {
    return Vec2{v.x / scalar, v.y / scalar};
}

double GmshMeshBuilder::Dot(const Vec2& a, const Vec2& b) {
    return a.x * b.x + a.y * b.y;
}

double GmshMeshBuilder::Norm(const Vec2& v) {
    return std::sqrt(Dot(v, v));
}

GmshMeshBuilder::Vec2 GmshMeshBuilder::Normalize(const Vec2& v) {
    const double norm = Norm(v);
    if (!(norm > 0.0)) {
        throw std::runtime_error("GmshMeshBuilder: cannot normalize zero vector");
    }
    return Divide(v, norm);
}

GmshMeshBuilder::Vec2 GmshMeshBuilder::ToVec2(const Node& node) {
    return Vec2{node.x, node.y};
}

// --- Geometry ---

GmshMeshBuilder::Vec2 GmshMeshBuilder::ComputePolygonCenter(
    const std::vector<std::size_t>& node_ids,
    const std::vector<Node>& nodes
) {
    Vec2 center{0.0, 0.0};
    for (const std::size_t node_id : node_ids) {
        center = Add(center, ToVec2(nodes[node_id]));
    }
    return Divide(center, static_cast<double>(node_ids.size()));
}

double GmshMeshBuilder::ComputeEdgeLength(
    const std::vector<std::size_t>& node_ids,
    const std::vector<Node>& nodes
) {
    return Norm(Subtract(ToVec2(nodes[node_ids[1]]), ToVec2(nodes[node_ids[0]])));
}

double GmshMeshBuilder::ComputeCellArea2D(
    const std::vector<std::size_t>& node_ids,
    const std::vector<Node>& nodes
) {
    double twice_area = 0.0;
    for (std::size_t i = 0; i < node_ids.size(); ++i) {
        const Node& a = nodes[node_ids[i]];
        const Node& b = nodes[node_ids[(i + 1) % node_ids.size()]];
        twice_area += a.x * b.y - b.x * a.y;
    }
    return 0.5 * std::abs(twice_area);
}

// --- Builder Methods ---

std::map<GmshMeshBuilder::FaceKey, GmshMeshBuilder::BoundaryFaceRecord>
GmshMeshBuilder::BuildBoundaryFaceMap(
    const int dim,
    const std::unordered_map<std::size_t, std::size_t>& gmsh_to_internal_node
) {
    std::map<FaceKey, BoundaryFaceRecord> boundary_faces;

    std::vector<std::pair<int, int>> physical_groups;
    gmsh::model::getPhysicalGroups(physical_groups, dim - 1); // 1D boundaries for 2D

    for (const auto& [physical_dim, physical_tag] : physical_groups) {
        if (physical_dim != dim - 1) continue;

        std::vector<int> entity_tags;
        gmsh::model::getEntitiesForPhysicalGroup(physical_dim, physical_tag, entity_tags);

        for (const int entity_tag : entity_tags) {
            std::vector<int> element_types;
            std::vector<std::vector<std::size_t>> element_tags, element_node_tags;

            gmsh::model::mesh::getElements(element_types, element_tags, element_node_tags, physical_dim, entity_tag);

            for (std::size_t block_id = 0; block_id < element_types.size(); ++block_id) {
                const int element_type = element_types[block_id];
                if (!IsSupportedBoundaryElementType(dim, element_type)) continue;

                const int node_count = GetExpectedNodeCountForElementType(element_type);
                const std::vector<std::size_t>& flat_nodes = element_node_tags[block_id];
                const std::size_t element_count = flat_nodes.size() / static_cast<std::size_t>(node_count);

                for (std::size_t element_id = 0; element_id < element_count; ++element_id) {
                    std::vector<std::size_t> face_nodes;
                    for (int local_id = 0; local_id < node_count; ++local_id) {
                        const std::size_t gmsh_node_tag = flat_nodes[element_id * node_count + local_id];
                        face_nodes.push_back(gmsh_to_internal_node.at(gmsh_node_tag));
                    }
                    boundary_faces[MakeFaceKey(face_nodes)] = BoundaryFaceRecord{physical_tag};
                }
            }
        }
    }
    return boundary_faces;
}

void GmshMeshBuilder::BuildNodes(
    Mesh& mesh,
    std::unordered_map<std::size_t, std::size_t>& gmsh_to_internal_node
) {
    std::vector<std::size_t> node_tags;
    std::vector<double> coords, params;

    gmsh::model::mesh::getNodes(node_tags, coords, params, -1, -1, false, false);

    auto& nodes = mesh.Nodes();
    nodes.reserve(node_tags.size());

    for (std::size_t i = 0; i < node_tags.size(); ++i) {
        Node node;
        node.id = nodes.size();
        node.x = coords[3 * i + 0]; // Gmsh всегда возвращает 3D координаты
        node.y = coords[3 * i + 1];

        gmsh_to_internal_node[node_tags[i]] = node.id;
        nodes.push_back(node);
    }
}

void GmshMeshBuilder::BuildCellsAndFaces(
    Mesh& mesh,
    const int dim,
    const std::unordered_map<std::size_t, std::size_t>& gmsh_to_internal_node,
    const std::map<FaceKey, BoundaryFaceRecord>& boundary_faces
) {
    std::vector<int> element_types;
    std::vector<std::vector<std::size_t>> element_tags, element_node_tags;

    gmsh::model::mesh::getElements(element_types, element_tags, element_node_tags, dim, -1);

    auto& cells = mesh.Cells();
    auto& faces = mesh.Faces();
    const auto& nodes = mesh.Nodes();

    std::map<FaceKey, std::size_t> face_key_to_face_id;

    for (std::size_t block_id = 0; block_id < element_types.size(); ++block_id) {
        const int element_type = element_types[block_id];
        if (!IsSupportedCellElementType(dim, element_type)) continue;

        const int node_count = GetExpectedNodeCountForElementType(element_type);
        const std::vector<std::size_t>& flat_nodes = element_node_tags[block_id];
        const std::size_t element_count = flat_nodes.size() / static_cast<std::size_t>(node_count);
        const std::vector<std::vector<int>> local_faces = GetLocalFacesForElementType(element_type);

        for (std::size_t element_id = 0; element_id < element_count; ++element_id) {
            Cell cell;
            cell.id = cells.size();
            cell.local_id = cells.size();
            cell.node_ids.reserve(static_cast<std::size_t>(node_count));

            for (int local_id = 0; local_id < node_count; ++local_id) {
                const std::size_t gmsh_node_tag = flat_nodes[element_id * node_count + local_id];
                cell.node_ids.push_back(gmsh_to_internal_node.at(gmsh_node_tag));
            }

            const Vec2 cell_center = ComputePolygonCenter(cell.node_ids, nodes);
            cell.center_x = cell_center.x;
            cell.center_y = cell_center.y;
            cell.volume = ComputeCellArea2D(cell.node_ids, nodes);

            cells.push_back(cell);

            for (const std::vector<int>& local_face_pattern : local_faces) {
                const std::vector<std::size_t> face_nodes = ConvertElementConnectivity(cell.node_ids, local_face_pattern);
                const FaceKey key = MakeFaceKey(face_nodes);
                const auto face_it = face_key_to_face_id.find(key);

                if (face_it == face_key_to_face_id.end()) {
                    Face face;
                    face.id = faces.size();
                    face.node_ids = face_nodes;
                    face.owner_cell_id = cell.local_id;
                    face.neighbor_cell_id = Face::k_invalid_cell_id;
                    face.kind = FaceKind::PhysicalBoundary;
                    face.remote_rank = -1;
                    face.remote_cell_id = Face::k_invalid_cell_id;

                    const Vec2 face_center = ComputePolygonCenter(face.node_ids, nodes);
                    face.center_x = face_center.x;
                    face.center_y = face_center.y;
                    face.measure = ComputeEdgeLength(face.node_ids, nodes);

                    const auto boundary_it = boundary_faces.find(key);
                    face.boundary_tag = (boundary_it != boundary_faces.end()) ? boundary_it->second.boundary_tag : -1;

                    faces.push_back(face);
                    face_key_to_face_id[key] = face.id;
                    cells[cell.local_id].face_ids.push_back(face.id);
                } else {
                    Face& face = faces[face_it->second];
                    face.neighbor_cell_id = cell.local_id;
                    face.kind = FaceKind::Interior;
                    face.boundary_tag = -1;
                    cells[cell.local_id].face_ids.push_back(face.id);
                }
            }
        }
    }
}

void GmshMeshBuilder::FinalizeFaceNormals(Mesh& mesh) {
    auto& faces = mesh.Faces();
    const auto& cells = mesh.Cells();
    const auto& nodes = mesh.Nodes();

    for (Face& face : faces) {
        const Vec2 owner_center{ cells[face.owner_cell_id].center_x, cells[face.owner_cell_id].center_y };

        Vec2 direction;
        if (face.IsInternal()) {
            const Vec2 neighbor_center{ cells[face.neighbor_cell_id].center_x, cells[face.neighbor_cell_id].center_y };
            direction = Subtract(neighbor_center, owner_center);
        } else {
            const Vec2 face_center{.x = face.center_x, .y = face.center_y};
            direction = Subtract(face_center, owner_center);
        }

        const Vec2 a = ToVec2(nodes[face.node_ids[0]]);
        const Vec2 b = ToVec2(nodes[face.node_ids[1]]);
        const Vec2 tangent = Subtract(b, a);

        Vec2 normal{.x = -tangent.y, .y = tangent.x};

        if (Dot(normal, direction) < 0.0) {
            normal = Multiply(-1.0, normal);
        }

        normal = Normalize(normal);
        face.normal_x = normal.x;
        face.normal_y = normal.y;
    }
}

Mesh GmshMeshBuilder::BuildFromCurrentModel(const int dim) {
    ValidateInputDimension(dim);
    EnsureGmshInitialized();

    // Передаем dim в Mesh только если его конструктор этого требует.
    // Если ты удалил размерность из Mesh, замени на `Mesh mesh;`
    Mesh mesh(dim);

    std::unordered_map<std::size_t, std::size_t> gmsh_to_internal_node;
    BuildNodes(mesh, gmsh_to_internal_node);

    const std::map<FaceKey, BoundaryFaceRecord> boundary_faces = BuildBoundaryFaceMap(dim, gmsh_to_internal_node);

    BuildCellsAndFaces(mesh, dim, gmsh_to_internal_node, boundary_faces);
    FinalizeFaceNormals(mesh);

    mesh.SetOwnedCellCount(mesh.GetCellCount());
    mesh.SetGhostCellCount(0);

    mesh.Validate();
    return mesh;
}

Mesh GmshMeshBuilder::BuildFromFile(const std::string& file_path, const int dim) {
    ValidateInputDimension(dim);

    const bool was_initialized = gmsh::isInitialized();
    if (!was_initialized) gmsh::initialize();

    GmshSessionGuard session_guard(!was_initialized);
    gmsh::clear();
    gmsh::open(file_path);

    return BuildFromCurrentModel(dim);
}

Mesh GmshMeshBuilder::BuildFromGeoFile(const std::string& file_path, const int dim) {
    ValidateInputDimension(dim);

    const bool was_initialized = gmsh::isInitialized();
    if (!was_initialized) gmsh::initialize();

    GmshSessionGuard session_guard(!was_initialized);
    gmsh::clear();
    gmsh::open(file_path);
    gmsh::model::mesh::generate(dim);

    return BuildFromCurrentModel(dim);
}