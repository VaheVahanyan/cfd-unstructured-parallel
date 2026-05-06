#include "reconstruction/P1Reconstruction.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "data/Workspace.hpp"
#include "geometry/Cell.hpp"
#include "geometry/Face.hpp"
#include "geometry/Mesh.hpp"

PrimitiveCell P1Reconstruction::LoadCellPrimitive(const Workspace& workspace, const std::size_t cell_id) const {
    const auto& W = workspace.W();
    PrimitiveCell state;
    state.rho = W(cell_id, Workspace::k_rho);
    state.u = W(cell_id, Workspace::k_u);
    state.v = W(cell_id, Workspace::k_v);
    state.P = W(cell_id, Workspace::k_p);
    return state;
}

P1Reconstruction::PrimitiveGradient P1Reconstruction::ComputeUnlimitedGradient(
    const Mesh& mesh,
    const Workspace& workspace,
    const Cell& cell
) const {
    PrimitiveGradient grad{};

    const std::size_t begin = mesh.GetCellNeighborBegin(cell.local_id);
    const std::size_t end = mesh.GetCellNeighborEnd(cell.local_id);

    if (begin == end) {
        return grad;
    }

    const PrimitiveCell wc = LoadCellPrimitive(workspace, cell.local_id);

    double a11 = 0.0;
    double a12 = 0.0;
    double a22 = 0.0;

    double b1_rho = 0.0;
    double b2_rho = 0.0;
    double b1_u = 0.0;
    double b2_u = 0.0;
    double b1_v = 0.0;
    double b2_v = 0.0;
    double b1_P = 0.0;
    double b2_P = 0.0;

    for (std::size_t offset = begin; offset < end; ++offset) {
        const std::size_t n_id = mesh.GetCellNeighborId(offset);
        const Cell& neighbor = mesh.GetCell(n_id);
        const PrimitiveCell wn = LoadCellPrimitive(workspace, neighbor.local_id);

        const double dx = neighbor.center_x - cell.center_x;
        const double dy = neighbor.center_y - cell.center_y;

        a11 += dx * dx;
        a12 += dx * dy;
        a22 += dy * dy;

        b1_rho += dx * (wn.rho - wc.rho);
        b2_rho += dy * (wn.rho - wc.rho);

        b1_u += dx * (wn.u - wc.u);
        b2_u += dy * (wn.u - wc.u);

        b1_v += dx * (wn.v - wc.v);
        b2_v += dy * (wn.v - wc.v);

        b1_P += dx * (wn.P - wc.P);
        b2_P += dy * (wn.P - wc.P);
    }

    const double det = a11 * a22 - a12 * a12;

    if (std::abs(det) > 1e-14) {
        grad.dx.rho = (b1_rho * a22 - b2_rho * a12) / det;
        grad.dy.rho = (-b1_rho * a12 + b2_rho * a11) / det;

        grad.dx.u = (b1_u * a22 - b2_u * a12) / det;
        grad.dy.u = (-b1_u * a12 + b2_u * a11) / det;

        grad.dx.v = (b1_v * a22 - b2_v * a12) / det;
        grad.dy.v = (-b1_v * a12 + b2_v * a11) / det;

        grad.dx.P = (b1_P * a22 - b2_P * a12) / det;
        grad.dy.P = (-b1_P * a12 + b2_P * a11) / det;
    }

    return grad;
}

double P1Reconstruction::ComputeBarthJespersenPhi(const double w_cell, const double w_min, const double w_max,
                                                  const double w_face_candidate) const {
    constexpr double eps = 1e-14;
    const double delta = w_face_candidate - w_cell;

    if (delta > eps) {
        return std::min(1.0, (w_max - w_cell) / delta);
    }
    if (delta < -eps) {
        return std::min(1.0, (w_min - w_cell) / delta);
    }
    return 1.0;
}

double P1Reconstruction::ComputeLimiter(const Mesh& mesh,
                                        const Workspace& workspace,
                                        const Cell& cell,
                                        const PrimitiveGradient& grad) const {
    PrimitiveCell w_min = LoadCellPrimitive(workspace, cell.local_id);
    PrimitiveCell w_max = w_min;
    const PrimitiveCell wc = w_min;

    const std::size_t begin = mesh.GetCellNeighborBegin(cell.local_id);
    const std::size_t end = mesh.GetCellNeighborEnd(cell.local_id);

    for (std::size_t offset = begin; offset < end; ++offset) {
        const std::size_t n_id = mesh.GetCellNeighborId(offset);
        const PrimitiveCell wn = LoadCellPrimitive(workspace, n_id);

        w_min.rho = std::min(w_min.rho, wn.rho);
        w_max.rho = std::max(w_max.rho, wn.rho);

        w_min.u = std::min(w_min.u, wn.u);
        w_max.u = std::max(w_max.u, wn.u);

        w_min.v = std::min(w_min.v, wn.v);
        w_max.v = std::max(w_max.v, wn.v);

        w_min.P = std::min(w_min.P, wn.P);
        w_max.P = std::max(w_max.P, wn.P);
    }

    double phi = 1.0;

    for (const std::size_t face_id : cell.face_ids) {
        const Face& face = mesh.GetFace(face_id);

        const double dx = face.center_x - cell.center_x;
        const double dy = face.center_y - cell.center_y;

        const PrimitiveCell w_face{
            wc.rho + grad.dx.rho * dx + grad.dy.rho * dy,
            wc.u + grad.dx.u * dx + grad.dy.u * dy,
            wc.v + grad.dx.v * dx + grad.dy.v * dy,
            wc.P + grad.dx.P * dx + grad.dy.P * dy
        };

        phi = std::min(phi, ComputeBarthJespersenPhi(wc.rho, w_min.rho, w_max.rho, w_face.rho));
        phi = std::min(phi, ComputeBarthJespersenPhi(wc.u, w_min.u, w_max.u, w_face.u));
        phi = std::min(phi, ComputeBarthJespersenPhi(wc.v, w_min.v, w_max.v, w_face.v));
        phi = std::min(phi, ComputeBarthJespersenPhi(wc.P, w_min.P, w_max.P, w_face.P));
    }

    return std::clamp(phi, 0.0, 1.0);
}

void P1Reconstruction::ComputeGradients(const Mesh& mesh, Workspace& workspace) const {
    workspace.ZeroGradients();

    auto& grad_rho = workspace.GradRho();
    auto& grad_u = workspace.GradU();
    auto& grad_v = workspace.GradV();
    auto& grad_p = workspace.GradP();

    const std::size_t cell_count = mesh.GetCellCount();

#pragma omp parallel for schedule(static) default(none) shared(mesh, workspace, grad_rho, grad_u, grad_v, grad_p, cell_count)
    for (std::size_t cell_id = 0; cell_id < cell_count; ++cell_id) {
        const Cell& cell = mesh.GetCell(cell_id);

        const PrimitiveGradient grad_unlimited =
            ComputeUnlimitedGradient(mesh, workspace, cell);

        const double phi =
            ComputeLimiter(mesh, workspace, cell, grad_unlimited);

        grad_rho(cell_id, 0) = phi * grad_unlimited.dx.rho;
        grad_rho(cell_id, 1) = phi * grad_unlimited.dy.rho;

        grad_u(cell_id, 0) = phi * grad_unlimited.dx.u;
        grad_u(cell_id, 1) = phi * grad_unlimited.dy.u;

        grad_v(cell_id, 0) = phi * grad_unlimited.dx.v;
        grad_v(cell_id, 1) = phi * grad_unlimited.dy.v;

        grad_p(cell_id, 0) = phi * grad_unlimited.dx.P;
        grad_p(cell_id, 1) = phi * grad_unlimited.dy.P;
    }
}

void P1Reconstruction::ReconstructInteriorFace(const Mesh& mesh, const Workspace& workspace, const Face& face,
                                               PrimitiveCell& owner_state, PrimitiveCell& neighbor_state) const {
    const Cell& owner = mesh.GetCell(face.owner_cell_id);
    const Cell& neighbor = mesh.GetCell(face.neighbor_cell_id);

    const PrimitiveCell wc_owner = LoadCellPrimitive(workspace, owner.local_id);
    const PrimitiveCell wc_neighbor = LoadCellPrimitive(workspace, neighbor.local_id);

    const auto& grad_rho = workspace.GradRho();
    const auto& grad_u = workspace.GradU();
    const auto& grad_v = workspace.GradV();
    const auto& grad_p = workspace.GradP();

    const double dx_o = face.center_x - owner.center_x;
    const double dy_o = face.center_y - owner.center_y;

    owner_state.rho = wc_owner.rho + grad_rho(owner.local_id, 0) * dx_o + grad_rho(owner.local_id, 1) * dy_o;
    owner_state.u = wc_owner.u + grad_u(owner.local_id, 0) * dx_o + grad_u(owner.local_id, 1) * dy_o;
    owner_state.v = wc_owner.v + grad_v(owner.local_id, 0) * dx_o + grad_v(owner.local_id, 1) * dy_o;
    owner_state.P = wc_owner.P + grad_p(owner.local_id, 0) * dx_o + grad_p(owner.local_id, 1) * dy_o;

    const double dx_n = face.center_x - neighbor.center_x;
    const double dy_n = face.center_y - neighbor.center_y;

    neighbor_state.rho = wc_neighbor.rho + grad_rho(neighbor.local_id, 0) * dx_n + grad_rho(neighbor.local_id, 1) *
        dy_n;
    neighbor_state.u = wc_neighbor.u + grad_u(neighbor.local_id, 0) * dx_n + grad_u(neighbor.local_id, 1) * dy_n;
    neighbor_state.v = wc_neighbor.v + grad_v(neighbor.local_id, 0) * dx_n + grad_v(neighbor.local_id, 1) * dy_n;
    neighbor_state.P = wc_neighbor.P + grad_p(neighbor.local_id, 0) * dx_n + grad_p(neighbor.local_id, 1) * dy_n;
}

void P1Reconstruction::ReconstructBoundaryFaceInterior(const Mesh& mesh, const Workspace& workspace, const Face& face,
                                                       PrimitiveCell& interior_state) const {
    const Cell& owner = mesh.GetCell(face.owner_cell_id);
    const PrimitiveCell wc_owner = LoadCellPrimitive(workspace, owner.local_id);

    const double dx_o = face.center_x - owner.center_x;
    const double dy_o = face.center_y - owner.center_y;

    interior_state.rho = wc_owner.rho + workspace.GradRho()(owner.local_id, 0) * dx_o + workspace.
        GradRho()(owner.local_id, 1) * dy_o;
    interior_state.u = wc_owner.u + workspace.GradU()(owner.local_id, 0) * dx_o + workspace.GradU()(owner.local_id, 1) *
        dy_o;
    interior_state.v = wc_owner.v + workspace.GradV()(owner.local_id, 0) * dx_o + workspace.GradV()(owner.local_id, 1) *
        dy_o;
    interior_state.P = wc_owner.P + workspace.GradP()(owner.local_id, 0) * dx_o + workspace.GradP()(owner.local_id, 1) *
        dy_o;
}
