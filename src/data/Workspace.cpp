#include "data/Workspace.hpp"

#include "geometry/Mesh.hpp"

void Workspace::ResizeFrom(const Mesh& mesh) {
    const std::size_t n_cells = mesh.GetCellCount();

    if (n_cells == n_cells_ && IsAllocated()) {
        return;
    }

    Allocate(n_cells);
}

// -------------------- cell-centered getters --------------------

xt::xtensor<double, 2>& Workspace::W() {
    return W_;
}

const xt::xtensor<double, 2>& Workspace::W() const {
    return W_;
}

xt::xtensor<double, 2>& Workspace::Rhs() {
    return rhs_;
}

const xt::xtensor<double, 2>& Workspace::Rhs() const {
    return rhs_;
}

// -------------------- cell-centered gradient getters --------------------

xt::xtensor<double, 2>& Workspace::GradRho() {
    return grad_rho_;
}

const xt::xtensor<double, 2>& Workspace::GradRho() const {
    return grad_rho_;
}

xt::xtensor<double, 2>& Workspace::GradU() {
    return grad_u_;
}

const xt::xtensor<double, 2>& Workspace::GradU() const {
    return grad_u_;
}

xt::xtensor<double, 2>& Workspace::GradV() {
    return grad_v_;
}

const xt::xtensor<double, 2>& Workspace::GradV() const {
    return grad_v_;
}

xt::xtensor<double, 2>& Workspace::GradP() {
    return grad_p_;
}

const xt::xtensor<double, 2>& Workspace::GradP() const {
    return grad_p_;
}

// -------------------- zero helpers --------------------

void Workspace::ZeroW() {
    W_.fill(0.0);
}

void Workspace::ZeroRhs() {
    rhs_.fill(0.0);
}

void Workspace::ZeroGradients() {
    grad_rho_.fill(0.0);
    grad_u_.fill(0.0);
    grad_v_.fill(0.0);
    grad_p_.fill(0.0);
}

void Workspace::ZeroAll() {
    ZeroW();
    ZeroRhs();
    ZeroGradients();
}

bool Workspace::IsAllocated() const {
    return
        W_.dimension() == 2 && rhs_.dimension() == 2 &&
        grad_rho_.dimension() == 2 && grad_u_.dimension() == 2 &&
        grad_v_.dimension() == 2 && grad_p_.dimension() == 2 &&
        W_.shape()[0] == n_cells_ &&
        W_.shape()[1] == k_nvar &&
        grad_rho_.shape()[0] == n_cells_ && grad_rho_.shape()[1] == 2;
}

std::size_t Workspace::GetCellCount() const {
    return n_cells_;
}

void Workspace::Allocate(const std::size_t n_cells) {
    n_cells_ = n_cells;

    const std::vector<std::size_t> var_shape = {n_cells_, k_nvar};
    const std::vector<std::size_t> grad_shape = {n_cells_, 2};

    W_ = xt::zeros<double>(var_shape);
    rhs_ = xt::zeros<double>(var_shape);

    grad_rho_ = xt::zeros<double>(grad_shape);
    grad_u_ = xt::zeros<double>(grad_shape);
    grad_v_ = xt::zeros<double>(grad_shape);
    grad_p_ = xt::zeros<double>(grad_shape);
}
