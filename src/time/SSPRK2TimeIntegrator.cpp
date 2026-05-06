#include "time/SSPRK2TimeIntegrator.hpp"

#include "data/DataLayer.hpp"
#include "data/Workspace.hpp"
#include "geometry/Mesh.hpp"
#include "parallel/HaloExchange.hpp"
#include "solver/PositivityLimiter.hpp"
#include "spatial/SpatialOperator.hpp"

void SSPRK2TimeIntegrator::Advance(DataLayer& layer,
                                   const Mesh& mesh,
                                   Workspace& workspace,
                                   const double dt,
                                   const double gamma,
                                   const SpatialOperator& op,
                                   const StateSynchronizer* halo_exchange) const {
    if (dt <= 0.0) {
        return;
    }

    workspace.ResizeFrom(mesh);

    auto& U = layer.U();
    const xt::xtensor<double, 2> U0 = U;

    const std::size_t n_owned = mesh.GetOwnedCellCount();

    if (halo_exchange) {
        halo_exchange->Synchronize(layer);
    }

    op.ComputeRHS(layer, mesh, workspace, gamma, dt);


    const auto& rhs_stage1 = workspace.Rhs();
#pragma omp parallel for default(none) shared(U, U0, rhs_stage1, n_owned) firstprivate(dt)
    for (std::size_t cell_id = 0; cell_id < n_owned; ++cell_id) {
        U(cell_id, 0) = U0(cell_id, 0) + dt * rhs_stage1(cell_id, 0);
        U(cell_id, 1) = U0(cell_id, 1) + dt * rhs_stage1(cell_id, 1);
        U(cell_id, 2) = U0(cell_id, 2) + dt * rhs_stage1(cell_id, 2);
        U(cell_id, 3) = U0(cell_id, 3) + dt * rhs_stage1(cell_id, 3);
    }

    PositivityLimiter::Apply(layer, mesh, gamma, rho_min_, p_min_);

    if (halo_exchange) {
        halo_exchange->Synchronize(layer);
    }

    op.ComputeRHS(layer, mesh, workspace, gamma, dt);

    const auto& rhs_stage2 = workspace.Rhs();
#pragma omp parallel for default(none) shared(U, U0, rhs_stage2, n_owned) firstprivate(dt)
    for (std::size_t cell_id = 0; cell_id < n_owned; ++cell_id) {
        U(cell_id, 0) = 0.5 * U0(cell_id, 0) + 0.5 * (U(cell_id, 0) + dt * rhs_stage2(cell_id, 0));
        U(cell_id, 1) = 0.5 * U0(cell_id, 1) + 0.5 * (U(cell_id, 1) + dt * rhs_stage2(cell_id, 1));
        U(cell_id, 2) = 0.5 * U0(cell_id, 2) + 0.5 * (U(cell_id, 2) + dt * rhs_stage2(cell_id, 2));
        U(cell_id, 3) = 0.5 * U0(cell_id, 3) + 0.5 * (U(cell_id, 3) + dt * rhs_stage2(cell_id, 3));
    }

    PositivityLimiter::Apply(layer, mesh, gamma, rho_min_, p_min_);
}
