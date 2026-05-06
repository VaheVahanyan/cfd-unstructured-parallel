#include "time/ForwardEulerTimeIntegrator.hpp"

#include "data/DataLayer.hpp"
#include "data/Workspace.hpp"
#include "geometry/Mesh.hpp"
#include "parallel/HaloExchange.hpp"
#include "solver/PositivityLimiter.hpp"
#include "spatial/SpatialOperator.hpp"

void ForwardEulerTimeIntegrator::Advance(DataLayer& layer,
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
    if (halo_exchange) {
        halo_exchange->Synchronize(layer);
    }
    op.ComputeRHS(layer, mesh, workspace, gamma, dt);

    auto& U = layer.U();
    const auto& rhs = workspace.Rhs();

    const std::size_t n_owned = mesh.GetOwnedCellCount();

#pragma omp parallel for default(none) shared(U, rhs, n_owned) firstprivate(dt)
    for (std::size_t i = 0; i < n_owned; ++i) {
        U(i, 0) += dt * rhs(i, 0);
        U(i, 1) += dt * rhs(i, 1);
        U(i, 2) += dt * rhs(i, 2);
        U(i, 3) += dt * rhs(i, 3);
    }

    PositivityLimiter::Apply(layer, mesh, gamma, rho_min_, p_min_);
}
