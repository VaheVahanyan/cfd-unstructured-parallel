#include "config/InitialConditionInitializer.hpp"

#include <algorithm>
#include <stdexcept>

#include "data/DataLayer.hpp"
#include "data/Variables.hpp"
#include "geometry/Cell.hpp"
#include "geometry/Mesh.hpp"

InitialConditionInitializer::InitialConditionInitializer(
    const Settings& settings,
    const InitialConditions& initial_conditions
) : settings_(settings),
    initial_conditions_(initial_conditions) {}

void InitialConditionInitializer::Apply(const Mesh& mesh, DataLayer& layer) const {
    if (!layer.IsAllocated()) {
        throw std::runtime_error("InitialConditionInitializer: DataLayer is not allocated");
    }

    if (layer.GetCellCount() != mesh.GetCellCount()) {
        throw std::runtime_error("InitialConditionInitializer: DataLayer cell count does not match mesh");
    }

    switch (initial_conditions_.type) {
    case InitialConditionType::StructuredRegions:
        ApplyStructuredRegions(mesh, layer);
        return;

    case InitialConditionType::Constant:
        ApplyConstant(mesh, layer);
        return;

    case InitialConditionType::RegionMarkers:
        throw std::runtime_error(
                                 "InitialConditionInitializer: region_markers initial condition is not implemented yet"
                                );
    }

    throw std::runtime_error("InitialConditionInitializer: unsupported initial condition type");
}

void InitialConditionInitializer::ApplyStructuredRegions(const Mesh& mesh, DataLayer& layer) const {
    if (!initial_conditions_.structured_regions.has_value()) {
        throw std::runtime_error(
                                 "InitialConditionInitializer: structured_regions data is missing"
                                );
    }

    const StructuredRegionInitialCondition& ic = *initial_conditions_.structured_regions;
    ValidateStructuredRegionShape(ic);

    auto& U = layer.U();

    for (std::size_t cell_id = 0; cell_id < mesh.GetCellCount(); ++cell_id) {
        const Cell& cell = mesh.GetCell(cell_id);

        const std::size_t ix = RegionIndex(cell.center_x, ic.interfaces_x);
        const std::size_t iy = (settings_.mesh.dim >= 2)
                                   ? RegionIndex(cell.center_y, ic.interfaces_y)
                                   : 0;

        PrimitiveCell primitive;
        primitive.rho = ic.rho.At(ix, iy);
        primitive.u = ic.u.At(ix, iy);
        primitive.v = (settings_.mesh.dim >= 2) ? ic.v.At(ix, iy) : 0.0;
        primitive.P = ic.p.At(ix, iy);

        const ConservativeCell conservative =
            ConservativeFromPrimitive(primitive, settings_.gamma);

        U(cell_id, DataLayer::k_rho) = conservative.rho;
        U(cell_id, DataLayer::k_rhoU) = conservative.rhoU;
        U(cell_id, DataLayer::k_rhoV) = conservative.rhoV;
        U(cell_id, DataLayer::k_E) = conservative.E;
    }
}

void InitialConditionInitializer::ApplyConstant(const Mesh& mesh, DataLayer& layer) const {
    (void)mesh;

    if (!initial_conditions_.constant.has_value()) {
        throw std::runtime_error("InitialConditionInitializer: constant initial condition is missing");
    }

    const ConstantInitialCondition& ic = *initial_conditions_.constant;

    PrimitiveCell primitive;
    primitive.rho = ic.rho;
    primitive.u = ic.u;
    primitive.v = (settings_.mesh.dim >= 2) ? ic.v : 0.0;
    primitive.P = ic.p;

    const ConservativeCell conservative =
        ConservativeFromPrimitive(primitive, settings_.gamma);

    auto& U = layer.U();

    for (std::size_t cell_id = 0; cell_id < layer.GetCellCount(); ++cell_id) {
        U(cell_id, DataLayer::k_rho) = conservative.rho;
        U(cell_id, DataLayer::k_rhoU) = conservative.rhoU;
        U(cell_id, DataLayer::k_rhoV) = conservative.rhoV;
        U(cell_id, DataLayer::k_E) = conservative.E;
    }
}

std::size_t InitialConditionInitializer::RegionIndex(const double coord,
                                                     const std::vector<double>& interfaces) {
    return static_cast<std::size_t>(
        std::upper_bound(interfaces.begin(), interfaces.end(), coord) - interfaces.begin()
    );
}

void InitialConditionInitializer::ValidateStructuredRegionShape(
    const StructuredRegionInitialCondition& ic
) {
    const std::size_t nx = ic.RegionCountX();
    const std::size_t ny = ic.RegionCountY();

    auto validate = [&](const Field2DValues& field, const char* name) {
        if (field.Nx() != nx) {
            throw std::runtime_error(
                                     std::string("InitialConditionInitializer: ") + name +
                                     " x-shape does not match interface count"
                                    );
        }
        if (field.Ny() != ny) {
            throw std::runtime_error(
                                     std::string("InitialConditionInitializer: ") + name +
                                     " y-shape does not match interface count"
                                    );
        }
    };

    validate(ic.rho, "rho");
    validate(ic.u, "u");
    validate(ic.v, "v");
    validate(ic.p, "p");
}
