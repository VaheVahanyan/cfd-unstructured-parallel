#ifndef INITIALCONDITIONS_HPP
#define INITIALCONDITIONS_HPP

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "config/Settings.hpp"

/**
 * @brief Supported initial-condition representations.
 */
enum class InitialConditionType {
    StructuredRegions,
    Constant,
    RegionMarkers
};

/**
 * @brief Structured region values stored in normalized 3D form.
 *
 * Layout:
 * - 1D: [nx][1]
 * - 2D: [nx][ny]
 */
struct Field2DValues final {
    std::vector<std::vector<double>> values;

    [[nodiscard]] std::size_t Nx() const {
        return values.size();
    }

    [[nodiscard]] std::size_t Ny() const {
        return values.empty() ? 0 : values[0].size();
    }

    [[nodiscard]] double At(std::size_t ix, std::size_t iy) const {
        return values.at(ix).at(iy);
    }

    [[nodiscard]] bool Empty() const {
        return values.empty();
    }
};

/**
 * @brief Region-based structured initial-condition description.
 *
 * Region counts:
 * - x regions = interfaces_x.size() + 1
 * - y regions = interfaces_y.size() + 1
 */
struct StructuredRegionInitialCondition final {
    std::vector<double> interfaces_x;
    std::vector<double> interfaces_y;

    Field2DValues rho;
    Field2DValues u;
    Field2DValues v;
    Field2DValues p;

    [[nodiscard]] std::size_t RegionCountX() const {
        return interfaces_x.size() + 1;
    }

    [[nodiscard]] std::size_t RegionCountY() const {
        return interfaces_y.size() + 1;
    }

};

/**
 * @brief Constant primitive initial state over the whole mesh.
 */
struct ConstantInitialCondition final {
    double rho = 0.0;
    double u = 0.0;
    double v = 0.0;
    double p = 0.0;
};

/**
 * @brief Initial-condition container for one case.
 *
 * At the current stage only one representation is expected to be populated
 * according to the selected type.
 */
struct InitialConditions final {
    InitialConditionType type = InitialConditionType::StructuredRegions;

    std::optional<StructuredRegionInitialCondition> structured_regions;
    std::optional<ConstantInitialCondition> constant;

    /**
     * @brief Case-local runtime overrides.
     */
    CaseSettings overrides;
};

#endif  // INITIALCONDITIONS_HPP
