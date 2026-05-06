#ifndef VARIABLES_HPP
#define VARIABLES_HPP

#include <cstddef>

/**
 * @file Variables.hpp
 * @brief Minimal variable utilities for 2D Euler equations (ideal gas).
 *
 * Conventions:
 *  - Conservative state: (rho, rhoU, rhoV, E)
 *  - Primitive state:    (rho, u, v, P)
 */

namespace var {
    static constexpr std::size_t rho = 0;
    static constexpr std::size_t u1 = 1;
    static constexpr std::size_t u2 = 2;
    static constexpr std::size_t p_or_E = 3;

    static constexpr std::size_t nvar = 4;
} // namespace var

/**
 * @brief Primitive state in one 2D cell.
 */
struct PrimitiveCell final {
    double rho = 0.0;
    double u = 0.0;
    double v = 0.0;
    double P = 0.0;
};

/**
 * @brief Conservative state in one 2D cell.
 * @details Order matches U = (rho, rhoU, rhoV, E).
 */
struct ConservativeCell final {
    double rho = 0.0;
    double rhoU = 0.0;
    double rhoV = 0.0;
    double E = 0.0;

    ConservativeCell& operator+=(const ConservativeCell& other);
    ConservativeCell& operator-=(const ConservativeCell& other);
};

[[nodiscard]] ConservativeCell operator+(ConservativeCell lhs, const ConservativeCell& rhs);
[[nodiscard]] ConservativeCell operator-(ConservativeCell lhs, const ConservativeCell& rhs);
[[nodiscard]] ConservativeCell operator*(double scalar, ConservativeCell value);
[[nodiscard]] ConservativeCell operator*(ConservativeCell value, double scalar);

/**
 * @brief Euler flux vector (4 components) in Cartesian form.
 */
struct FluxCell final {
    double mass = 0.0;
    double mom_x = 0.0;
    double mom_y = 0.0;
    double energy = 0.0;
};

/**
 * @brief One unit face normal in 2D.
 */
struct FaceNormal final {
    double x = 0.0;
    double y = 0.0;
};

[[nodiscard]] ConservativeCell ConservativeFromPrimitive(const PrimitiveCell& w, double gamma);

[[nodiscard]] PrimitiveCell PrimitiveFromConservativeCell(const ConservativeCell& U,
                                                          double gamma,
                                                          double rho_floor = 1e-14,
                                                          double p_floor = 1e-14);

void PrimitiveToConservative(const PrimitiveCell& w,
                             double gamma,
                             double& rho,
                             double& rhoU,
                             double& rhoV,
                             double& E);

[[nodiscard]] double SoundSpeed(const PrimitiveCell& w,
                                double gamma,
                                double p_floor = 1e-14,
                                double rho_floor = 1e-14);

[[nodiscard]] double KineticEnergyDensity(const PrimitiveCell& w);

[[nodiscard]] double TotalEnergyDensity(const PrimitiveCell& w, double gamma);

[[nodiscard]] bool IsUnitNormal(const FaceNormal& normal, double tolerance = 1e-10);

[[nodiscard]] double NormalVelocity(const PrimitiveCell& w, const FaceNormal& normal);

[[nodiscard]] ConservativeCell PhysicalFlux(const PrimitiveCell& w,
                                            double gamma,
                                            const FaceNormal& normal);

[[nodiscard]] FluxCell EulerFlux(const PrimitiveCell& w,
                                 double gamma,
                                 const FaceNormal& normal);

#endif  // VARIABLES_HPP
