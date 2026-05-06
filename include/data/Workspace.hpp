#ifndef WORKSPACE_HPP
#define WORKSPACE_HPP

#include <cstddef>
#include <vector>

#include <xtensor.hpp>

class Mesh;

/**
 * @brief Reusable scratch buffers for one unstructured single-process mesh in 2D.
 *
 * Cell-centered fields:
 *   W(cell, var)          : auxiliary primitive-like cell data (rho, u, v, P)
 *   rhs(cell, var)        : conservative RHS
 *
 * Cell-centered gradient fields (for Kolgan reconstruction, shape [n_cells, 2]):
 *   grad_rho(cell, dir)
 *   grad_u(cell, dir)
 *   grad_v(cell, dir)
 *   grad_p(cell, dir)
 */
class Workspace final {
public:
    static constexpr std::size_t k_nvar = 4;

    // W(cell,var) indices
    static constexpr std::size_t k_rho = 0;
    static constexpr std::size_t k_u = 1;
    static constexpr std::size_t k_v = 2;
    static constexpr std::size_t k_p = 3;

    Workspace() = default;

    /**
     * @brief Resize all buffers to match mesh cell and face counts.
     */
    void ResizeFrom(const Mesh& mesh);

    // -------------------- cell-centered arrays --------------------

    [[nodiscard]] xt::xtensor<double, 2>& W();
    [[nodiscard]] const xt::xtensor<double, 2>& W() const;

    [[nodiscard]] xt::xtensor<double, 2>& Rhs();
    [[nodiscard]] const xt::xtensor<double, 2>& Rhs() const;

    // -------------------- cell-centered gradient arrays --------------------

    [[nodiscard]] xt::xtensor<double, 2>& GradRho();
    [[nodiscard]] const xt::xtensor<double, 2>& GradRho() const;

    [[nodiscard]] xt::xtensor<double, 2>& GradU();
    [[nodiscard]] const xt::xtensor<double, 2>& GradU() const;

    [[nodiscard]] xt::xtensor<double, 2>& GradV();
    [[nodiscard]] const xt::xtensor<double, 2>& GradV() const;

    [[nodiscard]] xt::xtensor<double, 2>& GradP();
    [[nodiscard]] const xt::xtensor<double, 2>& GradP() const;

    // -------------------- zero helpers --------------------

    void ZeroW();
    void ZeroRhs();
    void ZeroGradients();
    void ZeroAll();

    [[nodiscard]] bool IsAllocated() const;

    [[nodiscard]] std::size_t GetCellCount() const;

private:
    std::size_t n_cells_ = 0;

    // cell-centered
    xt::xtensor<double, 2> W_;
    xt::xtensor<double, 2> rhs_;

    // cell-centered gradients (n_cells, 2)
    xt::xtensor<double, 2> grad_rho_;
    xt::xtensor<double, 2> grad_u_;
    xt::xtensor<double, 2> grad_v_;
    xt::xtensor<double, 2> grad_p_;

    void Allocate(std::size_t n_cells);
};

#endif  // WORKSPACE_HPP
