#ifndef DATALAYER_HPP
#define DATALAYER_HPP

#include <cstddef>
#include <stdexcept>

#include <xtensor.hpp>

class Mesh;

/**
 * @brief Storage owner for cell-centered conservative solution fields.
 *
 * Storage layout:
 *   U(cell, var), shape = (n_cells, 4)
 */
class DataLayer final {
public:
    static constexpr std::size_t k_nvar = 4;
    static constexpr std::size_t k_rho = 0;
    static constexpr std::size_t k_rhoU = 1;
    static constexpr std::size_t k_rhoV = 2;
    static constexpr std::size_t k_E = 3;

    DataLayer() = default;
    explicit DataLayer(std::size_t n_cells);

    /**
     * @brief Resize storage to the given number of cells.
     * @details Reallocates only if shape changed.
     */
    void Resize(std::size_t n_cells);

    /**
     * @brief Resize storage to match mesh cell count.
     */
    void ResizeFrom(const Mesh& mesh);

    /** @brief Conservative state array U(cell,var). Shape (n_cells, 4). */
    [[nodiscard]] xt::xtensor<double, 2>& U();
    [[nodiscard]] const xt::xtensor<double, 2>& U() const;

    /** @brief Number of cells in storage. */
    [[nodiscard]] std::size_t GetCellCount() const;

    /** @brief Check whether storage is allocated. */
    [[nodiscard]] bool IsAllocated() const;

private:
    std::size_t n_cells_ = 0;

    xt::xtensor<double, 2> U_;

    void Allocate(std::size_t n_cells);
};

#endif  // DATALAYER_HPP
