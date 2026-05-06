#ifndef MORTONORDER_HPP
#define MORTONORDER_HPP

#include <cstddef>
#include <vector>

class Mesh;

/**
 * @class MortonOrder
 * @brief Reorders mesh cells by 2D Morton Z-order for improved cache locality.
 *
 * The reorder keeps owned cells before ghost cells. Within each group, cells are
 * sorted by Morton code computed from cell centers. Cell local ids, face owner
 * ids, face neighbor ids, and cell connectivity are updated consistently.
 */
class MortonOrder final {
public:
    /**
     * @brief Apply 2D Morton cell ordering in-place.
     *
     * @param mesh Mesh to reorder.
     */
    static void ApplyToCells(Mesh& mesh);

private:
    struct Entry final {
        std::size_t old_id = 0;
        std::size_t new_id = 0;
        unsigned long long code = 0;
    };

    static std::vector<Entry> BuildEntries(const Mesh& mesh,
                                           std::size_t begin,
                                           std::size_t end);

    static unsigned long long Encode2D(unsigned int x, unsigned int y);
    static unsigned long long Part1By1(unsigned int value);

    static void ReorderRange(Mesh& mesh,
                             std::size_t begin,
                             std::size_t end,
                             std::vector<std::size_t>& old_to_new);
};

#endif
