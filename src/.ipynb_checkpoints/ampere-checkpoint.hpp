#ifndef HYBRIDIR_AMPERE_HPP
#define HYBRIDIR_AMPERE_HPP

#include "vecfield.hpp"

#include <cstddef>
#include <iostream>

template<std::size_t dimension>
class Ampere
{
public:
    Ampere(std::shared_ptr<GridLayout<dimension>> grid)
        : m_grid{grid}
    {
        if (!m_grid)
            throw std::runtime_error("GridLayout is null");
    }

    void operator()(VecField<dimension> const& B, VecField<dimension>& J)
    {
        auto const dx = m_grid->cell_size(Direction::X);

        if constexpr (dimension == 1)
        {
            // TODO your code here
            // Jy, Jz are primal
            for (auto ix = m_grid->primal_dom_start(Direction::X);
                 ix + 1 <= m_grid->primal_dom_end(Direction::X); ++ix)
            {
                auto const& Bx = B.x;
                auto const& By = B.y;
                auto const& Bz = B.z;

                auto& Jx = J.x;
                auto& Jy = J.y;
                auto& Jz = J.z;

                Jx(ix) = 0.0;
                Jy(ix) = - (Bz(ix + 1)- Bz(ix - 1))/(2*dx);
                Jz(ix) = (By(ix + 1)- By(ix - 1))/(2*dx);
            }
        }
        else
            throw std::runtime_error("Ampere not implemented for this dimension");
    }

private:
    std::shared_ptr<GridLayout<dimension>> m_grid;
};

#endif // HYBRIDIR_AMPERE_HPP
