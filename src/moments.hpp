#ifndef HYBIRT_MOMENTS_HPP
#define HYBIRT_MOMENTS_HPP

#include "field.hpp"
#include "vecfield.hpp"
#include "population.hpp"

#include <vector>


template<std::size_t dimension>
void total_density(std::vector<Population<dimension>> const& populations, Field<dimension>& N)
{
    for (auto ix = 0; ix < N.data().size(); ++ix)
    {
        N(ix) = 0;
    }
    for (auto const& pop : populations)
    {
        for (auto ix = 0; ix < N.data().size(); ++ix)
        {
            N(ix) += pop.density()(ix);
        }
    }
}

template<std::size_t dimension>
void bulk_velocity(std::vector<Population<dimension>> const& populations, Field<dimension> const& N,
                   VecField<dimension>& V)
{
    for (auto ix = 0; ix < N.data().size(); ++ix)
    {
        V.x(ix) = 0;
        V.y(ix) = 0;
        V.z(ix) = 0;
    }
    for (auto& pop : populations)
    {
        for (auto ix = 0; ix < N.data().size(); ++ix)
        {
            V.x(ix) += pop.flux().x(ix);
            V.y(ix) += pop.flux().y(ix);
            V.z(ix) += pop.flux().z(ix);
        }
    }
    // TODO calculate bulk velocity by dividing by density N
    for (auto ix = 0; ix < N.data().size(); ++ix)
    {
        V.x(ix) /= N(ix);
        V.y(ix) /= N(ix);
        V.z(ix) /= N(ix);
    }
}

#endif
