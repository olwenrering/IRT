#ifndef HYBIRT_BOUNDARY_CONDITION_HPP
#define HYBIRT_BOUNDARY_CONDITION_HPP

#include "field.hpp"
#include "vecfield.hpp"
#include "pusher.hpp"

#include <iostream>
#include <memory>
#include <string>
#include <stdexcept>
#include <cmath>


template<std::size_t dimension>
class BoundaryCondition
{
public:
    BoundaryCondition(std::shared_ptr<GridLayout<dimension>> const& grid)
        : m_grid(grid)
    {
        if (!m_grid)
            throw std::runtime_error("GridLayout is null");
    }

    virtual void fill(Field<dimension>& field) = 0;

    void fill(VecField<dimension>& vecfield)
    {
        fill(vecfield.x);
        fill(vecfield.y);
        fill(vecfield.z);
    }

    virtual void particles(std::vector<Particle<dimension>>& particles) = 0;

protected:
    std::shared_ptr<GridLayout<dimension>> m_grid;
};


template<std::size_t dimension>
class PeriodicBoundaryCondition : public BoundaryCondition<dimension>
{
public:
    PeriodicBoundaryCondition(std::shared_ptr<GridLayout<dimension>> const& grid)
        : BoundaryCondition<dimension>(grid)
    {
    }

    void fill(Field<dimension>& field) override
    {
        if constexpr (dimension == 1)
        {
            // left side
            auto gsi = this->m_grid->ghost_start(field.quantity(), Direction::X);
            auto dsi = this->m_grid->dom_start(field.quantity(), Direction::X);
            auto dei = this->m_grid->dom_end(field.quantity(), Direction::X);
            auto gei = this->m_grid->ghost_end(field.quantity(), Direction::X);

            auto const nbr_nodes = this->m_grid->nbr_cells(Direction::X);

            if (field.quantity() == Quantity::N or field.quantity() == Quantity::Vx
                or field.quantity() == Quantity::Vy or field.quantity() == Quantity::Vz)
            {
                for (auto ix_left = gsi; ix_left <= dsi; ++ix_left)
                {
                    auto const ix_right = ix_left + nbr_nodes;
                    field(ix_left) += field(ix_right);
                }

                for (auto ix_right = gei; ix_right > dei; --ix_right)
                {
                    auto const ix_left = ix_right - nbr_nodes;
                    field(ix_right) += field(ix_left);
                }
                field(dei) = field(dsi);
            }
            else
            {
                for (auto ix_left = gsi; ix_left < dsi; ++ix_left)
                {
                    auto const ix_right = ix_left + nbr_nodes;
                    field(ix_left)      = field(ix_right);
                }

                // std::cout << "Filling right side\n";
                for (auto ix_right = gei; ix_right > dei; --ix_right)
                {
                    auto const ix_left = ix_right - nbr_nodes;
                    field(ix_right)    = field(ix_left);
                }
            }
        }
    }

    void particles(std::vector<Particle<dimension>>& particles) override
    {
        if constexpr (dimension == 1)
        {
            for (auto& particle : particles)
            {
                double cell
                    = std::floor(particle.position[0] / this->m_grid->cell_size(Direction::X))
                      + this->m_grid->dual_dom_start(Direction::X);
                auto cell_save     = cell;
                auto position_save = particle.position[0];

                // particles left the right border injected on left side
                if (cell > this->m_grid->dual_dom_end(Direction::X))
                {
                    particle.position[0] -= this->m_grid->dom_size(Direction::X);
                }
                // particles left the left border injected on right side
                else if (cell < this->m_grid->dual_dom_start(Direction::X))
                {
                    // Wrap around to the right side
                    particle.position[0] += this->m_grid->dom_size(Direction::X);
                }

                if (particle.position[0] < 0.0
                    or particle.position[0] >= this->m_grid->dom_size(Direction::X))
                {
                    std::cout << "Particle position out of bounds after periodic BC: "
                              << particle.position[0] << " cell: " << cell
                              << " cell_save: " << cell_save << " position_save: " << position_save
                              << " dom_size: " << this->m_grid->dom_size(Direction::X) << "\n";
                    throw std::runtime_error("Particle position out of bounds after periodic BC");
                }
            }
        }
    }
};



template<std::size_t dimension>
class BoundaryConditionFactory
{
public:
    static std::unique_ptr<BoundaryCondition<dimension>>
    create(std::string const& type, std::shared_ptr<GridLayout<dimension>> grid)
    {
        if (type == "periodic")
            return std::make_unique<PeriodicBoundaryCondition<dimension>>(grid);
        // Add more boundary condition types as needed
        throw std::runtime_error("Unknown boundary condition type: " + type);
    }
};



#endif
