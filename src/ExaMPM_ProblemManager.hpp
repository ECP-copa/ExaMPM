/****************************************************************************
 * Copyright (c) 2018-2020 by the ExaMPM authors                            *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the ExaMPM library. ExaMPM is distributed under a   *
 * BSD 3-clause license. For the licensing terms see the LICENSE file in    *
 * the top-level directory.                                                 *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#ifndef EXAMPM_PROBLEMMANAGER_HPP
#define EXAMPM_PROBLEMMANAGER_HPP

#include <ExaMPM_ParticleInit.hpp>
#include <ExaMPM_Mesh.hpp>

#include <Cabana_Core.hpp>

#include <Cajita.hpp>

#include <memory>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
// Field locations
namespace Location
{
struct Cell {};
struct Node {};
struct Particle {};
} // end namespace Location

//---------------------------------------------------------------------------//
// Fields.
namespace Field
{
struct Velocity {};
struct Affine {};
struct Position {};
struct Mass {};
struct Volume {};
struct J {};
struct Momentum {};
struct Force {};
struct PositionCorrection {};
struct Density {};
struct Mark {};
} // end namespace Field.

//---------------------------------------------------------------------------//
template<class MemorySpace>
class ProblemManager
{
  public:

    using memory_space = MemorySpace;
    using execution_space = typename memory_space::execution_space;

    using particle_members = Cabana::MemberTypes<double[3][3],
                                                 double[3],
                                                 double[3],
                                                 double,
                                                 double,
                                                 double>;
    using particle_list = Cabana::AoSoA<particle_members,MemorySpace>;
    using particle_type = typename particle_list::tuple_type;

    using node_array = Cajita::Array<double,
                                     Cajita::Node,
                                     Cajita::UniformMesh<double>,
                                     MemorySpace>;

    using cell_array = Cajita::Array<double,
                                     Cajita::Cell,
                                     Cajita::UniformMesh<double>,
                                     MemorySpace>;

    using halo = Cajita::Halo<MemorySpace>;

    using mesh_type = Mesh<MemorySpace>;

    template<class InitFunc, class ExecutionSpace>
    ProblemManager( const ExecutionSpace& exec_space,
                    const std::shared_ptr<mesh_type>& mesh,
                    const InitFunc& create_functor,
                    const int particles_per_cell,
                    const double bulk_modulus,
                    const double rho,
                    const double gamma,
                    const double kappa )
        : _mesh( mesh )
        , _bulk_modulus( bulk_modulus )
        , _rho( rho )
        , _gamma( gamma )
        , _kappa( kappa )
        , _particles( "particles" )
    {
        initializeParticles(
            exec_space, *(_mesh->localGrid()), particles_per_cell,
            create_functor, _particles );

        auto node_vector_layout =
            Cajita::createArrayLayout( _mesh->localGrid(), 3, Cajita::Node() );
        auto node_scalar_layout =
            Cajita::createArrayLayout( _mesh->localGrid(), 1, Cajita::Node() );
        auto cell_scalar_layout =
            Cajita::createArrayLayout( _mesh->localGrid(), 1, Cajita::Cell() );

        _momentum =
            Cajita::createArray<double,MemorySpace>( "momentum", node_vector_layout );
        _mass =
            Cajita::createArray<double,MemorySpace>( "mass", node_scalar_layout );
        _force =
            Cajita::createArray<double,MemorySpace>( "force", node_vector_layout );
        _velocity =
            Cajita::createArray<double,MemorySpace>( "velocity", node_vector_layout );
        _position_correction =
            Cajita::createArray<double,MemorySpace>( "position_correction", node_vector_layout );
        _density =
            Cajita::createArray<double,MemorySpace>( "density", cell_scalar_layout );

        _mark =
            Cajita::createArray<double,MemorySpace>( "mark", cell_scalar_layout );

        _node_vector_halo =
            Cajita::createHalo<double,MemorySpace>( *node_vector_layout,
                                                    Cajita::FullHaloPattern() );
        _node_scalar_halo =
            Cajita::createHalo<double,MemorySpace>( *node_scalar_layout,
                                                    Cajita::FullHaloPattern() );
        _cell_scalar_halo =
            Cajita::createHalo<double,MemorySpace>( *cell_scalar_layout,
                                                    Cajita::FullHaloPattern() );
    }

    std::size_t numParticle() const
    {
        return _particles.size();
    }

    const std::shared_ptr<mesh_type>& mesh() const
    {
        return _mesh;
    }

    double bulkModulus() const
    {
        return _bulk_modulus;
    }

    double density() const
    {
        return _rho;
    }

    double gamma() const
    {
        return _gamma;
    }

    double kappa() const
    {
        return _kappa;
    }

    typename particle_list::template member_slice_type<0>
    get( Location::Particle, Field::Affine ) const
    {
        return Cabana::slice<0>( _particles, "affine" );
    }

    typename particle_list::template member_slice_type<1>
    get( Location::Particle, Field::Velocity ) const
    {
        return Cabana::slice<1>( _particles, "velocity" );
    }

    typename particle_list::template member_slice_type<2>
    get( Location::Particle, Field::Position ) const
    {
        return Cabana::slice<2>( _particles, "position" );
    }

    typename particle_list::template member_slice_type<3>
    get( Location::Particle, Field::Mass) const
    {
        return Cabana::slice<3>( _particles, "mass" );
    }

    typename particle_list::template member_slice_type<4>
    get( Location::Particle, Field::Volume ) const
    {
        return Cabana::slice<4>( _particles, "volume" );
    }

    typename particle_list::template member_slice_type<5>
    get( Location::Particle, Field::J ) const
    {
        return Cabana::slice<5>( _particles, "J" );
    }

    typename node_array::view_type
    get( Location::Node, Field::Momentum ) const
    {
        return _momentum->view();
    }

    typename node_array::view_type
    get( Location::Node, Field::Mass ) const
    {
        return _mass->view();
    }

    typename node_array::view_type
    get( Location::Node, Field::Force ) const
    {
        return _force->view();
    }

    typename node_array::view_type
    get( Location::Node, Field::Velocity ) const
    {
        return _velocity->view();
    }

    typename node_array::view_type
    get( Location::Node, Field::PositionCorrection ) const
    {
        return _position_correction->view();
    }

    typename cell_array::view_type
    get( Location::Cell, Field::Density ) const
    {
        return _density->view();
    }

    typename cell_array::view_type
    get( Location::Cell, Field::Mark ) const
    {
        return _mark->view();
    }

    void scatter( Location::Node, Field::Momentum ) const
    {
        _node_vector_halo->scatter( execution_space(),
                                    Cajita::ScatterReduce::Sum(),
                                    *_momentum );
    }

    void scatter( Location::Node, Field::Mass ) const
    {
        _node_scalar_halo->scatter( execution_space(),
                                    Cajita::ScatterReduce::Sum(),
                                    *_mass );
    }

    void scatter( Location::Node, Field::Force ) const
    {
        _node_vector_halo->scatter( execution_space(),
                                    Cajita::ScatterReduce::Sum(),
                                    *_force );
    }

    void scatter( Location::Node, Field::PositionCorrection ) const
    {
        _node_vector_halo->scatter( execution_space(),
                                    Cajita::ScatterReduce::Sum(),
                                    *_position_correction );
    }

    void scatter( Location::Cell, Field::Density ) const
    {
        _cell_scalar_halo->scatter( execution_space(),
                                    Cajita::ScatterReduce::Sum(),
                                    *_density );
    }

    void scatter( Location::Cell, Field::Mark ) const
    {
        _cell_scalar_halo->scatter( execution_space(),
                                    Cajita::ScatterReduce::Sum(),
                                    *_mark );
    }

    void gather( Location::Node, Field::Velocity ) const
    {
        _node_vector_halo->gather( execution_space(), *_velocity );
    }

    void gather( Location::Node, Field::PositionCorrection ) const
    {
        _node_vector_halo->gather( execution_space(), *_position_correction );
    }

    void communicateParticles( const int minimum_halo_width )
    {
        auto positions = get( Location::Particle(), Field::Position()) ;
        Cajita::particleGridMigrate(
            *(_mesh->localGrid()),
            positions,
            _particles,
            minimum_halo_width );
    }

  private:

    std::shared_ptr<mesh_type> _mesh;
    double _bulk_modulus;
    double _rho;
    double _gamma;
    double _kappa;
    particle_list _particles;
    std::shared_ptr<node_array> _momentum;
    std::shared_ptr<node_array> _mass;
    std::shared_ptr<node_array> _force;
    std::shared_ptr<node_array> _velocity;
    std::shared_ptr<node_array> _position_correction;
    std::shared_ptr<cell_array> _density;
    std::shared_ptr<cell_array> _mark;
    std::shared_ptr<halo> _node_vector_halo;
    std::shared_ptr<halo> _node_scalar_halo;
    std::shared_ptr<halo> _cell_scalar_halo;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_PROBLEMMANAGER_HPP
