/*
 *  volume_transmitter.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "exceptions.h"
#include "volume_transmitter.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "arraydatum.h"
#include "connector_base.h"
#include "spikecounter.h"

#include <numeric>

/* ----------------------------------------------------------------
 * Default constructor defining default parameters
 * ---------------------------------------------------------------- */

mynest::volume_transmitter::Parameters_::Parameters_()
  : deliver_interval_( 1 ) // in steps of mindelay
{
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
mynest::volume_transmitter::Parameters_::get( DictionaryDatum& d ) const
{
  def< nest::long_t >( d, "deliver_interval", deliver_interval_ );
}

void ::mynest::volume_transmitter::Parameters_::set( const DictionaryDatum& d )
{
  updateValue< nest::long_t >( d, "deliver_interval", deliver_interval_ );
}

/* ----------------------------------------------------------------
 * Default and copy constructor for volume transmitter
 * ---------------------------------------------------------------- */

mynest::volume_transmitter::volume_transmitter()
  : nest::Archiving_Node()
  , P_()
{
}

mynest::volume_transmitter::volume_transmitter( const volume_transmitter& n )
  : nest::Archiving_Node( n )
  , P_( n.P_ )
{
}

void
mynest::volume_transmitter::init_state_( const nest::Node& )
{
}

void
mynest::volume_transmitter::init_buffers_()
{
  B_.neuromodulatory_spikes_.clear();
  B_.spikecounter_.clear();
  //B_.spikecounter_.push_back(nest::spikecounter( 0.0, 0.0 ) ); // insert pseudo last dopa spike at t = 0.0
  nest::Archiving_Node::clear_history();
}

void
mynest::volume_transmitter::calibrate()
{
  // +1 as pseudo dopa spike at t_trig is inserted after trigger_update_weight
  B_.spikecounter_.reserve( nest::Scheduler::get_min_delay() * P_.deliver_interval_ + 1  );
}

void
mynest::volume_transmitter::update( const nest::Time& t_orig, const nest::long_t from, const nest::long_t to )
{

  //std::cout << "Updating volume transmitter state at time " << t_orig.get_ms() << " from " << from << " to " << to << std::endl;

  // spikes that arrive in this time slice are stored in spikecounter_
  nest::double_t t_spike;
  nest::double_t multiplicity;
  for ( nest::long_t lag = from; lag < to; ++lag )
  {
    multiplicity = B_.neuromodulatory_spikes_.get_value( lag );
    if ( multiplicity > 0 )
    {
      t_spike =
        nest::Time(
          nest::Time::step( network()->get_slice_origin().get_steps()
            + lag + 1 ) ).get_ms();
      B_.spikecounter_.push_back( nest::spikecounter( t_spike, multiplicity ) );
      //std::cout << "Adding spike at time " << t_spike << " with multiplicity " << multiplicity << std::endl;
    }
  }

  // all spikes stored in spikecounter_ are delivered to the target synapses
  if ( ( network()->get_slice_origin().get_steps() + to )
      % ( P_.deliver_interval_ * nest::Scheduler::get_min_delay() )
    == 0 )
  {
    nest::double_t t_trig =
      nest::Time(
        nest::Time::step( network()->get_slice_origin().get_steps()
          + to ) ).get_ms();

    if ( !B_.spikecounter_.empty() )
      network()->trigger_update_weight(
        get_gid(), B_.spikecounter_, t_trig );

    // clear spikecounter
    B_.spikecounter_.clear();

    // as with trigger_update_weight dopamine trace has been updated to t_trig,
    // insert pseudo last dopa spike at t_trig
    //B_.spikecounter_.push_back( spikecounter( t_trig, 0.0 ) );
  }
}

void
mynest::volume_transmitter::handle( nest::SpikeEvent& e )
{
  B_.neuromodulatory_spikes_.add_value(
    e.get_rel_delivery_steps( network()->get_slice_origin() ),
    static_cast< nest::double_t >( e.get_multiplicity() ) );
  //std::cout << "Handle spike event at time " << e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin() ) << " with multiplicity " << e.get_multiplicity() << std::endl;
}
