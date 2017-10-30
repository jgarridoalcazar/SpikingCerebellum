/*
 *  Archiving_Node_cs.h
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

/**
 * \file Archiving_Node_cs.h
 * Definition of Archiving_Node which is capable of
 * recording and managing a spike history.
 * \author Jesus Garrido
 * \date October 2017
 */

#ifndef ARCHIVING_NODE_CS_H
#define ARCHIVING_NODE_CS_H

#include "nest.h"
#include "archiving_node.h"
#include "dictdatum.h"
#include "nest_time.h"
#include "histentry_cs.h"
#include <deque>

#define DEBUG_ARCHIVER 1

namespace nest {

  namespace names
	{
    	// Neuron parameters
  }
}

namespace mynest {

/**
 * \class Archiving_Node
 * a node which archives spike history for the purposes of
 * timing dependent plasticity
 */
  class Archiving_Node_CS: public nest::Archiving_Node
{
   
 public:
  

  /**
   * \fn Archiving_Node()
   * Constructor.
   */
  Archiving_Node_CS();

  /**
   * \fn Archiving_Node_Sym()
   * Copy Constructor.
   */
  Archiving_Node_CS(const Archiving_Node_CS&);


  /**
   * \fn void get_CS_history(long_t t1, long_t t2, std::deque<Archiver::histentry>::iterator* start, std::deque<Archiver::histentry>::iterator* finish)
   * return the spike times (in steps) of spikes which occurred in the range (t1,t2].
   */
  void get_cs_history(double t1, double t2,
          std::deque<histentry_cs>::iterator* start,
    		  std::deque<histentry_cs>::iterator* finish);

    /**
     * Register a new incoming STDP connection.
     *
     * t_first_read: The newly registered synapse will read the history entries with t > t_first_read.
     */
    void register_stdp_connection_cs(double t_first_read);

    void get_status(DictionaryDatum & d) const;
    void set_status(const DictionaryDatum & d);

 protected:

  /**
   * \fn void set_spiketime(Time const & t_sp)
   * record spike history
   */
  void set_cs_spiketime(nest::Time const & t_sp, double offset=0);

  /**
   * \fn void clear_history()
   * clear spike history
   */
  void clear_history();


 private:

  // number of incoming connections from istdp connectors.
    // needed to determine, if every incoming connection has
    // read the spikehistory for a given point in time
    size_t n_incoming_cs_;

    // Accumulation variables
    
    // spiking history needed by stdp synapses
    std::deque<histentry_cs> history_cs_;

};
  
} // of namespace

#endif



