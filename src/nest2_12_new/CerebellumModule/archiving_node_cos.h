/*
 *  Archiving_Node_Cos.h
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
 * \file Archiving_Node_Cos.h
 * Definition of Archiving_Node which is capable of
 * recording and managing a spike history.
 * \author Jesus Garrido
 * \date October 2017
 */

#ifndef ARCHIVING_NODE_COS_H
#define ARCHIVING_NODE_COS_H

#include "nest.h"
#include "archiving_node.h"
#include "dictdatum.h"
#include "nest_time.h"
#include "histentry_cos.h"
#include <deque>

#define DEBUG_ARCHIVER 1

namespace nest {

  namespace names
	{
    // Neuron parameters
    extern const Name tau_cos;  
    extern const Name exponent;  
  }
}

namespace mynest {

/**
 * \class Archiving_Node
 * a node which archives spike history for the purposes of
 * cosine timing dependent plasticity
 */
  class Archiving_Node_Cos: public nest::Archiving_Node
{
   
 public:
  

  /**
   * \fn Archiving_Node()
   * Constructor.
   */
  Archiving_Node_Cos();

  /**
   * \fn Archiving_Node_Sym()
   * Copy Constructor.
   */
  Archiving_Node_Cos(const Archiving_Node_Cos&);


  /**
   * \fn void get_cos_history(long_t t1, long_t t2, std::deque<Archiver::histentry>::iterator* start, std::deque<Archiver::histentry>::iterator* finish)
   * return the spike times (in steps) of spikes which occurred in the range (t1,t2].
   */
  void get_cos_history(double t1, double t2,
          std::deque<histentry_cos>::iterator* start,
    		  std::deque<histentry_cos>::iterator* finish);


  /**
   * \fn void get_cos_value(double t, double cos2, double sin2, double cossin)
   * return the trace values at the specified time.
   */
  void get_cos_values(double t, double& cos2, double& sin2, double& cossin);

  /**
   * Register a new incoming STDP connection.
   *
   * t_first_read: The newly registered synapse will read the history entries with t > t_first_read.
   */
  void register_stdp_connection_cos(double t_first_read);

  void get_status(DictionaryDatum & d) const;
  void set_status(const DictionaryDatum & d);

 protected:

  /**
   * \fn void set_spiketime(Time const & t_sp)
   * record spike history
   */
  void set_cos_spiketime(nest::Time const & t_sp, double offset=0);

  /**
   * \fn void clear_history()
   * clear spike history
   */
  void clear_history();


 private:

  // number of incoming connections from istdp connectors.
    // needed to determine, if every incoming connection has
    // read the spikehistory for a given point in time
    size_t n_incoming_cos_;

    double tau_cos_;

    // Inverse of the learning rule tau
    double inv_tau_cos_;

    // Exponent of the cos function
    double exponent_;

    // Cos^2 accumulation variable
    double cos2_;

    // Sin^2 accumulation variable
    double sin2_;

    // Cos*Sin accumulation variable
    double cossin_;

    double last_cos_spike_;

    // spiking history needed by stdp synapses
    std::deque<histentry_cos> history_cos_;

    void evolve_cos_values( double ElapsedTime, 
                          double oldcos2, double oldsin2, double oldcossin,
                          double& cos2, double& sin2, double& cossin);

};
  
} // of namespace

#endif



