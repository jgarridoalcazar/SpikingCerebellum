/*
 *  Archiving_Node_cs.cpp
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
 * \file archiving_node_cs.cpp
 * Implementation of archiving_node to record and manage complex spike history
 * \author Jesus Garrido
 * \date october 2017
 */

#include "archiving_node_cs.h"
#include "dictutils.h"
#include <cmath>
#include <cstdlib>

namespace nest
{
	namespace names
  {

  }
}

namespace mynest {


  //member functions for Archiving_Node

Archiving_Node_CS::Archiving_Node_CS() :
    Node(),
		n_incoming_cs_(0),
    history_cs_()
		{
		}

Archiving_Node_CS::Archiving_Node_CS(const Archiving_Node_CS& n)
:Node(n),
n_incoming_cs_(n.n_incoming_cs_),
history_cs_()
  {}

void Archiving_Node_CS::register_stdp_connection_cs(double t_first_read){
  // Mark all entries in the deque, which we will not read in future as read by this input
  // input, so that we savely increment the incoming number of
  // connections afterwards without leaving spikes in the history.
  // For details see bug #218. MH 08-04-22
	for ( std::deque<histentry_cs>::iterator runner = history_cs_.begin();
    runner != history_cs_.end() && runner->t_ <= t_first_read;
    ++runner){
    (runner->access_counter_)++;
  }
  n_incoming_cs_++;
}

  void Archiving_Node_CS::get_cs_history(double t1, double t2,
  				   std::deque<histentry_cs>::iterator* start,
  				   std::deque<histentry_cs>::iterator* finish)
  {
    *finish = history_cs_.end();
    if (history_cs_.empty()){
      *start = *finish;
      return;
    } else {
      std::deque<mynest::histentry_cs>::iterator runner = history_cs_.begin();
      while ((runner != history_cs_.end()) && (runner->t_ <= t1)) ++runner;
      *start = runner;
      while ((runner != history_cs_.end()) && (runner->t_ <= t2)) {
        (runner->access_counter_)++;
        ++runner;
  	  }
      *finish = runner;
    }
  }

  void mynest::Archiving_Node_CS::set_cs_spiketime(nest::Time const & t_sp, double offset)
  {
    //Archiving_Node::set_spiketime(t_sp, offset);

    const double t_sp_ms = t_sp.get_ms() - offset;

    if (n_incoming_cs_){
      // prune all spikes from history which are no longer needed
      // except the penultimate one. we might still need it.
      while (history_cs_.size() > 1){
        if (history_cs_.front().access_counter_ >= n_incoming_cs_){
          history_cs_.pop_front();
        } else {
          break;
        }
      }  

      history_cs_.push_back( histentry_cs( t_sp_ms, 0) );
    }
  }


  void mynest::Archiving_Node_CS::get_status(DictionaryDatum & d) const
  {
	  //Archiving_Node::get_status(d);
    //def<double>(d, nest::names::t_spike, get_spiketime_ms());
  #ifdef DEBUG_ARCHIVER
    def<int>(d, nest::names::archiver_length, history_cs_.size());
  #endif
  }

  void mynest::Archiving_Node_CS::set_status(const DictionaryDatum & d)
  {
	  //Archiving_Node::set_status(d);
    // We need to preserve values in case invalid values are set
	  // check, if to clear spike history and K_minus
    bool clear = false;
    updateValue<bool>(d, nest::names::clear, clear);
    if ( clear )
    	clear_history();
  }

  void mynest::Archiving_Node_CS::clear_history()
  {
  	//Archiving_Node::clear_history();

  	history_cs_.clear();
  }

} // of namespace nest
