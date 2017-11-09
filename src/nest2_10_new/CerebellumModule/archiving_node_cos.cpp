/*
 *  Archiving_Node_cos.cpp
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
 * \file archiving_node_cos.cpp
 * Implementation of archiving_node to record and manage teaching signal history
 * \author Jesus Garrido
 * \date october 2017
 */

#include "archiving_node_cos.h"
#include "dictutils.h"
#include <cmath>
#include <cstdlib>

#include "ExponentialTable.h"
#include "TrigonometricTable.h"

namespace nest
{
	namespace names
  {
    const Name tau_cos("tau_cos");
    const Name exponent("exponent");
  }
}

namespace mynest {


  //member functions for Archiving_Node

  Archiving_Node_Cos::Archiving_Node_Cos() :
      Archiving_Node(),
      n_incoming_cos_(0),
      tau_cos_(1.0),
      inv_tau_cos_(1.0),
      exponent_(1.0),
      cos2_(0.0),
      sin2_(0.0),
      cossin_(0.0),
  		last_cos_spike_(-1.0),
      history_cos_()
  		{
  		}

  Archiving_Node_Cos::Archiving_Node_Cos(const Archiving_Node_Cos& n)
  :Archiving_Node(n),
  n_incoming_cos_(n.n_incoming_cos_),
  tau_cos_(n.tau_cos_),
  inv_tau_cos_(n.inv_tau_cos_),
  exponent_(n.exponent_),
  cos2_(n.cos2_),
  sin2_(n.sin2_),
  cossin_(n.cossin_),
  last_cos_spike_(n.last_cos_spike_),
  history_cos_()
    {}

  void Archiving_Node_Cos::register_stdp_connection_cos(double t_first_read){
    // Mark all entries in the deque, which we will not read in future as read by this input
    // input, so that we savely increment the incoming number of
    // connections afterwards without leaving spikes in the history.
    // For details see bug #218. MH 08-04-22
  	for ( std::deque<histentry_cos>::iterator runner = history_cos_.begin();
      runner != history_cos_.end() && runner->t_ <= t_first_read;
      ++runner){
      (runner->access_counter_)++;
    }
    n_incoming_cos_++;
  }

  void Archiving_Node_Cos::get_cos_values( double t,
                                    double& cos2,
                                    double& sin2,
                                    double& cossin ){
    // case when the neuron has not yet spiked. Evolved the state at the last
    // fired spike until the current time
    if ( history_cos_.empty() ) {
      cos2 = this->cos2_;
      sin2 = this->sin2_;
      cossin = this->cossin_;
      return;
    }
    
    // case
    int i = history_cos_.size() - 1;
    while ( i >= 0 ){
      if ( t > history_cos_[ i ].t_ ){
        this->evolve_cos_values(t - history_cos_[ i ].t_,
                            history_cos_[ i ].cos2_,
                            history_cos_[ i ].sin2_,
                            history_cos_[ i ].cossin_,
                            cos2, sin2, cossin);

        //std::cout << "Evolving postsynaptic trace from " << history_cos_[i].t_ << " to " << t << 
        //". Init values: " << history_cos_[ i ].cos2_ << " " << history_cos_[ i ].sin2_ << " " << history_cos_[ i ].cossin_ <<
        //". Final values: " << cos2 << " " << sin2 << " " << cossin << " Parameters: " << this->inv_tau_cos_ << " " << this->exponent_ << std::endl;

        return;
      }
      i--;
    }

    // we only get here if t< time of all spikes in history)
    // return 0.0 for both K values
    cos2 = 0.0;
    sin2 = 0.0;
    cossin = 0.0;
    return;
  }

  void Archiving_Node_Cos::evolve_cos_values( double ElapsedTime,
                                  double oldcos2, double oldsin2, double oldcossin,
                                  double& cos2, double& sin2, double& cossin){

    float ElapsedRelative = this->exponent_*ElapsedTime*this->inv_tau_cos_;
    float expon = ExponentialTable::GetResult(-ElapsedRelative);

    float ElapsedRelativeTrigonometric=ElapsedTime*this->inv_tau_cos_*1.5708f;


    int LUTindex=TrigonometricTable::CalculateOffsetPosition(ElapsedRelativeTrigonometric);
    LUTindex = TrigonometricTable::CalculateValidPosition(0,LUTindex);

    float SinVar = TrigonometricTable::GetElement(LUTindex);
    float CosVar = TrigonometricTable::GetElement(LUTindex+1);
    
    float auxCos2=CosVar*CosVar;
    float auxSin2=SinVar*SinVar;
    float auxCosSin=CosVar*SinVar;
  
    cos2 = expon*(oldcos2 * auxCos2 + oldsin2*auxSin2-2*oldcossin*auxCosSin);
    sin2 = expon*(oldsin2 * auxCos2 + oldcos2*auxSin2+2*oldcossin*auxCosSin);
    cossin = expon*(oldcossin *(auxCos2-auxSin2) + (oldcos2-oldsin2)*auxCosSin);

    return;
  }



  void Archiving_Node_Cos::get_cos_history(double t1, double t2,
  			   std::deque<histentry_cos>::iterator* start,
  			   std::deque<histentry_cos>::iterator* finish)
  {
    *finish = history_cos_.end();
    if (history_cos_.empty()){
      *start = *finish;
      return;
    } else {
      std::deque<mynest::histentry_cos>::iterator runner = history_cos_.begin();
      while ((runner != history_cos_.end()) && (runner->t_ <= t1)) ++runner;
      *start = runner;
      while ((runner != history_cos_.end()) && (runner->t_ <= t2)) {
        (runner->access_counter_)++;
        ++runner;
  	  }
      *finish = runner;
    }
  }

  void mynest::Archiving_Node_Cos::set_cos_spiketime(nest::Time const & t_sp, double offset)
  {
    const double t_sp_ms = t_sp.get_ms() - offset;

    if (n_incoming_cos_){
      // prune all spikes from history which are no longer needed
      // except the penultimate one. we might still need it.
      while (history_cos_.size() > 1){
        if (history_cos_.front().access_counter_ >= n_incoming_cos_){
          history_cos_.pop_front();
        } else {
          break;
        }
      }  

      this->evolve_cos_values( t_sp_ms - this->last_cos_spike_,
                                  this->cos2_, this->sin2_, this->cossin_,
                                  this->cos2_, this->sin2_, this->cossin_);

      this->cos2_ += 1.0;
      last_cos_spike_ = t_sp_ms;
      history_cos_.push_back( histentry_cos( last_cos_spike_, this->cos2_, this->sin2_, this->cossin_, 0 ) );
    } else {
      last_cos_spike_ = t_sp_ms;
    }
  }


  void mynest::Archiving_Node_Cos::get_status(DictionaryDatum & d) const
  {
	  Archiving_Node::get_status(d);

    def< double >( d, nest::names::tau_cos, this->tau_cos_ );
    def< double >( d, nest::names::exponent, this->exponent_ );
  #ifdef DEBUG_ARCHIVER
    def<int>(d, nest::names::archiver_length, history_cos_.size());
  #endif
  }

  void mynest::Archiving_Node_Cos::set_status(const DictionaryDatum & d)
  {
	  Archiving_Node::set_status(d);

    updateValue< double >( d, nest::names::tau_cos, this->tau_cos_ );
    updateValue< double >( d, nest::names::exponent, this->exponent_ );

    if ( this->tau_cos_ <= 0.0 )
    {
      throw nest::BadProperty( "All time constants must be strictly positive." );
    }

    this->inv_tau_cos_ = 1./this->tau_cos_;

    // We need to preserve values in case invalid values are set
	  // check, if to clear spike history and K_minus
    bool clear = false;
    updateValue<bool>(d, nest::names::clear, clear);
    if ( clear )
    	clear_history();
  }

  void mynest::Archiving_Node_Cos::clear_history()
  {
  	Archiving_Node::clear_history();

  	history_cos_.clear();
  }

} // of namespace nest
