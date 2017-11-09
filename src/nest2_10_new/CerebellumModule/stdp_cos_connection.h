/*
 *  stdp_cos_connection.h
 */

#ifndef STDP_COS_CONNECTION_H
#define STDP_COS_CONNECTION_H

/* BeginDocumentation

   Name: stdp_cos_synapse - Synapse type for DCN-like spike-timing dependent
   plasticity.

   Description:
   stdp_cos_synapse is a connection to create synapses with
   deep-cerebellar-nuclei-like spike-timing dependent plasticity (implemented
   based on [1]). Every time the complex spike generates an spike,
   this synapsis will produce LTD according to a temporal kernel based on
   f(t) = exp(-t/tau_c)*sin(-t/tau_c)^exponent, where t represents the time
   since the last presynaptic spike in the synapses.

   Parameters:
      A_plus    double - Amplitude of weight change for facilitation
      A_minus   double - Amplitude of weight change for depression
      Wmin      double - Minimal synaptic weight
      Wmax      double - Maximal synaptic weight
      exponent  unsigned int - Exponent of the sin function (integer between 1 and 20). The lower the exponent the wider the kernel function.
      tau_cos   double - Time constant of the learning rule (in ms)

  References:
   [1] Luque, N. R., Garrido, J. A., Naveros, F., Carrillo, R. R., D'Angelo, E., & Ros, E. (2016). 
   Distributed cerebellar motor learning: a spike-timing-dependent plasticity model. 
   Frontiers in computational neuroscience, 10. 


   Transmits: SpikeEvent

   Author: Jesus Garrido
   Remarks:
   - based on previous code for EDLUT simulator (https://github.com/EduardoRosLab/edlut)
   
   SeeAlso: iaf_cond_exp_cs
*/

#include "common_synapse_properties.h"

#include "connection.h"
#include "archiving_node_cos.h"

#include "ExponentialTable.h"
#include "TrigonometricTable.h"

#define A 1.0f/2.0f

namespace mynest
{

/**
 * Class representing an STDPCosConnection.
 */
template < typename targetidentifierT >
class STDPCosConnection : public nest::Connection< targetidentifierT >
{


public:
  typedef nest::CommonSynapseProperties CommonPropertiesType;
  typedef nest::Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  STDPCosConnection();

  /**
   * Copy constructor from a property object.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  STDPCosConnection( const STDPCosConnection& );

  // Explicitly declare all methods inherited from the dependent base ConnectionBase.
  // This avoids explicit name prefixes in all places these functions are used.
  // Since ConnectionBase depends on the template parameter, they are not automatically
  // found in the base class.
  using ConnectionBase::get_delay;
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;

  /**
   * Get all properties of this connection and put them into a dictionary.
   */
  void get_status( DictionaryDatum& d ) const;

  /**
   * Set properties of this connection from the values given in dictionary.
   */
  void set_status( const DictionaryDatum& d, nest::ConnectorModel& cm );

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   */
  void send( nest::Event& e, nest::thread t, double, const nest::CommonSynapseProperties& cp );

  class ConnTestDummyNode : public nest::ConnTestDummyNodeBase
  {
  public:
    // Ensure proper overriding of overloaded virtual functions.
    // Return values from functions are ignored.
    using ConnTestDummyNodeBase::handles_test_event;
    nest::port handles_test_event( nest::SpikeEvent&, nest::rport )
    {
      return nest::invalid_port_;
    }
  };

  /*
   * This function calls check_connection on the sender and checks if the receiver
   * accepts the event type and receptor type requested by the sender.
   * Node::check_connection() will either confirm the receiver port by returning
   * true or false if the connection should be ignored.
   *
   * \param s The source node
   * \param r The target node
   * \param receptor_type The ID of the requested receptor type
   * \param t_lastspike last spike produced by presynaptic neuron (in ms)
   */
  void
  check_connection( nest::Node& s,
    nest::Node& t,
    nest::rport receptor_type,
    double t_lastspike,
    const CommonPropertiesType& cp )
  {
    ConnTestDummyNode dummy_target;
    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );

    ((Archiving_Node_Cos *) (&t))->register_stdp_connection_cos( t_lastspike - get_delay() );
  }

  void
  set_weight( double w )
  {
    weight_ = w;
  }

private:
  // data members of each connection
  double weight_;

  double last_spike_weight_change_;
  
  // Cos^2 accumulation variable
  double cos2_;

  // Sin^2 accumulation variable
  double sin2_;

  // Cos*Sin accumulation variable
  double cossin_;

  double t_last_update_;

  // This vars could be common, but state_vars has to be resized according to Exponent.
  // Remaining vars depends on Exponent too.
  double exponent_;
  double inv_tau_;  

  // This vars could be also common, but this optimization might be performed as a future development.
  double A_plus_;
  double A_minus_;
  double Wmin_;
  double Wmax_;

  void evolve_cos_values( double ElapsedTime,
                          double oldcos2, double oldsin2, double oldcossin,
                          double& cos2, double& sin2, double& cossin);

  double check_weight_boundaries(double weight);
};

//
// Implementation of class STDPSinConnection.
//

template < typename targetidentifierT >
STDPCosConnection< targetidentifierT >::STDPCosConnection()
  : ConnectionBase(),
  weight_( 1.0 ),
  last_spike_weight_change_( 0.0 ),
  cos2_( 0.0 ),
  sin2_( 0.0 ),
  cossin_( 0.0 ),
  t_last_update_( 0.0 ),
  exponent_( 2 ),
  inv_tau_( 1.0 ),
  A_plus_( 1.0 ),
  A_minus_( 1.0 ),
  Wmin_( 0.0 ),
  Wmax_( 200.0 )
{
}

template < typename targetidentifierT >
STDPCosConnection< targetidentifierT >::STDPCosConnection( const STDPCosConnection& rhs )
  : ConnectionBase( rhs )
  , weight_( rhs.weight_ )
  , last_spike_weight_change_ ( rhs.last_spike_weight_change_ )
  , cos2_( rhs.cos2_ )
  , sin2_( rhs.sin2_ )
  , cossin_( rhs.cossin_ )
  , t_last_update_( rhs.t_last_update_ )
  , exponent_( rhs.exponent_ )
  , inv_tau_( rhs.inv_tau_ )
  , A_plus_( rhs.A_plus_ )
  , A_minus_( rhs.A_minus_ )
  , Wmin_( rhs.Wmin_ )
  , Wmax_( rhs.Wmax_ )
{
}

template < typename targetidentifierT >
void
STDPCosConnection< targetidentifierT >::get_status( DictionaryDatum& d ) const
{

  // base class properties, different for individual synapse
  ConnectionBase::get_status( d );
  def< double >( d, "A_plus", A_plus_ );
  def< double >( d, "A_minus", A_minus_ );
  def< double >( d, "Wmin", Wmin_ );
  def< double >( d, "Wmax", Wmax_ );
  def< double >( d, "tau_cos", 1./this->inv_tau_);
  def< double >( d, nest::names::weight, this->weight_ );

  def< double >( d, "exponent", this->exponent_ );
}

template < typename targetidentifierT >
void
STDPCosConnection< targetidentifierT >::set_status( const DictionaryDatum& d, nest::ConnectorModel& cm )
{
  // base class properties
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, nest::names::weight, weight_ );

  updateValue< double >( d, "A_plus", A_plus_ );
  updateValue< double >( d, "A_minus", A_minus_ );
  
  updateValue< double >( d, "Wmin", Wmin_ );
  updateValue< double >( d, "Wmax", Wmax_ );

  updateValue< double >( d, "exponent", exponent_ );
  double new_tau_cos;
  updateValue< double >( d, "tau_cos", new_tau_cos );
  this->inv_tau_ = 1./new_tau_cos;
}

template < typename targetidentifierT >
void STDPCosConnection< targetidentifierT >::evolve_cos_values( double ElapsedTime,
                                  double oldcos2,
                                  double oldsin2,
                                  double oldcossin,
                                  double& cos2,
                                  double& sin2,
                                  double& cossin){

    float ElapsedRelative = this->exponent_*ElapsedTime*this->inv_tau_;
    float expon = ExponentialTable::GetResult(-ElapsedRelative);

    float ElapsedRelativeTrigonometric=ElapsedTime*this->inv_tau_*1.5708f;


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

/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 * \param t_lastspike Time point of last spike emitted
 */
template < typename targetidentifierT >
inline void
STDPCosConnection< targetidentifierT >::send( nest::Event& e,
  nest::thread t,
  double,
  const nest::CommonSynapseProperties& cp )
{
  // t_lastspike_ = 0 initially

  //std::cout << "Synapse parameters: Aplus: " << this->A_plus_ << ". Aminus: " << this->A_minus_ << ". Wmin: " << this->Wmin_ << ". Wmax: " << this->Wmax_ << ". Expon: " << this->exponent_ << ". Tau: " << 1./this->inv_tau_ << ". Weight: " << this->weight_ << std::endl;
  
  nest::Node* target = get_target( t );

  // purely dendritic delay
  //float dendritic_delay = get_delay();
  double t_spike = e.get_stamp().get_ms();

  double new_cos2_, new_sin2_, new_cossin_;

  //std::cout << "Processing PF spike at time " << t_spike << "Last update: " << this->t_last_update_ << std::endl;

  this->weight_ += this->last_spike_weight_change_;

  // Check wether the weight stays within the boundaries
  this->weight_ = this->check_weight_boundaries(this->weight_);

  //std::cout << "Applying presynaptic spike weight change. New weight: " << this->weight_ << ".  Weight change: " << this->last_spike_weight_change_ << std::endl;

  //std::cout << "Sending spike in synapsis at time " << t_spike << std::endl;
  std::deque<mynest::histentry_cos>::iterator start;
  std::deque<mynest::histentry_cos>::iterator finish;
  //((mynest::Archiving_Node_Sym *)target_)->get_sym_history(t_lastspike - dendritic_delay, t_spike - dendritic_delay,&start, &finish);
  ((mynest::Archiving_Node_Cos *)target)->get_cos_history(this->t_last_update_, t_spike,&start, &finish);
  //weight change due to post-synaptic spikes since last pre-synaptic spike
  while (start != finish){

     // Evolve the state variables until the CS spike time
    this->evolve_cos_values( start->t_ - this->t_last_update_,
                              this->cos2_, this->sin2_, this->cossin_,
                              this->cos2_, this->sin2_, this->cossin_);

    this->t_last_update_ = start->t_;
    
    // Update the synaptic weight due to CS
    this->weight_ -= this->A_minus_*this->cos2_;
   
    // Check wether the weight stays within the boundaries
    this->weight_ = this->check_weight_boundaries(this->weight_);

    //std::cout << "Applying LTD with spike at time: " << start->t_ << ". New weight: " << this->weight_ << std::endl;

    ++start;
  }

  // Evolve the state variables until the presynaptic spike time
  this->evolve_cos_values( t_spike-this->t_last_update_,
                                  this->cos2_, this->sin2_, this->cossin_,
                                  this->cos2_, this->sin2_, this->cossin_);
  // Apply the effect of the incoming spike into the state variables
  this->cos2_ += 1.0;

  //std::cout << "Evolving state vars. cos2: " << this->cos2_ << " sin2: " << this->sin2_ << " cossin: " << this->cossin_ << std::endl;

  // Obtain weight change due to this spike (it will be applied when processing the
  // next presynaptic spike)
  ((mynest::Archiving_Node_Cos *)target)->get_cos_values( t_spike, new_cos2_, new_sin2_, new_cossin_);
  
  // Apply the LTD and LTP due to the previous presynaptic spike
  this->last_spike_weight_change_ = this->A_plus_ - this->A_minus_*new_cos2_;

  //std::cout << "Calculating next weight change: " << this->last_spike_weight_change_ << std::endl;

  this->t_last_update_ = t_spike;

  e.set_receiver( *target );
  e.set_weight( weight_ );
  e.set_delay( get_delay_steps() );
  e.set_rport( get_rport() );
  e();
}


template < typename targetidentifierT >
inline double STDPCosConnection< targetidentifierT >::check_weight_boundaries(double weight){
  if (weight_ > this->Wmax_){
    return this->Wmax_;
  } else if (weight < this->Wmin_) {
    return this->Wmin_;
  }
  
  return weight_;
}


} // of namespace nest

#endif // of #ifndef STDP_SIN_CONNECTION_H
