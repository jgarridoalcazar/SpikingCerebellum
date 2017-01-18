/*
 *  stdp_sin_connection.h
 */

#ifndef STDP_SIN_CONNECTION_H
#define STDP_SIN_CONNECTION_H

/* BeginDocumentation

   Name: stdp_sin_synapse - Synapse type for error-signal-driven spike-timing dependent
   plasticity.

   Description:
   stdp_sin_synapse is a connection to create synapses with
   error-signal-driven spike-timing dependent plasticity (implemented
   based on [1]). Every time the volume_transmitter process an spike,
   this synapsis will produce LTD according to a temporal kernel based on
   f(t) = exp(-t/tau_c)*sin(-t/tau_c)^exponent, where t represents the time
   since the last presynaptic spike in the synapses.

   Parameters:
     Common properties:
           exponent       unsigned int - Exponent of the sin function (integer between 1 and 20). The lower the exponent the wider the kernel function.
           peak    double - Time (in ms) of the peak of the kernel function (typically 100ms for the cerebellar parallel fibers).
           A_plus    double - Amplitude of weight change for facilitation
           A_minus   double - Amplitude of weight change for depression
           Wmin      double - Minimal synaptic weight
           Wmax      double - Maximal synaptic weight

     Individual properties:

   Remarks:
     The common properties can only be set by SetDefaults and apply to all synapses of
     the model.

   References:
   [1] Luque N. R., Garrido J. A., Carrillo R. R., Coenen O. J. Ros, E. (2011a). Cerebellarlike 
   corrective model inference engine for manipulation tasks. IEEE Trans. Syst. 
   Man Cybern. B Cybern. 41, 1299–1312. 10.1109/TSMCB.2011.2138693

   Transmits: SpikeEvent

   Author: Jesus Garrido
   Remarks:
   - based on previous code for EDLUT simulator (https://github.com/EduardoRosLab/edlut)
   
   SeeAlso: volume_transmitter
*/

#include "connection.h"
#include "volume_transmitter.h"
#include "spikecounter.h"
#include "numerics.h"

namespace nest
{

/**
 * Class containing the common properties for all synapses of type dopamine connection.
 */
class STDPSinCommonProperties : public CommonSynapseProperties
{
public:
  /**
   * Default constructor.
   * Sets all property values to defaults.
   */
  STDPSinCommonProperties();

  /**
   * Get all properties and put them into a dictionary.
   */
  void get_status( DictionaryDatum& d ) const;

  /**
   * Set properties from the values given in dictionary.
   */
  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  Node* get_node();

  double_t A_plus_;
  double_t A_minus_;
  double_t Peak_;
  int_t Exponent_;
  double_t inv_tau_
  double_t factor_;
  double_t Wmin_;
  double_t Wmax_;
};

/**
 * Class representing an STDPSinConnection with homogeneous parameters,
 * i.e. parameters are the same for all synapses.
 */
template < typename targetidentifierT >
class STDPSinConnection : public Connection< targetidentifierT >
{

public:
  typedef STDPSinCommonProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  STDPSinConnection();

  /**
   * Copy constructor from a property object.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  STDPSinConnection( const STDPSinConnection& );

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
  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   */
  void send( Event& e, thread t, double_t, const STDPSinCommonProperties& cp );

  void trigger_update_weight( thread t,
    const vector< spikecounter >& sin_spikes,
    double_t t_trig,
    const STDPSinCommonProperties& cp );

  class ConnTestDummyNode : public ConnTestDummyNodeBase
  {
  public:
    // Ensure proper overriding of overloaded virtual functions.
    // Return values from functions are ignored.
    using ConnTestDummyNodeBase::handles_test_event;
    port
    handles_test_event( SpikeEvent&, rport )
    {
      return invalid_port_;
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
  check_connection( Node& s,
    Node& t,
    rport receptor_type,
    double_t t_lastspike,
    const CommonPropertiesType& cp )
  {
    ConnTestDummyNode dummy_target;
    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );
  }

  void
  set_weight( double_t w )
  {
    weight_ = w;
  }

private:
  // data members of each connection
  double_t weight_;
  std::vector<double_t> state_vars_;
  double_t t_last_update_;
};

//
// Implementation of class STDPSinConnection.
//

template < typename targetidentifierT >
STDPSinConnection< targetidentifierT >::STDPSinConnection()
  : ConnectionBase()
  , weight_( 1.0 )
  , state_vars_( )
  , t_last_update_( 0.0 )
{
}

template < typename targetidentifierT >
STDPSinConnection< targetidentifierT >::STDPSinConnection( const STDPSinConnection& rhs )
  : ConnectionBase( rhs )
  , weight_( rhs.weight_ )
  , state_vars_( rhs.state_vars_ )
  , t_last_update_( rhs.t_last_update_ )
{
}

template < typename targetidentifierT >
void
STDPSinConnection< targetidentifierT >::get_status( DictionaryDatum& d ) const
{

  // base class properties, different for individual synapse
  ConnectionBase::get_status( d );
  def< double_t >( d, names::weight, weight_ );
}

template < typename targetidentifierT >
void
STDPSinConnection< targetidentifierT >::set_status( const DictionaryDatum& d, ConnectorModel& cm )
{
  // base class properties
  ConnectionBase::set_status( d, cm );
  updateValue< double_t >( d, names::weight, weight_ );
}

/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 * \param t_lastspike Time point of last spike emitted
 */
template < typename targetidentifierT >
inline void
STDPSinConnection< targetidentifierT >::send( Event& e,
  thread t,
  double_t,
  const STDPDopaCommonProperties& cp )
{
  // t_lastspike_ = 0 initially

  Node* target = get_target( t );

  // purely dendritic delay
  double_t dendritic_delay = get_delay();

  double_t t_spike = e.get_stamp().get_ms();

  // depression due to new pre-synaptic spike
  apply_presynaptic_spike_(cp);

  e.set_receiver( *target );
  e.set_weight( weight_ );
  e.set_delay( get_delay_steps() );
  e.set_rport( get_rport() );
  e();

  t_last_update_ = t_spike;
}

template < typename targetidentifierT >
inline void
STDPSinConnection< targetidentifierT >::trigger_update_weight( thread t,
  const vector< spikecounter >& sin_spikes,
  const double_t t_trig,
  const STDPSinCommonProperties& cp )
{
  apply_vt_spike_(sin_spikes, cp);
}

} // of namespace nest

#endif // of #ifndef STDP_SIN_CONNECTION_H
