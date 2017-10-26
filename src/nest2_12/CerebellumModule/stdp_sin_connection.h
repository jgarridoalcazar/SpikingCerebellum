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
           vt        long   - ID of volume_transmitter collecting the spikes from the pool of
                              dopamine releasing neurons and transmitting the spikes
                              to the synapse. A value of -1 indicates that no volume
                              transmitter has been assigned.
           A_plus    double - Amplitude of weight change for facilitation
           A_minus   double - Amplitude of weight change for depression
           Wmin      double - Minimal synaptic weight
           Wmax      double - Maximal synaptic weight
           

     Individual properties:
           exponent       unsigned int - Exponent of the sin function (integer between 1 and 20). The lower the exponent the wider the kernel function.
           peak    double - Time (in ms) of the peak of the kernel function (typically 100ms for the cerebellar parallel fibers).
           

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

#include "ExponentialTable.h"
#include "TrigonometricTable.h"

#define A 1.0f/2.0f

namespace mynest
{

/**
 * Class containing the common properties for all synapses of type dopamine connection.
 */
class STDPSinCommonProperties : public nest::CommonSynapseProperties
{
public:
  static float terms[11][11];

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
  void set_status( const DictionaryDatum& d, nest::ConnectorModel& cm );

  nest::Node* get_node();

  long get_vt_gid() const;

  mynest::volume_transmitter* vt_;
  float A_plus_;
  float A_minus_;
  float Wmin_;
  float Wmax_;
};

inline long
STDPSinCommonProperties::get_vt_gid() const
{
  if ( vt_ != 0 )
    return vt_->get_gid();
  else
    return -1;
}

/**
 * Class representing an STDPSinConnection with homogeneous parameters,
 * i.e. parameters are the same for all synapses.
 */
template < typename targetidentifierT >
class STDPSinConnection : public nest::Connection< targetidentifierT >
{

public:
  typedef STDPSinCommonProperties CommonPropertiesType;
  typedef nest::Connection< targetidentifierT > ConnectionBase;

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
  void set_status( const DictionaryDatum& d, nest::ConnectorModel& cm );

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   */
  void send( nest::Event& e, nest::thread t, float, const STDPSinCommonProperties& cp );

  void trigger_update_weight( nest::thread t,
    const std::vector< nest::spikecounter >& sin_spikes,
    float t_trig,
    const STDPSinCommonProperties& cp );

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
    float t_lastspike,
    const CommonPropertiesType& cp )
  {
    if ( cp.vt_ == 0 )
      throw nest::BadProperty( "No volume transmitter has been assigned to the STDPSin synapse." );
    
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
  float weight_;
  std::vector<float> state_vars_;
  float t_last_update_;

  // This vars could be common, but state_vars has to be resized according to Exponent.
  // Remaining vars depends on Exponent too.
  float Peak_;
  unsigned short int Exponent_;
  float inv_tau_;
  float factor_;
  float * TermPointer_;

  void update_synaptic_state(float t_spike);
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
  , Peak_( 100.0 )
  , Exponent_( 2 )
  , TermPointer_( 0 )
{
  this->state_vars_ = std::vector<float>(this->Exponent_+2);
  inv_tau_ = atan((float) this->Exponent_)/Peak_;
  factor_ = 1.0f/(exp(-atan((float)this->Exponent_))*pow(sin(atan((float)this->Exponent_)),(int) this->Exponent_));

  unsigned int ExponenLine = this->Exponent_/2;
  TermPointer_ = STDPSinCommonProperties::terms[ExponenLine]; 
}

template < typename targetidentifierT >
STDPSinConnection< targetidentifierT >::STDPSinConnection( const STDPSinConnection& rhs )
  : ConnectionBase( rhs )
  , weight_( rhs.weight_ )
  , state_vars_( rhs.state_vars_ )
  , t_last_update_( rhs.t_last_update_ )
  , Peak_( rhs.Peak_ )
  , Exponent_( rhs.Exponent_ )
  , inv_tau_( rhs.inv_tau_ )
  , factor_( rhs.factor_ )
  , TermPointer_( rhs.TermPointer_ )
{
}

template < typename targetidentifierT >
void
STDPSinConnection< targetidentifierT >::get_status( DictionaryDatum& d ) const
{

  // base class properties, different for individual synapse
  ConnectionBase::get_status( d );
  def< float >( d, nest::names::weight, this->weight_ );

  def< float >( d, "peak", this->Peak_ );
  def< long >( d, "exponent", this->Exponent_ );
}

template < typename targetidentifierT >
void
STDPSinConnection< targetidentifierT >::set_status( const DictionaryDatum& d, nest::ConnectorModel& cm )
{
  // base class properties
  ConnectionBase::set_status( d, cm );
  updateValue< float >( d, nest::names::weight, weight_ );

  long expon;
  if ( updateValue< long >( d, "exponent", expon ))
  { 
    if (expon<2 || expon>20 || expon%2!=0)
    {
      throw nest::BadProperty( "STDP sin exponent must be an even integer between 2 and 20" );
    }
    this->Exponent_ = (unsigned short int) expon;

    unsigned int ExponenLine = Exponent_/2;
    TermPointer_ = STDPSinCommonProperties::terms[ExponenLine]; 

    this->state_vars_.resize(this->Exponent_+2);
  }

  updateValue< float >( d, "peak", this->Peak_ );
  
  this->inv_tau_ = atan((float) this->Exponent_)/this->Peak_;
  this->factor_ = 1.0f/(exp(-atan((float)this->Exponent_))*pow(sin(atan((float)this->Exponent_)),(int) this->Exponent_)); 
}

template < typename targetidentifierT >
void STDPSinConnection< targetidentifierT >::update_synaptic_state(float t_spike){
  float OldExpon = this->state_vars_[1];

  float ElapsedTime = float(t_spike - this->t_last_update_);
  float ElapsedRelative = ElapsedTime*this->inv_tau_;

  float expon = ExponentialTable::GetResult(-ElapsedRelative);

  this->t_last_update_ = t_spike;
  
  float NewExpon = OldExpon * expon;
  float NewActivity =NewExpon*this->TermPointer_[0];

  int aux=TrigonometricTable::CalculateOffsetPosition(2*ElapsedRelative);
  int LUTindex=0;

  float SinVar, CosVar, OldVarCos, OldVarSin, NewVarCos, NewVarSin;
  int grade, offset;
  for (grade=2, offset=1; grade<=this->Exponent_; grade+=2, offset++){

    LUTindex =TrigonometricTable::CalculateValidPosition(LUTindex,aux);

    OldVarCos = this->state_vars_[grade];
    OldVarSin = this->state_vars_[grade + 1];

    SinVar = TrigonometricTable::GetElement(LUTindex);
    CosVar = TrigonometricTable::GetElement(LUTindex+1);

    NewVarCos = (OldVarCos*CosVar-OldVarSin*SinVar)*expon;
    NewVarSin = (OldVarSin*CosVar+OldVarCos*SinVar)*expon;

    NewActivity+= NewVarCos*this->TermPointer_[offset];

    this->state_vars_[grade] = NewVarCos;
    this->state_vars_[grade+1] = NewVarSin;
  }
  NewActivity*=this->factor_;
  this->state_vars_[0] = NewActivity;
  this->state_vars_[1] = NewExpon;
}

/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 * \param t_lastspike Time point of last spike emitted
 */
template < typename targetidentifierT >
inline void
STDPSinConnection< targetidentifierT >::send( nest::Event& e,
  nest::thread t,
  float,
  const STDPSinCommonProperties& cp )
{
  // t_lastspike_ = 0 initially

  nest::Node* target = get_target( t );

  // purely dendritic delay
  //float dendritic_delay = get_delay();

  float t_spike = e.get_stamp().get_ms();

  //std::cout << "Sending spike in synapsis at time " << t_spike << std::endl;

  // Increment the weight a fixed amount A_plus_
  this->weight_ += cp.A_plus_;
  if (this->weight_ > cp.Wmax_){
    this->weight_ = cp.Wmax_;
  }
  if (this->weight_ < cp.Wmin_){
    this->weight_ = cp.Wmin_;
  }
  
  // --------------------------------------------
  // Update the state vars until the current time
  this->update_synaptic_state(t_spike);

  // -----------------------------------
  // Add the effect of the presynaptic spike to the state vars
  this->state_vars_[1] += 1.0f;
  for (unsigned int grade=2; grade<=this->Exponent_; grade+=2){
    this->state_vars_[grade] += 1.0f;
  }

  e.set_receiver( *target );
  e.set_weight( weight_ );
  e.set_delay( get_delay_steps() );
  e.set_rport( get_rport() );
  e();

}

template < typename targetidentifierT >
inline void
STDPSinConnection< targetidentifierT >::trigger_update_weight( nest::thread t,
  const std::vector< nest::spikecounter >& sin_spikes,
  const float t_trig,
  const STDPSinCommonProperties& cp )
{
  //std:: cout << "Triggering update weight at time " << t_trig << std::endl;
  
  // Update the state and weight for each spike in coming from the volume transmitter
  for (std::vector< nest::spikecounter>::const_iterator it = sin_spikes.begin(); it!=sin_spikes.end(); ++it){
    float t_spike = it->spike_time_;

    this->update_synaptic_state(t_spike);

    // Update weight due to domaminergic input
    this->weight_ -= cp.A_minus_*this->state_vars_[0]*it->multiplicity_;
  }

  if (this->weight_ > cp.Wmax_){
    this->weight_ = cp.Wmax_;
  }
  if (this->weight_ < cp.Wmin_){
    this->weight_ = cp.Wmin_;
  }
}

} // of namespace nest

#endif // of #ifndef STDP_SIN_CONNECTION_H
