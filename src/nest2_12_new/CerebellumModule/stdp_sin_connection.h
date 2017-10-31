/*
 *  stdp_sin_connection.h
 */

#ifndef STDP_SIN_CONNECTION_H
#define STDP_SIN_CONNECTION_H

/* BeginDocumentation

   Name: stdp_sin_synapse - Synapse type for complex-spike-driven spike-timing dependent
   plasticity.

   Description:
   stdp_sin_synapse is a connection to create synapses with
   complex-spike-driven spike-timing dependent plasticity (implemented
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
      peak      double - Time (in ms) of the peak of the kernel function (typically 100ms for the cerebellar parallel fibers).
         

  References:
   [1] Luque N. R., Garrido J. A., Carrillo R. R., Coenen O. J. Ros, E. (2011a). Cerebellarlike 
   corrective model inference engine for manipulation tasks. IEEE Trans. Syst. 
   Man Cybern. B Cybern. 41, 1299–1312. 10.1109/TSMCB.2011.2138693

   Transmits: SpikeEvent

   Author: Jesus Garrido
   Remarks:
   - based on previous code for EDLUT simulator (https://github.com/EduardoRosLab/edlut)
   
   SeeAlso: iaf_cond_exp_cs
*/

#include "common_synapse_properties.h"

#include "connection.h"
#include "archiving_node_cs.h"

#include "ExponentialTable.h"
#include "TrigonometricTable.h"

#define A 1.0f/2.0f

namespace mynest
{
/**
 * Class representing an STDPSinConnection.
 */
template < typename targetidentifierT >
class STDPSinConnection : public nest::Connection< targetidentifierT >
{

private:
  static float terms[11][11];

public:
  typedef nest::CommonSynapseProperties CommonPropertiesType;
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
  void send( nest::Event& e, nest::thread t, float, const nest::CommonSynapseProperties& cp );

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
    ConnTestDummyNode dummy_target;
    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );

    ((Archiving_Node_CS *) (&t))->register_stdp_connection_cs( t_lastspike - get_delay() );
  }

  void
  set_weight( double_t w )
  {
    weight_ = w;
  }

private:
  // data members of each connection
  double weight_;
  std::vector<double> state_vars_;
  double t_last_update_;

  // This vars could be common, but state_vars has to be resized according to Exponent.
  // Remaining vars depends on Exponent too.
  double Peak_;
  unsigned short int Exponent_;
  double inv_tau_;
  double factor_;
  float * TermPointer_;
  
  // This vars could be also common, but this optimization might be performed as a future development.
  double A_plus_;
  double A_minus_;
  double Wmin_;
  double Wmax_;

  void apply_state_change(double new_time);

  double check_weight_boundaries(double weight);
};

template < typename targetidentifierT >
float STDPSinConnection< targetidentifierT >::terms[11][11] =
  {{1,0,0,0,0,0,0,0,0,0,0},
   {A,-A,0,0,0,0,0,0,0,0,0},
   {3.0f/2.0f*pow(A,2),-4.0f/2.0f*pow(A,2),1.0f/2.0f*pow(A,2),0,0,0,0,0,0,0,0},
   {10.0f/4.0f*pow(A,3),-15.0f/4.0f*pow(A,3),6.0f/4.0f*pow(A,3),-1.0f/4.0f*pow(A,3),0,0,0,0,0,0,0},
   {35.0f/8.0f*pow(A,4),-56.0f/8.0f*pow(A,4),28.0f/8.0f*pow(A,4),-8.0f/8.0f*pow(A,4),1.0f/8.0f*pow(A,4),0,0,0,0,0,0},
   {126.0f/16.0f*pow(A,5),-210.0f/16.0f*pow(A,5),120.0f/16.0f*pow(A,5),-45.0f/16.0f*pow(A,5),10.0f/16.0f*pow(A,5),-1.0f/16.0f*pow(A,5),0,0,0,0,0},
   {231.0f/16.0f*pow(A,6),-99.0f/4.0f*pow(A,6),495.0f/32.0f*pow(A,6),-55.0f/8.0f*pow(A,6),66.0f/32.0f*pow(A,6),-3.0f/8.0f*pow(A,6),1.0f/32.0f*pow(A,6),0,0,0,0},
   {429.0f/16.0f*pow(A,7),-3003.0f/64.0f*pow(A,7),1001.0f/32.0f*pow(A,7),-1001.0f/64.0f*pow(A,7),91.0f/16.0f*pow(A,7),-91.0f/64.0f*pow(A,7),7.0f/32.0f*pow(A,7),-1.0f/64.0f*pow(A,7),0,0,0},
   {6435.0f/128.0f*pow(A,8),-715.0f/8.0f*pow(A,8),1001.0f/16.0f*pow(A,8),-273.0f/8.0f*pow(A,8),455.0f/32.0f*pow(A,8),-35.0f/8.0f*pow(A,8),15.0f/16.0f*pow(A,8),-1.0f/8.0f*pow(A,8),1.0f/128.0f*pow(A,8),0,0},
   {12155.0f/128.0f*pow(A,9),-21879.0f/128.0f*pow(A,9),1989.0f/16.0f*pow(A,9),-4641.0f/64.0f*pow(A,9),1071.0f/32.0f*pow(A,9),-765.0f/64.0f*pow(A,9),51.0f/16.0f*pow(A,9),-153.0f/256.0f*pow(A,9),9.0f/128.0f*pow(A,9),-1.0f/256.0f*pow(A,9),0},
   {46189.0f/256.0f*pow(A,10),-20995.0f/64.0f*pow(A,10),62985.0f/256.0f*pow(A,10),-4845.0f/32.0f*pow(A,10),4845.0f/64.0f*pow(A,10),-969.0f/32.0f*pow(A,10),4845.0f/512.0f*pow(A,10),-285.0f/128.0f*pow(A,10),95.0f/256.0f*pow(A,10),-5.0f/128.0f*pow(A,10),1.0f/512.0f*pow(A,10)}};


//
// Implementation of class STDPSinConnection.
//

template < typename targetidentifierT >
STDPSinConnection< targetidentifierT >::STDPSinConnection()
  : ConnectionBase(),
  weight_( 1.0 ),
  state_vars_( ),
  t_last_update_( 0.0 ),
  Peak_( 100.0 ),
  Exponent_( 2 ),
  A_plus_( 1.0 ),
  A_minus_( 1.0 ),
  Wmin_( 0.0 ),
  Wmax_( 200.0 )
{
  this->state_vars_ = std::vector<double>(this->Exponent_+2);
  inv_tau_ = atan((float) this->Exponent_)/Peak_;
  factor_ = 1.0f/(exp(-atan((float)this->Exponent_))*pow(sin(atan((float)this->Exponent_)),(int) this->Exponent_));

  unsigned int ExponenLine = this->Exponent_/2;
  TermPointer_ = STDPSinConnection< targetidentifierT >::terms[ExponenLine]; 
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
  , A_plus_( rhs.A_plus_ )
  , A_minus_( rhs.A_minus_ )
  , Wmin_( rhs.Wmin_ )
  , Wmax_( rhs.Wmax_ )
{
  unsigned int ExponenLine = this->Exponent_/2;
  TermPointer_ = STDPSinConnection< targetidentifierT >::terms[ExponenLine];
}

template < typename targetidentifierT >
void
STDPSinConnection< targetidentifierT >::get_status( DictionaryDatum& d ) const
{

  // base class properties, different for individual synapse
  ConnectionBase::get_status( d );
  def< double >( d, "A_plus", A_plus_ );
  def< double >( d, "A_minus", A_minus_ );
  def< double >( d, "Wmin", Wmin_ );
  def< double >( d, "Wmax", Wmax_ );
  def< double >( d, nest::names::weight, this->weight_ );

  def< double >( d, "peak", this->Peak_ );
  def< double >( d, "exponent", this->Exponent_ );
}

template < typename targetidentifierT >
void
STDPSinConnection< targetidentifierT >::set_status( const DictionaryDatum& d, nest::ConnectorModel& cm )
{
  // base class properties
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, nest::names::weight, weight_ );

  updateValue< double >( d, "A_plus", A_plus_ );
  updateValue< double >( d, "A_minus", A_minus_ );
  
  updateValue< double >( d, "Wmin", Wmin_ );
  updateValue< double >( d, "Wmax", Wmax_ );

  long expon;
  if ( updateValue< double >( d, "exponent", expon ))
  { 
    double intpart;
    if (modf(expon, &intpart) != 0.0){
      throw nest::BadProperty( "STDP sin exponent must be an even integer between 2 and 20" ); 
    }

    if (expon<2 || expon>20 || expon%2!=0)
    {
      throw nest::BadProperty( "STDP sin exponent must be an even integer between 2 and 20" );
    }
    this->Exponent_ = (unsigned short int) expon;

    unsigned int ExponenLine = Exponent_/2;
    TermPointer_ = STDPSinConnection< targetidentifierT >::terms[ExponenLine]; 

    this->state_vars_.resize(this->Exponent_+2);
  }

  updateValue< double >( d, "peak", this->Peak_ );

  //std::cout << "Synapse parameters: Aplus: " << this->A_plus_ << ". Aminus: " << this->A_minus_ << ". Wmin: " << this->Wmin_ << ". Wmax: " << this->Wmax_ << ". Expon: " << this->Exponent_ << ". Peak: " << this->Peak_ << std::endl;
  
  this->inv_tau_ = atan((float) this->Exponent_)/this->Peak_;
  this->factor_ = 1.0f/(exp(-atan((float)this->Exponent_))*pow(sin(atan((float)this->Exponent_)),(int) this->Exponent_)); 
}

template < typename targetidentifierT >
void STDPSinConnection< targetidentifierT >::apply_state_change(double new_time){

  // Evolve all the state variables from last_cs_time until t_cs
  double OldExpon = this->state_vars_[1];

  double ElapsedTime = double(new_time - this->t_last_update_);
  double ElapsedRelative = ElapsedTime*this->inv_tau_;

  double expon = ExponentialTable::GetResult(-ElapsedRelative);

  this->t_last_update_ = new_time;
  
  double NewExpon = OldExpon * expon;
  double NewActivity =NewExpon*this->TermPointer_[0];

  int aux=TrigonometricTable::CalculateOffsetPosition(2*ElapsedRelative);
  int LUTindex=0;

  double SinVar, CosVar, OldVarCos, OldVarSin, NewVarCos, NewVarSin;
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
  const nest::CommonSynapseProperties& )
{
  // t_lastspike_ = 0 initially

  //std::cout << "Synapse parameters: Aplus: " << this->A_plus_ << ". Aminus: " << this->A_minus_ << ". Wmin: " << this->Wmin_ << ". Wmax: " << this->Wmax_ << ". Expon: " << this->Exponent_ << ". Peak: " << this->Peak_ << ". Weight: " << this->weight_ << std::endl;
  
  nest::Node* target = get_target( t );

  // purely dendritic delay
  //float dendritic_delay = get_delay();
  double t_spike = e.get_stamp().get_ms();

  //std::cout << "Processing PF spike at time " << t_spike << std::endl;

  if (this->t_last_update_>0.0){
    // Apply the LTP due to the previous presynaptic spike
    this->weight_ += this->A_plus_;

    // Check wether the weight stays within the boundaries
    this->weight_ = this->check_weight_boundaries(this->weight_);

    //std::cout << "Applying LTP. New weight: " << this->weight_ << ".  Aplus: " << this->A_plus_ << std::endl;

  }  
  
  //std::cout << "Sending spike in synapsis at time " << t_spike << std::endl;
  std::deque<mynest::histentry_cs>::iterator start;
  std::deque<mynest::histentry_cs>::iterator finish;
  //((mynest::Archiving_Node_Sym *)target_)->get_sym_history(t_lastspike - dendritic_delay, t_spike - dendritic_delay,&start, &finish);
  ((mynest::Archiving_Node_CS *)target)->get_cs_history(this->t_last_update_, t_spike,&start, &finish);
  //weight change due to post-synaptic spikes since last pre-synaptic spike
  while (start != finish){

     // Evolve the state variables until the CS spike time
     this->apply_state_change(start->t_);

     // Update the synaptic weight due to CS
     this->weight_ -= this->A_minus_*this->state_vars_[0];
   
     // Check wether the weight stays within the boundaries
     this->weight_ = this->check_weight_boundaries(this->weight_);

     //std::cout << "Applying LTD with spike at time: " << start->t_ << ". New weight: " << this->weight_ << std::endl;

     ++start;
  }

  // Evolve the state variables until the presynaptic spike time
  this->apply_state_change(t_spike);

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
inline double STDPSinConnection< targetidentifierT >::check_weight_boundaries(double weight){
  if (weight_ > this->Wmax_){
    return this->Wmax_;
  } else if (weight < this->Wmin_) {
    return this->Wmin_;
  }
  
  return weight_;
}

} // of namespace nest

#endif // of #ifndef STDP_SIN_CONNECTION_H
