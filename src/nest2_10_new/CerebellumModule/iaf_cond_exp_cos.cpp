/*
 *  iaf_cond_exp_cos.cpp
 *
 *  This file is based on the iaf_cond_exp cell model distributed with NEST.
 *  
 *  Modified by: Jesus Garrido (jgarridoalcazar at gmail.com) in 2017.
 */

#include "iaf_cond_exp_cos.h"

#ifdef HAVE_GSL

#include "exceptions.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include <limits>

#include "universal_data_logger_impl.h"
#include "event.h"

#include <iomanip>
#include <iostream>
#include <cstdio>

/* ---------------------------------------------------------------- 
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap<mynest::iaf_cond_exp_cos> mynest::iaf_cond_exp_cos::recordablesMap_;


namespace nest  // template specialization must be placed in namespace
{
  // Override the create() method with one call to RecordablesMap::insert_() 
  // for each quantity to be recorded.
  template <>
  void RecordablesMap<mynest::iaf_cond_exp_cos>::create()
  {
	  // use standard names whereever you can for consistency!
	  insert_(names::V_m,
		&mynest::iaf_cond_exp_cos::get_y_elem_<mynest::iaf_cond_exp_cos::State_::V_M>);
	  insert_(names::g_ex,
		&mynest::iaf_cond_exp_cos::get_y_elem_<mynest::iaf_cond_exp_cos::State_::G_EXC>);
	  insert_(names::g_in,
		&mynest::iaf_cond_exp_cos::get_y_elem_<mynest::iaf_cond_exp_cos::State_::G_INH>);
    insert_(names::g_ts,
    &mynest::iaf_cond_exp_cos::get_y_elem_<mynest::iaf_cond_exp_cos::State_::G_TS>);
  }
  
  namespace names
  {

  	  const Name tau_syn_ts("tau_syn_ts");
  	  const Name E_ts("e_ts");
      const Name g_ts("g_ts");
      const Name AMPA_R("AMPA");
      const Name GABA_R("GABA");
      const Name TEACHING_SIGNAL("TEACHING_SIGNAL");
  }
}

extern "C"
inline int mynest::iaf_cond_exp_cos_dynamics(double, const double y[], double f[], void* pnode)
{ 
  // a shorthand
  typedef mynest::iaf_cond_exp_cos::State_ S;

  // get access to node so we can almost work as in a member function
  assert(pnode);
  const mynest::iaf_cond_exp_cos& node =  *(reinterpret_cast<mynest::iaf_cond_exp_cos*>(pnode));

  // y[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.y[]. 

  // The following code is verbose for the sake of clarity. We assume that a
  // good compiler will optimize the verbosity away ...
  const double I_syn_exc = y[S::G_EXC] * (y[S::V_M] - node.P_.E_ex); 
  const double I_syn_inh = y[S::G_INH] * (y[S::V_M] - node.P_.E_in); 
  const double I_syn_ts = y[S::G_TS] * (y[S::V_M] - node.P_.E_ts); 
  const double I_L       = node.P_.g_L * ( y[S::V_M] - node.P_.E_L );
  const double I_total   = node.P_.I_e - I_syn_exc - I_syn_inh - I_syn_ts;
  
  //V dot
  f[0] = ( - I_L + node.B_.I_stim_ + I_total) / node.P_.C_m; // Vm diff. equation
  f[1] = -y[S::G_EXC] / node.P_.tau_synE; // Gexc diff. equation
  f[2] = -y[S::G_INH] / node.P_.tau_synI; // Ginh diff. equation
  f[3] = -y[S::G_TS] / node.P_.tau_synTS; // Gcs diff. equation

  return GSL_SUCCESS;
}

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
mynest::iaf_cond_exp_cos::Parameters_::Parameters_()
  : V_reset_   (-60.0    ),  // mV
  V_th_      (-55.0    ),  // mV
  t_ref_     (  2.0    ),  // ms
  g_L      ( 16.6667 ),  // nS
  C_m      (250.0    ),  // pF
  E_ex       (  0.0    ),  // mV
  E_in       (-85.0    ),  // mV
  E_ts       (  0.0    ),  // mV
  E_L        (-70.0    ),  // mV
	tau_synE   (  0.2    ),  // ms
  tau_synI   (  2.0    ),  // ms
  tau_synTS   (  5.0    ),  // ms
  I_e        (  0.0    )  // pA
{
}

mynest::iaf_cond_exp_cos::State_::State_(const Parameters_& p)
  : r_(0)
{
  y_[V_M] = p.E_L;
  y_[G_EXC] = y_[G_INH] = y_[G_TS] = 0;
}

mynest::iaf_cond_exp_cos::State_::State_(const State_& s)
  : r_(s.r_)
{
  for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
    y_[i] = s.y_[i];
}

mynest::iaf_cond_exp_cos::State_& mynest::iaf_cond_exp_cos::State_::operator=(const State_& s)
{
  assert(this != &s);  // would be bad logical error in program
  
  for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
    y_[i] = s.y_[i];
  r_ = s.r_;
  return *this;
}

/* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void mynest::iaf_cond_exp_cos::Parameters_::get(DictionaryDatum &d) const
{
	def<double>(d,nest::names::V_th,         V_th_);
	def<double>(d,nest::names::V_reset,      V_reset_);
	def<double>(d,nest::names::t_ref,        t_ref_);
	def<double>(d,nest::names::E_L,         E_L);
	def<double>(d,nest::names::g_L, 		 g_L);
	def<double>(d,nest::names::C_m,			 C_m);
	def<double>(d,nest::names::E_ex,         E_ex);
	def<double>(d,nest::names::E_in,         E_in);
  def<double>(d,nest::names::E_ts,         E_ts);
	def<double>(d,nest::names::tau_syn_ex,   tau_synE);
	def<double>(d,nest::names::tau_syn_in,   tau_synI);
  def<double>(d,nest::names::tau_syn_ts,   tau_synTS);
	def<double>(d,nest::names::I_e,          I_e);
}

void mynest::iaf_cond_exp_cos::Parameters_::set(const DictionaryDatum& d)
{
	// allow setting the membrane potential
	  updateValue<double>(d,nest::names::V_th,    V_th_);
	  updateValue<double>(d,nest::names::V_reset, V_reset_);
	  updateValue<double>(d,nest::names::t_ref,   t_ref_);
	  updateValue<double>(d,nest::names::E_L,     E_L);

	  updateValue<double>(d,nest::names::C_m, 	  C_m);
	  updateValue<double>(d,nest::names::g_L, 	  g_L);

	  updateValue<double>(d,nest::names::E_ex,    E_ex);
	  updateValue<double>(d,nest::names::E_in,    E_in);
    updateValue<double>(d,nest::names::E_ts,    E_ts);

	  updateValue<double>(d,nest::names::tau_syn_ex, tau_synE);
	  updateValue<double>(d,nest::names::tau_syn_in, tau_synI);
    updateValue<double>(d,nest::names::tau_syn_ts, tau_synTS);

	  updateValue<double>(d,nest::names::I_e,     I_e);

	// if ( V_reset_ >= V_th_ )
	//     throw nest::BadProperty("Reset potential must be smaller than threshold.");

	  if ( t_ref_ < 0 )
	    throw nest::BadProperty("Refractory time cannot be negative.");

	  if ( C_m <= 0 )
	      throw nest::BadProperty( "Capacitance must be strictly positive." );

	  if ( tau_synE <= 0 || tau_synI <= 0 || tau_synTS <= 0)
	    throw nest::BadProperty("All time constants must be strictly positive.");
}

void mynest::iaf_cond_exp_cos::State_::get(DictionaryDatum &d) const
{
	def<double>(d, nest::names::V_m, y_[V_M]); // Membrane potential
	def<double>(d, nest::names::g_ex, y_[G_EXC]);
  def<double>(d, nest::names::g_in, y_[G_INH]);
  def<double>(d, nest::names::g_ts, y_[G_TS]);
}

void mynest::iaf_cond_exp_cos::State_::set(const DictionaryDatum& d, const Parameters_&)
{
	updateValue<double>(d, nest::names::V_m, y_[V_M]);
	updateValue<double>(d, nest::names::g_ex, y_[G_EXC]);
  updateValue<double>(d, nest::names::g_in, y_[G_INH]);
  updateValue<double>(d, nest::names::g_ts, y_[G_TS]);
}

mynest::iaf_cond_exp_cos::Buffers_::Buffers_(mynest::iaf_cond_exp_cos& n)
  : logger_(n),
    s_(0),
    c_(0),
    e_(0)
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

mynest::iaf_cond_exp_cos::Buffers_::Buffers_(const Buffers_&, iaf_cond_exp_cos& n)
  : logger_(n),
    s_(0),
    c_(0),
    e_(0)
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ---------------------------------------------------------------- 
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

mynest::iaf_cond_exp_cos::iaf_cond_exp_cos()
  : mynest::Archiving_Node_Cos(),
    P_(), 
    S_(P_),
    B_(*this)
{
  recordablesMap_.create();
}

mynest::iaf_cond_exp_cos::iaf_cond_exp_cos(const iaf_cond_exp_cos& n)
  : mynest::Archiving_Node_Cos(n),
    P_(n.P_), 
    S_(n.S_),
    B_(n.B_, *this)
{
}

mynest::iaf_cond_exp_cos::~iaf_cond_exp_cos()
{
  // GSL structs may not have been allocated, so we need to protect destruction
  if ( B_.s_ ) gsl_odeiv_step_free(B_.s_);
  if ( B_.c_ ) gsl_odeiv_control_free(B_.c_);
  if ( B_.e_ ) gsl_odeiv_evolve_free(B_.e_);
}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void mynest::iaf_cond_exp_cos::init_state_(const Node& proto)
{
  const iaf_cond_exp_cos& pr = downcast<iaf_cond_exp_cos>(proto);
  S_ = pr.S_;
}

void mynest::iaf_cond_exp_cos::init_buffers_()
{
  B_.spike_exc_.clear();          // includes resize
  B_.spike_inh_.clear();          // includes resize
  B_.spike_ts_.clear();          // includes resize
  B_.currents_.clear();           // includes resize
  mynest::Archiving_Node_Cos::clear_history();

  B_.logger_.reset();

  B_.step_ = nest::Time::get_resolution().get_ms();
  B_.IntegrationStep_ = B_.step_;

  static const gsl_odeiv_step_type* T1 = gsl_odeiv_step_rkf45;
  
  if ( B_.s_ == 0 )
    B_.s_ = gsl_odeiv_step_alloc (T1, State_::STATE_VEC_SIZE);
  else 
    gsl_odeiv_step_reset(B_.s_);
    
  if ( B_.c_ == 0 )  
    B_.c_ = gsl_odeiv_control_y_new (1e-3, 0.0);
  else
    gsl_odeiv_control_init(B_.c_, 1e-3, 0.0, 1.0, 0.0);
    
  if ( B_.e_ == 0 )  
    B_.e_ = gsl_odeiv_evolve_alloc(State_::STATE_VEC_SIZE);
  else 
    gsl_odeiv_evolve_reset(B_.e_);
  
  B_.sys_.function  = iaf_cond_exp_cos_dynamics;
  B_.sys_.jacobian  = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params    = reinterpret_cast<void*>(this);

  B_.I_stim_ = 0.0;
}

void mynest::iaf_cond_exp_cos::calibrate()
{
  B_.logger_.init();  // ensures initialization in case mm connected after Simulate

  V_.RefractoryCounts_ = nest::Time(nest::Time::ms(P_.t_ref_)).get_steps();
  assert(V_.RefractoryCounts_ >= 0);  // since t_ref_ >= 0, this can only fail in error
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void mynest::iaf_cond_exp_cos::update(nest::Time const & origin, const long from, const long to)
{
   
  assert(to >= 0 && (nest::delay) from < nest::Scheduler::get_min_delay());
  assert(from < to);

  for ( long lag = from ; lag < to ; ++lag )
  {
    
    double t = 0.0;

    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t; this is of advantage
    // for a consistent and efficient integration across subsequent
    // simulation intervals
    while ( t < B_.step_ )
        {
          const int status = gsl_odeiv_evolve_apply(B_.e_, B_.c_, B_.s_,
    			   &B_.sys_,             // system of ODE
    			   &t,                   // from t
    			    B_.step_,            // to t <= step
    			   &B_.IntegrationStep_, // integration step size
    			    S_.y_); 	         // neuronal state
          if ( status != GSL_SUCCESS )
            throw nest::GSLSolverFailure(get_name(), status);
        }

        S_.y_[State_::G_EXC] += B_.spike_exc_.get_value(lag);
        S_.y_[State_::G_INH] += B_.spike_inh_.get_value(lag);
        S_.y_[State_::G_TS] += B_.spike_ts_.get_value(lag);

        // absolute refractory period
        if ( S_.r_ )
        {// neuron is absolute refractory
          --S_.r_;
          S_.y_[State_::V_M] = P_.V_reset_;
        }
        else
          // neuron is not absolute refractory
          if ( S_.y_[State_::V_M] >=  P_.V_th_)
    	    {
    	      S_.r_              = V_.RefractoryCounts_;
    	      S_.y_[State_::V_M] = P_.V_reset_;

            set_spiketime(nest::Time::step(origin.get_steps()+lag+1));

    	      nest::SpikeEvent se;
    	      network()->send(*this, se, lag);
    	    }
    
    // set new input current
    B_.I_stim_ = B_.currents_.get_value(lag);

    // log state data
    B_.logger_.record_data(origin.get_steps() + lag);

  }
}

void mynest::iaf_cond_exp_cos::handle(nest::SpikeEvent & e)
{
  assert(e.get_delay() > 0);

  assert(( e.get_rport() > INF_SPIKE_RECEPTOR ) && ( ( size_t ) e.get_rport() <= SUP_SPIKE_RECEPTOR ) );

  const long spike_time = e.get_rel_delivery_steps(network()->get_slice_origin());

  switch(e.get_rport()){
    case TEACHING_SIGNAL:
      B_.spike_ts_.add_value(spike_time, e.get_weight() * e.get_multiplicity() );
      set_cos_spiketime(nest::Time::step(network()->get_slice_origin().get_steps()+e.get_rel_delivery_steps(network()->get_slice_origin())));
      break;
    case AMPA:
      B_.spike_exc_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
          e.get_weight() * e.get_multiplicity() );
      break;
    case GABA:
      B_.spike_inh_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
          e.get_weight() * e.get_multiplicity() );  // ensure conductance is positive
      break;
    default:
      break;
    }
  

  
}

void mynest::iaf_cond_exp_cos::handle(nest::CurrentEvent& e)
{
  assert(e.get_delay() > 0);

  const double c=e.get_current();
  const double w=e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
		      w *c);
}

void mynest::iaf_cond_exp_cos::handle(nest::DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}

#endif //HAVE_GSL
