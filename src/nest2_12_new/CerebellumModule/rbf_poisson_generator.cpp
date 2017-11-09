/*
 *  rbf_poisson_generator.cpp
 *
 *  This file is based on the iaf_cond_exp cell model distributed with NEST.
 *  
 *  Modified by: Jesus Garrido (jgarridoalcazar at gmail.com) in 2017.
 */

#include "rbf_poisson_generator.h"

// Includes from nestkernel:
#include "event_delivery_manager_impl.h"
#include "exceptions.h"
#include "kernel_manager.h"

#include "universal_data_logger_impl.h"


// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"


/* ---------------------------------------------------------------- 
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap<mynest::rbf_poisson_generator> mynest::rbf_poisson_generator::recordablesMap_;


namespace nest  // template specialization must be placed in namespace
{
  // Override the create() method with one call to RecordablesMap::insert_() 
  // for each quantity to be recorded.
  template <>
  void RecordablesMap<mynest::rbf_poisson_generator>::create()
  {
    // use standard names whereever you can for consistency!
    insert_(names::rate,
    &mynest::rbf_poisson_generator::get_rate_);

    // use standard names whereever you can for consistency!
    insert_(names::I,
    &mynest::rbf_poisson_generator::get_current_);
  }
  
  namespace names
  {

  }
}


/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
mynest::rbf_poisson_generator::Parameters_::Parameters_()
   : min_rate_rbf_( 1.0 ) // Hz
   , max_rate_rbf_( 10.0 ) // Hz
   , mean_current_rbf_( 0.0 ) // nA
   , sigma_current_rbf_( 1.0 ) // nA
{
}

mynest::rbf_poisson_generator::State_::State_(const Parameters_& p)
  : rate_(5.0),
  input_current_(0.0)  
{
}

mynest::rbf_poisson_generator::State_::State_(const State_& s)
  : rate_(s.rate_),
  input_current_(s.input_current_)
{
}

/* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void mynest::cd_poisson_generator::Parameters_::get(DictionaryDatum &d) const
{
	def< double >( d, nest::names::min_rate, min_rate_ );
  def< double >( d, nest::names::max_rate, max_rate_ );
  def< double >( d, nest::names::min_current, min_current_ );
  def< double >( d, nest::names::max_current, max_current_ );
}

void mynest::cd_poisson_generator::Parameters_::set(const DictionaryDatum& d)
{
	def< double >( d, nest::names::min_rate_rbf, min_rate_rbf_ );
  def< double >( d, nest::names::max_rate_rbf, max_rate_rbf_ );
  def< double >( d, nest::names::mean_current_rbf, mean_current_rbf_ );
  def< double >( d, nest::names::sigma_current_rbf, sigma_current_rbf_ );
  if ( min_rate_rbf_ < 0 || max_rate_rbf_ < 0 || min_rate_rbf_ > max_rate_rbf_ )
  {
    throw nest::BadProperty( "The min_rate and max_rate parameters cannot be negative." );
  }
}

void mynest::cd_poisson_generator::State_::get(DictionaryDatum &d) const
{
  def<double>(d, nest::names::rate, rate_);
  def<double>(d, nest::names::I, input_current_);
}

void mynest::cd_poisson_generator::State_::set(const DictionaryDatum& d, const Parameters_&)
{
  updateValue<double>(d, nest::names::rate, rate_);
  updateValue<double>(d, nest::names::I, input_current_);
}

mynest::cd_poisson_generator::Buffers_::Buffers_(mynest::cd_poisson_generator& n)
  : logger_(n)
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

mynest::cd_poisson_generator::Buffers_::Buffers_(const Buffers_&, cd_poisson_generator& n)
  : logger_(n)
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}


/* ---------------------------------------------------------------- 
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

mynest::cd_poisson_generator::cd_poisson_generator()
  : Node()
  , P_()
  , S_(P_)
  , B_(*this)
{
  recordablesMap_.create();
}

mynest::cd_poisson_generator::cd_poisson_generator(const cd_poisson_generator& n)
  : Node( n )
  , P_( n.P_ )
  , S_( n.S_)
  , B_( n.B_, *this)
{
}

mynest::cd_poisson_generator::~cd_poisson_generator()
{
}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void mynest::cd_poisson_generator::init_state_(const Node& proto)
{
  const cd_poisson_generator& pr = downcast< cd_poisson_generator >( proto );

  S_ = pr.S_;
}

void mynest::cd_poisson_generator::init_buffers_()
{
  B_.currents_.clear();

  B_.logger_.reset();

  B_.step_ = nest::Time::get_resolution().get_ms();
}

void mynest::cd_poisson_generator::calibrate()
{
  B_.logger_.init();

  double gaussian = std::exp(-((S_.input_current_ - P_.mean_current_rbf_) * (S_.input_current_ - P_.mean_current_rbf_))/((2.0 * P_.sigma_current_rbf_ * P_.sigma_current_rbf_)));
  S_.rate_ = P_.min_rate_rbf_ + gaussian * (P_.max_rate_rbf_ - P_.min_rate_rbf_);

  // rate_ is in Hz, dt in ms, so we have to convert from s to ms
  V_.poisson_dev_.set_lambda(
    nest::Time::get_resolution().get_ms() * S_.rate_ * 1e-3 );
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void mynest::cd_poisson_generator::update(nest::Time const & T, const long from, const long to)
{
   
  assert(
    to >= 0 && ( nest::delay ) from < nest::kernel().connection_manager.get_min_delay() );
  assert( from < to );

  librandom::RngPtr rng = nest::kernel().rng_manager.get_rng( get_thread() );

  for ( long lag = from; lag < to; ++lag)
  {
    double new_current = B_.currents_.get_value(lag);

    // Update the firing rate only when the input current changes
    if (S_.input_current_!= new_current){

      S_.input_current_ = new_current;
      
      double gaussian = std::exp(-((S_.input_current_ - P_.mean_current_rbf_) * (S_.input_current_ - P_.mean_current_rbf_))/((2.0 * P_.sigma_current_rbf_ * P_.sigma_current_rbf_)));
      S_.rate_ = P_.min_rate_rbf_ + gaussian * (P_.max_rate_rbf_ - P_.min_rate_rbf_);

      // Lambda device frequency is updated only once per call to update.
      V_.poisson_dev_.set_lambda(
        nest::Time::get_resolution().get_ms() * S_.rate_ * 1e-3 );
    }

      
    if (S_.rate_ > 0.0){
      long n_spikes = V_.poisson_dev_.ldev( rng );

      if ( n_spikes > 0 ) // we must not send events with multiplicity 0
      {
        nest::SpikeEvent e;
        e.set_multiplicity( n_spikes );
        nest::kernel().event_delivery_manager.send( *this, e, lag );
      }
    }
 
    // log state data
    B_.logger_.record_data(T.get_steps() + lag);
  }


}

void mynest::cd_poisson_generator::handle(nest::CurrentEvent& e)
{
  assert(e.get_delay() > 0);

  const double c=e.get_current();
  const double w=e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(e.get_rel_delivery_steps(nest::kernel().simulation_manager.get_slice_origin()), 
          w *c);
}

void mynest::cd_poisson_generator::handle(nest::DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}

