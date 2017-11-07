/*
 *  rbf_poisson_generator.h
 *
 *  This file is based on the poisson generator model distributed with NEST.
 *  
 *  Modified by: Jesús Garrido (jgarridoalcazar at gmail.com) in 2017.
 */

#ifndef RBF_POISSON_GENERATOR_H
#define RBF_POISSON_GENERATOR_H

// Includes from librandom:
#include "poisson_randomdev.h"

// Includes from nestkernel:
#include "connection.h"
#include "event.h"
#include "nest.h"
#include "node.h"
#include "stimulating_device.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"

/* BeginDocumentation
Name: rbf_poisson_generator - simulate neuron firing with Poisson processes
                          statistics drive by input current.
Description:
  The rbf_poisson_generator simulates a neuron that is firing with Poisson
  statistics, i.e. exponentially distributed interspike intervals. Its firing
  rate is linearly calculated based on the total amount of input current. It will
  generate a _unique_ spike train for each of it's targets. If you do not want
  this behavior and need the same spike train for all targets, you have to use a
  parrot neuron inbetween the poisson generator and the targets.

Parameters:
   The following parameters appear in the element's status dictionary:

   origin   double - Time origin for device timer in ms
   start    double - begin of device application with resp. to origin in ms
   stop     double - end of device application with resp. to origin in ms
   min_rate double - Min firing rate in Hz
   max_rate double - Max firing rate in Hz
   min_current    double - The firing rate will be min_rate when the input
                      current is below this parameter.
   max_current    double - The firing rate will be max_rate when the input
                      current is above this parameter.

Sends: SpikeEvent

Remarks:
   A Poisson generator may, especially at high rates, emit more than one
   spike during a single time step. If this happens, the generator does
   not actually send out n spikes. Instead, it emits a single spike with
   n-fold synaptic weight for the sake of efficiency.

   The design decision to implement the Poisson generator as a device
   which sends spikes to all connected nodes on every time step and then
   discards the spikes that should not have happened generating random
   numbers at the recipient side via an event hook is twofold.

   On one hand, it leads to the saturation of the messaging network with
   an enormous amount of spikes, most of which will never get delivered
   and should not have been generated in the first place.

   On the other hand, a proper implementation of the Poisson generator
   needs to provide two basic features: (a) generated spike trains
   should be IID processes w.r.t. target neurons to which the generator
   is connected and (b) as long as virtual_num_proc is constant, each
   neuron should receive an identical Poisson spike train in order to
   guarantee reproducibility of the simulations across varying machine
   numbers.

   Therefore, first, as Network::get_network().send sends spikes to all the
   recipients, differentiation has to happen in the hook, second, the
   hook can use the RNG from the thread where the recipient neuron sits,
   which explains the current design of the generator. For details,
   refer to:

SeeAlso: poisson_generator, Device, parrot_neuron
*/

// Define name constants for state variables and parameters
namespace nest
{
	namespace names
	{
    	// Neuron parameters
    	extern const Name min_rate_rbf;  
    	extern const Name max_rate_rbf;  
      extern const Name mean_current_rbf;
      extern const Name sigma_current_rbf;
    }
}

namespace mynest
{
  /**
   * Function computing right-hand side of ODE for GSL solver.
   * @note Must be declared here so we can befriend it in class.
   * @note Must have C-linkage for passing to GSL. Internally, it is
   *       a first-class C++ function, but cannot be a member function
   *       because of the C-linkage.
   * @note No point in declaring it inline, since it is called
   *       through a function pointer.
   * @param void* Pointer to model neuron instance.
   */
  //extern "C"
  //int rbf_poisson_generator_dynamics (double, const double*, double*, void*);
  
  class rbf_poisson_generator : public nest::Node
  {
    
  public:        
    
    rbf_poisson_generator();
    rbf_poisson_generator(const rbf_poisson_generator&);
    ~rbf_poisson_generator();

  bool has_proxies() const
  {
    return false;
  }

  bool local_receiver() const {
    return true;
  }


    /**
     * Import sets of overloaded virtual functions.
     * We need to explicitly include sets of overloaded
     * virtual functions into the current scope.
     * According to the SUN C++ FAQ, this is the correct
     * way of doing things, although all other compilers
     * happily live without.
     */

    using Node::event_hook;
    using nest::Node::handles_test_event;
    using nest::Node::handle;

    nest::port send_test_event(nest::Node&, nest::rport, nest::synindex, bool);

    void handle(nest::CurrentEvent &);
    void handle(nest::DataLoggingRequest &); 
    
    nest::port handles_test_event(nest::CurrentEvent &, nest::rport);
    nest::port handles_test_event(nest::DataLoggingRequest &, nest::rport);

    
    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &);
    
  private:
    void init_state_(const Node& proto);
    void init_buffers_();
    void calibrate();
    
    void update(nest::Time const &, const long, const long);
    void event_hook( nest::DSSpikeEvent&);

    // END Boilerplate function declarations ----------------------------

    // Friends --------------------------------------------------------
    // The next two classes need to be friends to access the State_ class/member
    friend class nest::RecordablesMap<rbf_poisson_generator>;
    friend class nest::UniversalDataLogger<rbf_poisson_generator>;

  private:

    /**
      * Store independent parameters of the model.
      */
    struct Parameters_{
      double min_rate_rbf_;
      double max_rate_rbf_;
      double mean_current_rbf_;
      double sigma_current_rbf_;
      
      Parameters_(); //!< Sets default parameter values

      void get( DictionaryDatum& ) const; //!< Store current values in dictionary
      void set( const DictionaryDatum& ); //!< Set values from dicitonary
    };

  public:
    // ---------------------------------------------------------------- 

    /**
     * State variables of the model.
     * @note Copy constructor and assignment operator required because
     *       of C-style array.
     */
    struct State_ {

      double rate_; 

      State_(const Parameters_&);  //!< Default initialization
      State_(const State_&);
      
      void get(DictionaryDatum&) const;
      void set(const DictionaryDatum&, const Parameters_&);
    };    

    // ---------------------------------------------------------------- 

  private:

    /**
     * Buffers of the model.
     */
    struct Buffers_ {
      Buffers_(rbf_poisson_generator&);                   //!<Sets buffer pointers to 0
      Buffers_(const Buffers_&, rbf_poisson_generator&);  //!<Sets buffer pointers to 0

      //! Logger for all analog data
      nest::UniversalDataLogger<rbf_poisson_generator> logger_;

      /** buffers and sums up incoming spikes/currents */
      nest::RingBuffer currents_;

      // IntergrationStep_ should be reset with the neuron on ResetNetwork,
      // but remain unchanged during calibration. Since it is initialized with
      // step_, and the resolution cannot change after nodes have been created,
      // it is safe to place both here.
      double step_;           //!< step size in ms
    };

  // ------------------------------------------------------------

    struct Variables_
    {
      librandom::PoissonRandomDev poisson_dev_; //!< Random deviate generator
    };

    // Access functions for UniversalDataLogger -------------------------------
    
    //! Read out state vector elements, used by UniversalDataLogger
    double get_rate_() const { return S_.rate_; }

  // ------------------------------------------------------------

    Parameters_ P_;
    Variables_ V_;
    State_ S_;
    Buffers_ B_;

    //! Mapping of recordables names to access functions
    static nest::RecordablesMap<rbf_poisson_generator> recordablesMap_;

  };

  
  inline
  nest::port rbf_poisson_generator::send_test_event(nest::Node& target, nest::rport receptor_type, nest::synindex syn_id, bool dummy_target)
  {
  	if ( dummy_target )
    {
      nest::DSSpikeEvent e;
      e.set_sender( *this );
      return target.handles_test_event( e, receptor_type );
    }
    else
    {
      nest::SpikeEvent e;
      e.set_sender( *this );
      return target.handles_test_event( e, receptor_type );
    }
  }

  inline
  nest::port rbf_poisson_generator::handles_test_event(nest::CurrentEvent&, nest::rport receptor_type)
  {
    if (receptor_type != 0)
      throw nest::UnknownReceptorType(receptor_type, get_name());
    return 0;
  }

  inline
  nest::port rbf_poisson_generator::handles_test_event(nest::DataLoggingRequest& dlr, nest::rport receptor_type)
  {
    if (receptor_type != 0)
      throw nest::UnknownReceptorType(receptor_type, get_name());
    return B_.logger_.connect_logging_device(dlr, recordablesMap_);
  }

  inline
  void rbf_poisson_generator::get_status(DictionaryDatum &d) const
  {
    P_.get(d);
    S_.get(d);
  }

  inline
  void rbf_poisson_generator::set_status(const DictionaryDatum &d)
  {
    Parameters_ ptmp = P_;  // temporary copy in case of errors
    ptmp.set(d);                       // throws if BadProperty
    State_      stmp = S_;  // temporary copy in case of errors
    stmp.set(d, ptmp); 
    
    // if we get here, temporaries contain consistent set of properties
    P_ = ptmp;
    S_ = stmp;
    
  }
  
} // namespace

#endif //rbf_poisson_generator_H
