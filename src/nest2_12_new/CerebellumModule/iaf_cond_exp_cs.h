/*
 *  iaf_cond_exp_cs.h
 *
 *  This file is based on the iaf_cond_exp cell model distributed with NEST.
 *  
 *  Modified by: Jesús Garrido (jgarridoalcazar at gmail.com) in 2017.
 */

#ifndef IAF_COND_EXP_CS_H
#define IAF_COND_EXP_CS_H

#include "config.h"

#ifdef HAVE_GSL

#include "nest.h"
#include "event.h"
#include "archiving_node_cs.h"
#include "ring_buffer.h"
#include "connection.h"
#include "universal_data_logger.h"
#include "recordables_map.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

/* BeginDocumentation
Name: iaf_cond_exp_cs - Conductance based leaky integrate-and-fire neuron model with complex spike.

Description:
iaf_cond_exp_cs is an implementation of a spiking neuron using IAF dynamics with
conductance-based synapses and adaptive threshold. Incoming spike events induce a post-synaptic change
of conductance modelled by an exponential function. The exponential function 
is normalised such that an event of weight 1.0 results in a peak conductance of 1 nS.
The complex spike dynamics is implemented in a diferent port (number 1) and follows its own dynamics (tau
and E_rev constants). This neuron model includes a CS register that can be used for plasticity rules.

Parameters: 
The following parameters can be set in the status dictionary.

V_m        double - Membrane potential in mV 
E_L        double - Leak reversal potential in mV.
th_C       double - Threshold constant. Increment of the threshold voltage after each spike in mV.
t_ref      double - Duration of refractory period in ms. 
V_th       double - Spike threshold in mV.
V_reset    double - Reset potential of the membrane in mV.
E_ex       double - Excitatory reversal potential in mV.
E_in       double - Inhibitory reversal potential in mV.
E_cs       double - Complex spike reversal potential in mV.
C_m        double - Capacity of the membrane in pF
g_L        double - Leak conductance in nS;
tau_th	   double - Time constant of the threshold adaptation in ms.
tau_syn_ex double - Time constant of the excitatory synaptic exponential function in ms.
tau_syn_in double - Time constant of the inhibitory synaptic exponential function in ms.
tau_syn_cs double - Time constant of the complex spike synaptic exponential function in ms.
I_e        double - Constant external input current in pA.

Sends: SpikeEvent

Receives: SpikeEvent, CurrentEvent, DataLoggingRequest

References: 

Author: Jesus Garrido

SeeAlso: iaf_cond_exp
*/

// Define name constants for state variables and parameters
namespace nest
{
	namespace names
	{
    	// Neuron parameters
    	extern const Name tau_syn_cs;  //!<  Time constant of the complex spike synaptic exponential function in ms.
    	extern const Name E_cs;        //!<  Complex spike reversal potential in mV.
      extern const Name g_cs;
      extern const Name COMPLEX_SPIKE;
      extern const Name GABA;
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
  extern "C"
  int iaf_cond_exp_cs_dynamics (double, const double*, double*, void*);
  
  class iaf_cond_exp_cs : public mynest::Archiving_Node_CS
  {
    
  public:        
    
    iaf_cond_exp_cs();
    iaf_cond_exp_cs(const iaf_cond_exp_cs&);
    ~iaf_cond_exp_cs();

    /**
     * Import sets of overloaded virtual functions.
     * We need to explicitly include sets of overloaded
     * virtual functions into the current scope.
     * According to the SUN C++ FAQ, this is the correct
     * way of doing things, although all other compilers
     * happily live without.
     */

    using nest::Node::handles_test_event;
    using nest::Node::handle;

    nest::port send_test_event(nest::Node&, nest::rport, nest::synindex, bool);
    
    void handle(nest::SpikeEvent &);
    void handle(nest::CurrentEvent &);
    void handle(nest::DataLoggingRequest &); 
    
    nest::port handles_test_event(nest::SpikeEvent &, nest::rport);
    nest::port handles_test_event(nest::CurrentEvent &, nest::rport);
    nest::port handles_test_event(nest::DataLoggingRequest &, nest::rport);
    
    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &);
    
  private:
    void init_state_(const Node& proto);
    void init_buffers_();
    void calibrate();
    void update(nest::Time const &, const long, const long);

    // END Boilerplate function declarations ----------------------------

    // Friends --------------------------------------------------------

    // make dynamics function quasi-member
    friend int iaf_cond_exp_cs_dynamics(double, const double*, double*, void*);

    // The next two classes need to be friends to access the State_ class/member
    friend class nest::RecordablesMap<iaf_cond_exp_cs>;
    friend class nest::UniversalDataLogger<iaf_cond_exp_cs>;

  private:

    // ---------------------------------------------------------------- 
    enum SynapseTypes
  {
    INF_SPIKE_RECEPTOR = 0,
    AMPA,
    GABA,
    COMPLEX_SPIKE,
    SUP_SPIKE_RECEPTOR
  };

    //! Model parameters
	struct Parameters_ {
	  double V_reset_;    //!< Reset Potential in mV
	  double V_th_;    //!< Reset Potential in mV
    double t_ref_;      //!< Refractory period in ms
	  double g_L;      //!< Leak Conductance in nS
	  double C_m;      //!< Membrane Capacitance in pF
	  double E_ex;        //!< Excitatory reversal Potential in mV
	  double E_in;        //!< Inhibitory reversal Potential in mV
    double E_cs;        //!< CS reversal Potential in mV
	  double E_L;         //!< Leak reversal Potential (aka resting potential) in mV
	  double tau_synE;    //!< Synaptic Time Constant Excitatory Synapse in ms
	  double tau_synI;    //!< Synaptic Time Constant for Inhibitory Synapse in ms
    double tau_synCS;   //!< Synaptic Time Constant for CS Synapse in ms
	  double I_e;         //!< Constant Current in pA
	  
	  Parameters_();  //!< Sets default parameter values

	  void get(DictionaryDatum&) const;  //!< Store current values in dictionary
	  void set(const DictionaryDatum&);  //!< Set values from dicitonary
	};

  public:
    // ---------------------------------------------------------------- 

    /**
     * State variables of the model.
     * @note Copy constructor and assignment operator required because
     *       of C-style array.
     */
    struct State_ {

    	//! Symbolic indices to the elements of the state vector y
	  enum StateVecElems { V_M = 0,
			   G_EXC,     
			   G_INH,
         G_CS,
			   STATE_VEC_SIZE };

      double y_[STATE_VEC_SIZE];  //!< neuron state, must be C-array for GSL solver
      int    r_;                  //!< number of refractory steps remaining

      State_(const Parameters_&);  //!< Default initialization
      State_(const State_&);
      State_& operator=(const State_&);

      void get(DictionaryDatum&) const;
      void set(const DictionaryDatum&, const Parameters_&);
    };    

    // ---------------------------------------------------------------- 

  private:
    /**
     * Buffers of the model.
     */
    struct Buffers_ {
      Buffers_(iaf_cond_exp_cs&);                   //!<Sets buffer pointers to 0
      Buffers_(const Buffers_&, iaf_cond_exp_cs&);  //!<Sets buffer pointers to 0

      //! Logger for all analog data
      nest::UniversalDataLogger<iaf_cond_exp_cs> logger_;

      /** buffers and sums up incoming spikes/currents */
      nest::RingBuffer spike_exc_;
      nest::RingBuffer spike_inh_;
      nest::RingBuffer spike_cs_;
      nest::RingBuffer currents_;

      /** GSL ODE stuff */
      gsl_odeiv_step*    s_;    //!< stepping function
      gsl_odeiv_control* c_;    //!< adaptive stepsize control function
      gsl_odeiv_evolve*  e_;    //!< evolution function
      gsl_odeiv_system   sys_;  //!< struct describing system
      
      // IntergrationStep_ should be reset with the neuron on ResetNetwork,
      // but remain unchanged during calibration. Since it is initialized with
      // step_, and the resolution cannot change after nodes have been created,
      // it is safe to place both here.
      double step_;           //!< step size in ms
      double   IntegrationStep_;//!< current integration time step, updated by GSL

      /** 
       * Input current injected by CurrentEvent.
       * This variable is used to transport the current applied into the
       * _dynamics function computing the derivative of the state vector.
       * It must be a part of Buffers_, since it is initialized once before
       * the first simulation, but not modified before later Simulate calls.
       */
      double I_stim_;
    };

     // ---------------------------------------------------------------- 

     /**
      * Internal variables of the model.
      */
     struct Variables_ { 
    	int    RefractoryCounts_;
     };

    // Access functions for UniversalDataLogger -------------------------------
    
    //! Read out state vector elements, used by UniversalDataLogger
    template <State_::StateVecElems elem>
    double get_y_elem_() const { return S_.y_[elem]; }

    // ---------------------------------------------------------------- 

    Parameters_ P_;
    State_      S_;
    Variables_  V_;
    Buffers_    B_;

    //! Mapping of recordables names to access functions
    static nest::RecordablesMap<iaf_cond_exp_cs> recordablesMap_;
  };

  
  inline
  nest::port iaf_cond_exp_cs::send_test_event(nest::Node& target, nest::rport receptor_type, nest::synindex, bool)
  {
	nest::SpikeEvent e;
    e.set_sender(*this);
    return target.handles_test_event(e, receptor_type);
  }

  inline
  nest::port iaf_cond_exp_cs::handles_test_event(nest::SpikeEvent&, nest::rport receptor_type)
  {
    if (not( INF_SPIKE_RECEPTOR < receptor_type
         && receptor_type < SUP_SPIKE_RECEPTOR ))
      throw nest::UnknownReceptorType(receptor_type, get_name());
    return receptor_type;
  }
 
  inline
  nest::port iaf_cond_exp_cs::handles_test_event(nest::CurrentEvent&, nest::rport receptor_type)
  {
    if (receptor_type != 0)
      throw nest::UnknownReceptorType(receptor_type, get_name());
    return 0;
  }

  inline
  nest::port iaf_cond_exp_cs::handles_test_event(nest::DataLoggingRequest& dlr, nest::rport receptor_type)
  {
    if (receptor_type != 0)
      throw nest::UnknownReceptorType(receptor_type, get_name());
    return B_.logger_.connect_logging_device(dlr, recordablesMap_);
  }
 
  inline
  void iaf_cond_exp_cs::get_status(DictionaryDatum &d) const
  {
    P_.get(d);
    S_.get(d);
    Archiving_Node_CS::get_status(d);

    DictionaryDatum receptor_type = new Dictionary();

    ( *receptor_type )[ nest::names::AMPA ] = AMPA;
    ( *receptor_type )[ nest::names::GABA ] = GABA;
    ( *receptor_type )[ nest::names::COMPLEX_SPIKE ] = COMPLEX_SPIKE;
    
    ( *d )[ nest::names::receptor_types ] = receptor_type;

    (*d)[nest::names::recordables] = recordablesMap_.get_list();
  }

  inline
  void iaf_cond_exp_cs::set_status(const DictionaryDatum &d)
  {
    Parameters_ ptmp = P_;  // temporary copy in case of errors
    ptmp.set(d);                       // throws if BadProperty
    State_      stmp = S_;  // temporary copy in case of errors
    stmp.set(d, ptmp);                 // throws if BadProperty

    // We now know that (ptmp, stmp) are consistent. We do not 
    // write them back to (P_, S_) before we are also sure that 
    // the properties to be set in the parent class are internally 
    // consistent.
    Archiving_Node_CS::set_status(d);

    // if we get here, temporaries contain consistent set of properties
    P_ = ptmp;
    S_ = stmp;
  }
  
} // namespace

#endif //HAVE_GSL
#endif //IAF_COND_EXP_CS_H
