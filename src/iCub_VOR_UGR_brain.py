# -*- coding: utf-8 -*-
"""
This file contains the setup of the neuronal network running the Mouse experiment with neuronal image recognition
"""
# pragma: no cover

__author__ = 'Lazar Mateev, Georg Hinkel, Alina Roitberg'

### The following can be removed when PyNN 0.8 has been established or we have a more elegant
### solution
import pyNN.nest as sim
import nest
import numpy as np
import logging

logger = logging.getLogger(__name__)

# Synapsis parameters
gc_pc_weights = 1.35
mf_vn_weights = 0.005
pc_vn_weights = -0.002
io_pc_weights = 0.0

input_weights = 0.0005

# Network parameters
num_MF_neurons = 100
num_GC_neurons = 2000
num_PC_neurons = 200
num_VN_neurons = 200
num_IO_neurons = 200



def create_brain():
    """
    Initializes PyNN with the neuronal network that has to be simulated for the mouse experiment
    """

    # Following parameters were taken from the husky braitenberg brain experiment (braitenberg.py)

    SENSORPARAMS = {'cm': 0.25,
                    'v_rest': -60.5,
                    'tau_m': 10.0,
                    'e_rev_E': 0.0,
                    'e_rev_I': -75.0,
                    'v_reset': -60.5,
                    'v_thresh': -60.0,
                    'tau_refrac': 1.0,
                    'tau_syn_E': 0.5,
                    'tau_syn_I': 2.5}
    SENSORPARAMS2 = {'cm': 0.025,
                    'v_rest': -60.5,
                    'tau_m': 10.0,
                    'e_rev_E': 0.0,
                    'e_rev_I': -75.0,
                    'v_reset': -60.5,
                    'v_thresh': -60.0,
                    'tau_refrac': 10.0,
                    'tau_syn_E': 2.5,
                    'tau_syn_I': 2.5}  

    OLD_PC_PARAMS = {'cm': 0.314,
              'v_rest': -70.0,
              'tau_m': 26.17,
              'e_rev_E': 0.0,
              'e_rev_I': -75.0,
              'e_rev_CS': 0.0,
              'v_reset': -70.0,
              'v_thresh': -52.0,
              'tau_refrac': 2.0,
              'tau_syn_E': 0.85,
              'tau_syn_I': 5.45,
              'tau_syn_CS': 0.85}

    PC_PARAMS = {'C_m': 0.314,
    		  'g_L': 8.2,
              'E_L': -70.0,
              'E_ex': 0.0,
              'E_in': -75.0,
              'e_cs': 0.0,
              'V_reset': -70.0,
              'V_th': -52.0,
              't_ref': 1.0,
              'tau_syn_ex': 0.85,
              'tau_syn_in': 5.45,
              'tau_syn_cs': 0.85}
    
    VN_PARAMS = {'cm': 0.002,
              'v_rest': -70.0,
              'tau_m': 10.0,
              'e_rev_E': 0.0,
              'e_rev_I': -80.0,
              'v_reset': -70.5,
              'v_thresh': -40.0,
              'tau_refrac': 1.0,
              'tau_syn_E': 0.5,
              'tau_syn_I': 14.0}

    SYNAPSE_PARAMS = {'U': 1.0,
                      'tau_rec': 0.1,
                      'tau_facil': 0.1,
                      'weight': 1.5e-4,
                      'delay': 0.1}
    
    try:
      nest.Install('cerebellummodule')
    except nest.NESTError:
        pass        
    

    parrot_neuron = sim.native_cell_type('parrot_neuron')
    # Create MF population 
    MF_population = sim.Population(num_MF_neurons,parrot_neuron,{},label='MFLayer')
    
    # Create GrC population
    GC_population = sim.Population(num_GC_neurons,parrot_neuron,{},label='GCLayer')

    # Create IO population
    IO_neuron = sim.native_cell_type('cd_poisson_generator')
    IO_neuron.injectable = True
    IO_neuron.always_local = True
    IO_population = sim.Population(num_IO_neurons,IO_neuron(min_rate = 1.0, max_rate = 10.0, min_current = 0.0, max_current = 1.0),label='IOLayer')
    #IO_population = sim.Population(num_IO_neurons,pc_neuron(**PC_PARAMS),label='IOLayer')

    
    # Create PC population
    pc_neuron = sim.native_cell_type('iaf_cond_exp_cs')
    PC_population = sim.Population(num_PC_neurons,pc_neuron(**PC_PARAMS),label='PCLayer')
    
    # Create DCN population
    VN_population = sim.Population(num_VN_neurons,sim.IF_cond_exp(**VN_PARAMS),label='VNLayer')
    
    
    
    int_num_IO_agonist_neurons = int(num_IO_neurons/2)
    float_num_IO_agonist_neurons = float(int_num_IO_agonist_neurons)
    increment = 0.2
    for x in range(int_num_IO_agonist_neurons):
    	IO_population[x].min_current = x/float_num_IO_agonist_neurons - increment
    	IO_population[x].max_current = x/float_num_IO_agonist_neurons + increment
    	IO_population[x + int_num_IO_agonist_neurons].min_current = x/float_num_IO_agonist_neurons - increment
    	IO_population[x + int_num_IO_agonist_neurons].max_current = x/float_num_IO_agonist_neurons + increment




    Input_MF_population = sim.Population(1,parrot_neuron,{},label='InputMFLayer')
    Input_MF_connections = []
    
    for time in range(2, 1002, 2):
      neuron_id = ((np.sin(time*2*np.pi/1000)+1.0)/2.0*(num_MF_neurons-4)).astype(int)
      Input_MF_con = sim.Projection(Input_MF_population,
                      sim.PopulationView(MF_population, range(neuron_id, neuron_id+4)),
                                      sim.AllToAllConnector(),
                                      sim.StaticSynapse(delay=time, weight=input_weights))
      Input_MF_connections.append(Input_MF_con)


    Input_GC_population = sim.Population(1,parrot_neuron,{},label='InputGCLayer')
    Input_GC_connections = []
    
    RELATIVE_AMPLITUDE = 1.0

    for time in range(2, 1002, 2):
      neuron_id = ((np.sin(time*2*np.pi/1000)+1.0)*RELATIVE_AMPLITUDE/2.0*(num_GC_neurons/2.0-2)).astype(int)
      Input_GC_con = sim.Projection(Input_GC_population,
                      sim.PopulationView(GC_population, range(neuron_id, neuron_id+2)),
                                      sim.AllToAllConnector(),
                                      sim.StaticSynapse(delay=time, weight=input_weights))
      Input_GC_connections.append(Input_GC_con)

      neuron_id = (num_GC_neurons/2 + ((np.cos(time*2*np.pi/1000)+1.0)*RELATIVE_AMPLITUDE/2.0*(num_GC_neurons/2-2))).astype(int)
      Input_GC_con = sim.Projection(Input_GC_population,
                      sim.PopulationView(GC_population, range(neuron_id, neuron_id+2)),
                                      sim.AllToAllConnector(),
                                      sim.StaticSynapse(delay=time, weight=input_weights))
      Input_GC_connections.append(Input_GC_con)

    #neuron_id = 0
    #for time in range(2, 1002, 2):
    #  Input_GC_con = sim.Projection(Input_GC_population,
    #                  				  sim.PopulationView(GC_population, range(neuron_id, neuron_id+4)),
    #                                  sim.AllToAllConnector(),
    #                                  sim.StaticSynapse(delay=time, weight=input_weights))
    #  Input_GC_connections.append(Input_GC_con)
    #  neuron_id +=4

    
  
    
    # Create MF-VN connections
    mf_dcn_connections = sim.Projection(MF_population,
    VN_population,
    sim.AllToAllConnector(),
    sim.StaticSynapse(delay=1.0,weight=mf_vn_weights))
    
    # Create PC-VN connections
    pc_vn_connections = sim.Projection(PC_population,
                                       VN_population,
                                       sim.OneToOneConnector(),
                                       sim.StaticSynapse(delay=1.0, weight=pc_vn_weights))


    # Create IO-PC connections
    io_pc_connections = sim.Projection(IO_population,
                                       PC_population,
                                       sim.OneToOneConnector(),
                                       receptor_type='COMPLEX_SPIKE',
                                       synapse_type = sim.StaticSynapse(delay=1.0, weight=io_pc_weights))   
    
    
    
    # Create GC-PC connections
    stdp_syn = sim.native_synapse_type('stdp_sin_synapse')(**{'weight':gc_pc_weights,
                                    'delay':1.0,
                                                                'exponent': 20,
                                                                'peak': 100.0,
                                                                'A_plus': 0.01,
                                                                'A_minus': 0.06,
                                                                'Wmin': 0.0,
                                                                'Wmax': 3.5})
      
    gc_pc_connections = sim.Projection(GC_population,
                                       PC_population,
                                       sim.AllToAllConnector(),
                                       receptor_type='AMPA',
                                       synapse_type = stdp_syn)
  
    # Group all neural layers
    population = MF_population + GC_population + PC_population + VN_population + IO_population + Input_MF_population + Input_GC_population
    
    # Set Vm to resting potential
    sim.initialize(PC_population, V_m=PC_population.get('E_L'))
    sim.initialize(VN_population, v=VN_population.get('v_rest'))
    
    return population

sim.setup(timestep=0.5, threads=1)
circuit = create_brain()
sim_time = 10000
sim_step = 2
import time

el_sim_time = 0
el_ass_time = 0

start = time.time()
for cur_time in range(0,sim_time,sim_step):
  time1 = time.time()
  #for neu in poisson:
    #neu.rate = np.random.randint(6,10)
    #neu.rate = 1.0
  # pass
  time2 = time.time()
  sim.run(sim_step)
  time3 = time.time()
  el_sim_time += (time3-time2)
  el_ass_time += (time2-time1)
end = time.time()
print 'Elapsed time:',end-start,'Assignement time:',el_ass_time,'Simulation time:',el_sim_time
