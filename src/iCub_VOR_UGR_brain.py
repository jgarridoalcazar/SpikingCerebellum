# -*- coding: utf-8 -*-
"""
This file contains the setup of the neuronal network running the Mouse experiment with neuronal image recognition
"""
# pragma: no cover

__author__ = 'Lazar Mateev, Georg Hinkel, Alina Roitberg'

### The following can be removed when PyNN 0.8 has been established or we have a more elegant
### solution
import pyNN.nest as sim
import pyNN.random
import nest
import numpy as np
import logging

logger = logging.getLogger(__name__)

gc_pc_weights = 1.2
mf_vn_weights = 0.005 #cambiar
pc_vn_weights = 0.002 #cambiar
io_pc_weights = 0.0
mf_gc_weights = 0.00015
input_weights = 0.00025

# Network parameters
num_MF_neurons = 100
num_GC_neurons = 2000
num_PC_neurons = 1
num_VN_neurons = 1
num_IO_neurons = 1


sim.setup(timestep=0.5, threads=8)

def create_brain():
    """
    Initializes PyNN with the neuronal network that has to be simulated for the mouse experiment
    """

    # Following parameters were taken from the husky braitenberg brain experiment (braitenberg.py)

    GR_PARAMS = {'cm': 0.002,
                 'v_rest': -70.0,
                 'tau_m': 10.0,
                 'e_rev_E': 0.0,
                 'e_rev_I': -75.0,
                 'v_reset': -70.0,
                 'v_thresh': -40.0,
                 'tau_refrac': 10.0,
                 'tau_syn_E': 1.0,
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
    
    VN_PARAMS = {'C_m': 0.002,
              'g_L': 0.02, 
              'E_L': -70.0,
              'E_ex': 0.0,
              'E_in': -80.0,
              'e_ts': 0.0,
              'V_reset': -70.5,
              'V_th': -40.0,
              't_ref': 1.0,
              'tau_syn_ex': 0.5,
              'tau_syn_in': 14.0,
              'tau_syn_ts': 0.85,
              'tau_cos': 10.0,
              'exponent': 2.0}

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
    #MF_population = sim.Population(num_MF_neurons,parrot_neuron,{},label='MFLayer')
    pre_MF_neuron = sim.native_cell_type('rbf_poisson_generator')
    pre_MF_neuron.injectable = True
    #pre_MF_neuron.always_local = True
    pre_MF_population = sim.Population(num_MF_neurons,pre_MF_neuron(min_rate = 0.0, max_rate = 500.0, mean_current = 0.0, sigma_current = 0.01),label='inputMFLayer')


    num_pre_MF_pos = int(num_MF_neurons/2)
    for x in range(num_pre_MF_pos):
      pre_MF_population[x].mean_current = x/float(num_pre_MF_pos) + 0.01
      pre_MF_population[x + num_pre_MF_pos].mean_current = x/float(num_pre_MF_pos) + 0.01


    # Create MF population 
    MF_population = sim.Population(num_MF_neurons,parrot_neuron(),label='MFLayer')


    # Create input_MF-MF connections
    input_mf_mf_connections = sim.Projection(pre_MF_population,
                                       MF_population,
                                       sim.OneToOneConnector(),
                                       sim.StaticSynapse(delay=1.0, weight=0))



    
    
    # Create GrC population
    #GC_population = sim.Population(num_GC_neurons,parrot_neuron,{},label='GCLayer')
    GC_population = sim.Population(num_GC_neurons,sim.IF_cond_alpha(**GR_PARAMS),label='GCLayer')


    delay_distr = pyNN.random.RandomDistribution('uniform', (1.0, 10.0), rng=pyNN.random.NumpyRNG(seed=85524))
    weight_distr = pyNN.random.RandomDistribution('uniform', (mf_gc_weights * 0.5, mf_gc_weights * 1.5), rng=pyNN.random.NumpyRNG(seed=85524))


    float_num_MF_neurons = float (num_MF_neurons)
    for i in range (num_MF_neurons):
      GC_lower_index = int(round((i / float_num_MF_neurons) * num_GC_neurons)-40)
      GC_upper_index = GC_lower_index + 80
      if(GC_lower_index < 0):
        GC_lower_index = 0
        GC_upper_index = 80
        
      elif(GC_upper_index > num_GC_neurons):
        GC_lower_index = num_GC_neurons - 80
        GC_upper_index = num_GC_neurons

      MF_GC_connections = sim.Projection(sim.PopulationView(MF_population, range(i, i+1)),
                                      sim.PopulationView(GC_population, range(GC_lower_index, GC_upper_index)),
                                      sim.AllToAllConnector(),
                                      sim.StaticSynapse(delay=delay_distr, weight=weight_distr))





    
    # Create PC population
    pc_neuron = sim.native_cell_type('iaf_cond_exp_cs')
    PC_population = sim.Population(num_PC_neurons,pc_neuron(**PC_PARAMS),label='PCLayer')
    
    # Create DCN population
    vn_neuron = sim.native_cell_type('iaf_cond_exp_cos')
    VN_population = sim.Population(num_VN_neurons,vn_neuron(**VN_PARAMS),label='VNLayer')
    
    # Create IO population
    IO_population = sim.Population(num_IO_neurons,parrot_neuron,{},label='IOLayer')

    #Create pre IO population
    pre_IO_neuron = sim.native_cell_type('cd_poisson_generator')
    pre_IO_neuron.injectable = True
    #pre_IO_neuron.always_local = True
    pre_IO_population = sim.Population(num_IO_neurons,pre_IO_neuron(min_rate = 1.0, max_rate = 10.0, min_current = 0.0, max_current = 1.0),label='preIOLayer')

    
    int_num_IO_agonist_neurons = int(num_IO_neurons/2)
    float_num_IO_agonist_neurons = float(int_num_IO_agonist_neurons)
    increment = 1.0
    for x in range(int_num_IO_agonist_neurons):
      IO_population[x].min_current = x/float_num_IO_agonist_neurons 
      IO_population[x].max_current = x/float_num_IO_agonist_neurons + increment
      IO_population[x + int_num_IO_agonist_neurons].min_current = x/float_num_IO_agonist_neurons 
      IO_population[x + int_num_IO_agonist_neurons].max_current = x/float_num_IO_agonist_neurons + increment




  
  
    
    # Create MF-VN connections
    stdp_cos = sim.native_synapse_type('stdp_cos_synapse')(**{'weight':mf_vn_weights,
                                                              'delay':1.0,
                                                              'exponent': 2.0,
                                                              'tau_cos': 10.0,
                                                              'A_plus': 0.0,
                                                              'A_minus': 0.0,
                                                              'Wmin': 0.0,
                                                              'Wmax': 5.0})


    mf_vn_connections = sim.Projection(MF_population,
                                      VN_population,
                                      sim.AllToAllConnector(),
                                      receptor_type='AMPA',
                                      synapse_type = stdp_cos)
    
    # Create PC-VN connections
    pc_vn_connections = sim.Projection(PC_population,
                                       VN_population,
                                       sim.OneToOneConnector(),
                                       receptor_type='GABA',
                                       synapse_type = sim.StaticSynapse(delay=1.0, weight=pc_vn_weights))

    pc_vn_connections = sim.Projection(PC_population,
                                       VN_population,
                                       sim.OneToOneConnector(),
                                       receptor_type='TEACHING_SIGNAL',
                                       synapse_type = sim.StaticSynapse(delay=1.0, weight=0.0))


    # Create IO-PC connections
    io_pc_connections = sim.Projection(IO_population,
                                       PC_population,
                                       sim.OneToOneConnector(),
                                       receptor_type='COMPLEX_SPIKE',
                                       synapse_type = sim.StaticSynapse(delay=1.0, weight=io_pc_weights))   

    # Create pre IO-IO connections
    pre_io_io_connections = sim.Projection(pre_IO_population,
                                       IO_population,
                                       sim.OneToOneConnector(),
                                       synapse_type = sim.StaticSynapse(delay=10.0))  
    
    
    
    # Create GC-PC connections
    stdp_syn = sim.native_synapse_type('stdp_sin_synapse')(**{'weight':gc_pc_weights,
                                    'delay':1.0,
                                                                'exponent': 20,
                                                                'peak': 100.0,
                                                                'A_plus': 0.006,
                                                                'A_minus': 0.05,
                                                                'Wmin': 0.0,
                                                                'Wmax': 3.5})
      
    gc_pc_connections = sim.Projection(GC_population,
                                       PC_population,
                                       sim.AllToAllConnector(),
                                       receptor_type='AMPA',
                                       synapse_type = stdp_syn)
  
    # Group all neural layers
    population = MF_population + GC_population + PC_population + VN_population + IO_population + pre_MF_population + pre_IO_population
    
    # Set Vm to resting potential
    sim.initialize(PC_population, V_m=PC_population.get('E_L'))
    sim.initialize(VN_population, V_m=VN_population.get('E_L'))
    
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
