import pyNN.nest as sim
import nest
import numpy as np
import logging

sim.setup(timestep=0.5, threads=1, min_delay=0.5, max_delay=1000.0)

def create_brain():
    """
    Initializes PyNN with the neuronal network that has to be simulated for the mouse experiment
    """    
    try:
      nest.Install('cerebellummodule')
    except nest.NESTError:
        pass        
    

    global dc_generator
    # Create GrC population
    #GC_population = sim.Population(100,sim.IF_cond_exp(),label='GCLayer')
    
    # Create PC population
    pc_neuron = sim.native_cell_type('iaf_cond_exp_cs')
    PC_population = sim.Population(200,pc_neuron(),label='PCLayer')
    
    # Create PC population
    dc_generator = sim.DCSource(amplitude=0.5, start=20.0, stop=80.0)
    
    # Create IO population
    CD_neuron = sim.native_cell_type('cd_poisson_generator')
    CD_neuron.injectable = True
    CD_neuron.always_local = True
    CD_population = sim.Population(200,CD_neuron(min_rate = 1.0, max_rate = 10.0, min_current = 0.0, max_current = 1.0),label='IOLayer')

    dc_generator.inject_into(CD_population)

    # Create IO-PC connections
    io_pc_connections = sim.Projection(CD_population,
                                       PC_population,
                                       sim.OneToOneConnector(),
                                       receptor_type='COMPLEX_SPIKE',
                                       synapse_type = sim.StaticSynapse(delay=1.0, weight=1.0))   
    
    
    
    # Create GC-PC connections
    stdp_syn = sim.native_synapse_type('stdp_sin_synapse')(**{'weight':1.0,
                                                                'delay':1.0,
                                                                'exponent': 20,
                                                                'peak': 100.0,
                                                                'A_plus': 0.01,
                                                                'A_minus': 0.06,
                                                                'Wmin': 0.0,
                                                                'Wmax': 3.5})
      
    #gc_pc_connections = sim.Projection(GC_population,
    #                                   PC_population,
    #                                   sim.AllToAllConnector(),
    #                                   receptor_type='AMPA',
    #                                   synapse_type = stdp_syn)

    # Group all neural layers
    population = CD_population 
    
    return population


circuit = create_brain()

sim_time = 10000
sim_step = 2000
import time

el_sim_time = 0
el_ass_time = 0

start = time.time()
for cur_time in range(0,sim_time,sim_step):
  time1 = time.time()
  #for neu in poisson:
  dc_generator.amplitude = np.random.rand(1)[0]
    #neu.rate = 1.0
  # pass
  time2 = time.time()
  sim.run(sim_step)
  time3 = time.time()
  el_sim_time += (time3-time2)
  el_ass_time += (time2-time1)
end = time.time()
print 'Elapsed time:',end-start,'Assignement time:',el_ass_time,'Simulation time:',el_sim_time

