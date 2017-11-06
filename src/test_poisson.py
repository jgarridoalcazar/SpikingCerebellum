import nest
import time
import numpy
import matplotlib.pylab as plt

nest.set_verbosity('M_WARNING')

nest.Install('cerebellummodule')

num_neurons = 200

cur_generator = nest.Create('dc_generator', 1)

pop_poisson = nest.Create('cd_poisson_generator', num_neurons, params=
													{'min_rate': 3.0,
													'max_rate': 7.0,
													'min_current': -1.0,
													'max_current': 10.0})

pop_parrot = nest.Create('parrot_neuron', num_neurons)

spike_detector = nest.Create('spike_detector')

nest.Connect(cur_generator, pop_poisson, 'all_to_all')
nest.Connect(pop_poisson, pop_parrot, 'one_to_one')
nest.Connect(pop_poisson, spike_detector)

#nest.SetStatus(pop_poisson,{'rate': 10.0})
nest.SetStatus(cur_generator, {'amplitude': 5.0})

sim_time = 10000.0

sim_step = 2.0

el_ass_time = 0.0
el_sim_time = 0.0

for cur_time in numpy.arange(0.0, sim_time, sim_step):
	bef_ass = time.time()
	#nest.SetStatus(pop_poisson,{'rate': 10.0})
	for neu in pop_poisson:
		#nest.SetStatus([neu], {'rate': 10.0})
		pass
	bef_sim = time.time()
	nest.Simulate(sim_step)
	aft_sim = time.time()
	el_ass_time += (bef_sim-bef_ass)
	el_sim_time += (aft_sim-bef_sim)

neuron_act = nest.GetStatus(spike_detector, 'events') [0]
sp_id = neuron_act['senders']
sp_time = neuron_act['times']

print 'Simulation time:',el_sim_time,'. Assignement time:',el_ass_time
print 'Average frequency:',sp_id.size/(sim_time*num_neurons)*1.e3,'Hz'

plt.figure()
plt.plot(sp_time, sp_id, '.')
plt.show(block=True)
