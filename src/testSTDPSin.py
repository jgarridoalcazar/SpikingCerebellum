import nest
import numpy
import matplotlib.pylab as pylab

nest.Install('cerebellummodule')

nest.SetKernelStatus({'local_num_threads': 2})

num_neuron_pre = 1000
num_neuron_post = 1

SpGeneratorPre = nest.Create('spike_generator', num_neuron_pre)

NeuronPre = nest.Create('parrot_neuron', num_neuron_pre)
NeuronPost = nest.Create('iaf_neuron', num_neuron_post)

SpGeneratorVol = nest.Create('spike_generator', num_neuron_post)
Modulator = nest.Create('parrot_neuron', num_neuron_post)

Volumes = nest.Create('volume_transmitter_sin', num_neuron_post)
nest.SetStatus(Volumes,{'deliver_interval': 1})

SpDetector = nest.Create('spike_detector', 1)

nest.Connect(SpGeneratorPre, NeuronPre, 'one_to_one')
nest.Connect(SpGeneratorPre+NeuronPre+NeuronPost+SpGeneratorVol+Modulator, SpDetector, 'all_to_all')

# Frequency of presynaptic spikes (one per synapsis)
spike_freq_pre = 1000.0 
spike_times_pre = numpy.arange(400.0, 400.0+(1000.0/spike_freq_pre)*(num_neuron_pre), 1000./spike_freq_pre)

modulator_spike_time = [1000.0]
nest.SetStatus(SpGeneratorVol, {'spike_times': modulator_spike_time})
for id_neuron,neuron in enumerate(SpGeneratorPre):
	nest.SetStatus([neuron], {'spike_times': [spike_times_pre[id_neuron]]})

for neuron in range(num_neuron_post):
	synapse_name = 'STDP_' + str(neuron)
	nest.CopyModel('stdp_sin_synapse',synapse_name, {'weight': 1.0,
													'exponent': 20,
													'peak': 100.0})
	nest.SetDefaults(synapse_name, {'vt' : Volumes[neuron]})
	nest.SetDefaults(synapse_name, {'A_plus' : 0.05,
									'A_minus' : 0.2,
									'Wmin' : 0.00,
									'Wmax' : 2.00})
	nest.Connect([SpGeneratorVol[neuron]],[Modulator[neuron]])
	nest.Connect([Modulator[neuron]],[Volumes[neuron]])
	nest.Connect(NeuronPre,[NeuronPost[neuron]],'all_to_all',synapse_name)

connections = nest.GetConnections(source=NeuronPre, target=NeuronPost)
source = numpy.array(nest.GetStatus(connections, "source"))
weight_before = numpy.array(nest.GetStatus(connections, "weight"))

nest.Simulate(2000)

weight_after = numpy.array(nest.GetStatus(connections, "weight"))

neuron_act = nest.GetStatus(SpDetector, 'events') [0]
sp_id = neuron_act['senders']
sp_time = neuron_act['times']

pylab.figure()
pylab.plot(spike_times_pre, weight_after-weight_before)
pylab.xlabel('PF spike time (ms)')
pylab.ylabel('Diff. Weight (nS)')


pylab.figure()
pylab.vlines(sp_time, sp_id-0.5, sp_id+0.5)
pylab.xlabel('Time (ms)')
pylab.ylabel('Neuron')
pylab.show()


