import nest
import numpy
import matplotlib.pylab as pylab

nest.Install('cerebellummodule')

nest.SetKernelStatus({'local_num_threads': 1})

num_neuron_pre = 1000
num_neuron_post = 1

SpGeneratorPF = nest.Create('spike_generator', num_neuron_pre)

NeuronPF = nest.Create('parrot_neuron', num_neuron_pre)
NeuronPC = nest.Create('iaf_cond_exp_cs', num_neuron_post)

PCvoltmeter = nest.Create("voltmeter",1,{"withgid": True, "withtime": True})

nest.Connect(PCvoltmeter, NeuronPC)

SpGeneratorCF = nest.Create('spike_generator', num_neuron_post)
NeuronCF = nest.Create('parrot_neuron', num_neuron_post)

SpDetector = nest.Create('spike_detector', 1)

nest.Connect(SpGeneratorPF, NeuronPF, 'one_to_one')
nest.Connect(SpGeneratorPF+NeuronPF+NeuronPC+SpGeneratorCF+NeuronCF, SpDetector, 'all_to_all')

# Frequency of presynaptic spikes (one per synapsis)
spike_freq_pre = 1000.0 
spike_times_pre = numpy.arange(400.0, 400.0+(1000.0/spike_freq_pre)*(num_neuron_pre), 1000./spike_freq_pre)
test_spike_time = 1500.0

CF_spike_time = [1000.0]
nest.SetStatus(SpGeneratorCF, {'spike_times': CF_spike_time})
for id_neuron,neuron in enumerate(SpGeneratorPF):
	nest.SetStatus([neuron], {'spike_times': [spike_times_pre[id_neuron]]+[test_spike_time]})

nest.Connect(SpGeneratorCF, NeuronCF, 'one_to_one')
PCReceptor = {'Normal': 0, 'CS_receptor' : 1}
syn_dict_CFPC = {'model': 'static_synapse', 'weight': 10.0, 'delay':1.0, 'receptor_type':PCReceptor['CS_receptor']}
nest.Connect(NeuronCF, NeuronPC, 'one_to_one', syn_spec=syn_dict_CFPC)
syn_dict_PFPC = {'model': 'stdp_sin_synapse', 'weight': 1.0, 'delay':1.0, 'receptor_type':PCReceptor['Normal'], 
				'A_plus': 0.05, 'A_minus': 0.2, 'Wmin':0.00, 'Wmax':2.00,
				'exponent': 20, 'peak': 100.0}
nest.Connect(NeuronPF, NeuronPC, syn_spec=syn_dict_PFPC)

connections = nest.GetConnections(source=NeuronPF, target=NeuronPC)
source = numpy.array(nest.GetStatus(connections, "source"))
weight_before = numpy.array(nest.GetStatus(connections, "weight"))

nest.Simulate(1550)

# Disable the plasticity and apply 
#nest.SetStatus(connections, {'A_plus': 0.00, 'A_minus': 0.00})

weight_after = numpy.array(nest.GetStatus(connections, "weight"))

neuron_act = nest.GetStatus(SpDetector, 'events') [0]
sp_id = neuron_act['senders']
sp_time = neuron_act['times']

pylab.figure()
pylab.plot(spike_times_pre, weight_after-weight_before)
pylab.ylim(min(weight_after-weight_before)-0.05,max(weight_after-weight_before)+0.05)
pylab.xlabel('PF spike time (ms)')
pylab.ylabel('Weight Diff. (nS)')


pylab.figure()
pylab.vlines(sp_time, sp_id-0.5, sp_id+0.5)
pylab.xlabel('Time (ms)')
pylab.ylabel('Neuron')

pylab.show()
