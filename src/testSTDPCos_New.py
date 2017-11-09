import nest
import numpy
import matplotlib.pylab as pylab

nest.Install('cerebellummodule')

nest.SetKernelStatus({'local_num_threads': 1})

num_neuron_pre = 1000
num_neuron_post = 1

SpGeneratorMF = nest.Create('spike_generator', num_neuron_pre)

NeuronMF = nest.Create('parrot_neuron', num_neuron_pre)
NeuronDCN = nest.Create('iaf_cond_exp_cos', num_neuron_post,  
						params={"tau_cos": 100.0,
								"exponent": 2.0})

DCNvoltmeter = nest.Create("voltmeter",1,{"withgid": True, "withtime": True})

nest.Connect(DCNvoltmeter, NeuronDCN)

SpGeneratorPC = nest.Create('spike_generator', num_neuron_post)
NeuronPC = nest.Create('parrot_neuron', num_neuron_post)

SpDetector = nest.Create('spike_detector', 1)

nest.Connect(SpGeneratorMF, NeuronMF, 'one_to_one')
nest.Connect(SpGeneratorMF+NeuronMF+NeuronDCN+SpGeneratorPC+NeuronPC, SpDetector, 'all_to_all')

# Frequency of presynaptic spikes (one per synapsis)
spike_freq_pre = 1000.0 
spike_times_pre = numpy.arange(400.0, 400.0+(1000.0/spike_freq_pre)*(num_neuron_pre), 1000./spike_freq_pre)
test_spike_time = 1500.0

PC_spike_time = [1000.0]
nest.SetStatus(SpGeneratorPC, {'spike_times': PC_spike_time})
for id_neuron,neuron in enumerate(SpGeneratorMF):
	nest.SetStatus([neuron], {'spike_times': [spike_times_pre[id_neuron]]+[test_spike_time]})

nest.Connect(SpGeneratorPC, NeuronPC, 'one_to_one')
DCNReceptor = {'AMPA': 1, 'GABA': 2, 'TEACHING_SIGNAL' : 3}
syn_dict_PCDCN = {'model': 'static_synapse', 'weight': 10.0, 'delay':1.0, 'receptor_type':DCNReceptor['TEACHING_SIGNAL']}
nest.Connect(NeuronPC, NeuronDCN, 'one_to_one', syn_spec=syn_dict_PCDCN)
syn_dict_MFDCN = {'model': 'stdp_cos_synapse', 'weight': 1.0, 'delay':1.0, 'receptor_type':DCNReceptor['AMPA'], 
				'A_plus': 0.05, 'A_minus': 0.2, 'Wmin':0.00, 'Wmax':2.00,
				'exponent': 2.0, 'tau_cos': 100.0}
nest.Connect(NeuronMF, NeuronDCN, syn_spec=syn_dict_MFDCN)

connections = nest.GetConnections(source=NeuronMF, target=NeuronDCN)
source = numpy.array(nest.GetStatus(connections, "source"))
weight_before = numpy.array(nest.GetStatus(connections, "weight"))

#nest.Simulate(1550)
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
pylab.xlabel('MF-DCN spike time (ms)')
pylab.ylabel('Weight Diff. (nS)')


pylab.figure()
pylab.vlines(sp_time, sp_id-0.5, sp_id+0.5)
pylab.xlabel('Time (ms)')
pylab.ylabel('Neuron')

pylab.show()
