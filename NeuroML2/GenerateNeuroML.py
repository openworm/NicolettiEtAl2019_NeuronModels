from neuroml import NeuroMLDocument
from neuroml.utils import component_factory

from pyneuroml import pynml

c=1.2 

ek=-80
eca=60
eleak=-80
ena=30

colors = {'AWCon':'0 0 0.8', 'RMD':'0 0.8 0', 'GenericMuscleCell':'0.8 0 0'}

def create_channel_file(chan_id):

    chan_doc = NeuroMLDocument(id=chan_id, notes="A channel from Nicoletti et al. 2019")
    chan_fn = "%s.channel.nml"%chan_id

    channel = component_factory(
        "IonChannelHH", id=chan_id, conductance="10pS", notes="%s channel from Nicoletti et al. 2019"%chan_id
    )
    
    chan_doc.add(channel)
    chan_doc.validate(recursive=True)

    pynml.write_neuroml2_file(
        nml2_doc=chan_doc, nml2_file_name=chan_fn, validate=True
    )

    return chan_fn


def generate_nmllite(cell, duration=700, config='IClamp',parameters = None):

    from neuromllite import Cell, InputSource
    #from neuromllite.NetworkGenerator import *
    from neuromllite.utils import create_new_model
    
    reference = "%s_%s"%(config, cell)

    cell_id = '%s'%cell
    cell_nmll = Cell(id=cell_id, neuroml2_source_file='%s.cell.nml'%(cell))
 
    ################################################################################
    ###   Add some inputs
    
    if 'IClamp' in config:
        
        if not parameters:
            parameters = {}
            parameters['stim_amp'] = '10pA'
        
        input_source = InputSource(id='iclamp_0', 
                                   neuroml2_input='PulseGenerator', 
                                   parameters={'amplitude':'stim_amp', 'delay':'100ms', 'duration':'500ms'})
      
        
    else:

        if not parameters:
            parameters = {}
            parameters['average_rate'] = '100 Hz'
            parameters['number_per_cell'] = '10'
            
        input_source = InputSource(id='pfs0', 
                                   neuroml2_input='PoissonFiringSynapse', 
                                   parameters={'average_rate':'average_rate', 
                                               'synapse':syn_exc.id, 
                                               'spike_target':"./%s"%syn_exc.id})
                                               
    sim, net = create_new_model(reference,
                     duration, 
                     dt=0.025, # ms 
                     temperature=34, # degC
                     default_region='Worm',
                     parameters = parameters,
                     cell_for_default_population=cell_nmll,
                     color_for_default_population=colors[cell],
                     input_for_default_population=input_source)

    return sim, net




for cell_id in ['AWCon','RMD']:

    # Create the nml file and add the ion channels
    cell_doc = NeuroMLDocument(id=cell_id, notes="A cell from Nicoletti et al. 2019")
    cell_fn = "%s.cell.nml"%cell_id

    # Define a cell
    cell = cell_doc.add(
        "Cell", id=cell_id, notes="%s cell from Nicoletti et al. 2019"%cell_id
    )  
    diam = 10
    cell.add_segment(
        prox=[0, 0, 0, diam],
        dist=[0, 0, 0, diam],
        name="soma",
        parent=None,
        fraction_along=1.0,
        seg_type='soma'
    )

    
    # Leak channel
    cell.add_channel_density(
        cell_doc,
        cd_id="leak_chans",
        cond_density="3.0 S_per_m2",
        erev="%smV"%eleak,
        ion="non_specific",
        ion_channel="leak",
        ion_chan_def_file=create_channel_file('leak'),
    )

    cell.set_specific_capacitance("%s uF_per_cm2"%c)

    cell.add_membrane_property("SpikeThresh", value="0mV")
    cell.set_init_memb_potential("-75mV")

    # This value is not really used as it's a single comp cell model
    cell.set_resistivity("0.1 kohm_cm")

    cell.info(show_contents=True)

    cell_doc.validate(recursive=True)
    pynml.write_neuroml2_file(
        nml2_doc=cell_doc, nml2_file_name=cell_fn, validate=True
    )

    sim, net = generate_nmllite(cell_id, duration=700, config='IClamp', parameters = None)

    ################################################################################
    ###   Run in some simulators

    from neuromllite.NetworkGenerator import check_to_generate_or_run
    import sys

    check_to_generate_or_run(sys.argv, sim)