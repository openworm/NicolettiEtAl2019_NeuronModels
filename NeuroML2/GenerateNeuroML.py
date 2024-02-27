from neuroml import NeuroMLDocument
from neuroml.utils import component_factory

from pyneuroml import pynml

c=1.2 

ek=-80
eca=60
eleak=-80
ena=30

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

    # This value is not really used as it's a single comp cell model
    cell.set_resistivity("0.1 kohm_cm")

    cell.info(show_contents=True)

    cell_doc.validate(recursive=True)
    pynml.write_neuroml2_file(
        nml2_doc=cell_doc, nml2_file_name=cell_fn, validate=True
    )