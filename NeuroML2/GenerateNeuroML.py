from neuroml import NeuroMLDocument
from neuroml.utils import component_factory

from pyneuroml import pynml
from pyneuroml.xppaut import parse_script

from neuroml import GateHHRates


xpps = {"RMD": parse_script("../RMD.ode"), "AWCon": parse_script("../AWC.ode")}



colors = {"AWCon": "0 0 0.8", "RMD": "0 0.8 0", "GenericMuscleCell": "0.8 0 0"}


def create_channel_file(chan_id_in_cell, cell_id, xpp, species, gates={},extra_params=[]):
    chan_id = "%s_%s" % (cell_id, chan_id_in_cell)
    chan_doc = NeuroMLDocument(
        id=chan_id,
        notes="An ion channel from cell %s from Nicoletti et al. 2019" % cell_id,
    )

    chan_fn = "%s.channel.nml" % (chan_id)

    channel = component_factory(
        "IonChannelHH",
        id=chan_id,
        conductance="10pS",
        species=species,
        notes="%s channel from Nicoletti et al. 2019" % chan_id,
    )

    for g in gates:

        ss = component_factory("HHVariable", type="%s_%s_inf"%(chan_id_in_cell,g)) 
        tc = component_factory("HHTime", type="%s_%s_tau"%(chan_id_in_cell,g)) 

        gc = component_factory("GateHHTauInf", id=g, instances=gates[g][0], steady_state=ss, time_course=tc, validate=True) 
        channel.add(gc)
        
        ssct = component_factory("ComponentType", name=ss.type, extends="baseVoltageDepVariable")
        tcct = component_factory("ComponentType", name=tc.type, extends="baseVoltageDepTime")
        chan_doc.add(ssct)
        chan_doc.add(tcct)

        vscale = component_factory("Constant", name='VOLT_SCALE',  dimension="voltage", value="1 mV")
        ssct.add(vscale)
        tcct.add(vscale)
        tscale = component_factory("Constant", name='TIME_SCALE',  dimension="time", value="1 ms")
        tcct.add(tscale)


        inf_expr = xpp["derived_variables"][gates[g][1]]
        tau_expr = xpp["derived_variables"][gates[g][2]]

        potential_parameters = []
        for p in xpp["parameters"]:
            if chan_id_in_cell in p:
                potential_parameters.append(p)

        for p in potential_parameters+extra_params:
            const = component_factory(
                "Constant", name=p, dimension="none", value=str(xpp["parameters"][p])
            )
            if p in inf_expr:
                ssct.add(const)
            if p in tau_expr:
                tcct.add(const)
        

        d = component_factory("Dynamics")
        ssct.add(d)
        dvv = component_factory("DerivedVariable", name="V", dimension="none", value="(v) / VOLT_SCALE")
        d.add(dvv)
        

        import sympy
        from sympy.parsing.sympy_parser import parse_expr
        v,V = sympy.symbols('v V')

        s_expr = parse_expr(inf_expr, evaluate=False)
        s_expr = s_expr.subs(v,V)

        dv = component_factory("DerivedVariable", name="x", exposure="x", dimension="none", value=s_expr)
        d.add(dv)
        
        d = component_factory("Dynamics")
        tcct.add(d)
        d.add(dvv)

        s_expr = parse_expr(tau_expr, evaluate=False)
        s_expr = s_expr.subs(v,V)

        dv = component_factory("DerivedVariable", name="t", exposure="t", dimension="time", value='(%s)/1000'%s_expr)
        d.add(dv)
        
            


    chan_doc.add(channel)
    #chan_doc.validate(recursive=True)

    pynml.write_neuroml2_file(nml2_doc=chan_doc, nml2_file_name=chan_fn, validate=True)

    return chan_fn


def generate_nmllite(cell, duration=700, config="IClamp", parameters=None):
    from neuromllite import Cell, InputSource

    # from neuromllite.NetworkGenerator import *
    from neuromllite.utils import create_new_model

    reference = "%s_%s" % (config, cell)

    cell_id = "%s" % cell
    cell_nmll = Cell(id=cell_id, neuroml2_source_file="%s.cell.nml" % (cell))

    ################################################################################
    ###   Add some inputs

    if "IClamp" in config:
        if not parameters:
            parameters = {}
            parameters["stim_amp"] = "10pA"

        input_source = InputSource(
            id="iclamp_0",
            neuroml2_input="PulseGenerator",
            parameters={"amplitude": "stim_amp", "delay": "310ms", "duration": "500ms"},
        )

    else:
        if not parameters:
            parameters = {}
            parameters["average_rate"] = "100 Hz"
            parameters["number_per_cell"] = "10"

        input_source = InputSource(
            id="pfs0",
            neuroml2_input="PoissonFiringSynapse",
            parameters={
                "average_rate": "average_rate",
                "synapse": syn_exc.id,
                "spike_target": "./%s" % syn_exc.id,
            },
        )

    sim, net = create_new_model(
        reference,
        duration,
        dt=0.025,  # ms
        temperature=34,  # degC
        default_region="Worm",
        parameters=parameters,
        cell_for_default_population=cell_nmll,
        color_for_default_population=colors[cell],
        input_for_default_population=input_source,
    )

    return sim, net


def create_cells(channels_to_include):
    for cell_id in ["AWCon", "RMD"]:
        # Create the nml file and add the ion channels
        cell_doc = NeuroMLDocument(
            id=cell_id, notes="A cell from Nicoletti et al. 2019"
        )
        cell_fn = "%s.cell.nml" % cell_id

        # Define a cell
        cell = cell_doc.add(
            "Cell", id=cell_id, notes="%s cell from Nicoletti et al. 2019" % cell_id
        )
        diam = 17.841242  # Gives a convenient surface area of 1000.0 um^2
        cell.add_segment(
            prox=[0, 0, 0, diam],
            dist=[0, 0, 0, diam],
            name="soma",
            parent=None,
            fraction_along=1.0,
            seg_type="soma",
        )

        # Leak channel
        if 'leak' in channels_to_include:
            cell.add_channel_density(
                cell_doc,
                cd_id="leak_chans",
                cond_density="%s S_per_m2" % xpps[cell_id]["parameters"]["gleak"],
                erev="%smV" % xpps[cell_id]["parameters"]["eleak"],
                ion="non_specific",
                ion_channel="%s_leak" % cell_id,
                ion_chan_def_file=create_channel_file("leak", cell_id, xpps[cell_id], species='non_specific'),
            )

        # SHK1 CHANNELS
        if 'shak' in channels_to_include:
            chan_id = 'shak'
            ion = 'k'
            g_param = 'gshak'
            gates={'m':[1,'minf_shak','tm_shak'],'h':[1,'hinf_shak','th_shak']}
            extra_params = ['shiftV05']

            cell.add_channel_density(
                cell_doc,
                cd_id="%s_chans"%chan_id,
                cond_density="%s S_per_m2" % xpps[cell_id]["parameters"][g_param],
                erev="%smV" % xpps[cell_id]["parameters"]["e%s"%ion],
                ion=ion,
                ion_channel="%s_%s" % (cell_id,chan_id),
                ion_chan_def_file=create_channel_file(
                    chan_id,
                    cell_id,
                    xpps[cell_id],
                    species=ion,
                    gates=gates,
                    extra_params = extra_params
                ),
            )

        # IRK/Kir channel
        if 'kir' in channels_to_include:
            cell.add_channel_density(
                cell_doc,
                cd_id="kir_chans",
                cond_density="%s S_per_m2" % xpps[cell_id]["parameters"]["gkir"],
                erev="%smV" % xpps[cell_id]["parameters"]["ek"],
                ion="k",
                ion_channel="%s_kir" % cell_id,
                ion_chan_def_file=create_channel_file(
                    "kir",
                    cell_id,
                    xpps[cell_id],
                    species='k',
                    gates={'m':[1,'minf_kir','tm_kir']},
                ),
            )

        # CCA-1 channels
        if 'cca' in channels_to_include:
            chan_id = 'cca'
            ion = 'ca'

            cell.add_channel_density(
                cell_doc,
                cd_id="%s_chans"%chan_id,
                cond_density="%s S_per_m2" % xpps[cell_id]["parameters"]["gcca1"],
                erev="%smV" % xpps[cell_id]["parameters"]["e%s"%ion],
                ion=ion,
                ion_channel="%s_%s" % (cell_id,chan_id),
                ion_chan_def_file=create_channel_file(
                    chan_id,
                    cell_id,
                    xpps[cell_id],
                    species=ion,
                    gates={'m':[2,'minf_cca1','tm_cca1'],'h':[1,'hinf_cca1','th_cca1']},
                    extra_params = ['f3ca','f4ca']
                ),
            )

        cell.set_specific_capacitance("%s F_per_m2" % (float(xpps[cell_id]["parameters"]["c"]) * 1e-3))

        cell.add_membrane_property("SpikeThresh", value="0mV")
        cell.set_init_memb_potential("-70mV")

        # This value is not really used as it's a single comp cell model
        cell.set_resistivity("0.1 kohm_cm")

        cell.info(show_contents=True)

        cell_doc.validate(recursive=True)
        pynml.write_neuroml2_file(
            nml2_doc=cell_doc, nml2_file_name=cell_fn, validate=True
        )

        sim, net = generate_nmllite(
            cell_id, duration=1400, config="IClamp", parameters=None
        )

        ################################################################################
        ###   Run in some simulators

        from neuromllite.NetworkGenerator import check_to_generate_or_run
        import sys

        check_to_generate_or_run(sys.argv, sim)


if __name__ == "__main__":

    channels_to_include = ['leak']
    channels_to_include = ['leak','kir']
    channels_to_include = ['leak','kir','shak','cca']
    create_cells(channels_to_include)
