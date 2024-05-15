from neuroml import NeuroMLDocument
from neuroml.utils import component_factory

from pyneuroml import pynml
from pyneuroml.xppaut import parse_script
from pprint import pprint

from neuroml import GateHHRates
from neuroml import IncludeType


xpps = {"RMD": parse_script("../RMD.ode"), "AWCon": parse_script("../AWC.ode")}
# xpps = {"RMD": parse_script("../RMD.ode")}
# pprint(xpps)

colors = {"AWCon": "0 0 0.8", "RMD": "0 0.8 0", "GenericMuscleCell": "0.8 0 0"}


def _add_tc_ss(
    chan_id_in_cell,
    gate_id,
    inf_expr_param,
    tau_expr_param,
    chan_doc,
    xpp,
    extra_params,
    ca_conc_var=None,
):
    verbose = True
    ss = component_factory("HHVariable", type="%s_%s_inf" % (chan_id_in_cell, gate_id))
    tc = component_factory("HHTime", type="%s_%s_tau" % (chan_id_in_cell, gate_id))

    ssct = component_factory(
        "ComponentType",
        name=ss.type,
        extends="baseVoltageDepVariable"
        if not ca_conc_var
        else "baseVoltageConcDepVariable",
    )
    tcct = component_factory(
        "ComponentType", name=tc.type, extends="baseVoltageDepTime"
    )
    chan_doc.add(ssct)
    chan_doc.add(tcct)

    vscale = component_factory(
        "Constant", name="VOLT_SCALE", dimension="voltage", value="1 mV"
    )
    ssct.add(vscale)
    tcct.add(vscale)

    if ca_conc_var:
        cascale = component_factory(
            "Constant",
            name="CONC_SCALE",
            dimension="concentration",
            value="1 mM",
        )
        ssct.add(cascale)

    tscale = component_factory(
        "Constant", name="TIME_SCALE", dimension="time", value="1 ms"
    )
    tcct.add(tscale)

    inf_expr = xpp["derived_variables"][inf_expr_param]
    tau_expr = xpp["derived_variables"][tau_expr_param]

    if verbose:
        print(
            "For channel %s, gate %s, inf = [%s], tau = [%s]"
            % (chan_id_in_cell, gate_id, inf_expr, tau_expr)
        )

    potential_parameters = []
    for p in xpp["parameters"]:
        if chan_id_in_cell in p:
            potential_parameters.append(p)

    for p in potential_parameters + extra_params:
        const = component_factory(
            "Constant", name=p, dimension="none", value=str(xpp["parameters"][p])
        )
        if p in inf_expr:
            ssct.add(const)
        if p in tau_expr:
            tcct.add(const)

    d = component_factory("Dynamics")
    ssct.add(d)

    dvv = component_factory(
        "DerivedVariable", name="V", dimension="none", value="(v) / VOLT_SCALE"
    )
    d.add(dvv)

    if ca_conc_var:
        dvca = component_factory(
            "DerivedVariable",
            name=ca_conc_var,
            dimension="none",
            value="(caConc) / CONC_SCALE",
        )
        d.add(dvca)

    import sympy
    from sympy.parsing.sympy_parser import parse_expr

    v, V = sympy.symbols("v V")

    s_expr = parse_expr(inf_expr, evaluate=False)
    s_expr = s_expr.subs(v, V)

    dv = component_factory(
        "DerivedVariable", name="x", exposure="x", dimension="none", value=s_expr
    )
    d.add(dv)

    d = component_factory("Dynamics")
    tcct.add(d)
    d.add(dvv)

    if verbose:
        print("Converting %s" % tau_expr)
    s_expr = parse_expr(tau_expr.replace("^", "**"), evaluate=False)
    s_expr = str(s_expr.subs(v, V)).replace("**", "^")

    dv = component_factory(
        "DerivedVariable",
        name="t",
        exposure="t",
        dimension="time",
        value="(%s)* TIME_SCALE" % s_expr,
    )
    d.add(dv)

    return ss, tc


def create_channel_file(
    chan_id_in_cell, cell_id, xpp, species, gates={}, extra_params=[], ca_conc_var=None
):
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

    for gate_id in gates:
        if type(gates[gate_id]) == list:
            ss, tc = _add_tc_ss(
                chan_id_in_cell,
                gate_id,
                gates[gate_id][1],
                gates[gate_id][2],
                chan_doc,
                xpp,
                extra_params,
                ca_conc_var,
            )

            gc = component_factory(
                "GateHHUndetermined",
                type="gateHHtauInf",
                id=gate_id,
                instances=gates[gate_id][0],
                steady_state=ss,
                time_course=tc,
                validate=True,
            )
            channel.add(gc)

        elif type(gates[gate_id]) == dict:
            gc = component_factory(
                "GateHHUndetermined",
                type="gateFractional",
                id=gate_id,
                instances=1,
                validate=False,
            )

            for sub_gate in gates[gate_id]:
                # gates={'m':[3,'minf_shal','tm_shal'],'h':{'hf': [0.7,'hinf_shal','thf_shal'],'hs': [0.3,'hinf_shal','ths_shal']} }
                fract = gates[gate_id][sub_gate][0]

                ss, tc = _add_tc_ss(
                    chan_id_in_cell,
                    "%s_%s" % (gate_id, sub_gate),
                    gates[gate_id][sub_gate][1],
                    gates[gate_id][sub_gate][2],
                    chan_doc,
                    xpp,
                    extra_params,
                )

                sg = component_factory(
                    "GateFractionalSubgate",
                    id=sub_gate,
                    fractional_conductance=str(fract),
                    steady_state=ss,
                    time_course=tc,
                    validate=True,
                )
                gc.add(sg)
            channel.add(gc)

    chan_doc.add(channel)
    chan_doc.validate(recursive=True)

    pynml.write_neuroml2_file(nml2_doc=chan_doc, nml2_file_name=chan_fn, validate=True)

    return chan_fn


def generate_nmllite(
    cell,
    duration=700,
    config="IClamp",
    parameters=None,
    stim_delay=310,
    stim_duration=500,
    channels_to_include=[],
):
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
            parameters["stim_delay"] = "%sms" % stim_delay
            parameters["stim_duration"] = "%sms" % stim_duration

        input_source = InputSource(
            id="iclamp_0",
            neuroml2_input="PulseGenerator",
            parameters={
                "amplitude": "stim_amp",
                "delay": "stim_delay",
                "duration": "stim_duration",
            },
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
    sim.record_variables = {"caConc": {"all": "*"}}
    for c in channels_to_include:
        if c == "ca":
            c = "sk"

        if c is not "egl36" and cell is not "AWCon":
            sim.record_variables["biophys/membraneProperties/%s_chans/gDensity" % c] = {
                "all": "*"
            }
            sim.record_variables["biophys/membraneProperties/%s_chans/iDensity" % c] = {
                "all": "*"
            }
            if c is not "leak" and c is not "nca":
                sim.record_variables[
                    "biophys/membraneProperties/%s_chans/%s_%s/m/q" % (c, cell, c)
                ] = {"all": "*"}
    sim.to_json_file()

    return sim, net


def create_cells(channels_to_include, duration=700, stim_delay=310, stim_duration=500):
    for cell_id in xpps.keys():
        # Create the nml file and add the ion channels
        cell_doc = NeuroMLDocument(
            id=cell_id, notes="A cell from Nicoletti et al. 2019"
        )
        cell_fn = "%s.cell.nml" % cell_id

        # Define a cell
        cell = cell_doc.add(
            "Cell", id=cell_id, notes="%s cell from Nicoletti et al. 2019" % cell_id
        )
        diam = 1.7841242  # Gives a convenient surface area of ?? um^2
        length = 2.26  # Gives a convenient surface area of ?? um^2
        cell.add_segment(
            prox=[0, 0, 0, diam],
            dist=[0, length, 0, diam],
            name="soma",
            parent=None,
            fraction_along=1.0,
            seg_type="soma",
        )

        density_factor = 1000 / 12.66728074437459

        # Leak channel
        if "leak" in channels_to_include:
            cell.add_channel_density(
                cell_doc,
                cd_id="leak_chans",
                cond_density="%s S_per_m2"
                % (float(xpps[cell_id]["parameters"]["gleak"]) * density_factor),
                erev="%smV" % xpps[cell_id]["parameters"]["eleak"],
                ion="non_specific",
                ion_channel="%s_leak" % cell_id,
                ion_chan_def_file=create_channel_file(
                    "leak", cell_id, xpps[cell_id], species="non_specific"
                ),
            )

        # Nca channel
        if "nca" in channels_to_include:
            cell.add_channel_density(
                cell_doc,
                cd_id="nca_chans",
                cond_density="%s S_per_m2"
                % (float(xpps[cell_id]["parameters"]["gnca"]) * density_factor),
                erev="%smV" % xpps[cell_id]["parameters"]["ena"],
                ion="non_specific",
                ion_channel="%s_nca" % cell_id,
                ion_chan_def_file=create_channel_file(
                    "nca", cell_id, xpps[cell_id], species="non_specific"
                ),
            )

        # SHL-1 CHANNELS
        if "shal" in channels_to_include:
            chan_id = "shal"
            ion = "k"
            g_param = "gshal"
            gates = {
                "m": [3, "minf_shal", "tm_shal"],
                "h": {
                    "hf": [0.7, "hinf_shal", "thf_shal"],
                    "hs": [0.3, "hinf_shal", "ths_shal"],
                },
            }
            extra_params = []

            cell.add_channel_density(
                cell_doc,
                cd_id="%s_chans" % chan_id,
                cond_density="%s S_per_m2"
                % (float(xpps[cell_id]["parameters"][g_param]) * density_factor),
                erev="%smV" % xpps[cell_id]["parameters"]["e%s" % ion],
                ion=ion,
                ion_channel="%s_%s" % (cell_id, chan_id),
                ion_chan_def_file=create_channel_file(
                    chan_id,
                    cell_id,
                    xpps[cell_id],
                    species=ion,
                    gates=gates,
                    extra_params=extra_params,
                ),
            )

        # SHK1 CHANNELS
        if "shak" in channels_to_include:
            chan_id = "shak"
            ion = "k"
            g_param = "gshak"
            gates = {"m": [1, "minf_shak", "tm_shak"], "h": [1, "hinf_shak", "th_shak"]}
            extra_params = ["shiftV05"]

            cell.add_channel_density(
                cell_doc,
                cd_id="%s_chans" % chan_id,
                cond_density="%s S_per_m2"
                % (float(xpps[cell_id]["parameters"][g_param]) * density_factor),
                erev="%smV" % xpps[cell_id]["parameters"]["e%s" % ion],
                ion=ion,
                ion_channel="%s_%s" % (cell_id, chan_id),
                ion_chan_def_file=create_channel_file(
                    chan_id,
                    cell_id,
                    xpps[cell_id],
                    species=ion,
                    gates=gates,
                    extra_params=extra_params,
                ),
            )

        # EGL-36 CHANNELS
        if cell_id != "AWCon" and "egl36" in channels_to_include:
            chan_id = "egl36"
            ion = "k"
            g_param = "gegl36"

            gates = {
                "m": {
                    "m1": [
                        xpps[cell_id]["parameters"]["a1"],
                        "min1_egl36",
                        "tm1_egl36",
                    ],
                    "m2": [
                        xpps[cell_id]["parameters"]["a2"],
                        "min2_egl36",
                        "tm2_egl36",
                    ],
                    "m3": [
                        xpps[cell_id]["parameters"]["a3"],
                        "min3_egl36",
                        "tm3_egl36",
                    ],
                }
            }

            extra_params = []

            # print(xpps[cell_id]["parameters"].keys())

            cell.add_channel_density(
                cell_doc,
                cd_id="%s_chans" % chan_id,
                cond_density="%s S_per_m2"
                % (float(xpps[cell_id]["parameters"][g_param]) * density_factor),
                erev="%smV" % xpps[cell_id]["parameters"]["e%s" % ion],
                ion=ion,
                ion_channel="%s_%s" % (cell_id, chan_id),
                ion_chan_def_file=create_channel_file(
                    chan_id,
                    cell_id,
                    xpps[cell_id],
                    species=ion,
                    gates=gates,
                    extra_params=extra_params,
                ),
            )

        # IRK/Kir channel
        if "kir" in channels_to_include:
            cell.add_channel_density(
                cell_doc,
                cd_id="kir_chans",
                cond_density="%s S_per_m2"
                % (float(xpps[cell_id]["parameters"]["gkir"]) * density_factor),
                erev="%smV" % xpps[cell_id]["parameters"]["ek"],
                ion="k",
                ion_channel="%s_kir" % cell_id,
                ion_chan_def_file=create_channel_file(
                    "kir",
                    cell_id,
                    xpps[cell_id],
                    species="k",
                    gates={"m": [1, "minf_kir", "tm_kir"]},
                ),
            )

        # UNC-2 CHANNELS
        if "unc2" in channels_to_include:
            chan_id = "unc2"
            ion = "ca"
            g_param = "gunc2"
            gates = {"m": [1, "minf_unc2", "tm_unc2"], "h": [1, "hinf_unc2", "th_unc2"]}
            extra_params = ["stm2", "sth2"]

            cell.add_channel_density(
                cell_doc,
                cd_id="%s_chans" % chan_id,
                cond_density="%s S_per_m2"
                % (float(xpps[cell_id]["parameters"][g_param]) * density_factor),
                erev="%smV" % xpps[cell_id]["parameters"]["e%s" % ion],
                ion=ion,
                ion_channel="%s_%s" % (cell_id, chan_id),
                ion_chan_def_file=create_channel_file(
                    chan_id,
                    cell_id,
                    xpps[cell_id],
                    species=ion,
                    gates=gates,
                    extra_params=extra_params,
                ),
            )

        # EGL-19 CHANNELS
        if "egl19" in channels_to_include:
            chan_id = "egl19"
            ion = "ca"
            g_param = "gegl19"
            gates = {
                "m": [1, "minf_egl19", "tm_egl19"],
                "h": [1, "hinf_egl19", "ths_egl19"],
            }
            extra_params = [
                "stm19",
                "sth19",
                "pdg1",
                "pdg2",
                "pdg3",
                "pdg4",
                "pdg5",
                "pdg6",
                "pdg7",
                "stau19",
                "pds1",
                "pds2",
                "pds3",
                "pds4",
                "pds5",
                "pds6",
                "pds7",
                "pds8",
                "pds9",
                "pds10",
                "pds11",
                "shiftdps",
            ]

            cell.add_channel_density(
                cell_doc,
                cd_id="%s_chans" % chan_id,
                cond_density="%s S_per_m2"
                % (float(xpps[cell_id]["parameters"][g_param]) * density_factor),
                erev="%smV" % xpps[cell_id]["parameters"]["e%s" % ion],
                ion=ion,
                ion_channel="%s_%s" % (cell_id, chan_id),
                ion_chan_def_file=create_channel_file(
                    chan_id,
                    cell_id,
                    xpps[cell_id],
                    species=ion,
                    gates=gates,
                    extra_params=extra_params,
                ),
            )

        # CCA-1 channels
        if "cca" in channels_to_include:
            chan_id = "cca"
            ion = "ca"

            cell.add_channel_density(
                cell_doc,
                cd_id="%s_chans" % chan_id,
                cond_density="%s S_per_m2"
                % (float(xpps[cell_id]["parameters"]["gcca1"]) * density_factor),
                erev="%smV" % xpps[cell_id]["parameters"]["e%s" % ion],
                ion=ion,
                ion_channel="%s_%s" % (cell_id, chan_id),
                ion_chan_def_file=create_channel_file(
                    chan_id,
                    cell_id,
                    xpps[cell_id],
                    species=ion,
                    gates={
                        "m": [2, "minf_cca1", "tm_cca1"],
                        "h": [1, "hinf_cca1", "th_cca1"],
                    },
                    extra_params=["f3ca", "f4ca"],
                ),
            )

        # KCNL CHANNELS
        if "ca" in channels_to_include:
            chan_id = "sk"
            ion = "k"
            g_param = "gca"
            gates = {"m": [1, "minf_sk", "t_sk"]}
            extra_params = ["k_sk2"]

            cell.add_channel_density(
                cell_doc,
                cd_id="%s_chans" % chan_id,
                cond_density="%s S_per_m2"
                % (float(xpps[cell_id]["parameters"][g_param]) * density_factor),
                erev="%smV" % xpps[cell_id]["parameters"]["e%s" % ion],
                ion=ion,
                ion_channel="%s_%s" % (cell_id, chan_id),
                ion_chan_def_file=create_channel_file(
                    chan_id,
                    cell_id,
                    xpps[cell_id],
                    species=ion,
                    gates=gates,
                    extra_params=extra_params,
                    ca_conc_var="ca_intra1",
                ),
            )

        cell.set_specific_capacitance(
            "%s F_per_m2"
            % (float(xpps[cell_id]["parameters"]["c"]) / 12.66728074437459)
        )

        cell.add_membrane_property("SpikeThresh", value="0mV")
        cell.set_init_memb_potential("-70mV")

        # This value is not really used as it's a single comp cell model
        cell.set_resistivity("0.1 kohm_cm")

        cell_doc.includes.append(IncludeType(href="CaDynamics.nml"))
        # <species id="ca" ion="ca" concentrationModel="CaDynamics" initialConcentration="1e-4 mM" initialExtConcentration="2 mM"/>
        species = component_factory(
            "Species",
            id="ca",
            ion="ca",
            concentration_model="CaDynamics",
            initial_concentration="5e-5 mM",
            initial_ext_concentration="2 mM",
        )

        cell.biophysical_properties.intracellular_properties.add(species)

        cell.info(show_contents=True)

        cell_doc.validate(recursive=True)
        pynml.write_neuroml2_file(
            nml2_doc=cell_doc, nml2_file_name=cell_fn, validate=True
        )

        sim, net = generate_nmllite(
            cell_id,
            duration=duration,
            config="IClamp",
            parameters=None,
            stim_delay=stim_delay,
            stim_duration=stim_duration,
            channels_to_include=channels_to_include,
        )

        ################################################################################
        ###   Run in some simulators

        from neuromllite.NetworkGenerator import check_to_generate_or_run
        import sys

        check_to_generate_or_run(sys.argv, sim)


if __name__ == "__main__":
    channels_to_include = ["leak"]
    channels_to_include = ["leak", "nca"]
    channels_to_include = ["leak", "kir"]
    channels_to_include = ["leak", "kir", "shak", "cca"]
    channels_to_include = ["leak", "shal", "kir", "shak", "cca"]
    channels_to_include = ["leak", "shal", "egl36", "kir", "shak", "cca"]
    channels_to_include = ["leak", "shal", "egl36", "kir", "shak", "cca", "unc2"]
    channels_to_include = ["leak", "unc2"]
    channels_to_include = ["leak", "shal", "egl36", "kir", "shak", "cca", "unc2"]
    channels_to_include = ["leak", "kir"]
    channels_to_include = ["leak", "kir", "cca"]
    channels_to_include = ["leak", "kir", "cca", "ca"]
    channels_to_include = [
        "leak",
        "nca",
        "shal",
        "egl36",
        "kir",
        "shak",
        "cca",
        "unc2",
        "egl19",
        "ca",
    ]

    create_cells(channels_to_include, duration=1800, stim_delay=1310, stim_duration=500)
