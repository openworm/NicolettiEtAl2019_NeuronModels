from neuroml import NeuroMLDocument
from neuroml.utils import component_factory

from pyneuroml import pynml
from pyneuroml.xppaut import parse_script
from pprint import pprint

from neuroml import GateHHRates
from neuroml import IncludeType
import sympy
from sympy.parsing.sympy_parser import parse_expr


xpps = {"RMD": parse_script("../RMD.ode"), "AWCon": parse_script("../AWC.ode")}
# xpps = {"RMD": parse_script("../RMD.ode")}
# pprint(xpps)

colors = {"AWCon": "0 0 0.8", "RMD": "0 0.8 0", "GenericMuscleCell": "0.8 0 0"}
local_sympy_dict = {"beta": sympy.Symbol("beta")}


def replace_v(expr):
    v, V = sympy.symbols("v V")
    expr = expr.replace("^", "**")
    expr = parse_expr(expr, evaluate=False, local_dict=local_sympy_dict)
    expr = str(expr.subs(v, V))
    expr = expr.replace("**", "^")
    expr = expr.replace("Abs", "abs")

    return expr


def _add_tc_ss(
    chan_id_in_cell,
    gate_id,
    inf_expr_param,
    tau_expr_param,
    chan_doc,
    xpp,
    extra_params,
    ca_conc_var=None,
    add_all=False,
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

    tscale_ms = component_factory(
        "Constant", name="TIME_SCALE", dimension="time", value="1 ms"
    )
    tscale_s = component_factory(
        "Constant", name="TIME_SCALE_S", dimension="time", value="1 s"
    )
    tcct.add(tscale_ms)
    if add_all:
        ssct.add(tscale_ms)
        # tcct.add(tscale_s)
        # ssct.add(tscale_s)

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

    dss = component_factory("Dynamics")
    ssct.add(dss)

    dtc = component_factory("Dynamics")
    tcct.add(dtc)

    dvv = component_factory(
        "DerivedVariable", name="V", dimension="none", value="(v) / VOLT_SCALE"
    )
    dss.add(dvv)

    inf_tau_exprs = {inf_expr: ssct, tau_expr: tcct}
    inf_tau_all_deps = {inf_expr: inf_expr, tau_expr: tau_expr}

    for p in potential_parameters + extra_params:
        if p in xpp["parameters"]:
            val = xpp["parameters"][p]
            const = component_factory(
                "Constant", name=p, dimension="none", value=str(val)
            )
            for e in inf_tau_exprs:
                if p in inf_tau_all_deps[e] or add_all:
                    inf_tau_exprs[e].add(const)
                    inf_tau_all_deps[e] += "__%s" % p

        elif p in xpp["time_derivatives"]:
            val0 = replace_v(str(xpp["time_derivatives"][p]))

            td = component_factory(
                "TimeDerivative", variable=p, value="(%s) / %s" % (val0, tscale_ms.name)
            )
            sv = component_factory("StateVariable", name=p, dimension="none")
            for e in inf_tau_exprs:
                if p in inf_tau_all_deps[e]:
                    inf_tau_exprs[e].Dynamics[0].add(sv)
                    inf_tau_exprs[e].Dynamics[0].add(td)
                    inf_tau_all_deps[e] += "__%s" % p

        elif p in xpp["derived_variables"]:
            val = replace_v(str(xpp["derived_variables"][p]))
            dv = component_factory(
                "DerivedVariable", name=p, dimension="none", value=val
            )
            for e in inf_tau_exprs:
                if p in inf_tau_all_deps[e]:
                    inf_tau_exprs[e].Dynamics[0].add(dv)
                    inf_tau_all_deps[e] += "__%s" % p

            # print("====  %s" % dv)
            for pp in potential_parameters + extra_params:
                # print(pp)
                if pp in dv.value:
                    # print("---  %s" % pp)
                    if pp in xpp["parameters"]:
                        # print(">>>" + pp)
                        val2 = str(xpp["parameters"][pp])
                        const2 = component_factory(
                            "Constant", name=pp, dimension="none", value=val2
                        )
                        for e in inf_tau_exprs:
                            if pp in inf_tau_all_deps[e] or add_all:
                                inf_tau_exprs[e].add(const2)

                    if pp in xpp["derived_variables"]:
                        # print(".." + pp)
                        val3 = replace_v(str(xpp["derived_variables"][pp]))
                        dv3 = component_factory(
                            "DerivedVariable", name=pp, dimension="none", value=val3
                        )
                        for e in inf_tau_exprs:
                            if pp in inf_tau_all_deps[e] or add_all:
                                inf_tau_exprs[e].Dynamics[0].add(dv3)

    if ca_conc_var:
        dvca = component_factory(
            "DerivedVariable",
            name=ca_conc_var,
            dimension="none",
            value="(caConc) / CONC_SCALE",
        )
        dss.add(dvca)

    s_expr = replace_v(inf_expr)

    dv = component_factory(
        "DerivedVariable", name="x", exposure="x", dimension="none", value=s_expr
    )
    dss.add(dv)

    dtc.add(dvv)

    if verbose:
        print("Converting %s" % tau_expr)

    s_expr = replace_v(tau_expr)

    dv = component_factory(
        "DerivedVariable",
        name="t",
        exposure="t",
        dimension="time",
        value="(%s)* %s" % (s_expr, tscale_ms.name),
    )
    dtc.add(dv)

    return ss, tc


def create_channel_file(
    chan_id_in_cell,
    cell_id,
    xpp,
    species,
    gates={},
    extra_params=[],
    ca_conc_var=None,
    add_all=False,
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
                add_all,
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

        if c != "egl36" and cell != "AWCon":
            sim.record_variables["biophys/membraneProperties/%s_chans/gDensity" % c] = {
                "all": "*"
            }
            sim.record_variables["biophys/membraneProperties/%s_chans/iDensity" % c] = {
                "all": "*"
            }
            if c != "leak" and c != "nca":
                sim.record_variables[
                    "biophys/membraneProperties/%s_chans/%s_%s/m/q" % (c, cell, c)
                ] = {"all": "*"}
            if c != "leak" and c not in ["nca", "kir", "sk"]:
                sim.record_variables[
                    "biophys/membraneProperties/%s_chans/%s_%s/h/q" % (c, cell, c)
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

        # SLO1-UNC2 COMPLEX
        if "bk" in channels_to_include:
            chan_id = "bk"
            ion = "k"
            g_param = "gbk"
            gates = {"m": [1, "minf_bk", "tm_bkunc2"], "h": [1, "hinf_unc2", "th_unc2"]}
            extra_params = [
                "backgr",
                "cac_nano",
                "wom",
                "wyx",
                "kyx",
                "nyx",
                "wop",
                "wxy",
                "kxy",
                "nxy",
                "sth2",
                "tm_unc2",
                "pi",
                "r",
                "d",
                "F",
                "kb",
                "b",
                "gsc",
                "eca",
                "cao_nano",
                "kcm",
                "kom",
                "kop",
                "minf_unc2",
                "alpha",
                "beta",
                "stm2",
                "m_unc2",
            ]
            for p in xpps[cell_id]["parameters"]:
                if "unc2" in p:
                    extra_params.append(p)

            xpps[cell_id]["parameters"]["pi"] = 3.14159265359

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
                    add_all=True,
                ),
            )

        # SLO1-EGL19 COMPLEX
        if "slo1" in channels_to_include:
            chan_id = "slo1"
            ion = "k"
            g_param = "gslo1"
            gates = {
                "m": [1, "minf_slo1", "tm_slo1"],
                "h": [1, "hinf_egl19", "ths_egl19"],
            }
            extra_params = [
                "backgr",
                "cac_nano",
                "wom",
                "wyx",
                "kyx",
                "nyx",
                "wop",
                "wxy",
                "kxy",
                "nxy",
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
                "tm_egl19",
                "pi",
                "r",
                "d",
                "F",
                "kb",
                "b",
                "gsc",
                "eca",
                "cao_nano",
                "kcm2",
                "kom2",
                "kop2",
                "kop",
                "minf_egl19",
                "alpha1",
                "beta1",
                "stm2",
                "m_egl19",
            ]
            for p in xpps[cell_id]["parameters"]:
                if "egl19" in p:
                    extra_params.append(p)

            xpps[cell_id]["parameters"]["pi"] = 3.14159265359

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
                    add_all=True,
                ),
            )

        # SLO2-EGL19 COMPLEX
        if "bk2" in channels_to_include:
            chan_id = "bk2"
            ion = "k"
            g_param = "gbk2"
            gates = {
                "m": [1, "minf_bk2", "tm_bkegl19"],
                "h": [1, "hinf_egl19", "ths_egl19"],
            }
            extra_params = [
                "backgr",
                "cac_nano",
                "wom1",
                "wyx1",
                "kyx1",
                "nyx1",
                "wop1",
                "wxy1",
                "kxy1",
                "nxy1",
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
                "tm_egl19",
                "pi",
                "r",
                "d",
                "F",
                "kb",
                "b",
                "gsc",
                "eca",
                "cao_nano",
                "kcm1",
                "kom1",
                "kop1",
                "minf_egl19",
                "alpha1",
                "beta1",
                "stm2",
                "m_egl19",
            ]
            for p in xpps[cell_id]["parameters"]:
                if "egl19" in p:
                    extra_params.append(p)

            xpps[cell_id]["parameters"]["pi"] = 3.14159265359

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
                    add_all=True,
                ),
            )

        # SLO2-UNC2 COMPLEX
        if "slo2" in channels_to_include:
            chan_id = "slo2"
            ion = "k"
            g_param = "gslo2"
            gates = {"m": [1, "minf_slo2", "tm_slo2"], "h": [1, "hinf_unc2", "th_unc2"]}
            extra_params = [
                "backgr",
                "cac_nano",
                "wom1",
                "wyx1",
                "kyx1",
                "nyx1",
                "wop1",
                "wxy1",
                "kxy1",
                "nxy1",
                "sth2",
                "tm_unc2",
                "pi",
                "r",
                "d",
                "F",
                "kb",
                "b",
                "gsc",
                "eca",
                "cao_nano",
                "kcm3",
                "kom3",
                "kop3",
                "minf_unc2",
                "alpha",
                "beta",
                "stm2",
                "m_unc2",
            ]
            for p in xpps[cell_id]["parameters"]:
                if "unc2" in p:
                    extra_params.append(p)

            xpps[cell_id]["parameters"]["pi"] = 3.14159265359

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
                    add_all=True,
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
    channels_to_include = ["leak", "bk"]
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
    channels_to_include = ["leak", "unc2"]
    channels_to_include = ["leak", "unc2", "slo2"]
    channels_to_include = ["leak", "unc2", "bk"]
    channels_to_include = ["leak", "egl19", "slo1"]
    channels_to_include = ["leak", "egl19", "bk2"]
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
        "bk",
        "slo1",
        "bk2",
        "slo2",
    ]

    create_cells(channels_to_include, duration=1400, stim_delay=1310, stim_duration=50)
