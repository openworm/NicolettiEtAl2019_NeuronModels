from neuromllite import *
from neuromllite.NetworkGenerator import *
from neuromllite.utils import create_new_model
import sys


def generate(cell, config, parameters=None):
    reference = "%s_%s" % (config, cell)

    cell_id = "%s" % cell
    cell_nmll = Cell(id=cell_id, neuroml2_source_file="%s.cell.nml" % (cell))

    amps = [-4, 0, 4, 8, 12, 16, 20]
    stim_delay = 100
    stim_dur = 500
    post_stim = 100
    extra_init = 0

    if cell == "RMD":
        amps = [-2, 2, 6, 10]
        stim_delay = 10
        stim_dur = 50
        post_stim = 90
        extra_init = 0

    sim, net = create_new_model(
        reference,
        duration=extra_init + stim_delay + stim_dur + post_stim,
        dt=0.025,  # ms
        temperature=34,  # degC
        default_region="Worm",
        parameters=parameters,
        cell_for_default_population=cell_nmll,
    )

    net.populations[0].size = len(amps)

    for i in amps:
        ins = InputSource(
            id="iclamp_stim_%s" % str(i).replace("-", "min"),
            neuroml2_input="PulseGenerator",
            parameters={
                "amplitude": "%spA" % i,
                "delay": "%sms" % (extra_init + stim_delay),
                "duration": "%sms" % (stim_dur),
            },
        )
        net.input_sources.append(ins)
        net.inputs.append(
            Input(
                id="input_%s" % ins.id,
                input_source=ins.id,
                population=net.populations[0].id,
                cell_ids=[amps.index(i)],
            )
        )
    if config == "Fig7B":
        hyp_amp = -15
        hyp_delay = 50
        hyp_dur = 20
        ins = InputSource(
            id="iclamp_hyp",
            neuroml2_input="PulseGenerator",
            parameters={
                "amplitude": "%spA" % hyp_amp,
                "delay": "%sms" % (extra_init + stim_delay + stim_dur + hyp_delay),
                "duration": "%sms" % (hyp_dur),
            },
        )
        net.input_sources.append(ins)
        net.inputs.append(
            Input(
                id="input_%s" % ins.id,
                input_source=ins.id,
                population=net.populations[0].id,
                percentage=100,
            )
        )

    net.to_json_file()

    return sim, net


if __name__ == "__main__":
    sim, net = generate("AWCon", config="Fig4C")

    check_to_generate_or_run(sys.argv, sim)

    sim, net = generate("RMD", config="Fig7B")

    check_to_generate_or_run(sys.argv, sim)
