from pyneuroml.xppaut import parse_script, to_lems, to_xpp, run_xpp_file

from pprint import pprint

cell = "AWCon"
cell = "RMD"

parsed_data = parse_script("../%s.ode" % cell.replace("on", ""))

import sys

short_simulation = "-short" in sys.argv
if short_simulation:
    parsed_data = parse_script("../XPP_tests/%s_edited.ode" % cell.replace("on", ""))

pprint(parsed_data)

all_g = [
    "gleak",
    "gshal",
    "gkir",
    "gshak",
    "gegl36",
    "gunc2",
    "gegl19",
    "gcca1",
    "gslo1",
    "gbk",
    "gbk2",
    "gslo2",
    "gca",
    "gnca",
    "gkvs1",
    "gkqt3",
    "gegl2",
]


channels_to_include = ["leak", "nca"]
channels_to_include = ["leak", "shal", "egl36", "kir", "shak", "cca", "unc2"]
channels_to_include = ["leak", "unc2"]
channels_to_include = ["leak", "shal", "egl36", "kir", "shak", "cca", "unc2", "egl19"]
channels_to_include = ["leak", "kir"]
channels_to_include = ["leak", "kir", "cca"]
channels_to_include = ["leak", "kir", "cca", "ca"]
channels_to_include = ["leak", "bk"]
channels_to_include = ["leak", "unc2", "bk"]
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
channels_to_include = ["leak", "kir", "cca", "ca"]
channels_to_include = ["leak", "kvs1"]
channels_to_include = ["leak"]
channels_to_include = ["leak", "kqt3"]
channels_to_include = ["leak", "egl2"]
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
    "bk",
    "slo1",
    "bk2",
    "slo2",
    "kvs1",
    "kqt3",
    "egl2",
]

for c in channels_to_include:
    print("Including channel: %s" % c)
    if c == "cca":
        c = "cca1"
    all_g.remove("g%s" % c)

for p in all_g:
    parsed_data["parameters"][p] = 0

new_ode_file = "Test_%s.ode" % cell

# For testing, as some channel gates take quite some time to reach steady state
additional_transient_phase = 2000

parsed_data["settings"]["total"] = 400 + additional_transient_phase
parsed_data["settings"]["trans"] = 0
parsed_data["settings"]["dt"] = 0.005
parsed_data["parameters"]["ton"] = 310 + additional_transient_phase
parsed_data["parameters"]["toff"] = 360 + additional_transient_phase

if short_simulation:
    parsed_data["settings"]["total"] = 70
    parsed_data["parameters"]["ton"] = 10
    parsed_data["parameters"]["toff"] = 60

new_ode = to_xpp(parsed_data, new_ode_file)

mp_fig = "Membrane potentials"
ca_fig = "[Ca2+]"

print("Running XPP file for %s ms..." % parsed_data["settings"]["total"])

chans = [
    ["egl19", "m", "m_egl19"],
    ["egl19", "h", "hs_egl19"],
    ["bk2", "m", "mbk2"],
    ["unc2", "m", "m_unc2"],
    ["unc2", "h", "h_unc2"],
    ["bk", "m", "mbk"],
    ["slo2", "m", "mslo2"],
    ["slo1", "m", "mslo1"],
    ["kir", "m", "m_kir"],
    ["sk", "m", "m_sk"],
    ["shak", "m", "m_shak"],
    ["shak", "h", "h_shak"],
    ["cca", "m", "m_cca1"],
    ["cca", "h", "h_cca1"],
]
if cell == "AWCon":
    chans = [
        ["egl19", "m", "m_egl19"],
        ["egl19", "h", "hs_egl19"],
        ["bk2", "m", "mbk2"],
        ["unc2", "m", "m_unc2"],
        ["unc2", "h", "h_unc2"],
        ["bk", "m", "mbk"],
        ["slo2", "m", "mslo2"],
        ["slo1", "m", "mslo1"],
        ["kir", "m", "m_kir"],
        ["sk", "m", "m_sk"],
        ["shak", "m", "m_shak"],
        ["shak", "h", "h_shak"],
        ["cca", "m", "m_cca1"],
        ["cca", "h", "h_cca1"],
        ["egl2", "m", "m_egl2"],
        ["kqt3", "s", "s_kqt3"],
        ["kqt3", "w", "w_kqt3"],
    ]
"""
chans = [
["kir", "m", "m_kir"],
["cca", "m", "m_cca1"],
["cca", "h", "h_cca1"],]  """


plot_separately = {mp_fig: "v", ca_fig: "ca_intra1"}

for c in chans:
    print("Plotting %s" % c)
    plot_separately["var %s_%s %s" % (c[0], c[1], c[2])] = c[2]


axes = run_xpp_file(
    new_ode_file, plot=True, show_plot_already=False, plot_separately=plot_separately
)

from matplotlib import pyplot as plt

from pyneuroml.pynml import reload_standard_dat_file


data, indices = reload_standard_dat_file("Sim_IClamp_%s.pop_%s.v.dat" % (cell, cell))
t_ms = [t * 1000 for t in data["t"]]
v_mV = [v * 1000 for v in data[0]]

ax = axes[mp_fig]
ax.plot(t_ms, v_mV, label="NeuroML - v", linewidth=0.5)
ax.legend()


data, indices = reload_standard_dat_file("pop_%s_0.caConc.dat" % cell)
ca = [c * 1 for c in data[0]]

ax = axes[ca_fig]
ax.plot(t_ms, ca, label="NeuroML - [Ca2+]", linewidth=0.5)
ax.legend()

for c in chans:
    print("Plotting %s" % c)
    data, indices = reload_standard_dat_file(
        "pop_%s_0.biophys_membraneProperties_%s_chans_%s_%s_%s_q.dat"
        % (cell, c[0], cell, c[0], c[1])
    )

    vv = [v * 1 for v in data[0]]
    fig_name = "var %s_%s %s" % (c[0], c[1], c[2])
    ax = axes[fig_name]
    ax.plot(t_ms, vv, label="NeuroML %s" % fig_name, linewidth=0.5)
    ax.legend()


plt.show()
