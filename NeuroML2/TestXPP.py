from pyneuroml.xppaut import parse_script, to_lems, to_xpp, run_xpp_file

from pprint import pprint

cell = "RMD"

parsed_data = parse_script("../%s.ode" % cell)
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
    "gsc",
    "gnca",
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
    "bk",
]

for c in channels_to_include:
    print("Including channel: %s" % c)
    if c == "cca":
        c = "cca1"
    all_g.remove("g%s" % c)

for p in all_g:
    parsed_data["parameters"][p] = 0

new_ode_file = "Test_%s.ode" % cell
parsed_data["settings"]["total"] = 1800
parsed_data["parameters"]["ton"] = 1310
new_ode = to_xpp(parsed_data, new_ode_file)

mp_fig = "Membrane potentials"
ca_fig = "[Ca2+]"

print("Running XPP file for %s ms..." % parsed_data["settings"]["total"])

chans = [
    ["unc2", "m", "m_unc2"],
    ["unc2", "h", "h_unc2"],
    ["bk", "m", "mbk"],
    ["bk", "h", "h_unc2"],
]

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
