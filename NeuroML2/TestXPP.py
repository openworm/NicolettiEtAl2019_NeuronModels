from pyneuroml.xppaut import parse_script, to_lems, to_xpp, run_xpp_file

from pprint import pprint

parsed_data = parse_script("../RMD.ode")
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
]

for c in channels_to_include:
    print("Including channel: %s" % c)
    if c == "cca":
        c = "cca1"
    all_g.remove("g%s" % c)

for p in all_g:
    parsed_data["parameters"][p] = 0

new_ode_file = "Test_RMD.ode"
parsed_data["settings"]["total"] = 1800
parsed_data["parameters"]["ton"] = 1310
new_ode = to_xpp(parsed_data, new_ode_file)

mp_fig = "Membrane potentials"
ca_fig = "[Ca2+]"

print("Running XPP file for %s ms..." % parsed_data["settings"]["total"])
axes = run_xpp_file(
    new_ode_file,
    plot=True,
    show_plot_already=False,
    plot_separately={mp_fig: "v", ca_fig: "ca_intra1"},
)

from matplotlib import pyplot as plt

from pyneuroml.pynml import reload_standard_dat_file

data, indices = reload_standard_dat_file("Sim_IClamp_RMD.pop_RMD.v.dat")
t_ms = [t * 1000 for t in data["t"]]
v_mV = [v * 1000 for v in data[0]]

ax = axes[mp_fig]
ax.plot(t_ms, v_mV, label="NeuroML - v", linewidth=0.5)
ax.legend()

data, indices = reload_standard_dat_file("pop_AWCon_0.caConc.dat")
t_ms = [t * 1000 for t in data["t"]]
ca = [c * 1 for c in data[0]]

ax = axes[ca_fig]
ax.plot(t_ms, ca, label="NeuroML - [Ca2+]", linewidth=0.5)
ax.legend()

plt.show()
