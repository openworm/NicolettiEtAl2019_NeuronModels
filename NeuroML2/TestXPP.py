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
channels_to_include = ["leak", "nca", "shal", "egl36", "kir", "shak", "cca", "unc2", "egl19"]

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


run_xpp_file(new_ode_file, plot=True, show_plot_already=False)

from matplotlib import pyplot as plt

from pyneuroml.pynml import reload_standard_dat_file

data, indices = reload_standard_dat_file("Sim_IClamp_RMD.pop_RMD.v.dat")
t_ms = [t * 1000 for t in data["t"]]
v_mV = [v * 1000 for v in data[0]]
plt.plot(t_ms, v_mV, label="NeuroML - v", linewidth=0.5)
plt.legend()

plt.show()
