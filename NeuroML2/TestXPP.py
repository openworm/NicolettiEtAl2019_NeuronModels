
from pyneuroml.xppaut import parse_script, to_lems, to_xpp, run_xpp_file


parsed_data = parse_script("../RMD.ode")

all_g = ['gshal','gkir','gshak','gegl36','gunc2','gegl19','gcca1','gslo1','gbk','gbk2','gslo2','gca','gsc','gnca']

all_g.remove('gkir')

for p in all_g:
    parsed_data['parameters'][p] = 0

new_ode_file = 'Test_RMD.ode'
new_ode = to_xpp(parsed_data, new_ode_file)


run_xpp_file(new_ode_file, plot=True)
