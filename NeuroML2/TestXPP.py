
from pyneuroml.xppaut import parse_script, to_lems, to_xpp, run_xpp_file

from pprint import pprint 
parsed_data = parse_script("../RMD.ode")
pprint(parsed_data)

all_g = ['gleak','gshal','gkir','gshak','gegl36','gunc2','gegl19','gcca1','gslo1','gbk','gbk2','gslo2','gca','gsc','gnca']
  

channels_to_include = ['leak']
channels_to_include = ['leak','shal','egl36','kir','shak','cca','unc2']

for c in channels_to_include:
    print('Including channel: %s'%c)
    if c =='cca': c='cca1'
    all_g.remove('g%s'%c)

for p in all_g:
    parsed_data['parameters'][p] = 0

new_ode_file = 'Test_RMD.ode'
parsed_data['settings']['total'] = 800
new_ode = to_xpp(parsed_data, new_ode_file)


run_xpp_file(new_ode_file, plot=True)
