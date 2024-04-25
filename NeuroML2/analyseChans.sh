
pynml-channelanalysis RMD_shak.channel.nml  -iv -clampDe 100 -clampDu 700 -dur 800 -erev -80 -maxV 40 -minV -120  -clampBaseVoltage -80 -stepTargetVoltage 10
pynml-channelanalysis RMD_shal.channel.nml  -iv -clampDe 100 -clampDu 700 -dur 800 -erev -80 -maxV 40 -minV -120  -clampBaseVoltage -80 -stepTargetVoltage 10
pynml-channelanalysis RMD_egl36.channel.nml -iv -clampDe 500 -clampDu 6000 -dur 7000 -erev -80 -maxV 40 -minV -120  -clampBaseVoltage -80 -stepTargetVoltage 10
pynml-channelanalysis RMD_kir.channel.nml   -iv -clampDe 100 -clampDu 700 -dur 800 -erev -80 -maxV 40 -minV -120  -clampBaseVoltage -80 -stepTargetVoltage 10

pynml-channelanalysis RMD_cca.channel.nml   -iv -clampDe 100 -clampDu 200 -dur 400 -erev 60 -maxV 40 -minV -80  -clampBaseVoltage -80 -stepTargetVoltage 10
pynml-channelanalysis RMD_unc2.channel.nml  -iv -clampDe 100 -clampDu 200 -dur 400 -erev 60 -maxV 40 -minV -80  -clampBaseVoltage -80 -stepTargetVoltage 10