#!/bin/bash

#nohup root -l -q 'NeuRad_sim.C(10000, "sim500MeV10000evrc.root", "par500MeV10000evrc.root", 0.5)' &
#nohup root -l -q 'NeuRad_sim.C(10000, "sim400MeV10000evrc.root", "par400MeV10000evrc.root", 0.4)' &

nohup root -l -q 'NeuRad_digi.C(10000,"sim500MeV10000evrc.root", "par500MeV10000evrc.root", "digi500MeV10000evrc_onepe.root")' &

#nohup root -l -q 'NeuRad_digi_DT.C(10000,"sim500MeV10000evrc.root", "par500MeV10000evrc.root", "digi500MeV10000evrc_DT0001.root", 0.001)' &
#nohup root -l -q 'NeuRad_digi_DT.C(10000,"sim500MeV10000evrc.root", "par500MeV10000evrc.root", "digi500MeV10000evrc_DT001.root", 0.01)' &
#nohup root -l -q 'NeuRad_digi_DT.C(10000,"sim500MeV10000evrc.root", "par500MeV10000evrc.root", "digi500MeV10000evrc_DT01.root", 0.1)' &
#nohup root -l -q 'NeuRad_digi_DT.C(10000,"sim500MeV10000evrc.root", "par500MeV10000evrc.root", "digi500MeV10000evrc_DT1.root", 1.)' &
#nohup root -l -q 'NeuRad_digi.C(10000,"sim500MeV10000evrc.root", "par500MeV10000evrc.root", "digi500MeV10000evrc.root")' &
#nohup root -l -q 'NeuRad_sim.C(1000, "sim1.root", "par.root", 0.5)' &
#nohup root -l -q 'NeuRad_digi.C(1000,"sim.root", "par.root", "digi.root")' &

#nohup root -l -q 'NeuRad_digi.C(100000,"sim400MeV100000evrc.root", "par400MeV100000evrc.root", "digi400MeV100000evrc_woj.root")' &

#nohup root -l -q 'NeuRad_digi.C(100000,"sim500MeV100000evrc.root", "par500MeV100000evrc.root", "digi500MeV100000evrc_woj.root")' &