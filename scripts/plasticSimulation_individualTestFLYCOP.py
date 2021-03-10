#!/usr/bin/python3

############ FLYCOP ############
# Author: Beatriz García-Jiménez
# April 2018
################################

# Running an individual test, for a particular consortium configuration given by arguments
#cp -p -R ecoliLongTerm_TemplateOptimizeConsortiumV0 ecoliLongTerm_TestXX
#cd ecoliLongTerm_TestXX
#python3 ../../Scripts/ecoliLongTerm_individualTestFLYCOP.py -10 -16 -11 -12 -6 -16 'Yield'

import sys
import importlib
sys.path.append('../../Scripts')
import plasticDegradationSimulationEstrategy3 as flycop

biomass1 = float(sys.argv[1])
biomass2 = float(sys.argv[2])
gly=int(sys.argv[3])
tpha = int(sys.argv[4])
fitness = sys.argv[5]

flycop.Flycop(tpha,biomass1,biomass2,gly,fitness,'./IndividualRunsResults/',2)
