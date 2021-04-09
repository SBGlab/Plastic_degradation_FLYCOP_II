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
import plasticDegradationSimulationEstrategy3

biomass1 = float(sys.argv[1])
biomass2 = float(sys.argv[2])
gly1= float(sys.argv[3]) ## gly a gly1; int() a float()
gly2 = float(sys.argv[4]) ## tpha a a gly2; int() a float()
fitness = sys.argv[5]

plasticDegradationSimulationEstrategy3.Flycop(gly2,biomass1,biomass2,gly1,fitness,'./IndividualRunsResults/',2)
