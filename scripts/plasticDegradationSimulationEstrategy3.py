#!/usr/bin/python3

import cobra 
import pandas as pd
import re
import sys
import getopt
import os.path
import copy
import csv
import math
import cobra.flux_analysis.variability
import subprocess
import shutil, errno
import statistics
from cobra import Reaction
import cometspy as c
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'figure.max_open_warning': 0})

def initialize_models(file, mutations=False, percentage=None):
    for i in file:
        j=1
        #check if model is in format .json
        if model.endswith(".json"):
            cobra_model=cobra.io.json.load_json_model(model)
            cobra.io.sbml.write_cobra_model_to_sbml_file(cobra_model, "model_sbml_tmp")
            sbml_model=cobra.io.read_sbml_model("model_sbml_tmp.sbml") #check this
            #establish layout with a file with the model
            final_model=c.model(sbml_model)
            final_model.id='model_tmp'+str(j)
            final_model.write_comets_model()
            del(cobra_model, sbml_model, final_model)
        else:
            #Assuming the model is in .sbml format
            final_model=c.model(file)
            final_model.id='model_tmp'+str(j)
            final_model.write_comets_model()
            del(final_model)
            
        if mutations == True:
            reactions= final_model.get_reaction_names()
            total_reaction_changes= round(len(reactions)*percentage)
            for x in range(0,len(total_reaction_changes)):
                reaction_changed=random.choice(reactions)
                #deactivate that reaction
                final_model.change_bounds(reaction_changed, 0, 0)
            final_model.id='model_tmp'+str(j)
            final_model.write_comets_model()
            del(final_model)
        j = j+1
        
def model_modifications(tpha, model_id, gly=None):
    # Single GEMs parameter modifications

    # 1.1.- Establish modifications in model  
    model1=c.model('/home/iodmc/Documents/Plastic_degradation_FLYCOP_II/scripts/plasticDegradationEstrategy3-dPCA.cmd')
    #You can change the bounds of a reaction using the change_bounds(reaction name, lower bound, upper bound) method
    #You can create certain conditions for your simulation
    model1.change_bounds('EX_glycol_e', gly, 0)
    model1.change_bounds('EX_tpha_e',tpha,0)
    model1.change_bounds('EX_o2_e',-18.5,0)

    model1.id=model_id[0]
    model1.write_comets_model()
    del(model1)
     # 1.2.- Establish modifications in model 2
    model2=c.model('/home/iodmc/Documents/Plastic_degradation_FLYCOP_II/scripts/defaultModel-3.cmd')
    model2.change_bounds('EX_glycol_e', gly, 0)
    model2.change_bounds('EX_o2_e',-18.5,0)
    model2.change_bounds('EX_34dhbz_e',-5,0)
    
    model2.id=model_id[1]
    model2.write_comets_model()        
    del(model2)

    return (print('Modifications of the model done!'))
    
    
def initialize_layout(filename, biomass1, biomass2):
    #The layout is a file with the stablished format
    layout=c.layout(filename)
    layout.initial_pop=[[0.0,0.0,float(biomass1),float(biomass2)]]
    return layout

def initialize_params(package, globall):
    """Function to initialize the comets' params class
    it can be initialize by submitting two files, one for the package parameters
    and one for the global ones.
    If you don't submit a file, the params class will be initialize with the values stated below
    which have been tested in this simulation"""
    
    if package and globall is not None:
        params = c.params(global_params = globall, package_params= package)
    elif package is not None:
        params = c.params(package_params=package)
    elif globall is not None:
        params = c.params(global_params=globall)
        
    else:
        params = c.params()
        params.all_params['maxCycles']=322
        params.all_params['timeStep']=0.1
        params.all_params['spaceWidth']=0.05
        params.all_params['allowCellOverlap']= True
        params.all_params['deathRate']= 0.01
        params.all_params['numRunThreads']= 8
        params.all_params['maxSpaceBiomass']= 100
        params.all_params['defaultVmax']=20
        params.all_params['showCycleTime']=True
        params.all_params['writeTotalBiomassLog']=True
        params.all_params['writeMediaLog']=True
        params.all_params['writeFluxLog']=True
        params.all_params['useLogNameTimeStamp']=False
        params.all_params['FluxLogRate']=1
        params.all_params['MediaLogRate']=1
        params.all_params['exchangestyle']='Standard FBA'

        
    #check if param 'maxSpaceBiomass' has default value
    if (params.all_params['maxSpaceBiomass']== 0.1):
        print('The parameter "maxSpaceBiomass" is set to 0.1.\n'               'It may need to change if the mo used growths well.')

    return params

def make_df_and_graph(strains, metabolites, comets, layout, params):
    '''This function creates a figure and saves it to pdf format.
    It also creates the file biomass_vs_met.txt which contais the quantity
    of each strain and metabolite and has the following columns:
    time(h), strain1 ... strainX, met1 ... metX.'''
    file_name='_'.join(metabolites)
    df = comets.media #We get the media composition results'
    df2=comets.total_biomass #We get the biomass results
    columns=['cycle']
    for i in range(0,len(strains)):
        columns.append(strains[i])
    df2.columns=columns
    max_cycles=params.all_params['maxCycles']
    """For each metabolite a column with all zeros is added to the dataframe and each row that contains a value
     (metabolite concentration)is changed in the dataframe"""
    for d in metabolites:
        met =df[df['metabolite'] == d]

        temp=np.zeros(max_cycles+1) #Create an array with all zeros
        df2[d]=temp #We added it to the dataframe
        j=1
        while j < (max_cycles+1): #For each cycle
            if (met.cycle==j).any(): #If the row exists
                df2.loc[j-1,d] = float(met[met['cycle']==j]['conc_mmol']) #Its dataframe value is changed
            j+=1  
    np.savetxt(r'biomass_vs_'+file_name+'_template.txt', df2.values, fmt='%s',delimiter='\t') #The data is saved
    #---Starting the figure 
    plt.ioff()
    fig, ax = plt.subplots()

    ax.set_xlabel('time (h)')
    ax.set_ylabel('biomass (g/L)')
    c=['k', 'm', 'b', 'g', 'r']
    j=0
    for i in strains:
        ax.plot(df2['cycle']*0.1, df2[i], label=i, color=c[j])
        j+=1
    ax2 = ax.twinx()
    ax2.set_ylabel('metabolite conc (mM)')
    for m in metabolites:
        ax2.plot(df2['cycle']*0.1, df2[m], label=m)

    handles, labels = ax.get_legend_handles_labels()
    handle_list, label_list = [], []

    for handle, label in zip(handles, labels):
        if label not in label_list:
            handle_list.append(handle)
            label_list.append(label)
    handles, labels = ax2.get_legend_handles_labels()
    for handle, label in zip(handles, labels):
        if label not in label_list:
            handle_list.append(handle)
            label_list.append(label)
    plt.legend(handle_list, label_list)
    
    #saving the figure to a pdf
    plt.savefig('biomass_vs_'+file_name+'_template_plot.pdf')
    
    return df2

def end_simulation_cycle(df2,strains,params):
    """function that stablishes the endCyle after a certain stop condition
    If the initial biomass of a microorganism needs to be used, you can access it with
    float(df2.at[0,'Ecoli1']), where df2 is the dataframe created in make_df_and_graph() function.
    The following is a example of a stop condition in which the cycle is stopped when either:
    glucose and acetate are exhausted or when glucose is exhausted and the microorganisms
    can't use acetate to grow.
    The stop condition should be suitable for your simulation.
    """
    iniBiomass=float(df2.at[0,strains[0]])+float(df2.at[0, strains[1]])
    totTpha=float(df2.at[0,'tpha_e'])
    endCycle=params.all_params['maxCycles'] #total of cycles
        
       
    return iniBiomass, endCycle, totTpha

def biomass_yield(df2,endCycle, strains):
    # Function that compute final biomass as the maximum biomass of each strain.
    
    finalBiomass={}
    for i in range(1, len(strains)+1):
        finalBiomass['finalBiomass'+str(i)]=0
        count=0
        for row in df2.iterrows():
            if float(row[1][i])>finalBiomass['finalBiomass'+str(i)]:
                finalBiomass['finalBiomass'+str(i)]=float(row[1][i])
            if count > endCycle:
                break
            count+=1
        total_finalBiomass=sum(finalBiomass.values())
    return total_finalBiomass
    

def compute_fitness(biomassYieldNew, finalBiomass, iniBiomass, fitFunc, fitTime, comets,df2,params):

    if(fitFunc=='Yield'):
        phaProduction=float(df2.iloc[-2]['C80aPHA_e'])
        fitness= phaProduction
    elif(fitFunc=='MaxYield_MinTime'):
        fitness=0.5*(biomassYieldNew)+0.5*fitTime #Normalizing yield
    elif(fitFunc=='YieldNewScattered'): # (biomass^4)*10: To spread values from ~0.45-0.55 values to 0.5 to 1
        fitness=(biomassYieldNew**4)*10
    elif(fitFunc=='MaxYieldNewScattered_MinTime'):
        fitness=0.5*((biomassYieldNew**4)*10)+0.5*fitTime
    elif(fitFunc=='Biomass'):
        fitness=float(finalBiomass-iniBiomass)
    elif(fitFunc=='MaxBiomass_MinTime'):
        fitness=0.5*(float(finalBiomass-iniBiomass))+0.5*fitTime
    elif((fitFunc=='GR')or(fitFunc=='MaxGR_MinTime')):
        flux = comets.fluxes
        gr=[]
        for m in range(1,len(layout.models)+1):
            model_tmp = c.model("model_tmp"+str(m)+".cmd")
            objective = model_tmp.objective
            gr.append(flux.iat[47-m-1,int(objective+3)])
        fitGR=sum(gr)/len(gr)
    if(fitFunc=='GR'):
        fitness=fitGR
    elif(fitFunc=='MaxGR_MinTime'):
        fitness=0.5*fitGR+0.5*fitTime
    return fitness

def calculate_uptake(comets,layout, met):
    uptakeMet=[]
    flux=comets.fluxes
    #for i in range(1, len(layout.models)+1):
     #   j=0
        #Compute metabolite uptake
    model=layout.models[0]
    df3 = model.reactions
    numRxnMet=int(df3.loc[df3['REACTION_NAMES']=='EX_'+str(met)]['ID'])
    uptakeMet.append(float(flux.iat[0+2,numRxnMet+3]))
        #j+=1
    return uptakeMet


def Flycop(tpha,biomass1,biomass2,gly,fitFunc='Yield', dirPlot='', repeat=10, params_package=None, params_global=None):
    """Primarily function of the simulation
    You need to add to it the parameters to adjust.
    You need to fill these variables below(layout, model_id, model, wgMet, metabolites, strains) with your information
    layout= The layout for the comets simulation
    model_id = the names of the models in your layout
    model= the name of the model yours modifications are from
    wgMet=molecular weight of the metabolite to compute fitness. Must be divided by 1000.
    metabolites= list with the names of the metabolites to study
    strains = list with the names of the strains you're working with
    simulationID = will be part of the filename for the results 
    metUptake = the metabolite you want to calculate the uptake of"""
    layout_file='/home/iodmc/Documents/Plastic_degradation_FLYCOP_II/scripts/layout_3.txt'
    model_id=['noPCADeg_tmp','default3_tmp']
    wgMet=0.16411
    metabolites=['tpha_e','glycol_e','34dhbz_e','C80aPHA_e','o2_e']
    strains=['Pputida_no-pca', 'Pputida_default-3']
    metUptake = metabolites[0]
    
    if not(os.path.exists('IndividualRunsResults')):
        os.makedirs('IndividualRunsResults')
    if not(os.path.exists('/home/iodmc/Documents/Plastic_degradation_FLYCOP_II/scripts/noPCADeg_tmp.cmd')):
        
        model_modifications(tpha, model_id, gly)
        
    totfitness=0
    sumTotBiomass=0
    sumTotYield=0
    fitnessList=[]
    
    layout= initialize_layout(layout_file, biomass1, biomass2)
   
    params =initialize_params(params_package, params_global)
   
    #establish comets class
    comets = c.comets(layout, params)
    """You may need to adjust the comets.set_classpath object if
    your programs aren't in the default comets location.
    Changing using the comets method set_classpath('name of the program', 'path')"""
    comets.set_classpath('lang3', '/home/iodmc/comets/lib/commons-lang3-3.9/commons-lang3-3.9.jar')
    comets.set_classpath('gurobi', '/opt/gurobi911/linux64/lib/gurobi.jar')
    comets.set_classpath('bin', '/home/iodmc/comets/bin/comets_2.10.5.jar')

    # To repeat X times, due to random behaviour in COMETS:
    for i in range(repeat):
        comets.run(delete_files=True)
        
        #obtaining the results and writing them to csv
        comets.total_biomass.to_csv('Total_biomass_log_template.txt',sep='\t',index=False)
        comets.fluxes.to_csv('Flux_log_template.txt',sep='\t', index=False)
        comets.media.to_csv('Media_log_template.txt',sep='\t', index=False)        
        
        #make the graphic
        df2=make_df_and_graph(strains, metabolites[0:-1], comets, layout, params)
        
        # 7.- Compute fitness (measure to optimize):
        print('computing fitness...')
        
        # 7.1.- Determine endCycle: when glucose and acetate are exhausted
        iniBiomass, endCycle, totMet = end_simulation_cycle(df2, strains,params)
        
        # 7.2.- Compute first element fitness: maximize biomass yield
        finalBiomass =biomass_yield(df2,endCycle, strains)
        if wgMet is not None:
            
            biomassYieldNew=float((finalBiomass-iniBiomass)/(totMet*wgMet)) # molecular weigth of met per mmol
        
        
        # 7.3.- Compute second element fitnes: minimize time        
        fitTime=1-(float(endCycle)/float(params.all_params['maxCycles']))
        
        # 7.4.- Compute joint fitness, as a 50% each element.
        fitness= compute_fitness(biomassYieldNew, finalBiomass, iniBiomass, fitFunc, fitTime, comets,df2,params)
        
        # 7.5.- Calculate the uptake of metabolite
        uptakeMet =calculate_uptake(comets,layout,metUptake)
            
            
        print(" Total biomass: "+str(round(finalBiomass,6))+" infitness cycle "+str(endCycle)+". Biomass yield="+str(round(biomassYieldNew,6)))
   
        totfitness=totfitness+fitness
        fitnessList.append(fitness)
        sumTotBiomass=sumTotBiomass+finalBiomass
        sumTotYield=sumTotYield+biomassYieldNew
            
        # Copy individual solution
        x= '_'.join(metabolites[0:-1])
        y=str(uptakeMet[0])
        file='IndividualRunsResults/'+'biomass_vs_'+x+'_template_plot.pdf'        
        shutil.copy('biomass_vs_'+x+'_template_plot.pdf',file)
        if(dirPlot != ''):
            file2=dirPlot+'biomass_vs_tpha_'+'_run'+str(i)+'_'+str(round(fitness,6))+'_'+str(biomass1)+'_'+str(biomass2)+'_'+str(gly)+'_'+str(tpha)+'.pdf'
            shutil.copy(file,file2)
        file='IndividualRunsResults/'+'total_biomass_log_run'+str(i)+'.txt'
        shutil.move('Total_biomass_log_template.txt',file)
        file='IndividualRunsResults/'+'media_log_run'+str(i)+'.txt'
        shutil.move('Media_log_template.txt',file)
        file='IndividualRunsResults/'+'flux_log_run'+str(i)+'.txt'
        shutil.move('Flux_log_template.txt',file) 
            
    avgfitness=totfitness/repeat
    sdfitness=statistics.stdev(fitnessList)
    avgBiomass=sumTotBiomass/repeat
    avgYield=sumTotYield/repeat
    
    print("Fitness_function\tconfiguration\tfitness\tsd\tavg.Biomass\tavg.Yield\tendCycle")
                    
    print(fitFunc+"\t"+str(tpha)+','+str(gly)+','+str(biomass1)+','+str(biomass2)+','+"\t"+str(round(avgfitness,6))+"\t"+str(sdfitness)+"\t"+str(round(avgBiomass,6))+"\t"+str(round(avgYield,6))+"\t"+str(endCycle))
    with open(dirPlot+"configurationsResults"+fitFunc+".txt", "a") as myfile:
        myfile.write("Fitness_function\tconfiguration\tfitness\tsd\tavg.Biomass\tavg.Yield\tendCycle\n")
        myfile.write(fitFunc+"\t"+str(tpha)+','+str(gly)+','+str(biomass1)+','+str(biomass2)+','+"\t"+str(round(avgfitness,6))+"\t"+str(sdfitness)+"\t"+str(round(avgBiomass,6))+"\t"+str(round(avgYield,6))+"\t"+str(endCycle)+"\n")
  
    print("Avg.fitness(sd):\t"+str(avgfitness)+"\t"+str(sdfitness)+"\n")
    if(sdfitness>0.1):
        avgfitness=0.0
    #tmp models needs to be deleted
    for i in model_id:
        try:
            os.remove(str(i)+'.cmd')
        except FileNotFoundError:
            pass

    return avgfitness,sdfitness
