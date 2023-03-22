# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 15:12:46 2022

@author: Pierre Jouannais, Department of Planning, DCEA, Aalborg University
pijo@plan.aau.dk




"""


'''
Script which calls the functions to run simulations and produce the final results.
1) Initializes the variables needed for the simulations.
2) Calculate MonteCarlo LCA impacts for all activites in the foreground database
3) Simulate total LCAs in all scenarios

Execute the whole script to export the results of stochastic LCAs in all scenarios

Must be executed from the folder "Scripts Fish"

Takes a few minutes if several CPUS (remote servers).
Takes 3 hours on local computer for 1000 iterations
'''


size=2# How many MonteCarlo iterations?



import datetime
from time import *
import requests
import pickle
import cProfile
from scipy.integrate import odeint

import os
import sys


import pandas as pd
import decimal
from random import *
import pstats
from itertools import *
from math import*
import csv
import copy
import numpy as np
import seaborn as sns
import random


import bw2data
import bw2io
from bw2data.parameters import *
import brightway2 as bw


from SALib.test_functions import Ishigami
import math
from SALib.sample import saltelli
from SALib.sample import fast_sampler
from SALib.analyze import sobol
from SALib.analyze import fast
from SALib.sample import latin
from SALib.analyze import delta

import SALib

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits import mplot3d


#import Cultivation_simul_Night_Harvest_1_2nd_faster_noreinjection as cultsimul
import Functions_for_physical_and_biological_calculations_3rd as functions
import Main_functions_fish_scenarios_FCRs as mainfunc


import ray


import scipy
import numpy as np
from scipy import linalg

from ast import literal_eval



from itertools import compress





# Accessory functions

def createFolder(directory):
    '''Creates a folder/directory with the path given as input'''
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)


def export_pickle(var, name_var):
    '''Saves a pickle in the working driectory and
    saves the object in input across python sessions'''

    with open(name_var+'.pkl', 'wb') as pickle_file:
        pickle.dump(var, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)

def export_pickle_2(var, name_var, namefolder_in_root):
    '''Saves a pickle in the working directory and
    saves the object in input across python sessions'''

    path_object = "../"+namefolder_in_root+"/"+name_var+".pkl"
    with open(path_object, 'wb') as pickle_file:
        pickle.dump(var, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)



def importpickle(path):
    with open(path, 'rb') as pickle_load:
        obj = pickle.load(pickle_load)
    return obj    



"""1) INITIALIZE"""

# Uploading the necessary csv files



# Fish feed composition

fishfeed_table_withNP = pd.read_csv("../Data/Ingredient_composition_withNP.csv", sep=";",
                             header=0, encoding='unicode_escape', engine='python')

# Cleaning

fishfeed_table_withNP = fishfeed_table_withNP[0:8][[
    'Ingredient', 'kg.kg feed-1', 'Lipid', 'Protein', 'Carb', 'Ash', 'Water',"N","P"]]



# Calculating current fish feed biochmical profile

biochem_profile_feed_incumbent = functions.biochemprofile(fishfeed_table_withNP)


N_P_profile_feed_incumbent = [sum([fishfeed_table_withNP["N"][a] * fishfeed_table_withNP['kg.kg feed-1'][a] for a in range(len(fishfeed_table_withNP["N"]))]),
                              sum([fishfeed_table_withNP["P"][a] * fishfeed_table_withNP['kg.kg feed-1'][a] for a in range(len(fishfeed_table_withNP["P"]))])]
    
ingredient_profile_incumbent = list(fishfeed_table_withNP['kg.kg feed-1']) 



# Elemental composition of macronutrients

elemental_contents = pd.read_csv("../Data/elemental_contents.csv",
                                 sep=";",
                                 header=0,
                                 encoding='unicode_escape',
                                 engine='python')

elemental_contents = elemental_contents.iloc[:, 1:]



digestibility_list = [0.95,0.92,0.71,0.5,0.5] # lip, prot, carb, ash, phosphorus



# Import foreground Technosphere Fish Production


Techno_Matrix_Fish_loaded = pd.read_excel('../Data/Technosphere_Fish_17_03.xlsx',header=None)


Techno_Matrix_Fish_with_names =  np.array(Techno_Matrix_Fish_loaded)

Techno_Matrix_Fish_activities = Techno_Matrix_Fish_with_names[0,1:]

Techno_Matrix_Fish_products = Techno_Matrix_Fish_with_names[1:,0]

# Only amounts, no names
Techno_Matrix_Fish = Techno_Matrix_Fish_with_names[1:,1:]

# Quick test to see if the froeground technosphere looks ok
if Techno_Matrix_Fish.shape[0]!=Techno_Matrix_Fish.shape[1]:
    
     sys.exit("The Technosphere matrix is not square") 



# Keep the names and the positions of the activities which are ONLY connected to background activites
# This means excluding activites starting with "FG_"

filter_names_ok_activities_fish = [name[0:3] !="FG_" for name in Techno_Matrix_Fish_activities]
filter_names_exclude = [name[0:3] =="FG_" for name in Techno_Matrix_Fish_activities]


activities_fish_background = list(compress(Techno_Matrix_Fish_activities, filter_names_ok_activities_fish))

activities_fish_foreground = list(compress(Techno_Matrix_Fish_activities, filter_names_exclude))

positions_background = [np.where(Techno_Matrix_Fish_activities == act)[0][0] for act in activities_fish_background]

positions_foreground = [np.where(Techno_Matrix_Fish_activities == act)[0][0] for act in activities_fish_foreground]


# Collec the positions of some specific activities
index_growth_stages= positions_foreground[:-2]


names_growth_stages_no_solid_filter_laguna = ['FG_Seafarm_DK_1','FG_Seafarm_DK_2']

index_growth_stages_no_filter_laguna = [np.where(Techno_Matrix_Fish_activities == act)[0][0] for act in names_growth_stages_no_solid_filter_laguna]

index_roe = np.where(Techno_Matrix_Fish_activities == "Roe_prod")[0][0]

index_feed = np.where(Techno_Matrix_Fish_activities == "Fish_feed")[0][0]


index_dead_fish = np.where(Techno_Matrix_Fish_activities == "FG_Dead_fish_management_supply")[0][0]


index_Nemissions = np.where(Techno_Matrix_Fish_activities == "N_emissions_background")[0][0]

index_Pemissions = np.where(Techno_Matrix_Fish_activities == "P_emissions_background")[0][0]

index_sludge = np.where(Techno_Matrix_Fish_activities == "FG_Sludge management")[0][0]

index_growing_DK = np.where(Techno_Matrix_Fish_activities == "FG_Growing_out_DK")[0][0]

index_300g = np.where(Techno_Matrix_Fish_activities == "FG_Growing_out_IT_2_para")[0][0]



index_biogas_updgrade = np.where(Techno_Matrix_Fish_activities == "Heat market substitution Fish")[0][0]
index_N_substitution = np.where(Techno_Matrix_Fish_activities == "N source production fish")[0][0]
index_P_substitution = np.where(Techno_Matrix_Fish_activities == "P source production fish")[0][0]
index_heat_substitution = np.where(Techno_Matrix_Fish_activities == "Heat market substitution Fish")[0][0]

index_drug1=np.where(Techno_Matrix_Fish_activities == "Drug_1_prod")[0][0]
index_drug2=np.where(Techno_Matrix_Fish_activities == "Drug_2_prod")[0][0]
index_drug3=np.where(Techno_Matrix_Fish_activities == "Drug_3_prod")[0][0]
index_drug4=np.where(Techno_Matrix_Fish_activities == "Drug_4_prod")[0][0]



" Positions of the growth stages in the tehcnosphere matrix"
dict_correspondance_techno_growth_stagenames={0:"HAIT",
                                              1:"FRIT",
                                              5:"GOIT2_para",
                                              6:"GOIT1bis",
                                              7:"GODK",
                                              8:"SFDK1",
                                              9:"SFDK2"}
# Main production line
index_growth_stages_to_modif = [0, 1, 6, 7, 8, 9]



# Prepare demand vector for the LCA
demand_vector= np.zeros(Techno_Matrix_Fish.shape[0])

index_FU = np.where(Techno_Matrix_Fish_activities == "FG_Seafarm_DK_2")[0][0]

demand_vector[index_FU] = 1





                             
                             



def initial(Techno_Matrix_Fish,
            demand_vector,
            index_300g,
            index_growing_DK,
            index_roe,
            index_FU,
            index_feed,
            index_dead_fish):
    """Function to calculate the current economic indicators,
    before modification"""
    
    # Before modification
    HAIT_loss_lev=0  
    FRIT_loss_lev=0                        
    GOIT1_loss_lev=0                                                     
    GOIT2_loss_lev=0                                  
    GOIT1bis_loss_lev=0
    GODK_loss_lev=0
    SFDK1_loss_lev= 0                                
    SFDK2_loss_lev= 0                                 
    SFDK1_loss_red= 0 
    SFDK2_loss_red=0
    FRIT_loss_red=0
    GODK_loss_red=0
    GOIT2_loss_red =0   
    HAIT_loss_red=0
    GOIT1bis_loss_red=0     
    HAIT_para_loss_lev=0
    FRIT_para_loss_lev=0
    GOIT1_para_loss_lev=0
    GOIT2_para_loss_lev=0
    FRIT_micro_dose= 0 
    GODK_micro_dose= 0                                                                    
    GOIT2_micro_dose= 0                                
    GOIT1_micro_dose= 0                   
    GOIT1bis_micro_dose= 0                         
    HAIT_micro_dose= 0                                  
    SFDK2_micro_dose= 0                                  
    SFDK1_micro_dose= 0
    HAIT_FCR_red_ratio=1
    FRIT_FCR_red_ratio=1
    GOIT1_FCR_red_ratio=1
    GOIT2_FCR_red_ratio=1
    GOIT1bis_FCR_red_ratio=1
    GODK_FCR_red_ratio=1
    SFDK1_FCR_red_ratio=1
    SFDK2_FCR_red_ratio=1
    ratio_CO2_CH4_biogas_fish = 1
    CH4_volume_biogas_fish = 1
    
    # Just initialize with dummy values for these parameters here

    P_in_fish= 0.00415  # kg
    N_in_fish= 0.021627 # kg
    K_in_fish=0.004618313
    Mg_in_fish=0.00038425
    water_in_fish= 0.7
    gvs_gts_in_fish= 0.89
    
    gvs_gts_in_fish= 0.89
    
    
    
    
    biogenic_CO2_emitted = 10
    amount_of_combustion_natural_gas = 23
    m3_natural_gas_substituted = 9
    amount_of_upgrading_act = 2

    excretion_N_removal_efficiency = 1
    excretion_P_removal_efficiency = 1
    solid_filter_efficiency =1

    fertilizer_substi_digest_P=1
    fertilizer_substi_digest_N=1
    fertilizer_substi_digest_K=1
    fertilizer_substi_digest_Mg=1
    
    MJ_substituted_per_kilo_deadfish = 1
    
    
    # Initialize the numerical Tehcnosphere
    
    Techno_fish_num=np.zeros((Techno_Matrix_Fish.shape[0],Techno_Matrix_Fish.shape[0]))
    
    for i in range(Techno_fish_num.shape[0]):
    
        for j in range(Techno_fish_num.shape[1]):
            
            if type(Techno_Matrix_Fish[i,j])==str: # Then we evaluate the string
            
                Techno_fish_num[i,j] = eval(Techno_Matrix_Fish[i,j])
                
            else:        # Then we just take the float
                Techno_fish_num[i,j] = Techno_Matrix_Fish[i,j]
    


    # Calculate economic indicators
    
    [Tech_Matrix,
     eco_FCR_0, 
     bio_FCR_0,
     Loss_ratio_INC_0, 
     ratio_loss_biological_INC_0] = mainfunc.calculate_economic_indicators(Techno_fish_num,
                            index_300g,
                            index_growing_DK,
                            index_roe,
                            index_FU,
                            index_feed,
                            index_dead_fish) 
    


    return Tech_Matrix,eco_FCR_0,bio_FCR_0, Loss_ratio_INC_0, ratio_loss_biological_INC_0


[Tech_Matrix,ECO_FCR_0,bio_FCR_0, Loss_ratio_INC_0, ratio_loss_biological_INC_0]=initial(Techno_Matrix_Fish,
            demand_vector,
            index_300g,
            index_growing_DK,
            index_roe,
            index_FU,
            index_feed,
            index_dead_fish)

# Current dose for drugs

supply_ini=linalg.solve(Tech_Matrix, demand_vector)

ini_dose_drug1 = supply_ini[index_drug1]
ini_dose_drug2 = supply_ini[index_drug2]
ini_dose_drug3 = supply_ini[index_drug3]
ini_dose_drug4 = supply_ini[index_drug4]





""" Initialize brightway and calculate LCIA for 1 unit of each technosphere input """


"Fish_test_17_03"

bw.projects.set_current('Fish_project_check') 



#bw.databases
# Loading Ecoinvent
Ecoinvent = bw.Database('ecoinvent 3.8 conseq')

Ecoinvent.random()

# Fix overestimated uncertainty in ecoinvent for mercury flow in wheat production
act_problem=bw.Database("ecoinvent 3.8 conseq").get('28c279b19c7d0e0793614e3012843056')
for exc in list(act_problem.exchanges()):
    if exc["type"]=="biosphere" and exc["name"]=="Mercury":
        #print(exc['uncertainty type'])
        exc['uncertainty type']=0
        exc.save()
act_problem.save()   


# Loading foreground database

FISHMIC = bw.Database('AH_combi_1')



biosphere=bw.Database('biosphere3')


# Not used in the current version.
# Keep it anyways to keep the structure
electricity_low_voltage_input_per_m3_wastewater = 0.20571
electricity_high_voltage_input_per_m3_wastewater = 0.0086946




list_meth =[('ReCiPe Midpoint (H)', 'climate change', 'GWP100'), 
            ('ReCiPe Midpoint (H)', 'human toxicity', 'HTPinf'),
            ('ReCiPe Midpoint (H)', 'freshwater ecotoxicity', 'FETPinf'),
            ('ReCiPe Midpoint (H)', 'freshwater eutrophication', 'FEP'),
            ('ReCiPe Midpoint (H)', 'marine eutrophication', 'MEP'),
             ('TRACI', 'environmental impact', 'eutrophication'),
             ('CML 2001 (superseded)', 'eutrophication potential', 'generic'),
             ('ReCiPe Midpoint (H)', 'terrestrial ecotoxicity', 'TETPinf'),
             ('ReCiPe Midpoint (H)', 'fossil depletion', 'FDP'),
             ('ReCiPe Midpoint (H)', 'terrestrial acidification', 'TAP100'),
             ('ReCiPe Midpoint (H)', 'ozone depletion', 'ODPinf'),
             ('ReCiPe Midpoint (H)', 'particulate matter formation', 'PMFP')]


list_meth_code =  []

for meth in list_meth:
    list_meth_code.append(meth[-1])



list_cfs= [bw.Method((meth)).load() for meth in list_meth]

    









# Create dictionnaries with incumbent production outputs for the fish farm growth stages
# Needed to ease the calculations later

# the key is the name of the associated loss level parameter
Dict_incumbent_outputs_growth_stages_loss_level = {'HAIT_loss_lev': 2400,   
                                                   'FRIT_loss_lev': 12000,
                                                   'GOIT1_loss_lev': 64000, 
                                                   'GOIT2_loss_lev': 96000,
                                                   'GOIT1bis_loss_lev': 32000, 
                                                   'GODK_loss_lev': 319822+750100+194750, 
                                                   'SFDK1_loss_lev': 3512877, 
                                                   'SFDK2_loss_lev': 4934867}
                                 

# the key is the name of the associated loss reduction parameter

Dict_incumbent_outputs_growth_stages_loss_red = {'HAIT_loss_red': 2400,   
                                                   'FRIT_loss_red': 12000,
                                                   'GOIT1_loss_red': 64000, 
                                                   'GOIT2_loss_red': 96000,
                                                   'GOIT1bis_loss_red': 32000, 
                                                   'GODK_loss_red': 319822+750100+194750, 
                                                   'SFDK1_loss_red': 3512877, 
                                                   'SFDK2_loss_red': 4934867}




# Create dictionnaries with incumbent losses outputs for the fish farm growth stages


# the key is the name of the associated loss reduction parameter

Dict_incumbent_losses_growth_stages_loss_red = {'HAIT_loss_red': 400,   
                                                   'FRIT_loss_red': 0,
                                                   'GOIT1_loss_red': 16000, 
                                                   'GOIT2_loss_red': 0,
                                                   'GOIT1bis_loss_red': 0, 
                                                   'GODK_loss_red': 31053, 
                                                   'SFDK1_loss_red': 266783, 
                                                   'SFDK2_loss_red': 0}




# the key is the name of the associated loss level parameter
Dict_incumbent_losses_growth_stages_loss_level = {'HAIT_loss_lev': 400,   
                                                   'FRIT_loss_lev': 0,
                                                   'GOIT1_loss_lev': 16000, 
                                                   'GOIT2_loss_lev': 0,
                                                   'GOIT1bis_loss_lev': 0, 
                                                   'GODK_loss_lev': 31053, 
                                                   'SFDK1_loss_lev': 266783, 
                                                   'SFDK2_loss_lev': 0}
                                 

Dict_FCR_bio = {'HAIT_FCR': 2700/(2400+400),   
            'FRIT_FCR': 12400/(12000-2400),
            'GOIT1bis_FCR': 28000/(32000-12000),                  
            'GODK_FCR': 1135494/(319822+750100+194750+31053-162266), 
            'SFDK1_FCR': 2556661/(3512877+266783-1926316), 
            'SFDK2_FCR': 2941339/(4934867+710210-3512877)}
                               





"""2) Calculating MonteCarlo impacts for 1 unit of each foreground activity"""
     
#Collect only foreground activites that are only connected to background 



list_fish_FU = []
for name_act in activities_fish_background:
    
    for act in FISHMIC:
        
        if act["name"] == name_act:
            
            list_fish_FU.append({act:1})
            
drug_inputs_names=['Drug_1_prod','Drug_2_prod','Drug_3_prod','Drug_4_prod','Drug_5_prod']






# Collect Characterization matrixes in a dictionnary

C_matrixes ={}
# Initialize a LCA object, whatever object.
Lca=bw.LCA(list_fish_FU[2],('ReCiPe Midpoint (H)', 'climate change', 'GWP100'))
Lca.lci()
Lca.lcia()

for meth in list_meth:
    Lca.switch_method(meth)
    C_matrixes[meth]=Lca.characterization_matrix

    



# Initalize emplty list of arrays which will contain the MonteCarlo results
list_array_mc_sample_fish = [np.array([[0]*len(list_fish_FU)  ]*size,dtype="float32") for meth in range(len(list_meth))]

list_FU_combined_mc =  list_fish_FU 
list_FU_combined_names_mc = activities_fish_background 


#Initialize a MCLCA object with any FU
mc=bw.MonteCarloLCA(list_fish_FU[2],('ReCiPe Midpoint (I)', 'climate change', 'GWP20')) 


time1=time()
for it in range(size):
    
    print("iteration",it)
    next(mc)

    for i in range(0,len(list_FU_combined_mc)):
            
            
            mc.redo_lcia(list_FU_combined_mc[i])  # redo with new FU
            
            index_array_method=-1
            
            #print(i)
            for m in list_meth:
            
                #print("ok3",m)
    
                index_array_method+=1
            
                # This calculates the impact for 1 input, for 1 meth    
                list_array_mc_sample_fish[index_array_method][it,i]=(C_matrixes[m]*mc.inventory).sum()
                



# Concatenate, reorganize results
list_array_total_mc=[]
for index_meth in range(len(list_meth)):
    concat = np.concatenate([list_array_mc_sample_fish[index_meth]], axis= 1)
    list_array_total_mc.append(concat)
          

list_array_total_mc_sorted =[[[0 for a in range(list_array_total_mc[0].shape[1])] for meth in list_meth] for it in range(size)]


for it in range(size):
    row_to_add = [meth[it] for meth in list_array_total_mc]
    list_array_total_mc_sorted[it] = row_to_add
        

 
time2=time()

# Export as pickle if needed later

x = datetime.datetime.now()

month=str(x.month)
day=str(x.day)
microsec=str(x.strftime("%f"))


name_file_background ='montecarlo_background_fish'+"_"+month+"_"+day+"_"+microsec+"size="+str(size)


export_pickle_2(list_array_total_mc_sorted, name_file_background, "Background_mc_fish")










# Create the parameters dictionnaries corresponding to scenarios



# Description of the parameters are given in the appendix.
# Here values can be changed for parameters with unique values (no distributions)
# Values with distributions will be overwritten.



# Fish farm parameters

# Current primary data

Fish_farm_and_compound_effect_dict_INC= {'HAIT_loss_lev': 0,   
                              
                                 'HAIT_loss_red': 0,
                                 
                                 'FRIT_loss_lev': 0,
                                 
                                 'FRIT_loss_red': 0,
                                 
                                 'GOIT1_loss_lev': 0, 
                                 
                                 'GOIT1_loss_red': 0, 
                                 
                                 'GOIT2_loss_red': 0, 
                                 
                                 'GOIT2_loss_lev': 0,
                                 
                                 'GOIT1bis_loss_lev': 0, 
                                 
                                 'GOIT1bis_loss_red': 0, 
                                                                  
                                 'GODK_loss_lev': 0, 
                                 
                                 'GODK_loss_red': 0, 
                                                                  
                                 'SFDK1_loss_lev': 0, 
                                 
                                 'SFDK2_loss_lev': 0, 
                                 
                                 'SFDK1_loss_red': 0, 
                                 
                                 'SFDK2_loss_red': 0, 
                                 
                                 'FRIT_micro_dose': 0, 

                                 'GODK_micro_dose': 0,  
                                                                  
                                 'GOIT2_micro_dose': 0,
                                 
                                 'GOIT1_micro_dose': 0,
                   
                                 'GOIT1bis_micro_dose': 0,
                                 
                                 
                                 'HAIT_para_loss_lev': 0,   
                              
                                 
                                 'FRIT_para_loss_lev': 0,
                                 
                                 
                                 'GOIT1_para_loss_lev': 0, 
                                 
                                                                  
                                 'GOIT2_para_loss_lev': 0,
                                 

                                 
                                 'HAIT_micro_dose': 0, # kWh.kg dw -1
                                 
                                 'SFDK2_micro_dose': 0, # kWh.kg dw -1
                                 
                                 'SFDK1_micro_dose': 0,
                                 
                                 
                                 "HAIT_FCR_red_ratio_frac":0,
                                 "FRIT_FCR_red_ratio_frac":0,
                                 "GOIT1_FCR_red_ratio_frac":0,
                                 "GOIT2_FCR_red_ratio_frac":0,
                                 "GOIT1bis_FCR_red_ratio_frac":0,
                                 "GODK_FCR_red_ratio_frac":0,
                                 "SFDK1_FCR_red_ratio_frac":0,
                                 "SFDK2_FCR_red_ratio_frac":0,
                                 "ratio_CO2_CH4_biogas_fish":0.465, 
                                 "CH4_volume_biogas_fish":200, # ml
                                 "P_in_fish": 0.00415,  # kg
                                 "N_in_fish": 0.022, # kg 
                                  "K_in_fish":0.0046,
                                  "Mg_in_fish":0.00038,
                                  "water_in_fish": 0.7,
                                  "gvs_gts_in_fish": 0.89,
                                  "fraction_non_ingested" : 0.05,
                                  "excretion_N_removal_efficiency" :0.742095258,
                                  "excretion_P_removal_efficiency" : 0.656726614,
                                  "solid_filter_efficiency": 0.3485,
                                 
    
                                  "CH4_volume_sludge_and_manure": 350,
                                  "share_fish_sludge_in_substrate": 0.33,
                                  "ratio_CO2_CH4_biogas_sludge": 0.465,
                                  "N_manure": 0.05,
                                  "P_manure": 0.01,
                                  "fertilizer_substi_digest_N":0.9,
                                  "fertilizer_substi_digest_P":0.5,
                                  "fertilizer_substi_digest_K": 0.5,
                                  "fertilizer_substi_digest_Mg": 0.5,
                                  

                                  "fertilizer_substi_manure_field_N": 0.54,
                                  "fertilizer_substi_manure_field_P": 0.3}


# Physical parameters
Physicdict = {'CH4_LHV':50} # MJ.m-3




# Distribution dictionnaries
# each parameter is assigned a list containing  :
   # [Distribution,min,max,mode,sd]

# Distribution :
#   - 'unique' if no distribution. The value indicated indicated in
#     the normal static dictionnary will be considered.
#   - 'unif' for uniform, uses min and max
#   - 'triang, uses mim max and mode with mode as a fracion of max-min

#


# Current primary data
Fish_farm_and_compound_effect_dict_distributions_INC= {'HAIT_loss_lev': ['unique', [0, 0, 0.15, 0, 0]],   
                              
                                 'HAIT_loss_red': ['unique', [0, 0, 1, 0, 0]],
                                 
                                 'FRIT_loss_lev': ['unique', [0, 0, 0.15, 0, 0]],
                                 
                                 'FRIT_loss_red': ['unique', [0, 0, 1, 0, 0]],
                                 
                                 'GOIT1_loss_lev': ['unique', [0, 0, 0, 0, 0]],  # Never modified
                                 
                                 'GOIT1_loss_red': ['unique', [0, 0, 0, 0, 0]],  # Never modified
                                 
                                 'GOIT2_loss_red': ['unique', [0, 0, 0, 0, 0]],  # Never modified
                                 
                                 'GOIT2_loss_lev': ['unique', [0, 0, 0, 0, 0]],  # Never modified
                                 
                                 'GOIT1bis_loss_lev': ['unique', [0, 0, 0.15, 0, 0]],   # Never modified
                                 
                                 'GOIT1bis_loss_red': ['unique', [0, 0, 1, 0, 0]],  # Never modified
                                                                  
                                 'GODK_loss_lev': ['unique', [0, 0, 0.15, 0, 0]], 
                                 
                                 'GODK_loss_red': ['unique', [0, 0, 1, 0, 0]], 
                                 
                                 'SFDK1_loss_lev': ['unique', [0, 0, 0.15, 0, 0]], 
                                 
                                 'SFDK2_loss_lev': ['unique', [0, 0, 0.15, 0, 0]], 
                                 
                                 'SFDK1_loss_red': ['unique', [0, 0, 1, 0, 0]], 
                                 
                                 'SFDK2_loss_red': ['unique', [0, 0, 1, 0, 0]], 
                                 
                                 'FRIT_micro_dose': ['unique', [0, 0.00000017, 0.0047, 0, 0]], 

                                 'GODK_micro_dose': ['unique', [0, 0.00000017, 0.0047, 0, 0]],    
                                                                  
                                 'GOIT2_micro_dose': ['unique', [0, 0, 0, 0, 0]], # Always 0
                                 
                                 'GOIT1_micro_dose': ['unique', [0, 0, 0, 0, 0]],  # Always 0
                   
                                 'GOIT1bis_micro_dose': ['unique', [0, 0.00000017, 0.0047, 0, 0]], # Always 0
                                 
                                 
                                 'HAIT_para_loss_lev': ['unique', [0, 0, 0, 0, 0]],   
                              
                                 
                                 'FRIT_para_loss_lev': ['unique', [0, 0, 0, 0, 0]],
                                 
                                 
                                 'GOIT1_para_loss_lev': ['unique', [0, 0, 0, 0, 0]], 
                                 
                                                                  
                                 'GOIT2_para_loss_lev': ['unique', [0, 0, 0, 0, 0]],
                                 
                                 
                                 'HAIT_micro_dose': ['unique', [0, 0.00000017, 0.0047, 0, 0]], 
                                 
                                 'SFDK2_micro_dose': ['unique', [0, 0.00000017, 0.0047, 0, 0]], 
                                 
                                 'SFDK1_micro_dose': ['unique', [0, 0.00000017, 0.0047, 0, 0]], 
                                 
               
                                 "HAIT_FCR_red_ratio_frac":['unique', [0, 0, 1, 0, 0]],
                                 "FRIT_FCR_red_ratio_frac":['unique', [0, 0, 1, 0, 0]],
                                 "GOIT1_FCR_red_ratio_frac":['unique', [0,  0, 0, 0, 0]], # Never modified
                                 "GOIT2_FCR_red_ratio_frac":['unique', [0,  0, 0, 0, 0]],  # Never modified
                                 "GOIT1bis_FCR_red_ratio_frac":['unique', [0,  0, 1, 0, 0]],  # Never modified
                                 "GODK_FCR_red_ratio_frac":['unique', [0, 0, 1, 0, 0]],
                                 "SFDK1_FCR_red_ratio_frac":['unique', [0, 0, 1, 0, 0]],
                                 "SFDK2_FCR_red_ratio_frac":['unique', [0, 0, 1, 0, 0]],
                                 
                                 "ratio_CO2_CH4_biogas_fish":['unique', [0.465, 0.6, 0.9, 0, 0]], 
                                 "CH4_volume_biogas_fish": ['triang', [0, 390, 920, 550, 0]]  , # ml
                                 "P_in_fish":['unique', [0.00415, 0.0026, 0.005513, 0, 0]] ,
                                 "N_in_fish":['unique', [0.021627, 0, 0, 0, 0]]  ,
                                  "K_in_fish": ['unique', [0.004618313, 0.002954, 0.005513, 0, 0]]  ,
                                  "Mg_in_fish": ['unique', [0.00038425, 0.000301, 0.000431, 0, 0]],
                                  "water_in_fish": ['unique', [0.7, 0, 0, 0, 0]],
                                  "gvs_gts_in_fish": ['unique', [0.89, 0, 0, 0, 0]],
                                  "fraction_non_ingested" : ['unif', [0.05, 0.00, 0.05, 0, 0]],
                                  "excretion_N_removal_efficiency" :['unique', [0, 0, 0, 0, 0]],
                                  "excretion_P_removal_efficiency" : ['unique', [0, 0, 0, 0, 0]],
                                  "solid_filter_efficiency": ['unique', [0.3485, 0, 0, 0, 0]],
                                  "CH4_volume_sludge_and_manure": ['unif', [350, 300, 400, 0, 0]],
                                  "share_fish_sludge_in_substrate": ['unif', [0.33, 0.30, 0.40, 0, 0]],
                                  "ratio_CO2_CH4_biogas_sludge": ['unique', [0.465, 0, 0, 0, 0]],
                                  "N_manure": ['unique', [0.05, 0, 0, 0, 0]],
                                  "P_manure": ['unique', [0.01, 0, 0, 0, 0]],
                                  "fertilizer_substi_digest_N": ['unif', [0.9, 0.8, 1, 0, 0]],
                                  "fertilizer_substi_digest_P": ['unif', [0.6, 0.2, 0.8, 0, 0]],
                                  "fertilizer_substi_digest_K": ['unif', [0.6, 0.2, 0.8, 0, 0]],
                                  "fertilizer_substi_digest_Mg": ['unif', [0.6, 0.2, 0.8, 0, 0]],
                                  
                                  "fertilizer_substi_manure_field_N": ['unif', [0.5, 0.53, 0.55, 0, 0]],
                                  "fertilizer_substi_manure_field_P": ['unif', [0.28, 0.20, 0.385, 0, 0]]}

   
    

Physicdict_distributions = {
    
    'CH4_LHV':['unique', [50, 0, 0, 0, 0]]} # MJ.kg-1
    






# Modification of the dictionnaries specific to the tested scenarios

# To update current to 0 mortality and same FCR
dict_update_sce_A = {"HAIT_loss_red":1,
                     'FRIT_loss_red':1,
                    'GOIT1bis_loss_red':1,
                    'GODK_loss_red':1,
                    'SFDK1_loss_red':1,
                    'SFDK2_loss_red':1}
                    

# To update current to 0 mortality and best FCR

dict_update_sce_B = {"HAIT_loss_red":1,
                     'FRIT_loss_red':1,
                    'GOIT1bis_loss_red':1,
                    'GODK_loss_red':1,
                    'SFDK1_loss_red':1,
                    'SFDK2_loss_red':1,
                    "HAIT_FCR_red_ratio_frac":1,
                    "FRIT_FCR_red_ratio_frac":1,
                    "GOIT1bis_FCR_red_ratio_frac":1,
                    "GODK_FCR_red_ratio_frac":1,
                    "SFDK1_FCR_red_ratio_frac":1,
                    "SFDK2_FCR_red_ratio_frac":1}


# To update current to incumbent mortality and best FCR

dict_update_sce_C = {"HAIT_FCR_red_ratio_frac":1,
                    "FRIT_FCR_red_ratio_frac":1,
                    "GOIT1bis_FCR_red_ratio_frac":1,
                    "GODK_FCR_red_ratio_frac":1,
                    "SFDK1_FCR_red_ratio_frac":1,
                    "SFDK2_FCR_red_ratio_frac":1}



# OAT Loss level increase

dict_update_sce_stage_hatch = {"HAIT_loss_lev":0.15}

dict_update_sce_stage_fry = {"FRIT_loss_lev":0.15}

dict_update_sce_stage_goit = {"GOIT1bis_loss_lev":0.15}


dict_update_sce_stage_godk = {'GODK_loss_lev':0.15}


dict_update_sce_stage_sfdk1 = {'SFDK1_loss_lev':0.15}

dict_update_sce_stage_sfdk2 = {'SFDK2_loss_lev':0.15}


# Put all these scenarios in a list
list_scenarios = [dict_update_sce_A,
                  dict_update_sce_B,
                  dict_update_sce_C,
                  dict_update_sce_stage_hatch,
                  dict_update_sce_stage_fry,
                  dict_update_sce_stage_goit,
                  dict_update_sce_stage_godk,
                  dict_update_sce_stage_sfdk1,
                  dict_update_sce_stage_sfdk2]


scenario_names=  ["dict_update_sce_A",
                  "dict_update_sce_B",
                  "dict_update_sce_C",
                  "dict_update_sce_stage_hatch",
                  "dict_update_sce_stage_fry",
                  "dict_update_sce_stage_goit",
                  "dict_update_sce_stage_godk",
                  "dict_update_sce_stage_sfdk1",
                  "dict_update_sce_stage_sfdk2"]





"""3) SIMULATE"""

#Here we simulate the stochastic LCAs in the different scenarios, 
#using the same Montecarlo results for the background, that we have already simulated and stored (paired)


x = datetime.datetime.now()

month=str(x.month)
day=str(x.day)
microsec=str(x.strftime("%f"))
             
# current primary data   
 
res_INC,r_INC= mainfunc.simulations_fish_micro(Physicdict_distributions,
                           Fish_farm_and_compound_effect_dict_distributions_INC,
                           size, 
                           Physicdict,
                           Fish_farm_and_compound_effect_dict_INC,
                           fishfeed_table_withNP,
                           elemental_contents,
                           demand_vector,
                           Techno_Matrix_Fish,
                           list_FU_combined_names_mc,
                           list_array_total_mc_sorted,  ## Stored MC results for the background
                           activities_fish_background,
                           filter_names_ok_activities_fish,
                           list_meth,
                           Dict_incumbent_losses_growth_stages_loss_level,
                           Dict_incumbent_losses_growth_stages_loss_red,
                           Dict_incumbent_outputs_growth_stages_loss_red,
                           Dict_incumbent_outputs_growth_stages_loss_level,
                           biochem_profile_feed_incumbent,
                           N_P_profile_feed_incumbent,
                           ingredient_profile_incumbent,
                           index_growth_stages,
                           index_feed,
                           index_dead_fish,
                           digestibility_list,
                           index_Nemissions,
                           index_Pemissions,
                           index_growth_stages_no_filter_laguna,
                           index_sludge,
                           electricity_low_voltage_input_per_m3_wastewater,
                           electricity_high_voltage_input_per_m3_wastewater,
                           index_roe,
                           Dict_FCR_bio,
                           list_cfs,
                           ECO_FCR_0,
                           ratio_loss_biological_INC_0,
                           index_growing_DK,
                           index_300g,
                           drug_inputs_names,
                           index_FU,
                           dict_correspondance_techno_growth_stagenames,
                           index_growth_stages_to_modif,
                           bio_FCR_0,
                           index_biogas_updgrade,
                           index_N_substitution,
                           index_P_substitution,
                           index_heat_substitution)    


columns_to_keep =  [meth + "AH" for meth in list_meth_code]
indicators_names = ["Loss ratio",
                      "Economic FCR",
                      "Biological FCR"]


res_sub = res_INC[columns_to_keep+indicators_names]

new_columns_names = [col + "uu" for col in list(res_sub.columns)]

rees = pd.concat([res_INC,res_sub],axis=1)


# Now simulate all the scenarios with modifications of the current status

list_res = []
for scenario_index in range(len(list_scenarios)):
    
    scenario = list_scenarios[scenario_index]
    
    scenario_name =  scenario_names[scenario_index]
    
    Fish_farm_dict_sce = Fish_farm_and_compound_effect_dict_INC.copy()
    
    Fish_farm_dict_sce.update(scenario)
    
    
    np.random.seed(2)
    random.seed(2)
    
    time3=time()
    res,r= mainfunc.simulations_fish_micro(Physicdict_distributions,
                               Fish_farm_and_compound_effect_dict_distributions_INC,
                               size, 
                               Physicdict,
                               Fish_farm_dict_sce,
                               fishfeed_table_withNP,
                               elemental_contents,
                               demand_vector,
                               Techno_Matrix_Fish,
                               list_FU_combined_names_mc,
                               list_array_total_mc_sorted,  ## Stored MC results for the background
                               activities_fish_background,
                               filter_names_ok_activities_fish,
                               list_meth,
                               Dict_incumbent_losses_growth_stages_loss_level,
                               Dict_incumbent_losses_growth_stages_loss_red,
                               Dict_incumbent_outputs_growth_stages_loss_red,
                               Dict_incumbent_outputs_growth_stages_loss_level,
                               biochem_profile_feed_incumbent,
                               N_P_profile_feed_incumbent,
                               ingredient_profile_incumbent,
                               index_growth_stages,
                               index_feed,
                               index_dead_fish,
                               digestibility_list,
                               index_Nemissions,
                               index_Pemissions,
                               index_growth_stages_no_filter_laguna,
                               index_sludge,
                               electricity_low_voltage_input_per_m3_wastewater,
                               electricity_high_voltage_input_per_m3_wastewater,
                               index_roe,
                               Dict_FCR_bio,
                               list_cfs,
                               ECO_FCR_0,
                               ratio_loss_biological_INC_0,
                               index_growing_DK,
                               index_300g,
                               drug_inputs_names,
                               index_FU,
                               dict_correspondance_techno_growth_stagenames,
                               index_growth_stages_to_modif,
                               bio_FCR_0,
                               index_biogas_updgrade,
                               index_N_substitution,
                               index_P_substitution,
                               index_heat_substitution)
    
    
    
    timetot=time()-time3

    print("tiiiime",timetot)
    x = datetime.datetime.now()

    month=str(x.month)
    day=str(x.day)
    microsec=str(x.strftime("%f"))
                 
    
    
    name_file='result_fish'+"_"+month+"_"+day+"_"+microsec + "_" +scenario_name
    name_file_contribution='result_contribution_fisho'+"_"+month+"_"+day+"_"+microsec + "_" +scenario_name
    
    
    res.to_csv("../Outputs_Fish/"+str(name_file)+'.csv', sep=';', encoding='utf-8')
    
    res_sub = res[columns_to_keep+indicators_names]
    
    new_columns_names = [col + scenario_name  for col in list(res_sub.columns)]

    res_sub.columns = new_columns_names
    
    list_res.append(res_sub)



    export_pickle_2(r, name_file_contribution, "Outputs_Fish")


#res_INC + list_res
combined_res = pd.concat(list_res,axis = 1)
res_tot = pd.concat([res_INC,combined_res], axis = 1)

# Export total results

name_file_tot ='result_fish'+"_"+month+"_"+day+"_"+microsec + "_Total"

res_tot.to_csv("../Outputs_Fish/"+str(name_file_tot)+'.csv', sep=';', encoding='utf-8')



