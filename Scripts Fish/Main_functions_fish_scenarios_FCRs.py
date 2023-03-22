# -*- coding: utf-8 -*-
"""
Created on Sun Dec 18 12:50:34 2022


@author: Pierre Jouannais, Department of Planning, DCEA, Aalborg University
pijo@plan.aau.dk

"""
'''
Script containing the main functions perfoming the simulations 
of the fish farm LCAs


'''


import datetime
from time import *
import pickle

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
from SALib.analyze import delta

import SALib

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits import mplot3d

import shapefile as shp
import geopandas as gpd
import pysal as ps
from shapely.geometry import Polygon, mapping, Point

import Functions_for_physical_and_biological_calculations_3rd as functions


import Technosphere_matrix_modifications_fish_3rd as techno_modif

import ray


import scipy
import numpy as np
from scipy import linalg

from ast import literal_eval




"""Sampling Functions """


    
def sampling_func_total_montecarlo_fish(
                  Physicdict_distributions,
                  Fish_farm_and_compound_effect_dict_distributions,
                  size):
    '''Function which returns a random sample for the input space of 
    non-constant parameters 
    
    Inputs:
        
        # All parameters distributions dictionnaries : 

                  Physicdict_distributions
                  Fish_farm_and_compound_effect_dict_distributions
                  
        #size : Size of the sample 
    
    
    Outputs :
        
        -sample :  Generated sample. Array with 1 row = 1 combination of uncertain parameters
        -names_param :  List of names of the uncertain parameters
        -names_param_phy :  List of names of the uncertain parameters from Physicdict_distributions
        -names_param_fish_farm_compound_effect : List of names of the uncertain parameters from fish_farm_compound_effect_distributions
    '''
    

    # Creation of the problem

    # Creating a copy of the distribtutions dictionnaires containing only the
    # parameters with distributions(not the constant ones)

    # Physics

    Physicdict_distributions_input = Physicdict_distributions.copy()

    to_be_deleted_phy = []

    for param in Physicdict_distributions_input:
        if Physicdict_distributions_input[param][0] == 'unique' or Physicdict_distributions_input[param][0] == 'binary':
            to_be_deleted_phy.append(param)

    for a in to_be_deleted_phy:
        Physicdict_distributions_input.pop(a)

    # Fish farm

    Fish_farm_and_compound_effect_dict_distributions_input = Fish_farm_and_compound_effect_dict_distributions.copy()

    to_be_deleted_fish = []

    for param in Fish_farm_and_compound_effect_dict_distributions_input:
        
        if Fish_farm_and_compound_effect_dict_distributions_input[param][0] == 'unique' or Fish_farm_and_compound_effect_dict_distributions_input[param][0] == 'binary':
            
            to_be_deleted_fish.append(param)

    for a in to_be_deleted_fish:
        Fish_farm_and_compound_effect_dict_distributions_input.pop(a)
        
       

    # Collecting names, bounds , dists to create Saltelli problem
    names_param = []
    bounds = []
    dists = []
    

    #  Physics
    
    names_param_phy = []

    for param in Physicdict_distributions_input:
        
        distrib = Physicdict_distributions_input[param][0]
        
        dists.append(distrib)

        
        names_param_phy.append(param)

        if distrib == 'unif':  # then bounds = upper bound, lower bound

            bounds.append([Physicdict_distributions_input[param]
                           [1][1], Physicdict_distributions_input[param][1][2]])

        elif distrib == 'norm':  # then bounds = mean, sd

            bounds.append([Physicdict_distributions_input[param]
                           [1][3], Physicdict_distributions_input[param][1][4]])

        elif distrib == 'triang':  # then bounds = width, location of the peak in % of width, Assume lower bound = 0
            
            bounds.append([Physicdict_distributions_input[param][1][1],
                       Physicdict_distributions_input[param][1][3],
                       Physicdict_distributions_input[param][1][2]])
            




    # Fish farm 
    
    names_param_fish_farm_compound_effect = []

    for param in Fish_farm_and_compound_effect_dict_distributions_input:
        
        distrib = Fish_farm_and_compound_effect_dict_distributions_input[param][0]
        
        dists.append(distrib)

        
        names_param_fish_farm_compound_effect.append(param)

        if distrib == 'unif':  # then bounds = upper bound, lower bound

            bounds.append([Fish_farm_and_compound_effect_dict_distributions_input[param]
                           [1][1], Fish_farm_and_compound_effect_dict_distributions_input[param][1][2]])

        elif distrib == 'norm':  # then bounds = mean, sd

            bounds.append([Fish_farm_and_compound_effect_dict_distributions_input[param]
                           [1][3], Fish_farm_and_compound_effect_dict_distributions_input[param][1][4]])
            

        elif distrib == 'triang':  # then bounds = width, location of the peak in % of width, Assume lower bound = 0


            bounds.append([Fish_farm_and_compound_effect_dict_distributions_input[param][1][1],
                       Fish_farm_and_compound_effect_dict_distributions_input[param][1][3],
                       Fish_farm_and_compound_effect_dict_distributions_input[param][1][2]])            
                        
            


    names_param =  names_param_phy + names_param_fish_farm_compound_effect 
    
    
    # Create Empty array for sample
    

    sample_array = np.array([[0]*len(names_param)]*size,dtype="float32")

    def random_draw(dist_param,bounds_param):
        
        if dist_param == "unif":
            value = np.random.uniform(bounds_param[0], bounds_param[1], size=None)
            ##print("value",value)
        if dist_param == "norm":
            value = np.random.normal(bounds_param[0], bounds_param[1], size=None)
            ##print("value",value)
            
        if dist_param == "triang":
            value = np.random.triangular(bounds_param[0], bounds_param[1], bounds_param[2])
            ##print("value",value)
        return(value)    
            

    ##print(dists)
    for row_index in range(sample_array.shape[0]):

        
        sample_row =[random_draw(dist_param,bounds_param) for dist_param,bounds_param in zip(dists,bounds)]
        sample_array[row_index,:]=sample_row






    return sample_array, names_param, names_param_phy, names_param_fish_farm_compound_effect






""" Calculating functions"""




def simulations_fish_micro(Physicdict_distributions,
                           Fish_farm_and_compound_effect_dict_distributions,
                           size, 
                           Physicdict,
                           Fish_farm_and_compound_effect_dict,
                           fishfeed_table_withNP,
                           elemental_contents,
                           demand_vector,
                           Techno_Matrix_Fish,
                           list_FU_combined_names_mc,
                           list_array_total_mc_sorted,
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
                           index_heat_substitution):
    
    """Function which calls the other functions to simulate all the MonteCarlo 
    iterations and collect all results"""      
    
  
    

    # Generate sample of parameters
    
    timeA=time()

    
    output_sample_fish_micro = sampling_func_total_montecarlo_fish(Physicdict_distributions,
                  Fish_farm_and_compound_effect_dict_distributions,
                  size)
    
    # return sample_array, names_param, names_param_phy, names_param_fish_farm_compound_effect

    # Disaggregate results
    sample_total_input=output_sample_fish_micro[0]
    
    
    names_param_set=output_sample_fish_micro[1]    
   
    
    names_param_phy=output_sample_fish_micro[2]
    
    names_param_fish_farm_compound_effect=output_sample_fish_micro[3]
    
        
    

    
    time1=time()
    
    
    # Start Ray.
    ray.shutdown()
    
    ray.init()

    # Common variable across ray simulations
    constant_inputs = ray.put([Physicdict,
                             Fish_farm_and_compound_effect_dict,
                             fishfeed_table_withNP,
                             elemental_contents,
                             names_param_phy,
                             names_param_fish_farm_compound_effect,
                             names_param_set,
                             demand_vector,
                             Techno_Matrix_Fish,
                             biochem_profile_feed_incumbent,
                             index_growth_stages,
                             index_feed,
                             index_dead_fish,
                             index_Nemissions,
                             index_Pemissions,
                             index_growth_stages_no_filter_laguna,
                             index_sludge,
                             electricity_low_voltage_input_per_m3_wastewater,
                             electricity_high_voltage_input_per_m3_wastewater,
                             ingredient_profile_incumbent,
                             N_P_profile_feed_incumbent,
                             digestibility_list,
                             index_roe,
                             Dict_FCR_bio,
                             ECO_FCR_0,
                             ratio_loss_biological_INC_0,
                             index_growing_DK,
                             index_300g,
                             index_FU,
                             list_meth,
                             dict_correspondance_techno_growth_stagenames,
                             index_growth_stages_to_modif,
                             bio_FCR_0,
                             index_biogas_updgrade,
                             index_N_substitution,
                             index_P_substitution,
                             index_heat_substitution]) 

    # Calculate LCIs in parallel
    
    arrayresult_raw =ray.get([calculateLCI_1param_parallel.remote(constant_inputs,
                                                                  param_set) for param_set in sample_total_input])        


    time_B=time()
    time_LCI =time_B-timeA
    print("time_LCI",time_LCI)


    """ Calculate LCIAs in parallel by combining LCis with mc results for the background"""
        
    # Now calculate LCIAs in parallel
    
    constant_inputs_LCIA = ray.put([filter_names_ok_activities_fish,
                                    list_FU_combined_names_mc,
                                    activities_fish_background,
                                    list_meth,
                                    Dict_incumbent_losses_growth_stages_loss_level,
                                    Dict_incumbent_losses_growth_stages_loss_red,
                                    Dict_incumbent_outputs_growth_stages_loss_red,
                                    Dict_incumbent_outputs_growth_stages_loss_level,
                                    names_param_set,
                                    Fish_farm_and_compound_effect_dict,
                                    list_cfs])
    



    arrayresult_LCIA = ray.get([LCIA_parallel.remote(constant_inputs_LCIA,
                                                                  row_LCI,row_mc) for row_LCI,row_mc in zip(arrayresult_raw,list_array_total_mc_sorted)])
    
    
    ray.shutdown()
    #print("DONE LCIIIIAAA")
    
    print("Done with LCIAS")
    time_LCIA=time()-time_B
    print("time_LCIA",time_LCIA)
    

    # Rebuild a proper dataframe
    
    # Separate the contributions and the absolute results
    
    table_contrib=[pd.DataFrame(np.zeros([len(arrayresult_LCIA),len(activities_fish_background)]),columns=activities_fish_background)for meth in range(len(list_meth))]

    
    list_result_LCIA_without_contrib =[]


    print("Now sorting results")
    
    for row_index in range(len(arrayresult_LCIA)):
        
        row_LCIA =arrayresult_LCIA[row_index]  # Collect the row containing the LCIA
        
        ##print("len(row_LCIA)",len(row_LCIA)) # Should be 2 : [[absolute],[contrib]]
       # #print("len(row_LCIA[0])",len(row_LCIA[0])) 

        row_performance_indicators = arrayresult_raw[row_index][-1]  # Collect the row containing the performance indicators
        
        list_result_LCIA_without_contrib.append(row_LCIA[0]+row_performance_indicators)

        for index_meth in range(len(list_meth)):
            
            table_contrib[index_meth] = row_LCIA[-1][index_meth]            
 


    # Collect the rest before returning everything

        
    list_meth_short = [meth[-1] for meth in list_meth]
    
    names_total_AH = [meth+"AH" for meth in list_meth_short]
    

    name_indicators =["Loss ratio",
                      "Economic FCR",
                      "Biological FCR"]
 


    names_col_dataframe = names_param_set +  names_total_AH  + activities_fish_background  + name_indicators


    results_table_df = pd.DataFrame(np.array(list_result_LCIA_without_contrib), columns=names_col_dataframe)
    


    return results_table_df, table_contrib
    
    
    



@ray.remote
def calculateLCI_1param_parallel(constant_inputs,param_set) :   
        """Function which calculates one LCI and returns the supply vector for
        a set of input parameters.
        The function is built to be called in parallel with ray."""      
    

        print("1LCI")
    
        (Physicdict,
                             Fish_farm_and_compound_effect_dict,
                             fishfeed_table_withNP,
                             elemental_contents,
                             names_param_phy,
                             names_param_fish_farm_compound_effect,
                             names_param,
                             demand_vector,
                             Techno_Matrix_Fish,
                             biochem_profile_feed_incumbent,
                             index_growth_stages,
                             index_feed,
                             index_dead_fish,
                             index_Nemissions,
                             index_Pemissions,
                             index_growth_stages_no_filter_laguna,
                             index_sludge,
                             electricity_low_voltage_input_per_m3_wastewater,
                             electricity_high_voltage_input_per_m3_wastewater,
                             ingredient_profile_incumbent,
                             N_P_profile_feed_incumbent,
                             digestibility_list,
                             index_roe,
                             Dict_FCR_bio,
                             ECO_FCR_0,
                             ratio_loss_biological_INC_0,
                             index_growing_DK,
                             index_300g,
                             index_FU,
                             list_meth,
                             dict_correspondance_techno_growth_stagenames,
                             index_growth_stages_to_modif,
                             bio_FCR_0,
                             index_biogas_updgrade,
                             index_N_substitution,
                             index_P_substitution,
                             index_heat_substitution)=constant_inputs    
    

        # Update the dictionnaries whith the values of the sample

       

        #  Physicdict
        
        new_start = 0 
        for param in Physicdict: 
            # We browse the parameters to look for the variable ones
            # that need to be updated

            for index in range(new_start, new_start+len(names_param_phy)):


                # Looking for the corresponding parameter in the set
                if names_param[index] == param:

                    Physicdict[param] = param_set[index]

        new_start = new_start+len(names_param_phy)


        # Fish_farm_and_compound_effect_dict
        
        for param in Fish_farm_and_compound_effect_dict: 
            # We browse the parameters to look for the variable ones
            # that need to be updated

            for index in range(new_start, new_start+len(names_param_fish_farm_compound_effect)):

                # Looking for the corresponding parameter in the saltelli set
                if names_param[index] == param:

                    Fish_farm_and_compound_effect_dict[param] = param_set[index]

        new_start = new_start+len(names_param_fish_farm_compound_effect)
        
         
        
        """Calculate the LCI for the fish farms"""
        
        
        
  

        # Collect necessary values in the dict for intermediate calculations
        
        
        CH4_LHV = Physicdict["CH4_LHV"]
        

        HAIT_FCR_red_ratio_frac = Fish_farm_and_compound_effect_dict["HAIT_FCR_red_ratio_frac"]
        FRIT_FCR_red_ratio_frac = Fish_farm_and_compound_effect_dict["FRIT_FCR_red_ratio_frac"]
        GOIT1_FCR_red_ratio_frac = Fish_farm_and_compound_effect_dict["GOIT1_FCR_red_ratio_frac"]
        GOIT2_FCR_red_ratio_frac = Fish_farm_and_compound_effect_dict["GOIT2_FCR_red_ratio_frac"]
        GOIT1bis_FCR_red_ratio_frac = Fish_farm_and_compound_effect_dict["GOIT1bis_FCR_red_ratio_frac"]
        GODK_FCR_red_ratio_frac = Fish_farm_and_compound_effect_dict["GODK_FCR_red_ratio_frac"]
        SFDK1_FCR_red_ratio_frac = Fish_farm_and_compound_effect_dict["SFDK1_FCR_red_ratio_frac"]
        SFDK2_FCR_red_ratio_frac = Fish_farm_and_compound_effect_dict["SFDK2_FCR_red_ratio_frac"]


        ratio_CO2_CH4_biogas_fish = Fish_farm_and_compound_effect_dict["ratio_CO2_CH4_biogas_fish"] 
        CH4_volume_biogas_fish = Fish_farm_and_compound_effect_dict["CH4_volume_biogas_fish"]
        P_in_fish = Fish_farm_and_compound_effect_dict["P_in_fish"]

        Mg_in_fish =Fish_farm_and_compound_effect_dict["Mg_in_fish"]
        water_in_fish = Fish_farm_and_compound_effect_dict["water_in_fish"]
        gvs_gts_in_fish = Fish_farm_and_compound_effect_dict["gvs_gts_in_fish"]
        
        fraction_non_ingested = Fish_farm_and_compound_effect_dict["fraction_non_ingested"]


        # Fish waste management
        [MJ_substituted_per_kilo_deadfish,
         amount_of_upgrading_act, 
         biogenic_CO2_emitted] = functions.natural_gas_subst(
                    ratio_CO2_CH4_biogas_fish,
                    CH4_volume_biogas_fish,
                    CH4_LHV,
                    water_in_fish,
                    gvs_gts_in_fish)
        
                                                             
        # Calculate the minimum FCR with constant digestibility
        
        min_FCR = P_in_fish/(N_P_profile_feed_incumbent[1]*(digestibility_list[-1]-fraction_non_ingested))

        min_HAIT_FCR_red_ratio = min_FCR/Dict_FCR_bio['HAIT_FCR'] 
        min_FRIT_FCR_red_ratio = min_FCR/Dict_FCR_bio['FRIT_FCR'] 
        min_GOIT1bis_FCR_red_ratio = min_FCR/Dict_FCR_bio['GOIT1bis_FCR'] 
        min_GODK_FCR_red_ratio = min_FCR/Dict_FCR_bio['GODK_FCR'] 
        min_SFDK1_FCR_red_ratio = min_FCR/Dict_FCR_bio['SFDK1_FCR'] 
        min_SFDK2_FCR_red_ratio = min_FCR/Dict_FCR_bio['SFDK2_FCR'] 
        
        # Calculate the FCR reduction ratios according to the input parameters and the minimum FCRs
        HAIT_FCR_red_ratio =  (min_HAIT_FCR_red_ratio-1)* HAIT_FCR_red_ratio_frac + 1
        FRIT_FCR_red_ratio =  (min_FRIT_FCR_red_ratio-1)* FRIT_FCR_red_ratio_frac + 1
        GOIT1bis_FCR_red_ratio =  (min_GOIT1bis_FCR_red_ratio-1)* GOIT1bis_FCR_red_ratio_frac + 1
        GODK_FCR_red_ratio =  (min_GODK_FCR_red_ratio-1)* GODK_FCR_red_ratio_frac + 1
        SFDK1_FCR_red_ratio =  (min_SFDK1_FCR_red_ratio-1)* SFDK1_FCR_red_ratio_frac + 1
        SFDK2_FCR_red_ratio =  (min_SFDK2_FCR_red_ratio-1)* SFDK2_FCR_red_ratio_frac + 1
        

        # Store them in a new dictionnary
        dict_param_fish_farms_ready = Fish_farm_and_compound_effect_dict.copy() 
        
        dict_param_fish_farms_ready["HAIT_FCR_red_ratio"] = HAIT_FCR_red_ratio
        dict_param_fish_farms_ready["FRIT_FCR_red_ratio"]  = FRIT_FCR_red_ratio
        dict_param_fish_farms_ready["GOIT1bis_FCR_red_ratio"]  = GOIT1bis_FCR_red_ratio
        dict_param_fish_farms_ready["GODK_FCR_red_ratio"] = GODK_FCR_red_ratio
        dict_param_fish_farms_ready["SFDK1_FCR_red_ratio"]  = SFDK1_FCR_red_ratio
        dict_param_fish_farms_ready["SFDK2_FCR_red_ratio"] = SFDK2_FCR_red_ratio
        


        
        




        # Update technosphere with the new parameters
        
        Tech_num_AH_modif_with_excretion = techno_modif.calculate_fish_technosphere_matrix_scenario(Techno_Matrix_Fish,
                                            dict_param_fish_farms_ready,
                                            dict_correspondance_techno_growth_stagenames,
                                            index_growth_stages_to_modif,
                                            index_growth_stages,
                                            index_dead_fish,
                                            index_roe,
                                            MJ_substituted_per_kilo_deadfish,
                                            amount_of_upgrading_act,
                                            biogenic_CO2_emitted,
                                            N_P_profile_feed_incumbent,
                                            biochem_profile_feed_incumbent,
                                            fraction_non_ingested,
                                            digestibility_list,
                                            CH4_LHV,
                                            index_growth_stages_no_filter_laguna,
                                            index_feed,
                                            index_Nemissions,
                                            index_Pemissions,
                                            index_sludge,
                                            index_biogas_updgrade,
                                            index_N_substitution,
                                            index_P_substitution,
                                            index_heat_substitution,
                                            index_300g,
                                            index_growing_DK,
                                            index_FU)


        # Supply vector in the foreground           
        supply_vector_tot_AH = linalg.solve(Tech_num_AH_modif_with_excretion, demand_vector)
        
        
        # Activities in supply vector are ranked as in list_ acitviites -Names
        
        
        
        
                
        # Calculate Performance indicators
        
        [_,
        ECO_FCR_AH_1,
        bio_FCR_AH_1,
        Loss_fraction_AH_1,
        ratio_loss_biological_AH_1]=calculate_economic_indicators(Tech_num_AH_modif_with_excretion,
                                  index_300g,
                                  index_growing_DK,
                                  index_roe,
                                  index_FU,
                                  index_feed,
                                  index_dead_fish)  
                                                                  


        Loss_fraction_increase_AH_indic =   ratio_loss_biological_AH_1/ ratio_loss_biological_INC_0                                                             
        


        Economic_FCR_red_0_indic = (ECO_FCR_0-ECO_FCR_AH_1)/ECO_FCR_0
                
        Biological_FCR_red_0_indic = (bio_FCR_0-bio_FCR_AH_1)/bio_FCR_0
        

        
  

        performance_indicators= [ratio_loss_biological_AH_1,
                                ECO_FCR_AH_1,
                                bio_FCR_AH_1]
        
        return supply_vector_tot_AH, param_set, performance_indicators





    
@ray.remote
def LCIA_parallel(constant_inputs_LCIA, 
                  row_LCI, row_mc):
    
    [filter_names_ok_activities_fish,
     list_FU_combined_names_mc,
     activities_fish_background,
     list_meth,
    Dict_incumbent_losses_growth_stages_loss_level,
    Dict_incumbent_losses_growth_stages_loss_red,
    Dict_incumbent_outputs_growth_stages_loss_red,
    Dict_incumbent_outputs_growth_stages_loss_level,
    names_param_set,
    Fish_farm_and_compound_effect_dict,
    list_cfs]= constant_inputs_LCIA

    """Function which calculates one LCIA for
    a set of input parameters.
    The function is built to be called in parallel with ray."""        

    print("1LCIA")

    # Disaggregate inputs
    supply_vector_AH = row_LCI[0]
        
    param_set = row_LCI[1]
    

    # Keep only the LCI for activities directly connected to a background activity
    # Exclude the different growth stages
    LCI_fish_AH = list(compress(supply_vector_AH, filter_names_ok_activities_fish)) 
    
    # These activities are in the same order as their names in "activities_fish_background"
    
    
    # Put the mc results for 1 unit of each activity of the foreground
    # in a dictionnary with the name of the activity as the key
    Dict_mc_FU = {}
    
    # Prepare dictionnary such as { "Actname": [impact for this mc with meth1,impact for this mc with meth2],
    #                              "Actname2:"[impact for this mc with meth1,impact for this mc with meth2]]}
        

    for index_name_FU in range(len(row_mc[0])):  

        name_FU = list_FU_combined_names_mc[index_name_FU]
                
        Dict_mc_FU[name_FU] = [row_mc[meth_index][index_name_FU] for meth_index in range(len(row_mc))]
        
        

    # Put the paramters values in a dictionnary with paramter name as key
    
    Dict_param_set_values = {}
    
    for index_name_param in range(len(param_set)):  
        
        Dict_param_set_values[names_param_set[index_name_param]] = param_set[index_name_param]
    
    
    # Calculate LCIA
    LCIAs_AH=[]  
    
    for index_meth in range(len(list_meth)): # For each method
        
        for_1_meth_AH=[]
        for index_FU in range(len(LCI_fish_AH)):# For each foreground activity
            
            # Amount for this activity * Impact per unit of this activity
            
            for_1_meth_AH.append(LCI_fish_AH[index_FU] * Dict_mc_FU[activities_fish_background[index_FU]][index_meth])
           
            # For 1 meth [Total Impact due to first act, Total impact due to second act, ...]
            
        LCIAs_AH.append(for_1_meth_AH)
        
    # Total impact = SUM of impacts due to each activity in the foreground
    Total_impacts_per_kg_fish_AH = [np.sum(list_LCIAs_meth) for list_LCIAs_meth in LCIAs_AH]  # [impact meth1, impact meth 2 ---]
    
    
    # Calculate contributions
   
    list_contributions_AH=[]
    

    for index_meth in range(len(list_meth)):
        
        cont_for_1_meth_AH = [impact_input/Total_impacts_per_kg_fish_AH[index_meth] for impact_input in LCIAs_AH[index_meth]]
        
        
        list_contributions_AH.append(cont_for_1_meth_AH)
        
    
 
    row = list(param_set) + list(Total_impacts_per_kg_fish_AH) + list(LCI_fish_AH)  
  


    return [row, list_contributions_AH]














    


    

def calculate_economic_indicators(Tech_Matrix,
                                  index_300g,
                                  index_growing_DK,
                                  index_roe,
                                  index_FU,
                                  index_feed,
                                  index_dead_fish):
    
    """Function which calculates perfomance indicators associated with a tehcnosphere matrix
    Tech_Matrix = Numerical matrix"""
    
    Tech_Matrix_modif=Tech_Matrix.copy()
    
    # 1) Add the output of 300 g trout as an economic output
    Tech_Matrix_modif[index_FU,index_growing_DK] = Tech_Matrix_modif[index_FU,index_growing_DK] +Tech_Matrix_modif[index_300g,index_growing_DK]
    
    Tech_Matrix_modif[index_300g,index_growing_DK] = 0
    
    # 2) Add the roe output as a an economic output
    
    Tech_Matrix_modif[index_FU,index_FU] = Tech_Matrix_modif[index_FU,index_FU] + Tech_Matrix_modif[index_roe,index_FU]

    
    Tech_Matrix_modif[index_roe,index_FU] = 0
    
    
    demand_vector_unity = [0 for i in range(len(Tech_Matrix_modif))]
    demand_vector_unity[index_FU] = 1

    supply_vector_tot = linalg.solve(Tech_Matrix_modif, demand_vector_unity)
    
    eco_FCR =  supply_vector_tot[index_feed]    # Feed per kilo economic output = Sea raised trout(FU) + 300g trout + Roe
    Loss_ratio = supply_vector_tot[index_dead_fish]
    ratio_loss_biological = Loss_ratio/(Loss_ratio+1)
    
    # FCR biological 
    bio_FCR =  eco_FCR/(1+Loss_ratio) # Feed per kilo biological output = Sea raised trout(FU) + 300g trout + Roe +Dead
    

    return Tech_Matrix_modif,eco_FCR,bio_FCR, Loss_ratio, ratio_loss_biological

    

    