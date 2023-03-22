# -*- coding: utf-8 -*-
"""
Created on Mon May 31 15:46:57 2021

@author: Pierre Jouannais, Department of Planning, DCEA, Aalborg University
pijo@plan.aau.dk




"""


'''
Script containing the functions for anaerobic digestion and valorization of digestates
for dead fish and sludge

'''

from scipy.optimize import minimize
import pandas as pd
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

currentfolder=os.path.dirname(os.path.realpath(__file__))
os.chdir(currentfolder)


##

elemental_contents = pd.read_csv(
                    "../Data/elemental_contents.csv", sep=";",
                    header=0, encoding='unicode_escape',
                    engine='python')
elemental_contents = elemental_contents.iloc[:, 1:]
##









def natural_gas_subst(ratio_CO2_CH4_biogas_fish,
                    CH4_volume_biogas_fish,
                    CH4_LHV,
                    water_in_fish,
                    gvs_gts_in_fish):  # m-3

    """Returns the amount of substituted heat MJ and of biogas upgrading activity per
    kg of dead fish wet weight"""
    
    gvs_kg_fish = (1-water_in_fish) * gvs_gts_in_fish *1000
    
    V_CH4_biogas = gvs_kg_fish * CH4_volume_biogas_fish #ml
    
    V_CO2_biogas = V_CH4_biogas * ratio_CO2_CH4_biogas_fish #ml
    
    Volume_biogas =  V_CH4_biogas + V_CO2_biogas #ml
    

    amount_of_upgrading_act = (Volume_biogas*10**-6) #This is the amount of upgrading activity needed to upgrade the methane and make it substitutable to natural gas

    #print("amount_of_upgrading_act",amount_of_upgrading_act)

    

    MJ_substituted_per_kilo_deadfish = gvs_kg_fish * CH4_volume_biogas_fish * (0.66/1000) *CH4_LHV/1000# ml * g/ml * Mj.kg-1/1000

 
    g_CO2_biogas = 0.00187 * V_CO2_biogas # g.ml-1 * emitted as biogenic carbon




    biogenic_CO2_emitted = (g_CO2_biogas + V_CH4_biogas*(0.66/1000)/16 * 44)/1000 # kg

    
    return MJ_substituted_per_kilo_deadfish,  amount_of_upgrading_act,  biogenic_CO2_emitted
    





def fish_sludge_management(CH4_volume_sludge_and_manure,
                           share_fish_sludge_in_substrate,
                           CH4_LHV,
                           ratio_CO2_CH4_biogas_sludge,
                           N_in_sludge_dw,
                           P_in_sludge_dw,
                           N_manure,
                           P_manure,
                           fertilizer_substi_digest_N,
                           fertilizer_substi_digest_P,
                           fertilizer_substi_manure_field_N,
                           fertilizer_substi_manure_field_P):
    """Returns the amount of substituted heat MJ,of biogas upgrading activity,
    of N, P ferlitizers per kg dry weight fish sludge"""
    
    [MJ_substituted_per_kilo_dw_sludge, 
     amount_of_upgrading_act_sludge,  
     biogenic_CO2_emitted] = natural_gas_subst(ratio_CO2_CH4_biogas_sludge,
                    CH4_volume_sludge_and_manure,
                    CH4_LHV,
                    0,  # No water, already dry weight
                    1)   # Volatile_solid / total solid, already dry weight
                       
             
                              
                              
    needed_dm_manure = (1/share_fish_sludge_in_substrate-1)*61/130  # Modeled according to dry weight balances in Brod et al. 2017
                                    
   
    N_substitution_manure_co_sub = needed_dm_manure * N_manure * (fertilizer_substi_digest_N - fertilizer_substi_manure_field_N)
    
    P_substitution_manure_co_sub = needed_dm_manure * P_manure * (fertilizer_substi_digest_P - fertilizer_substi_manure_field_P)
    
    # Less of this substrate on the fields
   
    fertilizer_subst_N = N_in_sludge_dw * fertilizer_substi_digest_N + N_substitution_manure_co_sub
                           
    fertilizer_subst_P = P_in_sludge_dw * fertilizer_substi_digest_P + P_substitution_manure_co_sub
    
    return fertilizer_subst_N,fertilizer_subst_P,MJ_substituted_per_kilo_dw_sludge,amount_of_upgrading_act_sludge
    

    


def biochemprofile(fishfeed_table):
    """Returns the current biochemical/nutritional profile
        """
    #############
    # 1)Calculation current nutritional profile
    #############

    # Multiplying each ingredien's mass by its associated nutrient content
    incumbent_lipid = sum(
        fishfeed_table['kg.kg feed-1']*fishfeed_table['Lipid'])

    incumbent_prot = sum(
        fishfeed_table['kg.kg feed-1']*fishfeed_table['Protein'])

    incumbent_carb = sum(fishfeed_table['kg.kg feed-1']*fishfeed_table['Carb'])

    incumbent_ash = sum(fishfeed_table['kg.kg feed-1']*fishfeed_table['Ash'])

    incumbent_water = sum(
        fishfeed_table['kg.kg feed-1']*fishfeed_table['Water'])


    # Fish feed standard macronutrient profile Lip,prot,carb,ashes,water
    incumbentvect = [incumbent_lipid, incumbent_prot,
                     incumbent_carb, incumbent_ash, incumbent_water]

        

    return incumbentvect


