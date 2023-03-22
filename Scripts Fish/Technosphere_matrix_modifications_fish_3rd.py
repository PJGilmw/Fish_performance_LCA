# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 14:24:06 2022

@author: Pierre Jouannais, Department of Planning, DCEA, Aalborg University
pijo@plan.aau.dk

"""
'''
Script containing the functions which modify the growth stages, 
update the technosphere matrix and calculate the excretion.


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
import random

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits import mplot3d


import ray
import scipy
import numpy as np
from scipy import linalg

from ast import literal_eval
from itertools import compress

import Functions_for_physical_and_biological_calculations_3rd as functions





             
def calculate_A(FCR_red_ratio,Fish_dead_0,Total_output_fish,Fish_input_0):

    num = Fish_dead_0 + Total_output_fish-Fish_input_0 * (1-FCR_red_ratio) 
    A = num/FCR_red_ratio
    
    return A

def calculate_F_dead_1(A,Fish_dead_0,Fish_input_0,loss_lev,loss_red,Total_output_fish):

    Total_bio = Fish_dead_0 + Total_output_fish

    F_dead1 = A * (1-loss_red) * (Fish_dead_0 + loss_lev * Total_output_fish)/Total_bio
    return F_dead1 

def calculate_F_out_x_1(A,Fish_dead_0,Fish_input_0,loss_lev,loss_red,F_out_0,Total_output_fish):
    
    Total_bio = Fish_dead_0 + Total_output_fish

    F_outx1 = (A/Total_bio) *(F_out_0 *( 1-loss_lev + loss_red * loss_lev + loss_red * Fish_dead_0/Total_output_fish))
    
    return F_outx1







def calculate_fish_technosphere_matrix_scenario(Techno_Matrix_Fish,
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
                                            index_FU):
    
    '''
    Function which updates the tehcnosphere matrix according to the sampled parameters. 
    Modification of the fish flows, sludge and N,P excretion

    
      #Outputs : 
          
          Updated technosphere matrix
         
    '''
    
    # Collect Parameters
    
    P_in_fish=dict_param_fish_farms_ready["P_in_fish"]
    N_in_fish=dict_param_fish_farms_ready["N_in_fish"] 
    K_in_fish=dict_param_fish_farms_ready["K_in_fish"]
    Mg_in_fish=dict_param_fish_farms_ready["Mg_in_fish"]
    water_in_fish=dict_param_fish_farms_ready["water_in_fish"]
    gvs_gts_in_fish=dict_param_fish_farms_ready["gvs_gts_in_fish"] 
        

    excretion_N_removal_efficiency =dict_param_fish_farms_ready["excretion_N_removal_efficiency"]
    excretion_P_removal_efficiency =dict_param_fish_farms_ready["excretion_P_removal_efficiency"]
    solid_filter_efficiency =dict_param_fish_farms_ready["solid_filter_efficiency"]
    
   
    fraction_non_ingested = dict_param_fish_farms_ready["fraction_non_ingested"]
    
    CH4_volume_sludge_and_manure = dict_param_fish_farms_ready["CH4_volume_sludge_and_manure"]

    share_fish_sludge_in_substrate = dict_param_fish_farms_ready["share_fish_sludge_in_substrate"]
     
    ratio_CO2_CH4_biogas_sludge = dict_param_fish_farms_ready["ratio_CO2_CH4_biogas_sludge"]
    N_manure = dict_param_fish_farms_ready["N_manure"]
    P_manure= dict_param_fish_farms_ready["P_manure"]
    
    

    
    fertilizer_substi_digest_N = dict_param_fish_farms_ready["fertilizer_substi_digest_N"]
    fertilizer_substi_digest_P = dict_param_fish_farms_ready["fertilizer_substi_digest_P"]
    fertilizer_substi_digest_K = dict_param_fish_farms_ready["fertilizer_substi_digest_K"]
    fertilizer_substi_digest_Mg = dict_param_fish_farms_ready["fertilizer_substi_digest_Mg"]

    fertilizer_substi_manure_field_N = dict_param_fish_farms_ready["fertilizer_substi_manure_field_N"]
    fertilizer_substi_manure_field_P = dict_param_fish_farms_ready["fertilizer_substi_manure_field_P"]
    

    

    
    
    # Initialize the numerical Tehcnosphere
    
    Techno_fish_num_AH = np.zeros((Techno_Matrix_Fish.shape[0],Techno_Matrix_Fish.shape[0]))
    
    # Evaluate the values of paramters in the technosphere matrix
    
    for i in range(Techno_fish_num_AH.shape[0]):
    
        for j in range(Techno_fish_num_AH.shape[1]):
            
            if type(Techno_Matrix_Fish[i,j])==str: # Then we evaluate the string
            
                Techno_fish_num_AH[i,j] = eval(Techno_Matrix_Fish[i,j])
                
            else:        # Then we just take the float
                Techno_fish_num_AH[i,j] = Techno_Matrix_Fish[i,j]

    
    

    Tech_num_AH_modif = Techno_fish_num_AH.copy()
    
    

    # First modif 
    for column_index in index_growth_stages_to_modif:
        
        # Associate growth stages numbers and parameters
        
        
        # Code
        code_growth_stage = dict_correspondance_techno_growth_stagenames[column_index]

        # Collect the names of the parameters corresponding to the growth stage
        name_corresponding_FCR_red_ratio =  code_growth_stage+"_FCR_red_ratio"
        
        name_corresponding_loss_level =  code_growth_stage+"_loss_lev"
        
        name_corresponding_loss_red =  code_growth_stage+"_loss_red"
        
        
        # Collect Values parameters with their names
        FCR_red_ratio = dict_param_fish_farms_ready[name_corresponding_FCR_red_ratio]

        loss_lev = dict_param_fish_farms_ready[name_corresponding_loss_level]
        

        loss_red = dict_param_fish_farms_ready[name_corresponding_loss_red]
        
        
        
        # Collect incumbent values in the technosphere matrix
        
        # Incumbent dead
        Fish_dead_0 = -Tech_num_AH_modif[index_dead_fish,column_index]
        
        
        # Incumbent fish outputs and inputs to the growth stage
        dict_fish_outputs = {}
        Fish_input_0=0
        
        for index_input in index_growth_stages+[index_roe]: # The production outputs have a "+" and the fish inputs from other stages have a "-". So to callcualte the weight gain, we can just add all input values
           

            if Tech_num_AH_modif[index_input,column_index] < 0:  # then it's a fish input from a previous stage
                
                Fish_input_0 = -Tech_num_AH_modif[index_input,column_index]

            elif Tech_num_AH_modif[index_input,column_index] > 0:  # It's a fish or roe output
                
                dict_fish_outputs[index_input] = Tech_num_AH_modif[index_input,column_index]



        Total_output_fish_0 = sum(dict_fish_outputs.values()) 

        """ Modif the matrix"""


        # New total biomass
        A = calculate_A(FCR_red_ratio,Fish_dead_0,Total_output_fish_0,Fish_input_0)
        
        
    
        # New Dead
        dead_1 = calculate_F_dead_1(A,Fish_dead_0,Fish_input_0,loss_lev,loss_red,Total_output_fish_0)


        Tech_num_AH_modif[index_dead_fish,column_index]  =  -dead_1
        
        # Fish outputs
        summ=0 # Just for indication
   
        for index_input_2 in dict_fish_outputs:
            
            Fish_0 = dict_fish_outputs[index_input_2]
            Fish_1 = calculate_F_out_x_1(A,Fish_dead_0,Fish_input_0,loss_lev,loss_red,Fish_0,Total_output_fish_0)
            
            summ+=Fish_1
            
            
            Tech_num_AH_modif[index_input_2,column_index]  =  Fish_1
            
            
        
        feed_input = Tech_num_AH_modif[index_feed,column_index]
        new_bio_FCR_AH = -feed_input/(summ+dead_1-Fish_input_0)
        new_eco_FCR_AH = -feed_input/(summ-Fish_input_0)
        

        
    # Calculate excretion

    
    Tech_num_AH_modif_with_excretion = excretion(
              Tech_num_AH_modif,
              index_growth_stages,
              N_P_profile_feed_incumbent,
              P_in_fish,
              N_in_fish,
              biochem_profile_feed_incumbent,
              index_dead_fish,
              fraction_non_ingested,
              digestibility_list,
              excretion_P_removal_efficiency,
              excretion_N_removal_efficiency,
              solid_filter_efficiency,
              CH4_volume_sludge_and_manure,
              share_fish_sludge_in_substrate,
              CH4_LHV,
              ratio_CO2_CH4_biogas_sludge,
              N_manure,
              P_manure,
              fertilizer_substi_digest_N,
              fertilizer_substi_digest_P,
              fertilizer_substi_manure_field_N,
              fertilizer_substi_manure_field_P,
              index_growth_stages_no_filter_laguna,
              index_feed,
              index_roe,
              index_Nemissions,
              index_Pemissions,
              index_sludge,
              index_biogas_updgrade,
              index_N_substitution,
              index_P_substitution,
              index_heat_substitution) 
        

                                                                  
    return Tech_num_AH_modif_with_excretion
    



def calculate_economic_indicators(Tech_Matrix,
                                  index_300g,
                                  index_growing_DK,
                                  index_roe,
                                  index_FU,
                                  index_feed,
                                  index_dead_fish):
    
    """Tecg_Matrix = Numerical matrix"""
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

    





def excretion(Tech_num_AH_modif,
              index_growth_stages,
              N_P_profile_feed_incumbent,
              P_in_fish,
              N_in_fish,
              biochem_profile_feed_incumbent,
              index_dead_fish,
              fraction_non_ingested,
              digestibility_list,
              excretion_P_removal_efficiency,
              excretion_N_removal_efficiency,
              solid_filter_efficiency,
              CH4_volume_sludge_and_manure,
              share_fish_sludge_in_substrate,
              CH4_LHV,
              ratio_CO2_CH4_biogas_sludge,
              N_manure,
              P_manure,
              fertilizer_substi_digest_N,
              fertilizer_substi_digest_P,
              fertilizer_substi_manure_field_N,
              fertilizer_substi_manure_field_P,
              index_growth_stages_no_filter_laguna,
              index_feed,
              index_roe,
              index_Nemissions,
              index_Pemissions,
              index_sludge,
              index_biogas_updgrade,
              index_N_substitution,
              index_P_substitution,
              index_heat_substitution):
 
     '''
     Function which calculates the exrection according to the new fish flows in the tehcnosphere matrix
    
        
     #Outputs : 
              
     Updated technosphere matrix with excretion
             
     '''

    
     Tech_num_AH_modif_with_excretion = Tech_num_AH_modif.copy()

     for growth_stage_index in index_growth_stages:
         
         

         total_feed_AH = - Tech_num_AH_modif_with_excretion[index_feed,growth_stage_index]

         total_weight_gain_economic_AH =  0
         
         # Some growth stages have multiple output of fish, so we make sure to collect all outputs
    
         for index_input in index_growth_stages+[index_roe]: # The produciton outputs have a "+" and the fish inputs from other stages have a "-". So to callcualte the weight gain, we can just add all input values
             
             total_weight_gain_economic_AH += Tech_num_AH_modif_with_excretion[index_input,growth_stage_index]
            
         # Biological weight gain includes the losses ( - as there is a - in the matrix)
         total_weight_gain_biological_AH = total_weight_gain_economic_AH - Tech_num_AH_modif_with_excretion[index_dead_fish,growth_stage_index] 
         

      
         # Total N and P in feed
         """So far we assume no significant change in the fish feed compo"""
         

         N_in_feed_AH = total_feed_AH * N_P_profile_feed_incumbent[0]
         P_in_feed_AH = total_feed_AH * N_P_profile_feed_incumbent[1]

         #Total N and P in fish stock (biological = )
         N_stock_fish_AH = total_weight_gain_biological_AH * N_in_fish
         P_stock_fish_AH = total_weight_gain_biological_AH * P_in_fish
   
         # Feed non ingested 
         

         
         non_ingested_feed_AH = total_feed_AH * fraction_non_ingested
         

         """So far we assume the N and P content modification with microalgae is neglectable"""
         N_in_non_ingested_feed_AH = non_ingested_feed_AH * N_P_profile_feed_incumbent[0]
         P_in_non_ingested_feed_AH = non_ingested_feed_AH * N_P_profile_feed_incumbent[1]
    

         
         #Solid faeces 
         #Total Biochemical class * digestibility
         
         # keep only Lipid, Protein and add phosphorus
         
         biochem_profile_feed_clean_AH = biochem_profile_feed_incumbent[0:4] + [N_P_profile_feed_incumbent[1]]
         
     

         total_excreted_per_biochemical_class_AH = [biochem * (1-digest) * total_feed_AH for biochem,digest in zip(biochem_profile_feed_clean_AH,digestibility_list)]
    
     
         solid_faeces_AH = sum(total_excreted_per_biochemical_class_AH)
         
         
         
         N_in_solid_faeces_AH = total_excreted_per_biochemical_class_AH[1] * 0.16
        
         P_in_solid_faeces_AH = total_excreted_per_biochemical_class_AH[4]
         
         # Total solid
         
         Total_solid_AH =solid_faeces_AH + non_ingested_feed_AH

         # Liquid excretion
         # Input - N_ non ingested - N_fish - Nsolid
         
         N_liquid_AH = N_in_feed_AH - N_in_non_ingested_feed_AH - N_stock_fish_AH - N_in_solid_faeces_AH
       
        
         P_liquid_AH = P_in_feed_AH - P_in_non_ingested_feed_AH - P_stock_fish_AH - P_in_solid_faeces_AH
         

         
    
         # Filter and lagune efficiency
         N_liquid_after_treatment_AH = N_liquid_AH * (1-excretion_N_removal_efficiency)
         P_liquid_after_treatment_AH = P_liquid_AH * (1-excretion_P_removal_efficiency)
         
         Solid_faeces_after_filter_AH = solid_faeces_AH * (1-solid_filter_efficiency)
     
         Solid_faeces_sludge_AH = solid_faeces_AH * solid_filter_efficiency
         
         total_sludge_AH =  Solid_faeces_sludge_AH + non_ingested_feed_AH*solid_filter_efficiency
    
         #
         N_total_without_removal_AH = N_liquid_AH + N_in_solid_faeces_AH + N_in_non_ingested_feed_AH 
         P_total_without_removal_AH = P_liquid_AH + P_in_solid_faeces_AH + P_in_non_ingested_feed_AH 
     
         N_solid_and_liquid_after_removal_AH =  (N_in_solid_faeces_AH + N_in_non_ingested_feed_AH) * (1-solid_filter_efficiency)*(1-excretion_N_removal_efficiency) + N_liquid_after_treatment_AH 
         P_solid_and_liquid_after_removal_AH =  (P_in_solid_faeces_AH + P_in_non_ingested_feed_AH) * (1-solid_filter_efficiency)*(1-excretion_P_removal_efficiency) + P_liquid_after_treatment_AH 
     
         # Content N, P per kilo dm sludge
         N_in_sludge_dw = (N_in_solid_faeces_AH +N_in_non_ingested_feed_AH)*solid_filter_efficiency/total_sludge_AH
         P_in_sludge_dw = (P_in_solid_faeces_AH +P_in_non_ingested_feed_AH)*solid_filter_efficiency/total_sludge_AH
 
         


         # Here calculate input/output sludge treatment
         
         [fertilizer_subst_N_per_kilo_dw_sludge,
          fertilizer_subst_P_per_kilo_dw_sludge,
          MJ_substituted_per_kilo_dw_sludge,
          amount_of_upgrading_act_sludge]=functions.fish_sludge_management(CH4_volume_sludge_and_manure,
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
                           fertilizer_substi_manure_field_P)
                                                          
     # Update the slude management activity
         
         # Upgrading biogas
         Tech_num_AH_modif_with_excretion[index_biogas_updgrade,index_sludge] = -amount_of_upgrading_act_sludge
         
         # Fertilizer subst
         Tech_num_AH_modif_with_excretion[index_N_substitution,index_sludge] = fertilizer_subst_N_per_kilo_dw_sludge

         Tech_num_AH_modif_with_excretion[index_P_substitution,index_sludge] = fertilizer_subst_P_per_kilo_dw_sludge

         # Heat subst
         Tech_num_AH_modif_with_excretion[index_heat_substitution,index_sludge] = MJ_substituted_per_kilo_dw_sludge

                                                        


                                                                 
     # Update the rest of the technosphere matrix
         
         
         if growth_stage_index in index_growth_stages_no_filter_laguna:
             
                 Tech_num_AH_modif_with_excretion[index_Nemissions,growth_stage_index] = -N_total_without_removal_AH
                 Tech_num_AH_modif_with_excretion[index_Pemissions,growth_stage_index] = -P_total_without_removal_AH
                 Tech_num_AH_modif_with_excretion[index_sludge,growth_stage_index] = 0


     
         else :
                 Tech_num_AH_modif_with_excretion[index_Nemissions,growth_stage_index] = -N_solid_and_liquid_after_removal_AH
                 Tech_num_AH_modif_with_excretion[index_Pemissions,growth_stage_index] = -P_solid_and_liquid_after_removal_AH
                 Tech_num_AH_modif_with_excretion[index_sludge,growth_stage_index] = -total_sludge_AH

                 
     return Tech_num_AH_modif_with_excretion









    
    