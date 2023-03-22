# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 11:38:38 2023

@author: Pierre Jouannais, Department of Planning, DCEA, Aalborg University
pijo@plan.aau.dk




"""


'''
Script fo contribution analysis based on the outputs from the main simulation
'''








import datetime
from time import *
import requests
import pickle
import cProfile
from scipy.integrate import odeint

import os
import sys
# Set working directory to file location 
# (works only when executing the whole file and not only sections (Run Current cell))

# currentfolder=os.path.dirname(os.path.realpath(__file__))
# os.chdir(currentfolder)
# os.getcwd()
#os.chdir('/home/ubuntu/Work_folder/Code_third_paper/Scripts')
#os.chdir("/home/ubuntu/work_folder/Code_third_paper/Scripts")

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

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits import mplot3d




def importpickle(path):
    with open(path, 'rb') as pickle_load:
        obj = pickle.load(pickle_load)
    return obj    





list_meth =[('ReCiPe Midpoint (H)', 'climate change', 'GW100'), 
            ('ReCiPe Midpoint (H)', 'human toxicity', 'HTinf'),
            ('ReCiPe Midpoint (H)', 'freshwater ecotoxicity', 'FETinf'),
            ('ReCiPe Midpoint (H)', 'freshwater eutrophication', 'FE'),
            ('ReCiPe Midpoint (H)', 'marine eutrophication', 'ME'),
             ('TRACI', 'environmental impact', 'eutrophication'),
             ('CML 2001 (superseded)', 'eutrophication potential', 'Eutro.'),
             ('ReCiPe Midpoint (H)', 'terrestrial ecotoxicity', 'TETinf'),
             ('ReCiPe Midpoint (H)', 'fossil depletion', 'FD'),
             ('ReCiPe Midpoint (H)', 'terrestrial acidification', 'TA100'),
             ('ReCiPe Midpoint (H)', 'ozone depletion', 'ODinf'),
             ('ReCiPe Midpoint (H)', 'particulate matter formation', 'PMF')]


# Import results for contribution

list_contrib_A= importpickle("../Outputs_fish/result_contribution_fisho_2_23_842481_dict_update_sce_A.pkl")
 
list_contrib_fish_fry_95perc= importpickle("../Outputs_fish/result_contribution_fisho_2_14_476451_dict_update_sce_stage_fry.pkl")


list_contrib_sfdk2_loss= importpickle("../Outputs_fish/result_contribution_fisho_2_23_561982_dict_update_sce_stage_sfdk2.pkl")





dataframe_contrib= pd.DataFrame(list_contrib_A[0])
dataframe_contrib.columns=["Contribution"]
dataframe_contrib["Process"] = ['Fish_feed',
 'Water_fish_source',
 'Oxygen_fish_prod',
 'Dicopper_oxide_prod',
 'Copper_pyrithione_prod',
 'Diesel_trucks_use',
 'Diesel_boats_use',
 'Roe_prod',
 'Drug_1_prod',
 'Drug_2_prod',
 'Drug_3_prod',
 'Drug_4_prod',
 'Drug_5_prod',
 'Peracetic_acid_aqua_oxide_15_prod',
 'Formalin_24_prod',
 'Sodium_hydroxide_hatronlud_27_7_prod',
 'Sodium_chloride_prod',
 'Hydrochloric_acid_30_prod',
 'Chloramine_prod',
 'Quarternary_ammonium_salt_prod',
 'Electricity_fish_prod',
 'N_emissions_background',
 'P_emissions_background',
 'P source production fish',
 'N source production fish',
 'K source production fish',
 'Heat for AD fish',
 'Mg source production fish',
 'Electricity anaerobic digestion fish',
 'Biogenic_CO2_re-emission',
 'Heat market substitution Fish',
 'biogas purification to biomethane by amino washing']


ax_all = dataframe_contrib.plot.bar(x="Process", y="Contribution", ylim=(-0.5,1),rot=90,title='My Title')


ax_all.set_ylabel("Contribution")
ax_all.get_legend().remove()
ax_all.set_title().remove()

fig_all= ax_all.get_figure()

name_fig_all = "contrib_"+meth+"_"+ identif 
   
fig_all.savefig("../Outputs_fish/Contrib/"+name_fig_all,dpi=500,bbox_inches="tight")


identif = "2402.jpg"

meth_count= -1
for meth_contrib in list_contrib_A:
    
    meth_count+= 1
    
    dataframe_contrib= pd.DataFrame(meth_contrib)
    dataframe_contrib.columns=["Contribution"]
    dataframe_contrib["Process"] = list_FU_combined_names_mc
    
    meth = list_meth[meth_count][-1]

    if max(dataframe_contrib["Contribution"])>1:
        max_=max(dataframe_contrib["Contribution"])
    else:
        max_=1
        
    if min(dataframe_contrib["Contribution"])<0:  
        min_=min(dataframe_contrib["Contribution"])
    else:
        min_=0
    
    
    ax_all = dataframe_contrib.plot.bar(x="Process", y="Contribution", ylim=(min_,max_),rot=90,title=meth)

    
    ax_all.set_ylabel("Contribution")
    ax_all.get_legend().remove()

    fig_all= ax_all.get_figure()

    name_fig_all = "contrib_"+meth+"_"+ identif 
   
    fig_all.savefig("../Outputs_fish/Contrib/"+name_fig_all,dpi=500,bbox_inches="tight")

 
    
    
    
max(dataframe_contrib["Contribution"])
    

identif = "2402.jpg"

meth_count= -1
for meth_contrib in list_contrib_fish_fry_95perc:
    
    meth_count+= 1
    
    dataframe_contrib= pd.DataFrame(meth_contrib)
    dataframe_contrib.columns=["Contribution"]
    dataframe_contrib["Process"] = list_FU_combined_names_mc
    
    meth = list_meth[meth_count][-1]
    
    if max(dataframe_contrib["Contribution"])>1:
        max_=max(dataframe_contrib["Contribution"])
    else:
        max_=1
        
    if min(dataframe_contrib["Contribution"])<0:  
        min_=min(dataframe_contrib["Contribution"])
    else:
        min_=0
    
    ax_all = dataframe_contrib.plot.bar(x="Process", y="Contribution", ylim=(min_,max_),rot=90,title=meth)

    
    ax_all.set_ylabel("Contribution")
    ax_all.get_legend().remove()

    fig_all= ax_all.get_figure()

    name_fig_all = "contrib_fishfryloss95"+meth+"_"+ identif 
   
    fig_all.savefig("../Outputs_fish/Contrib/"+name_fig_all,dpi=500,bbox_inches="tight")

 







identif = "2402.jpg"

meth_count= -1
for meth_contrib in list_contrib_sfdk2_loss:
    
    meth_count+= 1
    
    dataframe_contrib= pd.DataFrame(meth_contrib)
    dataframe_contrib.columns=["Contribution"]
    dataframe_contrib["Process"] = list_FU_combined_names_mc
    
    meth = list_meth[meth_count][-1]
    
    
    if max(dataframe_contrib["Contribution"])>1:
        max_=max(dataframe_contrib["Contribution"])
    else:
        max_=1
        
    if min(dataframe_contrib["Contribution"])<0:  
        min_=min(dataframe_contrib["Contribution"])
    else:
        min_=0
    
    ax_all = dataframe_contrib.plot.bar(x="Process", y="Contribution", ylim=(min_,max_),rot=90,title=meth)

    
    ax_all.set_ylabel("Contribution")
    ax_all.get_legend().remove()

    fig_all= ax_all.get_figure()

    name_fig_all = "contrib_sfdk2_15perc"+meth+"_"+ identif 
   
    fig_all.savefig("../Outputs_fish/Contrib/"+name_fig_all,dpi=500,bbox_inches="tight")

 