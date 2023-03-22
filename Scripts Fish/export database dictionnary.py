# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 19:03:54 2021

@author: Pierre Jouannais
pijo@plan.aau.dk


"""

'''Script to export a BW database as a dictionnary in a
 json file and to import it from this file.
'''

import bw2data
import bw2io
from bw2data.parameters import *
import brightway2 as bw
from json import loads, dumps
from ast import literal_eval
import json

import ast






bw.projects.set_current('Fish_project_check') 


# Loading Ecoinvent
Ecoinvent = bw.Database('ecoinvent 3.8 conseq')
# Loading Database
FISHMIC = bw.Database('AH_combi_1')


MICAH = bw.Database('Micro_for_combi')

currentfolder = os.getcwd()

def remap_keys(mapping):
    
    return [{'key':k, 'value': v} for k, v in mapping.iteritems()]


def exchange_to_dict(exchange):
    
    dicti={'amount':exchange['amount'],
           'input':exchange['input'],
           'type':exchange['type']}
    
    return dicti


def export_as_dictionnary_json(bw_database,namedatabase):
    
    dict_database={}
    
    for act in bw_database:

        list_exchanges=[]
        for exc in list(act.exchanges()):
            list_exchanges.append(exchange_to_dict(exc))
        
        dicti_act={'name':act['name'],
                   'exchanges':list_exchanges,
                   'unit':act['unit'],
                   'location':act['location'],
                   }
        dict_database[(namedatabase,act['code'])]=dicti_act

    ready = json.dumps({str(k): v for k, v in dict_database.items()})
    
    # open file for writing, "w" 
    f = open(namedatabase+".json","w")
    
    # write json object to file
    f.write(ready)
    
    # close file
    f.close()





def import_database_data_from_json(jsonfilepath):
      
    with open(jsonfilepath+'.json') as json_file:
        data = json.load(json_file)
      
    database = {ast.literal_eval(k): v for k, v in data.items()}
 
       
    return database   




# Export the database


export_as_dictionnary_json(MICAH,'Micro_for_combi')

export_as_dictionnary_json(FISHMIC,'AH_combi_1')

