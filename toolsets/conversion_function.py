#!/usr/bin/env python
# coding: utf-8

# In[1]:


import requests
import numpy as np
import pandas as pd

def get_something_from_pubchem(inputt,content, something1):


    r = requests.get(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{inputt}/{content}/property/{something1}/JSON').json()
    if r =={'Fault': {'Code': 'PUGREST.NotFound',
  'Message': 'No CID found',
  'Details': ['No CID found that matches the given InChI key']}}:
#         print("i am in if loop")
        r = requests.get(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{inputt}/{content[0:14]}/property/{something1}/JSON').json()
    if r =={'Fault': {'Code': 'PUGREST.NotFound',
  'Message': 'No CID found',
  'Details': ['No CID found that matches the given InChI key(s)']}}:
        return(np.NaN)
    else:
        return r['PropertyTable']['Properties'][0][something1]


# formulas.append(cf.get_something_from_pubchem("inchikey", row['InChIKey'], "MolecularFormula"))

