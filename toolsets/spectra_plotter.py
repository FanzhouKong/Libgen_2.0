#!/usr/bin/env python
# coding: utf-8

# In[6]:


import pandas as pd
import numpy as np

from matplotlib import rcParams

import toolsets.spectra_operations as so

import matplotlib.pyplot as plt

import seaborn as sns

import plotly.express as px


# In[7]:



# In[169]:


# In[3]:


def head_to_tail_plot(msms_1, msms_2,pmz1=None, pmz2= None, color1 = None, color2 = None,lower=None, upper=None, identity = False, normalize = True, ifNIST = False,savepath = None, show= True):
    if ifNIST ==True:
        msms_1 = so.convert_nist_to_string(msms_1)
        msms_2 = so.convert_nist_to_string(msms_2)
    mass1, intensity1 = so.break_spectra(msms_1)
    mass1 = [float(x) for x in mass1]
    intensity1 = [float(x) for x in intensity1]
    d = {'m/z':mass1, 'intensity':intensity1}
    msms1 = pd.DataFrame(d)
    max_val = np.max(msms1['intensity'])
    if normalize:
        msms1["normalized_intensity"] = msms1['intensity'] / max_val * 100.0  # normalize intensity to percent
    else:
        msms1["normalized_intensity"]=msms1['intensity']

    mass2, intensity2 = so.break_spectra(msms_2)
    mass2 = [float(x) for x in mass2]
    intensity2 = [float(x) for x in intensity2]
    d = {'m/z':mass2, 'intensity':intensity2}

    msms2 = pd.DataFrame(d)
    if identity == True:
        max_val = np.max(msms1['intensity'])
    else:
        max_val = np.max(msms2['intensity'])
    if normalize:
        msms2["normalized_intensity"] = msms2['intensity'] / max_val * 100.0  # normalize intensity to percent
    else:
        msms2["normalized_intensity"]=msms2['intensity']


    msms2['inverse_normalized_intensity']=-msms2['normalized_intensity']
    
    fig = plt.figure(figsize = (7, 5))
    plt.subplots_adjust()
    ax = fig.add_subplot()
    for i in range(len(msms1['m/z'])):
        if color1 == None:
            plt.vlines(x = msms1["m/z"][i], ymin = 0, ymax = msms1["normalized_intensity"][i],color = 'blue')
        elif color1 != None:
            plt.vlines(x = msms1["m/z"][i], ymin = 0, ymax = msms1["normalized_intensity"][i],color = color1)
    if pmz1 != None:
        plt.vlines(x = pmz1, ymin = 0, ymax = 100,color = 'grey', linestyle='dashed')
    for i in range(len(msms2['m/z'])):
        if color2 ==None:
            plt.vlines(x = msms2["m/z"][i], ymax = 0, ymin = msms2["inverse_normalized_intensity"][i], color = 'r')
        elif color2 != None:
            plt.vlines(x = msms2["m/z"][i], ymax = 0, ymin = msms2["inverse_normalized_intensity"][i], color = color2)
    if pmz2 != None:
        plt.vlines(x = pmz2, ymin = -100, ymax = 0,color = 'grey', linestyle='dashed')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_xlabel(r"$m/z$")
    ax.set_ylabel(r"$Intensity\,[\%]$")
    plt.xticks(rotation='vertical')
    start, end = ax.get_xlim()
    start, end = ax.get_xlim()
    if(lower!=None and upper!= None):
        ax.set_xlim(lower, upper)
    ax.set_ylim(-100, 100)
    plt.axhline(y=0, color='black', linestyle='-')
    start, end = ax.get_ylim()
    plt.tight_layout()
    ax.set_facecolor("none")
    # ax.yaxis.set_ticks(np.arange(start, end + 1, 10))
    ax.grid(False)
    plt.grid(True, axis="y", color='black', linestyle=':', linewidth=0.1)
    ax.grid(False)
    if savepath != None:
        plt.savefig(savepath, dpi = 300,facecolor = 'none', edgecolor = 'none')
    if show==True:
        return(plt)
    else:
        return()


# In[10]:
from toolsets.API_gets import GNPS_get
from rdkit import Chem
import os
from toolsets.API_gets import pubchem_get
from mimas.external.features_by_alphapept.load_mzml_data import load_mzml_data
from tqdm import tqdm
from toolsets.file_io import find_files
from toolsets.spectra_operations import break_spectra, pack_spectra
import matplotlib.pyplot as plt
def molplot_from_inchikey(inchikey):
    smile_temp = pubchem_get(content=inchikey)
    mol = Chem.MolFromSmiles(smile_temp)
    mol
    return(mol)
def fast_EIC(mix,mzml_dir, precursormz,rt_start = -1, rt_end = -1, vlines_location = []):
    filename = find_files(mzml_dir, mix+'.mzml')
    mzml_file_path = os.path.join(mzml_dir, filename[0])
    test_mzml = load_mzml_data(mzml_file_path)
    EIC_intensity = []
    ms1_ppm=20
    for i in tqdm(range(0, len(test_mzml['rt_list_ms1']))):
        # ms1s.append(_extract_MS1(test_mzml, i))
        chopped_ms1_temp = chop_msms(_extract_MS1(test_mzml, i),precursormz*((1 - ms1_ppm * 1e-6)), precursormz*((1 + ms1_ppm * 1e-6)))
        mass_chopped, intensity_chopped = break_spectra(chopped_ms1_temp)
        EIC_intensity.append(np.sum(intensity_chopped) )
    fig, ax = plt.subplots(
    figsize = (9, 5)
                      )
    ax= sns.lineplot(test_mzml['rt_list_ms1'], EIC_intensity)
    if rt_start != -1 and rt_end != -1:
        ax.set_xlim(rt_start, rt_end)
        ax.set_ylim(0, np.max(EIC_intensity)+100)
    if len(vlines_location)>0:
        for position in vlines_location:
            plt.axvline(x = position, color = 'red')
    ax.grid(False)

def chop_msms(msms, lowest_allowed, highest_allowed):
    import bisect
    import toolsets.spectra_operations as so
    msms = so.sort_spectra(msms)
    mass, intensity = so.break_spectra(msms)
    upper_allowed = bisect.bisect_right(mass, highest_allowed)
    lower_allowed = bisect.bisect_left(mass, lowest_allowed)
    # pep_index = mass.index(parention)
    mass_choppped = mass[lower_allowed:upper_allowed]
    intensity_chopped = intensity[lower_allowed:upper_allowed]
    if len(mass_choppped)>0 and len(intensity_chopped)>0:
        return(so.pack_spectra(mass_choppped, intensity_chopped))
    else:
        return(so.pack_spectra([0],[0]))
def _extract_MS1(mzml, scan_number):
    from toolsets.spectra_operations import pack_spectra
    ms1_1 = pack_spectra(mzml['mass_list_ms1'][mzml['indices_ms1'][scan_number]:mzml['indices_ms1'][scan_number+1]],
                     mzml['int_list_ms1'][mzml['indices_ms1'][scan_number]:mzml['indices_ms1'][scan_number+1]])
    return (ms1_1)
def ms2_plot(msms_1, pmz1 = None, lower=None, upper=None, savepath = None, color = 'blue'):
    mass1, intensity1 = so.break_spectra(msms_1)
    mass1 = [float(x) for x in mass1]
    intensity1 = [float(x) for x in intensity1]
    d = {'m/z':mass1, 'intensity':intensity1}
    msms1 = pd.DataFrame(d)
    max_val = np.max(msms1['intensity'])
    msms1["normalized_intensity"] = msms1['intensity'] / max_val * 100.0  # normalize intensity to percent
    
    fig = plt.figure(figsize = (7, 5))
    plt.subplots_adjust()
    ax = fig.add_subplot()
    for i in range(len(msms1['m/z'])):
        plt.vlines(x = msms1["m/z"][i], ymin = 0, ymax = msms1["normalized_intensity"][i],color = color, linewidth=3)
    if pmz1 != None:
        plt.vlines(x = pmz1, ymin = 0, ymax = 100,color = 'grey', linestyle='dashed')
    # plt.legend()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_xlabel(r"$m/z$")
    ax.set_ylabel(r"$Intensity\,[\%]$")
    plt.xticks(rotation='vertical')
    start, end = ax.get_xlim()
    # start, end = ax.get_xlim(), 
    if(lower!=None and upper!= None):
        ax.set_xlim(lower, upper)
    ax.set_ylim(0, 100)
    plt.axhline(y=0, color='black', linestyle='-')
    start, end = ax.get_ylim()
    # ax.yaxis.set_ticks(np.arange(start, end + 1, 10))
    plt.grid(True, axis="y", color='black', linestyle=':', linewidth=0.1)
    ax.grid(False)
    ax.set_facecolor("none")
    if savepath != None:
        plt.savefig(savepath, dpi = 300,transparent=True,facecolor = 'none', edgecolor = 'none')

    return(plt)

def ms_figure_simple(msms1, savepath = None, color = 'blue'):
    import toolsets.spectra_operations as so
    # msms1 = (data_temp.iloc[0]['peaks'])
    mass, intensity =so.break_spectra(msms1)
    d = {'m/z':mass, 'intensity':intensity}
    msms1 = pd.DataFrame(d)
    max_val = np.max(msms1['intensity'])
    msms1["normalized_intensity"] = msms1['intensity'] / max_val * 100.0  # normalize intensity to percent
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize = (7, 5),facecolor="white")
    # fig.figure()
    plt.subplots_adjust()
    ax = fig.subplots()
    for i in range(len(msms1['m/z'])):
        plt.vlines(x = msms1["m/z"][i], ymin = 0, ymax = msms1["normalized_intensity"][i],color = color, linewidth = 3)
    # ax.set_xlim(150, 210)
    ax.set_ylim(0, 100)
    ax.set(xticklabels=[])
    ax.set(yticklabels=[])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['bottom'].set(color = 'black', linewidth = 2)
    ax.spines['left'].set_visible(False)
    fig.tight_layout()
    ax.grid(False)
    if savepath == None:
        return(plt)
    else:
        plt.savefig(savepath, dpi = 300,transparent=True, edge_color = 'none')

# In[17]:


print("i am ms2 spectra plotter, and I have been sideloaded successfully")
print("I have 2 functions, head to tail plot and ms2 plot")
