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


def head_to_tail_plot(msms_1, msms_2, lower=None, upper=None, identity = False, ifnormalized = True):
    mass1, intensity1 = so.break_spectra(msms_1)
    mass1 = [float(x) for x in mass1]
    intensity1 = [float(x) for x in intensity1]
    d = {'m/z':mass1, 'intensity':intensity1}
    msms1 = pd.DataFrame(d)
    max_val = np.max(msms1['intensity'])
    if ifnormalized:
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
    if ifnormalized:
        msms2["normalized_intensity"] = msms2['intensity'] / max_val * 100.0  # normalize intensity to percent
    else:
        msms2["normalized_intensity"]=msms2['intensity']


    msms2['inverse_normalized_intensity']=-msms2['normalized_intensity']
    
    fig = plt.figure(figsize = (7, 5))
    plt.subplots_adjust()
    ax = fig.add_subplot()
    for i in range(len(msms1['m/z'])):
        plt.vlines(x = msms1["m/z"][i], ymin = 0, ymax = msms1["normalized_intensity"][i],color = 'blue')
    for i in range(len(msms2['m/z'])):
        plt.vlines(x = msms2["m/z"][i], ymax = 0, ymin = msms2["inverse_normalized_intensity"][i], color = 'r')
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
    # ax.yaxis.set_ticks(np.arange(start, end + 1, 10))
    ax.grid(False)
    plt.grid(True, axis="y", color='black', linestyle=':', linewidth=0.1)
    return(plt)


# In[10]:


def ms2_plot(msms_1, lower=None, upper=None):
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
        plt.vlines(x = msms1["m/z"][i], ymin = 0, ymax = msms1["normalized_intensity"][i],color = 'blue')
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
    ax.set_ylim(0, 100)
    plt.axhline(y=0, color='black', linestyle='-')
    start, end = ax.get_ylim()
    # ax.yaxis.set_ticks(np.arange(start, end + 1, 10))
    plt.grid(True, axis="y", color='black', linestyle=':', linewidth=0.1)
    return(plt)


# In[17]:


print("i am ms2 spectra plotter, and I have been sideloaded successfully")
print("I have 2 functions, head to tail plot and ms2 plot")
