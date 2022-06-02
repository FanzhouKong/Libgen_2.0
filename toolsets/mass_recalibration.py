#!/usr/bin/env python
# coding: utf-8

# In[10]:


from sklearn.linear_model import LinearRegression
import numpy as np
import pandas as pd
import toolsets.spectra_operations as so

def mass_predicting(mass_measured, lm_temp):
    # predicting mass values using constructed linear model
    mass_cali = lm_temp.predict(np.array(mass_measured, dtype = float ).reshape(-1,1)) 
    mass_cali = [round(num, 6) for num in mass_cali]
#     mass_cali = [str(integer) for integer in mass_cali]
    return(mass_cali)
    

def fit_model(measured, actual):
    # fit the measured (averagemz) to actual value (precursormz) using a linear model
    if(len(measured)>1):
        x = measured.values.reshape(-1,1)
        y = actual.values
        lm = LinearRegression().fit(x,y)
        return(lm)
    else:
        # if there is only 1 compound in a mix, do not change it but just return the original value
        x = np.array([1,2,3]).reshape(-1,1)
        y = np.array([1,2,3])
        lm = LinearRegression().fit(x,y)
        return(lm)
    
def mass_recalibrate(lm, msms):
    # recalibrate a single mix
    mass, intensity = so.break_spectra(msms)
    mass_recali = mass_predicting(mass, lm)
    msms_recali = so.pack_spectra(mass_recali, intensity)
    return(msms_recali)
        

def data_recalibrate(data):
    # automatically recalibrate the whole dataset, based on all mixes
    data_recalibrated = pd.DataFrame()
    for i in data['mix_label'].unique():
        data_temp = data.loc[data['mix_label']==i]
        x_temp = data_temp['Average_mz']
        y_temp = data_temp['PRECURSORMZ']
        lm_temp = fit_model(x_temp, y_temp)
        msms_recalibrated = []
        for n in range(len(data_temp)):
            msms_recalibrated.append(mass_recalibrate(lm_temp, data_temp.iloc[n]['msms']))
        data_temp['msms_recalibrated']=msms_recalibrated
        data_recalibrated = pd.concat([data_recalibrated, data_temp], ignore_index = True, axis = 0)
    return(data_recalibrated)


def data_recalibrate_precursor(data):
    # automatically recalibrate the whole dataset, based on all mixes
    data_recalibrated = pd.DataFrame()
    diff_raw = pd.Series()
    diff_recalibrated = pd.Series()
    for i in data['mix_label'].unique():
        data_temp = data.loc[data['mix_label']==i]
        x_temp = data_temp['Average_mz']
        y_temp = data_temp['PRECURSORMZ']
        diff_raw = diff_raw.append((x_temp-y_temp))
        lm_temp = fit_model(x_temp, y_temp)
        y_pred = lm_temp.predict(np.array(x_temp, dtype = float ).reshape(-1,1))
        diff_recalibrated=diff_recalibrated.append((x_temp-y_pred))
    # return (, )
        msms_recalibrated = []
        for n in range(len(data_temp)):
            msms_recalibrated.append(mass_recalibrate(lm_temp, data_temp.iloc[n]['msms']))
        data_temp['msms_recalibrated']=msms_recalibrated
        data_recalibrated = pd.concat([data_recalibrated, data_temp], ignore_index = True, axis = 0)
        # break
    data_recalibrated['diff_raw']=diff_raw.tolist()
    data_recalibrated['diff_recalibrated']=diff_recalibrated.tolist()
    return(data_recalibrated)
    # data_recalibrated['diff_raw']=diff_raw
    # data_recalibrated['diff_recalibrated']=diff_recalibrated
    # return(data_recalibrated)


# In[2]:


print("I am mass recalibration, usage: mass_recalibrate(data)")
print("the data column must have columns of mix_label, Average_mz, PRECURSORMZ, msms")
print("the msms should in a string format, e.g. mass1\tintensity1\nmass2\tintensity2\n....")
print("the recalibrated column would be msms_recalibrated")

