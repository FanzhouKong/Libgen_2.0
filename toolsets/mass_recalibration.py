#!/usr/bin/env python
# coding: utf-8

# In[10]:


# from sklearn.linear_model import LinearRegression
from moepy import lowess
import numpy as np
import pandas as pd
from tqdm import tqdm
import toolsets.spectra_operations as so
from toolsets.search import string_search, num_search
def data_recalibrate(data, save_diff = False):
    # automatically recalibrate the whole dataset, based on all mixes
    data_recalibrated = pd.DataFrame()
    diff_raw = []
    diff_recalibrated = []
    # comments = []
    # parent_ion = []
    for mix in tqdm(data['reference_mix'].unique()):
        # poly = PolynomialFeatures(degree = 3)
        data_temp = string_search(data, "reference_mix", mix)
        # data_temp = data.loc[data['reference_Mix label']==i]
        x_temp = data_temp['precursor_mz'].values
        y_temp =  data_temp['reference_precursor_mz'].values
        diff_raw.extend((x_temp-y_temp).tolist())
        model_temp = fit_model(x_temp, y_temp)
        y_pred = model_temp.predict(x_temp.reshape(-1,1))
        if len(x_temp)<=1:
            
            data_temp['calibration']='not recalibrated_'
        else:
            # y_pred = model_temp.predict(poly.fit_transform(x_temp.reshape(-1,1)))
            data_temp['calibration']='recalibrated_'
        diff_recalibrated.extend((y_pred-y_temp).tolist())
        # parent_ion.extend(y_pred.tolist())
        msms_recalibrated = []
        for n in range(len(data_temp)):

            msms_recalibrated.append(mass_recalibrate(model_temp, data_temp.iloc[n]['peaks'] ))
        data_temp['peaks_recalibrated']=msms_recalibrated
        data_recalibrated = pd.concat([data_recalibrated, data_temp], ignore_index = True, axis = 0)
    if save_diff == True:
        data_recalibrated['diff_raw']=diff_raw
        data_recalibrated['diff_recalibrated']=diff_recalibrated
    # data_recalibrated['parent_ion']=parent_ion
    # data_recalibrated['comments']=comments
    data_recalibrated.reset_index(inplace=True, drop=True)
    return(data_recalibrated)




def mass_predicting(mass_measured, model):
    # predicting mass values using constructed linear model
    mass_cali = model.predict(np.array(mass_measured).reshape(-1,1))

    
    mass_cali = [round(num, 6) for num in mass_cali]
#     mass_cali = [str(integer) for integer in mass_cali]
    return(mass_cali)
    

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
def fit_model(x, y):
    # fit the measured (averagemz) to actual value (precursormz) using a linear model
    if(len(x)>1):
        pass
    else:
        # if there is only 1 compound in a mix, do not change it but just return the original value
        x = np.array([1,2,3])
        y = np.array([1,2,3])

        # return(model)
        # print("i am in else")
    # poly = PolynomialFeatures(degree = 3)
    # X_poly = poly.fit_transform(x.reshape(-1,1))
    # poly.fit(X_poly, y)

    # model = LinearRegression()
    # model.fit(X_poly, y)
    model = LinearRegression()
    model.fit(x.reshape(-1,1),y)
    return (model)
    
def mass_recalibrate(model, msms):
    # recalibrate a single mix
    mass, intensity = so.break_spectra(msms)
    mass_recali = mass_predicting(mass, model)
    msms_recali = so.pack_spectra(mass_recali, intensity)
    return(msms_recali)
        





