#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import toolsets.mass_to_formula as mtf
import re
import pandas as pd
import spectral_entropy as se
import itertools
import numpy as np
import scipy.stats
import os
import warnings
from tqdm import tqdm
import math
warnings.filterwarnings("ignore")
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import toolsets.denoising_related_functions as de
def sort_spectra(msms):
    mass, intensity = break_spectra(msms)
    mass_sorted, intensity_sorted = zip(*sorted(zip(mass, intensity)))
    return(pack_spectra(list(mass_sorted), list(intensity_sorted)))

def duplicate_handling(data, typeofmsms,mass_error = 0.01, ifppm = False, ifnormalize = True, method = 'weighedaverage'):
    data_unique = data.drop_duplicates(subset=['key'])
    consensus_msms = []
    for key in tqdm(data_unique['key'].unique()):
        data_temp =  data.loc[data['key']==key]
        if method =='weighedaverage':
            consensus_msms.append(weighted_average_spectra(data_temp, tolerance=mass_error, ppm = ifppm, typeofmsms =typeofmsms, ifnormalize=ifnormalize))
        elif method =='consensus':
            consensus_msms.append(make_consensus_spectra(data_temp, tolerance=mass_error, ppm = ifppm, typeofmsms =typeofmsms, ifnormalize=ifnormalize))
        else:
            print('please provide a valid duplicate handling method')

    data_unique['msms']=consensus_msms
    # print("i am done processing the concensus spectra")
    return(data_unique)
import toolsets.denoising_related_functions as de

def denoising(data, typeofmsms, mass_error = 0.01, ifppm = False):
    msms_consensus_denoised = []
    for index, row in tqdm(data.iterrows(), total = data.shape[0]):
    # break
        try:
            msms_consensus_denoised.append(de.denoise_blacklist(row, typeofmsms=typeofmsms, mass_error =mass_error, ppm = ifppm))
        except:
            msms_consensus_denoised.append(row['msms'])
    data['msms_denoised']=msms_consensus_denoised
    return(data)


def denoising_evaluation(data, msms1 = 'msms', msms2 = 'msms_denoised', min_explained_intensity = 80, allowed_max_unassigned_intensity = 30):
    explained_intensity = []
    max_unassigned_intensity = []
    for index, row in (data.iterrows()):
        explained_intensity.append(calculate_explained_intensity(row[msms1], row[msms2]))
        max_unassigned_intensity.append(identify_max_unassigned_intensity(row[msms1], row[msms2]))
    data['explained_intensity']=explained_intensity
    data['max_unassigned_intensity']=max_unassigned_intensity
    evaluations = []
    # print(min_explained_intensity/10)
    for index, row in (data.iterrows()):
        if row['explained_intensity']<min_explained_intensity/100 and row['max_unassigned_intensity']>allowed_max_unassigned_intensity:
            evaluations.append('flagged: poor quality')
        elif row['explained_intensity']<min_explained_intensity/100:
            evaluations.append('flagged:low assigned intensity')
        elif row['max_unassigned_intensity']>allowed_max_unassigned_intensity:
            evaluations.append('flagged: high unassigned intensity')
        else:
            evaluations.append('good quality')
    data['evaluations']=evaluations
    return(data)


def bin_spectra(msms, precursormz,tolerance = 0.01, ifnormalize = False, ppm = False):
    if ppm:
        if float(precursormz) <400:
            tol = 0.004
        else:
            tol = float(precursormz)*(tolerance/1E6)
    else:
        tol = tolerance
    # if ifnormalize == True:
    #     msms = normalize_spectra(msms)
    mass_temp, intensity_temp = break_spectra(msms)
    bin_left = pd.DataFrame({'mass': mass_temp, 'intensity': intensity_temp})
    mass_bin = []
    intensity_bin = []
    while(len(bin_left)>0):

        bin, bin_left = make_bin(bin_left, tol)

        # mass_bin.append(round(bin['mass'].loc[bin['intensity']== bin['intensity'].max()], 6))
        # max_index = bin['intensity'].idxmax()
        mass_bin.append(round(bin['mass'].median(), 6))
        # mass_bin.append(round(bin['mass'][max_index], 6))
        # intensity_bin.append(round(bin['intensity'].max(),6))
        intensity_bin.append(round(bin['intensity'].sum(),6))
    msms_bin = sort_spectra(pack_spectra(mass_bin, intensity_bin))
    if ifnormalize:
        return(normalize_spectra(msms_bin))
    else:
        return(msms_bin)


def make_bin(bin_left, tol):
    max_index = bin_left['intensity'].idxmax()
    bin = bin_left[bin_left['mass'].between(bin_left['mass'][max_index]-tol, bin_left['mass'][max_index]+tol, inclusive=True)]
    bin_left_return = bin_left[~bin_left['mass'].between(bin_left['mass'][max_index]-tol, bin_left['mass'][max_index]+tol, inclusive=True)]
    return(bin, bin_left_return)



def make_composite_spectra(data_subset, typeofmsms = 'msms', tolerance = 0.01, ifnormalize = False, ppm = False):
    precursormz = float(data_subset.iloc[0]['PRECURSORMZ'])
    if ppm:
        if precursormz <400:
            tol = 0.004
        else:
            tol = precursormz*(tolerance/1E6)
    else:
        tol = tolerance
    mass_com = []
    intensity_com = []
    for index, row in data_subset.iterrows():
        # msms_bin_temp = bin_spectra(row[typeofmsms],precursormz, tol,ifnormalize=ifnormalize, ppm =ppm)

        mass_temp, intensity_temp = break_spectra(normalize_spectra(row[typeofmsms]) )
        mass_com.extend(mass_temp)
        intensity_com.extend(intensity_temp)
    msms_com = sort_spectra(pack_spectra(mass_com, intensity_com))
    msms_com = bin_spectra(msms_com,precursormz, tol, ppm =ppm, ifnormalize=ifnormalize)
    return((msms_com))







def weighted_average_spectra(data_subset, typeofmsms = 'msms', tolerance = 0.01, ifnormalize = True, ppm = False, ifold = False):
    # if(len(data_subset)<2):
    #     print("you cannot make weighted average spectra using only 1 spectra")
    #     return(np.NAN)
    precursormz = float(data_subset.iloc[0]['PRECURSORMZ'])
    if ppm:
        if precursormz <400:
            tol = 0.004
        else:
            tol = precursormz*(tolerance/1E6)
    else:
        tol = tolerance
    # msms_com = []
    mass_com = []
    intensity_com = []
    ms1_intensity_com =[]
    for index, row in data_subset.iterrows():
        msms_bin_temp = bin_spectra(row[typeofmsms],precursormz, tol,ifnormalize=ifnormalize, ppm =ppm)
        mass_temp, intensity_temp = break_spectra(msms_bin_temp)
        mass_com.extend(mass_temp)
        intensity_com.extend(intensity_temp)
        ms1_intensity_com.extend([row['intensity']]*len(mass_temp))
    bin_left = pd.DataFrame({'mass': mass_com, 'intensity': intensity_com, 'ms1_intensity':ms1_intensity_com})
    # return(bin_left)
    mass_consensus = []
    intensity_consensus =[]

    while(len(bin_left)>0):
        bin, bin_left = make_bin(bin_left, tol)
        if ifold:
            sum =bin['ms1_intensity'].unique().sum()
        else:
            sum =bin['ms1_intensity'].sum()
        temp_mass = bin['mass']*bin['ms1_intensity']/sum
        temp_intensity =bin['intensity']*bin['ms1_intensity']/sum
        mass_consensus.append(round(temp_mass.sum(),6))
        intensity_consensus.append(round(temp_intensity.sum(),6))

    msms_consensus = sort_spectra(pack_spectra(mass_consensus, intensity_consensus))
    msms_consensus = bin_spectra(msms_consensus,precursormz,ifnormalize=ifnormalize)
    return(msms_consensus)



def make_consensus_spectra(data_subset, typeofmsms = 'msms', tolerance = 0.01, ifnormalize = False, ppm = False):
    # if(len(data_subset)<2):
    #     print("you cannot make weighted average spectra using only 1 spectra")
    #     return(np.NAN)
    precursormz = float(data_subset.iloc[0]['PRECURSORMZ'])

    mass_consensus = []
    intensity_consensus =[]
    mass_com =[]
    intensity_com=[]
    if ppm:
        if float(precursormz) <400:
            tol = 0.004
        else:
            tol = float(precursormz)*(tolerance/1E6)
    else:
        tol = tolerance
    # for msms in msms_s:
    #     msms_bin_temp = bin_spectra(msms,precursormz, tolerance=tol,ifnormalize=ifnormalize, ppm =ppm)
    #     mass_temp, intensity_temp = break_spectra(msms_bin_temp)
    #     mass_com.extend(mass_temp)
    #     intensity_com.extend(intensity_temp)
    for index, row in data_subset.iterrows():
        msms_bin_temp = bin_spectra(row[typeofmsms],precursormz, tol,ifnormalize=ifnormalize, ppm =ppm)
        mass_temp, intensity_temp = break_spectra(msms_bin_temp)
        mass_com.extend(mass_temp)
        intensity_com.extend(intensity_temp)
        # ms1_intensity_com.extend([row['intensity']]*len(mass_temp))
    bin_left = pd.DataFrame({'mass': mass_com, 'intensity': intensity_com})

    while(len(bin_left)>0):
        bin, bin_left = make_bin(bin_left, tol)
        mass_consensus.append(round(bin['mass'].median(), 6))
        intensity_consensus.append(round(bin['intensity'].median(),6))
    msms_consensus = sort_spectra(pack_spectra(mass_consensus, intensity_consensus))
    msms_consensus = normalize_spectra(msms_consensus)
    return(msms_consensus)






def denoising_by_threshold(msms, threshold = 1, need_normalized = False):
    mass_raw, intensity_raw = break_spectra(msms)
    if need_normalized == True:
        intensity_normalized = [x/max(intensity_raw) for x in intensity_raw]
    idx=[index for (index, number) in enumerate(intensity_normalized) if number > threshold]
    intensity_updated = [intensity_normalized[i] for i in idx]
    mass_updated = [mass_raw[i] for i in idx]
    msms_updated = pack_spectra(mass_updated, intensity_updated)
    return(msms_updated)




def break_spectra(spectra):
    split_msms = re.split('\t|\n',spectra)
    intensity = split_msms[1:][::2]
    mass = split_msms[::2]
    mass = [float(item) for item in mass]
    intensity = [float(item) for item in intensity]
    return(mass, intensity)

def pack_spectra(mass, intensity):
    intensity_return = [str(inten) + '\n' for (inten) in (intensity[:-1])]
    intensity_return.append(str(intensity[-1]))
    mass_cali_tab = [str(mas) + '\t' for (mas) in mass]
    list_temp = [None]*(len(mass_cali_tab)+len(intensity_return))
    list_temp[::2] = mass_cali_tab
    list_temp[1::2] = intensity_return
    list_temp = ''.join(list_temp)
    return(list_temp)


def convert_nist_to_string(msms):
    mass = []
    intensity = []
    for n in range(0, len(msms)):
        mass.append(msms[n][0])
        intensity.append(msms[n][1])
    return(pack_spectra(mass,intensity))


def convert_string_to_nist(msms):
    spec_raw = np.array([x.split('\t') for x in msms.split('\n')], dtype=np.float32)
    return(spec_raw)

def normalized_entropy(msms, order = 4):
    npeak = num_peaks(msms)
    return ((scipy.stats.entropy(convert_string_to_nist(msms)[:, 1])/math.log(npeak))**order)
    # return(((scipy.stats.entropy(convert_string_to_nist(msms)[:, 1]))/math.log(npeak))^order)

def entropy_similarity_default(msms1, msms2, typeofmsms='msms', threshold = 0.01):
    # return(se.similarity(convert_string_to_nist(msms1), convert_string_to_nist(msms2), 'entropy', ms2_da = 0.01, need_clean_spectra = True, need_normalize_result = True))
    return(se.similarity(convert_string_to_nist(msms1), convert_string_to_nist(msms2), 'entropy',
                                     ms2_da = threshold, need_clean_spectra = True, need_normalize_result = True))



def average_entropy_calculation(data_temp, typeofmsms = "msms"):
    if len(data_temp)==1:
        # print('you just cannot calculate entropy based on 1 spectrum!!!')
        return(np.NaN)
    else:
        entropy_temp = []
        combinations_object =itertools.combinations(data_temp[typeofmsms], 2)
        for n in combinations_object:
            entropy_temp.append(se.similarity(convert_string_to_nist(n[0]), convert_string_to_nist(n[1]), 'entropy', ms2_da = 0.01, need_clean_spectra = True, need_normalize_result = True))
        return(sum(entropy_temp)/len(entropy_temp))

def average_entropy_dataframe(data_pfp, typeofmsms):
    average_entropy = []
    for i in data_pfp['key'].unique():
        data_temp = data_pfp.loc[data_pfp['key']==i,:]
        entropy_temp = average_entropy_calculation(data_temp, typeofmsms)
        average_entropy.append([entropy_temp]*len(data_temp))
    flat = [x for sublist in average_entropy for x in sublist]
    data_pfp['Average_Entropy']= flat
    return(data_pfp)

# def duplicate_handling(data):
#     data_final = pd.DataFrame()
#     for i in data['key'].unique():
#         temp_df = data.loc[data.key == i,:]
#         if(len(temp_df) ==1):
#             data_final = pd.concat([data_final, temp_df], ignore_index = True, axis = 0)
#         else:
#             if(temp_df.iloc[0]['Average_Entropy'] >= 0.5):
#                 data_final = pd.concat([data_final, temp_df[temp_df['intensity'] == temp_df['intensity'].max()]], ignore_index = True, axis = 0)
#     return(data_final)




def num_peaks(msms):
    mass, intensity = break_spectra(msms)
    return(len(mass))

from operator import itemgetter

# still building
# def evaluate_spectra(msm1, msms2):
#     mass_raw, intensity_raw = break_spectra(msms1)
#     mass_dr, intensity_dr = break_spectra(msms2)
#     mass_raw_fl = [float(x) for x in mass_raw]
#     mass_dr_fl = [float(x) for x in mass_dr]
#     diff_index = [i for i, item in enumerate(mass_raw_fl) if item not in set(mass_dr_fl)]
#     mass_diff = list(itemgetter(*diff_index)(mass_raw))
#     intensity_diff = list(itemgetter(*diff_index)(intensity_raw))
#     rel_intensity_diff = [number / max([float(x) for x in intensity_raw])*100 for number in [float(y) for y in intensity_diff]]
#     rel_intensity_kept = [number / max([float(x) for x in intensity_raw])*100 for number in [float(y) for y in intensity_dr]]
#     return(mass_diff,rel_intensity_diff,rel_intensity_kept)

def normalize_spectra(msms):
    mass, intensity = break_spectra(msms)
    mass_fl = [float(x) for x in mass]
    # intensity_fl = [float(x) for x in intensity]
    intensity_rel = [number / max([float(x) for x in intensity])*100 for number in [float(y) for y in intensity]]
    intensity_rel = [round(number, 6) for number in intensity_rel]
    return(pack_spectra(mass_fl, intensity_rel))


# def comparing_spectrum(msms1, msms2):
#     if(num_peaks(msms1)<num_peaks(msms2)):
#         temp_msms = msms1
#         msms1 = msms2
#         msms2 = temp_msms
#     mass_raw, intensity_raw = break_spectra(msms1)
#     mass_dr, intensity_dr = break_spectra(msms2)
#     mass_raw_fl = [float(x) for x in mass_raw]
#     mass_dr_fl = [float(x) for x in mass_dr]
#     diff_index = [i for i, item in enumerate(mass_raw_fl) if item not in set(mass_dr_fl)]
#     mass_diff = list(itemgetter(*diff_index)(mass_raw))
#     intensity_diff = list(itemgetter(*diff_index)(intensity_raw))
#     rel_intensity_diff = [number / max([float(x) for x in intensity_raw])*100 for number in [float(y) for y in intensity_diff]]
#     rel_intensity_kept = [number / max([float(x) for x in intensity_raw])*100 for number in [float(y) for y in intensity_dr]]
#     return(mass_diff,rel_intensity_diff,rel_intensity_kept)




# below is some evaluation functions
def calculate_explained_intensity(msms1, msms2):
    if(num_peaks(msms1)<num_peaks(msms2)):
        temp_msms = msms1
        msms1 = msms2
        msms2 = temp_msms
    mass_raw, intensity_raw = break_spectra(msms1)
    mass_dr, intensity_dr = break_spectra(msms2)
    return(sum(intensity_dr)/sum(intensity_raw))

def identify_max_unassigned_intensity(msms1, msms2):
    if(num_peaks(msms1)<num_peaks(msms2)):
        temp_msms = msms1
        msms1 = msms2
        msms2 = temp_msms
    mass_raw, intensity_raw = break_spectra(msms1)
    mass_de, intensity_de = break_spectra(msms2)
    diff_index = [i for i, item in enumerate(mass_raw) if item not in set(mass_de)]
    # return(diff_index)
    if(len(diff_index)>1):
        # print("i am in wrong if")
        intensity_diff = list(itemgetter(*diff_index)(intensity_raw))
        return(max(intensity_diff))
    elif(len(diff_index)==1):
        # print("i am in right if")
        return(intensity_raw[diff_index[0]])
    else:
        return(0)
    #
    # return(mass_diff)







def export_library(data_dup,output_location, typeofmsms='msms', ifcollision_energy = False):
    entry = ''
    for index, row in data_dup.iterrows():
        entry = entry + 'Name: ' + row['NAME'] + '\n'
        entry = entry +'Spectrum_type: '+row['Spectrum_type']+ '\n'
        entry = entry + 'PrecursorMZ: ' + str(row['PRECURSORMZ']) + '\n'
        entry = entry + 'InChIKey: ' + row['InChIKey'] + '\n'
        entry = entry + 'Formula: ' + row['Formula'] + '\n'
        entry = entry + 'ExactMass: ' + str(row['ExactMass']) + '\n'
        entry = entry + 'Precursor_type: ' + row['Adduct'] + '\n'
        if ifcollision_energy:
            entry = entry + 'Collision_enerty: ' + str(row['Collision_energy']) + '\n'
        # entry = entry + 'RETENTIONTIME: ' + str(row['RETENTIONTIME']) + '\n'
        entry = entry+'Ion_mode: '+row['Ion_mode']+ '\n'
        entry = entry + 'Comment: ' + str(row['Comment']) + '\n'
        entry = entry + 'Num peaks: ' + str(num_peaks(row[typeofmsms])) + '\n'
        entry = entry + row[typeofmsms]
        # entry = entry +str(row['count'])
        entry = entry + '\n'
        entry = entry + '\n'

    #open text file
    text_file = open(output_location, "w",encoding='utf-8')
     
    #write string to file
    text_file.write(entry)
     
    #close file
    text_file.close()

def export_ms_sirius(row, output):
    # if row['Adduct'][-1]=='+':
    #         charge = '1+'
    #     else:
    #         charge = '1-'
    mass_1, intensity_1 = break_spectra(row['ms1'])
    pep_mass =de.find_parention(mass_1,intensity_1, row['PRECURSORMZ'])
    entry = ''
    entry = entry + '>compound '+str(row['NAME'])+'\n'
    entry = entry + '>parentmass '+str((pep_mass))+'\n'
    entry = entry + '>ionization '+str((row['Adduct']))+'\n'
    entry = entry +'\n'
    entry = entry + '>collision 35' + '\n'
    entry = entry+(row['msms'])+'\n'
    entry = entry +'\n'
    entry = entry +'\n'
    entry = entry + '>ms1peaks' + '\n'
    entry = entry+(row['ms1'])+'\n'
    entry = entry +'\n'
    text_file = open(output+'.ms', "w",encoding='utf-8')
    text_file.write(entry)
    text_file.close()



def export_mgf_sirius(inputfile, output):

    entry = ''
    for index, row in inputfile.iterrows():

        mass_1, intensity_1 = break_spectra(row['ms1'])
        if row['Adduct'][-1]=='+':
            charge = '1+'
        else:
            charge = '1-'
        pep_mass =de.find_parention(mass_1,intensity_1, row['PRECURSORMZ'])
        # output = os.path.join(output_dir, row['NAME']+'.mgf')
        # ms1
        entry = entry + 'BEGIN IONS'+'\n'
        entry = entry + 'PEPMASS='+str(pep_mass)+'\n'
        entry = entry + 'MSLEVEL=1'+'\n'
        entry = entry+'CHARGE=' + charge +'\n'
#         entry = entry+'Adduct=' +str(row['adduct']) +'\n'
        entry = entry+(row['ms1'])+'\n'
        entry = entry + 'END IONS'+'\n'
        entry = entry +'\n'
    #     ms2
        entry = entry + 'BEGIN IONS'+'\n'
        entry = entry + 'PEPMASS='+str(pep_mass)+'\n'
        entry = entry + 'MSLEVEL=2'+'\n'
        entry = entry+'CHARGE=' + charge +'\n'
        entry = entry+(row['msms'])+'\n'
        entry = entry + 'END IONS'+'\n'
    text_file = open(output+'.mgf', "w",encoding='utf-8')
    text_file.write(entry)
    text_file.close()
        # break

def export_mat(data,output_location, typeofms1='ms1',typeofmsms = 'msms', ifcollision_energy = True):
    entry = ''
    for index, row in data.iterrows():
        entry = entry + 'NAME: ' + row['NAME'] + '\n'
        entry = entry + 'PrecursorMZ: ' + str(row['PRECURSORMZ']) + '\n'
        entry = entry + 'PRECURSORTYPE: ' + row['Adduct'] + '\n'
        entry = entry + 'INSTRUMENTTYPE: ' +'\n' #empty line
        entry = entry + 'INSTRUMENT: ' +'\n'
        entry = entry + 'Authors: ' +'Arpana, Parker and Fanzhou'+'\n'
        entry = entry + 'License: '+'\n'
        entry = entry + 'FORMULA: ' + str(row['Formula']) + '\n'
        entry = entry + 'ONTOLOGY: ' +'\n'
        entry = entry + 'SMILES: ' +'\n'
        entry = entry + 'INCHIKEY: ' + row['InChIKey'] + '\n'
        entry = entry + 'INCHI: ' +'\n'
        entry = entry+'IONMODE: '+row['Ion_mode']+ '\n'
        if ifcollision_energy:
            entry = entry + 'Collision_enerty: ' + str(row['Collision_energy']) +'eV'+ '\n'
        entry = entry+'SPECTRUMTYPE: Centroid and composite'+ '\n'
        entry = entry + 'METABOLITENAME: ' + '\n'
        entry = entry + 'SCANNUMBER: Alignment ID '+str(row['Alignment_ID']) + '\n'
        entry = entry + 'RETENTIONTIME: ' + str(row['RETENTIONTIME']) + '\n'
        entry = entry + 'RETENTIONINDEX: N/A' +'\n'
        entry = entry + 'CCS: ' +'\n'
        entry = entry + 'INTENSITY: ' +str(row['intensity'])+'\n'
        entry = entry + '#Specific field for labeled experiment' +'\n'
        entry = entry + 'IsMarked: False' +'\n'
        entry = entry + 'Comment: ' + str(row['Comment']) + '\n'
        entry = entry + 'MSTYPE: MS1'+ '\n'
        entry = entry + 'Num Peaks: ' + str(num_peaks(row[typeofms1])) + '\n'
        mass1, intensity1= break_spectra(row[typeofms1])
        for i in range(0, len(mass1)):
            entry = entry + str(mass1[i]) +'\t'+ str(intensity1[i]) + '\t'+'"'+str(mass1[i])+'"'+'\n'
        entry = entry + 'MSTYPE: MS2'+ '\n'
        entry = entry + 'Num Peaks: ' + str(num_peaks(row[typeofmsms])) + '\n'
        mass2, intensity2= break_spectra(row[typeofmsms])
        for i in range(0, len(mass2)):
            entry = entry + str(mass2[i]) +'\t'+ str(intensity2[i]) + '\t'+'"'+str(mass2[i])+'"'+'\n'
        entry = entry + '\n'
        #
        #
        # entry = entry +'Spectrum_type: '+row['Spectrum_type']+ '\n'
        #
        # entry = entry + 'InChIKey: ' + row['InChIKey'] + '\n'
        #
        # entry = entry + 'ExactMass: ' + str(row['ExactMass']) + '\n'
        #
        #
        #
        # entry = entry + 'Num peaks: ' + str(num_peaks(row[typeofmsms])) + '\n'
        # entry = entry + row[typeofmsms]
        # # entry = entry +str(row['count'])

        # entry = entry + '\n'

    #open text file
    text_file = open(output_location, "w",encoding='utf-8')

    #write string to file
    text_file.write(entry)

    #close file
    text_file.close()
# print("i am spectra operation")

