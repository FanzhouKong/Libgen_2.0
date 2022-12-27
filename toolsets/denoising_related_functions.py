import os
import sys

# from multiprocess import Pool
import pandas as pd
# from multiprocessing import Pool
from ast import literal_eval
import bisect
import itertools
import time
import toolsets.mass_to_formula as mtf
import numpy as np
from toolsets.search import num_search, string_search
import numpy as np
from molmass import Formula

import toolsets.mass_to_formula as mtf

def denoise_blacklist(instance, typeofmsms = 'peaks_recalibrated', mass_error =0.02, ifppm = False):
    # print("i am in new")
    import toolsets.spectra_operations as so
    mass, intensity = so.break_spectra(instance[typeofmsms])
    parention = float(instance['reference_precursor_mz'])
    # pep_mass, pep_intensity = get_parentpeak(mass, intensity, parention)

    # mass_frag, intensity_frag, mass_precursor, intensity_precursor = remove_precursor(mass, intensity, parention)
    
    adduct = instance['reference_adduct']
    formula = prep_formula(instance['reference_smiles'], adduct)
    # if instance['reference_adduct']=='[M-H2O-H]-' or instance['reference_adduct']=='[M+H2O+H]+':
    #     formula = formula+'H2O'
    mass_d = []
    intensity_d = []
    threshold = so.set_tolerance(mass_error = mass_error, ifppm=ifppm, precursormz=parention)

    # print(threshold)
    for i in range(0,len(mass)):
        # print("i am evaluating", mass[i])
        nl_candidates = mtf.nl_to_formula(parention - mass[i], threshold, formula)
        if nl_candidates != [] and nl_candidates!=['']:
            for nl in nl_candidates:
                if evaluate_nl_blacklist(nl)==True:
                    mass_d.append(mass[i])
                    intensity_d.append(intensity[i])
                    break
    # if intensity_precursor[0] != -1:
    #     mass_d.extend(mass_precursor)
    #     intensity_d.extend(intensity_precursor)
    return(so.pack_spectra(mass_d, intensity_d))
def denoise_single(frag_ion, parention, smiles, mass_error = 0.02, adduct = '[M+H]+',ifppm = False, ifnl = True):
    import toolsets.spectra_operations as so
    parention = float(parention)
    threshold = so.set_tolerance(mass_error = mass_error, ifppm=ifppm, precursormz=parention)
    # formula= str(formula)
    formula = prep_formula(smiles, adduct)
    if ifnl ==True:
        nl_candidates = mtf.nl_to_formula(parention - frag_ion, threshold, formula)
        return(nl_candidates)
    if ifnl == False:
        nl_candidates = mtf.nl_to_formula(frag_ion, threshold, formula)
        return(nl_candidates)

def prep_formula(smiles, adduct = ['[M+H]+']):
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    benzene_pattern = Chem.MolFromSmiles('C1=CC=CC=C1')
    if mol.HasSubstructMatch(benzene_pattern) ==True:
        formula = CalcMolFormula(mol)
        formula = formula +'N2'
    else:
        formula = CalcMolFormula(mol)
    # formula = CalcMolFormula(mol)

    # formula= str(formula)
    if adduct=='[M+Na]+':
        formula = formula+'Na'
    if adduct=='[M+NH4+]':
        formula = formula+'NH3'
    if adduct=='[M+Cl]-':
        formula = formula+'Cl'
    if adduct=='[M+C2H4O2-H]-':
        formula = formula+'C2H4O2'
    if formula[-1]=='+' or formula[-1]=='-':
        formula = formula[:-1]
    return(formula)
def post_processing(data, ei_threshold = 90, denoised_column = 'peaks_denoised', high_quality = True):
    import toolsets.spectra_operations as so
    peaks_denoisied_normalized = []
    normalized_entropy = []
    spectrum_entropy = []
    for index, row in data.iterrows():
        peaks_denoisied_normalized.append(so.normalize_spectrum(row[denoised_column]))
        normalized_entropy.append(so.normalized_entropy(row[denoised_column], order=4))
        spectrum_entropy.append(so.spectral_entropy(row[denoised_column]))
    data['peaks_denoised_normalized']=peaks_denoisied_normalized
    data['spectrum_entropy']=spectrum_entropy
    data['normalized_entropy']=normalized_entropy
    data_good=num_search(data, 'ei', ei_threshold, direction=">", inclusion=True)
        # data=num_search(data, 'ei', ei_threshold, direction="<", inclusion=False)
    data_good =num_search(data_good, 'spectrum_entropy', 0.5, direction='>', inclusion=True)
    data_good.reset_index(drop = True, inplace= True)
    if high_quality ==True:
        return(data_good)
    else:
        data_bad = data[~data['c_id'].isin(data_good['c_id'])]
        data_bad =num_search(data_bad, 'spectrum_entropy', 0.5, direction='>', inclusion=True)
        data_bad.reset_index(drop = True, inplace= True)

        return(data_bad)
# def prep_formula(formula, adduct = ['[M+H]+']):
#     formula= str(formula)
#     if adduct=='[M+Na]+':
#         formula = formula+'Na'
#     if adduct=='[M+NH4+]':
#         formula = formula+'NH3'
#     if adduct=='[M+Cl]-':
#         formula = formula+'Cl'
#     if adduct=='[M+C2H4O2-H]-':
#         formula = formula+'C2H4O2'
#     if formula[-1]=='+' or formula[-1]=='-':
#         formula = formula[:-1]
#     return(formula)
# def find_parention(mass, intensity, precursormz, mass_error = 10, ifppm = True):
#     precursormz = float(precursormz)
#     tol = so.set_tolerance(mass_error, ifppm, precursormz)
#     # mass_raw, intensity_raw = so.break_spectra(msms)
#     # if precursormz <=200:
#     #     tolerance = 0.01
#     # if precursormz >200:
#     #     tolerance = precursormz*tolerance/1E6
#
#     # mass_sorted, intensity_sorted = zip(*sorted(zip(mass_raw, intensity_raw)))
#     lower_bound = precursormz-tol
#     upper_bound = precursormz+tol
#     # print(" i am here!")
#     lower_bound_i = bisect.bisect_left(mass, lower_bound)
#     # print("i have passed lower")
#     upper_bound_i = bisect.bisect_right(mass, upper_bound, lo=lower_bound_i)
#     mass_candi = mass[lower_bound_i:upper_bound_i]
#     intensity_candi = intensity[lower_bound_i:upper_bound_i]
#     # return(mass_candi)
#     if mass_candi ==[]:
#         return(precursormz)
#     else:
#         max_index = intensity_candi.index(max(intensity_candi))
#         return(mass_candi[max_index])

    # return (mass_candi, max_index)
# def get_parentpeak(mass, intensity, parention, mass_error = 0.01, ifppm = False):
#     parention = float(parention)
#     tolerance = so.set_tolerance(mass_error, ifppm, parention)
#     lower_bound = parention-tolerance
#     upper_bound = parention+tolerance

#     lower_bound_i = bisect.bisect_left(mass, lower_bound)
#     upper_bound_i = bisect.bisect_right(mass, upper_bound, lo=lower_bound_i)

#     mass_candi = mass[lower_bound_i:upper_bound_i]
#     # print("checkpoint!")
#     intensity_candi = intensity[lower_bound_i:upper_bound_i]
#     # return(mass_candi)
#     if mass_candi !=[]:
#         max_index = intensity_candi.index(max(intensity_candi))
#         return(mass_candi[max_index], intensity_candi[max_index])
#     else:
#         return(parention, 0)



def evaluate_nl_blacklist(nl):
    if len(Formula(nl).composition())==1:
        if Formula(nl).composition()[0][0] =='C' or Formula(nl).composition()[0][0] =='N':
            return(False)
        else:
            return(True)
    else:
        return(True)



# def evaluate_frag(instance):
#     allowed_formula = []
#     precursormz = float(spec['PrecursorMZ'])
#     mass, intensity = prep_msms(spec['spectrum'], ifnist=True)
#     parention = find_parention(mass, intensity, precursormz)
#     mass, intensity = truncate_msms(mass, intensity, parention)
#     for i in range(0, len(mass)):
#         allowed_formula_temp = mtf.nl_to_formula(parention - mass[i], 10, spec['Formula'])
#         if allowed_formula_temp!=[] and allowed_formula_temp!=['']:
#             allowed_formula.extend(allowed_formula_temp)
#             # return(allowed_formula)
#         # else:
#         #     allowed_formula.extend([-1])
#             # return(allowed_formula)
#     return (allowed_formula)
def find_losses_nist(spec):
    allowed_formula = []
    precursormz = float(spec['PrecursorMZ'])
    mass, intensity = prep_msms(spec['spectrum'], ifnist=True)
    parention = get_parentpeak(mass, intensity, precursormz)
    mass, intensity = truncate_msms(mass, intensity, parention)
    for i in range(0, len(mass)):
        allowed_formula_temp = mtf.nl_to_formula(parention - mass[i], 10, spec['Formula'])
        if allowed_formula_temp!=[] and allowed_formula_temp!=['']:
            allowed_formula.extend(allowed_formula_temp)
            # return(allowed_formula)
        # else:
        #     allowed_formula.extend([-1])
            # return(allowed_formula)
    return (allowed_formula)





    # pep_mass, pep_intensity = get_parentpeak(mass, intensity, precursormz)
    # mass, intensity = truncate_msms(mass, intensity, pep_mass)
    # if adduct=='[M+Na]+':
    #     formula = formula+'Na'
    # if adduct=='[M+NH4+]':
    #     formula = formula+'NH3'
    # if adduct=='[M+Cl]-':
    #     formula = formula+'Cl'
    # mass_d = []
    # intensity_d = []
    # threshold = so.set_tolerance(mass_error, ifppm = ifppm, precursormz=precursormz)
    #
    # # print(threshold)
    # for i in range(0,len(mass)):
    #     # print("i am evaluating", mass[i])
    #     nl_candidates = mtf.nl_to_formula(pep_mass - mass[i], threshold, formula)
    #     if nl_candidates != [] and nl_candidates!=['']:
    #         for nl in nl_candidates:
    #             if evaluate_nl_blacklist(nl)==True:
    #                 mass_d.append(mass[i])
    #                 intensity_d.append(intensity[i])
    #                 break
    # if pep_intensity != 0:
    #     mass_d.append(pep_mass)
    #     intensity_d.append(pep_intensity)
    # return(so.sort_spectra(so.pack_spectra(mass_d, intensity_d)) )