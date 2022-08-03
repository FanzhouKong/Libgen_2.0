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
import toolsets.spectra_operations as so
import numpy as np
from molmass import Formula

import toolsets.mass_to_formula as mtf
def prep_msms(msms_nist, ifnist = False):
    if ifnist ==True:
        msms = so.convert_nist_to_string(msms_nist)
    else:
        msms = msms_nist
    msms_sorted = so.sort_spectra(msms)

    return(so.break_spectra(msms_sorted))




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
def get_parentpeak(mass, intensity, parention, mass_error = 0.01, ifppm = False):
    parention = float(parention)
    tolerance = so.set_tolerance(mass_error, ifppm, parention)
    lower_bound = parention-tolerance
    upper_bound = parention+tolerance

    lower_bound_i = bisect.bisect_left(mass, lower_bound)
    upper_bound_i = bisect.bisect_right(mass, upper_bound, lo=lower_bound_i)

    mass_candi = mass[lower_bound_i:upper_bound_i]
    # print("checkpoint!")
    intensity_candi = intensity[lower_bound_i:upper_bound_i]
    # return(mass_candi)
    if mass_candi !=[]:
        max_index = intensity_candi.index(max(intensity_candi))
        return(mass_candi[max_index], intensity_candi[max_index])
    else:
        return(parention, 0)

def truncate_msms(mass, intensity, parention):
    upper_allowed = bisect.bisect_right(mass, parention)
    # pep_index = mass.index(parention)
    mass = mass[0:upper_allowed]
    intensity = intensity[0:upper_allowed]
    return(mass, intensity)

def evaluate_nl_blacklist(nl):
    if len(Formula(nl).composition())==1:
        if Formula(nl).composition()[0][0] =='C' or Formula(nl).composition()[0][0] =='N':
            return(False)
        else:
            return(True)
    else:
        return(True)

def get_unique_list(list1):
    list_set = set(list1)
    unique_list = (list(list_set))
    return(unique_list)

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

def denoise_blacklist(instance, typeofmsms = 'msms', mass_error =0.01, ifppm = False):
    mass, intensity = prep_msms(instance[typeofmsms], ifnist=False)
    parention = float(instance['parent_ion'])
    pep_mass, pep_intensity = get_parentpeak(mass, intensity, parention)
    mass, intensity = truncate_msms(mass, intensity, pep_mass)
    formula= instance['Formula']
    if instance['Adduct']=='[M+Na]+':
        formula = formula+'Na'
    if instance['Adduct']=='[M+NH4+]':
        formula = formula+'NH3'
    if instance['Adduct']=='[M+Cl]-':
        formula = formula+'Cl'
    mass_d = []
    intensity_d = []
    threshold = so.set_tolerance(mass_error = mass_error, ifppm=ifppm, precursormz=instance['parent_ion'])

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
    if pep_intensity != 0:
        mass_d.append(pep_mass)
        intensity_d.append(pep_intensity)
    return(so.pack_spectra(mass_d, intensity_d))


def denoise_blacklist_msms(msms, parent_ion, adduct,formula, mass_error = 0.01, ifppm = False):
    mass, intensity = prep_msms(msms, ifnist=False)
    pep_mass, pep_intensity = get_parentpeak(mass, intensity, parent_ion)
    mass, intensity = truncate_msms(mass, intensity, pep_mass)
    if adduct=='[M+Na]+':
        formula = formula+'Na'
    if adduct=='[M+NH4+]':
        formula = formula+'NH3'
    if adduct=='[M+Cl]-':
        formula = formula+'Cl'
    mass_d = []
    intensity_d = []
    threshold = so.set_tolerance(mass_error = mass_error, ifppm=ifppm, precursormz=parent_ion)

    # print(threshold)
    for i in range(0,len(mass)):
        # print("i am evaluating", mass[i])
        nl_candidates = mtf.nl_to_formula(parent_ion - mass[i], threshold, formula)
        if nl_candidates != [] and nl_candidates!=['']:
            for nl in nl_candidates:
                if evaluate_nl_blacklist(nl)==True:
                    mass_d.append(mass[i])
                    intensity_d.append(intensity[i])
                    break
    if pep_intensity != 0:
        mass_d.append(pep_mass)
        intensity_d.append(pep_intensity)
    return(so.pack_spectra(mass_d, intensity_d))
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