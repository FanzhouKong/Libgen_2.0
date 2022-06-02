from molmass import Formula
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
def calculate_precursormz(formula, adduct):
    proton = 1.00727646677
    Na_plus = 22.989218
    NH4_plus = 18.033823
    HacH_minus = 59.013851
    H2OH_minus = -19.01839
    FaH_minus = 44.998201
    Cl_minus = 34.969402
    if(adduct=='[M+NH4]+'):
        pmz = Formula(formula).isotope.mass+NH4_plus
    if(adduct=='[M+H]+'):
        pmz = Formula(formula).isotope.mass+proton
    if(adduct=='[M+Na]+'):
        pmz = Formula(formula).isotope.mass+Na_plus
    if(adduct=='[M-H]-'):
        pmz = Formula(formula).isotope.mass-proton
    if(adduct=='[M+C2H4O2-H]-'):
        pmz = Formula(formula).isotope.mass+HacH_minus
    if(adduct=='[M-H2O-H]-'):
        pmz = Formula(formula).isotope.mass+H2OH_minus
    if(adduct=='[M+FA-H]-'):
        pmz = Formula(formula).isotope.mass+FaH_minus
    if(adduct=='[M+Cl]-'):
        pmz = Formula(formula).isotope.mass+Cl_minus
    try:
        pmz
    except NameError:
        print("wrong adduct type is passed!")
        print(adduct)
        pmz = np.nan
    return(round(pmz,4))
def nl_calc(formula):
    pmz = Formula(formula).isotope.mass
    return(round(pmz,4))