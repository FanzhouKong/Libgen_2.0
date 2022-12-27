import pandas as pd
from molmass import Formula
from functools import reduce
def save_value_counts(data, column):
    vc = data[column].value_counts().rename_axis('unique_values').to_frame('counts')
    vc.index.name = 'unique_values'
    vc.reset_index(inplace=True)
    return(vc)

def get_mono_mass(data, formula_column = "Formula"):
    mass = []
    for index, row in data.iterrows():
        mass.append(Formula(row[formula_column]).isotope.mass)
    data.insert(data.columns.get_loc(formula_column)+1, "Monoisotopic mass", mass)
    return(data)
def get_unique_list(list1):
    list_set = set(list1)
    unique_list = (list(list_set))
    return(unique_list)
def find_common_items(columns):
    return (list(reduce(set.intersection, map(set, columns))))




