import pandas as pd
def check_for_missing_compounds(data, std_list):
    missing_compound = list(set(std_list['InChIKey']) - set(data['InChIKey']))
    missing_compounds_list = std_list.loc[std_list['InChIKey'].isin(missing_compound)]
    return(missing_compounds_list)