import pandas as pd
import toolsets.conversion_function as cf
from tqdm import tqdm
def pick_results(data_binbase ):
    score = []
    notes = []
    instrument_type = []
    comments = []
    splash = []
    precursormz = []
    normalized_entropy = []
    entropy = []
    library_inchi = []
    library_adduct = []
    retention_time = []
    # msms = []
    for key in tqdm(data_binbase['query_splash'].unique()):
        data_temp = data_binbase.loc[data_binbase['query_splash']==key]
        if len(data_temp['library_inchikey'].str[0:14].unique())==1:
            score.append(data_temp['score'].max())
            notes.append(data_temp.iloc[0]['query_notes'])
            instrument_type.append(data_temp.iloc[0]['query_instrument_type'])
            comments.append(data_temp.iloc[0]['query_comment'])
            splash.append(data_temp.iloc[0]['query_splash'])
            precursormz.append(data_temp.iloc[0]['query_precursor_mz'])
            retention_time.append(data_temp.iloc[0]['query_rt'])
            normalized_entropy.append(data_temp.iloc[0]['query_entropy_normalized'])
            library_inchi.append(data_temp.iloc[0]['library_inchikey'])
            library_adduct.append(data_temp.iloc[0]['library_adduct'])
            entropy.append(data_temp.iloc[0]['query_entropy'])
            # msms.append(data_temp.iloc[0]['query_'])


    result = pd.DataFrame()
    result['score']=score
    result['notes']=notes
    result['instrument_type']=instrument_type
    result['comments']=comments
    result['splash']=splash
    result['precursormz']=precursormz
    result['normalized_entropy']=normalized_entropy
    result['library_inchi']=library_inchi
    result['library_adduct']=library_adduct
    result['retention_time']=retention_time
    result['entropy']=entropy
    result_filtered = result.loc[result['entropy']>=0.5]

    smiles = []
    for index, row in tqdm(result_filtered.iterrows(), total = result_filtered.shape[0]):
        smiles.append(cf.get_something_from_pubchem("inchikey", row['library_inchi'], "CanonicalSMILES"))

    result_filtered['SMILES']=smiles
    result_filtered=result_filtered[~result_filtered['SMILES'].isnull()]
    return(result_filtered)