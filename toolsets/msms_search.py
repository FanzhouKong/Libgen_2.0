import toolsets.spectra_operations as so
import pandas as pd
import numpy as np
import multiprocessing as mp
# import yuanyue_code.msp_file as msp
# import spectral_entropy as se
from tqdm import tqdm
import warnings
from toolsets.spectra_operations import entropy_similarity_default
from toolsets.search import num_search
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import os
import sys
def mute():
    sys.stdout = open(os.devnull, 'w')

# def initialize_search_table():
#     isf_table = {}
#     for status in ['Unknown', 'Reference']:
#         for column in ['PRECURSORMZ', 'PRECURSORTYPE','RETENTIONTIME','Comment','peaks']:
#             isf_table[status+'_'+column] = []
#     isf_table['Compound_NAME']=[]
#     isf_table['Entropy_similarity']=[]
#     isf_table['Precursor_found'] =[]
#     return(isf_table)
def entropy_search_data(sample, library):
    search_result = pd.DataFrame()
    for index, row in tqdm(sample.iterrows(), total = len(sample)):
        result_temp = entropy_search_identity(row['peaks_cleaned'], row['PRECURSORMZ'], library)
        if result_temp is not np.NAN:

            result_temp['sample_id']=row['Comment']
            result_temp['sample_peaks']=row['peaks_cleaned']
            result_temp['sample_rt']=row['RETENTIONTIME']
            search_result = search_result.append(result_temp)
            # break
    search_result.reset_index(inplace = True, drop = True)
    return(search_result)
def entropy_search_identity(msms, precursormz,library, return_type = 'all'):
    library_temp = num_search(library, 'reference_precursor_mz',precursormz, direction='between', step = 0.01,inclusion = True)
    library_temp.reset_index(drop = True, inplace= True)
    # library_temp = library[(library['PRECURSORMZ'].between(precursormz-(precursormz*10/1E6), precursormz+(precursormz*10/1E6), inclusive=False))]
    if(len(library_temp)<1):
        return(np.NAN)
    entropy = []
    for index, row in library_temp.iterrows():
        entropy_temp = entropy_similarity_default(msms, row['peaks_denoised_normalized'], pmz1 = precursormz, pmz2 = row['reference_precursor_mz'], need_clean = True)
        # entropy_temp = se.similarity(so.convert_string_to_nist(msms), so.convert_string_to_nist(row['peaks_recalibrated_denoised']), 'entropy',
        #                              ms2_da = tolerence, need_clean_spectra = True, need_normalize_result = True)
        entropy.append(entropy_temp)
    # return(entropy_temp)
    if max(entropy)<0.75:
        return(np.NAN)
    if return_type=='all':
        indices = []

        for index, item in enumerate(entropy):
            if item >= 0.75:
                indices.append(index)
        search_result = pd.DataFrame()
        matched_entropy = []
        for index in indices:
            search_result = search_result.append(library_temp.iloc[index])
            matched_entropy.append(entropy[index])
        search_result['matched_entropy']=matched_entropy

        return(search_result)
    elif return_type == 'max':
        search_result = pd.DataFrame()
        index = np.argmax(entropy)
        search_result = search_result.append(library_temp.iloc[index])
        search_result['matched_entropy']=entropy[index]

        return(search_result)





def exact_lookup(instance, library, typeofmsms='msms', threshold = 0.01):
    entropy = []
    library_subset = library.loc[library['key']==instance['key']]
    for i in range(0, len(library_subset)):
        entropy.append(se.similarity(so.convert_string_to_nist(instance[typeofmsms]), so.convert_string_to_nist(library_subset.iloc[i]['msms']), 'entropy',
                                     ms2_da = threshold, need_clean_spectra = True, need_normalize_result = True))
    index_max = np.argmax(entropy)

    return(entropy[index_max], library_subset.iloc[index_max]['msms'])


print("i am msms_search!!!!!")



