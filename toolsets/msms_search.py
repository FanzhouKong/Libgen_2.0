import toolsets.spectra_operations as so
import pandas as pd
import numpy as np
import multiprocessing as mp
import yuanyue_code.msp_file as msp
import spectral_entropy as se
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import os
import sys
def mute():
    sys.stdout = open(os.devnull, 'w')



def entropy_search_fast(msms, precursormz,library, tolerence= 0.05):
    library_temp = library[(library['PRECURSORMZ'].between(precursormz-(precursormz*10/1E6), precursormz+(precursormz*10/1E6), inclusive=False))]
    entropy = []
    for index, row in library_temp.iterrows():
        entropy_temp = se.similarity(so.convert_string_to_nist(msms), so.convert_string_to_nist(row['msms']), 'entropy',
                                     ms2_da = tolerence, need_clean_spectra = True, need_normalize_result = True)
        entropy.append(entropy_temp)
    index_max = np.argmax(entropy)
    return(entropy[index_max], library_temp.iloc[index_max])





def exact_lookup(instance, library, typeofmsms='msms', threshold = 0.01):
    entropy = []
    library_subset = library.loc[library['key']==instance['key']]
    for i in range(0, len(library_subset)):
        entropy.append(se.similarity(so.convert_string_to_nist(instance[typeofmsms]), so.convert_string_to_nist(library_subset.iloc[i]['msms']), 'entropy',
                                     ms2_da = threshold, need_clean_spectra = True, need_normalize_result = True))
    index_max = np.argmax(entropy)

    return(entropy[index_max], library_subset.iloc[index_max]['msms'])


print("i am msms_search!!!!!")



