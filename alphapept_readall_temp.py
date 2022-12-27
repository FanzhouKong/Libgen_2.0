import sys
import re
import gc
import pandas as pd
# %load_ext autoreload
# %autoreload 2
import os
import itertools
import gc
import os
from tqdm import tqdm
from mimas.tools.spectral_file.extract_ms1_feature import process_mzml_file
import time
import fnmatch
import os
import sys
# mzml_dir = sys.argv[1]
mzml_folders = ['/Volumes/Samsung_T5/20220713_TS_NP_HCD-UVPD_sequential/neg_UVPD-MS2',
'/Volumes/Samsung_T5/20220713_TS_NP_HCD-UVPD_sequential/pos_UVPD-MS2',
'/Volumes/Samsung_T5/20220713_TS_NP_HCD-UVPD_sequential/pos_HCD-MS2']
# output_folder_path = (sys.argv[2])
# if os.path.exists(output_folder_path) == False:
#     os.mkdir(output_folder_path)
output_folders = ['/Users/fanzhoukong/Documents/GitHub/Libgen_data/UVPD/Features/UVPD_neg_MS2',
'/Users/fanzhoukong/Documents/GitHub/Libgen_data/UVPD/Features/UVPD_pos_MS2',
'/Users/fanzhoukong/Documents/GitHub/Libgen_data/UVPD/Features/HCD_pos_MS2']
import glob
for i in range(0,3):
    os.chdir(mzml_folders[i])
    bad_files = []
    for file in glob.glob("*.mzML"):
        # print(file[0:-5])
        # break
        try:
            ms2_selected = process_mzml_file(os.path.join(mzml_folders[i], file))
            ms2_selected.to_csv((os.path.join(output_folders[i], file[0:-5]+".csv")), index=False)
        except:
            bad_files.append(file)
    bad_files =  pd.DataFrame(bad_files)
    bad_files.to_csv(os.path.join(output_folders[i], "bad_files.csv"), index = False)
