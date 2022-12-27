import sys
import re
import pandas as pd
# %load_ext autoreload
# %autoreload 2
import os
import itertools
import os
from tqdm import tqdm
from mimas.tools.spectral_file.extract_ms1_feature import process_mzml_file
import time
import fnmatch
import os
import sys
mzml_dir = sys.argv[1]
output_folder_path = (sys.argv[2])
if os.path.exists(output_folder_path) == False:
    os.mkdir(output_folder_path)

import glob
os.chdir(mzml_dir)
bad_files = []
for file in glob.glob("*.mzML"):
    # print(file[0:-5])
    # break
    try:
        ms2_selected = process_mzml_file(os.path.join(mzml_dir, file))
        ms2_selected.to_csv((os.path.join(output_folder_path, file[0:-5]+".csv")), index=False)
    except:
        bad_files.append(file)
bad_files =  pd.DataFrame(bad_files)
bad_files.to_csv(os.path.join(output_folder_path, "bad_files.csv"), index = False)