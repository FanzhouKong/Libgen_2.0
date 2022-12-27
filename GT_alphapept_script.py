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
def find_files(base, pattern):
    '''Return list of files matching pattern in base folder.'''
    return [filename for filename in os.listdir(base) if re.search(pattern, filename, re.IGNORECASE)]

mzml_dir = sys.argv[1]
source_dir = "/Users/fanzhoukong/Documents/GitHub/Libgen_data/GT"
std_list = pd.read_csv(os.path.join(source_dir, "std_list.csv"))
std_list = std_list.applymap(lambda x: x.strip() if isinstance(x, str) else x)
enzyme_list = pd.read_csv(os.path.join(source_dir, "enzyme_list.csv"))
enzyme_list = enzyme_list.applymap(lambda x: x.strip() if isinstance(x, str) else x)
enzyme_list['List of enzymes']=enzyme_list['List of enzymes'].replace('_', '', regex=True)

for index, row in enzyme_list.iterrows():
    if row['List of enzymes'].startswith(str(1))==True:
        row['List of enzymes']=(row['List of enzymes'][1:])

output_folder_path = os.path.join(source_dir, "features_additional")
if os.path.exists(output_folder_path) == False:
        os.mkdir(output_folder_path)

missing_files = []
duplicated_files = []
susceptible_files = []
for mix in std_list['Mix'].unique():
    print("i am processing mix %s, out of %s" %(str(mix), str(len(std_list['Mix'].unique()))))
    print("progress of each enzyme in mix %s" %(str(mix)) )
    for enzyme in tqdm(enzyme_list['List of enzymes'].unique()):

        filename =enzyme+"_"+"Mix"+"_"+str(mix)
        found_file = find_files(mzml_dir, f"^{filename}.mzML")
        if len(found_file)==0:
            found_file = find_files(mzml_dir, f"^{filename}(.*).mzML$")

        if len(found_file)>1:
            duplicated_files.append(filename)
            # print(filename)
        if len(found_file)==0:
            # print("the file %s doesnt exist" %filename)
            missing_files.append(filename)
        if len(found_file)==1:
            try:
                ms2_selected = process_mzml_file(os.path.join(mzml_dir, found_file[0]))
                ms2_selected.to_csv((os.path.join(output_folder_path, filename+".csv")), index=False)
            except:
                susceptible_files.append(filename)
            gc.collect()
print("there are %s of files duplicated" % (str(len(duplicated_files))))
print("there are %s of files missing" % (str(len(missing_files))))
susceptible_files =  pd.DataFrame(susceptible_files)
susceptible_files.to_csv(os.path.join(source_dir, "susceptible_files"+".csv" ),index = False)
missing_files = pd.DataFrame(missing_files)
missing_files.to_csv(os.path.join(source_dir, "missing_files"+".csv"), index = False)