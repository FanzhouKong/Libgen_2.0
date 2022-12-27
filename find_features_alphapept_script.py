import os
import sys
from toolsets.features_by_alphapept import find_features_alphapept
import pandas as pd
std_list = pd.read_csv("/Users/fanzhoukong/Documents/GitHub/Libgen_data/EAD/standard_list_ori_values.csv")
mzml_folder = '/Users/fanzhoukong/Documents/GitHub/Libgen_data/EAD/mzml'
find_features_alphapept(std_list,mzml_folder, outputfoldername=sys.argv[1] )