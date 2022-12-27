import time
import os
import seaborn as sns
import numpy as np
import pandas as pd
import sys

from mimas.tools.spectral_file.extract_ms1_feature import process_mzml_file, extract_ms2_spectra, find_features
import logging
from toolsets.spectra_operations import entropy_similarity_default
from toolsets.adduct_calculator import complete_adducts, complete_formula
# from toolsets.API_gets import complete_smiles
from toolsets.search import string_search, num_search
import toolsets.spectra_operations as so
from tqdm import tqdm
from toolsets.features_by_alphapept import find_features_alphapept


#  usage argv[1]: std_list, in csv;
# argv[2]: feature folder (alphapepet)
# argv[3]: output filename
print("all packages imported")
std_list = pd.read_csv(sys.argv[1])
print("std list imported")
from toolsets.precursor_matching import precursor_matching_global
feature_dir = (sys.argv[2])
print("started initial matching...")
matched = precursor_matching_global(std_list, feature_dir, adducts = std_list.columns[-3:].values.tolist(),step=10, ppm = True)
print("initial matching done")
from toolsets.helpers import split_unique_duplicates
print("first round of unique keys confirming...")
matched_confirmed, matched_unconfirmed = split_unique_duplicates(matched)
print("first round of unique keys confirmed")
from toolsets.helpers import match_by_rt
print("second round of confirming via rt...")
matched_confirmed_rt, matched_unconfirmed_rt =match_by_rt(matched_confirmed, matched_unconfirmed)
print("second round of confirming via rt done")
from toolsets.helpers import deduplicate
print("start deduplicating of unconfirmed matches...")
deduplicated = deduplicate(matched_unconfirmed_rt)
print("deduplicating of unconfirmed matches done")
result = pd.concat([matched_confirmed_rt, deduplicated], axis=0)
result.reset_index(inplace=True, drop= True)
from toolsets.mass_recalibration import data_recalibrate
print("recalibrating started...")
data_r = data_recalibrate(result, save_diff=True)
print("recalibrating done")
from toolsets.spectra_operations import denoising
print("denoising started")
data_d = denoising(data_r, typeofmsms='peaks_recalibrated',mass_error = 0.02)
print("denoised done")
explained_intensity = []
unassigned_peak_intensity = []
for index, row in data_d.iterrows():
    explained_intensity.append(so.calculate_explained_intensity(row['peaks_recalibrated'], row['peaks_recalibrated_denoised'], row['reference_precursor_mz']))
    unassigned_peak_intensity.append(so.identify_max_unassigned_intensity(row['peaks_recalibrated'], row['peaks_recalibrated_denoised'], row['reference_precursor_mz']))
data_d['explained_intensity']=explained_intensity
data_d['max_unassigned_peak_intensity']=unassigned_peak_intensity
data_d.to_csv(sys.argv[3], index = False)


