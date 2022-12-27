
import sys
import os
import numpy as np
from tqdm import tqdm
import pandas as pd
import numpy as np
source_dir = sys.argv[1]
filename = sys.argv[2]
for_denoise = pd.read_csv(os.path.join(source_dir, filename))
from toolsets.denoising_related_functions import denoise_blacklist
from toolsets.spectra_operations import calculate_explained_intensity
msms_d = []
ei = []
for index, row in tqdm(for_denoise.iterrows(), total=len(for_denoise)) :
    try:
        msms_d_temp =(denoise_blacklist(row, typeofmsms = 'peaks_recalibrated', mass_error=0.05, ifppm=False))
        msms_d.append(msms_d_temp)
        ei.append(calculate_explained_intensity(msms_d_temp, row['peaks_recalibrated'],row['reference_precursor_mz'] ))
    except:
        msms_d.append(np.NAN)
        ei.append(np.NAN)
for_denoise['peaks_denoised'] = msms_d
for_denoise['ei'] = ei
for_denoise['c_id']=np.arange(len(for_denoise))
from toolsets.denoising_related_functions import post_processing
for_denoise.to_csv(os.path.join(source_dir, filename.split('.')[:-1][0]+'_denoised.csv'), index = False)
for_denoise_post = post_processing(for_denoise)
for_denoise_post.to_csv(os.path.join(source_dir, filename.split('.')[:-1][0]+'_denoised_post.csv'), index = False)