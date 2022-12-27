#!/usr/bin/env python
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
from mimas.external.alphapept import feature_finding
from mimas.external.alphapept.constants import averagine_aa, isotopes
from mimas.external.features_by_alphapept.load_mzml_data import load_mzml_data
from mimas.external.features_by_alphapept.ms_spectrum import MSSpectrum
from mimas.helper.arguments import Arguments
from mimas.spectra.similarity.tools import clean_spectrum
from mimas.file_io import spec_file


def main(parameters):
    filename_base = parameters.file_input.name
    filename_base = filename_base[:filename_base.rfind(".mzML")]
    process_mzml_file(parameters.file_input, parameters.path_output/filename_base)


def process_mzml_file(file_input, ifSciex = True, ifdebug = False):
    # path_output.mkdir(parents=True, exist_ok=True)
    # file_features_all = path_output / "features_all.csv"
    # file_features_mapped = path_output / "features_mapped.csv"
    # file_output_select_msp = path_output / "ms2_selected.msp"

    # Read mzML file
    try:
        ms_file = load_mzml_data(str(file_input), ifSciex=ifSciex)
    except Exception as e:
        # traceback.print_exc(file=open(path_output / "error_in_reading_mzml.txt", "w"))
        return 0
    if len(ms_file["ms_list_ms2"]) == 0:
        print("there is something going wrong")
        # open(path_output/"no_ms2", "wt").close()
        return 0

    # Process features
    features_all, features_mapped = find_features(ms_file)
    features_mapped.dropna(subset=["feature_idx"], inplace=True)
    features_mapped = features_mapped[abs(features_mapped["mz_offset"]) <= 0.5]#can be changed
    features_mapped["scan_number"] = features_mapped["query_idx"].map(lambda x: ms_file["scan_list_ms2"][x])
    features_mapped["charge"] = features_mapped["charge_matched"]*abs(features_mapped["charge"])/features_mapped["charge"]
    features_mapped.drop(["mass", "rt_matched", "rt_offset", "mass_matched", "mass_offset",
                         "mz_offset", "charge_matched", "charge_offset", "fwhm", "dist"], axis=1, inplace=True)
    if len(features_mapped) == 0:
        print("there is no mapped ms2")
        return 0
    if ifdebug == True:
        return(features_mapped)
    # features_all.to_csv(file_features_all, index=False)
    # features_mapped.to_csv(file_features_mapped, index=False)

    # debug purpose only
    # return(features_mapped)
     # Extract MS/MS spectra
    ms2_selected= extract_ms2_spectra(ms_file, features_mapped)

    # write_spectra_to_file(ms2_selected, file_output_select_msp)
    return ((ms2_selected)
            # , pd.DataFrame.from_dict(ms2_unmapped)
            )


# def write_spectra_to_file(ms2_collection, file_output):
#     with open(file_output, "wt") as fo:
#         for i in range(len(ms2_collection["scan_number"])):
#             spec = {
#                 "scan": ms2_collection["scan_number"][i],
#                 "precursormz": ms2_collection["precursor_mz"][i],
#                 "charge": ms2_collection["charge"][i],
#                 "Precursor_type": {1: "[M+H]+", -1: "[M-H]-", 2: "[M+2H]2+", -2: "[M-2H]2-"}[ms2_collection["charge"][i]],
#                 "ms1_intensity_ratio": ms2_collection["ms1_intensity_ratio"][i],
#                 "retention_time": ms2_collection["retention_time"][i],
#                 "peaks": ms2_collection["peaks"][i]
#             }
#             spec_file.write_one_spectrum(fo, spec, ".msp")


def extract_ms2_spectra(ms_file, features_mapped: pd.DataFrame):
    features_mapped["rt_offset"] = abs(features_mapped["rt"]-features_mapped["rt_apex"])
    features_mapped.sort_values(by=["rt_offset"], inplace=True)
    features_mapped_selected = features_mapped.drop_duplicates(subset=["feature_idx"], keep="first", inplace=False).copy()

    # Calculate MS1 intensity ratios
    ms1_ppm = 20
    df_ms1_idx_next = np.searchsorted(ms_file["rt_list_ms1"], features_mapped_selected["rt"], side='left')
    df_ms1_idx_pre = np.where(df_ms1_idx_next > 0, df_ms1_idx_next-1, df_ms1_idx_next)
    df_ms1_rt_next = ms_file["rt_list_ms1"][df_ms1_idx_next]
    df_ms1_rt_pre = ms_file["rt_list_ms1"][df_ms1_idx_pre]
    df_ms1_idx = np.where(abs(features_mapped_selected["rt"]-df_ms1_rt_next) < abs(features_mapped_selected["rt"]-df_ms1_rt_pre),
                          df_ms1_idx_next, df_ms1_idx_pre)
    features_mapped_selected["ms1_idx"] = df_ms1_idx
    for i, row in features_mapped_selected.iterrows():
        intensity_precursor_ion = _extract_ms1_intensity(
            ms_file, int(row["ms1_idx"]), row["mz"] * (1 - ms1_ppm * 1e-6), row["mz"] * (1 + ms1_ppm * 1e-6))
        intensity_select_window = _extract_ms1_intensity(
            ms_file, int(row["ms1_idx"]),  *(ms_file["select_windows_ms2"][int(row["query_idx"])])
        )
        if intensity_select_window == 0 or np.isnan(intensity_select_window) or np.isnan(intensity_precursor_ion):
            features_mapped_selected.loc[i, "ms1_intensity_ratio"] = 0
            features_mapped_selected.loc[i, "precursor_intensity"] = 0#i added this
        else:
            # features_mapped_selected.loc[i, "ms1_intensity_ratio"] = intensity_precursor_ion/intensity_select_window
            features_mapped_selected.loc[i, "ms1_intensity_ratio"] = 1
            features_mapped_selected.loc[i, "precursor_intensity"] = intensity_precursor_ion#i added this

    # Generate selected MS2 spectra
    mapped_ms2_idx = set(features_mapped["query_idx"].tolist())
    # selected_ms2_idx = features_mapped_selected.loc[features_mapped_selected["ms1_intensity_ratio"] >= 0.5, "query_idx"].tolist()
    selected_ms2_idx = features_mapped_selected.loc[:, "query_idx"].tolist()
    selected_ms2_idx.sort()
    idx_to_ms1_intensity_ratio = {int(row["query_idx"]): row["ms1_intensity_ratio"] 
                                  for i, row in features_mapped_selected.iterrows()}
    idx_to_ms1_precursor_intensity =  {int(row["query_idx"]): row["precursor_intensity"] 
                                  for i, row in features_mapped_selected.iterrows()}                             
    features_mapped_selected.sort_values(by=['query_idx'], ascending=True, inplace=True, ignore_index=True)
    selected_ms2 = {
        "scan_number": [ms_file["scan_list_ms2"][i] for i in selected_ms2_idx],
        "precursor_mz": [ms_file["mono_mzs2"][i] for i in selected_ms2_idx],
        "charge": [ms_file["charge2"][i] for i in selected_ms2_idx],
        "ms1_intensity_ratio": [idx_to_ms1_intensity_ratio[i] for i in selected_ms2_idx],
        "ms1_precursor_intensity":[idx_to_ms1_precursor_intensity[i] for i in selected_ms2_idx],
        "retention_time": [ms_file["rt_list_ms2"][i] for i in selected_ms2_idx],
        "peaks": [np.array([ms_file["mass_list_ms2"][i], ms_file["int_list_ms2"][i]]).T for i in selected_ms2_idx],
        # "ms1_index":features_mapped_selected['ms1_idx'],
        # "query_idx":features_mapped_selected['query_idx'],
    }

    unmapped_ms2 = {
        "scan_number": [ms_file["scan_list_ms2"][i] for i in range(len(ms_file["scan_list_ms2"])) if i not in mapped_ms2_idx],
        "precursor_mz": [ms_file["mono_mzs2"][i] for i in range(len(ms_file["mono_mzs2"])) if i not in mapped_ms2_idx],
        "charge": [ms_file["charge2"][i] for i in range(len(ms_file["charge2"])) if i not in mapped_ms2_idx],
        "peaks": [np.array([ms_file["mass_list_ms2"][i], ms_file["int_list_ms2"][i]]).T
                  for i in range(len(ms_file["mass_list_ms2"])) if i not in mapped_ms2_idx],
    }

    # Clean the MS/MS spectra
    # def clean_ms2(ms2_info):
    #     spec_all = ms2_info["peaks"]
    #     for i, (spec, precursor_mz, charge) in enumerate(zip(spec_all, ms2_info["precursor_mz"], ms2_info["charge"])):
    #         # spec_clean = clean_spectrum(
    #             # spec, max_mz=precursor_mz-1.5/abs(charge), noise_threshold=0.01, ms2_da=0.05)#this is original one
    #         spec_clean = 
    #         spec_all[i] = spec_clean
    #     ms2_info["peaks"] = spec_all
    #     for key, value in ms2_info.items():
    #         ms2_info[key] = [value[i] for i in range(len(spec_all)) if len(spec_all[i]) > 0]
    #     return ms2_info
    # selected_ms2_cleaned = clean_ms2(selected_ms2)
    # unmapped_ms2_cleaned = clean_ms2(unmapped_ms2)
    selected_ms2_cleaned = pd.DataFrame.from_dict(selected_ms2)
    unmapped_ms2_cleaned = pd.DataFrame.from_dict(unmapped_ms2)
    selected_ms2_cleaned['ms1_index']=features_mapped_selected['ms1_idx']
    selected_ms2_cleaned['query_idx']=features_mapped_selected['query_idx']

    return selected_ms2_cleaned


def _extract_ms1_intensity(ms_file, ms1_idx, mz_lower, mz_upper):
    # Find MS1 scan
    idx_start = ms_file["indices_ms1"][ms1_idx]
    idx_end = ms_file["indices_ms1"][ms1_idx+1]
    mz = ms_file["mass_list_ms1"][idx_start:idx_end]
    peaks_select = np.bitwise_and(mz >= mz_lower, mz <= mz_upper)
    return np.sum(ms_file["int_list_ms1"][idx_start:idx_end][peaks_select])


def find_features(query_data):
    f_settings = {'centroid_tol': 8,
                  'hill_check_large': 40,
                  'hill_length_min': 3,
                  'hill_nboot': 150,
                  'hill_nboot_max': 300,
                  # 'hill_smoothing': 1,
                  'hill_smoothing':5,
                  'hill_split_level': 1.3,
                  'iso_charge_max': 3,
                  'iso_charge_min': 1,
                  'iso_corr_min': 0.6,
                  'iso_mass_range': 5,
                  'iso_n_seeds': 100,
                  'iso_split_level': 1.3,
                  'map_mob_range': 0.3,
                  # 'map_mz_range': 0.5,
                  'map_mz_range': 1,
                  'map_n_neighbors': 5,
                  'map_rt_range': 0.5,
                  'max_gap': 2,
                  'search_unidentified': True}
    # print(f_settings['map_mz_range'])
    max_gap = f_settings['max_gap']
    centroid_tol = f_settings['centroid_tol']
    hill_split_level = f_settings['hill_split_level']
    iso_split_level = f_settings['iso_split_level']
    window = f_settings['hill_smoothing']
    hill_check_large = f_settings['hill_check_large']
    iso_charge_min = f_settings['iso_charge_min']
    iso_charge_max = f_settings['iso_charge_max']
    iso_n_seeds = f_settings['iso_n_seeds']
    hill_nboot_max = f_settings['hill_nboot_max']
    hill_nboot = f_settings['hill_nboot']
    iso_mass_range = f_settings['iso_mass_range']
    iso_corr_min = f_settings['iso_corr_min']

    query_data["int_list_ms1"] = np.array(query_data["int_list_ms1"]).astype(int)
    int_data = query_data['int_list_ms1']
    # logging.info('Feature finding on {}'.format(file_name))
    # logging.info(f'Hill extraction with centroid_tol {centroid_tol} and max_gap {max_gap}')

    hill_ptrs, hill_data, path_node_cnt, score_median, score_std = feature_finding.extract_hills(
        query_data, max_gap, centroid_tol)

    # logging.info(f'Number of hills {len(hill_ptrs):,}, len = {np.mean(path_node_cnt):.2f}')
    # logging.info(f'Repeating hill extraction with centroid_tol {score_median+score_std*3:.2f}')

    hill_ptrs, hill_data, path_node_cnt, score_median, score_std = feature_finding.extract_hills(
        query_data, max_gap, score_median + score_std * 3)

    # logging.info(f'Number of hills {len(hill_ptrs):,}, len = {np.mean(path_node_cnt):.2f}')

    hill_ptrs, hill_data = feature_finding.remove_duplicate_hills(hill_ptrs, hill_data, path_node_cnt)
    # logging.info(f'After duplicate removal of hills {len(hill_ptrs):,}')

    hill_ptrs = feature_finding.split_hills(hill_ptrs, hill_data, int_data, hill_split_level=hill_split_level,
                                            window=window)  # hill lenght is inthere already
    # logging.info(f'After split hill_ptrs {len(hill_ptrs):,}')

    hill_data, hill_ptrs = feature_finding.filter_hills(
        hill_data, hill_ptrs, int_data, hill_check_large=hill_check_large, window=window)

    # logging.info(f'After filter hill_ptrs {len(hill_ptrs):,}')

    stats, sortindex_, idxs_upper, scan_idx, hill_data, hill_ptrs = feature_finding.get_hill_data(
        query_data, hill_ptrs, hill_data, hill_nboot_max=hill_nboot_max, hill_nboot=hill_nboot)
    # logging.info('Extracting hill stats complete')

    pre_isotope_patterns = feature_finding.get_pre_isotope_patterns(
        stats, idxs_upper, sortindex_, hill_ptrs, hill_data, int_data, scan_idx, feature_finding.maximum_offset,
        iso_charge_min=iso_charge_min, iso_charge_max=iso_charge_max, iso_mass_range=iso_mass_range, cc_cutoff=iso_corr_min)
    # logging.info('Found {:,} pre isotope patterns.'.format(len(pre_isotope_patterns)))

    isotope_patterns, iso_idx, isotope_charges = feature_finding.get_isotope_patterns(
        pre_isotope_patterns, hill_ptrs, hill_data, int_data, scan_idx, stats, sortindex_, averagine_aa, isotopes,
        iso_charge_min=iso_charge_min, iso_charge_max=iso_charge_max, iso_mass_range=iso_mass_range,
        iso_n_seeds=iso_n_seeds, cc_cutoff=iso_corr_min, iso_split_level=iso_split_level, callback=None)
    # logging.info('Extracted {:,} isotope patterns.'.format(len(isotope_charges)))

    feature_table, lookup_idx = feature_finding.feature_finder_report(
        query_data, isotope_patterns, isotope_charges, iso_idx, stats, sortindex_, hill_ptrs, hill_data)
    # logging.info('Report complete.')

    # Calculate additional params
    feature_table['rt_length'] = feature_table['rt_end'] - feature_table['rt_start']
    feature_table['rt_right'] = feature_table['rt_end'] - feature_table['rt_apex']
    feature_table['rt_left'] = feature_table['rt_apex'] - feature_table['rt_start']
    feature_table['rt_tail'] = feature_table['rt_right'] / feature_table['rt_left']

    # logging.info('Matching features to query data.')

    if 'mono_mzs2' not in query_data.keys():
        # logging.info('No MS2-data to match.')
        features = pd.DataFrame()
    else:
        features = feature_finding.map_ms2(feature_table, query_data, **f_settings)

    return feature_table, features


if __name__ == '__main__':
    args = Arguments()
    para = {
        "file_input": r"data/from_shen-2022_03_03/mzml/NIH_Lip_Std_CSH_NEG_TissuePool_01.mzML.gz",
        "path_output": r"result/2022_03/0312_alphapept_vs_msdial/data/alphapept"
    }

    args.add_argument_from_dictionary(para)
    para = args.parse_args(print_parameter=False)
    main(para)
