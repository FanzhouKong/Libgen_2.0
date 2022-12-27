#!/usr/bin/python3
import pandas as pd
import logging
import os
import numpy as np
import h5py

from mimas.file_io import spec_file
from mimas.helper.arguments import Arguments
from alphapept.pyrawfilereader import RawFileReader
import alphapept.io
from alphapept.feature_finding import extract_hills, split_hills, filter_hills, get_hill_data, \
    get_pre_isotope_patterns, maximum_offset, get_isotope_patterns, feature_finder_report, extract_bruker, \
    convert_bruker, map_ms2

from alphapept.constants import averagine_aa, isotopes
from mimas.file_io.hdf5_file import AlphaPept

"""
This script will extract features from .raw file
"""


def main(para):
    alpha = AlphaPept()
    alpha.process_file_from_raw(file_raw=para["file_input"])
    features = alpha.get_feature_ms2()
    feature_table = alpha.get_feature_ms1()
    spectra = alpha.get_spectra()
    return


def find_features(file_input, file_hdf5=None, settings: dict = None):
    if file_hdf5 is None:
        file_hdf5 = os.path.splitext(file_input)[0] + ".ms_data.hdf"
    f_settings = {
        "max_gap": 2,
        "centroid_tol": 8,
        "hill_length_min": 3,
        "hill_split_level": 1.3,
        "iso_split_level": 1.3,
        "hill_smoothing": 1,
        "hill_check_large": 40,
        "iso_charge_min": 1,
        "iso_charge_max": 6,
        "iso_n_seeds": 100,
        "hill_nboot_max": 300,
        "hill_nboot": 150,
        "iso_mass_range": 5,
        "iso_corr_min": 0.6,
        "map_mz_range": 0.1,
        "map_rt_range": 0.5,
        "map_mob_range": 0.3,
        "map_n_neighbors": 5,
        "search_unidentified": True
    }
    if settings is not None:
        f_settings.update(settings)

    base, ext = os.path.splitext(file_input)

    if ext.lower() == '.raw':
        datatype = 'thermo'
    elif ext.lower() == '.d':
        datatype = 'bruker'
    elif ext.lower() == '.mzml':
        datatype = 'mzml'
    else:
        raise NotImplementedError('File extension {} not understood.'.format(ext))

    ms_file = alphapept.io.MS_Data_File(file_hdf5, is_new_file=True)
    ms_file.import_raw_DDA_data(file_input, n_most_abundant=-1)
    query_data = ms_file.read_DDA_query_data()

    # if 1:
    if datatype in ['thermo', 'mzml']:

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

        logging.info('Feature finding on {}'.format(file_input))

        logging.info(f'Hill extraction with centroid_tol {centroid_tol} and max_gap {max_gap}')

        hill_ptrs, hill_data, path_node_cnt, score_median, score_std = extract_hills(query_data, max_gap,
                                                                                     centroid_tol)
        logging.info(f'Number of hills {len(hill_ptrs):,}, len = {np.mean(path_node_cnt):.2f}')

        logging.info(f'Repeating hill extraction with centroid_tol {score_median + score_std * 3:.2f}')

        hill_ptrs, hill_data, path_node_cnt, score_median, score_std = extract_hills(query_data, max_gap,
                                                                                     score_median + score_std * 3)
        logging.info(f'Number of hills {len(hill_ptrs):,}, len = {np.mean(path_node_cnt):.2f}')

        int_data = np.array(query_data['int_list_ms1'])

        hill_ptrs = split_hills(hill_ptrs, hill_data, int_data, hill_split_level=hill_split_level,
                                window=window)  # hill lenght is inthere already
        logging.info(f'After split hill_ptrs {len(hill_ptrs):,}')

        hill_data, hill_ptrs = filter_hills(hill_data, hill_ptrs, int_data,
                                            hill_check_large=hill_check_large, window=window)

        logging.info(f'After filter hill_ptrs {len(hill_ptrs):,}')

        stats, sortindex_, idxs_upper, scan_idx, hill_data, hill_ptrs = get_hill_data(query_data, hill_ptrs,
                                                                                      hill_data,
                                                                                      hill_nboot_max=hill_nboot_max,
                                                                                      hill_nboot=hill_nboot)
        logging.info('Extracting hill stats complete')

        pre_isotope_patterns = get_pre_isotope_patterns(stats, idxs_upper, sortindex_, hill_ptrs, hill_data,
                                                        int_data, scan_idx, maximum_offset,
                                                        iso_charge_min=iso_charge_min,
                                                        iso_charge_max=iso_charge_max,
                                                        iso_mass_range=iso_mass_range,
                                                        cc_cutoff=iso_corr_min)
        logging.info('Found {:,} pre isotope patterns.'.format(len(pre_isotope_patterns)))

        isotope_patterns, iso_idx, isotope_charges = get_isotope_patterns(pre_isotope_patterns, hill_ptrs,
                                                                          hill_data, int_data, scan_idx,
                                                                          stats, sortindex_, averagine_aa,
                                                                          isotopes,
                                                                          iso_charge_min=iso_charge_min,
                                                                          iso_charge_max=iso_charge_max,
                                                                          iso_mass_range=iso_mass_range,
                                                                          iso_n_seeds=iso_n_seeds,
                                                                          cc_cutoff=iso_corr_min,
                                                                          iso_split_level=iso_split_level,
                                                                          callback=None)
        logging.info('Extracted {:,} isotope patterns.'.format(len(isotope_charges)))

        feature_table = feature_finder_report(query_data, isotope_patterns, isotope_charges, iso_idx, stats,
                                              sortindex_, hill_ptrs, hill_data)

        logging.info('Report complete.')

    elif datatype == 'bruker':
        logging.info('Feature finding on {}'.format(file_input))
        feature_path = extract_bruker(file_input)
        feature_table = convert_bruker(feature_path)
        logging.info('Bruker featurer finder complete. Extracted {:,} features.'.format(len(feature_table)))
    else:
        raise NotImplementedError()

    # Calculate additional params
    feature_table['rt_length'] = feature_table['rt_end'] - feature_table['rt_start']
    feature_table['rt_right'] = feature_table['rt_end'] - feature_table['rt_apex']
    feature_table['rt_left'] = feature_table['rt_apex'] - feature_table['rt_start']
    feature_table['rt_tail'] = feature_table['rt_right'] / feature_table['rt_left']

    logging.info('Matching features to query data.')
    features = map_ms2(feature_table, query_data, **f_settings)

    logging.info('Saving feature table.')
    ms_file.write(feature_table, dataset_name="feature_table")
    logging.info('Feature table saved to {}'.format(file_hdf5))

    logging.info('Saving features.')
    ms_file.write(features, dataset_name="features")
    logging.info(f'Feature finding of file {file_input} complete.')


if __name__ == '__main__':
    args = Arguments()
    para = {
        'file_input': r'D:\test\masssearch\2_A4_Mix_1.raw',
        'file_hdf5': r'D:\test\masssearch\2_A4_Mix_1.ms_data.hdf',
        'file_output': "/t/test.csv",
    }

    args.add_argument('-out', nargs='+')
    args.add_parameter(para)
    para = args.parse_args()
    main(para)
