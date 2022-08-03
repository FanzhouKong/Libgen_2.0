
import alphapept.feature_finding
import alphapept.performance
from alphapept.load_mzml_data import load_mzml_data
from alphapept.constants import averagine_aa, isotopes
import numpy as np

def find_features(ms_file):
    f_settings = {'centroid_tol': 8,
                  'hill_check_large': 40,
                  'hill_length_min': 3,
                  'hill_nboot': 150,
                  'hill_nboot_max': 300,
                  'hill_smoothing': 1,
                  'hill_split_level': 1.3,
                  'iso_charge_max': 3,
                  'iso_charge_min': 1,
                  'iso_corr_min': 0.6,
                  'iso_mass_range': 5,
                  'iso_n_seeds': 100,
                  'iso_split_level': 1.3,
                  'map_mob_range': 0.3,
                  'map_mz_range': 0.5,
                  'map_n_neighbors': 5,
                  'map_rt_range': 0.5,
                  'max_gap': 2,
                  'search_unidentified': True}
    alphapept.performance.set_worker_count(1)
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

    ms_file["int_list_ms1"] = np.array(ms_file["int_list_ms1"]).astype(int)
    int_data = ms_file['int_list_ms1']
    # logging.info('Feature finding on {}'.format(file_name))
    # logging.info(f'Hill extraction with centroid_tol {centroid_tol} and max_gap {max_gap}')

    hill_ptrs, hill_data, path_node_cnt, score_median, score_std = alphapept.feature_finding.extract_hills(
        ms_file, max_gap, centroid_tol)

    # logging.info(f'Number of hills {len(hill_ptrs):,}, len = {np.mean(path_node_cnt):.2f}')
    # logging.info(f'Repeating hill extraction with centroid_tol {score_median+score_std*3:.2f}')

    hill_ptrs, hill_data, path_node_cnt, score_median, score_std = alphapept.feature_finding.extract_hills(
        ms_file, max_gap, score_median + score_std * 3)

    # logging.info(f'Number of hills {len(hill_ptrs):,}, len = {np.mean(path_node_cnt):.2f}')

    hill_ptrs, hill_data = alphapept.feature_finding.remove_duplicate_hills(hill_ptrs, hill_data, path_node_cnt)
    # logging.info(f'After duplicate removal of hills {len(hill_ptrs):,}')

    hill_ptrs = alphapept.feature_finding.split_hills(
        hill_ptrs, hill_data, int_data, hill_split_level=hill_split_level, window=window)  # hill lenght is inthere already
    # logging.info(f'After split hill_ptrs {len(hill_ptrs):,}')

    hill_data, hill_ptrs = alphapept.feature_finding.filter_hills(
        hill_data, hill_ptrs, int_data, hill_check_large=hill_check_large, window=window)

    # logging.info(f'After filter hill_ptrs {len(hill_ptrs):,}')

    stats, sortindex_, idxs_upper, scan_idx, hill_data, hill_ptrs = alphapept.feature_finding.get_hill_data(
        ms_file, hill_ptrs, hill_data, hill_nboot_max=hill_nboot_max, hill_nboot=hill_nboot)
    # logging.info('Extracting hill stats complete')

    pre_isotope_patterns = alphapept.feature_finding.get_pre_isotope_patterns(
        stats, idxs_upper, sortindex_, hill_ptrs, hill_data, int_data, scan_idx, alphapept.feature_finding.maximum_offset,
        iso_charge_min=iso_charge_min, iso_charge_max=iso_charge_max, iso_mass_range=iso_mass_range, cc_cutoff=iso_corr_min)
    # logging.info('Found {:,} pre isotope patterns.'.format(len(pre_isotope_patterns)))

    isotope_patterns, iso_idx, isotope_charges = alphapept.feature_finding.get_isotope_patterns(
        pre_isotope_patterns, hill_ptrs, hill_data, int_data, scan_idx, stats, sortindex_, averagine_aa, isotopes,
        iso_charge_min=iso_charge_min, iso_charge_max=iso_charge_max, iso_mass_range=iso_mass_range,
        iso_n_seeds=iso_n_seeds, cc_cutoff=iso_corr_min, iso_split_level=iso_split_level, callback=None)
    # logging.info('Extracted {:,} isotope patterns.'.format(len(isotope_charges)))

    feature_table, lookup_idx = alphapept.feature_finding.feature_finder_report(
        ms_file, isotope_patterns, isotope_charges, iso_idx, stats, sortindex_, hill_ptrs, hill_data)
    # logging.info('Report complete.')

    # Calculate additional params
    feature_table['rt_length'] = feature_table['rt_end'] - feature_table['rt_start']
    feature_table['rt_right'] = feature_table['rt_end'] - feature_table['rt_apex']
    feature_table['rt_left'] = feature_table['rt_apex'] - feature_table['rt_start']
    feature_table['rt_tail'] = feature_table['rt_right'] / feature_table['rt_left']

    # logging.info('Matching features to query data.')

    if 'mono_mzs2' not in ms_file.keys():
        # logging.info('No MS2-data to match.')
        features = pd.DataFrame()
    else:
        features = alphapept.feature_finding.map_ms2(feature_table, ms_file, **f_settings)

    return feature_table, features
