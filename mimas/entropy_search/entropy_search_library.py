import numpy as np
from typing import Union
from pathlib import Path
from mimas.file_io import spec_file
from mimas.spectra.similarity import spectral_entropy, clean_spectrum, entropy_similarity, hybrid_similarity
from .spectra_library import SpectraLibrary


def entropy_search(search_type: str, file_search: Union[str, Path], spectral_library: SpectraLibrary, parameter: dict):
    """
    :param search_type: "identity", "shift", "open" or "hybrid".
    :param file_search: file name / file object of the search file.
    :param spectral_library: spectral library object.
    :param parameter: the following fields are required:
    {
        "score_min": -1,                # Minimum similarity score to be reported.
        "result_max": -1,               # Maximum number of results to be reported.
        "ms1_tolerance_in_ppm": 20,     # MS1 tolerance in ppm.
        "ms1_tolerance_in_da": 0.1,     # MS1 tolerance in Da. (optional)
        "ms2_tolerance_in_da": 0.02,    # MS2 tolerance in Da.
        "ms2_tolerance_in_ppm": 20,     # MS2 tolerance in ppm. (optional)
        "precursor_removal": 1.6,       # Ion removaled from precursor m/z before searching.
        "noise": 0.01,                  # Noise threshold. Noise is defined as the peak with the highest intensity multiplied by this value.

        "clean_spectra": True,          # Clean spectra before searching, need to be True.
    }

    """

    # Score
    if search_type == "identity":
        func_search = _identity_search_library
    elif search_type == "open":
        func_search = _open_search_library
    elif search_type == "hybrid":
        func_search = _hybrid_search_library
    elif search_type == "shift":
        func_search = _shift_search_library

    search_result = []
    query_spec_metadata = {}
    for spec_search_info in spec_file.read_one_spectrum(file_search, ms2_only=True):
        spec_query = _pre_process_query_spectrum(spec_search_info, parameter=parameter)

        if spec_query["precursor_mz"] > 0 and len(spec_query["peaks"]) > 0:
            scan_num, precursor_mz, charge, peaks = \
                spec_search_info["scan_number"], spec_search_info["precursor_mz"], spec_search_info["charge"], spec_search_info["peaks"]

            # Search this spectrum.
            cur_result = func_search(scan_num, precursor_mz, charge, peaks, spectral_library, parameter)
            if cur_result:
                if parameter["result_max"] > 0:
                    cur_result.sort(key=lambda x: x[-1], reverse=True)
                    cur_result = cur_result[:parameter["result_max"]]
                search_result += cur_result
                spec_query.pop("peaks")
                query_spec_metadata[spec_query["scan_number"]] = spec_query
            pass
        # break

    return search_result, query_spec_metadata


def _pre_process_query_spectrum(spec_info, parameter):
    spec_file.standardize_spectrum(spec_info, {
        "precursor_mz": [["precursormz"], -1, float],
        "adduct": [["precursortype", "precursor_type"], "", str],
        "name": [["title"], "", str],
        "id": [["spectrum_id", "db#", "spectrumid", "NISTNO"], None, None],
        "scan_number": [["_scan_number"], None, int],
        "rt": [["retention_time", "retentiontime"], None, float]
    })

    # If charge is given as paramter, use it. Otherwise, use the charge from adduct.
    charge = 1
    if "charge" in parameter and parameter["charge"]:
        charge = parameter["charge"]
    else:
        # If charge is written in the spectrum, use it.
        try:
            charge = int(spec_info["charge"])
        except:
            try:
                charge = {"+": 1, "-": -1}[spec_info["adduct"][-1]]
                if spec_info["adduct"][-2] != "]":
                    charge *= int(spec_info["adduct"][-2])
            except:
                print("Can't determine charge in spectrum {}, set it as +1.".format(spec_info["id"]))
    spec_info["charge"] = charge

    # Add spectral entropy
    peaks = clean_spectrum(
        spec_info["peaks"],
        max_mz=spec_info["precursor_mz"]-parameter["precursor_removal"],
        noise_threshold=parameter["noise"],
        max_peak_num=None,
        ms2_ppm=parameter.get("ms2_tolerance_in_ppm"),
        ms2_da=parameter.get("ms2_tolerance_in_da")
    )
    if "entropy" not in spec_info or (not isinstance(spec_info["entropy"], float)):
        spec_info["entropy"] = spectral_entropy(peaks, clean_spectra=False)
    if "entropy_normalized" not in spec_info or (not isinstance(spec_info["entropy_normalized"], float)):
        if len(peaks) > 1:
            spec_info["entropy_normalized"] = float(spec_info["entropy"])/np.log(len(peaks))
        else:
            spec_info["entropy_normalized"] = 0

    spec_info["peaks"] = peaks
    return spec_info


def _identity_search_library(scan_num, precursor_mz, charge, peaks, spec_library, parameter):
    spec_candidate_info_all = spec_library.get_candidate_spectra(
        precursor_mz=precursor_mz,
        charge=charge,
        delta_da=parameter.get("ms1_tolerance_in_da", None), delta_ppm=parameter.get("ms1_tolerance_in_ppm", None)
    )
    # Calculate similarity score
    result = []
    for spec_candidate_info in spec_candidate_info_all:
        similarity = entropy_similarity(peaks,
                                        spec_candidate_info[1],
                                        ms2_ppm=parameter.get("ms2_tolerance_in_ppm", None),
                                        ms2_da=parameter.get("ms2_tolerance_in_da", None),
                                        clean_spectra=False)
        if similarity >= parameter["score_min"]:
            result.append([scan_num, spec_candidate_info[-1], similarity])
    return result


def _hybrid_search_library(scan_num, precursor_mz, charge, peaks, spec_library, parameter):
    spec_candidate_info_all = spec_library.get_candidate_spectra_for_hybrid_search(
        precursor_mz=precursor_mz,
        charge=charge,
        peaks=peaks,
    )

    # Calculate similarity score
    result = []
    for spec_candidate_info in spec_candidate_info_all:
        similarity = hybrid_similarity.similarity(peaks, spec_candidate_info[1],
                                                  precursor_mz - spec_candidate_info[0],
                                                  ms2_ppm=parameter.get("ms2_tolerance_in_ppm", None),
                                                  ms2_da=parameter.get("ms2_tolerance_in_da", None),
                                                  clean_spectra=False)
        if similarity >= parameter["score_min"]:
            result.append([scan_num, spec_candidate_info[-1], similarity])
    return result


def _open_search_library(scan_num, precursor_mz, charge, peaks, spec_library, parameter):
    spec_candidate_info_all = spec_library.get_candidate_spectra_all(charge=charge)

    # Calculate similarity score
    result = []
    for spec_candidate_info in spec_candidate_info_all:
        similarity = entropy_similarity(spec_search,
                                        spec_candidate_info[1],
                                        ms2_ppm=parameter.get("ms2_tolerance_in_ppm", None),
                                        ms2_da=parameter.get("ms2_tolerance_in_da", None),
                                        clean_spectra=False)
        if similarity >= parameter["score_min"]:
            result.append([spec_search_info["scan_number"], spec_candidate_info[-1], similarity])
    return result


def _shift_search_library(scan_num, precursor_mz, charge, peaks, spec_library, parameter):
    spec_candidate_info_all_1 = spec_library.get_candidate_spectra(
        precursor_mz=precursor_mz,
        charge=charge,
        delta_da=parameter.get("ms1_tolerance_in_da", None), delta_ppm=parameter.get("ms1_tolerance_in_ppm", None)
    )
    spec_candidate_info_all_2 = spec_library.get_candidate_spectra(
        precursor_mz=precursor_mz - parameter["shift"],
        charge=charge,
        delta_da=parameter.get("ms1_tolerance_in_da", None), delta_ppm=parameter.get("ms1_tolerance_in_ppm", None)
    )
    spec_candidate_info_all = spec_candidate_info_all_1 + spec_candidate_info_all_2

    # Calculate similarity score
    result = []
    for spec_candidate_info in spec_candidate_info_all:
        similarity = hybrid_similarity.similarity(peaks_search, spec_candidate_info[1],
                                                  precursor_mz - spec_candidate_info[0],
                                                  ms2_ppm=parameter.get("ms2_tolerance_in_ppm", None),
                                                  ms2_da=parameter.get("ms2_tolerance_in_da", None),
                                                  clean_spectra=False)
        if similarity >= parameter["score_min"]:
            result.append([spec_search_info["scan_number"], spec_candidate_info[-1], similarity])
    return result
