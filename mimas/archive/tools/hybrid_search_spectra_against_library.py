import numpy as np
import pandas as pd

from database.spectra_collector import SpectraCollector
from file_io import spec_file
from helper import multiplecore
from helper.arguments import Arguments
from spectra.Daphnis import hybrid_similarity
from spectra.Daphnis.tools import clean_spectrum


def remove_noise(spec, noise, max_mz, ms2_da=None, ms2_ppm=None):
    spec = clean_spectrum(spec, max_mz=max_mz, ms2_da=ms2_da, ms2_ppm=ms2_ppm)
    if spec.shape[0] > 0:
        spec_max = np.max(spec[:, 1])
        spec = spec[spec[:, 1] > spec_max * noise]
        spec_sum = np.sum(spec[:, 1])
        spec[:, 1] /= spec_sum
    return spec


def remove_noise_in_spectral_library(spec_list, para):
    ms2_da = para.get("ms2_da", None)
    ms2_ppm = para.get("ms2_ppm", None)

    spec_list_new = multiplecore.run_multiple_process(
        remove_noise,
        all_para_individual=([(x["spectrum"], para["noise"],
                               x["precursormz"] - para["precursor_removal"],
                               ms2_da, ms2_ppm) for x in spec_list]),
        threads=para.get("threads", 1))

    for spec_info, spec in zip(spec_list, spec_list_new):
        spec_info["spectrum"] = spec
    return spec_list


def main(para):
    """
    :param para: the following fields are required:
            file_library: The pathname for all library files.
            file_search: The .msp file which contains spectra to search.
            path_output: The result file.
            ms2_da or ms2_ppm
            method: the name for similarity
    """

    # Read spectral library.
    spec_library = SpectraCollector()
    spec_library.read_spectra_file(para["file_library"])
    spec_library.spectra = remove_noise_in_spectral_library(spec_library.spectra, para)

    # Score
    print("Start score.")

    result = []
    for i, spec_search_info in enumerate(spec_file.read_one_spectrum(para["file_search"])):
        SpectraCollector().organize_spec_info(spec_search_info, i)
        result_cur = _search_library(spec_search_info, spec_library, para)
        result += result_cur
        if i > 10:
            break

    df = pd.DataFrame(result)
    df.to_csv(para["file_output"], index=False)


def _search_library(spec_search_info, spec_library, para):
    # print(spec_search_info["spectrum_id"])
    result = []

    # Add noise, the m/z for noise will be higher than precursor m/z.
    spec_search = remove_noise(
        spec_search_info["spectrum"], para["noise"],
        max_mz=float(spec_search_info["precursormz"]) - para["precursor_removal"],
        ms2_da=para.get("ms2_da", None), ms2_ppm=para.get("ms2_ppm", None))

    # Calculate similarity score
    for id, spec_library_info in enumerate(spec_library.spectra):
        if para['file_search'] == para['file_library'] and \
                spec_search_info["spectrum_id"] == spec_library_info["spectrum_id"]:
            continue

        similarity_score = hybrid_similarity.similarity(
            spec_search, spec_library_info["spectrum"],
            float(spec_library_info["precursormz"]) - float(spec_search_info["precursormz"]),
            ms2_ppm=para.get("ms2_ppm", None),
            ms2_da=para.get("ms2_da", None),
            need_clean_spectra=False)
        if similarity_score > 0:
            result.append([id, similarity_score])

    result.sort(key=lambda x: -x[1])
    if len(result) > para["output_number"]:
        score_cutoff = result[para["output_number"]][1]
    else:
        score_cutoff = 0
    result_final = []
    for r in result:
        if r[1] < score_cutoff:
            break
        spec_library_info = spec_library.spectra[r[0]]
        result_cur = {
            'spec_id_query': spec_search_info["spectrum_id"],
            'inchikey_query': spec_search_info["inchikey"],
            'spec_id_library': spec_library_info["spectrum_id"],
            'inchikey_library': spec_library_info["inchikey"],
            "ion_type": spec_library_info["precursortype"],
            "similarity_score": r[1]
        }
        result_final.append(result_cur)

    return result_final


if __name__ == '__main__':
    args = Arguments()
    para = {
        'threads': 22,
        "precursor_removal": 1.6,
        "ms2_da": 0.05,
        "ion_type": [],
        "score_method": [],
        "output_number": 20,

        "noise": 0.01,
        "file_search": "/home/yli/project/SimilarityConfidence/result/2020/10/1028_hybrid_search/library/alberto.msp",
        "file_library": "/home/yli/project/SimilarityConfidence/result/2020/10/1028_hybrid_search/library/nist20.msp",
        'file_output': "/home/yli/project/SimilarityConfidence/result/2020/10/1028_hybrid_search/test.csv",
    }

    args.add_parameter(para)

    para = args.parse_args()
    main(para)
