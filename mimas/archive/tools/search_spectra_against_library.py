import numpy as np
import pandas as pd

from database.spectra_collector import SpectraCollector
from file_io import spec_file
from helper import multiplecore
from helper.arguments import Arguments
from spectra.similarity import spectral_similarity_simple
from spectra.similarity.tools import clean_spectrum


def preprocess_peaks(precursor_mz, peaks, noise, precursor_removal=None, ms2_da=None, ms2_ppm=None):
    if precursor_removal is not None:
        peaks = clean_spectrum(peaks, max_mz=precursor_mz - precursor_removal, ms2_da=ms2_da, ms2_ppm=ms2_ppm)
    else:
        peaks = clean_spectrum(peaks, max_mz=None, ms2_da=ms2_da, ms2_ppm=ms2_ppm)

    if peaks.shape[0] > 0:
        spec_max = np.max(peaks[:, 1])
        peaks = peaks[peaks[:, 1] > spec_max * noise]
        spec_sum = np.sum(peaks[:, 1])
        peaks[:, 1] /= spec_sum
    return peaks


def main(para):
    untarget_identity_search(para)


def untarget_identity_search(para):
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
    spec_library.preprocess_peaks(preprocess_peaks, para["threads"],
                                  para["noise"], para["precursor_removal"],
                                  para.get("ms2_da", None), para.get("ms2_ppm", None))

    # Score
    print("Start score.")
    # """
    parallel = multiplecore.MPRunner(func_run=_search_library,
                                     para_share=(spec_library, para),
                                     copy_shared_para=True,
                                     threads=para.get("threads", 1))
    # """
    result_all = []

    for spec_search_info in spec_file.read_one_spectrum(para["file_search"]):
        SpectraCollector.standardize_query_spectrum(spec_search_info)
        spec_query = {"scan_number": spec_search_info["scan_number"],
                      "name": spec_search_info["name"],
                      "precursor_mz": spec_search_info["precursor_mz"],
                      "precursor_type": spec_search_info["precursor_type"],
                      "peaks": spec_search_info["peaks"]}
        # result_all.append(_search_library(spec_query, spec_library, para))
        parallel.add_parameter_for_job((spec_query,), debug=0)

    result_all = parallel.get_result()
    # result_all = [y for x in result_all for y in x]

    # 4. Output result
    df = pd.DataFrame.from_dict(result_all)
    df.to_csv(para["file_output"], index=False)


def _search_library(spec_search_info, spec_library, para):
    # print(spec_search_info["spectrum_id"])
    result = []

    precursor_mz, precursor_type, peaks = \
        spec_search_info["precursor_mz"], spec_search_info["precursor_type"], \
        np.asarray(spec_search_info.pop("peaks"), dtype=np.float32, order="C")

    spec_candidate_info_all = spec_library.get_candidate_spectra(
        precursor_mz=precursor_mz,
        adduct=precursor_type,
        delta_da=para["ms1_da"], delta_ppm=para["ms1_ppm"]
    )
    if len(spec_candidate_info_all) == 0:
        return result

    # Add noise, the m/z for noise will be higher than precursor m/z.
    spec_search = preprocess_peaks(precursor_mz, peaks, para["noise"], para["precursor_removal"],
                                   ms2_da=para["ms2_da"], ms2_ppm=para["ms2_ppm"])

    # Calculate similarity score
    for spec_candidate_info in spec_candidate_info_all:
        similarity = spectral_similarity_simple.entropy_similarity(spec_search,
                                                                   spec_candidate_info[1],
                                                                   ms2_ppm=para["ms2_ppm"],
                                                                   ms2_da=para["ms2_da"],
                                                                   need_clean_spectra=False)
        if similarity >= para["score_min"]:
            result_cur = {"score": similarity}
            for item in spec_candidate_info[-1]:
                result_cur["library_" + item] = spec_candidate_info[-1][item]
            for item in spec_search_info:
                result_cur["query_" + item] = spec_search_info[item]
            result.append(result_cur)

    return result


if __name__ == '__main__':
    args = Arguments()
    para = {
        'threads': 20,
        "precursor_removal": 1.6,
        "ms1_ppm": 10,
        "ms2_da": 0.05,
        "ion_type": ["[M+H]+", "[M-H]-"],
        "score_method": [],

        "noise": 0.01,
        'file_search': '/home/yli/project/SimilarityConfidence/result/2020/10/1016_similarity_score_v8/library/nist20.msp',
        "file_library": "/home/yli/project/SimilarityConfidence/result/2020/10/1016_similarity_score_v8/library/nist20.msp",
        'path_output': "/home/yli/project/SimilarityConfidence/result/2020/10/1016_similarity_score_v8/test/",
    }
    para = {
        'threads': 1,
        "ms1_ppm": 10,
        "ms2_da": 0.05,
        "score_min": 0.5,

        'file_search': r"D:\test\db\test.msp",
        "file_library": r"D:\test\db\mona.msp",
        'file_output': r"D:\test\db\test.csv",
    }

    args.add_parameter(para)

    para = args.parse_args()

    # Check parameter
    para["ms1_da"] = para.get("ms1_da", None)
    if para["ms1_da"] is not None:
        para["ms1_da"] = float(para["ms1_da"])
        para["ms1_ppm"] = None
    else:
        para["ms1_ppm"] = float(para.get("ms1_ppm", 10.))

    para["ms2_ppm"] = para.get("ms2_ppm", None)
    if para["ms2_ppm"] is not None:
        para["ms2_ppm"] = float(para["ms2_ppm"])
        para["ms2_da"] = None
    else:
        para["ms2_da"] = float(para.get("ms2_da", 0.05))

    para["precursor_removal"] = float(para.get("precursor_removal", 1.6))
    para["noise"] = float(para.get("noise", 0.01))
    para["score_method"] = para.get("score_method", "entropy")
    main(para)
