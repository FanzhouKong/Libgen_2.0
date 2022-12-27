from pathlib import Path
import pandas as pd
from mimas.helper.arguments import Arguments
from mimas.helper.fileio import md5_for_object, md5_for_file
from entropy_search_library import entropy_search
from spectra_library import SpectraLibrary


def get_spectra_library(all_file_library, parameter_for_search):
    """
    Get the spectra library
    """
    for file_name in all_file_library:
        file = Path(file_name)
        index_parameter = {
            "calibrate_mz": parameter_for_search["calibrate_mz"],
            "clean_spectra": parameter_for_search["clean_spectra"],
            "ms2_tolerance_in_da": parameter_for_search["ms2_tolerance_in_da"],
            "index_for_hybrid_search": True,
        }
        parameter_md5 = md5_for_object(index_parameter)[:4]
        file_index = file.with_suffix(f".{parameter_md5}.index")
        spectra_library = SpectraLibrary(parameter=index_parameter)
        try:
            spectra_library.read_index_from_file(file_index)
        except:
            print(f"Index file {file_index} not found, re-indexing...")
            spectra_library.build_index_from_spectral_library(name=file.stem, file_spectra=file)
            spectra_library.write_index_to_file(file_index)
        yield spectra_library


def add_metadata(result, query_metadata, library_metadata):
    info = {"score": result[2]}
    for key in query_metadata:
        info[f"query_{key}"] = query_metadata[key]
    for key in library_metadata:
        if key != "metadata":
            info[f"library_{key}"] = library_metadata[key]
    return info


def main(parameter):
    """
    Common parameter: {
        "threads": 1,
        "precursor_removal": 1.6,
        "ms2_da": 0.05,
        "noise": 0.01,
        "score_min": 0.5,
    }
    """
    # Prepare the search parameter
    parameter_for_search = {
        "score_min": -1,
        "result_max": -1,
        "precursor_removal": 1.6,
        "noise": 0.01,
        "unique_result": 0,

        "ms1_tolerance_in_ppm": 20,
        "ms2_tolerance_in_da": 0.02,

        "clean_spectra": True,
        "calibrate_mz": True,
    }
    for item in ["shift",
                 "score_min", "result_max", "precursor_removal", "noise", "unique_result"
                 "ms1_tolerance_in_ppm", "ms1_tolerance_in_da", "ms2_tolerance_in_ppm", "ms2_tolerance_in_da",
                 "clean_spectra"]:
        if item in parameter:
            parameter_for_search[item] = parameter[item]

    # Prepare the library
    if not isinstance(parameter.file_library, list):
        all_file_library = [parameter.file_library]
    else:
        all_file_library = parameter.file_library

    final_result = []
    for spectral_library in get_spectra_library(all_file_library, parameter_for_search):
        matching_score, query_metadata = entropy_search(
            search_type=parameter.method, file_search=parameter.file_search, spectral_library=spectral_library, parameter=parameter_for_search)
        for score in matching_score:
            result = add_metadata(score, query_metadata[score[0]], spectral_library.get_metadata(score[1]))
            result.update({"library_name": spectral_library.name})
            final_result.append(result)

    df = pd.DataFrame(final_result)
    df.to_csv(parameter.file_output, index=False)


if __name__ == '__main__':
    args = Arguments()
    para_shared = {
        "precursor_removal": 1.6,
        "score_min": 0.75,
        "noise": 0.01,
        "charge": -1,
        "unique_result": 0,

        "ms2_tolerance_in_da": 0.05,

        'file_search': r"result/2022_03/0312_alphapept_vs_msdial/data/alphapept/NIH_Lip_Std_CSH_NEG_TissuePool_01/ms2_selected.msp",
        "file_library": [
            r"/home/yli/project/Database/SpectralDatabase_v5.3/data/nist20/nist20-20220421.msp",
            r"/home/yli/project/Database/SpectralDatabase_v5.3/data/gnps/gnps-20220421.msp",
            r"/home/yli/project/Database/SpectralDatabase_v5.3/data/mona/mona-20220421.msp",
        ],
        'file_output': r"result/2022_03/0312_alphapept_vs_msdial/data/alphapept/NIH_Lip_Std_CSH_NEG_TissuePool_01/ms2_selected-identity.csv",
    }
    para_untarget_hybrid = {
        "method": "hybrid",
    }
    para_untarget_identity = {
        "method": "identity",
        "ms1_tolerance_in_ppm": 20,
    }
    para_target_shift = {
        "method": "shift",
        "shift": 162.0533,
    }
    para_untarget_open = {
        "method": "open",
    }
    para_test = {
        "file_library": r"result/2022_05/0506_hybrid_search_on_nist20_linear_spectra/data/spec_test_1.msp",
        "file_search": r"result/2022_05/0506_hybrid_search_on_nist20_linear_spectra/data/spec_test_1.msp",
        "file_output": r"result/2022_05/0506_hybrid_search_on_nist20_linear_spectra/data/spec_test_1-out.csv",
        "method": "hybrid",
    }

    para = para_shared
    para.update(para_untarget_identity)
    para.update(para_test)

    args.add_argument_from_dictionary(para)
    para = args.parse_args()
    main(para)
