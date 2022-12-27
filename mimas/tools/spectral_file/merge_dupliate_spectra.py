#!/usr/bin/env python3
import numpy as np
import tempfile
from pathlib import Path

from mimas.helper.arguments import Arguments
from mimas.file_io import spec_file
from mimas.spectra.similarity.tools import clean_spectrum, match_two_spectra
from mimas.spectra.similarity import entropy_similarity


def main(parameter):
    file_input = parameter.file_input
    file_output = tempfile.NamedTemporaryFile(delete=True).name+".msp"

    has_duplicate_spectra = True
    while has_duplicate_spectra:
        print("Write to file: "+file_output)
        has_duplicate_spectra = False
        with open(file_output, 'wt') as fo:
            all_good_spectra = []
            for spec in spec_file.read_one_spectrum(file_input):
                peaks = spec["peaks"]
                if len(peaks) > 0:
                    spec_file.standardize_spectrum(spec, {
                        "precursormz": [["precursor_mz"], None, float]
                    })
                    precursor_mz = spec["precursormz"]
                    if precursor_mz is not None:
                        spec["peaks"] = clean_spectrum(peaks, precursor_mz-1.6, noise_threshold=0.01, ms2_da=parameter.ms2_tolerance_da)
                        all_good_spectra.append(spec)
                        continue
                spec_file.write_one_spectrum(fo, spec, output_type=".msp")

            all_good_spectra.sort(key=lambda x: x["precursormz"])

            # Go through all spectra and merge duplicates
            all_spec_idx = {x["_scan_number"] for x in all_good_spectra}
            for spec_idx, spec in enumerate(all_good_spectra):
                if spec["_scan_number"] not in all_spec_idx:
                    continue
                # Find all spectra that are in the same precursor ion range
                precursor_mz_max = spec["precursormz"] * (1+2*parameter.ms1_tolerance_ppm*1e-6)
                similar_spectra = []
                for spec_j in all_good_spectra[spec_idx+1:]:
                    if spec_j["precursormz"] > precursor_mz_max:
                        break
                    similarity = entropy_similarity(spec["peaks"], spec_j["peaks"], ms2_da=parameter.ms2_tolerance_da, clean_spectra=False)
                    if similarity >= parameter.similarity_threshold:
                        similar_spectra.append(spec_j)

                # Merge all similar spectra
                if similar_spectra:
                    has_duplicate_spectra = True
                    for x in similar_spectra:
                        if x["_scan_number"] in all_spec_idx:
                            all_spec_idx.remove(x["_scan_number"])
                    total_spectra_number = len(similar_spectra)+1
                    peaks = np.copy(spec["peaks"])
                    peaks[:, 1] = peaks[:, 1]/total_spectra_number
                    for spec_j in similar_spectra:
                        spec_matched = match_two_spectra(peaks, spec_j["peaks"], ms2_da=parameter.ms2_tolerance_da)
                        spec_matched[:, 2] /= total_spectra_number
                        peaks = spec_matched[:, :2]
                        peaks[:, 1] += spec_matched[:, 2]
                        peaks = clean_spectrum(peaks, spec["precursormz"]-1.6, noise_threshold=0.01, ms2_da=parameter.ms2_tolerance_da)

                    # Generate new spectrum
                    similar_spectra.append(spec)
                    spec_merged = {"peaks": peaks}
                    for info in spec:
                        if info != "peaks":
                            value = {x[info] for x in similar_spectra if info in x}
                            if isinstance(spec[info], float):
                                spec_merged[info] = np.mean([float(x) for x in value])
                            else:
                                spec_merged[info] = "; ".join([str(x) for x in value])
                    spec = spec_merged
                else:
                    # Output the spectrum
                    all_spec_idx = all_spec_idx-{spec["_scan_number"]}

                # Output spec
                spec_file.write_one_spectrum(fo, spec, output_type=".msp")

        file_input = file_output
        file_output = tempfile.NamedTemporaryFile(delete=True).name+".msp"
    Path(file_input).rename(parameter.file_output)


if __name__ == "__main__":
    args = Arguments(description="Find neutral loss peak in spectra.")
    args.add_argument("-file_input", type=str, help="File of spectra to analysis.")
    args.add_argument("-file_output", type=str, help="File of output.")

    args.add_argument("-ms1_tolerance_ppm", type=float, help="Tolerance of MS2 spectrum in ppm.")
    args.add_argument("-ms2_tolerance_da", type=float, help="Tolerance of MS2 spectrum in Da.")
    args.add_argument("-similarity_threshold", type=float, help="Threshold of similarity.")

    para = {
        "file_input": r"/home/yli/project/Tethys/data/lipid_shen/unk_lipid_org_20211118/unk_lipid-pos.msp",
        # "file_input": r"/home/yli/project/Tethys/result/2022_02/0213_merge_similar_spectra/data/unk_lipid-neg-merged.msp",
        "file_output": r"/home/yli/project/Tethys/result/2022_02/0213_merge_similar_spectra/data/unk_lipid-pos-merged.msp",

        "ms1_tolerance_ppm": 10,
        "ms2_tolerance_da": 0.05,
        "similarity_threshold": 0.8,
    }
    args.add_argument_from_dictionary(para)
    para = args.parse_args()
    main(para)
