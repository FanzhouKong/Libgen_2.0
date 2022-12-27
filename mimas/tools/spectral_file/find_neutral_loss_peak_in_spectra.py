#!/usr/bin/env python3
import numpy as np
import pandas as pd
from mimas.helper.arguments import Arguments
from mimas.file_io import spec_file
from mimas.spectra.similarity.tools import clean_spectrum


def main(parameter):
    result = []
    for spec in spec_file.read_one_spectrum(parameter.file_input):
        peaks = spec["peaks"]
        if len(peaks) == 0:
            continue
        spec_file.standardize_spectrum(spec, {
            "precursormz": [["precursor_mz"], None, float]
        })
        precursor_mz = spec["precursormz"]
        if precursor_mz is None:
            continue
        peak_mz_min = precursor_mz+parameter.neutral_loss-parameter.tolerance_da
        peak_mz_max = precursor_mz+parameter.neutral_loss+parameter.tolerance_da

        select_peak_id = np.where(np.bitwise_and(peaks[:, 0] >= peak_mz_min, peaks[:, 0] <= peak_mz_max))[0]
        if len(select_peak_id) == 0:
            continue

        # Normalize spectrum
        peaks_intensity_sum = np.sum(peaks[:, 1])
        if peaks_intensity_sum > 0:
            peaks[:, 1] /= peaks_intensity_sum

        intensity_sum = peaks[select_peak_id, 1].sum()
        mz_mean= (peaks[select_peak_id, 0]*peaks[select_peak_id, 1]).sum()/intensity_sum

        spec["peaks"] = ";".join(["{:.6f}".format(p[0])+":"+"{:.6f}".format(p[1]) for p in peaks])
        spec["intensity"] = intensity_sum
        spec["neutral_loss"] = mz_mean-spec["precursormz"]
        result.append(spec)

    df = pd.DataFrame(result)
    df.to_csv(parameter.file_output, index=False)

if __name__ == "__main__":
    args = Arguments(description="Find neutral loss peak in spectra.")
    args.add_argument("-file_input", type=str, help="File of spectra to analysis.")
    args.add_argument("-file_output", type=str, help="File of output.")

    args.add_argument("-neutral_loss", type=float, help="Neutral loss mass.")
    args.add_argument("-tolerance_da", type=float, help="Tolerance of neutral loss mass in Da.")

    para = {
        # "file_input": r"/home/yli/project/Tethys/data/lipid_shen/unk_lipid_org_20211118/unk_lipid-neg.msp",
        # "file_output": r"/home/yli/project/Tethys/result/2022_02/0212_neutral_loss_search/data/unk_lipid-neg-result.csv",
        "file_input": r"/home/yli/project/Tethys/data/library/mona-20211108_clean.msp",
        "file_output": r"/home/yli/project/Tethys/result/2022_02/0212_neutral_loss_search/data/mona-57.0214.csv",

        "neutral_loss": -57.0214,
        "tolerance_da": 0.02
    }
    args.add_argument_from_dictionary(para)
    para = args.parse_args()
    main(para)
