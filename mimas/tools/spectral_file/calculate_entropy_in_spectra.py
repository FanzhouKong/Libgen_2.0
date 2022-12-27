import numpy as np
import pandas as pd
import typing

from file_io import spec_file
from helper.arguments import Arguments
import scipy.stats
from spectra.Daphnis.tools import clean_spectrum


def calculate_entropy(spectrum: typing.Union[list, np.ndarray]) -> float:
    entropy = scipy.stats.entropy([x[1] for x in spectrum])
    return entropy


def remove_noise(spec, noise, max_mz, ms2_da=None, ms2_ppm=None):
    spec = clean_spectrum(spec, max_mz=max_mz, ms2_da=ms2_da, ms2_ppm=ms2_ppm)
    if spec.shape[0] > 0:
        spec_max = np.max(spec[:, 1])
        spec = spec[spec[:, 1] > spec_max * noise]
        spec_sum = np.sum(spec[:, 1])
        spec[:, 1] /= spec_sum
    return spec


def main(para):
    result_all = []

    for spec_info in spec_file.read_one_spectrum(para["file_input"]):
        spec_info = spec_file.fill_spectrum_info_with_default_value(
            spec_info, ["spectrum_id", "precursor_mz"])
        spec = spec_info["spectrum"]
        entropy_no_noise = None
        peak_num_no_noise = None
        try:
            precursor_mz = float(spec_info["precursor_mz"])
            if precursor_mz > 0:
                spec_no_noise = remove_noise(spec, para["noise"],
                                             max_mz=precursor_mz - para["precursor_removal"],
                                             ms2_da=para.get("ms2_da", None), ms2_ppm=para.get("ms2_ppm", None))
                entropy_no_noise = calculate_entropy(spec_no_noise)
                peak_num_no_noise = len([x for x in spec_no_noise if x[1] > 0])
        except ValueError as e:
            print(e)

        result = {
            "spectrum_id": spec_info["spectrum_id"],
            'entropy': calculate_entropy(spec),
            'peak_num': len([x for x in spec if x[1] > 0]),

            'entropy_no_noise': entropy_no_noise,
            'peak_num_no_noise': peak_num_no_noise,

            'ce': spec_info.get("ce", ""),
            'mz': spec_info["precursor_mz"]
        }

        result_all.append(result)

    df = pd.DataFrame(result_all)

    with open("{path}/information.csv".format(path=para['path_output']), 'wt') as fo:
        df.to_csv(fo, index=False)


if __name__ == '__main__':
    args = Arguments()
    para = {
        "file_input": "/t/library/ALL_GNPS.mgf",
        'path_output': "/t/",

        "noise": 0.01,
        "ms2_da": 0.05,
        "precursor_removal": 1.6,
    }

    args.add_parameter(para)
    para = args.parse_args()
    main(para)
