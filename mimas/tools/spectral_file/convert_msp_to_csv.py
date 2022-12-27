#!/usr/bin/python3
import pandas as pd

from mimas.file_io import spec_file
from mimas.helper.arguments import Arguments
from collections.abc import Iterable


"""
This script will convert msp to csv file.
"""


def main(para):
    fi = para["file_input"]
    fo = open(para["file_output"], "wt")

    output_first_line = True
    for spec_info in spec_file.read_one_spectrum(fi, parse_comments=True):
        output_info = {}
        if output_first_line:
            output_item = [x for x in spec_info.keys()]

        for key in spec_info:
            value = spec_info[key]
            if key == "peaks":
                output_info[key] = ";".join([f"{x[0]}:{x[1]}" for x in value])
            elif isinstance(value, str):
                output_info[key] = value
            elif isinstance(value, Iterable):
                output_info[key] = ",".join(map(str, value))
            else:
                output_info[key] = str(value)

        for item in output_info:
            if isinstance(output_info[item], list):
                output_info[item] = ";".join(output_info[item])

        df = pd.DataFrame(output_info, index=[0])
        df.to_csv(fo, sep=",", index=False, header=output_first_line, columns=output_item,)
        if output_first_line:
            output_first_line = False

    fo.close()


if __name__ == '__main__':
    args = Arguments()
    para = {
        'file_input': 'result/2022_03/0312_alphapept_vs_msdial/data/alphapept/NIH_Lip_Std_CSH_NEG_TissuePool_01/ms2_selected.msp',
        'file_output': "result/2022_03/0312_alphapept_vs_msdial/data/alphapept/NIH_Lip_Std_CSH_NEG_TissuePool_01/ms2_selected.csv",
    }

    # args.add_argument('-out', nargs='+')
    args.add_argument_from_dictionary(para)
    para = args.parse_args()
    main(para)
