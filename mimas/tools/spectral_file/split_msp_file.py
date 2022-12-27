#!/usr/bin/python3
import pandas as pd

from file_io import spec_file
from helper.arguments import Arguments

"""
This script will split one .msp file into multiple .msp file.
"""


def main(para):
    fi = para["file_input"]
    file_num = 0
    spec_input = spec_file.read_one_spectrum(fi, include_raw=True)
    while True:
        fo = open(para["file_output_pre"] + str(file_num) + ".msp", "wt")
        spec_num = 0
        try:
            while True:
                spec = next(spec_input)
                fo.write(spec["raw"])
                spec_num += 1
                if spec_num == para["spec_num"]:
                    break
            fo.close()
            file_num += 1
        except StopIteration:
            break


if __name__ == '__main__':
    args = Arguments()
    para = {
        'file_input': r'd:/test/db/test.msp',
        'file_output_pre': "d:/test/db/test_sub_",
        "spec_num": 100,
    }

    args.add_parameter(para)
    para = args.parse_args()
    main(para)
