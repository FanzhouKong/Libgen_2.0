#!/usr/bin/python3
import sys
from mimas.file_io import spec_file
from mimas.helper import filename

"""
This script will extract spectra in a specific file which contains the text indicated.
"""


def main():
    if len(sys.argv) <= 1:
        print("""
Usage: python extract_spectra_from_msp.py filename_input filename_output item_to_extract.
        """)
        return
    filename_in = sys.argv[1]
    filename_out = sys.argv[2]
    items = " ".join(sys.argv[3:])

    fo = filename.smart_io(filename_out, "wt")

    for spec in spec_file.read_one_spectrum(filename_in, include_raw=1):
        if spec['raw'].find(items) > 0:
            fo.writelines(spec["raw"])


if __name__ == '__main__':
    main()
