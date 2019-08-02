"""
Convert mouse to human gene IDs for full fat tsv
example usage:
    python mouse_to_human_ortho.py -i /path/to/input -o /path/to/output
"""
import argparse
import sys


def parse_args():

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i',"--in_file", type=str, help=""
                        )

    parser.add_argument('-c',"--convert_file", type=str,
                        help="mouse to human orthology conversion, \
                              mouse=\"Gene name\" and \
                              human=\"Human gene name\"",
                        default="mouse_human_orthologs.txt"
                        )

    parser.add_argument('-o',"--out_file", type=str, help="",
                        )

    parser.add_argument('-s',"--sep", type=str,
                        help="in/out separator",
                        default="\t"
                        )

    parser.add_argument('-d',"--detectinputsep", type=str,
                        help="detect input separators if not same as output; \
                              note that this will be slower to read data.",
                        default=False
                        )

    opts = parser.parse_args()
    in_file, convert_file, out_file, separator, detect_input_separator = \
        opts.in_file, opts.convert_file, opts.out_file, opts.sep, \
        opts.detectinputsep

    return in_file, convert_file, out_file, separator, detect_input_separator


def read_data(in_file, convert_file, separator, detect_input_separator):
    import pandas as pd
    # if detect_input_separator is true then auto detect using python not c
    if detect_input_separator == True:
        data = pd.read_csv(in_file, sep=None)
        convert = pd.read_csv(convert_file, sep=None)
    else:
        data = pd.read_csv(in_file, sep=separator, index_col=0)
        convert = pd.read_csv(convert_file, sep=separator)
    return data, convert

def data_transform(data_obj, convert_obj):
    import pandas as pd
    # generate mouse to human conversion lookup dict
    mm = convert_obj["Gene name"]
    hh = convert_obj["Human gene name"]
    mm_to_hh = dict(zip(mm, hh))

    # initialize new column names
    data_mm = data_obj.columns.values[:]
    data_hh = ["NA"]*len(data_mm)

    # fill new column names where possible
    for i, datum_mm in enumerate(data_mm):
        if datum_mm in mm_to_hh:
            data_hh[i] = mm_to_hh[datum_mm]

    # output
    data_obj.columns = data_hh
    return data_obj

def write_data(out_file, data_obj, separator):
    import pandas as pd
    data_obj.to_csv(out_file, sep=separator)

def main():

    in_file, convert_file, out_file, separator, detect_input_separator = parse_args()

    data, convert = read_data(in_file, convert_file, separator, detect_input_separator)

    new_data = data_transform(data, convert)

    write_data(out_file, new_data, separator)

if __name__ == "__main__":
    sys.exit(main())
