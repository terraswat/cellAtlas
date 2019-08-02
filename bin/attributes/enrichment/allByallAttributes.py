import argparse
import sys
import pandas as pd
from pairwiseStats import allbyallStats
from data_type_dict import data_type_dictionary
def parse_args():

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i',"--in_file", type=str, help=""
                        )

    parser.add_argument('-o',"--out_file", type=str, help="",
                        default="attrXattr.tab"
                        )
    opts = parser.parse_args()
    in_file, out_file = opts.in_file, opts.out_file
    return in_file, out_file


def main():

    in_file, out_file = parse_args()

    attrdf = pd.read_table(in_file, index_col=0)
    data_type_dict = data_type_dictionary(attrdf)
    attrdf = allbyallStats(attrdf, data_type_dict)
    attrdf.to_csv(out_file, sep="\t")

if __name__ == "__main__":
    sys.exit(main())