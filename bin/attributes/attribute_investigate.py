"""
This is to help investigate which attributes we want for
attribute enrichment.

It makes a tab delimeted file with columns attribute name, dtype,
and number of values.

example usage:
    python3 attribute_investigate.py -i /path/to/scanpyObj -o /path/to/output.tab
"""
import argparse
import sys
import scanpy.api as sc
import pandas as pd

def parse_args():

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i',"--in_file", type=str, help=""
                        )

    parser.add_argument('-o',"--out_file", type=str, help="",
                        default="full_sim.euc.tab"
                        )

    opts = parser.parse_args()
    in_file, out_file = opts.in_file, opts.out_file

    return in_file, out_file


def read_data(in_file):
    """
    Returns a dataframe from a path to a scanpy object. 
    """
    return sc.read(in_file).obs


def data_transform(df):

    attr_dict = {
        "names": df.columns,
        "dtypes": df.dtypes,
        "n_unique": df.apply(lambda x: len(x.unique())), 
        "example": df.apply(lambda x: x.unique()[0]) 
    }
    return pd.DataFrame(attr_dict)


def write_data(out_file, df):
    df.to_csv(out_file, sep="\t", index=None)


def main():

    in_file, out_file = parse_args()

    data = read_data(in_file)

    new_data = data_transform(data)

    write_data(out_file, new_data)


if __name__ == "__main__":
    sys.exit(main())
