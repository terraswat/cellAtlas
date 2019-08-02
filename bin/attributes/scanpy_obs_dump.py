"""
This dumps the sample/cell annotations in out scanpy object to a tab delimeted file. 

The output has attributes in the columns and cell ids in the rows.

and number of values.

example usage:
    python3 scanpy_obs_dump.py -i /path/to/scanpyObj -o /path/to/output.tab
"""
import argparse
import sys
import scanpy.api as sc

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
    return df


def write_data(out_file, df):
    df.to_csv(out_file, sep="\t")


def main():

    in_file, out_file = parse_args()

    data = read_data(in_file)

    new_data = data_transform(data)

    write_data(out_file, new_data)


if __name__ == "__main__":
    sys.exit(main())
