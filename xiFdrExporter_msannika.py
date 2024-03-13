#!/usr/bin/env python3

# Exporter of MS Annika CSM Results to xiFDR input format
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import argparse
import pandas as pd

__version = "1.0.0"
__date = "20240313"

"""
DESCRIPTION:
A script to export MS Annika CSM results (.xlsx) to a xiFDR input file (.csv).
CSMs should be unfiltered, therefore include decoys and not be validated for any
FDR.
USAGE:
xiFdrExporter_msannika.py f [f]
                            [-o OUTPUT]
                            [-h]
                            [--version]
positional arguments:
  f                     Crosslink-Spectrum-Matches (CSMs) exported from
                        MS Annika in Microsoft Excel (.xlsx) format.
optional arguments:
  -o OUTPUT, --output OUTPUT
                        Prefix of the output file.
  -h, --help            show this help message and exit
  --version             show program's version number and exit
"""

# Exporter class with constructor that takes one MS Annika CSM result file as
# input. CSMs should not be in any way filtered and exported to Microsoft Excel
# .xlsx format from Proteome Discoverer.
class MSAnnika_Exporter:

    def __init__(self, input_file: str):
        self.input_file = input_file

    # static method to generate pandas dataframe of xiFDR export without class
    # instance. Takes the file name of the CSM file as input.
    @staticmethod
    def generate_df(input_file: str) -> pd.DataFrame:

        df = pd.read_excel(input_file)
        df.rename(columns = {"Spectrum File": "run",
                             "First Scan": "scan",
                             "Sequence A": "peptide1",
                             "Sequence B": "peptide2",
                             "Crosslinker Position A": "peptide link 1",
                             "Crosslinker Position B": "peptide link 2",
                             "Charge": "precursor charge",
                             "Combined Score": "score",
                             "Score Alpha": "peptide1 score",
                             "Score Beta": "peptide2 score",
                             "Accession A": "accession1",
                             "Accession B": "accession2",
                             "A in protein": "peptide position 1",
                             "B in protein": "peptide position 2"},
                  inplace = True,
                  errors = "raise")
        df["is decoy 1"] = df["Alpha T/D"].apply(lambda x: "false" if "t" in str(x).lower() else "true")
        df["is decoy 2"] = df["Beta T/D"].apply(lambda x: "false" if "t" in str(x).lower() else "true")
        df["peptide position 1"] = df["peptide position 1"].apply(lambda x: ";".join([str(int(y) + 1) for y in x.split(";")]))
        df["peptide position 2"] = df["peptide position 2"].apply(lambda x: ";".join([str(int(y) + 1) for y in x.split(";")]))

        return df

    # classmethod implementation of the static generate_df
    def __generate_csv_df(self) -> pd.DataFrame:
        return self.generate_df(self.input_file)

    # export function, takes one argument "output_file" which sets the prefix
    # of generated output file
    def export(self, output_file: str = None) -> pd.DataFrame:
        csv = self.__generate_csv_df()

        if output_file is None:
            output_file = ".".join(self.input_file.split(".")[:-1])

        csv.to_csv(output_file + "_xiFDR.csv", index = False)

        return csv

# initialize exporter and export xiFDR csv file
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(metavar = "f",
                        dest = "file",
                        help = "Name/Path of the MS Annika CSM result file (in .xlsx format) to process.",
                        type = str,
                        nargs = 1)
    parser.add_argument("-o", "--output",
                        dest = "output",
                        default = None,
                        help = "Prefix of the output file.",
                        type = str)
    parser.add_argument("--version",
                        action = "version",
                        version = __version)
    args = parser.parse_args()

    exporter = MSAnnika_Exporter(args.file[0])

    exporter.export(args.output)

if __name__ == "__main__":
    main()
