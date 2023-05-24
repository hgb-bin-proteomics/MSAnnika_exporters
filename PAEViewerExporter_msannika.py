#!/usr/bin/env python3

# Exporter of MS Annika Crosslink Results to PAE Viewer input format
# 2023 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import argparse
import pandas as pd

__version = "1.0.0"
__date = "20230523"

"""
DESCRIPTION:
A script to export MS Annika results from pyXlinkViewer to PAE Viewer input
files (CSV).
USAGE:
PAEViewerExporter_msannika.py f [f]
                                [-t DISTANCE]
                                [-o OUTPUT]
                                [-h]
                                [--version]
positional arguments:
  f                     Crosslinks exported from pyXlinkViewer in csv format.
optional arguments:
  -t DISTANCE, --threshold DISTANCE
                        threshold (float) that specifies if a crosslink
                        satisfies the crosslinker-specific distance constraint.
  -o OUTPUT, --output OUTPUT
                        Prefix of the output file.
  -h, --help            show this help message and exit
  --version             show program's version number and exit
"""

# Exporter class with constructor that takes one pyXlinkViewer CSV as input, and
# optional the threshold to be considered for constraint violation
class MSAnnika_Exporter:

    def __init__(self, input_file: str, threshold: float = None):
        self.input_file = input_file
        self.threshold = threshold

    def __get_constraint(self, row: pd.Series, threshold: float = None) -> bool:

        if threshold is not None:
            return float(row["CA distance"]) <= float(threshold)

        return str(row["Sat-Viol"]) == "S"

    # transform input pyXlinkViewer file to PAE Viewer input file
    def __generate_csv_df(self) -> pd.DataFrame:

        df = pd.read_csv(self.input_file)
        df_paeviewer = pd.DataFrame({"Protein1": df["Chain 1"],
                                     "SeqPos1": df["Residue 1"],
                                     "Protein2": df["Chain 2"],
                                     "SeqPos2": df["Residue 2"],
                                     "RestraintSatisfied": df.apply(lambda row: self.__get_constraint(row, threshold = self.threshold), axis = 1)})

        return df_paeviewer

    # export function, takes one argument "output_file" which sets the prefix
    # of generated output file
    def export(self, output_file: str = None) -> pd.DataFrame:
        csv = self.__generate_csv_df()

        if output_file == None:
            output_file = self.input_file.split(".")[0]

        csv.to_csv(output_file + "_PAEViewer.csv", index = False)

        return csv

# initialize exporter and export PAE Viewer file
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(metavar = "f",
                        dest = "file",
                        help = "Name/Path of the pyXlinkViewer result file to process.",
                        type = str,
                        nargs = 1)
    parser.add_argument("-t", "--threshold",
                        dest = "threshold",
                        default = None,
                        help = "Threshold that specifies if a crosslink satisfies the crosslinker-specific distance constraint.",
                        type = float)
    parser.add_argument("-o", "--output",
                        dest = "output",
                        default = None,
                        help = "Name of the output file.",
                        type = str)
    parser.add_argument("--version",
                        action = "version",
                        version = __version)
    args = parser.parse_args()

    exporter = MSAnnika_Exporter(args.file[0], args.threshold)

    exporter.export(args.output)

if __name__ == "__main__":
    main()
