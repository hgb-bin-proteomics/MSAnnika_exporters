#!/usr/bin/env python3

# Exporter of MS Annika Crosslink Results to xiNET Format
# 2022 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import argparse
import pandas as pd
from Bio import SeqIO
from typing import List

__version = "1.1.0"
__date = "20250121"

"""
DESCRIPTION:
A script to export MS Annika results to xiNET input files (CSV + FASTA).
USAGE:
xiNetExporter_msannika.py f [f ...]
                            [-fasta FASTA]
                            [-ignore IGNORE]
                            [-o OUTPUT]
                            [-h]
                            [--version]
positional arguments:
  f                     MS Annika crosslink result files in Microsoft Excel
                        format (.xlsx) to process. Decoys should be excluded!
required arguments:
  -fasta FASTAFILE, --fasta FASTAFILE
                        Fasta file used for crosslink search. Must contain
                        proteins identified in the MS Annika result files.
optional arguments:
  -ignore ACCESSION, --ignore ACCESSION
                        Protein accessions to be ignored. Crosslinks that only
                        link between ignored proteins will not be exported.
                        Supports input of multiple accessions.
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Prefix of the output files.
  --version             show program's version number and exit
"""

# Exporter class with constructor that takes one or multiple excel Crosslink
# result files as input, plus the fasta database and potentiol proteins to
# ignore
class MSAnnika_Exporter:

    def __init__(self, input_files: List[str], fasta_file: str, ignore_list: List[str]):
        self.input_files = input_files
        self.ignore_list = ignore_list

        # read fasta file
        sequences = dict()

        for entry in SeqIO.parse(fasta_file, "fasta"):
            identifier = str(entry.id).split("|")[1].strip()
            sequence = str(entry.seq)
            if identifier not in sequences:
                sequences[identifier] = {"sequence": sequence}
            else:
                print("WARNING: Identifier " + identifier + " is not unique!")

        self.database = sequences

    def __clean_sequence(self, sequence: str) -> str:
        return sequence.replace("[", "").replace("]", "").strip().upper()

    def __get_xl_position(self, sequence: str) -> int:
        pos = 1
        for AA in sequence:
            if AA == "[":
                return pos
            else:
                pos += 1

    def __get_proteins(self, accessions: str) -> str:
        accessions_list = accessions.split(";")
        accessions_output = ""
        for accession in accessions_list:
            if accession not in self.ignore_list:
                accessions_output = accessions_output + accession + ";"
        return accessions_output.rstrip(";")

    def __get_peptide_positions(self, sequence: str, identifiers_str: str) -> str:
        identifiers = identifiers_str.split(";")
        peptide_positions = ""
        for identifier in identifiers:
            protein = self.database[identifier]["sequence"]
            pos_in_protein = protein.find(self.__clean_sequence(sequence))
            peptide_positions = peptide_positions + str(pos_in_protein + 1) + ";"
        return peptide_positions.rstrip(";")

    def __generate_csv_df(self) -> pd.DataFrame:

        # columns
        Protein1 = []
        PepPos1 = []
        PepSeq1 = []
        LinkPos1 = []
        Protein2 = []
        PepPos2 = []
        PepSeq2 = []
        LinkPos2 = []
        Score = []
        Id = []

        xl_id = 1

        for input_file in self.input_files:
            df = pd.read_excel(input_file)

            for i, row in df.iterrows():
                if self.__get_proteins(row["Accession A"]) != "" and self.__get_proteins(row["Accession B"]) != "":
                    Protein1.append(self.__get_proteins(row["Accession A"]))
                    PepPos1.append(self.__get_peptide_positions(row["Sequence A"], self.__get_proteins(row["Accession A"])))
                    PepSeq1.append(self.__clean_sequence(row["Sequence A"]))
                    LinkPos1.append(self.__get_xl_position(row["Sequence A"]))
                    Protein2.append(self.__get_proteins(row["Accession B"]))
                    PepPos2.append(self.__get_peptide_positions(row["Sequence B"], self.__get_proteins(row["Accession B"])))
                    PepSeq2.append(self.__clean_sequence(row["Sequence B"]))
                    LinkPos2.append(self.__get_xl_position(row["Sequence B"]))
                    Score.append(row["Best CSM Score"])
                    Id.append(xl_id)
                    xl_id += 1

        result = pd.DataFrame({"Protein1": Protein1, "PepPos1": PepPos1, "PepSeq1": PepSeq1, "LinkPos1": LinkPos1,
                               "Protein2": Protein2, "PepPos2": PepPos2, "PepSeq2": PepSeq2, "LinkPos2": LinkPos2,
                               "Score": Score, "Id": Id})

        return result

    # export function, takes one argument "output_file" which sets the prefix
    # of generated output files
    def export(self, output_file = None, format = "xiNET") -> None:
        csv = self.__generate_csv_df()

        if output_file == None:
            output_file = self.input_files[0].split(".")[0]

        csv.to_csv(output_file + ".csv", index = False)

# initialize exporter and export xiNET files
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(metavar = "f",
                        dest = "files",
                        help = "Name/Path of the MS Annika result files to process.",
                        type = str,
                        nargs = "+")
    parser.add_argument("-fasta", "--fasta",
                        dest = "fasta",
                        default = None,
                        help = "Name/Path of the fasta file.",
                        required = True,
                        type = str)
    parser.add_argument("-ignore", "--ignore",
                        dest = "ignore",
                        default = [],
                        help = "Protein identifiers to ignore.",
                        type = str,
                        nargs = "*")
    parser.add_argument("-o", "--output",
                        dest = "output",
                        default = None,
                        help = "Name of the output files.",
                        type = str)
    parser.add_argument("--version",
                        action = "version",
                        version = __version)
    args = parser.parse_args()

    exporter = MSAnnika_Exporter(args.files, args.fasta, args.ignore)

    exporter.export(args.output)

if __name__ == "__main__":
    main()
