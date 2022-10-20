#!/usr/bin/env python3

# Exporter of MS Annika Crosslink Results to xiNET Format - PD Scripting Node
# 2022 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import sys
import json
import pandas as pd
from typing import List

__version = "1.0.0"
__date = "20221018"

"""
DESCRIPTION:
A Proteome Discoverer Scripting Node to export MS Annika results to xiNET
input files (CSV + FASTA).
USAGE:
xiNetExporter_msannikaPD.py NODEARGS
"""

# Exporter class with constructor that takes Crosslink pandas dataframe as
# input, plus the fasta database as string and potential proteins to
# ignore
class MSAnnika_Exporter:

    def __init__(self, input_files: pd.DataFrame, fasta_data: str, ignore_list: List[str]):
        self.input_files = input_files
        self.ignore_list = ignore_list

        # read fasta file
        sequences = dict()

        for entry in fasta_data.split(">"):
            lines = entry.split("\n")
            if lines[0].strip() != "":
                description = lines[0].strip()
                identifier = lines[0].strip().split("|")[1].strip()
                sequence = "".join([x.strip() for x in lines[1:]])
                if identifier not in sequences:
                    sequences[identifier] = {"description": description, "sequence": sequence}
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

        for df in [self.input_files]:
            for i, row in df.iterrows():
                if str(row["Decoy"]).lower() == "false" and str(row["Confidence"]) == "High":
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

    def __generate_fasta_str(self) -> str:
        fasta_str = ""
        for key in self.database:
            fasta_str = fasta_str + ">" + key + " " + self.database[key]["description"] + "\n" + self.database[key]["sequence"] + "\n"
        return fasta_str

    # export function, takes one argument "output_file" which sets the prefix
    # of generated output files
    def export(self, output_file = None, format = "xiNET") -> None:
        csv = self.__generate_csv_df()
        fasta = self.__generate_fasta_str()

        if output_file == None:
            output_file = self.input_files[0].split(".")[0]

        csv.to_csv(output_file + ".csv", index = False)
        with open(output_file + ".fasta", "w", encoding = "utf-8") as f:
            f.write(fasta)
            f.close()

# generate fasta file from Proteome Discoverer protein table
def generate_fasta_from_PD(proteins: pd.DataFrame) -> str:

    fasta_str = ""
    for i, row in proteins.iterrows():
        fasta_str = fasta_str + ">sp|" + row["Accession"] + "|" + row["Description"] + "\n"
        fasta_str = fasta_str + row["Sequence"] + "\n"

    return fasta_str

# loading data from Proteome Discoverer
def get_data_from_PD() -> dict:
    if len(sys.argv) != 2:
        raise Exception("Exactly one argument of type NODEARGS is required but was not given!")

    nodeargs = sys.argv[1]

    with open(nodeargs, "r", encoding = "utf-8") as f:
        nodeargs_data = json.load(f)
        f.close()

    result = dict()
    proteins = ""
    crosslinks = ""
    resultfile_path = nodeargs_data["ResultFilePath"]

    for table in nodeargs_data["Tables"]:
        if table["TableName"] == "Proteins":
            proteins = table["DataFile"]
        if table["TableName"] == "Crosslinks":
            crosslinks = table["DataFile"]

    if proteins == "" or crosslinks == "":
        raise Exception("Could not load data from Proteome Discoverer!")

    proteins_df = pd.read_csv(proteins, sep = "\t")
    crosslinks_df = pd.read_csv(crosslinks, sep = "\t", dtype = {"Accession A": str, "Accession B": str})

    result["files"] = crosslinks_df
    result["fasta"] = generate_fasta_from_PD(proteins_df)
    result["ignore"] = []
    result["output"] = ".".join(resultfile_path.split(".")[:-1]) + "_xiNET"

    return result

# initialize exporter and export xiNET files
def main() -> None:

    args = get_data_from_PD()

    exporter = MSAnnika_Exporter(args["files"], args["fasta"], args["ignore"])

    exporter.export(args["output"])

if __name__ == "__main__":
    main()
