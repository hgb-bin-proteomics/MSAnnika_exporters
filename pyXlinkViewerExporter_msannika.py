#!/usr/bin/env python3

# Exporter of MS Annika Crosslink Results to PyXlinkViewer for pyMOL
# 2022 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import re
import argparse
import pandas as pd
from typing import List
from typing import Union

# needs biopython: pip install biopython
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.Align import substitution_matrices

__version = "1.0.0"
__date = "20221028"

"""
DESCRIPTION:
A script to export MS Annika results to PyXlinkViewer format for visualizing
crosslinks in pyMOL.
USAGE:
pyXlinkViewerExporter_msannika.py f [f ...]
                                    [-pdb PDB_FILE]
                                    [-go GAP_OPEN_PENALTY]
                                    [-ge GAP_EXTENSION_PENALTY]
                                    [-si SEQUENCE_IDENTITY]
                                    [-tp TRUST_PDB]
                                    [-o OUTPUT]
                                    [-h]
                                    [--version]
positional arguments:
  f                     MS Annika crosslink result files in Microsoft Excel
                        format (.xlsx) to process.
required arguments:
  -pdb PDB_FILE, --pdb PDB_FILE
                        PDB file of the structure that crosslinks should be
                        exported to.
optional arguments:
  -go GAP_OPEN_PENALTY, --gap_open GAP_OPEN_PENALTY
                        Gap open penalty for sequence alignment.
                        Default: -10
  -ge GAP_EXTENSION_PENALTY, --gap_extension GAP_EXTENSION_PENALTY
                        Gap extension penalty for sequence alignment.
                        Default: -1
  -si SEQUENCE_IDENTITY, --sequence_identity SEQUENCE_IDENTITY
                        Sequence identity threshold in percent to consider two
                        aligned sequences as matching.
                        Default: 80
  -tp TRUST_PDB, --trust_pdb TRUST_PDB
                        Skip crosslinks that don't have crosslinks site in the
                        PDB sequence.
                        Default: True
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Prefix of the output files.
  --version             show program's version number and exit
"""

# Amino acid letter code translations
AA_translate = {"GLY": "G", "PRO": "P",
                "ALA": "A", "VAL": "V",
                "LEU": "L", "ILE": "I",
                "MET": "M", "CYS": "C",
                "PHE": "F", "TYR": "Y",
                "TRP": "W", "HIS": "H",
                "LYS": "K", "ARG": "R",
                "GLN": "Q", "ASN": "N",
                "GLU": "E", "ASP": "D",
                "SER": "S", "THR": "T"}

# BLOSUM62
Blosum62 = substitution_matrices.load("BLOSUM62")

# Exporter class with constructor that takes one or multiple excel Crosslink
# result files as input, the pdb structure of the crosslinked protein, gap
# penalties for alignment, a sequence identity threshold for matching
# crosslinked peptides to the PDB structure and specification if crosslinks
# that link to non crosslink reactive sites should be discarded
class MSAnnika_Exporter:

    def __init__(self, input_files: List[str], pdb_file: str,
                 gap_open: Union[int, float], gap_extension: Union[int, float],
                 sequence_identity: Union[int, float], trust_pdb: bool):

        pdb_data = self.__get_pdb_data(pdb_file)

        self.input_files = input_files
        self.sequence = pdb_data["sequence"]
        self.chains = pdb_data["chains"]
        self.residue_numbers = pdb_data["residue_numbers"]

        # alignment parameters
        self.gap_open = float(gap_open)
        self.gap_extension = float(gap_extension)
        self.sequence_identity = float(sequence_identity / 100)
        self.trust_pdb = bool(trust_pdb)

    # pdb file parser
    def __get_pdb_data(self, pdb_file: str) -> dict:

        sequence = []
        chains = []
        residue_numbers = dict()
        residue_numbers_lst = []

        with open(pdb_file, "r", encoding = "utf-8") as f:
            pdb_data = f.readlines()
            f.close()

        for line in pdb_data:
            if line.split()[0].strip().upper() == "ATOM":
                parts = line.split()
                if parts[3].strip() in AA_translate:
                    residue = AA_translate[parts[3].strip()]
                    chain = parts[4].strip()
                    residue_number = parts[5].strip()
                    if chain in residue_numbers:
                        if residue_number not in residue_numbers[chain]:
                            residue_numbers[chain].append(residue_number)
                            sequence.append(residue)
                            chains.append(chain)
                    else:
                        residue_numbers[chain] = [residue_number]
                        sequence.append(residue)
                        chains.append(chain)
                else:
                    print("WARNING: ", parts[3].strip(), " is not a supported amino acid.")

        for chain in sorted(residue_numbers.keys()):
            residue_numbers_lst = residue_numbers_lst + residue_numbers[chain]

        return {"sequence": "".join(sequence), "chains": "".join(chains), "residue_numbers": residue_numbers_lst}

    def __clean_sequence(self, sequence: str) -> str:
        return sequence.replace("[", "").replace("]", "").strip().upper()

    def __get_xl_position(self, sequence: str) -> int:
        pos = 0
        for AA in sequence:
            if AA == "[":
                return pos
            else:
                pos += 1

    def __get_xl_position_and_chain_in_protein(self, peptide_sequence: str) -> List[str]:
        pep_seq = self.__clean_sequence(peptide_sequence)
        pep_pos_in_proteins = [m.start() for m in re.finditer(pep_seq, self.sequence)]
        xl_pos_in_pep = self.__get_xl_position(peptide_sequence)
        if len(pep_pos_in_proteins) == 0:
            alignments = sorted(pairwise2.align.localds(Seq(self.sequence), Seq(pep_seq),
                                                        Blosum62, self.gap_open, self.gap_extension),
                                key = lambda alignment: alignment.score, reverse = True)
            if len(alignments) == 0:
                return []
            else:
                top_alignment = alignments[0]
                seqA = top_alignment.seqA[top_alignment.start:top_alignment.end]
                seqB = top_alignment.seqB[top_alignment.start:top_alignment.end]
                sequence_identity = self.__calculate_sequence_identity(seqA, seqB)
                if sequence_identity > self.sequence_identity:
                    xl_pos_in_alignment = xl_pos_in_pep
                    if len(pep_seq) != len(seqB):
                        xl_pos_in_alignment = self.__calculate_shifted_xl_pos(seqB, xl_pos_in_pep)
                    if self.trust_pdb:
                        if seqA[xl_pos_in_alignment] == seqB[xl_pos_in_alignment]:
                            pep_pos_in_protein = self.sequence.find(seqA)
                            xl_position = pep_pos_in_protein + xl_pos_in_alignment
                            xl_chain = self.chains[xl_position]
                            xl_residue = self.residue_numbers[xl_position]
                            return [str(xl_residue) + "|" + str(xl_chain) + "|"]
                        else:
                            return []
                    else:
                        pep_pos_in_protein = self.sequence.find(seqA)
                        xl_position = pep_pos_in_protein + xl_pos_in_alignment
                        xl_chain = self.chains[xl_position]
                        xl_residue = self.residue_numbers[xl_position]
                        return [str(xl_residue) + "|" + str(xl_chain) + "|"]
                else:
                    return []
        else:
            chain_residues = []
            for pep_pos_in_protein in pep_pos_in_proteins:
                xl_position = pep_pos_in_protein + xl_pos_in_pep
                xl_chain = self.chains[xl_position]
                xl_residue = self.residue_numbers[xl_position]
                chain_residues.append(str(xl_residue) + "|" + str(xl_chain) + "|")
            return chain_residues

    def __calculate_sequence_identity(self, seqA: str, seqB: str) -> float:
        ident = 0
        for i in range(len(seqA)):
            if seqA[i] == seqB[i]:
                ident += 1
        return float(ident/len(seqA))

    def __calculate_shifted_xl_pos(self, alignment: str, xl_pos_in_pep: int) -> int:
        new_xl_pos = xl_pos_in_pep
        if "-" not in alignment[:xl_pos_in_pep + 1]:
            return xl_pos_in_pep
        else:
            gaps = [m.start() for m in re.finditer("-", alignment)]
            curr_limit = xl_pos_in_pep + 1
            for gap in gaps:
                if gap < curr_limit:
                    curr_limit += 1
                    new_xl_pos += 1
            return new_xl_pos

    def __generate_output_string(self) -> str:

        if len(self.sequence) == len(self.chains) == len(self.residue_numbers):
            pass
        else:
            print(len(self.sequence))
            print(len(self.chains))
            print(len(self.residue_numbers))

            print(self.sequence)
            print(self.chains)
            print(self.residue_numbers)

            raise Exception("ERROR: Sequence, Chain and Residue Numbers are not matching! Exiting!")

        output_string = ""

        for input_file in self.input_files:
            df = pd.read_excel(input_file)

            for i, row in df.iterrows():
                links_a = self.__get_xl_position_and_chain_in_protein(row["Sequence A"])
                links_b = self.__get_xl_position_and_chain_in_protein(row["Sequence B"])
                if len(links_a) != 0 and len(links_b) != 0:
                    for link_a in links_a:
                        for link_b in links_b:
                            output_string = output_string + link_a + link_b + "\n"

        return output_string

    # export function, takes one argument "output_file" which sets the prefix
    # of generated output files
    def export(self, output_file = None, format = "PyXlinkViewer") -> None:
        output_string = self.__generate_output_string()

        if output_file == None:
            output_file = self.input_files[0].split(".")[0]

        with open(output_file + "_crosslinks.txt", "w", encoding = "utf-8") as f:
            f.write(output_string)
            f.close()

        parsed_pdb = ""
        for i, r in enumerate(self.residue_numbers):
            parsed_pdb = parsed_pdb + self.sequence[i] + " " + self.chains[i] + " " + r + "\n"

        with open(output_file + "_parsedPDB.txt", "w", encoding = "utf-8") as f:
            f.write(parsed_pdb)
            f.close()

        with open(output_file + "_sequence.fasta", "w", encoding = "utf-8") as f:
            f.write(">" + output_file + "\n" + self.sequence)
            f.close()

# initialize exporter and export PyXlinkViewer files
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(metavar = "f",
                        dest = "files",
                        help = "Name/Path of the MS Annika result files to process.",
                        type = str,
                        nargs = "+")
    parser.add_argument("-pdb", "--pdb",
                        dest = "pdb",
                        default = None,
                        help = "Name/Path of the pdb file.",
                        required = True,
                        type = str)
    parser.add_argument("-go", "--gap_open",
                        dest = "gap_open",
                        default = -10.0,
                        help = "Gap open penalty for alignment.",
                        type = float)
    parser.add_argument("-ge", "--gap_extension",
                        dest = "gap_extension",
                        default = -1.0,
                        help = "Gap extension penalty for alignment.",
                        type = float)
    parser.add_argument("-si", "--sequence_identity",
                        dest = "sequence_identity",
                        default = 80.0,
                        help = "Sequence identity threshold needed to match two sequences.",
                        type = float)
    parser.add_argument("-tp", "--trust_pdb",
                        dest = "trust_pdb",
                        default = True,
                        help = "Reject crosslinks that do not map to a reactive site on the PDB structure.",
                        type = bool)
    parser.add_argument("-o", "--output",
                        dest = "output",
                        default = None,
                        help = "Name of the output files.",
                        type = str)
    parser.add_argument("--version",
                        action = "version",
                        version = __version)
    args = parser.parse_args()

    exporter = MSAnnika_Exporter(args.files, args.pdb,
                                 args.gap_open, args.gap_extension,
                                 args.sequence_identity, args.trust_pdb)

    exporter.export(args.output)

if __name__ == "__main__":
    main()