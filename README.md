# MS Annika Exporters

Export [MS Annika](https://ms.imp.ac.at/index.php?action=ms-annika) crosslink results to different file formats.

## Requirements

Python 3.7+ installation with pandas, openpyxl and biopython to run the scripts or to use the Proteome Discoverer Scripting Nodes.
- Install [pandas](https://pandas.pydata.org/): `pip install pandas`
- Install [openpyxl](https://openpyxl.readthedocs.io/en/stable/): `pip install openpyxl`
- Install [biopython](https://biopython.org/): `pip install biopython`

Alternatively there are Windows binaries available in the [Releases](https://github.com/hgb-bin-proteomics/MSAnnika_exporters/releases) tab that don't require a python installation.

## Export to [xiNET](https://crosslinkviewer.org/)

```
EXPORTER DESCRIPTION:
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
                        format (.xlsx) to process.
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
```

Example usage:

```
python xiNetExporter_msannika.py "202001216_nsp8_trypsin_XL_REP1.xlsx" "202001216_nsp8_trypsin_XL_REP2.xlsx" "202001216_nsp8_trypsin_XL_REP3.xlsx" --fasta SARS-COV-2.fasta -o test --ignore P0DTC1 P0DTD1 P0DTC2
```

Or using the Windows binary:

```
xiNetExporter_msannika.exe "202001216_nsp8_trypsin_XL_REP1.xlsx" "202001216_nsp8_trypsin_XL_REP2.xlsx" "202001216_nsp8_trypsin_XL_REP3.xlsx" --fasta SARS-COV-2.fasta -o test --ignore P0DTC1 P0DTD1 P0DTC2
```

## Export to [xiVIEW](https://xiview.org/xiNET_website/index.php)

```
EXPORTER DESCRIPTION:
A script to export MS Annika results to xiVIEW input files (CSV + FASTA).
USAGE:
xiViewExporter_msannika.py f [f ...]
                             [-fasta FASTA]
                             [-ignore IGNORE]
                             [-o OUTPUT]
                             [-h]
                             [--version]
positional arguments:
  f                     MS Annika crosslink result files in Microsoft Excel
                        format (.xlsx) to process.
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
```

Example usage:

```
python xiViewExporter_msannika.py "202001216_nsp8_trypsin_XL_REP1.xlsx" "202001216_nsp8_trypsin_XL_REP2.xlsx" "202001216_nsp8_trypsin_XL_REP3.xlsx" --fasta SARS-COV-2.fasta -o test --ignore P0DTC1 P0DTD1 P0DTC2
```

Or using the Windows binary:

```
xiViewExporter_msannika.exe "202001216_nsp8_trypsin_XL_REP1.xlsx" "202001216_nsp8_trypsin_XL_REP2.xlsx" "202001216_nsp8_trypsin_XL_REP3.xlsx" --fasta SARS-COV-2.fasta -o test --ignore P0DTC1 P0DTD1 P0DTC2
```

## Export to [PyXlinkViewer for pyMOL](https://github.com/BobSchiffrin/PyXlinkViewer)

```
EXPORTER DESCRIPTION:
A script to export MS Annika results to PyXlinkViewer format for visualizing
crosslinks in pyMOL.
USAGE:
pyXlinkViewerExporter_msannika.py f [f ...]
                                    [-pdb PDB_FILE]
                                    [-go GAP_OPEN_PENALTY]
                                    [-ge GAP_EXTENSION_PENALTY]
                                    [-si SEQUENCE_IDENTITY]
                                    [-allowxlmismatch]
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
  -allowxlmismatch, --allowxlmismatch
                        Flag to report crosslinks that don't link to a crosslink
                        site in the PDB sequence.
                        Default: Do not report such crosslinks.
  -ic, --ignore_chains
                        Ignore specific chains in the PDB file.
                        Default: No chains are ignored.
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Prefix of the output files.
  --version             show program's version number and exit
```

Example usage:

```
python pyXlinkViewerExporter_msannika.py "202001216_nsp8_trypsin_XL_REP1.xlsx" "202001216_nsp8_trypsin_XL_REP2.xlsx" "202001216_nsp8_trypsin_XL_REP3.xlsx" --pdb 6yhu.pdb -o test
```

Or using the Windows binary:

```
pyXlinkViewerExporter_msannika.exe "202001216_nsp8_trypsin_XL_REP1.xlsx" "202001216_nsp8_trypsin_XL_REP2.xlsx" "202001216_nsp8_trypsin_XL_REP3.xlsx" --pdb 6yhu.pdb -o test
```

## Using exporters in Proteome Discoverer

To use the xiNET and xiVIEW exporters in Proteome Discoverer, add a "Scripting Node" from the "Post-Processing" tab in the "Workflow Nodes" window to your consensus workflow. You need to specify the following parameters in the Scripting Node:

- Path to Executable: Path of the python installation e.g. `C:\Users\Username\AppData\Local\Programs\Python\Python37\python.exe`
- Command Line Arguments: Path of the exporter script from `scripting_nodes` and `%NODEARGS%` e.g. `C:\Users\Username\Documents\PDScriptingNodes\xiViewExporter_msannikaPD.py %NODEARGS%`
- Requested Tables and Columns: Copy and paste the contents of `pd_tables.txt`

Re-running the consensus worklflow should create the xiNET/xiVIEW files in the study directory.

## Known Issues

[List of known issues](https://github.com/hgb-bin-proteomics/MSAnnika_exporters/issues)

## Citing

If you are using MS Annika please cite:
```
MS Annika: A New Cross-Linking Search Engine
Georg J. Pirklbauer, Christian E. Stieger, Manuel Matzinger, Stephan Winkler, Karl Mechtler, and Viktoria Dorfer
Journal of Proteome Research 2021 20 (5), 2560-2569
DOI: 10.1021/acs.jproteome.0c01000
```

## License

- [MIT](https://github.com/hgb-bin-proteomics/MSAnnika_exporters/blob/master/LICENSE)

## Contact

[micha.birklbauer@fh-hagenberg.at](mailto:micha.birklbauer@fh-hagenberg.at)
