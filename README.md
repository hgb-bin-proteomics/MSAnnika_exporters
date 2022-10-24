# MSAnnika_exporters

Export MS Annika crosslink results to different formats.

## Requirements

Python 3.7+ installation with pandas.

## Export to xiNET

Example usage:

```
python xiNetExporter_msannika.py "202001216_nsp8_trypsin_XL_REP1.xlsx" "202001216_nsp8_trypsin_XL_REP2.xlsx" "202001216_nsp8_trypsin_XL_REP3.xlsx" --fasta SARS-COV-2.fasta -o test --ignore P0DTC1 P0DTD1 P0DTC2
```

## Export to xiVIEW

Example usage:

```
python xiViewExporter_msannika.py "202001216_nsp8_trypsin_XL_REP1.xlsx" "202001216_nsp8_trypsin_XL_REP2.xlsx" "202001216_nsp8_trypsin_XL_REP3.xlsx" --fasta SARS-COV-2.fasta -o test --ignore P0DTC1 P0DTD1 P0DTC2
```

## Using exporters in Proteome Discoverer

To use the xiNET and xiVIEW exporters in Proteome Discoverer, add a "Scripting Node" from the "Post-Processing" tab in the "Workflow Nodes" window to your consensus workflow. You need to specify the following parameters in the Scripting Node:
- Path to Executable: Path of the python installation e.g. `C:\Users\Username\AppData\Local\Programs\Python\Python37\python.exe`
- Command Line Arguments: Path of the exporter script from `scripting_nodes` and `%NODEARGS%` e.g. `C:\Users\Username\Documents\PDScriptingNodes\xiViewExporter_msannikaPD.py %NODEARGS%`
- Requested Tables and Columns: Copy and paste the contents of `pd_tables.txt`
Re-running the consensus worklflow should create the xiNET/xiVIEW files in the study directory.

## Known Issues

[List of known issues](https://github.com/hgb-bin-proteomics/MSAnnika_exporters/issues)
