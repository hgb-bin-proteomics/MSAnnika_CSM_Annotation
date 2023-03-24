# MSAnnika CSM Annotation

Calculates fragment intensities and doublet intensities of cross-linked peptides from MS Annika CSMs.

## Usage

- Install requirements: `pip install -r requirements.txt`
- Export MS Annika CSMs from Proteome Discoverer to Microsoft Excel format. Filter out decoys beforehand. Optionally also filter out low confidence CSMs.
- [Optional] Export MS Annika Crosslink Doublets to text (tab seperated) if you want to annotate doublet peaks.
- Convert any RAW files to *.mgf format.
- Set your desired parameters in `get_intensities.py`.
- Run `python get_intensities.py`.
- If the script successfully finishes, there should be a new Excel file including intensities and annotations.

## Parameters

The following parameters should be set inside the script:

```python
##### PARAMETERS #####

SPECTRA_FILE = "20220215_Eclipse_LC6_PepMap50cm-cartridge_mainlib_DSSO_3CV_stepHCD_OT_001.mgf"
CSMS_FILE = "CSMs_unfiltered.xlsx"
DOUBLETS_FILE = "CrosslinkDoublets.txt"
MODIFICATIONS = \
    {"Oxidation": [15.994915],
     "Carbamidomethyl": [57.021464],
     "DSSO": [54.01056,  85.98264, 103.99320]}
ION_TYPES = ("b", "y")
MAX_CHARGE = 4
MATCH_TOLERANCE = 0.02
DOUBLET_IDENTIFICATION_MODE = "All"
# DOUBLET_IDENTIFICATION_MODE = "Evidence"
# DOUBLET_IDENTIFICATION_MODE = "Indication"

######################
```

The `DOUBLETS_FILE` parameter is optional, if a file is provided the doublet intensities will also be calculated, otherwise leave the field blank `DOUBLETS_FILE = None` and doublets will not be annotated and no doublet intensities will be calculated. `DOUBLET_IDENTIFICATION_MODE` can either be `"All"` which means all doublets will be considered or `"Evidence"` which means only doublets identified by MS Annika in evidence mode will be considered, or `"Indication"` which means only doublets identified by MS Annika in indication mode will be considered.

## License

- [MIT](https://github.com/hgb-bin-proteomics/MSAnnika_CSM_Annotation/blob/master/LICENSE)

## Contact

- [micha.birklbauer@fh-hagenberg.at](mailto:micha.birklbauer@fh-hagenberg.at)
