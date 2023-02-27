# MSAnnika CSM Annotation

Calculates fragment intensities of cross-linked peptides from MS Annika CSMs.

## Usage

- Install requirements: `pip install -r requirements.txt`
- Export MS Annika CSMs from Proteome Discoverer to Microsoft Excel format. Filter out decoys beforehand. Optionally also filter out low confidence CSMs.
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
MODIFICATIONS = \
    {"Oxidation": [15.994915],
     "Carbamidomethyl": [57.021464],
     "DSSO": [54.01056,  85.98264, 103.99320]}
ION_TYPES = ("b", "y")
MAX_CHARGE = 4
MATCH_TOLERANCE = 0.02

######################
```

## License

- [MIT](https://github.com/hgb-bin-proteomics/MSAnnika_CSM_Annotation/blob/master/LICENSE)

## Contact

- [micha.birklbauer@fh-hagenberg.at](mailto:micha.birklbauer@fh-hagenberg.at)
