#!/usr/bin/env python3

# FRAGMENT INTESITIES MS ANNIKA
# 2023 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

# REQUIREMENTS
# pip install pandas
# pip install openpyxl
# pip install tqdm
# pip install pyteomics

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

# import packages
import pandas as pd
from tqdm import tqdm
from pyteomics import mgf, mass

from typing import Dict
from typing import List
from typing import Tuple
from typing import Set
import warnings

# reading spectra
def read_spectra(filename: str) -> Dict[int, Dict]:
    """
    Returns a dictionary that maps scan numbers to spectra:
    Dict[int -> Dict["precursor" -> float
                     "charge" -> int
                     "peaks" -> Dict[m/z -> intensity]]
    """

    result_dict = dict()

    with mgf.read(filename) as reader:
        for spectrum in reader:
            scan_nr = int(spectrum["params"]["title"].split("scan=")[1].strip("\""))
            spectrum_dict = dict()
            spectrum_dict["precursor"] = spectrum["params"]["pepmass"]
            spectrum_dict["charge"] = spectrum["params"]["charge"]
            peaks = dict()
            for i, mz in enumerate(spectrum["m/z array"]):
                peaks[mz] = spectrum["intensity array"][i]
            spectrum_dict["peaks"] = peaks
            result_dict[scan_nr] = spectrum_dict
        reader.close()

    return result_dict

# generate a position to modification mass mapping
def generate_modifications_dict(peptide: str, modification_str: str) -> Dict[int, List[float]]:
    """
    Returns a mapping of peptide positions (0 based) to possible modification masses.
    modification_str is the modification string as returned by MS Annika e.g. K5(DSSO);M7(Oxidation)
    the modification in braces has to be defined in MODIFICATIONS
    """

    modifications_dict = dict()

    modifications = modification_str.split(";")
    for modification in modifications:
        # remove possible white spaces
        modification = modification.strip()
        # get modified amino acid and modification position
        aa_and_pos = modification.split("(")[0]
        # get modification type
        mod = modification.split("(")[1].rstrip(")")

        if aa_and_pos == "Nterm":
            pos = -1
        elif aa_and_pos == "Cterm":
            pos = len(peptide)
        else:
            pos = int(aa_and_pos[1:]) - 1
        if mod in MODIFICATIONS:
            modifications_dict[pos] = MODIFICATIONS[mod]
        else:
            warnings.warn("Modification '" + mod + "' not found!")

    return modifications_dict

# generate all theoretical fragments
# adapted from https://pyteomics.readthedocs.io/en/latest/examples/example_msms.html
def generate_theoretical_fragments(peptide: str, modifications: Dict[int, List[float]], ion_types: Tuple[str] = ("b", "y"), max_charge: int = 1) -> Dict[float, str]:
    """
    Generates a set of theoretical fragment ion masses of the specified peptide with the modifications.
    """

    fragments = dict()

    for i in range(1, len(peptide)):
        for ion_type in ion_types:
            for charge in range(1, max_charge + 1):
                if ion_type[0] in "abc":
                    frag_mass = mass.fast_mass(peptide[:i], ion_type = ion_type, charge = charge)
                    mass_possibilites = set()
                    for mod_pos in modifications.keys():
                        # if the modification is within the fragment:
                        if mod_pos < i:
                            # add the modification mass / charge if its a normal modification
                            if len(modifications[mod_pos]) == 1:
                                frag_mass += modifications[mod_pos][0] / charge
                            else:
                                # if it's a crosslinking modification we add the crosslinker fragment masses
                                # to a set of possible modification mass additions to generate a fragment mass
                                # for every crosslinker fragment
                                for modification in modifications[mod_pos]:
                                    mass_possibilites.add(modification / charge)
                    # we add all possible fragment masses including all crosslinker fragment possibilites
                    if len(mass_possibilites) == 0:
                        if frag_mass not in fragments:
                            fragments[frag_mass] = ion_type + str(i) + "+" + str(charge) + ": " + peptide[:i]
                    else:
                        for mass_possibility in mass_possibilites:
                            if frag_mass + mass_possibility not in fragments:
                                fragments[frag_mass + mass_possibility] = ion_type + str(i) + "+" + str(charge) + ": " + peptide[:i]
                else:
                    frag_mass = mass.fast_mass(peptide[i:], ion_type = ion_type, charge = charge)
                    mass_possibilites = set()
                    for mod_pos in modifications.keys():
                        if mod_pos >= i:
                            if len(modifications[mod_pos]) == 1:
                                frag_mass += modifications[mod_pos][0] / charge
                            else:
                                for modification in modifications[mod_pos]:
                                    mass_possibilites.add(modification / charge)
                    if len(mass_possibilites) == 0:
                        if frag_mass not in fragments:
                            fragments[frag_mass] = ion_type + str(len(peptide) - i) + "+" + str(charge) + ": " + peptide[i:]
                    else:
                        for mass_possibility in mass_possibilites:
                            if frag_mass + mass_possibility not in fragments:
                                fragments[frag_mass + mass_possibility] = ion_type + str(len(peptide) - i) + "+" + str(charge) + ": " + peptide[i:]

    return fragments

def get_intensities(row: pd.Series, alpha: bool, spectra: Dict[int, Dict]) -> Tuple[float, Dict[float, str], Dict[float, str]]:
    """
    Returns the sum of intensities of the matched fragment ions of an identified peptide, the matched ions and all potential theoretical ions.
    """

    scan_nr = row["First Scan"]

    if alpha:
        sequence = row["Sequence A"]
        modifications = row["Modifications A"]
    else:
        sequence = row["Sequence B"]
        modifications = row["Modifications B"]

    spectrum = spectra[scan_nr]

    modifications_processed = generate_modifications_dict(sequence, modifications)
    theoretical_fragments = generate_theoretical_fragments(sequence, modifications_processed, ion_types = ION_TYPES, max_charge = MAX_CHARGE)

    total_intensity = 0
    matched_fragments = dict()

    for peak_mz in spectrum["peaks"].keys():
        for fragment in theoretical_fragments.keys():
            if round(peak_mz, 4) < round(fragment + MATCH_TOLERANCE, 4) and round(peak_mz, 4) > round(fragment - MATCH_TOLERANCE, 4):
                total_intensity += spectrum["peaks"][peak_mz]
                matched_fragments[peak_mz] = theoretical_fragments[fragment]
                break

    return total_intensity, matched_fragments, theoretical_fragments

if __name__ == "__main__":

    print("INFO: Running CSMS annotation with input files:\nSpectra: " + SPECTRA_FILE + "\nCSMs: " + CSMS_FILE)
    print("INFO: Using the following modifications:")
    print(MODIFICATIONS)
    print("INFO: Using the following ion types:")
    print(ION_TYPES)
    print("INFO: Using the following charge states:")
    print([i for i in range(1, MAX_CHARGE + 1)])
    print("INFO: Using a match tolerance of: " + str(MATCH_TOLERANCE) + " Da")
    print("INFO: Starting annotation process...")

    print("INFO: Reading spectra...")
    spectra = read_spectra(SPECTRA_FILE)
    print("INFO: Done reading spectra!")

    print("INFO: Reading CSMs...")
    csms = pd.read_excel(CSMS_FILE)
    print("INFO: Done reading CSMs!")

    tqdm.pandas(desc = "INFO: Progress bar - Alpha Peptides:")
    csms["Fragment Intensities A"], csms["Matched Ions A"], csms["Theoretical Ions A"] = csms.progress_apply(lambda row: get_intensities(row, True, spectra), axis = 1)
    print("INFO: Done processing Alpha Peptides!")

    tqdm.pandas(desc = "INFO: Progress bar - Beta Peptides:")
    csms["Fragment Intensities B"], csms["Matched Ions B"], csms["Theoretical Ions B"] = csms.progress_apply(lambda row: get_intensities(row, False, spectra), axis = 1)
    print("INFO: Done processing Beta Peptides!")

    csms["Fragment Intensities Total"] = csms.apply(lambda row: row["Fragment Intensities A"] + row["Fragment Intensities B"], axis = 1)

    csms.to_excel(".".join(CSMS_FILE.split(".")[:-1]) + "_with_intensities.xlsx")
    print("SUCCESS: Output file generated as '" + ".".join(CSMS_FILE.split(".")[:-1]) + "_with_intensities.xlsx" + "'!")
