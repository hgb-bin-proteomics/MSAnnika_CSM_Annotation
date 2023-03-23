#!/usr/bin/env python3

# FRAGMENT INTESITIES MS ANNIKA
# 2023 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

# version tracking
__version = "1.1.0"
__date = "2023-03-23"

# REQUIREMENTS
# pip install pandas
# pip install openpyxl
# pip install tqdm
# pip install pyteomics

##### PARAMETERS #####

SPECTRA_FILE = "20220215_Eclipse_LC6_PepMap50cm-cartridge_mainlib_DSSO_3CV_stepHCD_OT_001.mgf"
CSMS_FILE = "CSMs_unfiltered.xlsx"
DOUBLETS_FILE = "CrosslinkDoublets.txt"
MODIFICATIONS = \
    {"Oxidation": [15.994915],
     "Carbamidomethyl": [57.021464],
     "DSSO": [54.01056, 85.98264, 103.99320]}
ION_TYPES = ("b", "y")
MAX_CHARGE = 4
MATCH_TOLERANCE = 0.02
DOUBLET_IDENTIFICATION_MODE = "All"
# DOUBLET_IDENTIFICATION_MODE = "Evidence"
# DOUBLET_IDENTIFICATION_MODE = "Indication"

######################

# import packages
import os
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
                     "charge"    -> int
                     "peaks"     -> Dict[m/z -> intensity]]
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

def read_doublets(filename: str) -> Dict[int, Dict]:
    """
    Returns a dictionary that maps scan numbers to doublets:
    Dict[int -> Dict["alpha_light" -> set(float)
                     "alpha_heavy" -> set(float)
                     "beta_light"  -> set(float)
                     "beta_heavy"  -> set(float)]
    """

    proton_mass = 1.007276466812
    result_dict = dict()
    doublet_df = pd.read_csv(filename, sep = "\t")

    for i, row in tqdm(doublet_df.iterrows(), total = doublet_df.shape[0], desc = "INFO: Progress bar - Reading doublet file"):

        if DOUBLET_IDENTIFICATION_MODE != "All":
            if DOUBLET_IDENTIFICATION_MODE not in row["Identification Mode"]:
                continue

        alpha_doublets = row["Complete Alpha Doublet"]
        uncharged_aL, uncharged_aH = [float(m.strip()) for m in alpha_doublets.split("|")]
        mz_aL = row["m/z Alpha Light"]
        mz_aH = row["m/z Alpha Heavy"]
        beta_doublets = row["Complete Beta Doublet"]
        uncharged_bL, uncharged_bH = [float(m.strip()) for m in beta_doublets.split("|")]
        mz_bL = row["m/z Beta Light"]
        mz_bH = row["m/z Beta Heavy"]
        scan_nr = row["Scan number"]

        # alpha light possible m/z
        aL_masses = set()
        if mz_aL != 0:
            aL_masses.add(mz_aL)

        for charge in range(1, MAX_CHARGE + 1):
            aL_masses.add((uncharged_aL + charge * proton_mass) / charge)

        # alpha heavy possible m/z
        aH_masses = set()
        if mz_aH != 0:
            aH_masses.add(mz_aH)

        for charge in range(1, MAX_CHARGE + 1):
            aH_masses.add((uncharged_aH + charge * proton_mass) / charge)

        # beta light possible m/z
        bL_masses = set()
        if mz_bL != 0:
            bL_masses.add(mz_bL)

        for charge in range(1, MAX_CHARGE + 1):
            bL_masses.add((uncharged_bL + charge * proton_mass) / charge)

        # beta heavy possible m/z
        bH_masses = set()
        if mz_bH != 0:
            bH_masses.add(mz_bH)

        for charge in range(1, MAX_CHARGE + 1):
            bH_masses.add((uncharged_bH + charge * proton_mass) / charge)

        if scan_nr not in result_dict:
            result_dict[scan_nr] = {"alpha_light": aL_masses,
                                    "alpha_heavy": aH_masses,
                                    "beta_light": bL_masses,
                                    "beta_heavy": bH_masses}
        else:
            result_dict[scan_nr]["alpha_light"].update(aL_masses)
            result_dict[scan_nr]["alpha_heavy"].update(aH_masses)
            result_dict[scan_nr]["beta_light"].update(bL_masses)
            result_dict[scan_nr]["beta_heavy"].update(bH_masses)

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

def get_doublets(row: pd.Series, spectra: Dict[int, Dict], doublets: Dict[int, Dict]) -> Tuple[float,
                                                                                               float,
                                                                                               float,
                                                                                               float,
                                                                                               List[Tuple],
                                                                                               List[Tuple],
                                                                                               List[Tuple],
                                                                                               List[Tuple]]:
    """
    Returns the sum of intensities of alpha light doublet peaks, alpha heavy doublet peaks, beta light doublet peaks, beta heavy doublet peaks (index 0 - 3).
    Returns the peaks (Tuple (m/z, intensity)) of alpha light doublet peaks, alpha heavy doublet peaks, beta light doublet peaks, beta heavy doublet peaks (index 4 - 7).
    """

    scan_nr = row["First Scan"]
    spectrum = spectra[scan_nr]
    aL_doublets = []
    aH_doublets = []
    bL_doublets = []
    bH_doublets = []

    if scan_nr not in doublets:
        return 0, 0, 0, 0, [], [], [], []
    else:
        for peak_mz in spectrum["peaks"].keys():

            for doublet_mz in doublets[scan_nr]["alpha_light"]:
                if round(peak_mz, 4) < round(doublet_mz + MATCH_TOLERANCE, 4) and round(peak_mz, 4) > round(doublet_mz - MATCH_TOLERANCE, 4):
                    aL_doublets.append((peak_mz, spectrum["peaks"][peak_mz]))
                    break

            for doublet_mz in doublets[scan_nr]["alpha_heavy"]:
                if round(peak_mz, 4) < round(doublet_mz + MATCH_TOLERANCE, 4) and round(peak_mz, 4) > round(doublet_mz - MATCH_TOLERANCE, 4):
                    aH_doublets.append((peak_mz, spectrum["peaks"][peak_mz]))
                    break

            for doublet_mz in doublets[scan_nr]["beta_light"]:
                if round(peak_mz, 4) < round(doublet_mz + MATCH_TOLERANCE, 4) and round(peak_mz, 4) > round(doublet_mz - MATCH_TOLERANCE, 4):
                    bL_doublets.append((peak_mz, spectrum["peaks"][peak_mz]))
                    break

            for doublet_mz in doublets[scan_nr]["beta_heavy"]:
                if round(peak_mz, 4) < round(doublet_mz + MATCH_TOLERANCE, 4) and round(peak_mz, 4) > round(doublet_mz - MATCH_TOLERANCE, 4):
                    bH_doublets.append((peak_mz, spectrum["peaks"][peak_mz]))
                    break

        return sum([p[1] for p in aL_doublets]), \
               sum([p[1] for p in aH_doublets]), \
               sum([p[1] for p in bL_doublets]), \
               sum([p[1] for p in bH_doublets]), \
               aL_doublets, aH_doublets, bL_doublets, bH_doublets

def get_total_fragment_intensity(row: pd.Series) -> float:
    """
    Returns the total intensity of fragment ions in the spectrum.
    """

    if row["Sequence A"].strip() == row["Sequence B"].strip():
        return (row["Fragment Intensities A (Sum)"] + row["Fragment Intensities B (Sum)"]) / 2
    else:
        return row["Fragment Intensities A (Sum)"] + row["Fragment Intensities B (Sum)"]

def get_total_doublet_intensity(row: pd.Series) -> float:
    """
    Returns the total intensity of doublet peaks in the spectrum.
    """

    if row["Sequence A"].strip() == row["Sequence B"].strip():
        return (row["Alpha Doublet Intensities Total"] + row["Beta Doublet Intensities Total"]) / 2
    else:
        return row["Alpha Doublet Intensities Total"] + row["Beta Doublet Intensities Total"]

def get_spectrum_intensity(row: pd.Series, spectra: Dict[int, Dict]) -> float:
    """
    Returns the total intensity of all peaks in a spectrum.
    """

    scan_nr = row["First Scan"]
    spectrum = spectra[scan_nr]

    total_intensity = 0
    for peak_mz in spectrum["peaks"]:
        total_intensity += spectrum["peaks"][peak_mz]

    return total_intensity

def main() -> pd.DataFrame:

    print("INFO: Running CSM annotation with input files:\nSpectra: " + SPECTRA_FILE + "\nCSMs: " + CSMS_FILE + "\nDoublets: " + DOUBLETS_FILE)
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
    print("INFO: Done reading CSMs! Starting fragment ion annotation...")

    tqdm.pandas(desc = "INFO: Progress bar - Annotating alpha peptide fragments")
    csms[["Fragment Intensities A (Sum)", "Matched Ions A", "Theoretical Ions A"]] = \
        csms.progress_apply(lambda row: get_intensities(row, True, spectra),
                            axis = 1,
                            result_type = "expand")
    print("INFO: Done processing alpha peptides!")

    tqdm.pandas(desc = "INFO: Progress bar - Annotating beta peptide fragments")
    csms[["Fragment Intensities B (Sum)", "Matched Ions B", "Theoretical Ions B"]] = \
        csms.progress_apply(lambda row: get_intensities(row, False, spectra),
                            axis = 1,
                            result_type = "expand")
    print("INFO: Done processing beta peptides!")

    csms["Fragment Intensities Total"] = csms.apply(lambda row: get_total_fragment_intensity(row), axis = 1)

    print("INFO: Done annotating fragment ions!")

    if DOUBLETS_FILE is not None and os.path.isfile(DOUBLETS_FILE):
        print("INFO: Doublet file was provided! Reading doublet file...")
        doublets = read_doublets(DOUBLETS_FILE)
        print("INFO: Done reading doublet file! Starting doublet annotation...")
        tqdm.pandas(desc = "INFO: Progress bar - Annotating doublets")
        csms[["Alpha Light Intensities (Sum)", "Alpha Heavy Intensities (Sum)", "Beta Light Intensities (Sum)", "Beta Heavy Intensities (Sum)",
              "Alpha Light Peaks", "Alpha Heavy Peaks", "Beta Light Peaks", "Beta Heavy Peaks"]] = \
              csms.progress_apply(lambda row: get_doublets(row, spectra, doublets),
                                  axis = 1,
                                  result_type = "expand")
        print("INFO: Done annotating doublets! Calculating intensities...")
        csms["Alpha Doublet Intensities Total"] = csms.apply(lambda row: row["Alpha Light Intensities (Sum)"] + row["Alpha Heavy Intensities (Sum)"], axis = 1)
        csms["Beta Doublet Intensities Total"] = csms.apply(lambda row: row["Beta Light Intensities (Sum)"] + row["Beta Heavy Intensities (Sum)"], axis = 1)
        csms["Doublet Intensities Total"] = csms.apply(lambda row: get_total_doublet_intensity(row), axis = 1)
        print("INFO: Done calculating doublet intensities!")
    else:
        print("INFO: Doublet file was not provided or not found! Skipping doublet annotation.")

    print("INFO: Calculating total spectrum intensities...")
    tqdm.pandas(desc = "INFO: Progress bar - Intensity per spectrum")
    csms["Total Intensity in Spectrum"] = csms.progress_apply(lambda row: get_spectrum_intensity(row, spectra), axis = 1)
    print("INFO: Done calculating total spectrum intensities!")

    csms.to_excel(".".join(CSMS_FILE.split(".")[:-1]) + "_with_intensities.xlsx")
    print("SUCCESS: Output file generated as '" + ".".join(CSMS_FILE.split(".")[:-1]) + "_with_intensities.xlsx" + "'!")

    return csms

if __name__ == "__main__":

        main()
