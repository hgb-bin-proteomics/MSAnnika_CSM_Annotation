#!/usr/bin/env python3

# FRAGMENT INTESITIES MS ANNIKA - TESTS
# 2022 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

def test_alpha_fragment_intensities():

    from get_intensities import main

    csms = main()

    assert csms.iloc[0]["Fragment Intensities A (Sum)"] > 30363 and csms.iloc[0]["Fragment Intensities A (Sum)"] < 30365

def test_beta_fragment_intensities():

    from get_intensities import main

    csms = main()

    assert csms.iloc[0]["Fragment Intensities B (Sum)"] > 30363 and csms.iloc[0]["Fragment Intensities B (Sum)"] < 30365

def test_alpha_doublet_intensities():

    from get_intensities import main

    csms = main()

    assert csms.iloc[0]["Alpha Doublet Intensities Total"] > 13651 and csms.iloc[0]["Alpha Doublet Intensities Total"] < 13653

def test_beta_doublet_intensities():

    from get_intensities import main

    csms = main()

    assert csms.iloc[0]["Beta Doublet Intensities Total"] > 13651 and csms.iloc[0]["Beta Doublet Intensities Total"] < 13653

def test_alpha_doublet_peaks():

    from get_intensities import main

    csms = main()

    assert csms.iloc[0]["Alpha Light Peaks"][0][0] > 403 and csms.iloc[0]["Alpha Light Peaks"][0][0] < 405 and \
           csms.iloc[0]["Alpha Heavy Peaks"][0][0] > 419 and csms.iloc[0]["Alpha Heavy Peaks"][0][0] < 421

def test_beta_doublet_peaks():

    from get_intensities import main

    csms = main()

    assert csms.iloc[0]["Beta Light Peaks"][0][0] > 403 and csms.iloc[0]["Beta Light Peaks"][0][0] < 405 and \
           csms.iloc[0]["Beta Heavy Peaks"][0][0] > 419 and csms.iloc[0]["Beta Heavy Peaks"][0][0] < 421
