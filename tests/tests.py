#!/usr/bin/env python3

# FRAGMENT INTESITIES MS ANNIKA - TESTS
# 2022 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

def test_intensities():

    from get_intensities import main

    csms = main()

    assert csms.iloc[0]["Fragment Intensities A (Sum)"] > 30363 and csms.iloc[0]["Fragment Intensities A (Sum)"] < 30365
