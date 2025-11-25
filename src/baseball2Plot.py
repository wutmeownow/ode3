#!/usr/bin/env python
## @package RKnplotDemo
# Example for making plots from our differential equation solver.
#
# Requires that baseball2 has been run first for all 4 pitches

import ROOT as r      # needed to load ROOT libraries
#import numpy as np    # only needed to use numpy arrays, see example below
import sys
from math import sqrt

# pdf name
pdf_name = "pitches.pdf"

# Create a canvas
c = r.TCanvas("c", "c", 800, 600)



# slider pitch file
tf1=r.TFile("slider.root") # open file for read access
if tf1.IsZombie():
    print("slider.root not found in current directory, run baseball2 first")
    sys.exit(1)

slider_zvsx=tf1.Get("g1")   # slider z vs x graph
slider_yvsx=tf1.Get("g2")   # slider y vs x graph

# plot slider graphs to first page
h1 = slider_zvsx.GetHistogram()
h1.GetXaxis().SetRangeUser(0, 60.)
h1.GetYaxis().SetRangeUser(-4, 2)
slider_zvsx.Draw("AC")
slider_yvsx.SetLineStyle(2)
slider_yvsx.Draw("L SAME")
c.Update()
# --- Start multi-page PDF ---
c.Print(pdf_name + "(")

# --- Close multi-page PDF ---
c.Print(pdf_name + ")")

print("Hit return to exit")
sys.stdin.readline()





