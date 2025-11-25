# plotting script for baseball2.cpp, plots all 4 pitches to the same pdf
# Requires that baseball2 has been run first for all 4 pitches

import ROOT as r      # needed to load ROOT libraries
#import numpy as np    # only needed to use numpy arrays, see example below
import sys

# pdf name
pdf_name = "pitches.pdf"
pitches = ["slider","curveball","screwball","fastball"]

# Create a canvas
c = r.TCanvas("c", "c", 800, 600)


for i in range(len(pitches)):
    # slider pitch file
    tf=r.TFile(f"{pitches[i]}.root") # open file for read access
    if tf.IsZombie():
        print(f"{pitches[i]}.root not found in current directory, run baseball2 first")
        sys.exit(1)

    zvsx=tf.Get("g1")   # z vs x graph
    yvsx=tf.Get("g2")   # y vs x graph

    # plot on canvas
    c.Clear()
    h = zvsx.GetHistogram()
    h.GetXaxis().SetRangeUser(0, 60.)
    h.GetYaxis().SetRangeUser(-4, 2)
    zvsx.SetTitle(pitches[i] + ";x (ft);y (ft) / z (ft)")
    zvsx.Draw("AC")
    yvsx.SetLineStyle(2)
    yvsx.Draw("L SAME")
    c.Update()

    if i == 0:
        # --- Start multi-page PDF ---
        c.Print(pdf_name + "(")
    elif i == len(pitches)-1:
        # --- Close multi-page PDF ---
        c.Print(pdf_name + ")")
    else:
        # save to next page
        c.Print(pdf_name)

print("Hit return to exit")
sys.stdin.readline()





