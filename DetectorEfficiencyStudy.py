import ROOT
from ROOT import TCanvas
from array import array
import numpy as np
import uproot

latex = ROOT . TLatex ()
latex . SetNDC ()
latex . SetTextSize (0.03)


LambdaFile = ROOT.TFile("Lambda_8Apr.root")
DstarFile = ROOT.TFile("Dstar_8Apr.root")
GammaFile = ROOT.TFile("Gamma_8Apr.root")

LambdaTree = LambdaFile.Get("Lambda")
DstarTree = DstarFile.Get("Dstar")
GammaTree = GammaFile.Get("Gamma")

Histogram2D_FirstElectron = ROOT.TH2D("Histogram2D_FirstElectron", "First Electron Energy vs Momentum;Energy (arb units);Momentum (GeV/c)", 100, 0,20,100,0, 5e6)
Histogram2D_SecondElectron = ROOT.TH2D("Histogram2D_SecondElectron", "Second Electron Energy vs Momentum;Energy (arb units);Momentum (GeV/c)", 100, 0,20,100,0, 5e6)
Histogram2D_SecondElectronNHitsUsed = ROOT.TH2D("Histogram2D_SecondElectronNHitsUsed", "Second Electron Energy vs NHitsUsed;Energy (arb units);NHitsUsed (arb units)", 100, 0,20,100,0, 5e6)
Histogram2D_FirstElectronNHitsUsed = ROOT.TH2D("Histogram2D_FirstElectronNHitsUsed", "First Electron Energy vs NHitsUsed;Energy (arb units);NHitsUsed (arb units)", 100, 0,20,100,0, 5e6)


InvM = ROOT.RooRealVar("InvM", "m(Lambda)", 1.11, 1.12, "GeV/c^{2}")
InvM_Dstar = ROOT.RooRealVar("deltaM", "m(D*)-m(D0)", 0.14,0.151, "GeV/c^{2}")


file2 = ROOT.TFile.Open("Gamma_8Apr.root")
preselTreeGamma = file2.Get("Gamma")
if not preselTreeGamma:
    raise RuntimeError("Failed to get 'Gamma' tree from Gamma_8Apr.root")
print('Files Loaded')

FirstElectronSVDdEdxList = np.empty(20, "float64")
FirstElectronSVDdEdx = 0

preselTreeGamma.Branch("FirstElectronSVDdEdxList",FirstElectronSVDdEdxList,"FirstElectronSVDdEdxList[20]/D")
preselTreeGamma.Branch("FirstElectronSVDdEdx",FirstElectronSVDdEdx,"FirstElectronSVDdEdx/D")

# for entry in preselTreeGamma:
#     print(entry.FirstElectronSVDdEdx)
#     print(entry.FirstElectronSVDdEdxList)
# #    for elem in entry.FirstElectronSVDdEdxList:
#  #       print(elem)

with uproot.open("Gamma_8Apr.root") as file:
    tree = file["Gamma"]
    FirstElectronSVDdEdxListArray = tree["FirstElectronSVDdEdxList"].array(library="np")

with uproot.open("Gamma_8Apr.root") as file:
    tree = file["Gamma"]
    FirstElectronSVDdEdxArray = tree["FirstElectronSVDdEdx"].array(library="np")

# event = 18
# bins_start = 2
# bins_fin = -3

# tester = FirstElectronSVDdEdxListArray[event]
# print(f'before sorting size: {tester.size}')
# tester = np.sort(tester)
# print(f'after sorting size: {tester.size}')

# if bins_fin == 100:
#     print(f'Sorted List Array: {tester[bins_start:]}')
#     print(f'Mean of Array: {np.mean(tester[bins_start:])}')
#     print(f'Expected Value: {FirstElectronSVDdEdxArray[event]}')
#     print(f'Difference of Values: {FirstElectronSVDdEdxArray[event] - np.mean(tester[bins_start:])}')    
# else:
#     print(f'Sorted List Array: {tester[bins_start:bins_fin]}')
#     print(f'Mean of Array: {np.mean(tester[bins_start:bins_fin])}')
#     print(f'Expected Value: {FirstElectronSVDdEdxArray[event]}')
#     print(f'Difference of Values: {FirstElectronSVDdEdxArray[event] - np.mean(tester[bins_start:bins_fin])}')
    

# print('Even ones')

event = 3826

tester = FirstElectronSVDdEdxListArray[event]

even_values = tester[::2]

even_sorted = np.sort(even_values)
even_without_max2 = even_sorted[:-2]

print(f'Even indices values: {even_values}')
print(f'Even indices sorted: {even_sorted}')
print(f'Even indices without max 2: {even_without_max2}')
print(f'Mean of even indices without max 2: {np.mean(even_without_max2)}')
print(f'Expected Value: {FirstElectronSVDdEdxArray[event]}')

# for i in range(2000,10000):
#     tester = FirstElectronSVDdEdxListArray[i]
#     even_values = tester[::2]
#     even_sorted = np.sort(even_values)
#     even_without_max2 = even_sorted[:-2]
#     Difference = np.mean(even_without_max2) - FirstElectronSVDdEdxArray[i]
#     # print(f'Difference between values: {Difference}')
#     if Difference != 0 and ~np.isnan(Difference):
#         print(f'Event {i}')


# Event 2271
# Event 2736
# Event 3443
# Event 3826
# Event 4556
# Event 4640
# Event 5203
# Event 6384
# Event 8813
# Event 9239
# Event 9556

# print(f'Size of event 11, {FirstElectronSVDdEdxListArray[11].size}')
# print(f'Size of event 12, {FirstElectronSVDdEdxListArray[12].size}')
# print(f'Size of event 13, {FirstElectronSVDdEdxListArray[13].size}')
# print(f'Size of event 14, {FirstElectronSVDdEdxListArray[14].size}')
# print(f'Size of event 15, {FirstElectronSVDdEdxListArray[15].size}')
# print(f'Size of event 16, {FirstElectronSVDdEdxListArray[16].size}')