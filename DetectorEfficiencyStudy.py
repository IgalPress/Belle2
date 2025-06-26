import ROOT
from ROOT import TCanvas
from array import array
import numpy as np
import uproot
import pandas as pd
# from scipy import stats

latex = ROOT . TLatex ()
latex . SetNDC ()
latex . SetTextSize (0.03)


LambdaFile = ROOT.TFile("Lambda_8Apr.root")
DstarFile = ROOT.TFile("Dstar_8Apr.root")
GammaFile = ROOT.TFile("Gamma_8Apr.root")

LambdaTree = LambdaFile.Get("Lambda")
DstarTree = DstarFile.Get("Dstar")
GammaTree = GammaFile.Get("Gamma")


file2 = ROOT.TFile.Open("Gamma_8Apr.root")
preselTreeGamma = file2.Get("Gamma")
if not preselTreeGamma:
    raise RuntimeError("Failed to get 'Gamma' tree from Gamma_8Apr.root")
print('Files Loaded')

FirstElectronSVDdEdxList = np.empty(20, "float64")
FirstElectronSVDdEdx = 0
FirstElectronMomentum = array('d', [0.])

preselTreeGamma.Branch("FirstElectronSVDdEdxList",FirstElectronSVDdEdxList,"FirstElectronSVDdEdxList[20]/D")
preselTreeGamma.Branch("FirstElectronSVDdEdx",FirstElectronSVDdEdx,"FirstElectronSVDdEdx/D")
preselTreeGamma.Branch("FirstElectronMomentum",FirstElectronMomentum,"FirstElectronMomentum/D")



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

with uproot.open("Gamma_8Apr.root") as file:
    tree = file["Gamma"]
    FirstElectronMomentumArray = tree["FirstElectronMomentum"].array(library="np")

def CMS_mean(arr):
    arr = np.array(arr)
    if len(arr) < 2:
        return np.nan
    else:
        return ((1/len(arr)) * np.sum((1.0/arr)**2))**(-1/2)

def ATLAS_mean(arr):
    arr = np.array(arr)
    if len(arr) <= 2:
        return (np.sum(arr)/len(arr))
    elif len(arr) == 3 or len(arr == 4):
        arr = arr[:-1]
        return (np.sum(arr)/len(arr))
    else:
        arr = arr[:-2]
        return (np.sum(arr)/len(arr))


def make_ALICE_weights(arr):
    n = len(arr)
    weights = np.ones(n)
    if n >= 2:
        weights[-2:] = 0.5
    return weights

def ALICE_mean(arr):
    arr = np.array(arr)
    if len(arr) in (3,4):
        return (arr[0] + arr[1])/2
    elif len(arr) < 2:
        return np.nan
    else:
        return sum(arr * make_ALICE_weights(arr)) / sum(make_ALICE_weights(arr))
    


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



event = 26955

tester = FirstElectronSVDdEdxListArray[event]


inf_indices = [
    i for i, row in enumerate(FirstElectronSVDdEdxListArray)
    if np.any(np.isinf(np.asarray(row, dtype=np.float64)))
]

mask = np.ones(len(FirstElectronSVDdEdxArray), dtype=bool)
mask[inf_indices] = False

filtered_FirstElectronSVDdEdxArray = FirstElectronSVDdEdxArray[mask]
filtered_FirstElectronSVDdEdxListArray = FirstElectronSVDdEdxListArray[mask]
filtered_FirstElectronMomentumArray = FirstElectronMomentumArray[mask]



cms_means = np.array([CMS_mean(row) for row in filtered_FirstElectronSVDdEdxListArray], dtype='float64')
alice_means = np.array([ALICE_mean(row) for row in filtered_FirstElectronSVDdEdxListArray], dtype='float64')
atlas_means = np.array([ATLAS_mean(row) for row in filtered_FirstElectronSVDdEdxListArray], dtype='float64')
harmonic_means = np.array([
    np.nan if len(row) < 2 else len(row) / np.sum(1.0 / np.array(row))
    for row in filtered_FirstElectronSVDdEdxListArray
], dtype='float64')

output_file = ROOT.TFile("Gamma_8Apr_withCMSmean.root", "UPDATE")
output_file.Delete("Gamma_Means;*")

new_tree = ROOT.TTree("Gamma_Means", "Gamma tree with multiple mean branches")

cms_mean_value = array('d', [0.])
alice_mean_value = array('d', [0.])
atlas_mean_value = array('d', [0.])
harmonic_mean_value = array('d', [0.])
filtered_FirstElectronSVDdEdxArray_value = array('d', [0.])
filtered_FirstElectronMomentumArray_value = array('d', [0.])

new_tree.Branch("CMS_mean", cms_mean_value, "CMS_mean/D")
new_tree.Branch("ALICE_mean", alice_mean_value, "ALICE_mean/D")
new_tree.Branch("ATLAS_mean", atlas_mean_value, "ATLAS_mean/D")
new_tree.Branch("Harmonic_mean", harmonic_mean_value, "Harmonic_mean/D")
new_tree.Branch("filtered_FirstElectronSVDdEdxArray", filtered_FirstElectronSVDdEdxArray, "filtered_FirstElectronSVDdEdxArray/D")
new_tree.Branch("filtered_FirstElectronMomentumArray", filtered_FirstElectronMomentumArray, "filtered_FirstElectronMomentumArray/D")

for cms, alice, atlas, harm, filtered_energy, filtered_momentum in zip(cms_means, alice_means, atlas_means, harmonic_means, filtered_FirstElectronSVDdEdxArray,filtered_FirstElectronMomentumArray):
    cms_mean_value[0] = cms
    alice_mean_value[0] = alice
    atlas_mean_value[0] = atlas
    harmonic_mean_value[0] = harm
    filtered_FirstElectronSVDdEdxArray_value[0] = filtered_energy
    filtered_FirstElectronMomentumArray_value[0] = filtered_momentum

    new_tree.Fill()

new_tree.Write()
output_file.Close()

# odd_values = tester[1::2]
# even_values = tester[::2]

# odd_sorted = np.sort(odd_values)
# odd_without_max2 = odd_sorted[:-2]
# even_sorted = np.sort(even_values)
# even_without_max2 = even_sorted[:-2]

# Harmonic_mean_even = len(even_without_max2)/(np.sum(1.0/even_without_max2))
# print(f'Harmonic Mean: {Harmonic_mean_even}')
# # print(f'Harm mean scipy: {stats.hmean(even_values)}')

# print(f'Full array: {FirstElectronSVDdEdxListArray[event]}')
# print(f'Even indices values: {even_values}')
# # print(f'Odd indices values: {odd_values}')
# print(f'Even indices sorted: {even_sorted}')
# # print(f'Even indices without max 2: {even_without_max2}')
# print(f'Mean of even indices without max 2: {np.mean(even_without_max2)}')
# # print(f'Mean of odd indices without max 2: {np.mean(odd_without_max2)}')
# print(f'Expected Value: {FirstElectronSVDdEdxArray[event]}')

# n = 0
# n_neg = 0
# for i in range(len(filtered_FirstElectronSVDdEdxListArray)):
#     tester = filtered_FirstElectronSVDdEdxListArray[i]
#     even_values = tester[::2]
#     even_sorted = np.sort(even_values)
#     even_without_max2 = even_sorted[:-2]
#     # Harmonic_mean_even = len(even_values)/(np.sum(1.0/even_values))
#     # Difference = Harmonic_mean_even - np.mean(even_without_max2)
#     Difference = CMS_mean(even_without_max2) - np.mean(even_without_max2)
#     # Difference = np.mean(even_without_max2) - FirstElectronSVDdEdxArray[i]
#     # print(f'Difference between values: {Difference}')
#     if Difference > 100000 and ~np.isnan(Difference):
#         # print(f'Event {i}')
#         # print(f'Even indices without max 2: {even_without_max2}')
#         # print(f'Harmonic mean: {Harmonic_mean_even}')
#         # print(f'Mean of even indices without max 2: {np.mean(even_without_max2)}')
#         # print(f'Difference {np.abs(Difference)}')
#         # print(f'Expected Value: {FirstElectronSVDdEdxArray[i]}')
#         n = n+1
#     elif Difference < -100000 and ~np.isnan(Difference):
#         n_neg = n_neg+1
#     if (i%50000 == 0):
#         print(f'Event number {i}')

# print(f'Number of times that ATLAS_mean_even - np.mean > 100 000: {n}')
# print(f'Number of times that ATLAS_mean_even - np.mean < -100 000: {n_neg}')
# print(f'Percent chance for this to happen: {((n+n_neg)/len(filtered_FirstElectronSVDdEdxListArray))*100}%')