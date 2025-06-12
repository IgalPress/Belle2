import ROOT
from ROOT import TCanvas
from array import array
import numpy as np

latex = ROOT . TLatex ()
latex . SetNDC ()
latex . SetTextSize (0.03)


LambdaFile = ROOT.TFile("Lambda_8Apr.root")
DstarFile = ROOT.TFile("Dstar_8Apr.root")
GammaFile = ROOT.TFile("Gamma_8Apr.root")

LambdaTree = LambdaFile.Get("Lambda")
DstarTree = DstarFile.Get("Dstar")
GammaTree = GammaFile.Get("Gamma")

Histogram2D_Kaon = ROOT.TH2D("Histogram2D_Kaon", "Kaon Energy vs Momentum;Energy (arb units);Momentum (GeV/c)",  100, 0,20,100,0, 5e6)
Histogram2D_Proton = ROOT.TH2D("Histogram2D_Proton", "Proton Energy vs Momentum;Energy (arb units);Momentum (GeV/c)",  100, 0,20,100,0, 5e6)
Histogram2D_Lambda = ROOT.TH2D("Histogram2D_Lambda", "Lambda Pion Energy vs Momentum;Energy (arb units);Momentum  (GeV/c)", 100, 0,20,100,0, 5e6)
Histogram2D_DstarPion = ROOT.TH2D("Histogram2D_DstarPion", "D* Pion Energy vs Momentum;Energy (arb units);Momentum (GeV/c)",  100, 0,20,100,0, 5e6)
Histogram2D_DstarSlowPion = ROOT.TH2D("Histogram2D_DstarSlowPion", "Slow Pion Energy vs Momentum;Energy (arb units);Momentum (GeV/c)", 100, 0,20,100,0, 5e6)
Histogram2D_FirstElectron = ROOT.TH2D("Histogram2D_FirstElectron", "First Electron Energy vs Momentum;Energy (arb units);Momentum (GeV/c)", 100, 0,20,100,0, 5e6)
Histogram2D_SecondElectron = ROOT.TH2D("Histogram2D_SecondElectron", "Second Electron Energy vs Momentum;Energy (arb units);Momentum (GeV/c)", 100, 0,20,100,0, 5e6)

InvM = ROOT.RooRealVar("InvM", "m(Lambda)", 1.11, 1.12, "GeV/c^{2}")
InvM_Dstar = ROOT.RooRealVar("deltaM", "m(D*)-m(D0)", 0.14,0.151, "GeV/c^{2}")

file0 = ROOT.TFile.Open("Lambda_8Apr.root")
preselTree = file0.Get("Lambda")
if not preselTree:
    raise RuntimeError("Failed to get 'Lambda' tree from Lambda_8Apr.root")

file1 = ROOT.TFile.Open("Dstar_8Apr.root")
preselTreeDstar = file1.Get("Dstar")
if not preselTreeDstar:
    raise RuntimeError("Failed to get 'Dstar' tree from Dstar_8Apr.root")
print('Files Loaded')

file2 = ROOT.TFile.Open("Gamma_8Apr.root")
preselTreeGamma = file2.Get("Gamma")
if not preselTreeGamma:
    raise RuntimeError("Failed to get 'Gamma' tree from Gamma_8Apr.root")
print('Files Loaded')

# ROOT.gStyle.SetStatW(0.35)
# ROOT.gStyle.SetStatH(0.35)
ROOT.gStyle.SetOptFit(1)

variables = ROOT.RooArgSet()
variables_Dstar = ROOT.RooArgSet()
variables_Gamma = ROOT.RooArgSet()
variables.add(InvM)
variables_Dstar.add(InvM_Dstar)

print('Data Loaded')

KaonMomentum = ROOT.RooRealVar("KaonMomentum", "", -1e8, 1e8)
KaonSVDdEdx = ROOT.RooRealVar("KaonSVDdEdx", "", -1e8, 1e8)
PionDMomentum = ROOT.RooRealVar("PionDMomentum", "", -1e8, 1e8)
PionDSVDdEdx = ROOT.RooRealVar("PionDSVDdEdx", "", -1e8, 1e8)
SlowPionMomentum = ROOT.RooRealVar("SlowPionMomentum", "", -1e8, 1e8)
SlowPionSVDdEdx = ROOT.RooRealVar("SlowPionSVDdEdx", "", -1e8, 1e8)
event_Dstar = ROOT.RooRealVar("event","",-1e8,1e8)

ProtonMomentum = ROOT.RooRealVar("ProtonMomentum", "momentum for Proton (GeV)", -1e8, 1e8)
ProtonSVDdEdxTrackMomentum = ROOT.RooRealVar("ProtonSVDdEdxTrackMomentum", "momentum for Proton (GeV), from the track", -1e8, 1e8)
ProtonSVDdEdx = ROOT.RooRealVar("ProtonSVDdEdx", "", -1e8, 1e8)
PionLambdaMomentum = ROOT.RooRealVar("PionLambdaMomentum", "momentum for pion (GeV)", -1e8, 1e8)
PionLambdaSVDdEdxTrackMomentum = ROOT.RooRealVar("PionLambdaSVDdEdxTrackMomentum", "momentum for pion (GeV), from the track", -1e8, 1e8)
PionLambdaSVDdEdx = ROOT.RooRealVar("PionLambdaSVDdEdx", "", -1e8, 1e8)
event_Lambda = ROOT.RooRealVar("event","",-1e8,1e8)

FirstElectronMomentum = ROOT.RooRealVar("FirstElectronMomentum", "", -1e8, 1e8)
FirstElectronSVDdEdx = ROOT.RooRealVar("FirstElectronSVDdEdx", "", -1e8, 1e8)
SecondElectronMomentum = ROOT.RooRealVar("SecondElectronMomentum", "", -1e8, 1e8)
SecondElectronSVDdEdx = ROOT.RooRealVar("SecondElectronSVDdEdx", "", -1e8, 1e8)


variables_Dstar.add(KaonMomentum)
variables_Dstar.add(PionDMomentum)
variables_Dstar.add(SlowPionMomentum)
variables_Dstar.add(PionDSVDdEdx)
variables_Dstar.add(SlowPionSVDdEdx)
variables_Dstar.add(KaonSVDdEdx)
variables_Dstar.add(event_Dstar)
variables.add(ProtonMomentum)
variables.add(ProtonSVDdEdxTrackMomentum)
variables.add(ProtonSVDdEdx)
variables.add(PionLambdaMomentum)
variables.add(PionLambdaSVDdEdxTrackMomentum)
variables.add(PionLambdaSVDdEdx)
variables.add(event_Lambda)
variables_Gamma.add(FirstElectronMomentum)
variables_Gamma.add(FirstElectronSVDdEdx)
variables_Gamma.add(SecondElectronMomentum)
variables_Gamma.add(SecondElectronSVDdEdx)


LambdaDataset = ROOT.RooDataSet("LambdaDataset", "LambdaDataset", preselTree, variables, "ProtonMomentum > 0.25")
DstarDataset = ROOT.RooDataSet("DstarDataset", "DstarDataset", preselTreeDstar, variables_Dstar)
GammaDataset = ROOT.RooDataSet("GammaDataset", "GammaDataset", preselTreeGamma, variables_Gamma)


DstarDataset.get().Print("v")



GaussMean = ROOT.RooRealVar("GaussMean", "GaussMean", 1.156, 1.114, 1.117)
GaussSigma1 = ROOT.RooRealVar("GaussSigma1", "GaussSigma1", 0.0008, 0.0003, 0.002)
GaussSigma2 = ROOT.RooRealVar("GaussSigma2", "GaussSigma2", 0.0012, 0.0003, 0.003)

LambdaGauss1 = ROOT.RooGaussian("LambdaGauss1", "LambdaGauss1", InvM, GaussMean, GaussSigma1)
LambdaGauss2 = ROOT.RooGaussian("LambdaGauss2", "LambdaGauss2", InvM, GaussMean, GaussSigma2)

fracGaussYield = ROOT.RooRealVar("fracGaussYield", "Fraction of two Gaussians", 0.5, 0, 1)
LambdaSignalPDF = ROOT.RooAddPdf("LambdaSignalPDF", "LambdaGauss1+LambdaGauss2",
                                 ROOT.RooArgList(LambdaGauss1, LambdaGauss2),
                                 fracGaussYield)

BkgPolyCoef0 = ROOT.RooRealVar("BkgPolyCoef0", "BkgPolyCoef0", 0.1, 0.0, 1.5)
BkgPolyCoef1 = ROOT.RooRealVar("BkgPolyCoef1", "BkgPolyCoef1", -0.5, -1.5, -1e-3)
LambdaBkgPDF = ROOT.RooChebychev("LambdaBkgPDF", "BkgPolyPDF", InvM,
                                 ROOT.RooArgList(BkgPolyCoef0, BkgPolyCoef1))
n_entries = preselTree.GetEntries()
nSignalLambda = ROOT.RooRealVar("nSignalLambda", "signal yield", 0.5 * n_entries, 0, n_entries)
nBkgLambda = ROOT.RooRealVar("nBkgLambda", "background yield", 0.5 * n_entries, 0, n_entries)


nSignalLambda = ROOT.RooRealVar("nSignalLambda","signal yield",0.5*preselTree.GetEntries(),0, preselTree.GetEntries())

totalPDFLambda = ROOT.RooAddPdf("totalPDFLambda","total PDF",ROOT.RooArgList(LambdaSignalPDF, LambdaBkgPDF),ROOT.RooArgList(nSignalLambda, nBkgLambda))

sPlotDatasetLambda = ROOT.RooStats.SPlot("sData", "An SPlot", LambdaDataset, totalPDFLambda, ROOT.RooArgList(nSignalLambda, nBkgLambda))

LambdaDatasetSWeighted = ROOT.RooDataSet(LambdaDataset.GetName(), LambdaDataset.GetTitle(), LambdaDataset, LambdaDataset.get()) 

treeLambda_sw = LambdaDatasetSWeighted.GetClonedTree()
treeLambda_sw.SetName("treeLambda_sw") 

################################

GaussMean_D = ROOT.RooRealVar("GaussMean_D", "GaussMean", 0.145, 0.14, 0.15)
GaussSigma1_D = ROOT.RooRealVar("GaussSigma1_D", "GaussSigma1", 0.01, 0.0001, 0.03)
GaussSigma2_D = ROOT.RooRealVar("GaussSigma2_D", "GaussSigma2", 0.015, 0.0001, 0.05)
fracGaussYield_D = ROOT.RooRealVar("fracGaussYield_D", "frac", 0.5, 0.0, 1.0)

Gauss1_D = ROOT.RooGaussian("Gauss1_D", "G1", InvM_Dstar, GaussMean_D, GaussSigma1_D)
Gauss2_D = ROOT.RooGaussian("Gauss2_D", "G2", InvM_Dstar, GaussMean_D, GaussSigma2_D)


fracGaussYield_D = ROOT.RooRealVar("fracGaussYield_D", "Fraction of two Gaussians", 0.5, 0, 1)


BkgPolyCoef0_D = ROOT.RooRealVar("BkgPolyCoef0_D", "BkgPolyCoef0", 0.1, 0.0, 1.5)
BkgPolyCoef1_D = ROOT.RooRealVar("BkgPolyCoef1_D", "BkgPolyCoef1", -0.5, -1.5, -1e-3)
DstarBkgPDF = ROOT.RooChebychev("DstarBkgPDF", "BkgPolyPDF_D", InvM_Dstar,
                                 ROOT.RooArgList(BkgPolyCoef0_D, BkgPolyCoef1_D))

n_entries_Dstar = preselTreeDstar.GetEntries()
nSignalDstar = ROOT.RooRealVar("nSignalDstar", "signal yield", 0.5 * n_entries_Dstar, 0, n_entries_Dstar)
nBkgDstar = ROOT.RooRealVar("nBkgDstar", "background yield", 0.5 * n_entries_Dstar, 0, n_entries_Dstar)

DstarSignalPDF = ROOT.RooAddPdf("DstarSignalPDF", "Signal", ROOT.RooArgList(Gauss1_D, Gauss2_D), ROOT.RooArgList(fracGaussYield_D))

nSignalDstar = ROOT.RooRealVar("nSignalDstar","signal yield",0.5*preselTreeDstar.GetEntries(),0, preselTreeDstar.GetEntries())

totalPDFDstar = ROOT.RooAddPdf("totalPDFDstar","total PDF",ROOT.RooArgList(DstarSignalPDF, DstarBkgPDF),ROOT.RooArgList(nSignalDstar, nBkgDstar))



###################################################################

DstarFitFrame = InvM_Dstar.frame(ROOT.RooFit.Title("D* Mass Fit"))

DstarFitResult = totalPDFDstar.fitTo(DstarDataset, ROOT.RooFit.Save(True)) 
DstarDataset.plotOn(DstarFitFrame, ROOT.RooFit.Name("DstarData"), ROOT.RooFit.MarkerColor(ROOT.kBlack), ROOT.RooFit.MarkerStyle(20), ROOT.RooFit.MarkerSize(0.5))
DstarFitFrame.SetTitle("D* Mass Fit")
totalPDFDstar.plotOn(DstarFitFrame, ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#4575b4")))
totalPDFDstar.paramOn(DstarFitFrame, ROOT.RooFit.Layout(0.63, 0.96, 0.93), ROOT.RooFit.Format("NEU", ROOT.RooFit.AutoPrecision(2)))
DstarFitFrame.getAttText().SetTextSize(0.03)

totalPDFDstar.plotOn(DstarFitFrame, ROOT.RooFit.Components("DstarSignalPDF"), ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#d73027")))
totalPDFDstar.plotOn(DstarFitFrame, ROOT.RooFit.Components("DstarBkgPDF"), ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#fc8d59")))
totalPDFDstar.plotOn(DstarFitFrame, ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#4575b4")))

DstarFitFrame.GetXaxis().SetTitle("#Deltam [GeV/c^{2}]")

canvDstar = ROOT.TCanvas("canvDstar", "canvDstar")
canvDstar.cd()
DstarFitFrame.Draw()

m_isMakePlots = True  
canvDstar.Print("SVDdEdxCalibrationFitDstar.pdf")

###################################################################





sPlotDatasetDstar = ROOT.RooStats.SPlot("sDataDstar", "An SPlot", DstarDataset, totalPDFDstar, ROOT.RooArgList(nSignalDstar, nBkgDstar))

DstarDatasetSWeighted = ROOT.RooDataSet(DstarDataset.GetName(), DstarDataset.GetTitle(), DstarDataset, DstarDataset.get()) 

treeDstar_sw = DstarDatasetSWeighted.GetClonedTree()
treeDstar_sw.SetName("treeDstar_sw") 

GammaDatasetSWeighted = ROOT.RooDataSet(GammaDataset.GetName(), GammaDataset.GetTitle(), GammaDataset, GammaDataset.get()) 

treeGamma_sw = GammaDatasetSWeighted.GetClonedTree()
treeGamma_sw.SetName("treeGamma_sw") 

#########################################################################################################################################



# Set estimate for the number of entries
treeGamma_sw.SetEstimate(treeGamma_sw.GetEntries() + 1)

# Draw the variable to fill the internal array (no plot, just data)
treeGamma_sw.Draw("FirstElectronMomentum/0.0005110", "", "goff")

# Retrieve the array of values
vXP = treeGamma_sw.GetV1()

# Number of bins
m_numPBins = 100

kdBinsP = ROOT.TKDTreeBinning(treeGamma_sw.GetEntries(), 1, vXP, m_numPBins)
binsMinEdgesP_orig = kdBinsP.SortOneDimBinEdges()

binsMinEdgesP = [binsMinEdgesP_orig[i] for i in range(m_numPBins + 2)]

# Fix the first and last bin edges
binsMinEdgesP[0] = 0.1
binsMinEdgesP[m_numPBins + 1] = 50.

# Print bin edges for checking
for i in range(m_numPBins + 2):
    print(binsMinEdgesP[i])

# Define the histogram with custom binning
m_numDEdxBins = 200  
m_dedxCutoff = 5e6   

# Convert to C array for ROOT
binsMinEdgesP_arr = array('d', binsMinEdgesP)

hGammaFirstElectron_bg = ROOT.TH2F(
    "hGammaFirstElectron",
    "hGammaFirstElectron",
    m_numPBins,
    binsMinEdgesP_arr,
    m_numDEdxBins,
    0,
    m_dedxCutoff
)

treeGamma_sw.Draw(
    "FirstElectronSVDdEdx:FirstElectronMomentum/0.0005110>>hGammaFirstElectron",
    "(FirstElectronSVDdEdx>0)",
    "goff"
)

cFirstElectronCustom = ROOT.TCanvas("cFirstElectronCustom", "Custom Binning FirstElectron", 800, 600)
hGammaFirstElectron_bg.Draw("COLZ")
cFirstElectronCustom.SaveAs("FirstElectronCustomBinning.png")



cFirstElectronCustomProjectionY = ROOT.TCanvas("cFirstElectronCustomProjectionY", "Custom Binning FirstElectron 2D ProjectionY", 800, 600)

cFirstElectronCustomProjectionY.SetRightMargin(0.15)
cFirstElectronCustomProjectionY.SetLeftMargin(0.12)
# cFirstElectronCustomProjectionY.SetLogy(True)
# cFirstElectronCustomProjectionY.SetLogx(True)


CustomProjectionY_FirstElectron = hGammaFirstElectron_bg.ProjectionY("CustomProjectionY_FirstElectron",80,150)
CustomProjectionY_FirstElectron.SetLineColor(ROOT.kRed)
CustomProjectionY_FirstElectron.SetLineWidth(2)
CustomProjectionY_FirstElectron.SetTitle("First Electron ProjectionY")
CustomProjectionY_FirstElectron.GetYaxis().SetTitle("Mean Energy Loss")
CustomProjectionY_FirstElectron.GetXaxis().SetTitle("P/M (GeV/C)")
CustomProjectionY_FirstElectron.Draw("E")



gaus2Landau = ROOT.TF1(
    "gaus2Landau",
    "[0]*TMath::Gaus(x, [1], [1]*[2]*[5]) + [3]*TMath::Gaus(x, [1], [1]*[4]*[5])*TMath::Landau(x, [1], [1]*[5])",
    1.e5, 1.5e6
)

gaus2Landau.SetLineColor(ROOT.kBlack)
# Parameters: gaus1_amp, gaus1_mean, gaus1_sigma, gaus2_amp, gaus2_mean, gaus2_sigma, landau_amp, landau_mp, landau_width
gaus2Landau.SetParameters(9000, 650e3, 1, 4000, 1, 0.1)
gaus2Landau.SetNpx(1000)


gaus2Landau.SetParLimits(5,0,1)
gaus2Landau.SetParLimits(2,1,5)
gaus2Landau.SetParLimits(4,1,5)

CustomProjectionY_FirstElectron.Fit("gaus2Landau", "R")


cFirstElectronCustomProjectionY.SaveAs("CustomBinningFirstElectron2DHistogram_ProjectionY.png")


##########################################################################################################################################


# Set estimate for the number of entries
treeGamma_sw.SetEstimate(treeGamma_sw.GetEntries() + 1)

# Draw the variable to fill the internal array (no plot, just data)
treeGamma_sw.Draw("SecondElectronMomentum/0.0005110", "", "goff")

# Retrieve the array of values
vXP = treeGamma_sw.GetV1()

# Number of bins
m_numPBins = 100

kdBinsP = ROOT.TKDTreeBinning(treeGamma_sw.GetEntries(), 1, vXP, m_numPBins)
binsMinEdgesP_orig = kdBinsP.SortOneDimBinEdges()

binsMinEdgesP = [binsMinEdgesP_orig[i] for i in range(m_numPBins + 2)]

# Fix the first and last bin edges
binsMinEdgesP[0] = 0.1
binsMinEdgesP[m_numPBins + 1] = 50.

# Print bin edges for checking
for i in range(m_numPBins + 2):
    print(binsMinEdgesP[i])


# Convert to C array for ROOT
binsMinEdgesP_arr = array('d', binsMinEdgesP)

hGammaSecondElectron_bg = ROOT.TH2F(
    "hGammaSecondElectron",
    "hGammaSecondElectron",
    m_numPBins,
    binsMinEdgesP_arr,
    m_numDEdxBins,
    0,
    m_dedxCutoff
)

treeGamma_sw.Draw(
    "SecondElectronSVDdEdx:SecondElectronMomentum/0.0005110>>hGammaSecondElectron",
    "(SecondElectronSVDdEdx>0)",
    "goff"
)

cSecondElectronCustom = ROOT.TCanvas("cSecondElectronCustom", "Custom Binning SecondElectron", 800, 600)
hGammaSecondElectron_bg.Draw("COLZ")
cSecondElectronCustom.SaveAs("SecondElectronCustomBinning.png")



cSecondElectronCustomProjectionY = ROOT.TCanvas("cSecondElectronCustomProjectionY", "Custom Binning SecondElectron 2D ProjectionY", 800, 600)

cSecondElectronCustomProjectionY.SetRightMargin(0.15)
cSecondElectronCustomProjectionY.SetLeftMargin(0.12)
# cSecondElectronCustomProjectionY.SetLogy(True)
# cSecondElectronCustomProjectionY.SetLogx(True)


CustomProjectionY_SecondElectron = hGammaSecondElectron_bg.ProjectionY("CustomProjectionY_SecondElectron",80,150)
CustomProjectionY_SecondElectron.SetLineColor(ROOT.kRed)
CustomProjectionY_SecondElectron.SetLineWidth(2)
CustomProjectionY_SecondElectron.SetTitle("Second Electron ProjectionY")
CustomProjectionY_SecondElectron.GetYaxis().SetTitle("Mean Energy Loss")
CustomProjectionY_SecondElectron.GetXaxis().SetTitle("P/M (GeV/C)")
CustomProjectionY_SecondElectron.Draw("E")



gaus2Landau = ROOT.TF1(
    "gaus2Landau",
    "[0]*TMath::Gaus(x, [1], [1]*[2]*[5]) + [3]*TMath::Gaus(x, [1], [1]*[4]*[5])*TMath::Landau(x, [1], [1]*[5])",
    1.e5, 1.5e6
)

gaus2Landau.SetLineColor(ROOT.kBlack)
# Parameters: gaus1_amp, gaus1_mean, gaus1_sigma, gaus2_amp, gaus2_sigma, landau_width
gaus2Landau.SetParameters(20000, 650e3, 1, 4000, 1, 0.1)
gaus2Landau.SetNpx(1000)


gaus2Landau.SetParLimits(5,0,1)
gaus2Landau.SetParLimits(2,1,5)
gaus2Landau.SetParLimits(4,1,5)

CustomProjectionY_SecondElectron.Fit("gaus2Landau", "R")


cSecondElectronCustomProjectionY.SaveAs("CustomBinningSecondElectron2DHistogram_ProjectionY.png")



######################################################################################################################################


# Set estimate for the number of entries
treeDstar_sw.SetEstimate(treeDstar_sw.GetEntries() + 1)

# Draw the variable to fill the internal array (no plot, just data)
treeDstar_sw.Draw("KaonMomentum/0.4937", "", "goff")

# Retrieve the array of values
vXP = treeDstar_sw.GetV1()

# Number of bins
m_numPBins = 100

kdBinsP = ROOT.TKDTreeBinning(treeDstar_sw.GetEntries(), 1, vXP, m_numPBins)
binsMinEdgesP_orig = kdBinsP.SortOneDimBinEdges()

binsMinEdgesP = [binsMinEdgesP_orig[i] for i in range(m_numPBins + 2)]

# Fix the first and last bin edges
binsMinEdgesP[0] = 0.1
binsMinEdgesP[m_numPBins + 1] = 50.

# Print bin edges for checking
for i in range(m_numPBins + 2):
    print(binsMinEdgesP[i])

# Define the histogram with custom binning
m_numDEdxBins = 200  # or your preferred value
m_dedxCutoff = 5e6   # or your preferred value

# Convert to C array for ROOT
binsMinEdgesP_arr = array('d', binsMinEdgesP)

hDstarK_bg = ROOT.TH2F(
    "hDstarK",
    "hDstarK",
    m_numPBins,
    binsMinEdgesP_arr,
    m_numDEdxBins,
    0,
    m_dedxCutoff
)

treeDstar_sw.Draw(
    "KaonSVDdEdx:KaonMomentum/0.4937>>hDstarK",
    "nSignalDstar_sw * (KaonSVDdEdx>0)",
    "goff"
)

cKaonCustom = ROOT.TCanvas("cKaonCustom", "Custom Binning Kaon", 800, 600)
hDstarK_bg.Draw("COLZ")
cKaonCustom.SaveAs("KaonCustomBinning.png")



cKaonCustomProjectionY = ROOT.TCanvas("cKaonCustomProjectionY", "Custom Binning Kaon 2D ProjectionY", 800, 600)

cKaonCustomProjectionY.SetRightMargin(0.15)
cKaonCustomProjectionY.SetLeftMargin(0.12)
# cKaonCustomProjectionY.SetLogy(True)
# cKaonCustomProjectionY.SetLogx(True)


CustomProjectionY_Kaon = hDstarK_bg.ProjectionY("CustomProjectionY_Kaon",80,150)
CustomProjectionY_Kaon.SetLineColor(ROOT.kRed)
CustomProjectionY_Kaon.SetLineWidth(2)
CustomProjectionY_Kaon.SetTitle("Kaon Pion ProjectionY")
CustomProjectionY_Kaon.GetYaxis().SetTitle("Mean Energy Loss")
CustomProjectionY_Kaon.GetXaxis().SetTitle("P/M (GeV/C)")
ROOT.gStyle.SetStatW(0.35)
ROOT.gStyle.SetStatH(0.35)
CustomProjectionY_Kaon.Draw("E")


ROOT.gStyle.SetOptFit(1)

gaus2Landau = ROOT.TF1(
    "gaus2Landau",
    "[0]*TMath::Gaus(x, [1], [1]*[2]*[5]) + [3]*TMath::Gaus(x, [1], [1]*[4]*[5])*TMath::Landau(x, [1], [1]*[5])",
    2.e5, 1.5e6
)

gaus2Landau.SetLineColor(ROOT.kBlack)
# Parameters: gaus1_amp, gaus1_mean, gaus1_sigma, gaus2_amp, gaus2_mean, gaus2_sigma, landau_amp, landau_mp, landau_width
gaus2Landau.SetParameters(8000, 650e3, 1, 4000, 1, 0.1)
gaus2Landau.SetNpx(1000)


gaus2Landau.SetParLimits(5,0,1)
gaus2Landau.SetParLimits(2,1,5)
gaus2Landau.SetParLimits(4,1,5)

CustomProjectionY_Kaon.Fit("gaus2Landau", "R")

mean1  = gaus2Landau.GetParameter(1)
sigma1 = gaus2Landau.GetParameter(2)
sigma2 = gaus2Landau.GetParameter(4)


cKaonCustomProjectionY.SaveAs("CustomBinningKaon2DHistogram_ProjectionY.png")



#############################################################################################################################
# Set estimate for the number of entries
treeLambda_sw.SetEstimate(treeLambda_sw.GetEntries() + 1)

# Draw the variable to fill the internal array (no plot, just data)
treeLambda_sw.Draw("PionLambdaMomentum/0.1396", "", "goff")

# Retrieve the array of values
vXP = treeLambda_sw.GetV1()

# Number of bins
m_numPBins = 100

kdBinsP = ROOT.TKDTreeBinning(treeLambda_sw.GetEntries(), 1, vXP, m_numPBins)
binsMinEdgesP_orig = kdBinsP.SortOneDimBinEdges()

binsMinEdgesP = [binsMinEdgesP_orig[i] for i in range(m_numPBins + 2)]

# Fix the first and last bin edges
binsMinEdgesP[0] = 0.1
binsMinEdgesP[m_numPBins + 1] = 50.

# Print bin edges for checking
for i in range(m_numPBins + 2):
    print(binsMinEdgesP[i])

# Convert to C array for ROOT
binsMinEdgesP_arr = array('d', binsMinEdgesP)

hPionLambda_bg = ROOT.TH2F(
    "hPionLambda",
    "hPionLambda",
    m_numPBins,
    binsMinEdgesP_arr,
    m_numDEdxBins,
    0,
    m_dedxCutoff
)

treeLambda_sw.Draw(
    "PionLambdaSVDdEdx:PionLambdaMomentum/0.1396>>hPionLambda",
    "nSignalLambda_sw * (PionLambdaSVDdEdx>0)",
    "goff"
)

cPionLambdaCustom = ROOT.TCanvas("cPionLambdaCustom", "Custom Binning PionLambda 2D", 800, 600)
hPionLambda_bg.Draw("COLZ")
cPionLambdaCustom.SaveAs("PionLambdaCustomBinning.png")



########################################################################################################
# Set estimate for the number of entries
treeDstar_sw.SetEstimate(treeDstar_sw.GetEntries() + 1)

# Draw the variable to fill the internal array (no plot, just data)
treeDstar_sw.Draw("PionDMomentum/0.1396", "", "goff")

# Retrieve the array of values
vXP = treeDstar_sw.GetV1()

# Number of bins
m_numPBins = 100

kdBinsP = ROOT.TKDTreeBinning(treeDstar_sw.GetEntries(), 1, vXP, m_numPBins)
binsMinEdgesP_orig = kdBinsP.SortOneDimBinEdges()

binsMinEdgesP = [binsMinEdgesP_orig[i] for i in range(m_numPBins + 2)]

# Fix the first and last bin edges
binsMinEdgesP[0] = 0.1
binsMinEdgesP[m_numPBins + 1] = 50.

# Print bin edges for checking
for i in range(m_numPBins + 2):
    print(binsMinEdgesP[i])

# Convert to C array for ROOT
binsMinEdgesP_arr = array('d', binsMinEdgesP)

hPionD_bg = ROOT.TH2F(
    "hPionD",
    "hPionD",
    m_numPBins,
    binsMinEdgesP_arr,
    m_numDEdxBins,
    0,
    m_dedxCutoff
)

treeDstar_sw.Draw(
    "PionDSVDdEdx:PionDMomentum/0.1396>>hPionD",
    "nSignalDstar_sw * (PionDSVDdEdx>0)",
    "goff"
)

cPionDCustom = ROOT.TCanvas("cPionDstar", "Custom Binning PionD* 2D", 800, 600)
hPionD_bg.Draw("COLZ")
cPionDCustom.SaveAs("PionDCustomBinning.png")

######################################################################################################################################



# Set estimate for the number of entries
treeDstar_sw.SetEstimate(treeDstar_sw.GetEntries() + 1)

# Draw the variable to fill the internal array (no plot, just data)
treeDstar_sw.Draw("SlowPionMomentum/0.1396", "", "goff")

# Retrieve the array of values
vXP = treeDstar_sw.GetV1()

# Number of bins
m_numPBins = 100

kdBinsP = ROOT.TKDTreeBinning(treeDstar_sw.GetEntries(), 1, vXP, m_numPBins)
binsMinEdgesP_orig = kdBinsP.SortOneDimBinEdges()

binsMinEdgesP = [binsMinEdgesP_orig[i] for i in range(m_numPBins + 2)]

# Fix the first and last bin edges
binsMinEdgesP[0] = 0.1
binsMinEdgesP[m_numPBins + 1] = 50.

# Print bin edges for checking
for i in range(m_numPBins + 2):
    print(binsMinEdgesP[i])

binsMinEdgesP_arr = array('d', binsMinEdgesP)

hSlowPion_bg = ROOT.TH2F(
    "hSlowPion",
    "hSlowPion",
    m_numPBins,
    binsMinEdgesP_arr,
    m_numDEdxBins,
    0,
    m_dedxCutoff
)

treeDstar_sw.Draw(
    "SlowPionSVDdEdx:SlowPionMomentum/0.1396>>hSlowPion",
    "nSignalDstar_sw * (SlowPionSVDdEdx>0)",
    "goff"
)

cSlowPionCustom = ROOT.TCanvas("cSlowPionstar", "Custom Binning Slow Pion 2D", 800, 600)
hSlowPion_bg.Draw("COLZ")
cSlowPionCustom.SaveAs("SlowPionCustomBinning.png")

#############################################################################

legend = ROOT.TLegend(0.7 ,0.6 ,0.85 ,0.75)

treeDstar_sw.SetEstimate(treeDstar_sw.GetEntries() + 1)

treeDstar_sw.Draw("SlowPionMomentum/0.140*(event%2==0) + PionDMomentum/0.140*(event%2==1)", "", "goff")
n_entries = treeDstar_sw.GetSelectedRows()
vXP = treeDstar_sw.GetV1()

xp_values = [vXP[i] for i in range(n_entries)]

bins_input = array('d', xp_values)


m_numPBins = 100
kdBinsP = ROOT.TKDTreeBinning(n_entries, 1, bins_input, m_numPBins)
binsMinEdgesP_orig = kdBinsP.SortOneDimBinEdges()

binsMinEdgesP = [binsMinEdgesP_orig[i] for i in range(m_numPBins + 2)]
binsMinEdgesP[0] = 0.1
binsMinEdgesP[-1] = 25.0
binsMinEdgesP_arr = array('d', binsMinEdgesP)

hPionLambda_bg = ROOT.TH2F("hPionLambda_bg", "PionLambda", m_numPBins, binsMinEdgesP_arr, m_numDEdxBins, 0, m_dedxCutoff)
hPionD_bg = ROOT.TH2F("hPionD_bg", "PionD", m_numPBins, binsMinEdgesP_arr, m_numDEdxBins, 0, m_dedxCutoff)
hSlowPion_bg = ROOT.TH2F("hSlowPion_bg", "SlowPion", m_numPBins, binsMinEdgesP_arr, m_numDEdxBins, 0, m_dedxCutoff)

hSlowPion_bg.Sumw2()
hPionD_bg.Sumw2()
hPionLambda_bg.Sumw2()

treeLambda_sw.Draw("PionLambdaSVDdEdx:PionLambdaMomentum/0.1396>>hPionLambda_bg", "nSignalLambda_sw*(PionLambdaSVDdEdx>0)", "goff")
treeDstar_sw.Draw("PionDSVDdEdx:PionDMomentum/0.1396>>hPionD_bg", "nSignalDstar_sw*(PionDSVDdEdx>0)", "goff")
treeDstar_sw.Draw("SlowPionSVDdEdx:SlowPionMomentum/0.1396>>hSlowPion_bg", "nSignalDstar_sw*(SlowPionSVDdEdx>0)", "goff")

hTotal = hPionLambda_bg.Clone("hTotal")
hTotal.Add(hPionD_bg)
hTotal.Add(hSlowPion_bg)


cCombined = ROOT.TCanvas("cCombined", "Combined 2D Histogram", 800, 600)
cCombined.SetRightMargin(0.15)
hTotal.SetMinimum(0.)
cCombined.SetLeftMargin(0.12) 
hTotal.SetTitle("Combined Energy vs Momentum Histogram")
hTotal.Draw("COLZ")
cCombined.SaveAs("Combined2DHistogram.png")



cCombinedCustomProjectionY = ROOT.TCanvas("cCombinedCustomProjectionY", "Custom Binning Combined 2D ProjectionY", 800, 600)

cCombinedCustomProjectionY.SetRightMargin(0.15)
cCombinedCustomProjectionY.SetLeftMargin(0.12)
# cCombinedCustomProjectionY.SetLogy(True)
# cCombinedCustomProjectionY.SetLogx(True)


CustomProjectionY_Combined = hTotal.ProjectionY("CustomProjectionY_Combined",80,150)
CustomProjectionY_Combined.SetLineColor(ROOT.kRed)
CustomProjectionY_Combined.SetLineWidth(2)
CustomProjectionY_Combined.SetTitle("Combined Pion ProjectionY")
CustomProjectionY_Combined.GetYaxis().SetTitle("Mean Energy Loss")
CustomProjectionY_Combined.GetXaxis().SetTitle("P/M (GeV/C)")
ROOT.gStyle.SetStatW(0.35)
ROOT.gStyle.SetStatH(0.35)
CustomProjectionY_Combined.Draw("E")


ROOT.gStyle.SetOptFit(1)

gaus2Landau = ROOT.TF1(
    "gaus2Landau",
    "[0]*TMath::Gaus(x, [1], [1]*[2]*[5]) + [3]*TMath::Gaus(x, [1], [1]*[4]*[5])*TMath::Landau(x, [1], [1]*[5])",
    2.e5, 1.5e6
)

gaus2Landau.SetLineColor(ROOT.kBlack)
# Parameters: gaus1_amp, gaus1_mean, gaus1_sigma, gaus2_amp, gaus2_mean, gaus2_sigma, landau_amp, landau_mp, landau_width
gaus2Landau.SetParameters(8000, 650e3, 1, 4000, 1, 0.1)
gaus2Landau.SetNpx(1000)


gaus2Landau.SetParLimits(5,0,1)
gaus2Landau.SetParLimits(2,1,5)
gaus2Landau.SetParLimits(4,1,5)

CustomProjectionY_Combined.Fit("gaus2Landau", "R")

mean1  = gaus2Landau.GetParameter(1)
sigma1 = gaus2Landau.GetParameter(2)
sigma2 = gaus2Landau.GetParameter(4)


cCombinedCustomProjectionY.SaveAs("CustomBinningCombined2DHistogram_ProjectionY.png")
