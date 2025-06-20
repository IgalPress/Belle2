import ROOT
from ROOT import TCanvas
from array import array
import numpy as np

latex = ROOT . TLatex ()
latex . SetNDC ()
latex . SetTextSize (0.03)


LambdaFile = ROOT.TFile("Lambda_8Apr.root")
DstarFile = ROOT.TFile("Dstar_8Apr.root")

LambdaTree = LambdaFile.Get("Lambda")
DstarTree = DstarFile.Get("Dstar")

# HistogramLambdaPionEnergy = ROOT.TH1D('','',100,0,5e6)
# HistogramDstarPionEnergy = ROOT.TH1D('','',100,0,5e6)
# HistogramDstarSlowPionEnergy = ROOT.TH1D('','',100,0,5e6)

# HistogramLambdaPionMomentum = ROOT.TH1D('','',100,0,4)
# HistogramDstarPionMomentum = ROOT.TH1D('','',100,0,4)
# HistogramDstarSlowPionMomentum = ROOT.TH1D('','',100,0,4)

Histogram2D_Kaon = ROOT.TH2D("Histogram2D_Kaon", "Kaon Energy vs Momentum;Energy (arb units);Momentum (GeV/c)",  100, 0,20,100,0, 5e6)
Histogram2D_Proton = ROOT.TH2D("Histogram2D_Proton", "Proton Energy vs Momentum;Energy (arb units);Momentum (GeV/c)",  100, 0,20,100,0, 5e6)
Histogram2D_Lambda = ROOT.TH2D("Histogram2D_Lambda", "Lambda Pion Energy vs Momentum;Energy (arb units);Momentum  (GeV/c)", 100, 0,20,100,0, 5e6)
Histogram2D_DstarPion = ROOT.TH2D("Histogram2D_DstarPion", "D* Pion Energy vs Momentum;Energy (arb units);Momentum (GeV/c)",  100, 0,20,100,0, 5e6)
Histogram2D_DstarSlowPion = ROOT.TH2D("Histogram2D_DstarSlowPion", "Slow Pion Energy vs Momentum;Energy (arb units);Momentum (GeV/c)", 100, 0,20,100,0, 5e6)
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

variables = ROOT.RooArgSet()
variables_Dstar = ROOT.RooArgSet()
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

LambdaDataset = ROOT.RooDataSet("LambdaDataset", "LambdaDataset", preselTree, variables, "ProtonMomentum > 0.25")
DstarDataset = ROOT.RooDataSet("DstarDataset", "DstarDataset", preselTreeDstar, variables_Dstar)

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

cCustomKaonProfileX = ROOT.TCanvas("cCustomKaonProfileX", "CustomBinning ProfileX of Kaons", 800, 600)
cCustomKaonProfileX.SetRightMargin(0.15)
cCustomKaonProfileX.SetLeftMargin(0.12)
cCustomKaonProfileX.SetLogy(True)
cCustomKaonProfileX.SetLogx(True)

CustomprofileX_Kaon = hDstarK_bg.ProfileX("CustomprofileX_Kaon")
CustomprofileX_Kaon.SetLineColor(ROOT.kRed)
CustomprofileX_Kaon.SetLineWidth(2)
CustomprofileX_Kaon.SetTitle("Custom Binning ProfileX")
CustomprofileX_Kaon.GetYaxis().SetTitle("Mean Energy Loss")
CustomprofileX_Kaon.GetXaxis().SetTitle("P/M (GeV/C)")

# max_error = 1e4 

# for bin in range(1, profileX_Kaon.GetNbinsX() + 1):
#     if profileX_Kaon.GetBinError(bin) > max_error:
#         profileX_Kaon.SetBinContent(bin, 0)
#         profileX_Kaon.SetBinError(bin, 0)

CustomprofileX_Kaon.Draw("E")
CustomprofileX_Kaon.SetStats(False)

FitKaon = ROOT.TF1(
    "FitKaon",
    "[0] + [1]/pow((x + [2]), [3]) + [4]*x",
    0.4, 24
)
FitKaon.SetParameters(100, 5000, 0.01, 12, 10)  # [6] is the slope of the rising tail

CustomprofileX_Kaon.Fit(FitKaon, "", "", 0.4, 8)

FitKaon.SetLineColor(ROOT.kBlue)
FitKaon.SetLineWidth(2)
FitKaon.Draw("SAME")

# Get parameter values
p0 = FitKaon.GetParameter(0)
p1 = FitKaon.GetParameter(1)
p2 = FitKaon.GetParameter(2)
p3 = FitKaon.GetParameter(3)
p4 = FitKaon.GetParameter(4)


# # Format the formula with parameter values
formula_str = (
    f"Fit: {p0:.2e} + {p1:.2e}/(x + {p2:.2e})^{p3:.2e} + {p4:.2e}*x"
)

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextSize(0.03)
latex.DrawLatex(0.15, 0.85, formula_str)


cCustomKaonProfileX.SaveAs("CustomKaon2DHistogram_ProfileX.png")

#############################################################################################################################
# Set estimate for the number of entries
treeLambda_sw.SetEstimate(treeLambda_sw.GetEntries() + 1)

# Draw the variable to fill the internal array (no plot, just data)
treeLambda_sw.Draw("ProtonMomentum/0.938", "", "goff")

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

hLambdaP_bg = ROOT.TH2F(
    "hLambdaP_bg",
    "hLambdaP_bg",
    m_numPBins,
    binsMinEdgesP_arr,
    m_numDEdxBins,
    0,
    m_dedxCutoff
)

treeLambda_sw.Draw(
    "ProtonSVDdEdx:ProtonMomentum/0.938>>hLambdaP_bg",
    "nSignalLambda_sw * (ProtonSVDdEdx>0)",
    "goff"
)

cProtonCustom = ROOT.TCanvas("cProtonCustom", "Custom Binning Proton 2D", 800, 600)
hLambdaP_bg.Draw("COLZ")
cProtonCustom.SaveAs("ProtonCustomBinning.png")


cProtonCustomProfileX = ROOT.TCanvas("cProtonCustomProfileX", "Custom Binning Proton 2D ProfileX", 800, 600)

cProtonCustomProfileX.SetRightMargin(0.15)
cProtonCustomProfileX.SetLeftMargin(0.12)
cProtonCustomProfileX.SetLogy(True)
cProtonCustomProfileX.SetLogx(True)

CustomprofileX_Proton = hLambdaP_bg.ProfileX("CustomprofileX_Proton")
CustomprofileX_Proton.SetLineColor(ROOT.kRed)
CustomprofileX_Proton.SetLineWidth(2)
CustomprofileX_Proton.SetTitle("ProtonProfileX")
CustomprofileX_Proton.GetYaxis().SetTitle("Mean Energy Loss")
CustomprofileX_Proton.GetXaxis().SetTitle("P/M (GeV/C)")

CustomprofileX_Proton.Draw("E")
CustomprofileX_Proton.SetStats(False)

CustomFitProton = ROOT.TF1(
    "CustomFitProton",
    "[0] + [1]/pow((x + [2]), [3]) + [4]*x",
    0.25, 24
)
CustomFitProton.SetParameters(100, 5000, 0.01, 12, 10)  # [6] is the slope of the rising tail

CustomprofileX_Proton.Fit(CustomFitProton, "", "", 0.25, 3)

CustomFitProton.SetLineColor(ROOT.kBlue)
CustomFitProton.SetLineWidth(2)
CustomFitProton.Draw("SAME")

# Get parameter values
p0 = CustomFitProton.GetParameter(0)
p1 = CustomFitProton.GetParameter(1)
p2 = CustomFitProton.GetParameter(2)
p3 = CustomFitProton.GetParameter(3)
p4 = CustomFitProton.GetParameter(4)


# # Format the formula with parameter values
formula_str = (
    f"Fit: {p0:.2e} + {p1:.2e}/(x + {p2:.2e})^{p3:.2e} + {p4:.2e}*x"
)

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextSize(0.03)
latex.DrawLatex(0.15, 0.85, formula_str)

cProtonCustomProfileX.SaveAs("CustomBinningProton2DHistogram_ProfileX.png")


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

# Define the histogram with custom binning
# m_numDEdxBins = 500  # or your preferred value
# m_dedxCutoff = 5e6   # or your preferred value

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

# Define the histogram with custom binning
# m_numDEdxBins = 500  
# m_dedxCutoff = 5e6 

# Convert to C array for ROOT
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

# m_numDEdxBins = 500
# m_dedxCutoff = 5e6

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
# hTotal.SetStats(False)
hTotal.SetTitle("Combined Energy vs Momentum Histogram")
hTotal.Draw("COLZ")
cCombined.SaveAs("Combined2DHistogram.png")




cCombinedCustomProjectionY = ROOT.TCanvas("cCombinedCustomProjectionY", "Custom Binning Combined 2D ProjectionY", 800, 600)

cCombinedCustomProjectionY.SetRightMargin(0.15)
cCombinedCustomProjectionY.SetLeftMargin(0.12)
# cCombinedCustomProjectionY.SetLogy(True)
# cCombinedCustomProjectionY.SetLogx(True)


CustomProjectionY_Combined = hTotal.ProjectionY("CustomProjectionY_Combined",50,100)
CustomProjectionY_Combined.SetLineColor(ROOT.kRed)
CustomProjectionY_Combined.SetLineWidth(2)
CustomProjectionY_Combined.SetTitle("Combined Pion ProjectionY")
CustomProjectionY_Combined.GetYaxis().SetTitle("Mean Energy Loss")
CustomProjectionY_Combined.GetXaxis().SetTitle("P/M (GeV/C)")
ROOT.gStyle.SetStatW(0.5)
ROOT.gStyle.SetStatH(0.5)
CustomProjectionY_Combined.Draw("E")
# CustomProjectionY_Combined.SetStats(False)


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

print("Gaussian 1: mean =", mean1, ", sigma =", sigma1)
print("Gaussian 2: sigma =", sigma2)

# # Format the formula with parameter values
# formula_str = (
#     f"Gauss1: Mean = {mean1:.2e} Sigma = {sigma1:.2e} \n Gauss2: Mean = {mean2:.2e} Sigma = {sigma2:.2e}"
# )

# latex = ROOT.TLatex()
# latex.SetNDC()
# latex.SetTextSize(0.03)
# latex.DrawLatex(0.15, 0.85, formula_str)


cCombinedCustomProjectionY.SaveAs("CustomBinningCombined2DHistogram_ProjectionY.png")




cCombinedCustomProfileX = ROOT.TCanvas("cCombinedCustomProfileX", "Custom Binning Combined 2D ProfileX", 800, 600)

cCombinedCustomProfileX.SetRightMargin(0.15)
cCombinedCustomProfileX.SetLeftMargin(0.12)
cCombinedCustomProfileX.SetLogy(True)
cCombinedCustomProfileX.SetLogx(True)

CustomprofileX_Combined = hTotal.ProfileX("CustomprofileX_Combined")
CustomprofileX_Combined.SetLineColor(ROOT.kRed)
CustomprofileX_Combined.SetLineWidth(2)
CustomprofileX_Combined.SetTitle("Combined Pion ProfileX")
CustomprofileX_Combined.GetYaxis().SetTitle("Mean Energy Loss")
CustomprofileX_Combined.GetXaxis().SetTitle("P/M (GeV/C)")

CustomprofileX_Combined.Draw("E")
CustomprofileX_Combined.SetStats(False)

CustomFitCombined = ROOT.TF1(
    "CustomFitCombined",
    "[0] + [1]/pow((x + [2]), [3]) + [4]*x",
    0.4, 24
)
CustomFitCombined.SetParameters(100, 5000, 0.01, 12, 10)  # [6] is the slope of the rising tail

CustomprofileX_Combined.Fit(CustomFitCombined, "", "", 0.4, 24)

CustomFitCombined.SetLineColor(ROOT.kBlue)
CustomFitCombined.SetLineWidth(2)
CustomFitCombined.Draw("SAME")

# Get parameter values
p0 = CustomFitCombined.GetParameter(0)
p1 = CustomFitCombined.GetParameter(1)
p2 = CustomFitCombined.GetParameter(2)
p3 = CustomFitCombined.GetParameter(3)
p4 = CustomFitCombined.GetParameter(4)


# # Format the formula with parameter values
formula_str = (
    f"Fit: {p0:.2e} + {p1:.2e}/(x + {p2:.2e})^{p3:.2e} + {p4:.2e}*x"
)

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextSize(0.03)
latex.DrawLatex(0.15, 0.85, formula_str)

cCombinedCustomProfileX.SaveAs("CustomBinningCombined2DHistogram_ProfileX.png")


#####################################################################

cProfileFormulae = ROOT.TCanvas("cProfileFormulae", "ProfileX of Pions, Kaons, and Protons", 800, 600)
cProfileFormulae.SetRightMargin(0.15)
cProfileFormulae.SetLeftMargin(0.12)
cProfileFormulae.SetLogy(True)
cProfileFormulae.SetLogx(True)

CustomFitCombined.Draw()
FitKaon.Draw('Same')
CustomFitProton.Draw('Same')
CustomFitCombined.SetLineColor(ROOT.kRed)
CustomFitCombined.SetLineWidth(2)
FitKaon.SetLineColor(ROOT.kBlue)
FitKaon.SetLineWidth(2)
CustomFitProton.SetLineColor(ROOT.kBlack)
CustomFitProton.SetLineWidth(2)


legend = ROOT.TLegend(0.7 ,0.6 ,0.85 ,0.75)
legend.AddEntry(CustomFitCombined,"Fit Pions")
legend.AddEntry(FitKaon,"Fit Kaons")
legend.AddEntry(CustomFitProton,"Fit Protons")
legend.SetLineWidth(0)
legend.Draw("Same")

CustomFitCombined.GetYaxis().SetRangeUser(1e6, 1e7)
CustomFitCombined.SetTitle("Momentum vs Mean Energy Loss")
CustomFitCombined.GetYaxis().SetTitle("Mean Energy Loss")
CustomFitCombined.GetXaxis().SetTitle("P/M (GeV/C)")

cProfileFormulae.SaveAs("CustomCombinedFormulae.png")

########################################################################################

cDataCustom = ROOT.TCanvas("cDataCustom", "Custom Binning ProfileX: Pions, Kaons, Protons", 800, 600)
cDataCustom.SetRightMargin(0.15)
cDataCustom.SetLeftMargin(0.12)
cDataCustom.SetLogy(True)
cDataCustom.SetLogx(True)

CustomprofileX_Combined.SetLineColor(ROOT.kRed)
CustomprofileX_Combined.SetLineWidth(2)
CustomprofileX_Combined.SetTitle("ProfileX: Pions (red), Kaons (blue), Protons (black)")
CustomprofileX_Combined.GetYaxis().SetTitle("Mean Energy Loss")
CustomprofileX_Combined.GetXaxis().SetTitle("P/M (GeV/C)")

CustomprofileX_Combined.Draw("E")
CustomprofileX_Kaon.SetLineColor(ROOT.kBlue)
CustomprofileX_Kaon.SetLineWidth(2)
CustomprofileX_Kaon.Draw("E SAME")
CustomprofileX_Proton.SetLineColor(ROOT.kBlack)
CustomprofileX_Proton.SetLineWidth(2)
CustomprofileX_Proton.Draw("E SAME")

legend = ROOT.TLegend(0.7, 0.6, 0.85, 0.75)
legend.AddEntry(CustomprofileX_Combined, "Pions")
legend.AddEntry(CustomprofileX_Kaon, "Kaons")
legend.AddEntry(CustomprofileX_Proton, "Protons")
legend.SetLineWidth(0)
legend.Draw("Same")

cDataCustom.SaveAs("DataCustomBinning.png")
