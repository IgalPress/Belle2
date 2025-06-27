import ROOT
from ROOT import TCanvas
from array import array
import numpy as np
# import uproot

latex = ROOT . TLatex ()
latex . SetNDC ()
latex . SetTextSize (0.03)


LambdaFile = ROOT.TFile("Lambda_8Apr.root")
DstarFile = ROOT.TFile("Dstar_8Apr.root")
GammaFile = ROOT.TFile("Gamma_8Apr.root")
GammaMeansFile = ROOT.TFile("Gamma_8Apr_withCMSmean.root")

LambdaTree = LambdaFile.Get("Lambda")
DstarTree = DstarFile.Get("Dstar")
GammaTree = GammaFile.Get("Gamma")
GammaMeansTree = GammaMeansFile.Get("Gamma_Means;2")

Histogram2D_Kaon = ROOT.TH2D("Histogram2D_Kaon", "Kaon Energy vs Momentum;Energy (arb units);Momentum (GeV/c)",  100, 0,20,100,0, 5e6)
Histogram2D_Proton = ROOT.TH2D("Histogram2D_Proton", "Proton Energy vs Momentum;Energy (arb units);Momentum (GeV/c)",  100, 0,20,100,0, 5e6)
Histogram2D_Lambda = ROOT.TH2D("Histogram2D_Lambda", "Lambda Pion Energy vs Momentum;Energy (arb units);Momentum  (GeV/c)", 100, 0,20,100,0, 5e6)
Histogram2D_DstarPion = ROOT.TH2D("Histogram2D_DstarPion", "D* Pion Energy vs Momentum;Energy (arb units);Momentum (GeV/c)",  100, 0,20,100,0, 5e6)
Histogram2D_DstarSlowPion = ROOT.TH2D("Histogram2D_DstarSlowPion", "Slow Pion Energy vs Momentum;Energy (arb units);Momentum (GeV/c)", 100, 0,20,100,0, 5e6)
Histogram2D_FirstElectron = ROOT.TH2D("Histogram2D_FirstElectron", "First Electron Energy vs Momentum;Energy (arb units);Momentum (GeV/c)", 100, 0,20,100,0, 5e6)
Histogram2D_SecondElectron = ROOT.TH2D("Histogram2D_SecondElectron", "Second Electron Energy vs Momentum;Energy (arb units);Momentum (GeV/c)", 100, 0,20,100,0, 5e6)
Histogram2D_SecondElectronNHitsUsed = ROOT.TH2D("Histogram2D_SecondElectronNHitsUsed", "Second Electron Energy vs NHitsUsed;Energy (arb units);NHitsUsed (arb units)", 100, 0,20,100,0, 5e6)
Histogram2D_FirstElectronNHitsUsed = ROOT.TH2D("Histogram2D_FirstElectronNHitsUsed", "First Electron Energy vs NHitsUsed;Energy (arb units);NHitsUsed (arb units)", 100, 0,20,100,0, 5e6)


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

file3 = ROOT.TFile.Open("Gamma_8Apr_withCMSmean.root")
preselTreeGamma_Means = file3.Get("Gamma_Means;2")
if not preselTreeGamma_Means:
    raise RuntimeError("Failed to get 'Gamma_Means;2' tree from Gamma_8Apr_withCMSmean.root")
print('Files Loaded')




FirstElectronSVDdEdxList = np.empty(20, "float64")
# FirstElectronSVDdEdx = 0

preselTreeGamma.Branch("FirstElectronSVDdEdxList",FirstElectronSVDdEdxList,"FirstElectronSVDdEdxList[20]/D")
# preselTreeGamma.Branch("FirstElectronSVDdEdx",FirstElectronSVDdEdx,"FirstElectronSVDdEdx/D")




ROOT.gStyle.SetOptFit(1)

variables = ROOT.RooArgSet()
variables_Dstar = ROOT.RooArgSet()
variables_Gamma = ROOT.RooArgSet()
variables_GammaMean = ROOT.RooArgSet()
variables.add(InvM)
variables_Dstar.add(InvM_Dstar)


KaonMomentum = ROOT.RooRealVar("KaonMomentum", "", -1e8, 1e8)
KaonSVDdEdx = ROOT.RooRealVar("KaonSVDdEdx", "", -1e8, 1e8)
KaonSVDdEdxTrackNHitsUsed = ROOT.RooRealVar("KaonSVDdEdxTrackNHitsUsed","", -1e8, 1e8)
PionDMomentum = ROOT.RooRealVar("PionDMomentum", "", -1e8, 1e8)
PionDSVDdEdx = ROOT.RooRealVar("PionDSVDdEdx", "", -1e8, 1e8)
PionDSVDdEdxTrackNHitsUsed = ROOT.RooRealVar("PionDSVDdEdxTrackNHitsUsed", "", -1e8, 1e8)
SlowPionMomentum = ROOT.RooRealVar("SlowPionMomentum", "", -1e8, 1e8)
SlowPionSVDdEdx = ROOT.RooRealVar("SlowPionSVDdEdx", "", -1e8, 1e8)
SlowPionSVDdEdxTrackNHitsUsed = ROOT.RooRealVar("SlowPionSVDdEdxTrackNHitsUsed", "", -1e8, 1e8)
event_Dstar = ROOT.RooRealVar("event","",-1e8,1e8)

ProtonMomentum = ROOT.RooRealVar("ProtonMomentum", "momentum for Proton (GeV)", -1e8, 1e8)
ProtonSVDdEdxTrackMomentum = ROOT.RooRealVar("ProtonSVDdEdxTrackMomentum", "momentum for Proton (GeV), from the track", -1e8, 1e8)
ProtonSVDdEdx = ROOT.RooRealVar("ProtonSVDdEdx", "", -1e8, 1e8)
PionLambdaMomentum = ROOT.RooRealVar("PionLambdaMomentum", "momentum for pion (GeV)", -1e8, 1e8)
PionLambdaSVDdEdxTrackMomentum = ROOT.RooRealVar("PionLambdaSVDdEdxTrackMomentum", "momentum for pion (GeV), from the track", -1e8, 1e8)
PionLambdaSVDdEdx = ROOT.RooRealVar("PionLambdaSVDdEdx", "", -1e8, 1e8)
PionLambdaSVDdEdxTrackNHitsUsed = ROOT.RooRealVar("PionLambdaSVDdEdxTrackNHitsUsed", "", -1e8, 1e8)
event_Lambda = ROOT.RooRealVar("event","",-1e8,1e8)

FirstElectronMomentum = ROOT.RooRealVar("FirstElectronMomentum", "", -1e8, 1e8)
FirstElectronSVDdEdx = ROOT.RooRealVar("FirstElectronSVDdEdx", "", -1e8, 1e8)
FirstElectronSVDdEdxTrackNHitsUsed = ROOT.RooRealVar("FirstElectronSVDdEdxTrackNHitsUsed", "", -1e8, 1e8)
FirstElectronSVDdEdxList = ROOT.RooRealVar("FirstElectronSVDdEdxList", "", -1e8, 1e8)
SecondElectronMomentum = ROOT.RooRealVar("SecondElectronMomentum", "", -1e8, 1e8)
SecondElectronSVDdEdx = ROOT.RooRealVar("SecondElectronSVDdEdx", "", -1e8, 1e8)
SecondElectronSVDdEdxTrackNHitsUsed = ROOT.RooRealVar("SecondElectronSVDdEdxTrackNHitsUsed", "", -1e8, 1e8)

CMS_mean = ROOT.RooRealVar("CMS_mean","",-1e8,1e8)
ALICE_mean = ROOT.RooRealVar("ALICE_mean","",-1e8,1e8)
ATLAS_mean = ROOT.RooRealVar("ATLAS_mean","",-1e8,1e8)
Harmonic_mean = ROOT.RooRealVar("Harmonic_mean","",-1e8,1e8)
filtered_FirstElectronSVDdEdxArray = ROOT.RooRealVar("filtered_FirstElectronSVDdEdxArray","",-1e8,1e8)
filtered_FirstElectronMomentumArray = ROOT.RooRealVar("filtered_FirstElectronMomentumArray","",-1e8,1e8)

variables_Dstar.add(KaonMomentum)
variables_Dstar.add(KaonSVDdEdxTrackNHitsUsed)
variables_Dstar.add(PionDMomentum)
variables_Dstar.add(SlowPionMomentum)
variables_Dstar.add(PionDSVDdEdx)
variables_Dstar.add(PionDSVDdEdxTrackNHitsUsed)
variables_Dstar.add(SlowPionSVDdEdx)
variables_Dstar.add(SlowPionSVDdEdxTrackNHitsUsed)
variables_Dstar.add(KaonSVDdEdx)
variables_Dstar.add(event_Dstar)
variables.add(ProtonMomentum)
variables.add(ProtonSVDdEdxTrackMomentum)
variables.add(ProtonSVDdEdx)
variables.add(PionLambdaMomentum)
variables.add(PionLambdaSVDdEdxTrackMomentum)
variables.add(PionLambdaSVDdEdx)
variables.add(PionLambdaSVDdEdxTrackNHitsUsed)
variables.add(event_Lambda)
variables_Gamma.add(FirstElectronMomentum)
variables_Gamma.add(FirstElectronSVDdEdx)
variables_Gamma.add(FirstElectronSVDdEdxTrackNHitsUsed)
variables_Gamma.add(FirstElectronSVDdEdxList)
variables_Gamma.add(SecondElectronMomentum)
variables_Gamma.add(SecondElectronSVDdEdx)
variables_Gamma.add(SecondElectronSVDdEdxTrackNHitsUsed)

variables_GammaMean.add(CMS_mean)
variables_GammaMean.add(ALICE_mean)
variables_GammaMean.add(ATLAS_mean)
variables_GammaMean.add(Harmonic_mean)
variables_GammaMean.add(filtered_FirstElectronSVDdEdxArray)
variables_GammaMean.add(filtered_FirstElectronMomentumArray)

LambdaDataset = ROOT.RooDataSet("LambdaDataset", "LambdaDataset", preselTree, variables, "ProtonMomentum > 0.25")
DstarDataset = ROOT.RooDataSet("DstarDataset", "DstarDataset", preselTreeDstar, variables_Dstar)
GammaDataset = ROOT.RooDataSet("GammaDataset", "GammaDataset", preselTreeGamma, variables_Gamma)
GammaMeansDataset = ROOT.RooDataSet("GammaMeansDataset", "GammaMeansDataset", preselTreeGamma_Means, variables_GammaMean)


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

GammaMeansDatasetSWeighted = ROOT.RooDataSet(GammaMeansDataset.GetName(), GammaMeansDataset.GetTitle(), GammaMeansDataset, GammaMeansDataset.get()) 

treeGammaMeans_sw = GammaMeansDatasetSWeighted.GetClonedTree()
treeGammaMeans_sw.SetName("treeGammaMeans_sw") 

bins_start = 0
bins_fin = 200

m_numDEdxBins = 200  
m_dedxCutoff = 5e6   

hCMS_mean_bg1 = ROOT.TH1F(
    "hCMS_mean",
    "hCMS_mean",
    m_numDEdxBins,
    0,
    m_dedxCutoff
)


preselTreeGamma_Means.Draw(
    "CMS_mean>>hCMS_mean",
    "(filtered_FirstElectronSVDdEdxArray>0)",
    "goff"
)
#NOTE THAT IF YOU HAVE THE FILTERED_ > 0 CUT IT WILL DROP THE TOTAL NUMBER OF EVENTS BY A BIT.



cCMS_meanCustom = ROOT.TCanvas("cCMS_meanCustom", "Custom Binning CMS_mean", 800, 600)
hCMS_mean_bg1.Draw("COLZ")
cCMS_meanCustom.SaveAs("MeanChecking/CMS_Mean_dEdx.png")



treeGammaMeans_sw.SetEstimate(treeGammaMeans_sw.GetEntries() + 1)


treeGammaMeans_sw.Draw("filtered_FirstElectronMomentumArray/0.0005110", "", "goff")


vXP = treeGammaMeans_sw.GetV1()


m_numPBins = 100

kdBinsP = ROOT.TKDTreeBinning(treeGammaMeans_sw.GetEntries(), 1, vXP, m_numPBins)
binsMinEdgesP_orig = kdBinsP.SortOneDimBinEdges()

binsMinEdgesP = [binsMinEdgesP_orig[i] for i in range(m_numPBins + 2)]

binsMinEdgesP[0] = 0.1
binsMinEdgesP[m_numPBins + 1] = 50.



binsMinEdgesP_arr = array('d', binsMinEdgesP)

hCMS_meanVsMom_bg0 = ROOT.TH2F(
    "hCMS_meanVsMomCustom",
    "hCMS_meanVsMomCustom",
    m_numPBins,
    binsMinEdgesP_arr,
    m_numDEdxBins,
    0,
    m_dedxCutoff
)

preselTreeGamma_Means.Draw(
    "CMS_mean:filtered_FirstElectronMomentumArray/0.0005110>>hCMS_meanVsMomCustom",
    "(filtered_FirstElectronSVDdEdxArray>0)",
    "goff"
)



cCMS_meanVsMom = ROOT.TCanvas("cCMS_meanVsMom", "Custom Binning FirstElectron", 800, 600)
hCMS_meanVsMom_bg0.Draw("COLZ")
cCMS_meanVsMom.SaveAs("MeanChecking/FirstElectronCMS_meanVsMomentumCustom.png")




cCMS_meanVsMomCustomProjectionY = ROOT.TCanvas("cCMS_meanVsMomCustomProjectionY", "Custom Binning CMS_mean 2D ProjectionY", 800, 600)



cCMS_meanVsMomCustomProjectionY.SetRightMargin(0.15)
cCMS_meanVsMomCustomProjectionY.SetLeftMargin(0.12)
# cFirstElectronCustomProjectionY.SetLogy(True)
# cFirstElectronCustomProjectionY.SetLogx(True)


CustomProjectionY_CMS_Mean = hCMS_meanVsMom_bg0.ProjectionY("CustomProjectionY_CMS_Mean", bins_start, bins_fin)
CustomProjectionY_CMS_Mean.SetLineColor(ROOT.kRed)
CustomProjectionY_CMS_Mean.SetLineWidth(2)
CustomProjectionY_CMS_Mean.SetTitle("CMS ProjectionY")
CustomProjectionY_CMS_Mean.GetYaxis().SetTitle("Mean Energy Loss")
CustomProjectionY_CMS_Mean.GetXaxis().SetTitle("P/M (GeV/C)")
CustomProjectionY_CMS_Mean.Draw("E1 SAME")

CustomProjectionY_CMS_Mean.SetStats(False)


gaus2Landau = ROOT.TF1(
    "gaus2Landau",
    "[0]*TMath::Gaus(x, [1], [1]*[2]*[5]) + [3]*TMath::Gaus(x, [1], [1]*[4]*[5]) + [6]*TMath::Landau(x, [1], [1]*[5])",
    1.e5, 1.5e6
)

gaus2Landau.SetLineColor(ROOT.kBlack)
# Parameters: gaus1_amp, gaus1_mean, gaus1_sigma, gaus2_amp, gaus2_sigma, landau_width, Landau_Amplitude
gaus2Landau.SetParameters(20000, 650e3, 1, 4000, 1, 0.1, 1e5)
gaus2Landau.SetNpx(1000)

gaus2Landau.FixParameter(4,2.1)
gaus2Landau.FixParameter(2,2.7)
gaus2Landau.SetParLimits(5,0,1)
gaus2Landau.SetParLimits(2,1,5)
gaus2Landau.SetParLimits(4,1,5)

CustomProjectionY_CMS_Mean.Fit("gaus2Landau", "R SAME")
LandauWidth0 = gaus2Landau.GetParameter(5)
LandauWidthErr0 = gaus2Landau.GetParError(5)

latex.DrawLatex(0.4,0.8, "L. width NHitsUsed>0: %.4f +- %.4f" % (LandauWidth0, LandauWidthErr0))
# p5FirstElectron = gaus2Landau.GetParameter(5)
# p5FirstElectronErr = gaus2Landau.GetParError(5)
# p1FirstElectron = gaus2Landau.GetParameter(1)


cCMS_meanVsMomCustomProjectionY.SaveAs("MeanChecking/CMS_meanProjectionY.png")


################################################################



hATLAS_mean_bg1 = ROOT.TH1F(
    "hATLAS_mean",
    "hATLAS_mean",
    m_numDEdxBins,
    0,
    m_dedxCutoff
)


preselTreeGamma_Means.Draw(
    "ATLAS_mean>>hATLAS_mean",
    "(filtered_FirstElectronSVDdEdxArray>0)",
    "goff"
)
#NOTE THAT IF YOU HAVE THE FILTERED_ > 0 CUT IT WILL DROP THE TOTAL NUMBER OF EVENTS BY A BIT.



cATLAS_meanCustom = ROOT.TCanvas("cATLAS_meanCustom", "Custom Binning ATLAS_mean", 800, 600)
hATLAS_mean_bg1.Draw("COLZ")
cATLAS_meanCustom.SaveAs("MeanChecking/ATLAS_Mean_dEdx.png")



treeGammaMeans_sw.SetEstimate(treeGammaMeans_sw.GetEntries() + 1)


treeGammaMeans_sw.Draw("filtered_FirstElectronMomentumArray/0.0005110", "", "goff")


vXP = treeGammaMeans_sw.GetV1()


m_numPBins = 100

kdBinsP = ROOT.TKDTreeBinning(treeGammaMeans_sw.GetEntries(), 1, vXP, m_numPBins)
binsMinEdgesP_orig = kdBinsP.SortOneDimBinEdges()

binsMinEdgesP = [binsMinEdgesP_orig[i] for i in range(m_numPBins + 2)]

binsMinEdgesP[0] = 0.1
binsMinEdgesP[m_numPBins + 1] = 50.



binsMinEdgesP_arr = array('d', binsMinEdgesP)

hATLAS_meanVsMom_bg0 = ROOT.TH2F(
    "hATLAS_meanVsMomCustom",
    "hATLAS_meanVsMomCustom",
    m_numPBins,
    binsMinEdgesP_arr,
    m_numDEdxBins,
    0,
    m_dedxCutoff
)

preselTreeGamma_Means.Draw(
    "ATLAS_mean:filtered_FirstElectronMomentumArray/0.0005110>>hATLAS_meanVsMomCustom",
    "(filtered_FirstElectronSVDdEdxArray>0)",
    "goff"
)



cATLAS_meanVsMom = ROOT.TCanvas("cATLAS_meanVsMom", "Custom Binning FirstElectron", 800, 600)
hATLAS_meanVsMom_bg0.Draw("COLZ")
cATLAS_meanVsMom.SaveAs("MeanChecking/FirstElectronATLAS_meanVsMomentumCustom.png")




cATLAS_meanVsMomCustomProjectionY = ROOT.TCanvas("cATLAS_meanVsMomCustomProjectionY", "Custom Binning ATLAS_mean 2D ProjectionY", 800, 600)



cATLAS_meanVsMomCustomProjectionY.SetRightMargin(0.15)
cATLAS_meanVsMomCustomProjectionY.SetLeftMargin(0.12)
# cFirstElectronCustomProjectionY.SetLogy(True)
# cFirstElectronCustomProjectionY.SetLogx(True)


CustomProjectionY_ATLAS_Mean = hATLAS_meanVsMom_bg0.ProjectionY("CustomProjectionY_ATLAS_Mean", bins_start, bins_fin)
CustomProjectionY_ATLAS_Mean.SetLineColor(ROOT.kRed)
CustomProjectionY_ATLAS_Mean.SetLineWidth(2)
CustomProjectionY_ATLAS_Mean.SetTitle("ATLAS ProjectionY")
CustomProjectionY_ATLAS_Mean.GetYaxis().SetTitle("Mean Energy Loss")
CustomProjectionY_ATLAS_Mean.GetXaxis().SetTitle("P/M (GeV/C)")
CustomProjectionY_ATLAS_Mean.Draw("E1 SAME")

CustomProjectionY_ATLAS_Mean.SetStats(False)


gaus2Landau = ROOT.TF1(
    "gaus2Landau",
    "[0]*TMath::Gaus(x, [1], [1]*[2]*[5]) + [3]*TMath::Gaus(x, [1], [1]*[4]*[5]) + [6]*TMath::Landau(x, [1], [1]*[5])",
    1.e5, 1.5e6
)

gaus2Landau.SetLineColor(ROOT.kBlack)
# Parameters: gaus1_amp, gaus1_mean, gaus1_sigma, gaus2_amp, gaus2_sigma, landau_width, Landau_Amplitude
gaus2Landau.SetParameters(20000, 650e3, 1, 4000, 1, 0.1, 1e5)
gaus2Landau.SetNpx(1000)

gaus2Landau.FixParameter(4,2.1)
gaus2Landau.FixParameter(2,2.7)
gaus2Landau.SetParLimits(5,0,1)
gaus2Landau.SetParLimits(2,1,5)
gaus2Landau.SetParLimits(4,1,5)

CustomProjectionY_ATLAS_Mean.Fit("gaus2Landau", "R SAME")
LandauWidth0 = gaus2Landau.GetParameter(5)
LandauWidthErr0 = gaus2Landau.GetParError(5)

latex.DrawLatex(0.4,0.8, "L. width NHitsUsed>0: %.4f +- %.4f" % (LandauWidth0, LandauWidthErr0))
# p5FirstElectron = gaus2Landau.GetParameter(5)
# p5FirstElectronErr = gaus2Landau.GetParError(5)
# p1FirstElectron = gaus2Landau.GetParameter(1)


cATLAS_meanVsMomCustomProjectionY.SaveAs("MeanChecking/ATLAS_meanProjectionY.png")




####################################################################################



hALICE_mean_bg1 = ROOT.TH1F(
    "hALICE_mean",
    "hALICE_mean",
    m_numDEdxBins,
    0,
    m_dedxCutoff
)


preselTreeGamma_Means.Draw(
    "ALICE_mean>>hALICE_mean",
    "(filtered_FirstElectronSVDdEdxArray>0)",
    "goff"
)
#NOTE THAT IF YOU HAVE THE FILTERED_ > 0 CUT IT WILL DROP THE TOTAL NUMBER OF EVENTS BY A BIT.



cALICE_meanCustom = ROOT.TCanvas("cALICE_meanCustom", "Custom Binning ALICE_mean", 800, 600)
hALICE_mean_bg1.Draw("COLZ")
cALICE_meanCustom.SaveAs("MeanChecking/ALICE_Mean_dEdx.png")



treeGammaMeans_sw.SetEstimate(treeGammaMeans_sw.GetEntries() + 1)


treeGammaMeans_sw.Draw("filtered_FirstElectronMomentumArray/0.0005110", "", "goff")


vXP = treeGammaMeans_sw.GetV1()


m_numPBins = 100

kdBinsP = ROOT.TKDTreeBinning(treeGammaMeans_sw.GetEntries(), 1, vXP, m_numPBins)
binsMinEdgesP_orig = kdBinsP.SortOneDimBinEdges()

binsMinEdgesP = [binsMinEdgesP_orig[i] for i in range(m_numPBins + 2)]

binsMinEdgesP[0] = 0.1
binsMinEdgesP[m_numPBins + 1] = 50.



binsMinEdgesP_arr = array('d', binsMinEdgesP)

hALICE_meanVsMom_bg0 = ROOT.TH2F(
    "hALICE_meanVsMomCustom",
    "hALICE_meanVsMomCustom",
    m_numPBins,
    binsMinEdgesP_arr,
    m_numDEdxBins,
    0,
    m_dedxCutoff
)

preselTreeGamma_Means.Draw(
    "ALICE_mean:filtered_FirstElectronMomentumArray/0.0005110>>hALICE_meanVsMomCustom",
    "(filtered_FirstElectronSVDdEdxArray>0)",
    "goff"
)



cALICE_meanVsMom = ROOT.TCanvas("cALICE_meanVsMom", "Custom Binning FirstElectron", 800, 600)
hALICE_meanVsMom_bg0.Draw("COLZ")
cALICE_meanVsMom.SaveAs("MeanChecking/FirstElectronALICE_meanVsMomentumCustom.png")




cALICE_meanVsMomCustomProjectionY = ROOT.TCanvas("cALICE_meanVsMomCustomProjectionY", "Custom Binning ALICE_mean 2D ProjectionY", 800, 600)



cALICE_meanVsMomCustomProjectionY.SetRightMargin(0.15)
cALICE_meanVsMomCustomProjectionY.SetLeftMargin(0.12)
# cFirstElectronCustomProjectionY.SetLogy(True)
# cFirstElectronCustomProjectionY.SetLogx(True)


CustomProjectionY_ALICE_Mean = hALICE_meanVsMom_bg0.ProjectionY("CustomProjectionY_ALICE_Mean", bins_start, bins_fin)
CustomProjectionY_ALICE_Mean.SetLineColor(ROOT.kRed)
CustomProjectionY_ALICE_Mean.SetLineWidth(2)
CustomProjectionY_ALICE_Mean.SetTitle("ALICE ProjectionY")
CustomProjectionY_ALICE_Mean.GetYaxis().SetTitle("Mean Energy Loss")
CustomProjectionY_ALICE_Mean.GetXaxis().SetTitle("P/M (GeV/C)")
CustomProjectionY_ALICE_Mean.Draw("E1 SAME")

CustomProjectionY_ALICE_Mean.SetStats(False)


gaus2Landau = ROOT.TF1(
    "gaus2Landau",
    "[0]*TMath::Gaus(x, [1], [1]*[2]*[5]) + [3]*TMath::Gaus(x, [1], [1]*[4]*[5]) + [6]*TMath::Landau(x, [1], [1]*[5])",
    1.e5, 1.5e6
)

gaus2Landau.SetLineColor(ROOT.kBlack)
# Parameters: gaus1_amp, gaus1_mean, gaus1_sigma, gaus2_amp, gaus2_sigma, landau_width, Landau_Amplitude
gaus2Landau.SetParameters(20000, 650e3, 1, 4000, 1, 0.1, 1e5)
gaus2Landau.SetNpx(1000)

gaus2Landau.FixParameter(4,2.1)
gaus2Landau.FixParameter(2,2.7)
gaus2Landau.SetParLimits(5,0,1)
gaus2Landau.SetParLimits(2,1,5)
gaus2Landau.SetParLimits(4,1,5)

CustomProjectionY_ALICE_Mean.Fit("gaus2Landau", "R SAME")
LandauWidth0 = gaus2Landau.GetParameter(5)
LandauWidthErr0 = gaus2Landau.GetParError(5)

latex.DrawLatex(0.4,0.8, "L. width NHitsUsed>0: %.4f +- %.4f" % (LandauWidth0, LandauWidthErr0))
# p5FirstElectron = gaus2Landau.GetParameter(5)
# p5FirstElectronErr = gaus2Landau.GetParError(5)
# p1FirstElectron = gaus2Landau.GetParameter(1)


cALICE_meanVsMomCustomProjectionY.SaveAs("MeanChecking/ALICE_meanProjectionY.png")


###################################################################################



hHarmonic_mean_bg1 = ROOT.TH1F(
    "hHarmonic_mean",
    "hHarmonic_mean",
    m_numDEdxBins,
    0,
    m_dedxCutoff
)


preselTreeGamma_Means.Draw(
    "Harmonic_mean>>hHarmonic_mean",
    "(filtered_FirstElectronSVDdEdxArray>0)",
    "goff"
)
#NOTE THAT IF YOU HAVE THE FILTERED_ > 0 CUT IT WILL DROP THE TOTAL NUMBER OF EVENTS BY A BIT.



cHarmonic_meanCustom = ROOT.TCanvas("cHarmonic_meanCustom", "Custom Binning Harmonic_mean", 800, 600)
hHarmonic_mean_bg1.Draw("COLZ")
cHarmonic_meanCustom.SaveAs("MeanChecking/Harmonic_Mean_dEdx.png")



treeGammaMeans_sw.SetEstimate(treeGammaMeans_sw.GetEntries() + 1)


treeGammaMeans_sw.Draw("filtered_FirstElectronMomentumArray/0.0005110", "", "goff")


vXP = treeGammaMeans_sw.GetV1()


m_numPBins = 100

kdBinsP = ROOT.TKDTreeBinning(treeGammaMeans_sw.GetEntries(), 1, vXP, m_numPBins)
binsMinEdgesP_orig = kdBinsP.SortOneDimBinEdges()

binsMinEdgesP = [binsMinEdgesP_orig[i] for i in range(m_numPBins + 2)]

binsMinEdgesP[0] = 0.1
binsMinEdgesP[m_numPBins + 1] = 50.



binsMinEdgesP_arr = array('d', binsMinEdgesP)

hHarmonic_meanVsMom_bg0 = ROOT.TH2F(
    "hHarmonic_meanVsMomCustom",
    "hHarmonic_meanVsMomCustom",
    m_numPBins,
    binsMinEdgesP_arr,
    m_numDEdxBins,
    0,
    m_dedxCutoff
)

preselTreeGamma_Means.Draw(
    "Harmonic_mean:filtered_FirstElectronMomentumArray/0.0005110>>hHarmonic_meanVsMomCustom",
    "(filtered_FirstElectronSVDdEdxArray>0)",
    "goff"
)



cHarmonic_meanVsMom = ROOT.TCanvas("cHarmonic_meanVsMom", "Custom Binning FirstElectron", 800, 600)
hHarmonic_meanVsMom_bg0.Draw("COLZ")
cHarmonic_meanVsMom.SaveAs("MeanChecking/FirstElectronHarmonic_meanVsMomentumCustom.png")




cHarmonic_meanVsMomCustomProjectionY = ROOT.TCanvas("cHarmonic_meanVsMomCustomProjectionY", "Custom Binning Harmonic_mean 2D ProjectionY", 800, 600)



cHarmonic_meanVsMomCustomProjectionY.SetRightMargin(0.15)
cHarmonic_meanVsMomCustomProjectionY.SetLeftMargin(0.12)
# cFirstElectronCustomProjectionY.SetLogy(True)
# cFirstElectronCustomProjectionY.SetLogx(True)


CustomProjectionY_Harmonic_Mean = hHarmonic_meanVsMom_bg0.ProjectionY("CustomProjectionY_Harmonic_Mean", bins_start, bins_fin)
CustomProjectionY_Harmonic_Mean.SetLineColor(ROOT.kRed)
CustomProjectionY_Harmonic_Mean.SetLineWidth(2)
CustomProjectionY_Harmonic_Mean.SetTitle("Harmonic ProjectionY")
CustomProjectionY_Harmonic_Mean.GetYaxis().SetTitle("Mean Energy Loss")
CustomProjectionY_Harmonic_Mean.GetXaxis().SetTitle("P/M (GeV/C)")
CustomProjectionY_Harmonic_Mean.Draw("E1 SAME")

CustomProjectionY_Harmonic_Mean.SetStats(False)


gaus2Landau = ROOT.TF1(
    "gaus2Landau",
    "[0]*TMath::Gaus(x, [1], [1]*[2]*[5]) + [3]*TMath::Gaus(x, [1], [1]*[4]*[5]) + [6]*TMath::Landau(x, [1], [1]*[5])",
    1.e5, 1.5e6
)

gaus2Landau.SetLineColor(ROOT.kBlack)
# Parameters: gaus1_amp, gaus1_mean, gaus1_sigma, gaus2_amp, gaus2_sigma, landau_width, Landau_Amplitude
gaus2Landau.SetParameters(20000, 650e3, 1, 4000, 1, 0.1, 1e5)
gaus2Landau.SetNpx(1000)

gaus2Landau.FixParameter(4,2.1)
gaus2Landau.FixParameter(2,2.7)
gaus2Landau.SetParLimits(5,0,1)
gaus2Landau.SetParLimits(2,1,5)
gaus2Landau.SetParLimits(4,1,5)

CustomProjectionY_Harmonic_Mean.Fit("gaus2Landau", "R SAME")
LandauWidth0 = gaus2Landau.GetParameter(5)
LandauWidthErr0 = gaus2Landau.GetParError(5)

latex.DrawLatex(0.4,0.8, "L. width NHitsUsed>0: %.4f +- %.4f" % (LandauWidth0, LandauWidthErr0))
# p5FirstElectron = gaus2Landau.GetParameter(5)
# p5FirstElectronErr = gaus2Landau.GetParError(5)
# p1FirstElectron = gaus2Landau.GetParameter(1)


cHarmonic_meanVsMomCustomProjectionY.SaveAs("MeanChecking/Harmonic_meanProjectionY.png")
