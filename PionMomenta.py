import ROOT
from ROOT import TCanvas

latex = ROOT . TLatex ()
latex . SetNDC ()
latex . SetTextSize (0.03)


LambdaFile = ROOT.TFile("Lambda_8Apr.root")
DstarFile = ROOT.TFile("Dstar_8Apr.root")
####ppasdpfapsdfpasdfpasdfp
LambdaTree = LambdaFile.Get("Lambda")
DstarTree = DstarFile.Get("Dstar")

# HistogramLambdaPionEnergy = ROOT.TH1D('','',100,0,5e6)
# HistogramDstarPionEnergy = ROOT.TH1D('','',100,0,5e6)
# HistogramDstarSlowPionEnergy = ROOT.TH1D('','',100,0,5e6)

# HistogramLambdaPionMomentum = ROOT.TH1D('','',100,0,4)
# HistogramDstarPionMomentum = ROOT.TH1D('','',100,0,4)
# HistogramDstarSlowPionMomentum = ROOT.TH1D('','',100,0,4)

Histogram2D_Kaon = ROOT.TH2D("Histogram2D_Kaon", "Kaon Energy vs Momentum;Energy (arb units);Momentum (GeV/c)",  500, 0,20,500,0, 5e6)
Histogram2D_Proton = ROOT.TH2D("Histogram2D_Proton", "Proton Energy vs Momentum;Energy (arb units);Momentum (GeV/c)",  500, 0,20,500,0, 5e6)
Histogram2D_Lambda = ROOT.TH2D("Histogram2D_Lambda", "Lambda Pion Energy vs Momentum;Energy (arb units);Momentum  (GeV/c)", 500, 0,20,500,0, 5e6)
Histogram2D_DstarPion = ROOT.TH2D("Histogram2D_DstarPion", "D* Pion Energy vs Momentum;Energy (arb units);Momentum (GeV/c)",  500, 0,20,500,0, 5e6)
Histogram2D_DstarSlowPion = ROOT.TH2D("Histogram2D_DstarSlowPion", "Slow Pion Energy vs Momentum;Energy (arb units);Momentum (GeV/c)", 500, 0,20,500,0, 5e6)
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

ProtonMomentum = ROOT.RooRealVar("ProtonMomentum", "momentum for Proton (GeV)", -1e8, 1e8)
ProtonSVDdEdxTrackMomentum = ROOT.RooRealVar("ProtonSVDdEdxTrackMomentum", "momentum for Proton (GeV), from the track", -1e8, 1e8)
ProtonSVDdEdx = ROOT.RooRealVar("ProtonSVDdEdx", "", -1e8, 1e8)
PionLambdaMomentum = ROOT.RooRealVar("PionLambdaMomentum", "momentum for pion (GeV)", -1e8, 1e8)
PionLambdaSVDdEdxTrackMomentum = ROOT.RooRealVar("PionLambdaSVDdEdxTrackMomentum", "momentum for pion (GeV), from the track", -1e8, 1e8)
PionLambdaSVDdEdx = ROOT.RooRealVar("PionLambdaSVDdEdx", "", -1e8, 1e8)

variables_Dstar.add(KaonMomentum)
variables_Dstar.add(PionDMomentum)
variables_Dstar.add(SlowPionMomentum)
variables_Dstar.add(PionDSVDdEdx)
variables_Dstar.add(SlowPionSVDdEdx)
variables_Dstar.add(KaonSVDdEdx)
variables.add(ProtonMomentum)
variables.add(ProtonSVDdEdxTrackMomentum)
variables.add(ProtonSVDdEdx)
variables.add(PionLambdaMomentum)
variables.add(PionLambdaSVDdEdxTrackMomentum)
variables.add(PionLambdaSVDdEdx)

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



c4 = ROOT.TCanvas("c4", "2D Plot", 800, 600)
c4.SetRightMargin(0.15)
c4.SetLeftMargin(0.12)   


treeDstar_sw.Draw("PionDSVDdEdx:PionDMomentum/0.1396 >> Histogram2D_DstarPion", "nSignalDstar_sw* (PionDSVDdEdx>0)", "COLZ")
Histogram2D_DstarPion.SetStats(False)
Histogram2D_DstarPion.SetMinimum(0.)
Histogram2D_DstarPion.SetYTitle("Energy Loss (arb)")
Histogram2D_DstarPion.SetXTitle("P/M (GeV/c)")
c4.Update()
c4.SaveAs("PionD.png")

c8 = ROOT.TCanvas("c8", "ProjectionY of PionD Hist", 800, 600)
c8.SetRightMargin(0.15)
c8.SetLeftMargin(0.12)
ProjectionY_DstarPion = Histogram2D_DstarPion.ProjectionY("ProjectionY_DstarPion")

max_error = 1e3 

for bin in range(1, ProjectionY_DstarPion.GetNbinsX() + 1):
    if ProjectionY_DstarPion.GetBinError(bin) > max_error:
        ProjectionY_DstarPion.SetBinContent(bin, 0)
        ProjectionY_DstarPion.SetBinError(bin, 0)

ProjectionY_DstarPion.SetLineColor(ROOT.kRed)
ProjectionY_DstarPion.SetLineWidth(2)
ProjectionY_DstarPion.SetTitle("ProjectionY: Normalized Event Number vs Mean Energy Loss")
ProjectionY_DstarPion.GetYaxis().SetTitle("Event #")
ProjectionY_DstarPion.GetXaxis().SetTitle("Mean Energy Loss (arb)")
ProjectionY_DstarPion.Draw("E1")
ProjectionY_DstarPion.SetStats(False)

c8.SaveAs("PionDProjectionY.png")

######################################################################################################################################

c9 = ROOT.TCanvas("c9", "2D Plot", 800, 600)
c9.SetRightMargin(0.15)
c9.SetLeftMargin(0.12)   

treeDstar_sw.Draw("SlowPionSVDdEdx:SlowPionMomentum/0.1396 >> Histogram2D_DstarSlowPion", "nSignalDstar_sw*(SlowPionSVDdEdx>0)", "COLZ")
Histogram2D_DstarSlowPion.SetStats(False)
Histogram2D_DstarSlowPion.SetMinimum(0.)
Histogram2D_DstarSlowPion.SetYTitle("Energy Loss (arb)")
Histogram2D_DstarSlowPion.SetXTitle("P/M (GeV/c)")
c9.Update()
c9.SaveAs("SlowPion.png")


######################################################################################################################################

cKaon = ROOT.TCanvas("cKaon", "2D Plot", 800, 600)
cKaon.SetRightMargin(0.15)
cKaon.SetLeftMargin(0.12)   

treeDstar_sw.Draw("KaonSVDdEdx:KaonMomentum/0.4937 >> Histogram2D_Kaon", "nSignalDstar_sw*(KaonSVDdEdx>0)", "COLZ")
Histogram2D_Kaon.SetStats(False)
Histogram2D_Kaon.SetMinimum(0.)
Histogram2D_Kaon.SetYTitle("Energy Loss (arb)")
Histogram2D_Kaon.SetXTitle("P/M (GeV/c)")

cKaon.Update()
cKaon.SaveAs("Kaon.png")

cKaonProfileX = ROOT.TCanvas("cKaonProfileX", "ProfileX of Kaon 2D Histogram", 800, 600)
cKaonProfileX.SetRightMargin(0.15)
cKaonProfileX.SetLeftMargin(0.12)
cKaonProfileX.SetLogy(True)

profileX_Kaon = Histogram2D_Kaon.ProfileX("profileX_Kaon")
profileX_Kaon.SetLineColor(ROOT.kRed)
profileX_Kaon.SetLineWidth(2)
profileX_Kaon.SetTitle("ProfileX: Momentum vs Mean Energy Loss")
profileX_Kaon.GetYaxis().SetTitle("Mean Energy Loss")
profileX_Kaon.GetXaxis().SetTitle("P/M (GeV/C)")

max_error = 1e4 

for bin in range(1, profileX_Kaon.GetNbinsX() + 1):
    if profileX_Kaon.GetBinError(bin) > max_error:
        profileX_Kaon.SetBinContent(bin, 0)
        profileX_Kaon.SetBinError(bin, 0)

profileX_Kaon.Draw("E")
profileX_Kaon.SetStats(False)

FitKaon = ROOT.TF1(
    "FitKaon",
    "[0] + [1]/pow((x + [2]), [3]) + [4]*x",
    0.05, 8
)
FitKaon.SetParameters(100, 5000, 0.01, 12, 10)  # [6] is the slope of the rising tail

profileX_Kaon.Fit(FitKaon, "", "", 0.05, 8)

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


cKaonProfileX.SaveAs("Kaon2DHistogram_ProfileX.png")

#############################################################################################################################
cProton = ROOT.TCanvas("cProton", "2D Plot", 800, 600)
cProton.SetRightMargin(0.15)
cProton.SetLeftMargin(0.12)   

treeLambda_sw.Draw("ProtonSVDdEdx:ProtonMomentum/0.938 >> Histogram2D_Proton", "nSignalLambda_sw*(ProtonSVDdEdx>0)", "COLZ")
Histogram2D_Proton.SetStats(False)
Histogram2D_Proton.SetMinimum(0.)
Histogram2D_Proton.SetYTitle("Energy Loss (arb)")
Histogram2D_Proton.SetXTitle("P/M (GeV/c)")
cProton.Update()
cProton.SaveAs("Proton.png")

cProtonProfileX = ROOT.TCanvas("cProtonProfileX", "ProfileX of Proton 2D Histogram", 800, 600)
cProtonProfileX.SetRightMargin(0.15)
cProtonProfileX.SetLeftMargin(0.12)
cProtonProfileX.SetLogy(True)

profileX_Proton = Histogram2D_Proton.ProfileX("profileX_Proton")
profileX_Proton.SetLineColor(ROOT.kRed)
profileX_Proton.SetLineWidth(2)
profileX_Proton.SetTitle("ProfileX: Momentum vs Mean Energy Loss")
profileX_Proton.GetYaxis().SetTitle("Mean Energy Loss")
profileX_Proton.GetXaxis().SetTitle("P/M (GeV/C)")

max_error = 1e4 

for bin in range(1, profileX_Proton.GetNbinsX() + 1):
    if profileX_Proton.GetBinError(bin) > max_error:
        profileX_Proton.SetBinContent(bin, 0)
        profileX_Proton.SetBinError(bin, 0)

profileX_Proton.Draw("E")
profileX_Proton.SetStats(False)


FitProton = ROOT.TF1(
    "FitProton",
    "[0] + [1]/pow((x + [2]), [3]) + [4]*x",
    0.05, 3
)
FitProton.SetParameters(100, 5000, 0.01, 12, 10)  # [6] is the slope of the rising tail

profileX_Proton.Fit(FitProton, "", "", 0.05, 3)

FitProton.SetLineColor(ROOT.kBlue)
FitProton.SetLineWidth(2)
FitProton.Draw("SAME")

# Get parameter values
p0 = FitProton.GetParameter(0)
p1 = FitProton.GetParameter(1)
p2 = FitProton.GetParameter(2)
p3 = FitProton.GetParameter(3)
p4 = FitProton.GetParameter(4)


# # Format the formula with parameter values
formula_str = (
    f"Fit: {p0:.2e} + {p1:.2e}/(x + {p2:.2e})^{p3:.2e} + {p4:.2e}*x"
)

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextSize(0.03)
latex.DrawLatex(0.15, 0.85, formula_str)


cProtonProfileX.SaveAs("Proton2DHistogram_ProfileX.png")

#############################################################################################################################
c10 = ROOT.TCanvas("c4", "2D Plot", 800, 600)
c10.SetRightMargin(0.15)
c10.SetLeftMargin(0.12)   

treeLambda_sw.Draw("PionLambdaSVDdEdx:PionLambdaMomentum/0.1396 >> Histogram2D_Lambda", "nSignalLambda_sw*(PionLambdaSVDdEdx>0)", "COLZ")
# LambdaTree.Draw("PionLambdaMomentum:PionLambdaSVDdEdx >> h2d(100, 0, 5e6, 100, 0, 4)", "", "COLZ")
Histogram2D_Lambda.SetStats(False)
Histogram2D_Lambda.SetMinimum(0.)
Histogram2D_Lambda.SetYTitle("Energy Loss (arb)")
Histogram2D_Lambda.SetXTitle("P/M (GeV/c)")
c10.Update()
c10.SaveAs("PionLambda.png")

#############################################################################

legend = ROOT.TLegend(0.7 ,0.6 ,0.85 ,0.75)
# legend.AddEntry(HistogramLambdaPionEnergy,"Lambda")
# legend.AddEntry(HistogramDstarPionEnergy,"Dstar")
# legend.AddEntry(gaussFitDstar,"Fit - Dstar")
# legend.AddEntry(gaussFitLambda,"Fit - Lambda")
# legend.SetLineWidth(0)
# legend.Draw("same")

# print('Plot Customization Completed')

# c.Draw()
# c.SaveAs("LambdaDstarEnergy.png")
# # c.Print("plots.pdf")

Histogram2D_Lambda.Add(Histogram2D_DstarPion)
Histogram2D_Lambda.Add(Histogram2D_DstarSlowPion)

# Draw combined histogram
c5 = ROOT.TCanvas("c5", "Combined 2D Histogram", 800, 600)
c5.SetRightMargin(0.15)
Histogram2D_Lambda.SetMinimum(0.)
c5.SetLeftMargin(0.12) 
Histogram2D_Lambda.SetStats(False)
Histogram2D_Lambda.SetTitle("Combined Energy vs Momentum Histogram")
Histogram2D_Lambda.Draw("COLZ")
c5.SaveAs("Combined2DHistogram.png")

##########################################

c6 = ROOT.TCanvas("c6", "ProjectionY of Combined 2D Histogram", 800, 600)
c6.SetRightMargin(0.15)
c6.SetLeftMargin(0.12)
ProjectionY = Histogram2D_Lambda.ProjectionY("ProjectionY",15*25,-1)
ProjectionY.SetLineColor(ROOT.kRed)
ProjectionY.SetLineWidth(2)
ProjectionY.SetTitle("ProjectionY: Event Number vs Mean Energy Loss")
ProjectionY.GetYaxis().SetTitle("Event #")
ProjectionY.GetXaxis().SetTitle("Mean Energy Loss (arb)")
ProjectionY.Draw("E1")
ProjectionY.GetXaxis().SetRangeUser(0, 2500e3)
ProjectionY.SetStats(False)

c6.SaveAs("Combined2DHistogram_ProjectionY.png")

########################################################

c7 = ROOT.TCanvas("c7", "ProfileX of Combined 2D Histogram", 800, 600)
c7.SetRightMargin(0.15)
c7.SetLeftMargin(0.12)
c7.SetLogy(True)
profileX = Histogram2D_Lambda.ProfileX("profileX")
profileX.SetLineColor(ROOT.kRed)
profileX.SetLineWidth(2)
profileX.SetTitle("ProfileX: Momentum vs Mean Energy Loss")
profileX.GetYaxis().SetTitle("Mean Energy Loss")
profileX.GetXaxis().SetTitle("P/M (GeV/C)")


max_error = 1e4 

for bin in range(1, profileX.GetNbinsX() + 1):
    if profileX.GetBinError(bin) > max_error:
        profileX.SetBinContent(bin, 0)
        profileX.SetBinError(bin, 0)

profileX.Draw("E")
profileX.SetStats(False)

betterFit = ROOT.TF1(
    "betterFit",
    "[0] + [1]/pow((x + [2]), [3]) + [4]*x",
    0.05, 20
)
betterFit.SetParameters(100, 5000, 0.01, 12, 10)  # [6] is the slope of the rising tail

profileX.Fit(betterFit, "", "", 0.05, 20)

betterFit.SetLineColor(ROOT.kBlue)
betterFit.SetLineWidth(2)
betterFit.Draw("SAME")

# Get parameter values
p0 = betterFit.GetParameter(0)
p1 = betterFit.GetParameter(1)
p2 = betterFit.GetParameter(2)
p3 = betterFit.GetParameter(3)
p4 = betterFit.GetParameter(4)

formula_str = (
    f"Fit: {p0:.2e} + {p1:.2e}/(x + {p2:.2e})^{p3:.2e} + {p4:.2e}*x"
)

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextSize(0.03)
latex.DrawLatex(0.15, 0.85, formula_str)

c7.SaveAs("Combined2DHistogram_ProfileX.png")


#####################################################################

cProfileFormulae = ROOT.TCanvas("cProfileFormulae", "ProfileX of Pions, Kaons, and Protons", 800, 600)
cProfileFormulae.SetRightMargin(0.15)
cProfileFormulae.SetLeftMargin(0.12)
cProfileFormulae.SetLogy(True)

betterFit.Draw()
FitKaon.Draw('Same')
FitProton.Draw('Same')
betterFit.SetLineColor(ROOT.kRed)
betterFit.SetLineWidth(2)
FitKaon.SetLineColor(ROOT.kBlue)
FitKaon.SetLineWidth(2)
FitProton.SetLineColor(ROOT.kBlack)
FitProton.SetLineWidth(2)


legend = ROOT.TLegend(0.7 ,0.6 ,0.85 ,0.75)
legend.AddEntry(betterFit,"Fit Pions")
legend.AddEntry(FitKaon,"Fit Kaons")
legend.AddEntry(FitProton,"Fit Protons")
legend.SetLineWidth(0)
legend.Draw("Same")

betterFit.GetYaxis().SetRangeUser(5e5, 7e7)
betterFit.SetTitle("Momentum vs Mean Energy Loss")
betterFit.GetYaxis().SetTitle("Mean Energy Loss")
betterFit.GetXaxis().SetTitle("P/M (GeV/C)")

cProfileFormulae.SaveAs("CombinedFormulae.png")
# gaussFitProfileY = ROOT.TF1("gaussfitProfileY","gaus",100e3 ,1000e3)

# profileY.Fit(gaussFitProfileY, "R 0")

# meanProfileY = gaussFitProfileY.GetParameter(1)/1e3
# meanProfileYerr = gaussFitProfileY.GetParError(1)/1e3
# sigmaProfileY = gaussFitProfileY.GetParameter(2)/1e3
# sigmaProfileYerr = gaussFitProfileY.GetParError(2)/1e3
# TextHeight = 1
# latex.DrawText(0.4, TextHeight-0.3, f" MeanProfileY = {meanProfileY:.1f} +/- {meanProfileYerr:.2f} e3")
# latex.DrawText(0.4, TextHeight-0.35, f" SigmaProfileY = {sigmaProfileY:.1f} +/- {sigmaProfileYerr:.2f} e3")

# gaussFitProfileY.SetLineColor(ROOT.kBlue)
# gaussFitProfileY.SetLineWidth(2)
# gaussFitProfileY.Draw("SAME")

# c6.SaveAs("Combined2DHistogram_ProfileY.png")









# print('Starting pion momentum plot')

# c2 = ROOT.TCanvas()
# c2.cd()
# # c2.SetLogy(True)

# # LambdaTree.Draw("HistogramDstarPionMomentum:HistogramDstarPionEnergy>>canvas")

# HistogramDstarPionMomentum.SetLineColor(ROOT.kBlack)
# HistogramLambdaPionMomentum.SetLineColor(ROOT.kRed)
# HistogramDstarSlowPionMomentum.SetLineColor(ROOT.kBlue)

# HistogramLambdaPionMomentum.Draw("Same")
# HistogramDstarPionMomentum.Draw("Same")
# HistogramDstarSlowPionMomentum.Draw("Same")

# # legend.AddEntry(HistogramDstarPionMomentum, "Dstar")
# # legend.AddEntry(HistogramLambdaPionMomentum,"Lambda")

# # gaussFitLambdaPionMomentum = ROOT.TF1("gaussfitLambdaPionMomentum","gaus",0 ,1)
# # gaussFitDstarPionMomentum = ROOT.TF1("gaussfitDstarPionMomentum","gaus",0.25 ,3.5)

# # HistogramLambdaPionMomentum.Fit(gaussFitLambdaPionMomentum, "R 0")  
# # HistogramDstarPionMomentum.Fit(gaussFitDstarPionMomentum, "R same 0") 

# # meanLambdaPionMomentum = gaussFitLambdaPionMomentum . GetParameter (1)
# # widthLambdaPionMomentum = gaussFitLambdaPionMomentum . GetParameter (0)

# # meanDstarPionMomentum = gaussFitDstarPionMomentum.GetParameter(1)
# # widthDstarPionMomentum = gaussFitDstarPionMomentum.GetParameter(0)

# # TextHeight = 1
# # latex . DrawText (0.45 ,TextHeight-0.3 , " MeanDstar = %.1f "%( meanDstarPionMomentum ))
# # # latex . DrawText (0.45 ,TextHeight-0.35 , " WidthDstar = %.1f "%( widthDstarPionMomentum ))

# # latex . DrawText (0.45 ,TextHeight-0.35 , " MeanLambda = %.1f "%( meanLambdaPionMomentum ))
# # # latex . DrawText (0.45 ,TextHeight-0.45 , " WidthLambda = %.1f "%( widthLambdaPionMomentum ))

# legend.AddEntry(HistogramDstarSlowPionMomentum, "D^{*+}->D^{0}\pi^{+}")
# legend.AddEntry(HistogramDstarPionMomentum, "D^{+}->\pi^{+}K^{-}")
# legend.AddEntry(HistogramLambdaPionMomentum, "\Lambda->P\pi^{-}")
# legend.SetLineWidth(0)
# legend.Draw("same")

# HistogramLambdaPionMomentum.SetTitle("Momentum Comparison of pions;Momentum;Events")

# c2.SaveAs("PionMomenta.png")

# Set estimate for the number of entries
treeLambda_sw.SetEstimate(treeLambda_sw.GetEntries() + 1)

# Draw the variable to fill the internal array (no plot, just data)
treeLambda_sw.Draw("ProtonMomentum/0.938", "", "goff")

# Retrieve the array of values
vXP = treeLambda_sw.GetV1()

# Number of bins
m_numPBins = 100

# Calculate bin edges using TKDTreeBinning
kdBinsP = ROOT.TKDTreeBinning(treeLambda_sw.GetEntries(), 1, vXP, m_numPBins)
binsMinEdgesP_orig = kdBinsP.SortOneDimBinEdges()
binsMinEdgesP = ROOT.std.vector('double')(binsMinEdgesP_orig, binsMinEdgesP_orig + m_numPBins + 2)

# Fix the first and last bin edges
binsMinEdgesP[0] = 0.1
binsMinEdgesP[m_numPBins + 1] = 50.

# Print bin edges for checking
for i in range(m_numPBins + 2):
    print(binsMinEdgesP[i])

# Define the histogram with custom binning
m_numDEdxBins = 500  # or your preferred value
m_dedxCutoff = 5e6   # or your preferred value

hLambdaP_bg = ROOT.TH2F(
    "hist_d1_2212_trunc_bg",
    "hist_d1_2212_trunc_bg",
    m_numPBins,
    binsMinEdgesP.data(),
    m_numDEdxBins,
    0,
    m_dedxCutoff
)

# Fill the histogram using the custom binning
treeLambda_sw.Draw(
    "ProtonSVDdEdx:ProtonMomentum/0.938>>hist_d1_2212_trunc_bg",
    "nSignalLambda_sw * (ProtonSVDdEdx>0)",
    "goff"
)

# Now you can draw or save the histogram as usual
cProtonCustom = ROOT.TCanvas("cProtonCustom", "Custom Binning Proton 2D", 800, 600)
hLambdaP_bg.Draw("COLZ")
cProtonCustom.SaveAs("ProtonCustomBinning.png")