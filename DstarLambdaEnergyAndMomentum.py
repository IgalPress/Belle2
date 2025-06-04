import ROOT

latex = ROOT . TLatex ()
latex . SetNDC ()
latex . SetTextSize (0.03)


LambdaFile = ROOT.TFile("Lambda_8Apr.root")
DstarFile = ROOT.TFile("Dstar_8Apr.root")


LambdaTree = LambdaFile.Get("Lambda")
HistogramLambdaPionEnergy = ROOT.TH1D('','',100,0,5e6)

DstarTree = DstarFile.Get("Dstar")
HistogramDstarPionEnergy = ROOT.TH1D('','',100,0,5e6)


HistogramLambdaPionMomentum = ROOT.TH1D('','',100,0,8)

HistogramDstarPionMomentum = ROOT.TH1D('','',100,0,8)

print('Files Loaded')

for i in range(LambdaTree.GetEntries()):
    LambdaTree.GetEntry(i)
    HistogramLambdaPionEnergy.Fill(LambdaTree.PionLambdaSVDdEdx)
    HistogramLambdaPionMomentum.Fill(LambdaTree.PionLambdaMomentum)

for i in range(DstarTree.GetEntries()):
    DstarTree.GetEntry(i)
    HistogramDstarPionEnergy.Fill(DstarTree.PionDSVDdEdx)
    HistogramDstarPionMomentum.Fill(DstarTree.PionDMomentum)

print('Data Loaded')

c = ROOT.TCanvas()
c.cd()
# c.SetLogy(True)
# c.Print("plots.pdf [")

print('Fitting')


gaussFitLambda = ROOT.TF1("gaussfitLambda","gaus",300.e3 ,1000.e3 )
gaussFitDstar = ROOT.TF1("gaussfitDstar","gaus",300.e3 ,1000.e3)


HistogramLambdaPionEnergy.SetLineColor(ROOT.kBlue)
HistogramDstarPionEnergy.SetLineColor(ROOT.kRed)

HistogramLambdaPionEnergy.Draw()
HistogramDstarPionEnergy.Draw("same")

HistogramLambdaPionEnergy.Fit(gaussFitLambda, "R 0")  
HistogramDstarPionEnergy.Fit(gaussFitDstar, "R same 0") 

# gaussFitLambda.SetLineColor(ROOT.kBlue)
# gaussFitDstar.SetLineColor(ROOT.kRed)
# gaussFitLambda.Draw("same")
# gaussFitDstar.Draw("same")

print('Histograms Drawn')

meanLambda = gaussFitLambda . GetParameter (1)/1000
widthLambda = gaussFitLambda . GetParameter (0)/1000

meanDstar = gaussFitDstar.GetParameter(1)/1000
widthDstar = gaussFitDstar.GetParameter(0)/1000

TextHeight = 1
latex . DrawText (0.4 ,TextHeight-0.3 , " MeanDstar = %.1f GeV"%( meanDstar ))
latex . DrawText (0.4 ,TextHeight-0.35 , " WidthDstar = %.1f GeV"%( widthDstar ))

latex . DrawText (0.4 ,TextHeight-0.4 , " MeanLambda = %.1f GeV"%( meanLambda ))
latex . DrawText (0.4 ,TextHeight-0.45 , " WidthLambda = %.1f GeV"%( widthLambda ))

HistogramLambdaPionEnergy.SetLineColor(ROOT.kRed)
HistogramDstarPionEnergy.SetLineColor(ROOT.kBlack)
HistogramLambdaPionEnergy.SetTitle("dE/dx Comparison of pions;dE/dx (ADC);Events")


HistogramLambdaPionEnergy.SetStats(0)
HistogramDstarPionEnergy.SetStats(0)

legend = ROOT.TLegend(0.7 ,0.6 ,0.85 ,0.75)
legend.AddEntry(HistogramLambdaPionEnergy,"Lambda")
legend.AddEntry(HistogramDstarPionEnergy,"Dstar")
# legend.AddEntry(gaussFitDstar,"Fit - Dstar")
# legend.AddEntry(gaussFitLambda,"Fit - Lambda")
legend.SetLineWidth(0)
legend.Draw("same")

print('Plot Customization Completed')

c.Draw()
c.SaveAs("LambdaDstarEnergy.png")
# c.Print("plots.pdf")

print('Starting pion momentum plot')

c2 = ROOT.TCanvas()
c2.cd()
# c2.SetLogy(True)


HistogramDstarPionMomentum.SetLineColor(ROOT.kBlack)
HistogramLambdaPionMomentum.SetLineColor(ROOT.kRed)

HistogramLambdaPionMomentum.Draw("Same")
HistogramDstarPionMomentum.Draw("Same")

# legend.AddEntry(HistogramDstarPionMomentum, "Dstar")
# legend.AddEntry(HistogramLambdaPionMomentum,"Lambda")

gaussFitLambdaPionMomentum = ROOT.TF1("gaussfitLambdaPionMomentum","gaus",0 ,1)
gaussFitDstarPionMomentum = ROOT.TF1("gaussfitDstarPionMomentum","gaus",0.25 ,3.5)

HistogramLambdaPionMomentum.Fit(gaussFitLambdaPionMomentum, "R 0")  
HistogramDstarPionMomentum.Fit(gaussFitDstarPionMomentum, "R same 0") 

meanLambdaPionMomentum = gaussFitLambdaPionMomentum . GetParameter (1)
widthLambdaPionMomentum = gaussFitLambdaPionMomentum . GetParameter (0)

meanDstarPionMomentum = gaussFitDstarPionMomentum.GetParameter(1)
widthDstarPionMomentum = gaussFitDstarPionMomentum.GetParameter(0)

TextHeight = 1
latex . DrawText (0.45 ,TextHeight-0.3 , " MeanDstar = %.1f "%( meanDstarPionMomentum ))
# latex . DrawText (0.45 ,TextHeight-0.35 , " WidthDstar = %.1f "%( widthDstarPionMomentum ))

latex . DrawText (0.45 ,TextHeight-0.35 , " MeanLambda = %.1f "%( meanLambdaPionMomentum ))
# latex . DrawText (0.45 ,TextHeight-0.45 , " WidthLambda = %.1f "%( widthLambdaPionMomentum ))

legend.SetLineWidth(0)
legend.Draw("same")

HistogramLambdaPionMomentum.SetTitle("Momentum Comparison of pions;Momentum;Events")

c2.SaveAs("LambdaDstarMomentum.png")