#include "ROOT/RDataFrame.hxx"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooChebychev.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TAttText.h"
#include "RooStats/SPlot.h"
#include <iterator>
#include <RooFitResult.h>
#include <RooLinkedList.h>
#include <memory>
#include <tuple>
#include <RooArgList.h>
using namespace RooFit; 
using namespace RooStats;


bool m_isMakePlots = 1;     
   int m_numDEdxBins = 100;                                                 /**< the number of dEdx bins for the payloads */
    int m_numPBins = 69;                                                     /**< the number of momentum bins for the payloads */
    double m_dedxCutoff = 5.e6;                                              /**< the upper edge of the dEdx binning for the payloads */

    /**
     * build the binning scheme
     */
    std::vector<double> CreatePBinningScheme()
    {
      std::vector<double> pbins;
      pbins.reserve(m_numPBins + 1);
      pbins.push_back(0.0);
      pbins.push_back(0.05);

      for (int iBin = 2; iBin <= m_numPBins; iBin++) {
        if (iBin <= 19)
          pbins.push_back(0.025 + 0.025 * iBin);
        else if (iBin <= 59)
          pbins.push_back(pbins.at(19) + 0.05 * (iBin - 19));
        else
          pbins.push_back(pbins.at(59) + 0.3 * (iBin - 59));
      }

      return pbins;
    }
  // std::vector<double> CreatePBinningScheme()
  // {
  //     std::vector<double> pbins;

  //     double p = 0.0;

  //     // Fine resolution from 0.0 to 0.5 (0.0125 steps)
  //     while (p < 0.5) {
  //         pbins.push_back(p);
  //         p += 0.0125;
  //     }

  //     // Medium resolution from 0.5 to 2.5 (0.025 steps)
  //     while (p < 2.5) {
  //         pbins.push_back(p);
  //         p += 0.025;
  //     }

  //     // Finer resolution from 2.5 to 3.0 (0.01 steps)
  //     while (p < 3.0) {
  //         pbins.push_back(p);
  //         p += 0.01;
  //     }

  //     // Final edge
  //     pbins.push_back(3.0);

  //     return pbins;
  // }





void LambdaMassFit()
{

  // gROOT->SetBatch(true); // to not show plots live
  // RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

TFile *_file0 = TFile::Open("Lambda_8Apr.root");
TTree* preselTree = (TTree*) _file0->Get("Lambda");

// define the fit variable
  //RooRealVar InvM("InvM", "m(Lambda)",0.1,1,"GeV/c^{2}");
  RooRealVar InvM("InvM", "m(Lambda)",0.1, 1.11, 1.12, "GeV/c^{2}");

// define other useful variables to be saved
  RooRealVar ProtonMomentum("ProtonMomentum", "momentum for Proton (GeV)", -1.e8, 1.e8);
  RooRealVar ProtonSVDdEdxTrackMomentum("ProtonSVDdEdxTrackMomentum", "momentum for Proton (GeV), from the track", -1.e8, 1.e8);
  RooRealVar ProtonSVDdEdx("ProtonSVDdEdx", "", -1.e8, 1.e8);
  RooRealVar PionLambdaMomentum("PionLambdaMomentum", "momentum for pion (GeV)", -1.e8, 1.e8);
  RooRealVar PionLambdaSVDdEdxTrackMomentum("PionLambdaSVDdEdxTrackMomentum", "momentum for pion (GeV), from the track", -1.e8, 1.e8);
  RooRealVar PionLambdaSVDdEdx("PionLambdaSVDdEdx", "", -1.e8, 1.e8);
  //RooRealVar SlowPionMomentum("SlowPionMomentum", "momentum for slow pion (GeV)", -1.e8, 1.e8);
  //RooRealVar SlowPionSVDdEdxTrackMomentum("SlowPionSVDdEdxTrackMomentum", "momentum for slow pion (GeV), from the track", -1.e8, 1.e8);
  //RooRealVar SlowPionSVDdEdx("SlowPionSVDdEdx", "", -1.e8, 1.e8);


// variables to be saved at the end
  auto variables = new RooArgSet();
  variables->add(InvM);
  variables->add(ProtonMomentum);
  variables->add(ProtonSVDdEdxTrackMomentum);
  variables->add(ProtonSVDdEdx);
  variables->add(PionLambdaMomentum);
  variables->add(PionLambdaSVDdEdxTrackMomentum);
  variables->add(PionLambdaSVDdEdx);
  //variables->add(SlowPionMomentum);
  //variables->add(SlowPionSVDdEdxTrackMomentum);
  //variables->add(SlowPionSVDdEdx);

// dataset to be fitted
  // RooDataSet* LambdaDataset = new RooDataSet("LambdaDataset", "LambdaDataset", preselTree, *variables);
  RooDataSet* LambdaDataset = new RooDataSet("LambdaDataset", "LambdaDataset", preselTree, *variables, "ProtonMomentum > 0.25");

  if (LambdaDataset->sumEntries() == 0) {
    //B2FATAL("The Lambda dataset is empty, stopping here");
  }

  RooPlot* LambdaFitFrame = LambdaDataset->plotOn(InvM.frame());

// define the fit model
// signal model is a bunch of Gaussians 
  RooRealVar GaussMean("GaussMean", "GaussMean", 0.1156,1.114,1.117);//, 0.140, 0.150);
  RooRealVar GaussSigma1("GaussSigma1", "GaussSigma1", 0.0008,0.0003,0.002);// 1.e-4, 1.0);
  RooGaussian LambdaGauss1("LambdaGauss1", "LambdaGauss1", InvM, GaussMean, GaussSigma1);
  RooRealVar GaussSigma2("GaussSigma2", "GaussSigma2", 0.0012,0.0003,0.003);// 1.e-4, 1.0);
  RooGaussian LambdaGauss2("LambdaGauss2", "LambdaGauss2", InvM, GaussMean, GaussSigma2);
  RooRealVar fracGaussYield("fracGaussYield", "Fraction of two Gaussians", 0.5, 0, 1);
  RooAddPdf LambdaSignalPDF("LambdaSignalPDF", "LambdaGauss1+LambdaGauss2", RooArgList(LambdaGauss1, LambdaGauss2), fracGaussYield);

// background model 
  RooRealVar BkgPolyCoef0("BkgPolyCoef0", "BkgPolyCoef0", 0.1, 0., 1.5);
  RooRealVar BkgPolyCoef1("BkgPolyCoef1", "BkgPolyCoef1", -0.5, -1.5, -1.e-3);
  RooChebychev LambdaBkgPDF("LambdaBkgPDF", "BkgPolyPDF", InvM, RooArgList(BkgPolyCoef0, BkgPolyCoef1));

// try this as well later 
  // RooRealVar dm0Bkg("dm0Bkg", "dm0", 0.13957018, 0.130, 0.140);
  // RooRealVar aBkg("aBkg", "a", -0.0784, -0.08, 3.0);
  // RooRealVar bBkg("bBkg", "b", -0.444713, -0.5, 0.4);
  // RooRealVar cBkg("cBkg", "c", 0.3);
  // RooDstD0BG LambdaBkgPDF("LambdaBkgPDF", "LambdaBkgPDF", InvM, dm0Bkg, cBkg, aBkg, bBkg);

// yields of signal and background
  RooRealVar nSignalLambda("nSignalLambda", "signal yield", 0.5 * preselTree->GetEntries(), 0, preselTree->GetEntries());
  RooRealVar nBkgLambda("nBkgLambda", "background yield", 0.5 * preselTree->GetEntries(), 0, preselTree->GetEntries());

  // the total signal + background model
  RooAddPdf totalPDFLambda("totalPDFLambda", "totalPDFLambda pdf", RooArgList(LambdaSignalPDF, LambdaBkgPDF),
                          RooArgList(nSignalLambda, nBkgLambda));

  //B2INFO("Lambda: Start fitting...");
  RooFitResult* LambdaFitResult = totalPDFLambda.fitTo(*LambdaDataset, Save(kTRUE), NumCPU(8));//, PrintLevel(-1));

// optional paranoia checks

  // int status = LambdaFitResult->status();
  // int covqual = LambdaFitResult->covQual();
  // double diff = nSignalLambda.getValV() + nBkgLambda.getValV() - LambdaDataset->sumEntries();

  // //B2INFO("Lambda: Fit status: " << status << "; covariance quality: " << covqual);
  // // if the fit is not healthy, try again once before giving up, with a slightly different setup:
  // if ((status > 0) || (TMath::Abs(diff) > 1.) || (nSignalLambda.getError() < sqrt(nSignalLambda.getValV()))
  //     || (nSignalLambda.getError() > (nSignalLambda.getValV()))) {

  //   LambdaFitResult = totalPDFLambda.fitTo(*LambdaDataset, Save(), Strategy(2), Offset(1));
  //   status = LambdaFitResult->status();
  //   covqual = LambdaFitResult->covQual();
  //   diff = nSignalLambda.getValV() + nBkgLambda.getValV() - LambdaDataset->sumEntries();
  // }

  // if ((status > 0) || (TMath::Abs(diff) > 1.) || (nSignalLambda.getError() < sqrt(nSignalLambda.getValV()))
  //     || (nSignalLambda.getError() > (nSignalLambda.getValV()))) {
  //   //B2WARNING("Lambda: Fit problem: fit status " << status << "; sum of component yields minus the dataset yield is " << diff <<
  //           //   "; signal yield is " << nSignalLambda.getValV() << ", while its uncertainty is " << nSignalLambda.getError());
  // }
  // if (covqual < 2) {
  //   //B2INFO("Lambda: Fit warning: covariance quality " << covqual);
  // }

  totalPDFLambda.plotOn(LambdaFitFrame, LineColor(TColor::GetColor("#4575b4")));

  // double chisquare = LambdaFitFrame->chiSquare();
  //B2INFO("Lambda: Fit chi2 = " << chisquare);
  totalPDFLambda.paramOn(LambdaFitFrame, Layout(0.63, 0.96, 0.93), Format("NEU", AutoPrecision(2)));
  LambdaFitFrame->getAttText()->SetTextSize(0.03);

  totalPDFLambda.plotOn(LambdaFitFrame, Components("LambdaSignalPDF"), LineColor(TColor::GetColor("#d73027")));
  totalPDFLambda.plotOn(LambdaFitFrame, Components("LambdaBkgPDF"), LineColor(TColor::GetColor("#fc8d59")));
  totalPDFLambda.plotOn(LambdaFitFrame, LineColor(TColor::GetColor("#4575b4")));

  LambdaFitFrame->GetXaxis()->SetTitle("#InvM [GeV/c^{2}]");
  TCanvas* canvLambda = new TCanvas("canvLambda", "canvLambda");
  canvLambda->cd();

  LambdaFitFrame->Draw();

  if (m_isMakePlots) {
    canvLambda->Print("SVDdEdxCalibrationFitLambda.pdf");
  }

  /////////////////// SPlot ///////////////////////////////////////////////////////////

  RooStats::SPlot* sPlotDatasetLambda = new RooStats::SPlot("sData", "An SPlot", *LambdaDataset, &totalPDFLambda,
                                                           RooArgList(nSignalLambda, nBkgLambda));

  for (int iEvt = 0; iEvt < 5; iEvt++) {
    if (TMath::Abs(sPlotDatasetLambda->GetSWeight(iEvt, "nSignalLambda") + sPlotDatasetLambda->GetSWeight(iEvt, "nBkgLambda") - 1) > 5.e-3)
      std::cout<<"Lambda: sPlot error: sum of weights not equal to 1"<<std::endl;
  }

  RooDataSet* LambdaDatasetSWeighted = new RooDataSet(LambdaDataset->GetName(), LambdaDataset->GetTitle(), LambdaDataset,
                                                     *LambdaDataset->get());

  TTree* treeLambda_sw = LambdaDatasetSWeighted->GetClonedTree();
  treeLambda_sw->SetName("treeLambda_sw");

  TTree* treePion_sw = LambdaDatasetSWeighted->GetClonedTree();
  treePion_sw->SetName("treePion_sw");
  //B2INFO("Lambda: sPlot done. Proceed to histogramming");

  std::vector<double> pbins = CreatePBinningScheme();

  TCanvas* canvLambda2 = new TCanvas("canvLambda2", "canvLambda2");
  TCanvas* canvLambda3 = new TCanvas("canvLambda3", "canvLambda3");
  // the Proton distribution
  TH2F* hLambdaP = new TH2F("hist_d1_321_trunc", "hist_d1_321_trunc", m_numPBins, pbins.data(),
                           m_numDEdxBins, 0, m_dedxCutoff);

  TH2F* hLambdaPion = new TH2F("hist_d2_321_trunc", "hist_d2_321_trunc", m_numPBins, pbins.data(),
                           m_numDEdxBins, 0, m_dedxCutoff);

  // treeLambda_sw->Draw("ProtonSVDdEdx:ProtonMomentum>>hist_d1_321_trunc", "nSignalLambda_sw * (ProtonSVDdEdx>0)", "");
  treeLambda_sw->Draw("ProtonSVDdEdx:ProtonMomentum>>hist_d1_321_trunc", "nSignalLambda_sw * (ProtonSVDdEdx>0 && ProtonMomentum>0.25)", "");
  treePion_sw->Draw("PionLambdaSVDdEdx:PionLambdaMomentum>>hist_d2_321_trunc", "nSignalLambda_sw * (PionLambdaSVDdEdx>0)", "");
  hLambdaPion->SetAxisRange(0,1,"X");

  canvLambda3->Print("PionDistribution.pdf");

  if (m_isMakePlots) {
      
    canvLambda2->Print("ProtonDistribution.pdf");

  }


  // Profile: Average dEdx vs Momentum (X: Momentum, Y: Avg dEdx)
  TProfile* ProtonProfileX = hLambdaP->ProfileX("ProtonProfileX");

  TProfile* ProtonProfileY = hLambdaP->ProfileY("ProtonProfileY");
  

  // Projection: Total dEdx count vs Momentum
  TH1D* ProtonProjectionX = hLambdaP->ProjectionX("ProtonProjectionX", 30,50);

  // Projection: Total count of entries per dEdx bin
  TH1D* ProtonProjectionY = hLambdaP->ProjectionY("ProtonProjectionY", 30,50);
  // gPad->Modified(); gPad->Update();
  //TH1D* ProtonProfile = (TH1D*)hLambdaP->ProfileX("ProtonProfile"); 
  // gPad->Modified(); gPad->Update();
  // Optionally draw or save them
  TCanvas* c1 = new TCanvas("c1", "Profiles",800,400);
  c1->Divide(2,1);

  ProtonProfileX->SetTitle("ProtonProfileX (p>0.25)");
  c1->cd(1); ProtonProfileX->Draw("E1");
  ProtonProfileX->GetXaxis()->SetRangeUser(0, 1);

  ProtonProfileY->SetTitle("ProtonProfileY (p>0.25)");
  c1->cd(2); ProtonProfileY->Draw("E1");
  ProtonProfileY->GetYaxis()->SetRangeUser(0, 2.5);

  TCanvas* c2 = new TCanvas("c2", "Projections",800,400);
  c2->Divide(2,1);

  ProtonProjectionX->SetTitle("ProtonProjectionX (p>0.25)");
  c2->cd(1); ProtonProjectionX->Draw("E1");
  ProtonProjectionX->GetXaxis()->SetRangeUser(0, 1);

  ProtonProjectionY->SetTitle("ProtonProjectionY (p>0.25)");
  c2->cd(2); ProtonProjectionY->Draw("E1");
  ProtonProjectionY->GetXaxis()->SetRangeUser(0, 5e6);

  //TCanvas* c3 = new TCanvas("c3","Tracks");
  //c3->cd(1); hLambdaP->Draw("colz");

  if (m_isMakePlots) {
      c1->Print("ProtonProfiles.pdf");
      c2->Print("ProtonProjections.pdf");
  }


  // Determine bin index for p > 2.5
  int binAbove2p5 = -1;
  for (size_t i = 0; i < pbins.size(); ++i) {
      if (pbins[i] >= 2.5) {
          binAbove2p5 = i + 1; // ROOT bins start at 1
          break;
      }
  }
  int binMax = hLambdaP->GetNbinsX();

  // ProjectionY: for flat region p > 2.5 GeV
  TH1D* ProtonProjectionY_HighP = hLambdaP->ProjectionY("ProtonProjectionY_HighP", binAbove2p5, binMax);

  TCanvas* c3 = new TCanvas("c3", "HighP Projection", 600, 400);
  ProtonProjectionY_HighP->SetTitle("dE/dx Resolution in Flat Region (p > 2.5 GeV)");
  ProtonProjectionY_HighP->GetXaxis()->SetTitle("dE/dx");
  ProtonProjectionY_HighP->Draw("E1");

  if (m_isMakePlots) {
      c3->Print("ProtonProjectionY_HighP.pdf");
  }


}



