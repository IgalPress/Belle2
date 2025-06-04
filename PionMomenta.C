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





void PionMomenta()
{

  // gROOT->SetBatch(true); // to not show plots live
  // RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

TFile *_file0 = TFile::Open("Lambda_8Apr.root");
TTree* preselTree = (TTree*) _file0->Get("Lambda");

TFile *_file1 = TFile::Open("Dstar_8Apr.root");
TTree* preselTreeDstar = (TTree*) _file1->Get("Dstar");

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
  RooRealVar SlowPionMomentum("SlowPionMomentum", "momentum for slow pion (GeV)", -1.e8, 1.e8);
  RooRealVar SlowPionSVDdEdxTrackMomentum("SlowPionSVDdEdxTrackMomentum", "momentum for slow pion (GeV), from the track", -1.e8, 1.e8);
  RooRealVar SlowPionSVDdEdx("SlowPionSVDdEdx", "", -1.e8, 1.e8);


// variables to be saved at the end
  auto variables = new RooArgSet();
  auto variables_Dstar = new RooArgSet();
  variables->add(InvM);
  variables->add(ProtonMomentum);
  variables->add(ProtonSVDdEdxTrackMomentum);
  variables->add(ProtonSVDdEdx);
  variables->add(PionLambdaMomentum);
  variables->add(PionLambdaSVDdEdxTrackMomentum);
  variables->add(PionLambdaSVDdEdx);
  variables_Dstar->add(SlowPionMomentum);
  variables_Dstar->add(SlowPionSVDdEdxTrackMomentum);
  variables_Dstar->add(SlowPionSVDdEdx);

// dataset to be fitted
  // RooDataSet* LambdaDataset = new RooDataSet("LambdaDataset", "LambdaDataset", preselTree, *variables);
  RooDataSet* LambdaDataset = new RooDataSet("LambdaDataset", "LambdaDataset", preselTree, *variables, "ProtonMomentum > 0.25");
  RooDataSet* DstarDataset = new RooDataSet("DstarDataset","DstarDataset",preselTreeDstar,*variables_Dstar, "");
  // if (LambdaDataset->sumEntries() == 0) {
  //   //B2FATAL("The Lambda dataset is empty, stopping here");
  // }

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


// yields of signal and background
  RooRealVar nSignalLambda("nSignalLambda", "signal yield", 0.5 * preselTree->GetEntries(), 0, preselTree->GetEntries());
  RooRealVar nBkgLambda("nBkgLambda", "background yield", 0.5 * preselTree->GetEntries(), 0, preselTree->GetEntries());

  // the total signal + background model
  RooAddPdf totalPDFLambda("totalPDFLambda", "totalPDFLambda pdf", RooArgList(LambdaSignalPDF, LambdaBkgPDF),
                          RooArgList(nSignalLambda, nBkgLambda));

  //B2INFO("Lambda: Start fitting...");
  RooFitResult* LambdaFitResult = totalPDFLambda.fitTo(*LambdaDataset, Save(kTRUE), NumCPU(8));//, PrintLevel(-1));


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

  //B2INFO("Lambda: sPlot done. Proceed to histogramming");

  std::vector<double> pbins = CreatePBinningScheme();

  TCanvas* canvLambda3 = new TCanvas("canvLambda3", "canvLambda3");

  // the Proton distribution

  TH2F* hLambdaPion = new TH2F("hist_d2_321_trunc", "hist_d2_321_trunc", m_numPBins, pbins.data(),
                           m_numDEdxBins, 0, m_dedxCutoff);

  treeLambda_sw->Draw("PionLambdaSVDdEdx:PionLambdaMomentum>>hist_d2_321_trunc", "nSignalLambda_sw * (PionLambdaSVDdEdx>0)", "");
  hLambdaPion->SetAxisRange(0,1,"X");

  

}



