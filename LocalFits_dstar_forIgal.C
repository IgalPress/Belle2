#include "ROOT/RDataFrame.hxx"
#include <iterator>
using namespace RooFit; 
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




void DstarMassFit()
{

  // gROOT->SetBatch(true); // to not show plots live
  // RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

TFile *_file0 = TFile::Open("Dstar_8Apr.root");
TTree* preselTree = (TTree*) _file0->Get("Dstar");

// define the fit variable
  RooRealVar deltaM("deltaM", "m(D*)-m(D^{0})", 0.139545, 0.151, "GeV/c^{2}");

// define other useful variables to be saved
  RooRealVar KaonMomentum("KaonMomentum", "momentum for Kaon (GeV)", -1.e8, 1.e8);
  RooRealVar KaonSVDdEdxTrackMomentum("KaonSVDdEdxTrackMomentum", "momentum for Kaon (GeV), from the track", -1.e8, 1.e8);
  RooRealVar KaonSVDdEdx("KaonSVDdEdx", "", -1.e8, 1.e8);
  RooRealVar PionDMomentum("PionDMomentum", "momentum for pion (GeV)", -1.e8, 1.e8);
  RooRealVar PionDSVDdEdxTrackMomentum("PionDSVDdEdxTrackMomentum", "momentum for pion (GeV), from the track", -1.e8, 1.e8);
  RooRealVar PionDSVDdEdx("PionDSVDdEdx", "", -1.e8, 1.e8);
  RooRealVar SlowPionMomentum("SlowPionMomentum", "momentum for slow pion (GeV)", -1.e8, 1.e8);
  RooRealVar SlowPionSVDdEdxTrackMomentum("SlowPionSVDdEdxTrackMomentum", "momentum for slow pion (GeV), from the track", -1.e8, 1.e8);
  RooRealVar SlowPionSVDdEdx("SlowPionSVDdEdx", "", -1.e8, 1.e8);


// variables to be saved at the end
  auto variables = new RooArgSet();
  variables->add(deltaM);
  variables->add(KaonMomentum);
  variables->add(KaonSVDdEdxTrackMomentum);
  variables->add(KaonSVDdEdx);
  variables->add(PionDMomentum);
  variables->add(PionDSVDdEdxTrackMomentum);
  variables->add(PionDSVDdEdx);
  variables->add(SlowPionMomentum);
  variables->add(SlowPionSVDdEdxTrackMomentum);
  variables->add(SlowPionSVDdEdx);

// dataset to be fitted
  RooDataSet* DstarDataset = new RooDataSet("DstarDataset", "DstarDataset", preselTree, *variables);

  if (DstarDataset->sumEntries() == 0) {
    //B2FATAL("The Dstar dataset is empty, stopping here");
  }

  RooPlot* DstarFitFrame = DstarDataset->plotOn(deltaM.frame());

// define the fit model
// signal model is a bunch of Gaussians 
  RooRealVar GaussMean("GaussMean", "GaussMean", 0.145, 0.140, 0.150);
  RooRealVar GaussSigma1("GaussSigma1", "GaussSigma1", 0.01, 1.e-4, 1.0);
  RooGaussian DstarGauss1("DstarGauss1", "DstarGauss1", deltaM, GaussMean, GaussSigma1);
  RooRealVar GaussSigma2("GaussSigma2", "GaussSigma2", 0.001, 1.e-4, 1.0);
  RooGaussian DstarGauss2("DstarGauss2", "DstarGauss2", deltaM, GaussMean, GaussSigma2);
  RooRealVar fracGaussYield("fracGaussYield", "Fraction of two Gaussians", 0.75, 0.0, 1.0);
  RooAddPdf DstarSignalPDF("DstarSignalPDF", "DstarGauss1+DstarGauss2", RooArgList(DstarGauss1, DstarGauss2), fracGaussYield);

// background model 
  RooRealVar BkgPolyCoef0("BkgPolyCoef0", "BkgPolyCoef0", 0.1, 0., 1.5);
  RooRealVar BkgPolyCoef1("BkgPolyCoef1", "BkgPolyCoef1", -0.5, -1.5, -1.e-3);
  RooChebychev DstarBkgPDF("DstarBkgPDF", "BkgPolyPDF", deltaM, RooArgList(BkgPolyCoef0, BkgPolyCoef1));

// try this as well later 
  // RooRealVar dm0Bkg("dm0Bkg", "dm0", 0.13957018, 0.130, 0.140);
  // RooRealVar aBkg("aBkg", "a", -0.0784, -0.08, 3.0);
  // RooRealVar bBkg("bBkg", "b", -0.444713, -0.5, 0.4);
  // RooRealVar cBkg("cBkg", "c", 0.3);
  // RooDstD0BG DstarBkgPDF("DstarBkgPDF", "DstarBkgPDF", deltaM, dm0Bkg, cBkg, aBkg, bBkg);

// yields of signal and background
  RooRealVar nSignalDstar("nSignalDstar", "signal yield", 0.5 * preselTree->GetEntries(), 0, preselTree->GetEntries());
  RooRealVar nBkgDstar("nBkgDstar", "background yield", 0.5 * preselTree->GetEntries(), 0, preselTree->GetEntries());

  // the total signal + background model
  RooAddPdf totalPDFDstar("totalPDFDstar", "totalPDFDstar pdf", RooArgList(DstarSignalPDF, DstarBkgPDF),
                          RooArgList(nSignalDstar, nBkgDstar));

  //B2INFO("Dstar: Start fitting...");
  RooFitResult* DstarFitResult = totalPDFDstar.fitTo(*DstarDataset, Save(kTRUE));//, PrintLevel(-1));

// optional paranoia checks

  // int status = DstarFitResult->status();
  // int covqual = DstarFitResult->covQual();
  // double diff = nSignalDstar.getValV() + nBkgDstar.getValV() - DstarDataset->sumEntries();

  // //B2INFO("Dstar: Fit status: " << status << "; covariance quality: " << covqual);
  // // if the fit is not healthy, try again once before giving up, with a slightly different setup:
  // if ((status > 0) || (TMath::Abs(diff) > 1.) || (nSignalDstar.getError() < sqrt(nSignalDstar.getValV()))
  //     || (nSignalDstar.getError() > (nSignalDstar.getValV()))) {

  //   DstarFitResult = totalPDFDstar.fitTo(*DstarDataset, Save(), Strategy(2), Offset(1));
  //   status = DstarFitResult->status();
  //   covqual = DstarFitResult->covQual();
  //   diff = nSignalDstar.getValV() + nBkgDstar.getValV() - DstarDataset->sumEntries();
  // }

  // if ((status > 0) || (TMath::Abs(diff) > 1.) || (nSignalDstar.getError() < sqrt(nSignalDstar.getValV()))
  //     || (nSignalDstar.getError() > (nSignalDstar.getValV()))) {
  //   //B2WARNING("Dstar: Fit problem: fit status " << status << "; sum of component yields minus the dataset yield is " << diff <<
  //           //   "; signal yield is " << nSignalDstar.getValV() << ", while its uncertainty is " << nSignalDstar.getError());
  // }
  // if (covqual < 2) {
  //   //B2INFO("Dstar: Fit warning: covariance quality " << covqual);
  // }

  totalPDFDstar.plotOn(DstarFitFrame, LineColor(TColor::GetColor("#4575b4")));

  // double chisquare = DstarFitFrame->chiSquare();
  //B2INFO("Dstar: Fit chi2 = " << chisquare);
  totalPDFDstar.paramOn(DstarFitFrame, Layout(0.63, 0.96, 0.93), Format("NEU", AutoPrecision(2)));
  DstarFitFrame->getAttText()->SetTextSize(0.03);

  totalPDFDstar.plotOn(DstarFitFrame, Components("DstarSignalPDF"), LineColor(TColor::GetColor("#d73027")));
  totalPDFDstar.plotOn(DstarFitFrame, Components("DstarBkgPDF"), LineColor(TColor::GetColor("#fc8d59")));
  totalPDFDstar.plotOn(DstarFitFrame, LineColor(TColor::GetColor("#4575b4")));

  DstarFitFrame->GetXaxis()->SetTitle("#Deltam [GeV/c^{2}]");
  TCanvas* canvDstar = new TCanvas("canvDstar", "canvDstar");
  canvDstar->cd();

  DstarFitFrame->Draw();

  if (m_isMakePlots) {
    canvDstar->Print("SVDdEdxCalibrationFitDstar.pdf");
  }

  /////////////////// SPlot ///////////////////////////////////////////////////////////

  RooStats::SPlot* sPlotDatasetDstar = new RooStats::SPlot("sData", "An SPlot", *DstarDataset, &totalPDFDstar,
                                                           RooArgList(nSignalDstar, nBkgDstar));

  for (int iEvt = 0; iEvt < 5; iEvt++) {
    if (TMath::Abs(sPlotDatasetDstar->GetSWeight(iEvt, "nSignalDstar") + sPlotDatasetDstar->GetSWeight(iEvt, "nBkgDstar") - 1) > 5.e-3)
      std::cout<<"Dstar: sPlot error: sum of weights not equal to 1"<<std::endl;
  }

  RooDataSet* DstarDatasetSWeighted = new RooDataSet(DstarDataset->GetName(), DstarDataset->GetTitle(), DstarDataset,
                                                     *DstarDataset->get());

  TTree* treeDstar_sw = DstarDatasetSWeighted->GetClonedTree();
  treeDstar_sw->SetName("treeDstar_sw");

  //B2INFO("Dstar: sPlot done. Proceed to histogramming");

  std::vector<double> pbins = CreatePBinningScheme();

  TCanvas* canvDstar2 = new TCanvas("canvDstar2", "canvDstar");
  // the kaon distribution
  TH2F* hDstarK = new TH2F("hist_d1_321_trunc", "hist_d1_321_trunc", m_numPBins, pbins.data(),
                           m_numDEdxBins, 0, m_dedxCutoff);


  treeDstar_sw->Draw("KaonSVDdEdx:KaonMomentum>>hist_d1_321_trunc", "nSignalDstar_sw * (KaonSVDdEdx>0)", "");

  if (m_isMakePlots) {
      
    canvDstar2->Print("KaonDistribution.pdf");

  }



}



