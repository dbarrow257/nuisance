// Copyright 2016 L. Pickering, P Stowell, R. Terri, C. Wilkinson, C. Wret

/*******************************************************************************
*    This file is part of NUISANCE.
*
*    NUISANCE is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    NUISANCE is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with NUISANCE.  If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************/

#include "MicroBooNE_CC1mu2p0pi_XSec_1D_nu.h"
#include "MicroBooNE_SignalDef.h"
#include "TH2D.h"
#include <cmath>

//********************************************************************
MicroBooNE_CC1mu2p0pi_XSec_1D_nu::MicroBooNE_CC1mu2p0pi_XSec_1D_nu(nuiskey samplekey) {
//********************************************************************
  fSettings = LoadSampleSettings(samplekey);
  std::string name = fSettings.GetS("name");

  if (!name.compare("MicroBooNE_CC1mu2p0pi_XSec_1DDeltaPT_nu")) {
    fDist = kDeltaPT;
    fSettings.SetXTitle("");
    fSettings.SetYTitle("");
  }
  else {
    assert(false);
  }

  // Sample overview ---------------------------------------------------
  std::string descrip = name + " sample.\n" \
                        "Target: Ar\n" \
                        "Flux: BNB numu\n" \
                        "Signal: CC1mu2p0pi\n";

  fSettings.SetDescription(descrip);
  fSettings.SetTitle(name);
  fSettings.SetAllowedTypes("FULL,DIAG/FREE,SHAPE,FIX/SYSTCOV/STATCOV", "FIX/FULL");
  fSettings.SetEnuRange(0.0, 20.0);
  fSettings.DefineAllowedTargets("Ar");
  fSettings.DefineAllowedSpecies("numu");
  FinaliseSampleSettings();

  /*
  std::string DataHistName = "";
  std::string CovMatName = "";
  std::string ACMatName = "";

  // Load data --------------------------------------------------------- 
  std::string inputFile = FitPar::GetDataBase()+DataFileName;
  SetDataFromRootFile(inputFile, DataHistName);
  ScaleData(1E-38);

  fScaleFactor = GetEventHistogram()->Integral("width") / (double(fNEvents) * TotalIntegratedFlux()); // Standard differential cross section per nucleon 
  fScaleFactor *= 1E-38; // Convert units

  SetCovarFromRootFile(inputFile, CovMatName);

  // Load smearing matrix ---------------------------------------------------------
  // Set up the additional smearing matrix Ac
  // All the MC predictions need to be multiplied by Ac to move to the regularized phase space

  TFile* inputRootFile = TFile::Open(inputFile.c_str());
  assert(inputRootFile && inputRootFile->IsOpen());
  TH2D* hSmearMat_TH2 = (TH2D*)inputRootFile->Get(ACMatName.c_str());
  assert(hSmearMat_TH2);

  int nrows = hSmearMat_TH2->GetNbinsX();
  int ncols = hSmearMat_TH2->GetNbinsY();
  fSmearingMatrix = new TMatrixD(nrows, ncols);
  for (int i=0; i<nrows; i++) {
    for (int j=0; j<ncols; j++) {
      (*fSmearingMatrix)(i,j) = hSmearMat_TH2->GetBinContent(i+1, j+1);
    }
  }

  inputRootFile->Close();
  */

  // Final setup ------------------------------------------------------
  FinaliseMeasurement();
};


bool MicroBooNE_CC1mu2p0pi_XSec_1D_nu::isSignal(FitEvent* event) {
  //return SignalDef::MicroBooNE::isCC1ENp(event, EnuMin, EnuMax);
  throw;
};


void MicroBooNE_CC1mu2p0pi_XSec_1D_nu::FillEventVariables(FitEvent* event) {
  fXVar = -999.;

  if (fDist == kDeltaPT) {
    throw;
  }

  throw;
}

void MicroBooNE_CC1mu2p0pi_XSec_1D_nu::ConvertEventRates() {
  // Do standard conversion
  Measurement1D::ConvertEventRates();

  /*
  // Apply Weiner-SVD additional smearing Ac - Needs to be applied on 'event rate' units then converted back to 'xsec units'
  int nBins = fMCHist->GetNbinsX();
  double Flux_CV = 6699174026.68;
  double nTargets = 4.240685683288815e+31;

  // First convert to event rate units
  TVectorD MC_EVRUnits(nBins);
  for (int iBin=0;iBin<nBins;iBin++) {
    MC_EVRUnits(iBin) = fMCHist->GetBinContent(iBin+1) * nTargets * Flux_CV * fMCHist->GetXaxis()->GetBinWidth(iBin+1);
  }

  // Apply smearing
  TVectorD SmearedMC_EVRUnits = (*fSmearingMatrix) * MC_EVRUnits;

  // Convert back to xsec units
  TVectorD SmearedMC_XSecUnits(nBins);
  for (int iBin=0;iBin<nBins;iBin++) {
    SmearedMC_XSecUnits(iBin) = SmearedMC_EVRUnits(iBin) / (nTargets * Flux_CV * fMCHist->GetXaxis()->GetBinWidth(iBin+1));
  }

  // Then copy results back to histogram
  for (int iBin=0;iBin<nBins;iBin++) {
    fMCHist->SetBinContent(iBin+1, SmearedMC_XSecUnits(iBin));
  }
  */
}

