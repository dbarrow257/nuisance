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

#include "MicroBooNE_SignalDef.h"
#include "FitUtils.h"

namespace SignalDef {
  namespace MicroBooNE {

bool isCC1MuNp(FitEvent* event, double EnuMin, double EnuMax) {
  // Check CC inclusive
  if (!SignalDef::isCCINC(event, 14, EnuMin, EnuMax)) return false;

  // Veto events which don't have exactly 1 FS muon
  if (event->NumFSMuon() != 1) return false;

  // Veto events with FS mesons
  if (event->NumFSPions() != 0) return false;

  // Veto events with FS electrons
  if (event->NumFSElectron() != 0) return false;

  // Veto events with FS photons
  if (event->NumFSPhoton() != 0) return false;

  // Muon momentum above threshold
  if (event->GetHMFSParticle(13)->fP.Vect().Mag() < 100) return false;

  // Leading proton within momentum range
  if (event->NumFSParticle(2212) == 0) return false;
  double plead = event->GetHMFSParticle(2212)->fP.Vect().Mag();
  if (plead > 300 && plead < 1200) return true;

  return false;
}

bool isCC1ENp(FitEvent* event, double EnuMin, double EnuMax) {
  bool Verbose = false;

  // ==============================================================================================================================
  // Check CC inclusive
  if (!SignalDef::isCCINC(event, 12, EnuMin, EnuMax)) return false;

  // ==============================================================================================================================
  // Veto events which don't have exactly 1 FS electron
  /*
  if (event->NumFSElectron() != 1) {
    if (Verbose) {std::cout << "DB: ElectronCut failed" << std::endl;}
    return false;
  }
  */

  // ==============================================================================================================================
  // Veto events which don't have 1 or more FS protons that have kinetic energy > 40MeV
  std::vector<FitParticle*> ProtonParticles = event->GetAllFSProton();

  double ProtonKEThreshold = 40.0;
  uint nProtonsWithKEAboveThreshold = 0;
  for (uint i=0;i<ProtonParticles.size();i++) {
    if (ProtonParticles[i]->KE()>=ProtonKEThreshold) {
      nProtonsWithKEAboveThreshold += 1;
    }
  }
  if (nProtonsWithKEAboveThreshold == 0) return false;

  // ==============================================================================================================================
  // Veto events with any charged pions that have kinetic energy > 40MeV
  std::vector<FitParticle*> PiPlusParticles = event->GetAllFSPiPlus();
  std::vector<FitParticle*> PiMinusParticles = event->GetAllFSPiMinus();

  std::vector<FitParticle*> ChargedPionParticles;
  for (uint i=0;i<PiPlusParticles.size();i++) {ChargedPionParticles.push_back(PiPlusParticles[i]);}
  for (uint i=0;i<PiMinusParticles.size();i++) {ChargedPionParticles.push_back(PiMinusParticles[i]);}

  double ChargedPionKEThreshold = 40.0;
  uint nChargedPionsWithKEAboveThreshold = 0;
  for (uint i=0;i<ChargedPionParticles.size();i++) {
    if (ChargedPionParticles[i]->KE()>=ChargedPionKEThreshold) {
      nChargedPionsWithKEAboveThreshold += 1;
    }
  }
  if (nChargedPionsWithKEAboveThreshold != 0) return false;

  // ==============================================================================================================================
  // Veto events with any neutral pions
  std::vector<FitParticle*> PiZeroParticles = event->GetAllFSPiZero();
  uint nNeutralPions = PiZeroParticles.size();
  if (nNeutralPions != 0) return false;

  // ==============================================================================================================================
  // Events which have not yet failed the selection are defined as CC1ENP (N>=1)
  return true;
}

bool isCC1mu2p0pi(FitEvent* event, double EnuMin, double EnuMax) {

  // ==============================================================================================================================
  // Check CC inclusive
  if (!SignalDef::isCCINC(event, 14, EnuMin, EnuMax)) return false;

  // ==============================================================================================================================
  // Veto events which don't have exactly one muon with momentum 0.1 < X (GeV) < 1.2

  if (event->NumFSMuon() != 1) return false;
  double MuonMomentum = event->GetHMFSParticle(13)->fP.Vect().Mag();
  if (!(MuonMomentum > 100. && MuonMomentum < 1200.)) return false;

  // ==============================================================================================================================
  // Veto events which don't have exactly 2 FS protons with momentum 0.3 > X(GeV) > 1.0
  std::vector<FitParticle*> ProtonParticles = event->GetAllFSProton();

  double ProtonMomMinThreshold = 300.; //MeV
  double ProtonMomMaxThreshold = 1000.0; //MeV

  uint nProtons = 0;
  for (uint i=0;i<ProtonParticles.size();i++) {
    double ProtonMom = ProtonParticles[i]->fP.Vect().Mag();
    if (ProtonMom > ProtonMomMinThreshold && ProtonMom < ProtonMomMaxThreshold) {
      nProtons += 1;
    }
  }
  if (nProtons != 2) return false;

  // ==============================================================================================================================
  // Veto events with any charged pions with momentum > 65MeV
  std::vector<FitParticle*> PiPlusParticles = event->GetAllFSPiPlus();
  std::vector<FitParticle*> PiMinusParticles = event->GetAllFSPiMinus();

  std::vector<FitParticle*> ChargedPionParticles;
  for (uint i=0;i<PiPlusParticles.size();i++) {ChargedPionParticles.push_back(PiPlusParticles[i]);}
  for (uint i=0;i<PiMinusParticles.size();i++) {ChargedPionParticles.push_back(PiMinusParticles[i]);}

  double ChargedPionMomThreshold = 65.0;
  uint nChargedPionsWithMomAboveThreshold = 0;
  for (uint i=0;i<ChargedPionParticles.size();i++) {
    if (ChargedPionParticles[i]->fP.Vect().Mag() > ChargedPionMomThreshold) {
      nChargedPionsWithMomAboveThreshold += 1;
    }
  }
  if (nChargedPionsWithMomAboveThreshold != 0) return false;

  // ==============================================================================================================================
  // Veto events with any neutral pions
  std::vector<FitParticle*> PiZeroParticles = event->GetAllFSPiZero();
  uint nNeutralPions = PiZeroParticles.size();
  if (nNeutralPions != 0) return false;

  // ==============================================================================================================================
  // Events which have not yet failed the selection are defined as CC1mu2p0pi
  return true;
}

  }  // namespace MicroBooNE
}  // namespace SignalDef

