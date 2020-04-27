#include "TMath.h"
#include "TF1.h"
#include "TTree.h"

#include "rootSupport.hpp"



//////////////////////////////////////////////////////////
/// Model functions
//////////////////////////////////////////////////////////

Double_t fNoise(Double_t x, Double_t *par){
  // Parameter Description:
  // p0 - the noise line
  Double_t value = par[0];
  return value;
}

Double_t fOnePeak(Double_t x, Double_t *par){
  // Parameter Description:
  // p0 - the position of the impact
  // p1 - decay rate
  Double_t value = (x>par[0]?1:0)*TMath::Exp(-par[1]*(x-par[0]));
  return value;
}

Double_t fTwoPeak(Double_t x, Double_t *par){
  Double_t arg = x;
  Double_t value = 0.;
  // Parameter Description:
  // p0 - the position of the first impact
  // p1 - the position of the second impact
  // p2 - time between impact and peak of second event
  // p3 - decay rate for first event
  // p4 - amplitude of second event
  // p5 - decay rate for second event
  // the quiet
  if (arg <= par[0]) return value;
  // first decay
  if ((arg > par[0]) && (arg <= par[1])) {
    value = TMath::Exp(-par[3]*(arg-par[0]));
    return value;
  }
  // second rise
  if ((arg>par[1]) && (arg<=(par[1]+par[2]))){
    // slope = dy/dx, so par[2] > 0
    value = (par[4]-TMath::Exp(-par[3]*(par[1]-par[0])))/par[2]*(arg-par[1]) + TMath::Exp(-par[3]*(par[1]-par[0]));
    return value;
  }
  // second decay
  if (arg>(par[1]+par[2])){
    value = par[4]*TMath::Exp(-par[5]*(arg-(par[1]+par[2])));
    return value;
  }
  return value;
}



//////////////////////////////////////////////////////////
/// Filtering Procedures
//////////////////////////////////////////////////////////


void findEventPairs(TTree *dev, TTree *pos_ch, TTree *neg_ch){

  int posAmp=0, negAmp=0;
  long iat = 0;
  long pos_iat=0, neg_iat=0;
  long impactPosition=0,decayPosition=0;
  long relativePosition=0;
  long posNewTimeTag=0, posOldTimeTag=0;
  long negNewTimeTag=0, negOldTimeTag=0;


  dev->Branch("iat",&iat,"iat/L");
  dev->Branch("relativePosition",&relativePosition,"relativePosition/I");
  dev->Branch("impactPosition",&impactPosition,"impactPosition/I");
  dev->Branch("decayPosition",&decayPosition,"decayPosition/I");
  dev->Branch("posAmp",&posAmp,"posAmp/I");
  dev->Branch("negAmp",&negAmp,"negAmp/I");

  pos_ch->SetBranchAddress("iat",&pos_iat);
  pos_ch->SetBranchAddress("newTimeTag",&posNewTimeTag);
  pos_ch->SetBranchAddress("oldTimeTag",&posOldTimeTag);
  pos_ch->SetBranchAddress("posAmp",&posAmp);

  neg_ch->SetBranchAddress("iat",&neg_iat);
  neg_ch->SetBranchAddress("newTimeTag",&negNewTimeTag);
  neg_ch->SetBranchAddress("oldTimeTag",&negOldTimeTag);
  neg_ch->SetBranchAddress("negAmp",&negAmp);

  int nPos = pos_ch->GetEntries();
  int nNeg = neg_ch->GetEntries();
  for(int i=0;i<nPos;i++){
    pos_ch->GetEntry(i);
    for(int j=0;j<nNeg;j++){
      neg_ch->GetEntry(j);
      impactPosition = posOldTimeTag - negOldTimeTag;
      decayPosition = posNewTimeTag - negNewTimeTag;
      relativePosition = impactPosition - decayPosition;
      // This position cut filter was based on Justins talk on 2019-05-20
      // The events seemed to fall within a +/- 5 cm range
      if (std::abs(impactPosition) < 20 && std::abs(decayPosition) < 20){
        iat = (pos_iat + neg_iat)/2; // average the interrArrival times
        dev->Fill();
      }
    }
  }
  return;
}
