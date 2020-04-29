#include "TMath.h"
#include "TF1.h"
#include "TTree.h"

#include "rootSupport.hpp"



//////////////////////////////////////////////////////////
/// Model functions
//////////////////////////////////////////////////////////

// Note: When making root functions, the x-value must be passed as a pointer,
// otherwise TF1.h is not a happy camper
Double_t fNoise(Double_t *x, Double_t *par){
  // Parameter Description:
  // p0 - the noise line
  Double_t value = par[0];
  return value;
}

Double_t fOnePeak(Double_t *x, Double_t *par){
  // Parameter Description:
  // p0 - the position of the impact
  // p1 - decay rate
  Double_t value = (x[0]>par[0]?1:0)*TMath::Exp(-par[1]*(x[0]-par[0]));
  return value;
}

Double_t fTwoPeak(Double_t *x, Double_t *par){
  Double_t arg = x[0];
  Double_t value = 0.;
  // Parameter Description:
  // p0 - the position of the first impact
  // p1 - the time between the first and second impact
  // p2 - time between impact and peak of second event
  // p3 - decay rate for first event
  // p4 - amplitude of second event
  // p5 - decay rate for second event
  // the quiet
  if (arg <= par[0]) return value;
  // first decay
  if ((arg > par[0]) && (arg <= (par[0]+par[1]))) {
    value = TMath::Exp(-par[3]*(arg-par[0]));
    return value;
  }
  // second rise
  if ((arg>(par[0]+par[1])) && (arg<=((par[0]+par[1])+par[2]))){
    // slope = dy/dx, so par[2] > 0
    value = (par[4]-TMath::Exp(-par[3]*((par[0]+par[1])-par[0])))/par[2]*(arg-(par[0]+par[1])) + TMath::Exp(-par[3]*((par[0]+par[1])-par[0]));
    return value;
  }
  // second decay
  if (arg>((par[0]+par[1])+par[2])){
    value = par[4]*TMath::Exp(-par[5]*(arg-((par[0]+par[1])+par[2])));
    return value;
  }
  return value;
}



//////////////////////////////////////////////////////////
/// Filtering Procedures
//////////////////////////////////////////////////////////


void findEventPairs(TTree *dev, TTree *pos_ch, TTree *neg_ch){

  Int_t posAmp=0, negAmp=0;
  Long64_t iat = 0;
  Long64_t pos_iat=0, neg_iat=0;
  Long64_t impactPosition=0,decayPosition=0;
  Long64_t relativePosition=0;
  Long64_t posNewTimeTag=0, posOldTimeTag=0;
  Long64_t negNewTimeTag=0, negOldTimeTag=0;

  dev->Branch("iat",&iat,"iat/L");
  dev->Branch("timeTag",&posNewTimeTag,"timeTag/L"); // still record the time of events
  dev->Branch("relativePosition",&relativePosition,"relativePosition/I");
  dev->Branch("impactPosition",&impactPosition,"impactPosition/I");
  dev->Branch("decayPosition",&decayPosition,"decayPosition/I");
  dev->Branch("posAmp",&posAmp,"posAmp/I");
  dev->Branch("negAmp",&negAmp,"negAmp/I");

  pos_ch->SetBranchAddress("iat",&pos_iat);
  pos_ch->SetBranchAddress("newTimeTag",&posNewTimeTag);
  pos_ch->SetBranchAddress("oldTimeTag",&posOldTimeTag);
  pos_ch->SetBranchAddress("amplitude",&posAmp);

  neg_ch->SetBranchAddress("iat",&neg_iat);
  neg_ch->SetBranchAddress("newTimeTag",&negNewTimeTag);
  neg_ch->SetBranchAddress("oldTimeTag",&negOldTimeTag);
  neg_ch->SetBranchAddress("amplitude",&negAmp);

  int nPos = pos_ch->GetEntries();
  int nNeg = neg_ch->GetEntries();

  // For each event in one of the channels, we are looking for an event on the other
  // channel within a physical time window (see comment below)

  for(int i=0;i<nPos;i++){
    pos_ch->GetEntry(i);
    for(int j=0;j<nNeg;j++){
      neg_ch->GetEntry(j);
      impactPosition = posOldTimeTag - negOldTimeTag;
      decayPosition = posNewTimeTag - negNewTimeTag;
      relativePosition = impactPosition - decayPosition;
      // This position cut filter was based on Justins talk on 2019-05-20
      // The events seemed to fall within a +/- 5 ns range
      // Also, looking a the positions for event with a long interArrivalTime
      // the locations are restricted to the +/- 10 ns
      if (std::abs(impactPosition) < 10 && std::abs(decayPosition) < 10){
        iat = (pos_iat + neg_iat)/2; // average the interrArrival times
        dev->Fill();
        continue;
      }
    }
  }
  return;
}


void iatByPosition(TTree *dev, TTree *iatByPosition){

  Int_t posAmp=0, negAmp=0;
  Long64_t iat = 0;
  Int_t impactPosition=0,decayPosition=0;
  Int_t relativePosition=0;
  Long64_t timeTag =0;
  int position;

  Long64_t eventTimeByPosition[17];
  for (int i=0; i<17;i++){
    eventTimeByPosition[i] = 0;
  }

  dev->SetBranchAddress("timeTag",&timeTag); // still record the time of events
  dev->SetBranchAddress("relativePosition",&relativePosition);
  dev->SetBranchAddress("impactPosition",&impactPosition);
  dev->SetBranchAddress("decayPosition",&decayPosition);
  dev->SetBranchAddress("posAmp",&posAmp);
  dev->SetBranchAddress("negAmp",&negAmp);

  iatByPosition->Branch("iatByPosition",&iat,"iatByPosition/L");

  int nEntries = dev->GetEntries();
  for(int i=0;i<nEntries;i++){
    dev->GetEntry(i);
    position = decayPosition + 8;
    if (position < 0 || position > 16){
      continue;
    }
    iat = timeTag - eventTimeByPosition[position];
    eventTimeByPosition[position] = timeTag;
    iatByPosition->Fill();
  }

  return;
}
