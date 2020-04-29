
#include <array>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <chrono>
#include <list>
#include <pthread.h>
//ROOT libraries, linked by g++ -o targetfile codeFile.C `root-config --cflags --libs`
//
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TString.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TGraph.h>

// Libraries from scratch
#include "globalConstants_snspi.h"
#include "rootSupport.hpp"
// Record length is now globally defined

int main(int argc, char **argv){

  //////////////////////////////
  // Grab the files
  //////////////////////////////
  std::stringstream datum;
  if (argc<2){
    std::cout << "Need to input some file...\n";
    return 1;
  }
  datum.str("");
  datum << argv[1];
  std::ifstream caenTextFile;
  TString fileName= datum.str() + ".root"; // Weird workaround I came up with a while ago to get file names readable by ROOT
  std::cout << "Making " << fileName << '\n';
  TFile *rootFile = TFile::Open(fileName,"RECREATE");
  rootFile->cd();

  //////////////////////////////
  // TTrees for processing
  //////////////////////////////
  TTree *channel0 = new TTree("channel0","Channel 0 data");
  TTree *channel1 = new TTree("channel1","Channel 1 data");
  TTree *channel2 = new TTree("channel2","Channel 2 data");
  TTree *channel3 = new TTree("channel3","Channel 3 data");
  TTree *active_ch = nullptr;

  // Need the * in TTree declaration since we are working with pointers to TTrees (tested in Cling)
  std::list<TTree*> channel_list{channel0,channel1,channel2,channel3};

  caenTextFile.open("output.txt");
  std::string dataline; // Where the output.txt lines are dumped into for processing

  //////////////////////////////
  // Variables for data
  //////////////////////////////
  int channel = 0;
  long timeTag = 0;
  long oldTimeTag = 0, newTimeTag =0;
  long iat = 0;
  int q_short = 0, q_long = 0, baseline = 0, amplitude=0, overshoot=0, n_center=0;
  //int sample[recordLength];
  int waveform[recordLength];
  double xaxis[recordLength];
  double normed_data[recordLength];

  for (int i=0;i<recordLength;i++){
    xaxis[i] = (double)i;
  }

  //////////////////////////////////
  // Define the branches of trees
  //////////////////////////////////
  for(auto &channel:channel_list){
    channel->Branch("q_short",&q_short,"q_short/I");
    channel->Branch("q_long",&q_long,"q_long/I");
    channel->Branch("baseline",&baseline,"baseline/I");
    channel->Branch("amplitude",&amplitude,"amplitude/I");
    channel->Branch("overshoot",&overshoot,"overshoot/I");
    channel->Branch("n_center",&n_center,"n_center/I");
    channel->Branch("iat",&iat,"iat/L");
    channel->Branch("oldTimeTag",&oldTimeTag,"oldTimeTag/L");
    channel->Branch("newTimeTag",&newTimeTag,"newTimeTag/L");
    channel->Branch("waveform",waveform,"waveform[240]/I");
  }


  ///////////////////////////////
  /// Functions for filtering
  ///////////////////////////////

  
  TF1 *noise = new TF1("noise",fNoise,0,240,1); // last number is the number of parameters
  TF1 *one_peak = new TF1("one_peak",fOnePeak,0,240,2);
  TF1 *two_peak = new TF1("two_peak",fTwoPeak,0,240,6);
  std::cout << "Functions made\n";
  //////////////////////////////
  /// Parameters for filtering
  //////////////////////////////

  double gamma =0.;
  double chi_noise=0.;
  double chi_onePeak=0;
  double chi_twoPeak=0.;

  //////////////////////////////
  /// Process the text file
  //////////////////////////////
  auto start = std::chrono::system_clock::now();
  while(std::getline(caenTextFile,dataline)){
    // This next keeps track of what channel we are looking at
    if (dataline.find("# CHANNEL") != std::string::npos){
      channel = dataline[10] - '0';
      if (channel == 0){
        active_ch = channel0;
      }
      if (channel == 1){
        active_ch = channel1;
      }
      if (channel == 2){
        active_ch = channel2;
      }
      if (channel == 3){
        active_ch = channel3;
      }
      timeTag = 0;
      newTimeTag = 0; // reset iat system
      oldTimeTag = 0;
      continue;
    }
    else if (dataline.find("# time_tag") != std::string::npos){
      continue;
    }

    std::stringstream ss(dataline); // converts dataline into a format amenable for processing

    // Before any sort of filtering, load the info from the stream
    ss >> timeTag >> q_short >> q_long >> baseline >> amplitude >> overshoot >> n_center;
    // The positive polarity channes have overshoot for their amplitudes
    if (channel == 0 || channel == 2){
      amplitude = overshoot;
    }
    for (int i=0; i<recordLength; i++){
      ss >> waveform[i];
      normed_data[i] = std::abs((double)waveform[i]-(double)baseline);
      normed_data[i] = normed_data[i]/(double)amplitude;
    }

    TGraph *g = new TGraph(recordLength,xaxis,normed_data); // Generate a TGraph object to call fit method
    //////////////////////////////
    /// Check for one peak
    //////////////////////////////

    // Set parameter limits, the function arguments go as SetParLimits(parameter #, lower_limit, upper_limit)
    one_peak->SetParLimits(0,0,240);
    one_peak->SetParLimits(1,0,1);

    // Initial parameter values in order of by parameter number SetParameters(par0,par1,par2,...)
    one_peak->SetParameters(10,1./150.);
    g->Fit("one_peak","Q"); // Q for quiet, don't need the output
    gamma = g->GetFunction("one_peak")->GetParameter(1); // grabs the decay rate
    chi_onePeak = g->GetFunction("one_peak")->GetChisquare();

    //////////////////////////////////////////////////////////////////////////////
    /// If one peak, then fill the tree and move to next iteration of while loop
    //////////////////////////////////////////////////////////////////////////////
    // The bounds were determined by experimentation
    // Fits of one_peak to noise and double peak events yield bad gamma values outside these bounds

    if ((gamma < 1./50.) && (gamma>1./180.)){
      oldTimeTag = newTimeTag;
      newTimeTag = timeTag;
      iat = newTimeTag - oldTimeTag;
      active_ch -> Fill();
      continue;
    }

    //////////////////////////////////////////////////////////////////////////////
    /// If not one peak, then see whether noise or two peak fits better
    //////////////////////////////////////////////////////////////////////////////
    noise->SetParLimits(0,0,1);
    g->Fit("noise","Q");
    chi_noise = g->GetFunction("noise")->GetChisquare();

    two_peak->SetParLimits(0,0,50);
    two_peak->SetParLimits(1,10,200);
    two_peak->SetParLimits(2,1,20);
    two_peak->SetParLimits(3,1/240.,1/50.);
    two_peak->SetParLimits(4,0,1);
    two_peak->SetParLimits(5,1/240.,1/50.);
    two_peak->SetParameters(10,50,1,1/150.,0.5,1/150.);
    g->Fit("two_peak","MQ"); // M option is for "improved", Minuit throws more function calls
    chi_twoPeak = g->GetFunction("two_peak")->GetChisquare();

    // The above order is nice because it lets one grab the info for the double peak should it be required

    //////////////////////////////////////////////////////////////////////////////
    /// If noise, throw the event away
    //////////////////////////////////////////////////////////////////////////////

    if(chi_noise <= chi_twoPeak || chi_twoPeak > 10.){
      // don't update any time tags and don't fill
      continue;
    }

    //////////////////////////////////////////////////////////////////////////////
    /// If two peaks, then grab data for both events
    //////////////////////////////////////////////////////////////////////////////
    if(chi_noise > chi_twoPeak){
      // log the first peak
      oldTimeTag = newTimeTag;
      newTimeTag = timeTag;
      iat = newTimeTag - oldTimeTag;
      active_ch->Fill();
      // grab info for second peak
      oldTimeTag = newTimeTag;
      newTimeTag = newTimeTag + (int)(g->GetFunction("two_peak")->GetParameter(1));
      iat = (int)(g->GetFunction("two_peak")->GetParameter(1));
      amplitude = (int)(amplitude*g->GetFunction("two_peak")->GetParameter(4));
      std::cout << chi_onePeak << " " << chi_twoPeak << '\n';
      g->Write();
      active_ch->Fill();
    }
  }

  caenTextFile.close();
  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start);
  std::cout << "Filtering Completed in : " << elapsed.count() << '\n';

  ///////////////////////////////////////////////////////////////////////////////////
  // Now we will take the collected data and attempt to extract the real events
  // using positional filtering. Given that the nanowires pixels are 100 um x 100 um
  // and each nanowire is 120 nm wide with 80 nm spacing, there should be roughly
  // 5 cm total length of wire for a pulse to travel at 0.02c. This should mean that
  // the maximum travel time in the wire is ~10 ns.
  ///////////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////
  // Declare the devA and devB trees
  //////////////////////////////////////////

  TTree *devA = new TTree("devA","devA");
  TTree *devB = new TTree("devB","devB");


  ///////////////////////////////////////////////////////////////
  // Find event pairs so that we may now the location of events
  ///////////////////////////////////////////////////////////////
  start = std::chrono::system_clock::now();
  findEventPairs(devA,channel0,channel1);
  findEventPairs(devB,channel2,channel3);

  TTree *positionDevA = new TTree("positionDevA","positionDevA");
  TTree *positionDevB = new TTree("positionDevB","positionDevB");

  iatByPosition(devA,positionDevA);
  iatByPosition(devB,positionDevB);

  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start);

  

  channel0->Write();
  channel1->Write();
  channel2->Write();
  channel3->Write();
  devA->Write();
  devB->Write();
  positionDevA->Write();
  positionDevB->Write();
  rootFile->Close();
  std::cout << "Event Pairing Completed in:" << elapsed.count() << '\n';
  return 0;
}
