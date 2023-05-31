/*

Mackenzie Scott and Mike Lisa -- May 2021

The purpose of this program is to convert FPGA text files into root objects. 

Compile with g++ 'root-config --cflags' -o ConvertNolanFrames.exe ConvertNolanFrames.cpp 'root-config --libs' -Wno-deprecated -lfftw3 
Execute with path/file on the execution line

*/

#include <iostream>
#include <fstream>

using namespace std;

#include <TGraph.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1D.h>
#include <TProfile2D.h>
#include <TDatime.h>
#include <TTree.h>
#include <TObjString.h>

int DecodeHeader(ifstream* fs);

//All of the variables we need for the header
int NumTimeBuckets,version,year,month,day,hour,minute,sec,T1,T2,corrDelay,scalingBit,nWindows;
float twindowSec,preampGain,AdcGain,samplingNs;
string Date, Time, AVEADC1, AVEADC2, text;
TString* sourceStr;
TDatime* datime;
TTree* headTree;

int main(int argc, char **argv){
  
  //Set Graph Styles
  gStyle->SetOptStat(0); //YES, THIS IS REQUIED EVEN THOUGH WE ARE NOT GRAPHING!!!!!!!!

  //Need Nolan's File
  if (argc < 2){
    cout << "You have to include the name of Nolans file\n";
    return 1;
  }
  
  //Can't Read Nolan's File
  string FileName = argv[1];
  ifstream* NolanDataFile = new ifstream(FileName); 
  if (!NolanDataFile->good()) {cout << "Not good!!\n\n"; return 1;}
  
  //Read in Header info
  cout << "First the header\n";
  if (DecodeHeader(NolanDataFile)!=0){
    cout << "Something wrong with the header - I quit\n";
    return 1;
  }
  
  //Buckets and Frames 
  double bucketMin = -NumTimeBuckets/2 - 0.5;
  double bucketMax =  NumTimeBuckets/2 - 0.5;  // yes, I want to SUBTRACT 0.5 in both cases
  double frameMin  =  0.5; 
  double frameMax  =  nWindows + 0.5; 
  
  //Nanoseconds and Seconds
  double minTime   = bucketMin * samplingNs;
  double maxTime   = bucketMax * samplingNs;
  double startTime = 1.0;
  double endTime   = nWindows * twindowSec + 1.0;
 
  //Define Histograms 
  TH1D* ADC1 = new TH1D("singles1","singles1", nWindows, frameMin, frameMax);  // we just fill once, with average.  but it keeps the average!
  TH1D* ADC2 = new TH1D("singles2","singles2", nWindows, frameMin, frameMax);
  TProfile2D* CF12 = new TProfile2D("CF12","CF12",NumTimeBuckets, minTime, maxTime, nWindows, startTime, endTime, "");
  
  cout << "Now Let's bring in the data\n";
  
  double timeinc = startTime;
  
  // Data Loop
  for(int i = 1; i<=nWindows; ++i){
      
    //Read in Times and ADCs from the run
    *NolanDataFile >> Date;  
    *NolanDataFile >> Time;  
    *NolanDataFile >> AVEADC1;
    *NolanDataFile >> AVEADC2;

    //Fill ADC Histograms 
    ADC1->Fill(i, atof(AVEADC1.c_str()));
    ADC2->Fill(i, atof(AVEADC2.c_str()));
    
    string CF; //Read in values as strings 
    double xval =  minTime;
    
    for (int j = 0; j < NumTimeBuckets; ++j){
      *NolanDataFile >> CF;
      double valNum = atof(CF.c_str());
      CF12->Fill(xval, timeinc, valNum);
      xval += samplingNs;
    }
    timeinc += twindowSec;
  } 
  
  // -------------------------- Make Header TTree ----------------------------------
  // stupid that you have to make a new TTree every time, but otherwise it doesn't get put (fully) into the TFile...
  headTree = new TTree("Header","Header Data");
  headTree->Branch("VersionBr",&version,"Version/I");
  headTree->Branch("SourceBr", "TString", &sourceStr);
  headTree->Branch("T1Br",&T1,"T1/I");
  headTree->Branch("T2Br",&T2,"T2/I");
  headTree->Branch("corrDelayBr",&corrDelay,"CorrDelay/I");
  headTree->Branch("tWindowSecBr",&twindowSec,"tWindowSec/F");
  headTree->Branch("NumWindowsBr",&nWindows,"NumWindows/I");
  headTree->Branch("samplingNsBr",&samplingNs,"SampleNS/F");
  headTree->Branch("scalingBitBr",&scalingBit,"ScalingBit/I");
  headTree->Branch("preampGainBr",&preampGain,"preampGain/F");
  headTree->Branch("AdcGainBr",&AdcGain,"AdcGain/F");
  
  TDatime* dtm = new TDatime(year,month,day,hour,minute,sec);
  dtm->Print();
  
  headTree->Branch("LocalTimeBr","TDatime",&dtm);
  headTree->Fill();
  // -------------------------- Done Header TTree ----------------------------------

  //Save everything to the root file 
  TFile  FrameFile("ZippedFrames.root","RECREATE");
  headTree->Write("Header");
  ADC1->Write("singles1");
  ADC2->Write("singles2");
  CF12->Write("CF12");
  FrameFile.Close();
  
  return 0;
}


// I will do something meaningful with this thing when I get the time.
int DecodeHeader(ifstream* fs){  
  string text, xx;
  for (int i=0; i<15; i++){
    getline(*fs,text);
    cout << i << " : " << text << endl;
    xx = "VERSION,v";        if (text.find(xx)!=string::npos) version        = atoi((text.substr(xx.size(),text.size())).data());
    xx = "SOURCE,";          if (text.find(xx)!=string::npos) sourceStr      = new TString(text.substr(xx.size(),text.size()).data());
    xx = "DATE,";            if (text.find(xx)!=string::npos) year           = atoi((text.substr(5,4)).data()); 
    xx = "DATE,";            if (text.find(xx)!=string::npos) month          = atoi((text.substr(10,2)).data());  
    xx = "DATE,";            if (text.find(xx)!=string::npos) day            = atoi((text.substr(13,2)).data());
    xx = "LOCAL_TIME,";      if (text.find(xx)!=string::npos) hour           = atoi((text.substr(11,2)).data());  
    xx = "LOCAL_TIME,";      if (text.find(xx)!=string::npos) minute         = atoi((text.substr(14,2)).data()); 
    xx = "LOCAL_TIME,";      if (text.find(xx)!=string::npos) sec            = atoi((text.substr(17,2)).data());
    xx = "TEL_X,";           if (text.find(xx)!=string::npos) T1             = atoi((text.substr(xx.size(),text.size())).data());
    xx = "TEL_Y,";           if (text.find(xx)!=string::npos) T2             = atoi((text.substr(xx.size(),text.size())).data());
    xx = "CORR_DELAY,";      if (text.find(xx)!=string::npos) corrDelay      = atoi((text.substr(xx.size(),text.size())).data());
    xx = "N_LAGS,";          if (text.find(xx)!=string::npos) NumTimeBuckets = atoi((text.substr(xx.size(),text.size())).data());
    xx = "TWINDOW_S,";       if (text.find(xx)!=string::npos) twindowSec     = atof((text.substr(xx.size(),text.size())).data());
    xx = "N_WINDOWS,";       if (text.find(xx)!=string::npos) nWindows       = atof((text.substr(xx.size(),text.size())).data());
    xx = "SAMPLING_NS,";     if (text.find(xx)!=string::npos) samplingNs     = atof((text.substr(xx.size(),text.size())).data());
    xx = "SCALING_BIT,";     if (text.find(xx)!=string::npos) scalingBit     = atoi((text.substr(xx.size(),text.size())).data());
    xx = "PREAMP_GAIN_V/A,"; if (text.find(xx)!=string::npos) preampGain     = atof((text.substr(xx.size(),text.size())).data());
    xx = "ADC_GAIN_C/V,";    if (text.find(xx)!=string::npos) AdcGain        = atof((text.substr(xx.size(),text.size())).data());
  }

  return 0;
}

  
