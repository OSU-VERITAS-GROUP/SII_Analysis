// josie 15 aug 2023 (5 oct 2023, 25 oct 2023)
// adapted from Mackenzie's CorrelationAnalysis.C
// NEW version of analysis script for the new versii root file format - processes all 6+ correlations at once, produces ONE analysis.root file per run

// run this macro with:      root VersiiAnalysis.C\(\"ZippedFrames.root\"\)

bool offData = false;
 
//function skeletons
double DoFFT (TProfile2D * CF, int times);
TProfile2D * NoiseRemoveExactFreq(TProfile2D *CF, double FFTFreq, int times);

void VersiiAnalysis(TString filename){
    
  TFile* zippedFile = new TFile(filename.Data(),"READONLY");
  
  gStyle->SetOptStat(0); // do not autofill stat box (allows info from fit later)
  gStyle->SetOptFit(1); // put fit result in stat box
  
  // get info from the MAIN header first
  TTree* header;                    zippedFile->GetObject("Header",header);
  TString* source;                  header->SetBranchAddress("SourceBr",&source);
  int version, nfiles, npairs;      header->SetBranchAddress("VersionBr",&version);   
  
  // ============================= get necessary info from headers and run python script to calculate delays ===================================================
  TString runpy;
  TString py = "python3 CalcDelays.py"
  TString length = "0";
  TString loctime = "0";
  
  if (version > 0){  // FPGA correlation = positive version
      frames = "0";
      width = "0";
      source = "0";
      loctime = "0";
  }
  else if (version < 0){  // OSU correlation = negative version (versii = -3,-4)
      
      double longestRun(0.0);
      TDatime* ltime = new TDatime;  ltime->Set(2030,1,1,1,1,1); // set to some time in the future, so that all times in runs would be less
      TDatime* timetemp;
      
      // get cf for all possible pairs and continue if file exists
      TProfile2D* cfT0T1 = new TProfile2D;      TTree* headT0T1 = new TTree;
      TProfile2D* cfT0T2 = new TProfile2D;      TTree* headT0T2 = new TTree;
      TProfile2D* cfT0T3 = new TProfile2D;      TTree* headT0T3 = new TTree;
      TProfile2D* cfT0T4 = new TProfile2D;      TTree* headT0T4 = new TTree;
      TProfile2D* cfT1T2 = new TProfile2D;      TTree* headT1T2 = new TTree;
      TProfile2D* cfT1T3 = new TProfile2D;      TTree* headT1T3 = new TTree;
      TProfile2D* cfT1T4 = new TProfile2D;      TTree* headT1T4 = new TTree;
      TProfile2D* cfT2T3 = new TProfile2D;      TTree* headT2T3 = new TTree;
      TProfile2D* cfT2T4 = new TProfile2D;      TTree* headT2T4 = new TTree;
      TProfile2D* cfT3T4 = new TProfile2D;      TTree* headT3T4 = new TTree;
      
      zippedFile->GetDirectory("T0T1")->GetObject("CFT0T1",cfT0T1);     zippedFile->GetDirectory("T0T1")->GetObject("PairHeader",headT0T1);
      zippedFile->GetDirectory("T0T2")->GetObject("CFT0T2",cfT0T2);     zippedFile->GetDirectory("T0T2")->GetObject("PairHeader",headT0T2);
      zippedFile->GetDirectory("T0T3")->GetObject("CFT0T3",cfT0T3);     zippedFile->GetDirectory("T0T3")->GetObject("PairHeader",headT0T3);
      zippedFile->GetDirectory("T0T4")->GetObject("CFT0T4",cfT0T4);     zippedFile->GetDirectory("T0T4")->GetObject("PairHeader",headT0T4);
      zippedFile->GetDirectory("T1T2")->GetObject("CFT1T2",cfT1T2);     zippedFile->GetDirectory("T1T2")->GetObject("PairHeader",headT1T2);
      zippedFile->GetDirectory("T1T3")->GetObject("CFT1T3",cfT1T3);     zippedFile->GetDirectory("T1T3")->GetObject("PairHeader",headT1T3);
      zippedFile->GetDirectory("T1T4")->GetObject("CFT1T4",cfT1T4);     zippedFile->GetDirectory("T1T4")->GetObject("PairHeader",headT1T4);
      zippedFile->GetDirectory("T2T3")->GetObject("CFT2T3",cfT2T3);     zippedFile->GetDirectory("T2T3")->GetObject("PairHeader",headT2T3);
      zippedFile->GetDirectory("T2T4")->GetObject("CFT2T4",cfT2T4);     zippedFile->GetDirectory("T2T4")->GetObject("PairHeader",headT2T4);
      zippedFile->GetDirectory("T3T4")->GetObject("CFT3T4",cfT3T4);     zippedFile->GetDirectory("T3T4")->GetObject("PairHeader",headT3T4);
      
      if(cfT0T2){   if(cfT0T2->GetYaxis()->GetXmax() > longestRun){ longestRun = cfT0T2->GetYaxis()->GetXmax(); }
                    headT0T2->SetBranchAddress("LocalTimeBr",&timetemp);     if(timetemp->GetTime() < ltime->GetTime()){ ltime = timetemp;}}
      if(cfT0T3){   if(cfT0T3->GetYaxis()->GetXmax() > longestRun){ longestRun = cfT0T3->GetYaxis()->GetXmax(); }
                    headT0T3->SetBranchAddress("LocalTimeBr",&timetemp);     if(timetemp->GetTime() < ltime->GetTime()){ ltime = timetemp;}}
      if(cfT0T4){   if(cfT0T4->GetYaxis()->GetXmax() > longestRun){ longestRun = cfT0T4->GetYaxis()->GetXmax(); }
                    headT0T4->SetBranchAddress("LocalTimeBr",&timetemp);     if(timetemp->GetTime() < ltime->GetTime()){ ltime = timetemp;}}
      if(cfT1T2){   if(cfT1T2->GetYaxis()->GetXmax() > longestRun){ longestRun = cfT1T2->GetYaxis()->GetXmax(); }
                    headT1T2->SetBranchAddress("LocalTimeBr",&timetemp);     if(timetemp->GetTime() < ltime->GetTime()){ ltime = timetemp;}}
      if(cfT1T3){   if(cfT1T3->GetYaxis()->GetXmax() > longestRun){ longestRun = cfT1T3->GetYaxis()->GetXmax(); }
                    headT1T3->SetBranchAddress("LocalTimeBr",&timetemp);     if(timetemp->GetTime() < ltime->GetTime()){ ltime = timetemp;}}
      if(cfT1T4){   if(cfT1T4->GetYaxis()->GetXmax() > longestRun){ longestRun = cfT1T4->GetYaxis()->GetXmax(); }
                    headT1T4->SetBranchAddress("LocalTimeBr",&timetemp);     if(timetemp->GetTime() < ltime->GetTime()){ ltime = timetemp;}}
      if(cfT2T3){   if(cfT2T3->GetYaxis()->GetXmax() > longestRun){ longestRun = cfT2T3->GetYaxis()->GetXmax(); }
                    headT2T3->SetBranchAddress("LocalTimeBr",&timetemp);     if(timetemp->GetTime() < ltime->GetTime()){ ltime = timetemp;}}
      if(cfT2T4){   if(cfT2T4->GetYaxis()->GetXmax() > longestRun){ longestRun = cfT2T4->GetYaxis()->GetXmax(); }
                    headT2T4->SetBranchAddress("LocalTimeBr",&timetemp);     if(timetemp->GetTime() < ltime->GetTime()){ ltime = timetemp;}}
      if(cfT3T4){   if(cfT3T4->GetYaxis()->GetXmax() > longestRun){ longestRun = cfT3T4->GetYaxis()->GetXmax(); }
                    headT3T4->SetBranchAddress("LocalTimeBr",&timetemp);     if(timetemp->GetTime() < ltime->GetTime()){ ltime = timetemp;}}
      
      length = to_string(longestRun);
      loctime = ltime->AsSQLString();
  }
  
  runpy = py + " " + length + " " + source + " " + loctime;
  gSystem->Exec(runpy);
  cout << "================> Running python script to calculate delays <=====================" << endl;    
    
} // end of main macro
