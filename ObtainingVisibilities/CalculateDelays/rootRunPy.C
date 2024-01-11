// josie 18 dec 2023
// run python script to calculate baselines, OPDs, uv coords, WITHOUT running the analysis

void rootRunPy(TString filename){

  // opening input "zipped frames" root file
    TFile* zippedFile = new TFile(filename.Data(),"READONLY");

    gStyle->SetOptStat(0); // do not autofill stat box (allows info from fit later)
    gStyle->SetOptFit(1); // put fit result in stat box

    // get info from the MAIN header first
    TTree* header = new TTree;        zippedFile->GetObject("Header",header);
    TString* source = new TString;    header->SetBranchAddress("SourceBr",&source);         //cout << "star: " << source->Data() << endl;
    int version, nfiles, numpairs;    header->SetBranchAddress("VersionBr",&version);       //header->SetBranchAddress("NfilesBr",&nfiles);   header->SetBranchAddress("NpairsBr",&numpairs)
;
    header->GetEvent(0); // you MUST have this to retrieve correct values from headers

  // =================== get necessary info from headers and run python script to calculate delays ========================================
  TString runpy;
  //TString py = "python3 delays/CalcMvt.py";
  TString py = "python3 delays/CalcMvtUVonly.py";
  TString length = "0";
  TString loctime = "0";
  TString nframes = "0";
  TString frameSize = "0";

  cout << "run taken on " << *(source) << "    version " << version << endl;

  if(version != -4){ cout << "wrong version, bye!" << endl; return;}

  double longestRun(0.0);
  int mostFrames(0);
  double frameWidth(0.0);
  TDatime* ltime = new TDatime;  ltime->Set(2030,1,1,23,59,59); // init to some time far in the future, so that all times in runs would be less
  TDatime* timetemp = new TDatime;

  header->SetBranchAddress("NfilesBr",&nfiles);     header->SetBranchAddress("NpairsBr",&numpairs);
  const int npairs = numpairs;

  // loop over all directories to find longest run, earliest start time, etc to send to python
  TList* keyListPre = gDirectory->GetListOfKeys();
  for (int ip=0; ip<keyListPre->GetSize(); ip++){
    TKey* key = (TKey*)keyListPre->At(ip);
    TString keyName = key->GetName();
    TString className = key->GetClassName();
    if ((keyName.Data()[0]=='T')&&(className=="TDirectoryFile")){ // if we find directory that starts w "T" - take it as a pair
      TTree* pairhead = new TTree;
      zippedFile->GetDirectory(keyName)->GetObject("PairHeader",pairhead);
      TProfile2D* cftemp = new TProfile2D;      zippedFile->GetDirectory(keyName)->GetObject(Form("CF%s", keyName.Data()), cftemp);
      if(cftemp->GetYaxis()->GetXmax() > longestRun){ longestRun = cftemp->GetYaxis()->GetXmax(); }
      if(cftemp->GetYaxis()->GetNbins() > mostFrames){ mostFrames = cftemp->GetYaxis()->GetNbins();  frameWidth = cftemp->GetYaxis()->GetBinWidth(5); }
      pairhead->SetBranchAddress("LocalTimeBr",&timetemp);   pairhead->GetEvent(0);
      if(timetemp->GetTime() < ltime->GetTime()){ ltime = timetemp; }
    }
  }

  // store length of longest run, earliest start time (local) as strings to send to python
  length = to_string(longestRun);   //cout << "to str longest run " << to_string(longestRun) << "  longest run " << longestRun << endl;
  loctime = ltime->AsSQLString();   //cout << "ltime num " << ltime << "  ltime string " << loctime << endl;
  nframes = to_string(mostFrames);  frameSize = to_string(frameWidth);      //cout << "most frames " << nframes << "   width " << frameSize << endl;

  cout << "checking python parameters " << length << "  " << nframes << "  " << frameSize << "  " << *(source) << "  " << loctime << endl << endl;

  // run python script to calculate delays/baselines via bash from within this macro
  cout << "================> Running python script to calculate delays <=====================" << endl;
  runpy = py + " " + length + " " + nframes + " " + frameSize + " " + *(source) + " " + loctime; // note: source is a pointer so must be dereferenced to add to other strings
  gSystem->Exec(runpy);
  cout << "===============================> Python is done! <================================" << endl;

}
