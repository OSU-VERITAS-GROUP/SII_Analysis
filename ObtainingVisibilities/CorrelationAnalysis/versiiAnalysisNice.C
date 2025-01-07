// josie 25 june 2024
// full analysis script - copied from versiiAnalysisFixOffRuns.C - 25 june 2024, add in john's frame by frame noise removal method and identifying noise in power spectra
// run with: root versiiAnalysisFixOffRuns.C\(\"eta\ Tau_y2023m12d28h23m54s52_ZippedFrames.root\",1\)
// 18 june 2024 - josie fixed saving the header to the analysis file
// 23 july 2024 - josie adding removal of frames with bad telescope tracking or dip in ADCs
// 6 nov 2024 - josie fixing normalization (after noise removal and dividing, not subtracting) and making functions etc nicer, removing all version -3 artifacts

// function headers
TNtuple* GeometryNtuple(TString keyName, int N);
double DoFFT (TProfile2D * CF, int times);
TProfile2D * NoiseRemoveExactFreq(TProfile2D *CF, double FFTFreq, int times);
TProfile2D* subtractOne(TProfile2D* cf);
//TProfile2D* forcedNormSubtOne(TProfile2D* cf, int one = 1);
TProfile2D* forcedNormSubtOne(TProfile2D* cf, int one = 1, TH1D* key = 0);
TH1D* identifyNoise(TProfile2D* powerSpectrum1, TProfile2D* powerSpectrum2, int noiseBin); // function to return a hist with noisy frames (filled w 1 or 0 for noise/no noise)
TProfile2D* removeFrames(TProfile2D* cf, TH1D* noisekey); 

void versiiAnalysisNice(TString filename, int off = 1, int rtf = 0){ // 1=incl off data, 0=no off data, int rtf for remove tracking frames, 0 = no, 1 = yes

    // define final frame size after rebinning
    int rebin(1);
    
    // opening input "zipped frames" root file
    TFile* zippedFile = new TFile(filename.Data(),"READONLY");
    
    gStyle->SetOptStat(0); // do not autofill stat box (allows info from fit later)
    gStyle->SetOptFit(1); // put fit result in stat box
    
    // get info from the MAIN header first
    TTree* header = new TTree;              zippedFile->GetObject("Header",header);
    TString* source = new TString;          header->SetBranchAddress("SourceBr",&source);         //cout << "star: " << source->Data() << endl;
    int version, nfiles, numpairs;          header->SetBranchAddress("VersionBr",&version);       header->SetBranchAddress("NfilesBr",&nfiles);   header->SetBranchAddress("NpairsBr",&numpairs);
    TString* filestr = new TString;         header->SetBranchAddress("FileListBr",&filestr);
    header->GetEvent(0); // you MUST have this to retrieve correct values from headers
  
  // =================== get necessary info from headers and run python script to calculate delays ========================================
  TString runpy;
  TString py = "python3 /users/PAS1977/jrose/macros/calcDelays/CalcMvt_hourangle.py"; //CalcMvtUVonly.py";
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
  
  header->SetBranchAddress("NfilesBr",&nfiles);     header->SetBranchAddress("NpairsBr",&numpairs); // remove this??
  const int npairs = numpairs;
 
  // loop over all directories to find longest run, earliest start time, etc to send to python
  TList* keyListPre = gDirectory->GetListOfKeys();
  for (int ip=0; ip<keyListPre->GetSize(); ip++){
    TKey* key = (TKey*)keyListPre->At(ip);
    TString keyName = key->GetName();
    TString className = key->GetClassName();
    if ((keyName.Data()[0]=='T')&&(className=="TDirectoryFile")){ // if we find directory that start w "T" - take it as a pair
      TTree* pairhead = new TTree;
      zippedFile->GetDirectory(keyName)->GetObject("PairHeader",pairhead);
      TProfile2D* cftemp = new TProfile2D;      
      zippedFile->GetDirectory(keyName)->GetObject(Form("CF%s", keyName.Data()), cftemp);
      if(cftemp->GetYaxis()->GetXmax() > longestRun){ longestRun = cftemp->GetYaxis()->GetXmax(); }
      if(cftemp->GetYaxis()->GetNbins() > mostFrames){ mostFrames = cftemp->GetYaxis()->GetNbins();  frameWidth = cftemp->GetYaxis()->GetBinWidth(5); }
      pairhead->SetBranchAddress("LocalTimeBr",&timetemp);   pairhead->GetEvent(0);
      if(timetemp->GetTime() < ltime->GetTime()){ ltime = timetemp; }
    }
  }
  // divide by rebinning factor to equal final number of frames at OPD shift
  mostFrames = mostFrames/rebin; //16;
  frameWidth = frameWidth*double(rebin); //16.0;
  
  // for the runs we have a 10 MHz file, we will use that timestamp instead for the start time
  
  cout << "tdatime first time is " << ltime->AsSQLString() << endl;
  int pyyeari(0), pymonthi(0), pydayi(0), pyhouri(0), pymini(0), pyseci(0);
  string pyyear, pymonth, pyday, pyhour, pymin, pysec;
  double tugtimedif(0.0);
  pyyeari = ltime->GetYear();         pyyear = to_string(pyyeari);        
  pymonthi = ltime->GetMonth();       pymonth = to_string(pymonthi);      if(pymonthi < 10){ pymonth = "0" + pymonth; }
  pydayi = ltime->GetDay();           pyday = to_string(pydayi);          if(pydayi < 10){ pyday = "0" + pyday; }
  pyhouri = ltime->GetHour();         pyhour = to_string(pyhouri);        if(pyhouri < 10){ pyhour = "0" + pyhour; }
  pymini = ltime->GetMinute();        pymin = to_string(pymini);          if(pymini < 10){ pymin = "0" + pymin; }
  pyseci = ltime->GetSecond();        pysec = to_string(pyseci);          if(pyseci < 10){ pysec = "0" + pysec; }
  
  cout << "pydate: " << pyyear << " " << pymonth << " " << pyday << " " << pyhour << " " << pymin << " " << pysec << endl;
   
  string line2, name10mhz;
  string files10mhz[28];
  double shift12(0.0), shift13(0.0), shift14(0.0), shift23(0.0), shift24(0.0), shift34(0.0);
  double shifts12[28], shifts13[28], shifts14[28], shifts23[28], shifts24[28], shifts34[28];
  int itug(0);
  
  ifstream tug10file;     tug10file.open("10MHzShifts.txt");    // this is a good candidate to turn into a function - return loctime
  if(tug10file){ 
      while(getline(tug10file, line2)){
          tug10file >> name10mhz >> shift12 >> shift13 >> shift14 >> shift23 >> shift24 >> shift34;
          files10mhz[itug] = name10mhz;
          shifts12[itug] = shift12;
          shifts13[itug] = shift13;
          shifts14[itug] = shift14;
          shifts23[itug] = shift23;
          shifts24[itug] = shift24;
          shifts34[itug] = shift34;
          itug++;
          
          //cout << name10mhz << endl;
      }
      tug10file.close();
      
      string tugyear, tugmonth, tugday, tughour, tugmin, tugsec;
      
      if(pyyear == "2024" && pymonth == "02"){
          // split 10 mhz file to get date info (and calc differences)
          //"T1_ON_y2024m02d18h19m32s48_10MHzSyncPulse.txt"
          for (int itf=0; itf<28; itf++){
              string tugfname;    tugfname = files10mhz[itf];
              tugyear = tugfname.substr(8,4);
              tugmonth = tugfname.substr(13,2);
              tugday = tugfname.substr(16,2);
              tughour = tugfname.substr(19,2);
              tugmin = tugfname.substr(22,2);
              tugsec = tugfname.substr(25,2);
              
              if(pyday == tugday && pyhour == tughour){ // do not need to compare year, month here bc it's all feb 2024 - but in future maybe for the peak shift
              // there is only one case where we have 2 runs in the same hour for gam cas - so we will treat it uniquely bc this doesn't need to be so generalized
                  if(pyday == "21" && pyhour == "20"){
                      if(stoi(pymin) < 40 && stoi(tugmin) < 40){ break; }
                      else if(stoi(pymin) > 40 && stoi(tugmin) > 40){ break; }
                  }
                  else{ break; }
              }
          }
          cout << tugyear << " " << tugmonth << " " << tugday << " " << tughour << " " << tugmin << " " << tugsec << endl;
          //2024-02-24 04:13:29
          loctime = tugyear + "-" + tugmonth + "-" + tugday + " " + tughour + ":" + tugmin + ":" + tugsec;
      }
      else{ loctime = ltime->AsSQLString(); }
  }
  else{ loctime = ltime->AsSQLString(); }
    
  // store length of longest run, earliest start time (local) as strings to send to python    
  length = to_string(longestRun);   //cout << "to str longest run " << to_string(longestRun) << "  longest run " << longestRun << endl;
  //loctime = ltime->AsSQLString();   //cout << "ltime num " << ltime << "  ltime string " << loctime << endl; 
  nframes = to_string(mostFrames);  frameSize = to_string(frameWidth);      //cout << "most frames " << nframes << "   width " << frameSize << endl; 
  
  cout << "checking python parameters " << length << "  " << nframes << "  " << frameSize << "  " << *(source) << "  " << loctime << endl << endl; 
  
  // run python script to calculate delays/baselines via bash from within this macro 
  cout << "================> Running python script to calculate delays <=====================" << endl; 
  runpy = py + " " + length + " " + nframes + " " + frameSize + " " + *(source) + " " + loctime; // note: source is a pointer so must be dereferenced to add to other strings
  gSystem->Exec(runpy);
  cout << "===============================> Python is done! <================================" << endl;

  TCanvas* c1 = new TCanvas;

  // make root file to save objs
  TString outName = filename.ReplaceAll("ZippedFrames","analysis");
  TFile* outfile = new TFile(Form("%s",outName.Data()),"RECREATE");
  cout << "saving everything to root file named: " << outName.Data() << endl;
  
  header->CloneTree()->Write("header");
  
  // ===================================================== start of the main analysis ===========================================================
  
  // loop over each directory we find in the root file (taken from AnalyzeVersiiFiles.C)
  zippedFile->cd();
  TList* keyList = gDirectory->GetListOfKeys();
  TString pairsList[npairs];
  
  for (int ip=0; ip<keyList->GetSize(); ip++){
      TKey* key = (TKey*)keyList->At(ip);
      TString keyName = key->GetName();
      TString className = key->GetClassName();
      if ((keyName.Data()[0]=='T')&&(className=="TDirectoryFile")){ // if we find directory that start w "T" - take it as a pair
          cout << "now starting pair " << keyName.Data() << endl;
          TString tee1 = keyName[1];    int t1 = atoi(tee1.Data()); 
          TString tee2 = keyName[3];    int t2 = atoi(tee2.Data());
          pairsList[ip] = keyName;
          
          // get pair header tree, for now all we need is the correlation delay (num buckets offset) but can print the rest as a check
          TTree* pairhead = new TTree;
          zippedFile->GetDirectory(keyName)->GetObject("PairHeader",pairhead);
          int nbo;   pairhead->SetBranchAddress("corrDelayBr",&nbo);   pairhead->GetEvent(0);
          
          // get original correlation function and draw
          TProfile2D* cfOrig = new TProfile2D;      
          zippedFile->GetDirectory(keyName)->GetObject(Form("CF%s", keyName.Data()), cfOrig);
          cfOrig->SetTitle(Form("raw correlation function %s;relative time (ns); frames (s)", keyName.Data()));
          outfile->cd();
          gDirectory->mkdir(Form("%s",keyName.Data()));
          outfile->cd(Form("%s",keyName.Data()));
          
          pairhead->CloneTree()->Write("pairheader");
          
          // get info from python-generated text files and create TNtuple -- note: super simple to draw baselines, etc w TNtuple --> PairGeometryInfo->Draw("v:u")
          TNtuple* geom = GeometryNtuple(keyName, int(mostFrames));
          geom->Write();
          
          cfOrig->Write("originalCF");
          
          // ------------------------------------------------------------- off data ---------------------------------------------------------------
          
          // off ADCs before and after run for both telescopes, and increments for change in ADC 
          double beforeADC1(0.0), afterADC1(0.0), beforeADC2(0.0), afterADC2(0.0); 
          double beforeADC1temp(0.0), afterADC1temp(0.0), beforeADC2temp(0.0), afterADC2temp(0.0); // fix issue with time matching correct but wrong ADC value - variable modified later by getevent
          double dOffadc1, dOffadc2;
          
          if (off == 1){
              
              // get start and end times of run and convert to UTC time
              TDatime* timeStart = new TDatime;    TDatime* timeEnd = new TDatime;
              //int timeStartConv, timeEndConv;
              int runlength = int(cfOrig->GetYaxis()->GetXmax());
              int addMin(0), addHr(0), addDay(0);
              double setSec(0.0), setMin(0.0), setHr(0.0);
              pairhead->SetBranchAddress("LocalTimeBr",&timeStart);    pairhead->GetEvent(0);
              if((timeStart->GetSecond() + (runlength%3600%60)) > 60.0){ addMin = 1; setSec = (timeStart->GetSecond() + (runlength%3600%60)) - 60.0; } // adding this to fix the weird time overflow problems
              else{ setSec = timeStart->GetSecond() + (runlength%3600%60); }
              if((timeStart->GetMinute() + (runlength%3600/60) + addMin) > 60.0){ addHr = 1; setMin = (timeStart->GetMinute() + (runlength%3600/60) + addMin) - 60.0; } 
              else{ setMin = timeStart->GetMinute() + (runlength%3600/60) + addMin; }
              if((timeStart->GetHour() + (runlength/3600) + addHr) > 24.0){ addDay = 1; setHr = (timeStart->GetHour() + (runlength/3600) + addHr) - 24.0; }
              else{ setHr = timeStart->GetHour() + (runlength/3600) + addHr; }
              timeEnd->Set(timeStart->GetYear(), timeStart->GetMonth(), timeStart->GetDay() + addDay, setHr, setMin, setSec);
              //timeEndOLD->Set(timeStart->GetYear(), timeStart->GetMonth(), timeStart->GetDay(), timeStart->GetHour() + (runlength/3600), timeStart->GetMinute() + (runlength%3600/60), timeStart->GetSecond() + (runlength%3600%60));
              cout << "start time: " << timeStart->AsSQLString() << "    end time: " << timeEnd->AsSQLString() << "   length: " << runlength << endl;
              //timeStartConv = timeStart->Convert(); 
              //timeEndConv = timeEnd->Convert(); //timeStartConv + runlength;
              
              // variables to hold off data info while looping to match correct off runs
              TTree* offbranch1 = new TTree;        zippedFile->GetObject(Form("OFF/OffDataT%s", tee1.Data()), offbranch1);  // get branches for both telescopes in pair
              TTree* offbranch2 = new TTree;        zippedFile->GetObject(Form("OFF/OffDataT%s", tee2.Data()), offbranch2);
              TDatime* tempofftime1 = new TDatime;      TDatime* tempofftime2 = new TDatime; // temps to hold each off run in loop, before and after to hold times closest to start/end
            
              // initialize before and after times to long/before after run - must initialize with a value so comparison works
              TDatime* beforeRun1 = new TDatime;    beforeRun1->Set(timeStart->GetYear(), timeStart->GetMonth(), timeStart->GetDay(), timeStart->GetHour() - 1.0, timeStart->GetMinute(), timeStart->GetSecond());
              TDatime* afterRun1 = new TDatime;     afterRun1->Set(timeStart->GetYear(), timeStart->GetMonth(), timeStart->GetDay(), timeStart->GetHour() + 5.0, timeStart->GetMinute(), timeStart->GetSecond());
              TDatime* beforeRun2 = new TDatime;    beforeRun2->Set(timeStart->GetYear(), timeStart->GetMonth(), timeStart->GetDay(), timeStart->GetHour() - 1.0, timeStart->GetMinute(), timeStart->GetSecond());
              TDatime* afterRun2 = new TDatime;     afterRun2->Set(timeStart->GetYear(), timeStart->GetMonth(), timeStart->GetDay(), timeStart->GetHour() + 5.0, timeStart->GetMinute(), timeStart->GetSecond());
              
              // loop through off run times to find closest off runs before start and after end of run - must check both telescopes bc they may not have same number of off runs
              for (int io=0; io<offbranch1->GetEntries(); io++){
                  offbranch1->SetBranchAddress("LocalTimBr",&tempofftime1);  offbranch1->GetEvent(io); // yes it should be LocalTim, typo in root file
                  if((tempofftime1->Convert() < timeStart->Convert()) && (timeStart->Convert() - tempofftime1->Convert() < timeStart->Convert() - beforeRun1->Convert())){ 
                      beforeRun1->Set(tempofftime1->GetYear(), tempofftime1->GetMonth(), tempofftime1->GetDay(), tempofftime1->GetHour(), tempofftime1->GetMinute(), tempofftime1->GetSecond());
                      offbranch1->SetBranchAddress("AveAdcBr",&beforeADC1temp);     offbranch1->GetEvent(io);   beforeADC1 = beforeADC1temp;
                      //cout << "ck   beforeADC1: " << beforeADC1 << "   event " << io << endl;
                  }
                  if((tempofftime1->Convert() > timeEnd->Convert()) && (tempofftime1->Convert() - timeEnd->Convert() < afterRun1->Convert() - timeEnd->Convert())){ 
                      afterRun1->Set(tempofftime1->GetYear(), tempofftime1->GetMonth(), tempofftime1->GetDay(), tempofftime1->GetHour(), tempofftime1->GetMinute(), tempofftime1->GetSecond());
                      offbranch1->SetBranchAddress("AveAdcBr",&afterADC1temp);      offbranch1->GetEvent(io);   afterADC1 = afterADC1temp;
                      //cout << "ck   afterADC1: " << afterADC1 << "   event " << io << endl;
                  }
              }
              // and repeat for the second telescope - just in case one would cut out at the beginning/end of the night
              for (int io=0; io<offbranch2->GetEntries(); io++){
                  offbranch2->SetBranchAddress("LocalTimBr",&tempofftime2);  offbranch2->GetEvent(io); // yes it should be LocalTim, typo in root file
                  if((tempofftime2->Convert() < timeStart->Convert()) && (timeStart->Convert() - tempofftime2->Convert() < timeStart->Convert() - beforeRun2->Convert())){ 
                      beforeRun2->Set(tempofftime2->GetYear(), tempofftime2->GetMonth(), tempofftime2->GetDay(), tempofftime2->GetHour(), tempofftime2->GetMinute(), tempofftime2->GetSecond());
                      offbranch2->SetBranchAddress("AveAdcBr",&beforeADC2temp);     offbranch2->GetEvent(io);   beforeADC2 = beforeADC2temp;
                      //cout << "ck   beforeADC2: " << beforeADC2 << "   event " << io << endl;
                  }
                  if((tempofftime2->Convert() > timeEnd->Convert()) && (tempofftime2->Convert() - timeEnd->Convert() < afterRun2->Convert() - timeEnd->Convert())){ 
                      afterRun2->Set(tempofftime2->GetYear(), tempofftime2->GetMonth(), tempofftime2->GetDay(), tempofftime2->GetHour(), tempofftime2->GetMinute(), tempofftime2->GetSecond());
                      offbranch2->SetBranchAddress("AveAdcBr",&afterADC2temp);      offbranch2->GetEvent(io);   afterADC2 = afterADC2temp;
                      //cout << "ck   afterADC2: " << afterADC2 << "   event " << io << endl;
                  }
              }
              // checking in case of errors
              if(abs(beforeADC1) > 50){ cout << "ERROR!!! before ADC for " << tee1.Data() << " = " << beforeADC1 << endl;   beforeADC1 = 0.0; }   
              if(abs(afterADC1) > 50){ cout << "ERROR!!! after ADC for " << tee1.Data() << " = " << afterADC1 << endl;      afterADC1 = 0.0; } // add/change something so if there is an "error" off run first, after an on run, we get the next time a few secs later
              if(abs(beforeADC2) > 50){ cout << "ERROR!!! before ADC for " << tee2.Data() << " = " << beforeADC2 << endl;   beforeADC2 = 0.0; }         
              if(abs(afterADC2) > 50){ cout << "ERROR!!! after ADC for " << tee2.Data() << " = " << afterADC2 << endl;      afterADC2 = 0.0; }
              
              cout << "check off times   before: " << beforeRun1->AsSQLString() << "   after: " << afterRun1->AsSQLString() << "    2nd tel    before: " << beforeRun2->AsSQLString() << "   after: " << afterRun2->AsSQLString() << endl;
              cout << "ADCs    before: " << beforeADC1 << "   after: " << afterADC1 << "     2nd tel    before: " << beforeADC2 << "   after: " << afterADC2 << endl;
              
              // check - put out a warning if the start or end times are more than 5 mins different
              if(abs(int(beforeRun1->Convert()) - int(beforeRun2->Convert())) > 300){ cout << endl << "big discrepancy between telescopes in off runs before! check for missing files" << endl; } // must explicitly type ints or abs() gets confused
              if(abs(int(afterRun1->Convert()) - int(afterRun2->Convert())) > 300){ cout << endl << "big discrepancy between telescopes in off runs after! check for missing files" << endl; }
              
              // I used to do - if off run is missing, skip the off subtraction --> changing to if off run is missing, assume ADC constant 
              //if(beforeADC1 == 0.0 || afterADC1 == 0.0){ cout << "missing off file - skipping off subtraction for 1st telescope" << endl;  beforeADC1 = 0;     afterADC1 = 0;}
              //if(beforeADC2 == 0.0 || afterADC2 == 0.0){ cout << "missing off file - skipping off subtraction for 2nd telescope" << endl;  beforeADC2 = 0;     afterADC2 = 0;}
              if(beforeADC1 == 0.0 && afterADC1 != 0.0){ cout << "missing off file - before for 1st telescope, assuming constant background" << endl;  beforeADC1 = afterADC1;}
              else if(afterADC1 == 0.0 && beforeADC1 != 0.0){ cout << "missing off file - after for 1st telescope, assuming constant background" << endl;  afterADC1 = beforeADC1;}
              if(beforeADC2 == 0.0 && afterADC2 != 0.0){ cout << "missing off file - before for 2nd telescope, assuming constant background" << endl;  beforeADC2 = afterADC2;}
              else if(afterADC2 == 0.0 && beforeADC2 != 0.0){ cout << "missing off file - after for 2nd telescope, assuming constant background" << endl;  afterADC2 = beforeADC2;}
              
              // ok - now we need to calculate what the off ADCs would be at the beginning and end of the run - instead of assuming before ADC = start ADC and after ADC = end ADC - should be a small time difference and very small effect but worth doing it right
              // first - we'll declare histograms to hold each bin of off data - we use histograms bc it IS binned/equal finite intervals and it is easier to get bin content 
              double ninc1 = (afterRun1->Convert() - beforeRun1->Convert());     // /cfOrig->GetYaxis()->GetBinWidth(1); // should this be +15 secs for midpt? (half of 30 sec off run) doesnt matter for length but DOES matter for dif btwn start/end time
              double ninc2 = (afterRun2->Convert() - beforeRun2->Convert());     // /cfOrig->GetYaxis()->GetBinWidth(1); // total length of time / size of units (frames) = number of units -> just doing in seconds!! NOT units!
              TH1D* fullADC1 = new TH1D("fullADC1","off to off ADC T1", ninc1+1, -0.5, ninc1+0.5); // check that this is really the correct logic for making the CENTER of the bins the correct values!!
              TH1D* fullADC2 = new TH1D("fullADC2","off to off ADC T2", ninc2+1, -0.5, ninc2+0.5);
              double dOffFulladc1(0.0), dOffFulladc2(0.0); //thisOffadc1(0.0), thisOffadc2(0.0);
              dOffFulladc1 = (afterADC1 - beforeADC1)/(afterRun1->Convert() - beforeRun1->Convert());
              dOffFulladc2 = (afterADC2 - beforeADC2)/(afterRun2->Convert() - beforeRun2->Convert());
              
              cout << "convert times, before run 1: " << beforeRun1->AsSQLString() << "  " << beforeRun1->Convert() << "     after run 1: " << afterRun1->AsSQLString() << "  " << afterRun1->Convert() << endl;
              cout << "before run 2: " << beforeRun2->AsSQLString() << "  " << beforeRun2->Convert() << "      after run 2: " << afterRun2->AsSQLString() << "  " << afterRun2->Convert() << endl;
              //cout << "ninc1: " << ninc1 << "   ninc2: " << ninc2 << endl;
              
              // fill hist w adc values for full run length
              for(int io=1; io<=ninc1+1; io++){   fullADC1->SetBinContent(io, beforeADC1+(dOffFulladc1*(io-1))); }
              for(int io=1; io<=ninc2+1; io++){   fullADC2->SetBinContent(io, beforeADC2+(dOffFulladc2*(io-1))); }
              
              // get adc values from bin with time matching start/end times to get start/end of run adc values (down to second, not FRAME (~2s) increments)
              double runStartADC1(0.0), runEndADC1(0.0), runStartADC2(0.0), runEndADC2(0.0);
              runStartADC1 = fullADC1->GetBinContent((timeStart->Convert() - beforeRun1->Convert()) +1); 
              //cout << endl << "start time: " << timeStart->Convert() << "   matched bin: " << double(fullADC1->GetBinCenter((timeStart->Convert()-beforeRun1->Convert())+1)) << endl;
              runEndADC1 = fullADC1->GetBinContent(fullADC1->GetXaxis()->GetNbins() - (afterRun1->Convert() - timeEnd->Convert()));
              //cout << "end time: " << timeEnd->Convert() << "   matched bin: " << double(fullADC1->GetBinCenter(fullADC1->GetXaxis()->GetNbins() - (afterRun1->Convert() - timeEnd->Convert()))) << endl;
              runStartADC2 = fullADC2->GetBinContent((timeStart->Convert() - beforeRun2->Convert()) +1);
              runEndADC2 = fullADC2->GetBinContent(fullADC2->GetXaxis()->GetNbins() - (afterRun2->Convert() - timeEnd->Convert()));
              
              cout << endl << endl << "off before " << beforeADC1 << "  run start " << runStartADC1 << "  run end " << runEndADC1 << "   off after " << afterADC1 << endl << endl;
              
              // get avg change in ADC over run for subtraction
              dOffadc1 = (runEndADC1 - runStartADC1)/double(cfOrig->GetYaxis()->GetNbins());       //(afterADC1 - beforeADC1)/double(cfOrig->GetYaxis()->GetNbins()); OLD method - approx start and end of run as SAME as off runs
              dOffadc2 = (runEndADC2 - runStartADC2)/double(cfOrig->GetYaxis()->GetNbins());       //(afterADC2 - beforeADC2)/double(cfOrig->GetYaxis()->GetNbins());
              cout << "d off ADC1 " << dOffadc1 << "   d off ADC2 " << dOffadc2 << endl;
              
              TH1D* storeADC1 = new TH1D("storeADC1","off adc tel 1",2,0,2);
              TH1D* storeADC2 = new TH1D("storeADC2","off adc tel 2",2,0,2);
              storeADC1->SetBinContent(1,beforeADC1);       storeADC1->SetBinContent(2,afterADC1);
              storeADC2->SetBinContent(1,beforeADC2);       storeADC2->SetBinContent(2,afterADC2);
              storeADC1->Write("offADCs1");
              storeADC2->Write("offADCs2");
          }

	  // ----------------------------------------------------------- normalization -------------------------------------------------------------
	  
	  // get ADCs and normalize correlation function by their product (loop through each frame to normalize)
	  TH2D* ADC1 = new TH2D;        TH2D* ADC2 = new TH2D;   
	  zippedFile->GetDirectory("Singles")->GetObject(Form("singlesT%s", tee1.Data()), ADC1);          zippedFile->GetDirectory("Singles")->GetObject(Form("singlesT%s", tee2.Data()), ADC2); 
	  TH1D* avgADC1 = new TH1D("avgADC1",Form("avg ADC T%s", tee1.Data()), ADC1->GetYaxis()->GetNbins(), ADC1->GetYaxis()->GetXmin(), ADC1->GetYaxis()->GetXmax());     
	  TH1D* avgADC2 = new TH1D("avgADC2",Form("avg ADC T%s", tee2.Data()), ADC2->GetYaxis()->GetNbins(), ADC2->GetYaxis()->GetXmin(), ADC2->GetYaxis()->GetXmax());
	  TProfile2D* cfNorm = new TProfile2D("cfNorm","normalized correlation function %s;relative time (ns);time in run (s)",cfOrig->GetXaxis()->GetNbins(),cfOrig->GetXaxis()->GetXmin(),cfOrig->GetXaxis()->GetXmax(),cfOrig->GetYaxis()->GetNbins(),cfOrig->GetYaxis()->GetXmin(),cfOrig->GetYaxis()->GetXmax());
	  
	  
	  // divide correlation function by product of ADCs at each frame
	  if(off == 0){
	      for (int iy=1; iy<=cfOrig->GetYaxis()->GetNbins(); iy++){
	          TH1D* adc1y = ADC1->ProjectionX("first ADC proj",iy,iy);     
	          for(int ix=0; ix<5; ix++){    adc1y->SetBinContent(adc1y->GetXaxis()->GetNbins()-ix, 0); } // handle overflow positive bins by setting the last 5 to zero (maybe better way?)
	          double t1mean = adc1y->GetMean();
	          avgADC1->SetBinContent(iy,t1mean);
	          TH1D* adc2y = ADC2->ProjectionX("second ADC proj",iy,iy);  
	          for(int ix=0; ix<5; ix++){    adc2y->SetBinContent(adc2y->GetXaxis()->GetNbins()-ix ,0); }
	          double t2mean = adc2y->GetMean();
	          avgADC2->SetBinContent(iy,t2mean);
	          
	          for (int ix=1; ix<=cfOrig->GetXaxis()->GetNbins(); ix++){
	              cfNorm->Fill(cfOrig->GetXaxis()->GetBinCenter(ix), cfOrig->GetYaxis()->GetBinCenter(iy), cfOrig->GetBinContent(ix,iy)/(t1mean*t2mean));
	          }
	          delete adc1y, adc2y;
	      }
	  }
	  // for off data, divide by product of ADCs with off subtraction at each frame
	  else if(off == 1){ 
	      double offadc1 = beforeADC1;    double offadc2 = beforeADC2;
	      for (int iy=1; iy<=cfOrig->GetYaxis()->GetNbins(); iy++){
	          TH1D* adc1y = ADC1->ProjectionX("first ADC proj",iy,iy);     
	          for(int ix=0; ix<5; ix++){    adc1y->SetBinContent(adc1y->GetXaxis()->GetNbins()-ix, 0); }  
	          double t1mean = adc1y->GetMean();
	          avgADC1->SetBinContent(iy,t1mean);
	          TH1D* adc2y = ADC2->ProjectionX("second ADC proj",iy,iy); 
	          for(int ix=0; ix<5; ix++){    adc2y->SetBinContent(adc2y->GetXaxis()->GetNbins()-ix, 0); }  
	          double t2mean = adc2y->GetMean();
	          avgADC2->SetBinContent(iy,t2mean);
	          
	          for (int ix=1; ix<=cfOrig->GetXaxis()->GetNbins(); ix++){
	              cfNorm->Fill(cfOrig->GetXaxis()->GetBinCenter(ix), cfOrig->GetYaxis()->GetBinCenter(iy), (cfOrig->GetBinContent(ix,iy) - (t1mean*offadc2) - (t2mean*offadc1) + (offadc1*offadc2))/((t1mean-offadc1)*(t2mean-offadc2)));
	          }
	          offadc1 += dOffadc1;   offadc2 += dOffadc2;
	          //cout << "check " << offadc1 << "  " << offadc2 << endl;
	          delete adc1y, adc2y;
	      }
	  }
	  cfNorm->Write("normalizedCF");
	  avgADC1->Write(Form("avgADCT%swhole", tee1.Data()));
	  avgADC2->Write(Form("avgADCT%swhole", tee2.Data()));
	  
	  // josie trying to understand what is going wrong with the normalization
	  /*TH1D* avgslice = new TH1D("avgslice","avg of each frame after ADC normalization",cfNorm->GetYaxis()->GetNbins(),0,cfNorm->GetYaxis()->GetNbins());
	  for(int ib=1; ib<=cfNorm->GetYaxis()->GetNbins(); ib++){
	      double sumslice(0);
	      for(int ix=1; ix<=cfNorm->GetXaxis()->GetNbins(); ix++){
	          sumslice += cfNorm->GetBinContent(ix,ib);
	      }
	      avgslice->SetBinContent(ib, sumslice/double(cfNorm->GetXaxis()->GetNbins()));
	  }
	  avgslice->Write("normslices");*/
	  
	  // -------------------------------------------------------------- noise removal -----------------------------------------------------------

	  // rebin for better noise removal
	  //cfNorm->RebinY(16); // ~2s frames -> ~32s frames
	  
	  cfNorm->RebinY(rebin); // testing this out w peak widths
	  
	  // and save after rebinning to compare
	  /*TH1D* avgslice3 = new TH1D("avgslice3","avg of each frame after ADC normalization",cfNorm->GetYaxis()->GetNbins(),0,cfNorm->GetYaxis()->GetNbins());
	  for(int ib=1; ib<=cfNorm->GetYaxis()->GetNbins(); ib++){
	      double sumslice(0);
	      for(int ix=1; ix<=cfNorm->GetXaxis()->GetNbins(); ix++){
	          sumslice += cfNorm->GetBinContent(ix,ib);
	      }
	      avgslice3->SetBinContent(ib, sumslice/double(cfNorm->GetXaxis()->GetNbins()));
	  }
	  avgslice3->Write("rebinslices");*/

	  // FFTs and noise removal
	  cfNorm->Draw("COLZ");
	  
	  double fft = DoFFT(cfNorm, 1);
	  
	  //cfNorm = subtractOne(cfNorm); // subtract one so cf is centered around zero
	  //cfNorm->Write("cfMinusOne");
	  
	  cfNorm = forcedNormSubtOne(cfNorm, 1); // normalize so that off-peak avg is truly zero
	  cfNorm->Write("cfForcedNorm");
	  
	  // identify 80 mhz noisy frames
	  TProfile2D* PS1 = new TProfile2D;     zippedFile->GetDirectory("Singles")->GetObject(Form("PowerSpectrumT%s", tee1.Data()), PS1);
	  TProfile2D* PS2 = new TProfile2D;     zippedFile->GetDirectory("Singles")->GetObject(Form("PowerSpectrumT%s", tee2.Data()), PS2);
	  //TH1D* frameKey1 = identifyNoise(PS1,651);
	  //TH1D* frameKey2 = identifyNoise(PS2,651);
	  TH1D* frameKey = identifyNoise(PS1, PS2, 651); // ok this is where we should go next - combine the power spectra for a higher noise threshold
	  /*TH1D* frameKeyBoth = new TH1D("frameKeyBoth","pair frame key", frameKey1->GetXaxis()->GetNbins(),0,frameKey1->GetXaxis()->GetNbins());
	  for(int iff=1; iff<=frameKey1->GetXaxis()->GetNbins(); iff++){
	      if(frameKey1->GetBinContent(iff) > 0.99 || frameKey2->GetBinContent(iff) > 0.99){     frameKeyBoth->SetBinContent(iff,1); }
	      else { frameKeyBoth->SetBinContent(iff, 0); }
	  }*/
	  frameKey->Write("frameKeyNoiseOnly");
	  
	  // get keys for this run/pair of when tracking is bad, from dips in avg ADCs
	  if(rtf == 1){
	      string namedate = filename.Data();
	      namedate = namedate.substr(namedate.find("_y")+1, namedate.length()); 
	      namedate = namedate.substr(0, namedate.find("_")); 
	      TFile* trackingfile = new TFile("/fs/ess/PAS1977/StellarAnalysis/gamCas/badTrackingKeys.root","READONLY");
	      TH1D* trackingKey1 = new TH1D;    trackingfile->GetDirectory(namedate.c_str())->GetObject(Form("key%s",tee1.Data()), trackingKey1);
	      TH1D* trackingKey2 = new TH1D;    trackingfile->GetDirectory(namedate.c_str())->GetObject(Form("key%s",tee2.Data()), trackingKey2);
	      TH1D* trackingKey = new TH1D("trackingKey","bad tracking key",trackingKey1->GetXaxis()->GetNbins(),0,trackingKey1->GetXaxis()->GetNbins());
	      for(int it=1; it<=trackingKey1->GetXaxis()->GetNbins(); it++){
	          if(trackingKey1->GetBinContent(it) > 0.5 || trackingKey2->GetBinContent(it) > 0.5){ 
	              trackingKey->SetBinContent(it,1); 
	              frameKey->Fill(it, 1);//SetBinContent(it, frameKey->GetBinContent(it) + 1);
	              avgADC1->SetBinContent(it,0);
	              avgADC2->SetBinContent(it,0);
	          }
	      }
	  }
	  
	  outfile->cd(Form("%s",keyName.Data()));
	  frameKey->Write("frameKey");
	  avgADC1->Write(Form("avgADCT%s", tee1.Data()));
	  avgADC2->Write(Form("avgADCT%s", tee2.Data()));
	  
	  for(int ix=1; ix<=frameKey->GetXaxis()->GetNbins(); ix++){
	      if(frameKey->GetBinContent(ix) > 0.99){
	          avgADC1->SetBinContent(ix,0);
	          avgADC2->SetBinContent(ix,0);
	      }
	  }
	  avgADC1->Write(Form("avgADCT%srem", tee1.Data()));
	  avgADC2->Write(Form("avgADCT%srem", tee2.Data()));
	  
	  frameKey->Rebin(rebin); // for regular hist and NOT a profile, rebinning adds counts -> threshold = rebin (rebin*1)
	  /*if(rebin > 1){ // normalize frame key if needed
	      for(int i=1; i<=frameKey->GetXaxis()->GetNbins(); i++){   frameKey->SetBinContent(i, frameKey->GetBinContent(i));///double(rebin)); }
	  }*/
	  
	  //cfNorm = NoiseRemoveExactFreq(cfNorm, 0.4991, 1); // 79 MHz
	  
	  // remove 80 mhz noise - setting all frames to 0
	  /*TProfile2D* cfNormTemp = new TProfile2D;      cfNormTemp = (TProfile2D*)cfNorm->Clone();
	  cfNorm->Reset();
	  //("cfNoiseRem","CF removed 80 MHz frames",cfNorm->GetXaxis()->GetNbins(),cfNorm->GetXaxis()->GetXmin(),cfNorm->GetXaxis()->GetXmax(),cfNorm->GetYaxis()->GetNbins(),cfNorm->GetYaxis()->GetXmin(),cfNorm->GetYaxis()->GetXmax());
	  for(int iy=1; iy<=cfNormTemp->GetYaxis()->GetNbins(); iy++){
	      for(int ix=1; ix<=cfNormTemp->GetXaxis()->GetNbins(); ix++){
	          if(frameKeyBoth->GetBinContent(iy) < 0.99){  cfNorm->Fill(cfNormTemp->GetXaxis()->GetBinCenter(ix), cfNormTemp->GetYaxis()->GetBinCenter(iy), cfNormTemp->GetBinContent(ix,iy)); }
	          else { cfNorm->Fill(cfNormTemp->GetXaxis()->GetBinCenter(ix), cfNormTemp->GetYaxis()->GetBinCenter(iy), 0); }
	      }
	  }*/
	  
	  cfNorm = removeFrames(cfNorm, frameKey);
	  cfNorm->Write("cfRemoveFrames");
	  
	  /*// checking normalization - after forced norm
	  TH1D* avgslice2 = new TH1D("avgslice2","avg of each frame after forced normalization",cfNorm->GetYaxis()->GetNbins(),0,cfNorm->GetYaxis()->GetNbins());
	  for(int ib=1; ib<=cfNorm->GetYaxis()->GetNbins(); ib++){
	      double sumslice(0);
	      for(int ix=1; ix<=cfNorm->GetXaxis()->GetNbins(); ix++){
	          sumslice += cfNorm->GetBinContent(ix,ib);
	      }
	      avgslice2->SetBinContent(ib, sumslice/double(cfNorm->GetXaxis()->GetNbins()));
	  }
	  avgslice2->Write("forcedslices");*/
	  
	  /*fft = DoFFT(cfNorm, 2); 
	  cfNorm = NoiseRemoveExactFreq(cfNorm, 0.0628, 2); // 10 MHz*/
	  fft = DoFFT(cfNorm, 2); 
	  cfNorm = NoiseRemoveExactFreq(cfNorm, 0.6754, 2); // 107 MHz
	  fft = DoFFT(cfNorm, 3); 
	  //cfNorm = NoiseRemoveExactFreq(cfNorm, 0.6280, 4);
	  //cfNorm = NoiseRemoveExactFreq(cfNorm, 0.5887, 3); // 93.69 MHz
	  //fft = DoFFT(cfNorm, 4);
	  //cfNorm = NoiseRemoveExactFreq(cfNorm, 0.6255, 4); // 99.55 MHz
	  //fft = DoFFT(cfNorm, 5);
	  //cfNorm = NoiseRemoveExactFreq(cfNorm, 0.5726, 6); // 91.125 MHz - no difference (for this data)
	  //fft = DoFFT(cfNorm, 5);
	  //cfNorm = NoiseRemoveExactFreq(cfNorm, 0.5963, 5); // 94.91 MHz
	  //fft = DoFFT(cfNorm, 6);*/
	  
	  cfNorm = forcedNormSubtOne(cfNorm, 0, frameKey);
	  cfNorm->Write("cfForcedNorm2");

	  // draw the heatmap
	  cfNorm->SetTitle("heatmap;relative time (ns);time in run (s)");
	  cfNorm->Write("heatmap");
	  
	  //--------------------------------------------------------- OPD shift and time average -----------------------------------------------------

	  // read in delays and baselines 
	  double cableDelays[5]    = {676.8, 676.8, 585.0, 955.0, 1063.7};
	  double refDelay;   refDelay = -1*(nbo*4.0 + cableDelays[t1] - cableDelays[t2]); // seg fault if not defined separately
	  TGraph* opd;                     TGraph* baseline;
	  TGraph* opdGraph = new TGraph;   TGraph* baselineGraph = new TGraph;
	  TSpline3* opdSpline;             TSpline3* baseSpline;
	  TProfile* cfFlat = new TProfile;

	  // skip opd shift, baselines, etc for T0T1 
	  if(t1 == 0 && t2 == 1){
	    cfFlat = cfNorm->ProfileX("",1,-1); 
	  }
	  else {
	    geom->Draw("OPD:time");
	    opd = new TGraph(geom->GetSelectedRows(), geom->GetV1(), geom->GetV2());
	    //opd->Write("opd");
	    
	    for(int r=0; r<opd->GetN(); r++){ opd->SetPointX(r,opd->GetPointX(r) - refDelay);}
	    opd->Write("opd");

	    // turn into splines and evaluate at each knot to make new TGraph
	    //opdSpline = new TSpline3("opdSpline",opd);
	    //opdGraph->Write("opdgraph");

	    // make new histogram to hold opd shifted cf and fill by opd shift
	    TProfile2D* cfShift = new TProfile2D("cfShift","OPD shifted correlation function;relative time (ns); time in run (s)", cfNorm->GetXaxis()->GetNbins()*4.0, cfNorm->GetXaxis()->GetXmin(),
					       cfNorm->GetXaxis()->GetXmax(), cfNorm->GetYaxis()->GetNbins(), cfNorm->GetYaxis()->GetXmin(), cfNorm->GetYaxis()->GetXmax());
	    for (int iy=1; iy<=cfNorm->GetYaxis()->GetNbins(); iy++){
	        for (int ix=1; ix<=cfNorm->GetXaxis()->GetNbins(); ix++){
	            double xshift = cfNorm->GetXaxis()->GetBinCenter(ix) - opd->GetPointX(iy-1); // iy - 1 ??? -- 2 sept 24 john pointed out I should check this
	            cfShift->Fill(xshift, cfNorm->GetYaxis()->GetBinCenter(iy), cfNorm->GetBinContent(ix,iy));
	        }
	    }
	    cfShift = removeFrames(cfShift, frameKey);
	    
	    // project shifted correlation function
	    cfFlat = cfShift->ProfileX("",1,-1);
	    cfShift->Write("shiftedCF");
	  }
	  
	  cfFlat->SetTitle("projected correlation function;relative time(ns);g^2");
	  cfFlat->Write("projectedCF");
	  
	  // ------------------------------------------------------- HBT peak fit -----------------------------------------------------------
	  
	  // now fit HBT peak
	  // define fit parameters
	  double reltimePar(-5);
	  double rtWindow(10.0); // HALF width of search window for hbt peak
	  double sigmaMin(2.5), sigmaMax(5.5); // ns, min and max width of hbt peak
	  //if((t1==1 && t2==2)){reltimePar = -20;} if((t1==2 && t2==1)){reltimePar =  20;} // need to adjust these??
	  //if((t1==1 && t2==3)){reltimePar = -10;} if((t1==3 && t2==1)){reltimePar =  10;}
	  //if((t1==1 && t2==4)){reltimePar =  2;}  if((t1==4 && t2==1)){reltimePar = -2;}
	  //if((t1==2 && t2==3)){reltimePar =  10;} if((t1==3 && t2==2)){reltimePar = -10;}
	  //if((t1==2 && t2==4)){reltimePar =  18;} if((t1==4 && t2==2)){reltimePar = -18;}
	  //if((t1==3 && t2==4)){reltimePar =  8;}  if((t1==4 && t2==3)){reltimePar = -8;}
	  if((t1==1 && t2==2)){reltimePar = -5;} if((t1==2 && t2==1)){reltimePar =  20;} // need to adjust these??
	  if((t1==1 && t2==3)){reltimePar = -10;} if((t1==3 && t2==1)){reltimePar =  10;}
	  if((t1==1 && t2==4)){reltimePar =  2;}  if((t1==4 && t2==1)){reltimePar = -2;}
	  if((t1==2 && t2==3)){reltimePar =  10;} if((t1==3 && t2==2)){reltimePar = -10;}
	  if((t1==2 && t2==4)){reltimePar =  18;} if((t1==4 && t2==2)){reltimePar = -18;}
	  if((t1==3 && t2==4)){reltimePar =  8;}  if((t1==4 && t2==3)){reltimePar = -8;}
	  
	  TF1* hbtfit = new TF1("hbtfit","([0]*exp(-pow((x-[1])/[2],2)/2.0))/([2]*sqrt(2*pi))",-256,256);
	  hbtfit->SetParName(0,"area");   hbtfit->SetParameter(0,0.0);
	  hbtfit->SetParName(1,"#tau_{0}"); hbtfit->SetParameter(1,reltimePar);  hbtfit->SetParLimits(1,reltimePar-rtWindow, reltimePar+rtWindow);
	  hbtfit->SetParName(2,"#sigma");  hbtfit->SetParameter(2,4.0);  hbtfit->SetParLimits(2,sigmaMin, sigmaMax);
	  
	  gStyle->SetOptStat(0);  gStyle->SetOptFit(1);
	  cfFlat->Fit("hbtfit");
	  gPad->Update(); gPad->Modified(); // get fit results to show up in stats box
	  cout << t1 << t2 << "   area: " << hbtfit->GetParameter(0) << "   tau: " << hbtfit->GetParameter(1) << "   sigma: " << hbtfit->GetParameter(2) << endl;
	  cfFlat->Write("fittedCF");
	  outfile->cd();
  
      }// end of loop over pairs
  }// end of loop over keys
  
  // ========================================================== save ADCs and power spectra =============================================================

  // now save all ADCs and power spectra to the output file
  outfile->mkdir("Singles");   outfile->cd("Singles");
  
  TList* singlesList = zippedFile->GetDirectory("Singles")->GetListOfKeys();
  TH2D* adc = new TH2D;
  TProfile2D* powspec = new TProfile2D;
  
  for (int ia=0; ia<singlesList->GetSize(); ia++){
      TKey* key = (TKey*)singlesList->At(ia);
      TString keyName = key->GetName();
      TString className = key->GetClassName();
      if((keyName.Data()[0] == 's') && (className == "TH2D")){
          zippedFile->GetDirectory("Singles")->GetObject(keyName.Data(),adc);     adc->SetTitle(Form("ADC T%c",keyName.Data()[8]));
          adc->Write(Form("ADCT%c", keyName.Data()[8]));
      }
      else if((keyName.Data()[0] == 'P') && (className == "TProfile2D")){
          zippedFile->GetDirectory("Singles")->GetObject(keyName.Data(),powspec);  powspec->SetTitle(Form("power spectrum T%c", keyName.Data()[14]));
          powspec->Write(Form("PowerSpectrumT%c", keyName.Data()[14]));
      }
  }
  
  // save off runs to output file also
  outfile->mkdir("OFF");    outfile->cd("OFF");
  TList* offList = zippedFile->GetDirectory("OFF")->GetListOfKeys();
  TTree* offTree;
  
  for(int is=0; is<offList->GetSize(); is++){
      TKey* key = (TKey*)offList->At(is);
          TString keyName = key->GetName();
          zippedFile->GetDirectory("OFF")->GetObject(keyName.Data(),offTree);
          outfile->cd("OFF");
          offTree->CloneTree()->Write();
  }
  
  outfile->cd();
  outfile->Close();
      
}// end of macro!


// ========================================================== FUNCTIONS ============================================================================

// FFT function

double DoFFT (TProfile2D * CF, int times){
  //cout << "Do an FFT of the data\n\n";
  int nPowerBins = CF->GetXaxis()->GetNbins()/2;
  double binsize = CF->GetXaxis()->GetBinWidth(1);   // multiply by 4.0e-9 to make seconds 
  double nyquistFrequency_RadiansPerNs = 2.0*TMath::Pi()*(0.5/binsize);
 //double nyquistFrequency_Hertz        = nyquistFrequency_RadianPerSec / (2.0*TMath::Pi());  //to make hertz
  
  TH2D * FourierTransform = new TH2D(Form("FFT%d", times) , "FFT", nPowerBins, 0.0, nyquistFrequency_RadiansPerNs, CF->GetYaxis()->GetNbins(), CF->GetYaxis()->GetXmin(),CF->GetYaxis()->GetXmax());
  
  for(int TimeSlice = 1; TimeSlice<=CF->GetYaxis()->GetNbins(); TimeSlice++){
    TH1D * CFy = CF->ProjectionX("temp", TimeSlice, TimeSlice);
    int nbins = CF->GetXaxis()->GetNbins();
    TVirtualFFT * myFFT = TVirtualFFT::FFT(1, &nbins, "DHT");
    TH1D * Results = (TH1D *)CFy->FFT(0, "MAG");
    for(int i=2; i<=nPowerBins; i++){     // note that I am NOT INCLUDING BIN #1, which has the DC component
      double W = 2.0*Results->GetBinContent(i)/double(Results->GetXaxis()->GetNbins());
      FourierTransform->SetBinContent(i,TimeSlice,W);
    }
    delete Results;
  }
 
  TH1D * FFTCombine = FourierTransform->ProjectionX("", 1, -1, "");
  FFTCombine->SetTitle("Fourier Transform; Frequency (Radians per ns); Magnitude");
  FFTCombine->Write(Form("FFTProjection%d",times));
  
  return FFTCombine->GetBinCenter(FFTCombine->GetMaximumBin());
}

// ==================================================================================================================================================

// noise removal function

TProfile2D * NoiseRemoveExactFreq(TProfile2D *CF, double FFTFreq, int times){
    
  TProfile2D* CFMinusNoise = (TProfile2D*)CF->Clone("CFRemoveNoise"); 
  CFMinusNoise->Reset();
  
  //Define graphs to save parameters to
  TGraph * HeightGraph = new TGraph(); HeightGraph->SetTitle("DC Component; Time (Seconds); Normalized Counts");
  TGraph * AmpGraph    = new TGraph(); AmpGraph->SetTitle("Amplitude of Noise; Time (Seconds); Amplitude");
  TGraph * FreqGraph   = new TGraph(); FreqGraph->SetTitle("Frequency of Noise; Time (Seconds); Frequency (radians per ns)");
  TGraph * PhaseGraph  = new TGraph(); PhaseGraph->SetTitle("Phase of Noise; Time(Seconds); Phase (Radians)");
 
  for(int TimeSlice = 1; TimeSlice<=CF->GetYaxis()->GetNbins() ; TimeSlice++){ 
    TH1D * CFy = CF->ProjectionX("temp", TimeSlice, TimeSlice);
    double height = 0.0;
    //if (times == 1){
    //    CFy->Fit("pol0", "S0Q");
    //    height = CFy->GetFunction("pol0")->GetParameter(0);
    //}
    
    double AmpCos(0), CosNorm(0), AmpSin(0), SinNorm(0);
    for(int ix = 1; ix <= CFy->GetXaxis()->GetNbins(); ++ix){
      if(CFy->GetBinCenter(ix) < -100 || CFy->GetBinCenter(ix) > 100){ 
        AmpCos += cos(FFTFreq*CFy->GetBinCenter(ix))*(CFy->GetBinContent(ix));// -height);
        CosNorm += pow(cos(FFTFreq*CFy->GetBinCenter(ix)), 2);
        AmpSin += sin(FFTFreq*CFy->GetBinCenter(ix))*(CFy->GetBinContent(ix)); //-height);
        SinNorm += pow(sin(FFTFreq*CFy->GetBinCenter(ix)), 2);
      }
    }
    
    AmpCos /= CosNorm; AmpSin /= SinNorm;

    double Amp = sqrt(pow(AmpCos,2)+pow(AmpSin,2));
    double Phi = atan2(AmpSin, AmpCos);
    
    HeightGraph->AddPoint(CF->GetYaxis()->GetBinCenter(TimeSlice), height);
    AmpGraph->AddPoint(CF->GetYaxis()->GetBinCenter(TimeSlice), Amp);
    FreqGraph->AddPoint(CF->GetYaxis()->GetBinCenter(TimeSlice), 0.0);
    PhaseGraph->AddPoint(CF->GetYaxis()->GetBinCenter(TimeSlice), Phi);
    
    for (int ix=1; ix<=CF->GetXaxis()->GetNbins(); ix++){
      CFMinusNoise->Fill(CF->GetXaxis()->GetBinCenter(ix),CF->GetYaxis()->GetBinCenter(TimeSlice), CF->GetBinContent(ix, TimeSlice) - (Amp*cos(FFTFreq*CF->GetXaxis()->GetBinCenter(ix)-Phi))); //height+Amp*cos..
    }
  }
  
  //Save Parameter Graphs to a canvas and write to output file 
  TCanvas * AllParams = new TCanvas(Form("NoiseParameters%d", times), "AllParams", 1600, 800);
  AllParams->Divide(2,2);
  AllParams->cd(1); HeightGraph->Draw("A*");
  AllParams->cd(2); AmpGraph->Draw("A*");
  AllParams->cd(3); FreqGraph->Draw("A*");
  AllParams->cd(4); PhaseGraph->Draw("A*");
  AllParams->Write(Form("ParameterGraphs%d", times));
  
  
  //Write and return Correlation Functions without noise 
  CFMinusNoise->SetTitle("Normalized Correlation Function Without Noise; Relative Time (ns); Time");
  CFMinusNoise->Write(Form("CFMinusNoise%d", times));
  
  return CFMinusNoise;
}

// ==========================================================================================================================

// forced normalization / subtract 1 function --> NO LONGER DOING FORCED NORM HERE - josie 6 nov 2024

TProfile2D* subtractOne(TProfile2D* cf){
    
    // new cf to fill
    TProfile2D* cfSubtOne = new TProfile2D("cfSubtOne","cf subtract one;relative time (ns);time in runs (s)", cf->GetXaxis()->GetNbins(), cf->GetXaxis()->GetXmin(), cf->GetXaxis()->GetXmax(), 
                            cf->GetYaxis()->GetNbins(), cf->GetYaxis()->GetXmin(), cf->GetYaxis()->GetXmax());
    
    for(int iy=1; iy<=cf->GetYaxis()->GetNbins(); iy++){
        for(int ix=1; ix<=cf->GetXaxis()->GetNbins(); ix++){
            cfSubtOne->Fill(cf->GetXaxis()->GetBinCenter(ix), cf->GetYaxis()->GetBinCenter(iy), cf->GetBinContent(ix,iy) - 1.0);
        }
    }
                            
    
    // trying both pol0 fit and simple avg to compare
    /*for(int iy=1; iy<=cf->GetYaxis()->GetNbins(); iy++){
        double pol0fit(0), regavg(0), sum(0);
        int navg(0);
        TH1D* tempslice = new TH1D;
        tempslice = cf->ProjectionX("",iy,iy);
        //if(times == 1){tempslice->Fit("pol0", "S0Q"); pol0fit = tempslice->GetFunction("pol0")->GetParameter(0); }
        for(int ix=1; ix<=tempslice->GetXaxis()->GetNbins(); ix++){
            //if(cf->GetXaxis()->GetBinCenter(ix) < -50 || cf->GetXaxis()->GetBinCenter(ix) > 50){ sum += tempslice->GetBinContent(ix); navg++;} // previous - as of 3 nov 2024
            if((cf->GetXaxis()->GetBinCenter(ix) > -230 && cf->GetXaxis()->GetBinCenter(ix) < -30) || (cf->GetXaxis()->GetBinCenter(ix) > 30 && cf->GetXaxis()->GetBinCenter(ix) < 230)){ sum += tempslice->GetBinContent(ix); navg++;} 
        }
        regavg = sum/((double)navg);
        
        //if(times == 1){ subtpol0->SetBinContent(iy,pol0fit); }
        //subtreg->SetBinContent(iy,regavg);
        
        for(int ix=1; ix<=cf->GetXaxis()->GetNbins(); ix++){
            cfSubtOne->Fill(cf->GetXaxis()->GetBinCenter(ix), cf->GetYaxis()->GetBinCenter(iy), cf->GetBinContent(ix,iy) - regavg);
        }
    } */
    /*if(times == 1){
        subtpol0->Write("subtpol01");
        subtreg->Write("subtreg1");
    }
    if(times == 2){
        //subtpol0->Write("subtpol02");
        subtreg->Write("subtreg2");
    }*/
    
    return cfSubtOne;
}

// ===========================================================================================================================

// forced normalization function plus subtract one

TProfile2D* forcedNormSubtOne(TProfile2D* cf, int one = 1){ // one = 1 if cf centered at 1 (needs -1), one = 0 if cf centered at 0 
    
    // new profile to fill with normalized cf
    TProfile2D* cfForcedNorm = new TProfile2D("cfForcedNorm","cf forced normalization;relative time (ns);time in runs (s)", 
        cf->GetXaxis()->GetNbins(), cf->GetXaxis()->GetXmin(), cf->GetXaxis()->GetXmax(), cf->GetYaxis()->GetNbins(), cf->GetYaxis()->GetXmin(), cf->GetYaxis()->GetXmax());
        
    // loop over each frame to calculate avg off-peak, then DIVIDE all bins in frame to normalize, center cf at zero
    for(int iy=1; iy<=cf->GetYaxis()->GetNbins(); iy++){
        double sum(0), rowavg(0);
        int navg(0);
        for(int ix=1; ix<=cf->GetXaxis()->GetNbins(); ix++){
            if(one == 1){
                if((cf->GetXaxis()->GetBinCenter(ix) > -230 && cf->GetXaxis()->GetBinCenter(ix) < -50) || (cf->GetXaxis()->GetBinCenter(ix) > 50 && cf->GetXaxis()->GetBinCenter(ix) < 230)){ sum += cf->GetBinContent(ix,iy); navg++;}
            }
            else if(one == 0){ 
                if((cf->GetXaxis()->GetBinCenter(ix) > -230 && cf->GetXaxis()->GetBinCenter(ix) < -50) || (cf->GetXaxis()->GetBinCenter(ix) > 50 && cf->GetXaxis()->GetBinCenter(ix) < 230)){ sum += 1.0+cf->GetBinContent(ix,iy); navg++;}
            }
        }
        rowavg = sum/navg;
        //if (abs(rowavg - 1.0) > 1e-8){ // check that row isn't empty - else we ignore it
        for(int ix=1; ix<=cf->GetXaxis()->GetNbins(); ix++){
            cfForcedNorm->Fill(cf->GetXaxis()->GetBinCenter(ix), cf->GetYaxis()->GetBinCenter(iy), (cf->GetBinContent(ix,iy)/rowavg) - 1.0);
        }
    }
    return cfForcedNorm;
}


TProfile2D* forcedNormSubtOne(TProfile2D* cf, int one = 1, TH1D* key = 0){ // this one uses the bad frames key - checks if frames are zero
    
    // new profile to fill with normalized cf
    TProfile2D* cfForcedNorm = new TProfile2D("cfForcedNorm","cf forced normalization;relative time (ns);time in runs (s)", 
        cf->GetXaxis()->GetNbins(), cf->GetXaxis()->GetXmin(), cf->GetXaxis()->GetXmax(), cf->GetYaxis()->GetNbins(), cf->GetYaxis()->GetXmin(), cf->GetYaxis()->GetXmax());
        
    // loop over each frame to calculate avg off-peak, then DIVIDE all bins in frame to normalize, center cf at zero
    for(int iy=1; iy<=cf->GetYaxis()->GetNbins(); iy++){
        double sum(0), rowavg(0);
        int navg(0);
        for(int ix=1; ix<=cf->GetXaxis()->GetNbins(); ix++){
            if(one == 1){
                if((cf->GetXaxis()->GetBinCenter(ix) > -230 && cf->GetXaxis()->GetBinCenter(ix) < -50) || (cf->GetXaxis()->GetBinCenter(ix) > 50 && cf->GetXaxis()->GetBinCenter(ix) < 230)){ sum += cf->GetBinContent(ix,iy); navg++;}
                //if(cf->GetXaxis()->GetBinCenter(ix) < -50 || cf->GetXaxis()->GetBinCenter(ix) < 50){ sum += cf->GetBinContent(ix,iy); navg++;}
            }
            else if(one == 0){ 
                if((cf->GetXaxis()->GetBinCenter(ix) > -230 && cf->GetXaxis()->GetBinCenter(ix) < -50) || (cf->GetXaxis()->GetBinCenter(ix) > 50 && cf->GetXaxis()->GetBinCenter(ix) < 230)){ sum += 1.0+cf->GetBinContent(ix,iy); navg++;}
                //if(cf->GetXaxis()->GetBinCenter(ix) < -50 || cf->GetXaxis()->GetBinCenter(ix) < 50){ sum += 1.0 + cf->GetBinContent(ix,iy); navg++;}
            }
        }
        rowavg = sum/navg;
        // only fill bins that are nonzero w good frames
        if(key != 0){
            if(key->GetBinContent(iy) < 0.5){
                for(int ix=1; ix<=cf->GetXaxis()->GetNbins(); ix++){
                    cfForcedNorm->Fill(cf->GetXaxis()->GetBinCenter(ix), cf->GetYaxis()->GetBinCenter(iy), (cf->GetBinContent(ix,iy)+1.0/rowavg) - 1.0);
                }
            }
        }
        else{
            for(int ix=1; ix<=cf->GetXaxis()->GetNbins(); ix++){
                cfForcedNorm->Fill(cf->GetXaxis()->GetBinCenter(ix), cf->GetYaxis()->GetBinCenter(iy), (cf->GetBinContent(ix,iy)/rowavg) - 1.0);
            }
        }
    }
    return cfForcedNorm;
}

// ===========================================================================================================================

// function to make TNtuple from python calculations stored in txt files

TNtuple* GeometryNtuple(TString keyName, int N){
  //TNtuple* geom = new TNtuple("PairGeometryInfo",Form("%s geometry",keyName.Data()),"time:baseline:OPD:u:v:w");
  TNtuple* geom = new TNtuple("PairGeometryInfo",Form("%s geometry",keyName.Data()),"frame:time:u:v:baseline:OPD:hourAngle");
  float vals[7], junk;
  //TString junk; // trying w junk as a string
  
  if (keyName == "T0T1"){                    // T0T1 is a special case
    for (int i=0; i<7; i++) vals[i]=0.0;
    for (int i=0; i<N; i++){
      geom->Fill(vals);
    }
    return geom;
  }

  // also..... T0Tx where x is "anything" is a special case..... :-(
  TString nameBase = keyName;
  nameBase.ReplaceAll("0","1");   // there is no "T0T2.delay" file :-(

  // from the file, to reference
  //frame#   time in run (s)   u coord   v coord   baseline(m)   opd(ns)
  //1   789646049.184000   -67.690173   -50.053073   84.183744   178.278883

  ifstream pyinfoFile;
  pyinfoFile.open(Form("pyinfo%s.txt", nameBase.Data()));       if (!pyinfoFile.is_open()){cout << "Problem opening python file!!!\n";  return 0;}
  //if (pyinfoFile.is_open()){ cout << "file opened ok! " << nameBase.Data() << endl;}

  pyinfoFile >> junk >> junk >> junk >> junk >> junk >> junk >> junk; 
  for (int i=0; i<N; i++){
    pyinfoFile >> vals[0] >> vals[1] >> vals[2] >> vals[3] >> vals[4] >> vals[5] >> vals[6]; // frame, time, u, v, baseline, opd, hour angle
    //if(i < 10){ cout << vals[0] << " " << vals[1] << " " << vals[2] << " " << vals[3] << " " << vals[4] << " " << vals[5] << endl;}
    geom->Fill(vals);
  }
  pyinfoFile.close();
  return geom;
}

// ===========================================================================================================================

// function to identify noisy frames to be removed

TH1D* identifyNoise(TProfile2D* powerSpectrum1, TProfile2D* powerSpectrum2, int noiseBin){  // input power spectrum and bin with noise (80 mhz = bin 651)

    // projected profile hists in clock time (y)
    TProfile* psproj1 = new TProfile;       psproj1 = powerSpectrum1->ProfileY("",noiseBin,noiseBin);
    TProfile* psproj2 = new TProfile;       psproj2 = powerSpectrum2->ProfileY("",noiseBin,noiseBin);
    
    // projected hist, "normalized" by subtracting avg of neighbors
    TH1D* localnorm = new TH1D("localnorm","normalized projection of noise in time;time (s);mag of freq",psproj1->GetXaxis()->GetNbins(),0,psproj1->GetXaxis()->GetXmax());
    for(int i=1; i<=psproj1->GetXaxis()->GetNbins(); i++){ 
        localnorm->SetBinContent(i,((psproj1->GetBinContent(i) - ((powerSpectrum1->GetBinContent(noiseBin-1, i)+(powerSpectrum1->GetBinContent(noiseBin+1, i)))/2.0)) 
            + (psproj2->GetBinContent(i) - ((powerSpectrum2->GetBinContent(noiseBin-1, i)+(powerSpectrum2->GetBinContent(noiseBin+1, i)))/2.0))));
    }
    localnorm->Write("localnormps");
    // hists to hold dist of values in projected plots - to calculate threshold of noise > normal fluctuations
    TH1D* dist = new TH1D("dist","dist to calc threshold",400,localnorm->GetMinimum(), localnorm->GetMaximum());
    for(int i=1; i<=localnorm->GetXaxis()->GetNbins(); i++){ 
        dist->Fill(localnorm->GetBinContent(i)); 
    }
    
    TF1* gausfull = new TF1("gausfull","[0]*TMath::Exp(-pow(x-[2], 2.0)/(2.0*pow([1],2.0)))",-100,100);
    gausfull->SetParName(0, "amp");         gausfull->SetParameter(0,100);
    gausfull->SetParName(1, "#sigma");      gausfull->SetParameter(1,0.00001);
    gausfull->SetParName(2, "mean");        gausfull->SetParameter(2,0);
    dist->Fit(gausfull,"Q");
    
    TF1* gauslim = new TF1("gauslim","[0]*TMath::Exp(-pow(x-[2], 2.0)/(2.0*pow([1],2.0)))",-gausfull->GetParameter(1),gausfull->GetParameter(1));
    gauslim->SetParName(0, "amp");          gauslim->SetParameter(0,100);
    gauslim->SetParName(1, "#sigma");       gauslim->SetParameter(1,0.00001);
    gauslim->SetParName(2, "mean");         gauslim->SetParameter(2,0);
    dist->Fit(gauslim,"Q");
    
    double threshold = abs(gauslim->GetParameter(1))*3.0;   // noise threshold is 3 sigma from 2nd fit - in theory should only be positive but we have at least one where it was neg
    dist->Write("tempnoisedist");
    
    // fill "key" hist to say good/bad frame
    int nbad(0);
    TH1D* frameKey = new TH1D("frameKey","noisy and good frames key",psproj1->GetXaxis()->GetNbins(),0,psproj1->GetXaxis()->GetNbins());
    for(int i=1; i<=psproj1->GetXaxis()->GetNbins(); i++){
        if(localnorm->GetBinContent(i) > threshold){ frameKey->SetBinContent(i,1); nbad++; }
        else { frameKey->SetBinContent(i,0); }
    }
    
    cout << "noisy frames found in " << (nbad/double(psproj1->GetXaxis()->GetNbins()))*100 << "% of the data" << endl;
    if((nbad/double(psproj1->GetXaxis()->GetNbins())) < 0.20 || (nbad/double(psproj1->GetXaxis()->GetNbins())) > 0.40){ cout << endl << endl << "CHECK BAD FRAME IDENTIFICATION!!!" << endl << endl << endl; }
    return frameKey;
}

// ===========================================================================================================================

// function to remove the bad frames, after they have been identified

TProfile2D* removeFrames(TProfile2D* cf, TH1D* noisekey){
    
    TProfile2D* cfnew = new TProfile2D("cfnew","cf removed frames",cf->GetXaxis()->GetNbins(),cf->GetXaxis()->GetXmin(),cf->GetXaxis()->GetXmax(),cf->GetYaxis()->GetNbins(),cf->GetYaxis()->GetXmin(),cf->GetYaxis()->GetXmax());      
    for(int i=1; i<=cf->GetYaxis()->GetNbins(); i++){
        if(noisekey->GetBinContent(i) < 0.99){
            for(int ix=1; ix<=cf->GetXaxis()->GetNbins(); ix++){
                if(cf->GetBinEntries(cf->FindBin(cf->GetXaxis()->GetBinCenter(ix),cf->GetYaxis()->GetBinCenter(i))) > 0.5){
                    cfnew->Fill(cf->GetXaxis()->GetBinCenter(ix), cf->GetYaxis()->GetBinCenter(i), cf->GetBinContent(ix,i)); 
                }
            //else { cfnew->Fill(cf->GetXaxis()->GetBinCenter(ix), cf->GetYaxis()->GetBinCenter(i), 0); }
            }
        }
    }
    return cfnew;
}
