// josie 13 nov 2023
// new full analysis script, processes as many telescope pairs as are in the file from versii, only takes version -4 for now
// with edits from dr lisa to add TNtuple, and josie edits to read in correct uvw info since i fixed the python
// 21 jan 2024 adding off data -- NOTE: takes a SECOND input, T/F off data

TNtuple* GeometryNtuple(TString keyName, int N);
double DoFFT (TProfile2D * CF, int times);
TProfile2D * NoiseRemoveExactFreq(TProfile2D *CF, double FFTFreq, int times);

void versiiAnalysisOff(TString filename){//, int off = 1){ // bool doesn't work - let's try an int
int off = 1; // idk man
    
    // opening input "zipped frames" root file
    TFile* zippedFile = new TFile(filename.Data(),"READONLY");
    
    gStyle->SetOptStat(0); // do not autofill stat box (allows info from fit later)
    gStyle->SetOptFit(1); // put fit result in stat box
    
    // get info from the MAIN header first
    TTree* header = new TTree;                    zippedFile->GetObject("Header",header);
    TString* source = new TString;                  header->SetBranchAddress("SourceBr",&source);         //cout << "star: " << source->Data() << endl;
    int version, nfiles, numpairs;    header->SetBranchAddress("VersionBr",&version);       //header->SetBranchAddress("NfilesBr",&nfiles);   header->SetBranchAddress("NpairsBr",&numpairs);
    header->GetEvent(0); // you MUST have this to retrieve correct values from headers
  
  // =================== get necessary info from headers and run python script to calculate delays ========================================
  TString runpy;
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
    if ((keyName.Data()[0]=='T')&&(className=="TDirectoryFile")){ // if we find directory that start w "T" - take it as a pair
      TTree* pairhead = new TTree;
      zippedFile->GetDirectory(keyName)->GetObject("PairHeader",pairhead);
      TProfile2D* cftemp = new TProfile2D;      zippedFile->GetDirectory(keyName)->GetObject(Form("CF%s", keyName.Data()), cftemp);
      if(cftemp->GetYaxis()->GetXmax() > longestRun){ longestRun = cftemp->GetYaxis()->GetXmax(); }
      if(cftemp->GetYaxis()->GetNbins() > mostFrames){ mostFrames = cftemp->GetYaxis()->GetNbins();  frameWidth = cftemp->GetYaxis()->GetBinWidth(5); }
      pairhead->SetBranchAddress("LocalTimeBr",&timetemp);   pairhead->GetEvent(0);
      if(timetemp->GetTime() < ltime->GetTime()){ ltime = timetemp; }
    }
  }
  // divide by rebinning factor to equal final number of frames at OPD shift
  mostFrames = mostFrames/16;
  frameWidth = frameWidth*16.0;
    
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

  TCanvas* c1 = new TCanvas;

  // make root file to save objs
  TFile* outfile = new TFile("analysisFull.root","RECREATE");

  // ===================================================== start of the main analysis ===========================================================
  
  // loop over each directory we find in the root file (taken from AnalyzeVersiiFiles.C)
  zippedFile->cd();
  TList* keyList = gDirectory->GetListOfKeys();
  for (int ip=0; ip<keyList->GetSize(); ip++){
      TKey* key = (TKey*)keyList->At(ip);
      TString keyName = key->GetName();
      TString className = key->GetClassName();
      if ((keyName.Data()[0]=='T')&&(className=="TDirectoryFile")){ // if we find directory that start w "T" - take it as a pair
          cout << "now starting pair " << keyName.Data() << endl;
          TString tee1 = keyName[1];    int t1 = atoi(tee1.Data()); 
          TString tee2 = keyName[3];    int t2 = atoi(tee2.Data());
          
          // get pair header tree, for now all we need is the correlation delay (num buckets offset) but can print the rest as a check
          TTree* pairhead = new TTree;
          zippedFile->GetDirectory(keyName)->GetObject("PairHeader",pairhead);
          int nbo;   pairhead->SetBranchAddress("corrDelayBr",&nbo);   pairhead->GetEvent(0);
          
          // get original correlation function and draw
          TProfile2D* cfOrig = new TProfile2D;      zippedFile->GetDirectory(keyName)->GetObject(Form("CF%s", keyName.Data()), cfOrig);
          cfOrig->SetTitle(Form("raw correlation function %s;relative time (ns); frames (s)", keyName.Data()));
          outfile->cd();
          gDirectory->mkdir(Form("%s",keyName.Data()));
          outfile->cd(Form("%s",keyName.Data()));
          
          // get info from python-generated text files and create TNtuple -- note: super simple to draw baselines, etc w TNtuple --> PairGeometryInfo->Draw("v:u")
          TNtuple* geom = GeometryNtuple(keyName, int(mostFrames));
          geom->Write();
          
          cfOrig->Write("originalCF");
          
          // ------------------------------------------------------------- off data ---------------------------------------------------------------
          
          // off ADCs before and after run for both telescopes, and increments for change in ADC 
          double beforeADC1(0.0), afterADC1(0.0), beforeADC2(0.0), afterADC2(0.0); 
          double dOffadc1, dOffadc2;
          
          if (off == 1){
              
              // get start and end times of run and convert to UTC time
              TDatime* timeStart = new TDatime;    TDatime* timeEnd = new TDatime;
              //int timeStartConv, timeEndConv;
              int runlength = int(cfOrig->GetYaxis()->GetXmax());
              pairhead->SetBranchAddress("LocalTimeBr",&timeStart);    pairhead->GetEvent(0);
              timeEnd->Set(timeStart->GetYear(), timeStart->GetMonth(), timeStart->GetDay(), timeStart->GetHour() + (runlength/3600), timeStart->GetMinute() + (runlength%3600/60), timeStart->GetSecond() + (runlength%3600%60));
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
                      offbranch1->SetBranchAddress("AveAdcBr",&beforeADC1);     offbranch1->GetEvent(io);
                  }
                  if((tempofftime1->Convert() > timeEnd->Convert()) && (tempofftime1->Convert() - timeEnd->Convert() < afterRun1->Convert() - timeEnd->Convert())){ 
                      afterRun1->Set(tempofftime1->GetYear(), tempofftime1->GetMonth(), tempofftime1->GetDay(), tempofftime1->GetHour(), tempofftime1->GetMinute(), tempofftime1->GetSecond());
                      offbranch1->SetBranchAddress("AveAdcBr",&afterADC1);      offbranch1->GetEvent(io);
                  }
              }
              // and repeat for the second telescope - just in case one would cut out at the beginning/end of the night
              for (int io=0; io<offbranch2->GetEntries(); io++){
                  offbranch2->SetBranchAddress("LocalTimBr",&tempofftime2);  offbranch2->GetEvent(io); // yes it should be LocalTim, typo in root file
                  if((tempofftime2->Convert() < timeStart->Convert()) && (timeStart->Convert() - tempofftime2->Convert() < timeStart->Convert() - beforeRun2->Convert())){ 
                      beforeRun2->Set(tempofftime2->GetYear(), tempofftime2->GetMonth(), tempofftime2->GetDay(), tempofftime2->GetHour(), tempofftime2->GetMinute(), tempofftime2->GetSecond());
                      offbranch2->SetBranchAddress("AveAdcBr",&beforeADC2);     offbranch2->GetEvent(io);
                  }
                  if((tempofftime2->Convert() > timeEnd->Convert()) && (tempofftime2->Convert() - timeEnd->Convert() < afterRun2->Convert() - timeEnd->Convert())){ 
                      afterRun2->Set(tempofftime2->GetYear(), tempofftime2->GetMonth(), tempofftime2->GetDay(), tempofftime2->GetHour(), tempofftime2->GetMinute(), tempofftime2->GetSecond());
                      offbranch2->SetBranchAddress("AveAdcBr",&afterADC2);      offbranch2->GetEvent(io);
                  }
              }
              cout << "check off times   before: " << beforeRun1->AsSQLString() << "   after: " << afterRun1->AsSQLString() << "    2nd tel    before: " << beforeRun2->AsSQLString() << "   after: " << afterRun2->AsSQLString() << endl;
              cout << "ADCs    before: " << beforeADC1 << "   after: " << afterADC1 << "     2nd tel    before: " << beforeADC2 << "   after: " << afterADC2 << endl;
              
              // check - put out a warning if the start or end times are more than 5 mins different
              if(abs(int(beforeRun1->Convert()) - int(beforeRun2->Convert())) > 300){ cout << endl << "big discrepancy between telescopes in off runs before! check for missing files" << endl; } // must explicitly type ints or abs() gets confused
              if(abs(int(afterRun1->Convert()) - int(afterRun2->Convert())) > 300){ cout << endl << "big discrepancy between telescopes in off runs after! check for missing files" << endl; }
              
              // get avg change in ADC over run for subtraction
              dOffadc1 = (afterRun1 - beforeRun1)/double(cfOrig->GetYaxis()->GetNbins());
              dOffadc2 = (afterRun2 - beforeRun2)/double(cfOrig->GetYaxis()->GetNbins());
          }

	  // ----------------------------------------------------------- normalization -------------------------------------------------------------
	  
	  // get ADCs and normalize correlation function by their product (loop through each frame to normalize)
	  TH2D* ADC1 = new TH2D;   zippedFile->GetDirectory("Singles")->GetObject(Form("singlesT%s", tee1.Data()), ADC1);
	  TH2D* ADC2 = new TH2D;   zippedFile->GetDirectory("Singles")->GetObject(Form("singlesT%s", tee2.Data()), ADC2);
	  TProfile2D* cfNorm = new TProfile2D("cfNorm","normalized correlation function %s;relative time (ns);time in run (s)",cfOrig->GetXaxis()->GetNbins(),cfOrig->GetXaxis()->GetXmin(),cfOrig->GetXaxis()->GetXmax(),cfOrig->GetYaxis()->GetNbins(),cfOrig->GetYaxis()->GetXmin(),cfOrig->GetYaxis()->GetXmax());
	  
	  // divide correlation function by product of ADCs at each frame
	  if(off == 0){
	      for (int iy=1; iy<=cfOrig->GetYaxis()->GetNbins(); iy++){
	          TH1D* adc1y = ADC1->ProjectionX("first ADC proj",iy,iy);     double t1mean = adc1y->GetMean(); // think of how to check for last bin overflow and adjust mean!!!!
	          TH1D* adc2y = ADC2->ProjectionX("second ADC proj",iy,iy);    double t2mean = adc2y->GetMean();
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
	          TH1D* adc1y = ADC1->ProjectionX("first ADC proj",iy,iy);     double t1mean = adc1y->GetMean();
	          TH1D* adc2y = ADC2->ProjectionX("second ADC proj",iy,iy);    double t2mean = adc2y->GetMean();
	          for (int ix=1; ix<=cfOrig->GetXaxis()->GetNbins(); ix++){
	              cfNorm->Fill(cfOrig->GetXaxis()->GetBinCenter(ix), cfOrig->GetYaxis()->GetBinCenter(iy), (cfOrig->GetBinContent(ix,iy) - (t1mean*offadc2) - (t2mean*offadc1) + (offadc1*offadc2))/((t1mean-offadc1)*(t2mean-offadc2)));
	          }
	          offadc1 += dOffadc1;   offadc2 += dOffadc2;
	          delete adc1y, adc2y;
	      }
	  }
	  cfNorm->Write("normalizedCF");
	  
	  // -------------------------------------------------------------- noise removal -----------------------------------------------------------

	  // rebin for better noise removal
	  cfNorm->RebinY(16); // ~2s frames -> ~32s frames

	  // FFTs and noise removal
	  cfNorm->Draw("COLZ");
	  double fft = DoFFT(cfNorm, 1);
	  cfNorm = NoiseRemoveExactFreq(cfNorm, 0.4991, 1);
	  cfNorm = NoiseRemoveExactFreq(cfNorm, 0.0628, 2);
	  cfNorm = NoiseRemoveExactFreq(cfNorm, 0.6754, 3);
	  fft = DoFFT(cfNorm, 2); 

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
	          double xshift = cfNorm->GetXaxis()->GetBinCenter(ix) - opd->GetPointX(iy);
	          cfShift->Fill(xshift, cfNorm->GetYaxis()->GetBinCenter(iy), cfNorm->GetBinContent(ix,iy));
	      }
	    }
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
    if (times == 1){
        CFy->Fit("pol0", "S0Q");
        height = CFy->GetFunction("pol0")->GetParameter(0);
    }
    
    double AmpCos(0), CosNorm(0), AmpSin(0), SinNorm(0);
    for(int ix = 1; ix <= CFy->GetXaxis()->GetNbins(); ++ix){
      if(CFy->GetBinCenter(ix) < -100 || CFy->GetBinCenter(ix) > 100){ 
        AmpCos += cos(FFTFreq*CFy->GetBinCenter(ix))*(CFy->GetBinContent(ix)-height);
        CosNorm += pow(cos(FFTFreq*CFy->GetBinCenter(ix)), 2);
        AmpSin += sin(FFTFreq*CFy->GetBinCenter(ix))*(CFy->GetBinContent(ix)-height);
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
      CFMinusNoise->Fill(CF->GetXaxis()->GetBinCenter(ix),CF->GetYaxis()->GetBinCenter(TimeSlice), CF->GetBinContent(ix, TimeSlice) - (height+Amp*cos(FFTFreq*CF->GetXaxis()->GetBinCenter(ix)-Phi)));
    }
  }
  
  //Save Parameter Graphs to a canvas and write to output file 
  /*  TCanvas * AllParams = new TCanvas(Form("NoiseParameters%d", times), "AllParams", 1600, 800);
  AllParams->Divide(2,2);
  AllParams->cd(1); HeightGraph->Draw("A*");
  AllParams->cd(2); AmpGraph->Draw("A*");
  AllParams->cd(3); FreqGraph->Draw("A*");
  AllParams->cd(4); PhaseGraph->Draw("A*");
  AllParams->Write(Form("ParameterGraphs%d", times));
  */
  
  //Write and return Correlation Functions without noise 
  CFMinusNoise->SetTitle("Normalized Correlation Function Without Noise; Relative Time (ns); Time");
  CFMinusNoise->Write(Form("CFMinusNoise%d", times));
  
  return CFMinusNoise;
}

// ===========================================================================================================================

// function to make TNtuple from python calculations stored in txt files

TNtuple* GeometryNtuple(TString keyName, int N){
  //TNtuple* geom = new TNtuple("PairGeometryInfo",Form("%s geometry",keyName.Data()),"time:baseline:OPD:u:v:w");
  TNtuple* geom = new TNtuple("PairGeometryInfo",Form("%s geometry",keyName.Data()),"frame:time:u:v:baseline:OPD");
  float vals[6], junk;
  //TString junk; // trying w junk as a string
  
  if (keyName == "T0T1"){                    // T0T1 is a special case
    for (int i=0; i<6; i++) vals[i]=0.0;
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
  pyinfoFile.open(Form("delays/pyinfo%s.txt", nameBase.Data()));       if (!pyinfoFile.is_open()){cout << "Problem opening python file!!!\n";  return 0;}
  //if (pyinfoFile.is_open()){ cout << "file opened ok! " << nameBase.Data() << endl;}

  pyinfoFile >> junk >> junk >> junk >> junk >> junk >> junk; 
  for (int i=0; i<N; i++){
    pyinfoFile >> vals[0] >> vals[1] >> vals[2] >> vals[3] >> vals[4] >> vals[5]; // frame, time, u, v, baseline, opd
    //if(i < 10){ cout << vals[0] << " " << vals[1] << " " << vals[2] << " " << vals[3] << " " << vals[4] << " " << vals[5] << endl;}
    geom->Fill(vals);
  }
  pyinfoFile.close();
  return geom;
}
  
  // below is the part dr lisa wrote
  /*
  ifstream delayFile, baselineFile, baselineCoordFile;
  delayFile.open(Form("delays/%s.delay",nameBase.Data()));                   if (!delayFile.is_open()){cout << "Problem opening delay file!!!\n";  return 0;}
  baselineFile.open(Form("delays/%s.baseline",nameBase.Data()));             if (!baselineFile.is_open()){cout << "Problem opening baseline file!!!\n";  return 0;}
  baselineCoordFile.open(Form("delays/%s.baselinecoord",nameBase.Data()));   if (!baselineCoordFile.is_open()){cout << "Problem opening baseline coordinate file!!!\n";  return 0;}
  for (int i=0; i<N; i++){
    baselineFile      >> vals[0] >> vals[1];                        // time, baseline
    delayFile         >> junk    >> vals[2];                        // OPD
    baselineCoordFile >> junk    >> vals[3] >> vals[4] >> vals[5];  // u, v, w
    geom->Fill(vals);
  }
  delayFile.close();
  baselineFile.close();
  baselineCoordFile.close();
  return geom;
  */
