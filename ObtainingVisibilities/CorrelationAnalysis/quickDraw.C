// josie 13 nov 2023
// abbreviated version of the full analysis to draw objects and print to pdf, check data quickly and easily without saving to root file
// run with: root quickDraw.C\(\"file\ name_ZippedFrames.root\"\)
// must run with python in delays directory

double DoFFT (TProfile2D * CF, int times);
TProfile2D * NoiseRemoveExactFreq(TProfile2D *CF, double FFTFreq, int times);

void quickDraw(TString filename){
    
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
  TString py = "python3 delays/CalcMvtUpdated.py"; // <-- edit path if python script is located somewhere else!
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
  TString pairsList[npairs];    int iplist(0);
 
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
    
  // store length of longest run, earliest start time (local) as strings to send to python    
  length = to_string(longestRun);   //cout << "to str longest run " << to_string(longestRun) << "  longest run " << longestRun << endl;
  loctime = ltime->AsSQLString();   //cout << "ltime num " << ltime << "  ltime string " << loctime << endl; 
  nframes = to_string(mostFrames);  frameSize = to_string(frameWidth);      //cout << "most frames " << nframes << "   width " << frameSize << endl; 
  
  cout << "checking python parameters " << length << "  " << nframes << "  " << frameSize << "  " << *(source) << "  " << loctime << endl << endl; 
  
  // run python script to calculate delays/baselines via bash from within this macro 
  cout << "================> Running python script to calculate delays <=====================" << endl;
  //runpy = py + " " + length + " " + *(source) + " " + loctime; // note: source is a pointer so must be dereferenced to add to other strings
  runpy = py + " " + length + " " + nframes + " " + frameSize + " " + *(source) + " " + loctime;
  gSystem->Exec(runpy);
  cout << "===============================> Python is done! <================================" << endl;

  TCanvas* c1 = new TCanvas;
  c1->Divide(2,2);
  //TString pdfName = filename.ReplaceAll(".root",".pdf");
  //cout << "output pdf is named: " << pdfName.Data() << endl;
  c1->Print("checkfiles.pdf[");
  //c1->SaveAs(Form("%s(",pdfName.Data()));
  
  // ===================================================== start of the main analysis ===========================================================
  
  // loop over each directory we find in the root file (taken from AnalyzeVersiiFiles.C)
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
	  cfOrig->SetTitle("raw correlation function;relative time (ns); frames (s)");
	  c1->cd(1);	  cfOrig->Draw("COLZ");
	  
	  // get ADCs and normalize correlation function by their product (loop through each frame to normalize)
	  TH2D* ADC1 = new TH2D;   zippedFile->GetDirectory("Singles")->GetObject(Form("singlesT%s", tee1.Data()), ADC1);
	  TH2D* ADC2 = new TH2D;   zippedFile->GetDirectory("Singles")->GetObject(Form("singlesT%s", tee2.Data()), ADC2);
	  TProfile2D* cfNorm = new TProfile2D("cfNorm","normalized correlation function;relative time (ns);time in run (s)",cfOrig->GetXaxis()->GetNbins(),cfOrig->GetXaxis()->GetXmin(),
					      cfOrig->GetXaxis()->GetXmax(),cfOrig->GetYaxis()->GetNbins(),cfOrig->GetYaxis()->GetXmin(),cfOrig->GetYaxis()->GetXmax());
	  for (int iy=1; iy<=cfOrig->GetYaxis()->GetNbins(); iy++){
	    TH1D* adc1y = ADC1->ProjectionX("first ADC proj",iy,iy);     double t1mean = adc1y->GetMean(); // think of how to check for last bin overflow and adjust mean!!!!
	    TH1D* adc2y = ADC2->ProjectionX("second ADC proj",iy,iy);    double t2mean = adc2y->GetMean();

	    for (int ix=1; ix<=cfOrig->GetXaxis()->GetNbins(); ix++){
              cfNorm->Fill(cfOrig->GetXaxis()->GetBinCenter(ix), cfOrig->GetYaxis()->GetBinCenter(iy), cfOrig->GetBinContent(ix,iy)/(t1mean*t2mean));
	    }
	    delete adc1y, adc2y;
	  }
	  c1->cd(2);	  cfNorm->Draw("COLZ");

	  // rebin for better noise removal
	  int rebinby = cfNorm->GetYaxis()->GetNbins()/300; // may need to adjust this more 
	  cfNorm->RebinY(rebinby);

	  // FFTs and noise removal
	  cfNorm->Draw("COLZ");
	  double fft = DoFFT(cfNorm, 1);
	  cfNorm = NoiseRemoveExactFreq(cfNorm, 0.4991, 1);
	  cfNorm = NoiseRemoveExactFreq(cfNorm, 0.0628, 2);
	  cfNorm = NoiseRemoveExactFreq(cfNorm, 0.6754, 3);
	  fft = DoFFT(cfNorm, 2); 

	  // draw the heatmap
	  cfNorm->SetTitle("heatmap;relative time (ns);time in run (s)");
	  c1->cd(3);   cfNorm->Draw("COLZ");

	  // read in delays and baselines 
	  double cableDelays[5]    = {676.8, 676.8, 585.0, 955.0, 1063.7};
	  double refDelay;   refDelay = -1*(nbo*4.0 + cableDelays[t2] - cableDelays[t1]); // seg fault if not defined separately
	  TGraph* opd;                     TGraph* baseline;
	  TGraph* opdGraph = new TGraph;   TGraph* baselineGraph = new TGraph;
	  TSpline3* opdSpline;             TSpline3* baseSpline;
	  TProfile* cfFlat = new TProfile;

	  // skip opd shift, baselines, etc for T0T1 
	  if(t1 == 0 && t2 == 1){
	    //cfFlat = cfNorm->ProfileX("",1,-1); // this is giving the early exit error and I can't figure out why!! 
	  }
	  else {
	    // **** for now - reading in opd and baseline the old way so I have a working version, but will update when I do the ttree ****
	    if(t1 == 0){
	      opd = new TGraph(Form("delays/T1T%s.delay", tee2.Data()));
	      baseline = new TGraph(Form("delays/T1T%s.baseline", tee2.Data()));
	    }
	    else{
	      opd = new TGraph(Form("delays/%s.delay", keyName.Data()));
	      baseline = new TGraph(Form("delays/%s.baseline", keyName.Data()));
	    }
	    for(int r=0; r<opd->GetN(); r++){ opd->SetPointY(r,-opd->GetPointY(r));}
	    for(int r=0; r<baseline->GetN(); r++){ baseline->SetPointY(r,-baseline->GetPointY(r));}

	    // turn into splines and evaluate at each knot to make new TGraph *** this may be eliminated now that python calculates opd and baseline for each frame, but will keep for now
	    opdSpline = new TSpline3("opdSpline",opd);
	    baseSpline = new TSpline3("baseSpline",baseline);
	    for(int y=1; y<=cfNorm->GetYaxis()->GetNbins(); y++){
	      double evalOPD = opdSpline->Eval(opd->GetPointX(0) + cfNorm->GetYaxis()->GetBinCenter(y)) - refDelay;
	      opdGraph->AddPoint(evalOPD, cfNorm->GetYaxis()->GetBinCenter(y));
	      double evalBase = baseSpline->Eval(baseline->GetPointX(0) + cfNorm->GetYaxis()->GetBinCenter(y));
	      baselineGraph->AddPoint(cfNorm->GetYaxis()->GetBinCenter(y), evalBase);
	    }

	    // make new histogram to hold opd shifted cf and fill by opd shift
	    TProfile2D* cfShift = new TProfile2D("cfShift","OPD shifted correlation function;relative time (ns); time in run (s)", cfNorm->GetXaxis()->GetNbins()*4.0, cfNorm->GetXaxis()->GetXmin(),
					       cfNorm->GetXaxis()->GetXmax(), cfNorm->GetYaxis()->GetNbins(), cfNorm->GetYaxis()->GetXmin(), cfNorm->GetYaxis()->GetXmax());
	    for (int iy=1; iy<=cfNorm->GetYaxis()->GetNbins(); iy++){
	      for (int ix=1; ix<=cfNorm->GetXaxis()->GetNbins(); ix++){
		double xshift = cfNorm->GetXaxis()->GetBinCenter(ix) - opdGraph->GetPointX(iy);
		cfShift->Fill(xshift, cfNorm->GetYaxis()->GetBinCenter(iy), cfNorm->GetBinContent(ix,iy));
	      }
	    }
	    // project shifted correlation function
	    cfFlat = cfShift->ProfileX("",1,-1);
	  }

	  cfFlat->SetTitle("projected correlation function;relative time(ns);g^2");
	  c1->cd(4); cfFlat->Draw();
	  c1->Print("checkfiles.pdf");
  
      }// end of loop over pairs
  }// end of loop over keys

  c1->Print("checkfiles.pdf]");     
      
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
  //FFTCombine->Write(Form("FFTProjection%d",times));
  
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
  //CFMinusNoise->Write(Form("CFMinusNoise%d", times));
  
  return CFMinusNoise;
}
