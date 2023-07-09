/*

Mackenzie Scott and Michael Lisa -- Begining of Summer 2022

The aim of this program is to analyze correlations created by both the FPGAs and C++ Software. 

Compile/Execute as a root Macro with: root CorrealtionAnalysis.C

It must be able to properly read in ZippedFrames.root as well as the optical path delay

*/


// Nolan does not follow a consistent convention for the telescope
// Looking at data y2021m02d23h20m44s45, in his header information, he identifies TEL_X and TEL_Y
// In T1T2, he has TEL_X=2, TEL_Y=1
// In T1T3, he has TEL_X=3, TEL_Y=1
// In T1T4, he has TEL_X=4, TEL_Y=1
// In T2T3, he has TEL_X=3, TEL_Y=2
// In T2T4, he has TEL_X=4, TEL_Y=2
// In T3T4, he has TEL_X=3, TEL_Y=4  <-------- why did you do this different?????
// ..... and for all I know, he might switch them up at random in other runs, just for kicks.
// /sigh.  Okay, well then I have to account for it.


// ADD Is Not Properly implemented for TProfile2D. 
// One MUST define a new histogram which is sad for the error bars. 
// One can also not SET the error bars for a Tprofile2D which would remedy the sitution 


//THIS IS THE ONLY THING YOU NEED TO CHANGE
//If you want to include OFF Data make true
//If runs are shorter than 30 minutes then the whole section zero might need to be edited
bool offData = true;

//Prototypes 
double DoFFT (TProfile2D * CF, int times);
TProfile2D * NoiseRemoveExactFreq(TProfile2D *CF, double FFTFreq, int times);

//TProfile2D * RemoveOne(TProfile2D * CF);
//TProfile2D * NoiseRemove(TProfile2D *CF, double FFTFreq, int times);

void CorrelationAnalysis(){
    
  //Set Graph Style
  gStyle->SetOptStat(0);
  
  //Values for the headers 
  int tx,ty,version,corrDelay;  // The two telescopes, version and Correlation Delay
  int NumberOfFrames(0);        // Number of Frames. How long is the run?
  float dtWindow;               // one second?  one-eighth of a second?
  float dtBucket;               // so far it's always 4 ns.
  TDatime* whenTaken(0);        // setting pointer to zero before using it in SetBranchAddress is completely necessary!
  TString * source;             // What star?
  
  //Define Objects that will read in Objects from Root File 
  TFile *CorrelationFile = new TFile("ZippedFrames.root","READONLY");
  TTree *header;
  TH1D *ADC1N(0), *ADC2N(0);
  TH2D *ADC1(0),  *ADC2(0);
  TProfile2D *CF(0);
  
  //Read in Header Data
  CorrelationFile->GetObject("Header",header);
  header->SetBranchAddress("VersionBr",&version);           // Version of code: Positive Incrament for Nolan, Negative Incrament for OSU
  header->SetBranchAddress("SourceBr",&source);             // String containing star name 
  header->SetBranchAddress("T1Br",&tx);                     // First Telescope
  header->SetBranchAddress("T2Br",&ty);                     // Second Telescope
  header->SetBranchAddress("corrDelayBr",&corrDelay);       // Correlator Delay used whe making the correlations 
  header->SetBranchAddress("tWindowSecBr",&dtWindow);       // Length of the "Frames"
  header->SetBranchAddress("LocalTimeBr",&whenTaken);       // When the run began
  header->SetBranchAddress("samplingNsBr",&dtBucket);       // Sampling rate 4ns
  header->SetBranchAddress("NumWindowsBr",&NumberOfFrames); // Number of Frames
  header->GetEntry(0);
  
  //Clone Tree to Save in Analysis File
  TTree* HeaderTreeCopy = header->CloneTree();
  HeaderTreeCopy->SetDirectory(0);
  
  //Output all the run information
  cout << "This correlation between T" << tx << " and T" << ty 
       << " on source " << source->Data() 
       << " taken " << whenTaken->AsString() 
       << " with " << NumberOfFrames << " frames of " << dtWindow << " seconds duration\n\n";
  
  //Now we have this "T0" to worry about.... It is located at T1
  if (0==tx) tx=1; if (0==ty) ty=1;
  
  //Read in CF
  CorrelationFile->GetObject("CF12",CF);
  
  //Create a Root File for Writing
  TFile* outputFile = new TFile("Analysis.root","RECREATE");
  HeaderTreeCopy->Write("Header");
  CF->SetTitle("Correlation Function; Relative Time (ns); Time into Observation Run (sec)");
  CF->Write("CFOriginal");
  
  //----------------------------------------------------------------------------
  //------------------- Zero: If OFF data. Read it in --------------------------
  
  //Initialize Background ADCs
  double BackADC1(0), BackADC2(0), dBackADC1(0), dBackADC2(0);
  
  if(offData == true){
    //Convert Times into just seconds so we can compare them
    int BeginingTime = whenTaken->Convert();
    int EndTime  = BeginingTime + (int)(NumberOfFrames *dtWindow);

    //Read in TextFile with all ADC values 
    std::ifstream myfile;
    myfile.open("../../../ADCValues/ADCValues.txt",std::ifstream::in);

    //Use these values to match lines of text file with this run
    double BackADC1Begin(0), BackADC2Begin(0), BackADC1End(0), BackADC2End(0), ADCCheck(0);
    int yearCheck(0), monthCheck(0), dayCheck(0), hourCheck(0), minuteCheck(0), teleCheck(0); 
    string line;
    
    //Read in the ADC values from the file and match them to this file
    while (getline (myfile,line)){
      
      //Read in data from ADC File 
      myfile >> yearCheck >> monthCheck >> dayCheck >> hourCheck >> minuteCheck >> teleCheck >> ADCCheck;
    
      //OSU saves time slighlty differently. Convert it
      if (version < 0){
          if(monthCheck == 3 && dayCheck == 13 && hourCheck == 4){hourCheck -=1;} //Daylight Savings Time
          hourCheck +=7;
          if (hourCheck > 24){hourCheck -= 24; ++dayCheck;}
      }
      
      //Save temporary value from ADC file as TDatime
      TDatime * temporary = new TDatime(yearCheck, monthCheck, dayCheck, hourCheck, minuteCheck, 0);
      
      //Compare values and match but only if different by less than 600 seconds
      if(abs(BeginingTime-(int)temporary->Convert()) < 750 && tx == teleCheck){BackADC1Begin = ADCCheck;} 
      if(abs(BeginingTime-(int)temporary->Convert()) < 750 && ty == teleCheck){BackADC2Begin = ADCCheck;} 
      if(abs((int)temporary->Convert()-EndTime)      < 750 && tx == teleCheck){BackADC1End = ADCCheck;} 
      if(abs((int)temporary->Convert()-EndTime)      < 750 && ty == teleCheck){BackADC2End = ADCCheck;} 
    }
    
    //Output to terminal so it can be double checked if they fault 
    cout << "Beg Background ADC Values: " << BackADC1Begin << " and " << BackADC2Begin << endl;
    cout << "End Backgroung ADC Values: " << BackADC1End << " and " << BackADC2End << endl << endl;
    
    //Because it does fault. Input by hand when it does. Doesn't work as part of a batch job.
    //if(BackADC1Begin == 0){cout << "Help I couldn't find BackADC1Begin. What is it? "; cin >> BackADC1Begin;}
    //if(BackADC2Begin == 0){cout << "Help I couldn't find BackADC2Begin. What is it? "; cin >> BackADC2Begin;}
    //if(BackADC1End == 0)  {cout << "Help I couldn't find BackADC1End. What is it? "; cin >> BackADC1End;}
    //if(BackADC2End == 0)  {cout << "Help I couldn't find BackADC2End. What is it? "; cin >> BackADC2End;}
    
    //Where the background ADC starts
    BackADC1 = BackADC1Begin;
    BackADC2 = BackADC2Begin;
    
    //How much the background ADC changes from frame to frame 
    dBackADC1 = (BackADC1End - BackADC1Begin)/CF->GetYaxis()->GetNbins();
    dBackADC2 = (BackADC2End - BackADC2Begin)/CF->GetYaxis()->GetNbins();
  }
  
  //Create CFTemp to hold new Correlations with the backgrounds removed
  TProfile2D* CFTemp = (TProfile2D*)CF->Clone("TempForOFF");  // I "clone" so that it has the same binning and axis limits
  CFTemp->Reset();
  
  //----------------------------------------------------------------------------
  //--------------- First: Normalize the Correlation Function ------------------

  //Create a Denominator to hold the Normalization factors for the CFs
  TProfile2D* denom = (TProfile2D*)CF->Clone("denominator");  
  denom->Reset();
    
  //Fill Denominator for Nolan's frames. Positive version means nolan  
  if(version > 0){
    //Read in Objects. Check they are nonzero. Save them to the out file
    CorrelationFile->GetObject("singles1", ADC1N);
    CorrelationFile->GetObject("singles2", ADC2N);
    if ((CF==0)||(ADC1N==0)||(ADC2N==0)) {cout << "Hey missing histogram! I quit!\n"; return;}
    ADC1N->Write("ADC1"); ADC2N->Write("ADC2");
    
    //Fill Numerator and Denominator
    cout << "Normalization \n\n";
    for (int TimeSlice=1; TimeSlice<= CF->GetYaxis()->GetNbins(); TimeSlice++){
      for (int ix=1; ix<=denom->GetXaxis()->GetNbins(); ix++){
        denom->Fill(CF->GetXaxis()->GetBinCenter(ix), CF->GetYaxis()->GetBinCenter(TimeSlice),(ADC1N->GetBinContent(TimeSlice)-BackADC1)*(ADC2N->GetBinContent(TimeSlice)-BackADC2));
        CFTemp->Fill(CFTemp->GetXaxis()->GetBinCenter(ix), CFTemp->GetYaxis()->GetBinCenter(TimeSlice), (CF->GetBinContent(ix,TimeSlice) - (ADC1N->GetBinContent(TimeSlice)*BackADC2) - (ADC2N->GetBinContent(TimeSlice)*BackADC1) + (BackADC1*BackADC2)));
      }
      BackADC1 += dBackADC1;
      BackADC2 += dBackADC2;
    }
  }

  //Fill Denominator OSU frames. Negative Version Means OSU
  if(version < 0){ 
    //Read in Objects. Check they are nonzero. Save them to the out filet
    CorrelationFile->GetObject("singles1",ADC1);
    CorrelationFile->GetObject("singles2",ADC2);
    if ((CF==0)||(ADC1==0)||(ADC2==0)) {cout << "Hey missing histogram! I quit!\n"; return;}
    if (ADC1->GetYaxis()->GetNbins() != CF->GetYaxis()->GetNbins()) {cout << "Illogical!\n"; return;}
    ADC1->Write("ADC1"); ADC2->Write("ADC2");
    
    //Fill Numerator and Denominator
    cout << "Normalization \n\n";
    for (int TimeSlice=1; TimeSlice<= CF->GetYaxis()->GetNbins(); TimeSlice++){
      TH1D* ADC1y = ADC1->ProjectionX("temp1",TimeSlice,TimeSlice);
      TH1D* ADC2y = ADC2->ProjectionX("temp2",TimeSlice,TimeSlice);
      double ADC1Mean(ADC1y->GetMean()), ADC2Mean(ADC2y->GetMean());
      for (int ix=1; ix<=CF->GetXaxis()->GetNbins(); ix++){
        denom->Fill(CF->GetXaxis()->GetBinCenter(ix),CF->GetYaxis()->GetBinCenter(TimeSlice), (ADC1Mean-BackADC1) * (ADC2Mean-BackADC2)); 
        CFTemp->Fill(CFTemp->GetXaxis()->GetBinCenter(ix), CFTemp->GetYaxis()->GetBinCenter(TimeSlice), (CF->GetBinContent(ix,TimeSlice) - (ADC1Mean*BackADC2) - (ADC2Mean*BackADC1) + (BackADC1*BackADC2)));
      }
      BackADC1 += dBackADC1;
      BackADC2 += dBackADC2;
      delete ADC1y, ADC2y;
    }
  }  
  
  //Reset CF and Save the new numerator
  CF->Reset();
  CF = (TProfile2D*)CFTemp->Clone("CFNew");
  
  //Normalize!!
  CF->Divide(denom);
  CF->SetTitle("Normalized Correlation Function; Relative Time (ns); Time into Observation Run (sec)");
  CF->Write("CFNormalized");

  //----------------------------------------------------------------------------
  //-------- Second: Take FFT and Remove Noise or Normalize to Zero ------------
  
  //Make the Frames Course so we can fit the noise better
  CF->RebinY(256);
  
  //Do FFTs and Remove Noise 
  double FFTFreq = DoFFT(CF, 1);  //Put How many times you call it at the end. Start with 1.
  CF = NoiseRemoveExactFreq(CF, 0.4991, 1); //80 MHz //Put how many times you call it at the end. Start with 1. 
  CF = NoiseRemoveExactFreq(CF, 0.0628, 2); //10 MHz
  CF = NoiseRemoveExactFreq(CF, 0.6754, 3); //107MHz
  FFTFreq = DoFFT(CF, 2);
  
  //---------------> or
  
  //If you want to simply remove 1 from the correlation rather than remove noise
  //CF = RemoveOne(CF);
  
  //----------------------------------------------------------------------------
  //----------------------- Third: Heat Map ------------------------------------
  
  //Output heatmap just for fun
  cout << "Generating HeatMap \n\n";
  TProfile2D* HeatMap = (TProfile2D*)CF->Clone("HeatMap");
  HeatMap->RebinY(2);
  HeatMap->SetTitle("Heat Map; Relative Time (ns); Time into Observation Run (sec)");
  HeatMap->Write("HeatMap"); 
  
  // ---------------------------------------------------------------------------
  //----------------- Finally: Now Calculate Delays ----------------------------
  
  //Create Finely Binned Graph 
  TProfile * FineShiftShort; 
  TProfile2D * OnlyShift;
  
  //Save Baseline and OPDs to Graphs
  TGraph* delays, *baseline;
  TGraph * OPDGraph  = new TGraph(); OPDGraph->SetLineWidth(4);
  TGraph * BaselineGraph = new TGraph(); 
  TSpline3 * DelaySpline, *BaselineSpline; 
  
  //Find the reference delay
  double CableDelays[5]    = {676.8, 676.8, 585.0, 955.0, 1063.7};
  double ReferenceDelay = (corrDelay*dtBucket + CableDelays[tx] - CableDelays[ty])*-1; 
  
  cout << "ReferenceDelay" << ReferenceDelay << endl;
  
  //If pair is beam split no delays required
  if(tx==ty){goto cont;}

  cout << "Now Time Delay Shift \n\n";
  
  //Read in delays and baselines
  if (tx>ty){
    delays = new TGraph(Form("../delays/T%dT%d.delay",ty,tx));
    baseline = new TGraph(Form("../delays/T%dT%d.baseline",ty,tx));
    //Flip if telescopes are in opposite order 
    for (int kkk=0; kkk < delays->GetN(); kkk++){delays->SetPointY(kkk,-delays->GetPointY(kkk));} 
    for (int kkk=0; kkk < baseline->GetN(); kkk++){baseline->SetPointY(kkk,-baseline->GetPointY(kkk));} 
  } // note order of tx,ty
 
  else {
    delays = new TGraph(Form("../delays/T%dT%d.delay",tx,ty));
    baseline = new TGraph(Form("../delays/T%dT%d.baseline",tx,ty));
  } // note order of tx,ty
  
  //Make those points into splines 
  DelaySpline = new TSpline3("DelaySpline",delays);
  BaselineSpline = new TSpline3("BaselineSpline",baseline);
  
  //Get Knots for Delays 
  for(int TimeSlice = 1; TimeSlice<=CF->GetYaxis()->GetNbins(); TimeSlice++){
    double OPDelay = DelaySpline->Eval(delays->GetPointX(0) + CF->GetYaxis()->GetBinCenter(TimeSlice)) - ReferenceDelay;
    double BaselineEval  = BaselineSpline->Eval(baseline->GetPointX(0) + CF->GetYaxis()->GetBinCenter(TimeSlice)); 
    OPDGraph->AddPoint(OPDelay, CF->GetYaxis()->GetBinCenter(TimeSlice));
    BaselineGraph->AddPoint(CF->GetYaxis()->GetBinCenter(TimeSlice), BaselineEval);
  }
  
  //Write Basline and OPD to output file
  BaselineGraph->Write("BaselineGraph");
  OPDGraph->Write("OPDGraph");
  
  //Define new graphs with their style
  gStyle->SetOptFit(1);
  FineShiftShort = new TProfile("Shift", "Time Delay Shifted Peak", 2*CF->GetXaxis()->GetNbins(), CF->GetXaxis()->GetXmin()/2, CF->GetXaxis()->GetXmax()/2, "");
  OnlyShift      = new TProfile2D("shifting", "Time Delay Shifted non-Projected", 4*CF->GetXaxis()->GetNbins(), CF->GetXaxis()->GetXmin(), CF->GetXaxis()->GetXmax(), CF->GetYaxis()->GetNbins(), CF->GetYaxis()->GetXmin(),CF->GetYaxis()->GetXmax());
  
  //Fill the new Graphs 
  for (int frame=1; frame<=CF->GetYaxis()->GetNbins(); frame++){
    for (int itau=1; itau<=CF->GetXaxis()->GetNbins(); itau++){
      double shifted_time = CF->GetXaxis()->GetBinCenter(itau) - OPDGraph->GetPointX(frame);
      FineShiftShort->Fill(shifted_time,CF->GetBinContent(itau,frame));
      OnlyShift->Fill(shifted_time, CF->GetYaxis()->GetBinCenter(frame), CF->GetBinContent(itau,frame));
    }
  }
  
  //Write Full TProfile2D to File 
  OnlyShift->Write("JustShift");
  
  //Skip here if pair is the same telescope
  cont:
  if(tx == ty){
    gStyle->SetOptFit(1);
    FineShiftShort = CF->ProfileX("HBTIGuess", 1, CF->GetYaxis()->GetNbins());
  }
  
  //Define Fit Parameters/Windows
  double relTimePar(-5), relTimeWindow(20);
  if((tx==1 && ty==2)){relTimePar = -20;} if((tx==2 && ty==1)){relTimePar =  20;}
  if((tx==1 && ty==3)){relTimePar = -10;} if((tx==3 && ty==1)){relTimePar =  10;}
  if((tx==1 && ty==4)){relTimePar =  2;}  if((tx==4 && ty==1)){relTimePar = -2;}
  if((tx==2 && ty==3)){relTimePar =  10;} if((tx==3 && ty==2)){relTimePar = -10;}
  if((tx==2 && ty==4)){relTimePar =  18;} if((tx==4 && ty==2)){relTimePar = -18;}
  if((tx==3 && ty==4)){relTimePar =  8;}  if((tx==4 && ty==3)){relTimePar = -8;}
  double sigMin(0), sigMax(7);
  relTimePar = -5;
  
  
  TCanvas * Peak = new TCanvas("c1", "c1", 0, 0, 1600, 1200);
  TF1* PeakFit = new TF1("HBTpeak","([0]*exp(-pow((x-[1])/[2],2)/2.0))/([2]*sqrt(2*pi))",-300,300); 
  PeakFit->SetParName(0,"Area");          PeakFit->SetParameter(0,0.0);     
  PeakFit->SetParName(1,"#tau_{0} (ns)"); PeakFit->SetParameter(1,relTimePar); PeakFit->SetParLimits(1,relTimePar-relTimeWindow, relTimePar+relTimeWindow); 
  PeakFit->SetParName(2,"#sigma (ns)");   PeakFit->SetParameter(2,4.5); PeakFit->SetParLimits(2, sigMin,  sigMax);
  FineShiftShort->Fit(PeakFit);
  FineShiftShort->SetTitle("; Relative Time (ns); g^{2}(#tau)-1");
  FineShiftShort->GetYaxis()->SetTitleOffset(1.05);
  

  
  //Write and Draw Peak With Fit 
  Peak->SaveAs("HBTPicture.pdf");
  FineShiftShort->SetTitle("Normalized and Time Shifted Correlation Function; Relative Time (ns); Magnitude");
  FineShiftShort->Write("HBTPeakWithFit");
  
  outputFile->Close();
  //gApplication->Terminate();

  
  //The End...Finally
  
}
 
//------------------------------- FFT ------------------------------------------
double DoFFT (TProfile2D * CF, int times){
  cout << "Do an FFT of the data\n\n";
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

//-------------------------- Remove Noise --------------------------------------
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





/*
//---------------------- Only Subtract One -----------------------------------
TProfile2D * RemoveOne(TProfile2D * CF){
  cout << "Subtract 1 From Correlations \n\n";
  
  TProfile2D* CFMinusOne = (TProfile2D*)CF->Clone("CFMinusOne");  // I "clone" so that it has the same binning and axis limits
  CFMinusOne ->Reset();
  
  for(int TimeSlice = 1; TimeSlice<=CF->GetYaxis()->GetNbins(); TimeSlice++){
    TH1D * CFy = CF->ProjectionX("temp", TimeSlice, TimeSlice);
    CFy->Fit("pol0", "S0QW");
    for (int ix=1; ix<=CF->GetXaxis()->GetNbins(); ix++){
      CFMinusOne->Fill(CFMinusOne->GetXaxis()->GetBinCenter(ix),CFMinusOne->GetYaxis()->GetBinCenter(TimeSlice), CF->GetBinContent(ix, TimeSlice) -CFy->GetFunction("pol0")->GetParameter(0)); 
    }
  }
  
  CFMinusOne->SetTitle("Correlation Function Minus DC Component; Relative Time (ns); Time into Observation Run (sec)");
  CFMinusOne->Write("CFMinusOne");
  
  return CFMinusOne;
}
*/


/*

//---------------------------- Remove Noise ------------------------------------
TProfile2D * NoiseRemove(TProfile2D *CF, double Freq, int times){
  cout << "Remove noise WITH the DC component\n\n";
  
  //Chi Square Check at Zero
  TF1 *NoiseAmpZero = new TF1("NR","[0]+[1]*cos([2]*x+[3])");
  if(times == 1){NoiseAmpZero->SetParameter(0,1.0);} else{NoiseAmpZero->SetParameter(0, 0.0);}
  NoiseAmpZero->FixParameter(1,0.0);
  NoiseAmpZero->SetParameter(2,0.499);    NoiseAmpZero->SetParLimits(2,Freq-0.001,Freq+0.001);
  NoiseAmpZero->SetParameter(3,0.0);      NoiseAmpZero->SetParLimits(3,-0.1,2.0*TMath::Pi()+0.1);
  
  //Fit to see if there is noise 
  TF1 *NoiseRemoval = new TF1("NR","[0]+[1]*cos([2]*x+[3])");
  if(times == 1){NoiseRemoval->SetParameter(0,1.0);} else{NoiseRemoval->SetParameter(0, 0.0);}
  NoiseRemoval->SetParameter(1,0.0);
  NoiseRemoval->SetParameter(2,0.499);    NoiseRemoval->SetParLimits(2,Freq-0.001,Freq+0.001);
  NoiseRemoval->SetParameter(3,0.0);      NoiseRemoval->SetParLimits(3,-0.1,2.0*TMath::Pi()+0.1);
  
  //Holds just the values of the noise
  TProfile2D* noise = (TProfile2D*)CF->Clone(Form("noise1%d", times));  
  noise->Reset();
  
  //Define graphs to save parameters to
  TGraph * HeightGraph = new TGraph(); HeightGraph->SetTitle("DC Component; Time (Seconds); Normalized Counts");
  TGraph * AmpGraph    = new TGraph(); AmpGraph->SetTitle("Amplitude of Noise; Time (Seconds); Amplitude");
  TGraph * FreqGraph   = new TGraph(); FreqGraph->SetTitle("Frequency of Noise; Time (Seconds); Frequency (radians per ns)");
  TGraph * PhaseGraph  = new TGraph(); PhaseGraph->SetTitle("Phase of Noise; Time(Seconds); Phase (Radians)");

  for(int TimeSlice = 1; TimeSlice<=CF->GetYaxis()->GetNbins() ; TimeSlice++){ 
    TH1D * CFy = CF->ProjectionX("temp", TimeSlice, TimeSlice);
    
    //Exclude parts from the fit 
    for(int i = 1; i < CFy->GetXaxis()->GetNbins();  ++i){CFy->SetBinError(i, 0.0);} //Set all Errors to Zero
    for(int i = CFy->FindBin(-100); i < CFy->FindBin(100); ++i){CFy->SetBinContent(i, 0.0);} //Set Peak Region to Zero 
    for(int i = CFy->FindBin(200); i < CFy->FindBin(CFy->GetXaxis()->GetNbins()); ++i){CFy->SetBinContent(i, 0.0);} //Set End bins to zero for glitched correlation smearing 
    
    //Fit Noise and Zero
    TFitResultPtr s = CFy->Fit(NoiseAmpZero, "S0QW"); int status2 = s;
    TFitResultPtr r = CFy->Fit(NoiseRemoval, "S0QW");int status = r;
    
    //Change things around to get postive amplitude and postive phase
    if ((status2!=0)||(NoiseAmpZero->GetParameter(1)<0.0)){
      NoiseAmpZero->SetParameter(1,-NoiseAmpZero->GetParameter(1));
      double temp = NoiseAmpZero->GetParameter(3);
      temp = (NoiseAmpZero->GetParameter(3)<TMath::Pi())?NoiseAmpZero->GetParameter(3)+TMath::Pi():NoiseAmpZero->GetParameter(3)-TMath::Pi();
      NoiseAmpZero->SetParameter(3,temp);
      r = CFy->Fit(NoiseAmpZero,"S0QW");
    }
    if (NoiseAmpZero->GetParameter(3)<0.0){
      NoiseAmpZero->SetParameter(3.0,NoiseAmpZero->GetParameter(3)+2.0*TMath::Pi());
      r = CFy->Fit(NoiseAmpZero,"S0QW");
    }
    
    //Change things around to get postive amplitude and postive phase
    if ((status!=0)||(NoiseRemoval->GetParameter(1)<0.0)){
      NoiseRemoval->SetParameter(1,-NoiseRemoval->GetParameter(1));
      double temp = NoiseRemoval->GetParameter(3);
      temp = (NoiseRemoval->GetParameter(3)<TMath::Pi())?NoiseRemoval->GetParameter(3)+TMath::Pi():NoiseRemoval->GetParameter(3)-TMath::Pi();
      NoiseRemoval->SetParameter(3,temp);
      r = CFy->Fit(NoiseRemoval,"S0QW");
    }
    if (NoiseRemoval->GetParameter(3)<0.0){
      NoiseRemoval->SetParameter(3.0,NoiseRemoval->GetParameter(3)+2.0*TMath::Pi());
      r = CFy->Fit(NoiseRemoval,"S0QW");
    }
    
    //Save either noise removed frame or just plain old normalized function 
    if (status == 0 && NoiseRemoval->GetChisquare()< NoiseAmpZero->GetChisquare()){
        for (int ix=1; ix<=noise->GetXaxis()->GetNbins(); ix++){
            noise->Fill(noise->GetXaxis()->GetBinCenter(ix),noise->GetYaxis()->GetBinCenter(TimeSlice), NoiseRemoval->Eval(noise->GetXaxis()->GetBinCenter(ix)));
        }
        //Save points to the Parameters graphs 
        HeightGraph->AddPoint(CF->GetYaxis()->GetBinCenter(TimeSlice), NoiseRemoval->GetParameter(0));
        AmpGraph->AddPoint(CF->GetYaxis()->GetBinCenter(TimeSlice), NoiseRemoval->GetParameter(1));
        FreqGraph->AddPoint(CF->GetYaxis()->GetBinCenter(TimeSlice), NoiseRemoval->GetParameter(2));
        PhaseGraph->AddPoint(CF->GetYaxis()->GetBinCenter(TimeSlice), NoiseRemoval->GetParameter(3));
    }
    else{
      double height;
      
      //If this is the first time through subtract 1 to get to 0
      if(times == 1){
          CFy->Fit("pol0", "S0QW");
          height = CFy->GetFunction("pol0")->GetParameter(0);
      } 
      else{height = 0.0;}
      
      for (int ix=1; ix<=noise->GetXaxis()->GetNbins(); ix++){
        noise->Fill(noise->GetXaxis()->GetBinCenter(ix),noise->GetYaxis()->GetBinCenter(TimeSlice), height);
      }
      //Save points to the Parameters graphs 
      HeightGraph->AddPoint(CF->GetYaxis()->GetBinCenter(TimeSlice), height);
      AmpGraph->AddPoint(CF->GetYaxis()->GetBinCenter(TimeSlice), 0.0);
      FreqGraph->AddPoint(CF->GetYaxis()->GetBinCenter(TimeSlice), 0.0);
      PhaseGraph->AddPoint(CF->GetYaxis()->GetBinCenter(TimeSlice), 0.0);
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
  
  //Clone empty CF to save new noise subtracted values to
  TProfile2D* CFMinusNoise = (TProfile2D*)CF->Clone("CFRemoveNoise"); 
  CFMinusNoise->Reset();
  
  //Fill new correlation now with noise subtracted values
  for(int TimeSlice = 1; TimeSlice<=CF->GetYaxis()->GetNbins(); TimeSlice++){
    for (int ix=1; ix<=CF->GetXaxis()->GetNbins(); ix++){
      CFMinusNoise->Fill(CF->GetXaxis()->GetBinCenter(ix),CF->GetYaxis()->GetBinCenter(TimeSlice), CF->GetBinContent(ix, TimeSlice) - noise->GetBinContent(ix, TimeSlice));
    }
  }
  
  //Write and return Correlation Functions without noise 
  CFMinusNoise->SetTitle("Normalized Correlation Function Without Noise; Relative Time (ns); Time");
  CFMinusNoise->Write(Form("CFMinusNoise%d", times));
  return CFMinusNoise;
}
*/

