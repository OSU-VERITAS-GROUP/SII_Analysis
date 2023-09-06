// uniform disk fit gathers all the HBT peaks and baselines, performs quality cuts of HBT fits, and plots and fits the visibility curve 
// ... to measure the angular diameter of the star

// josie added printing OPD and baseline TGraphs and splines - 7 july 2023

//Define Functions
int Setup(int irun, TString& DirName, int& HumanFlag, int& splitNum);
int RelativeTimeOffset(int tx, int ty);
int RelativeTimeOffsetOlder(int tx, int ty);
TF1 * HBTFit(TProfile * HBT, double relTimePar);
TF1 * HBTFitAmp(TH1D * Noise, double relTimePar);
double CalcSNRRatio(TProfile * HBT, TF1 *PeakFit);
void DrawGraph(TGraph * tgtemp, TF1 * Bess, TH1D * UnifDisk, const char * title);

//Initial Fit Parameters
double relTimeWindow(10);  //Maybe Change to 15!!
double sigMin(3), sigMax(6);
double ExpectPeaks(90); 

void UniformDiskFit(){
  //Just Make This Very Large
  int pairNum(4000);  
   
  //Read in from Analysis File
  TTree *header;
  TProfile2D *CFShifted;
  TGraph *Baselines;
  TH1D *ADC1N, *ADC2N;
  TH2D *ADC1, *ADC2;
  int tx, ty, version, NumberOfFrames(0), nx(4), ny(3);
  TDatime * ObsTime = 0;
  float dtWindow;
    
  //Pass to SetUp function
  TString DirName;
  int splitNum(0), HumanFlag(0);
  int nOnCan(0), nOnCanGood(0);
    
  //Create the Canvas for all of the HBT Peaks 
  TCanvas * ManyGraphs = new TCanvas("Many Graphs", "Postage Stamps", 1200, 900);
  ManyGraphs->Divide(nx,ny);
  ManyGraphs->Print("HBTStamps.pdf[");   // note bracket
    
  //Create a Canvas for only the "Good" HBT Peaks   
  TCanvas * ManyGraphsGood = new TCanvas("Many Graphs Good", "Postage Stamps Good", 1200, 900);
  ManyGraphsGood->Divide(nx,ny);
  ManyGraphsGood->Print("HBTStampsGood.pdf["); 

  //Create a Canvas for the entire realtive time window of all the "Good" HBT Peaks 
  TCanvas * FullRelativeTimeGraphsGood = new TCanvas("Full Relative Time", "Full Relative Time ", 1200, 900);
  FullRelativeTimeGraphsGood->Print("FullRelativeTimeHBTStamps.pdf[");   // note bracket
    
  //Plot Style for the graphs 
  gStyle->SetOptFit(1);           
  gStyle->SetOptStat(0);
  gStyle->SetStatX(0.9);          gStyle->SetStatY(0.9);
  gStyle->SetLabelFont(42,"X");   gStyle->SetLabelSize(0.04, "X");
  gStyle->SetLabelFont(42,"Y");   gStyle->SetLabelSize(0.04, "Y");
  gStyle->SetLegendFont(42);      gStyle->SetLegendTextSize(0.04);
  gStyle->SetStatFont(42);        gStyle->SetStatFontSize(0.04);
  gStyle->SetTitleFont(42,"X");   gStyle->SetTitleSize(0.04, "X");
  gStyle->SetTitleFont(42,"Y");   gStyle->SetTitleSize(0.04, "Y");
  gStyle->SetTitleFont(42,"T");   gStyle->SetTitleSize(0.07, "T");
  gStyle->SetStatFormat("6.2g");
  
  //Define Graphs to store the Data in 
  TGraphAsymmErrors* tgGood = new TGraphAsymmErrors(); tgGood->SetMarkerStyle(20); tgGood->SetMarkerColor(4); //Graph of all Points that make the cuts
  TGraphAsymmErrors* tgGoodSingle = new TGraphAsymmErrors(); tgGoodSingle->SetMarkerStyle(24); tgGoodSingle->SetMarkerColor(4); //Graph of all Points that make the cuts
  TGraphAsymmErrors* tgCut  = new TGraphAsymmErrors(); tgCut->SetMarkerStyle(20); tgCut->SetMarkerColor(13); //Graph of points that don't make Muck cut
  TGraphAsymmErrors* tgSNR  = new TGraphAsymmErrors(); tgSNR->SetMarkerStyle(20); tgSNR->SetMarkerColor(5);
  //TGraphAsymmErrors* tgbg = new TGraphAsymmErrors(); tgbg->SetMarkerStyle(20); tgbg->SetMarkerColor(17); //keep track of background points RMS scatter
    
  //Graph Points by telescope pair
  TGraphErrorAsymms* tg12 = new TGraphAsymmErrors(); tg12->SetMarkerStyle(20); tg12->SetMarkerColor(1);
  TGraphAsymmErrors* tg13 = new TGraphAsymmErrors(); tg13->SetMarkerStyle(20); tg13->SetMarkerColor(2);
  TGraphAsymmErrors* tg14 = new TGraphAsymmErrors(); tg14->SetMarkerStyle(20); tg14->SetMarkerColor(3);
  TGraphAsymmErrors* tg23 = new TGraphAsymmErrors(); tg23->SetMarkerStyle(20); tg23->SetMarkerColor(4);
  TGraphAsymmErrors* tg24 = new TGraphAsymmErrors(); tg24->SetMarkerStyle(20); tg24->SetMarkerColor(6);
  TGraphAsymmErrors* tg34 = new TGraphAsymmErrors(); tg34->SetMarkerStyle(20); tg34->SetMarkerColor(7);
    
  //Graph Points by month
  TGraphAsymmErrors * tgmonth1 = new TGraphAsymmErrors(); tgmonth1->SetMarkerStyle(20); tgmonth1->SetMarkerColor(1);
  TGraphAsymmErrors * tgmonth2 = new TGraphAsymmErrors(); tgmonth2->SetMarkerStyle(20); tgmonth2->SetMarkerColor(2);
  TGraphAsymmErrors * tgmonth3 = new TGraphAsymmErrors(); tgmonth3->SetMarkerStyle(20); tgmonth3->SetMarkerColor(3);
  TGraphAsymmErrors * tgmonth4 = new TGraphAsymmErrors(); tgmonth4->SetMarkerStyle(20); tgmonth4->SetMarkerColor(4);
    
  //Distribution histograms of background points
  TH1D * BackgroundPointProjectionCut = new TH1D("Cut Background Projection","Cut Background Projection", 50, -10e-6, 10e-6);
  TH1D * BackgroundPointProjectionAll = new TH1D("All Background Projection","All Background Projection", 50, -10e-6, 10e-6);
  TH1D * DataPointProjection = new TH1D("DataPointProjection","Data Point Projection", 72, -3.8e-6, 3.4e-6);
  TH1D * GoodPointsSNR = new TH1D("GoodPointsSNR","GoodPointsSNR", 24, 0, 6);
  TH1D * CutPointsSNR = new TH1D("CutPointsSNR","CutPointsSNR", 24, 0, 6);
    
  //Iterate over these for the different graphs 
  int nptsGood(0), nptsbg(0), nptsCut(0), nptsSNR(0);
  int nptsT1T2(0), nptsT1T3(0), nptsT1T4(0), nptsT2T3(0), nptsT2T4(0), nptsT3T4(0);
  int nptsmonth1(0), nptsmonth2(0), nptsmonth3(0), nptsmonth4(0);
    
  //Muck Parameter 
  double RMSMuck(0);
    
  //Clump Long Baselines 
  double SumAveVisClump1(0), SumAveVisClump2(0), SumAveVisClump3(0);
  double SumErrVisClump1(0), SumErrVisClump2(0), SumErrVisClump3(0);
  double SumAveBaseClump1(0), SumAveBaseClump2(0), SumAveBaseClump3(0);
  double SumErrBaseClump1(0), SumErrBaseClump2(0), SumErrBaseClump3(0);
  double nptsClump1(0), nptsClump2(0), nptsClump3(0);
    
//------------------------------------------------------------------------------
//             Iterate over all Directories to find RMS Muck
//------------------------------------------------------------------------------
 
  for(int idir = 0; idir < pairNum; idir++){
    cout << "Working on RMS Muck for " << idir << endl;
    
    //Read in directory tree
    if (Setup(idir, DirName, HumanFlag, splitNum)==0) break;
        
    //Read in information from output file for each run 
    TFile * CorrInfo = new TFile(Form("%s/Analysis.root", DirName.Data()), "READONLY");
    CorrInfo->GetObject("Header",header);
    header->SetBranchAddress("VersionBr",&version);
    header->SetBranchAddress("T1Br",&tx);       
    header->SetBranchAddress("T2Br",&ty);
    header->SetBranchAddress("tWindowSecBr",&dtWindow);       
    header->SetBranchAddress("LocalTimeBr" ,&ObsTime);
    header->SetBranchAddress("NumWindowsBr",&NumberOfFrames); // Number of Frames
    header->GetEntry(0);
    CorrInfo->GetObject("BaselineGraph", Baselines);
    CorrInfo->GetObject("JustShift", CFShifted);
        
    //Read in ADCs and make sure nothing is empty 
    if(version > 0){
      CorrInfo->GetObject("ADC1", ADC1N); 
      CorrInfo->GetObject("ADC2", ADC2N);
      if(Baselines == 0 || CFShifted == 0 || ADC1N == 0 || ADC2N == 0){cout << "Oh no we have a problem!! " << DirName << endl;}
    }
    if(version < 0){
      CorrInfo->GetObject("ADC1", ADC1);
      CorrInfo->GetObject("ADC2", ADC2);
      if(Baselines == 0 || CFShifted == 0 || ADC1 == 0 || ADC2 == 0){cout << "Oh no we have a problem!! " << DirName << endl;};
    }
        
    //Get Relative Time Offset for a Given Pair of Telescopes
    double relTimePar = RelativeTimeOffset(tx, ty); 

    //Figure out how to split the run
    int start(1), end(CFShifted->GetYaxis()->GetNbins()/splitNum);
    for(int sectionloop = 0; sectionloop < splitNum; ++sectionloop){
      TProfile * HBTProjbg = CFShifted->ProfileX("HBTProjbg", start, end);
    
      //Calculate Weighted Baselines
      double sumAveNum(0), sumAveDen(0), sumRMSTerm1Num(0);
      if(version > 0){
        for(int TimeSlice = start; TimeSlice <= end; TimeSlice++){
          sumAveNum += ADC1N->GetBinContent(TimeSlice)*ADC2N->GetBinContent(TimeSlice)*Baselines->GetPointY(TimeSlice);
          sumAveDen += ADC1N->GetBinContent(TimeSlice)*ADC2N->GetBinContent(TimeSlice);
          sumRMSTerm1Num += ADC1N->GetBinContent(TimeSlice)*ADC2N->GetBinContent(TimeSlice)*Baselines->GetPointY(TimeSlice)*Baselines->GetPointY(TimeSlice);
        }
      }
      if(version < 0){
        for(int TimeSlice = start; TimeSlice <= end; TimeSlice++){
          TH1D* ADC1y = ADC1->ProjectionX("temp1",TimeSlice,TimeSlice);
          TH1D* ADC2y = ADC2->ProjectionX("temp2",TimeSlice,TimeSlice);
          sumAveNum += ADC1y->GetMean()*ADC2y->GetMean()*Baselines->GetPointY(TimeSlice);
          sumAveDen += ADC1y->GetMean()*ADC2y->GetMean();
          sumRMSTerm1Num += ADC1y->GetMean()*ADC2y->GetMean()*Baselines->GetPointY(TimeSlice)*Baselines->GetPointY(TimeSlice);
        }
      }
            
      //Calculate the points and their errors
      //double BasePoint  = abs(sumAveNum/sumAveDen);
            
            
      for(int check = 14; check < HBTProjbg->GetXaxis()->GetNbins(); check +=3){
        if(HBTProjbg->GetBinCenter(check) > -100 && HBTProjbg->GetBinCenter(check) <  100){continue;}
        double relTimeParbg = HBTProjbg->GetBinCenter(check);
        TF1* PeakFitbg = new TF1("HBTpeak","([0]*exp(-pow((x-[1])/[2],2)/2.0))",-300,300); 
        PeakFitbg->SetParName(0,"Amp");           PeakFitbg->SetParameter(0,0.0);     
        PeakFitbg->SetParName(1,"#tau_{0} (ns)"); PeakFitbg->FixParameter(1,relTimeParbg); 
        PeakFitbg->SetParName(2,"#sigma (ns)");   PeakFitbg->SetParameter(2,4.5);        PeakFitbg->SetParLimits(2, sigMin,  sigMax);
        TFitResultPtr fitResultbg = HBTProjbg->Fit(PeakFitbg, "SQN");
        double PointValbg = PeakFitbg->GetParameter(0)*PeakFitbg->GetParameter(2)*sqrt(2*TMath::Pi());
        if(abs(PointValbg) > 10e-6){continue;}
        BackgroundPointProjectionAll->Fill(PointValbg);
        if(abs(PeakFitbg->GetParameter(2)- sigMax) > 0.01 && abs(PeakFitbg->GetParameter(2)- sigMin) < 0.01 && (abs(PeakFitbg->GetParameter(1)-relTimePar-relTimeWindow)) > 0.01 && (abs(PeakFitbg->GetParameter(1) - relTimePar+relTimeWindow)) > 0.01){
          //tgbg->AddPoint(BasePoint, PointValbg); 
          ++nptsbg; 
          RMSMuck += pow(PointValbg, 2);
          BackgroundPointProjectionCut->Fill(PointValbg);
        }
      }
    }
    CorrInfo->Close();
  }
    
    //Find RMS and plot on Graph
    RMSMuck = sqrt(RMSMuck/nptsbg);
    
    //Define Muck Boxes
    TBox * MuckBounds = new TBox(-0.05,-RMSMuck,200,RMSMuck);            MuckBounds->SetFillColorAlpha(15, 0.65);
    TBox * MuckBoundsOutter = new TBox(-0.05,-2*RMSMuck,200,2*RMSMuck);  MuckBoundsOutter->SetFillColorAlpha(15, 0.35);
    
//------------------------------------------------------------------------------
//             Iterate over all Directories to Get HBT Peaks
//------------------------------------------------------------------------------

  //Define file to write to throuhout loop
  ofstream GoodDirectories, ZeroDirectories, RelTauPars, flags;
  GoodDirectories.open ("GoodDirectories.txt");
  ZeroDirectories.open ("ZeroDirectories.txt");
  RelTauPars.open ("RelativeTauList.txt");
  flags.open("checkFlags.txt");
  flags << "txty  flag  par0   par1   par2" << endl;
    
    
  //Loop over all the runs/pairs 
  for(int idir = 0; idir < pairNum; idir++){

    cout << "Analyzing Data for " << idir << endl;
        
    //Read in directory tree 
    if (Setup(idir, DirName, HumanFlag, splitNum)==0) break;

    //Make Set HumanFlags Have Overriding Importance 
    double tempImport(0); 
    if(HumanFlag !=0){tempImport = HumanFlag;}
        
    //Read in information from output file for each run 
    TFile * CorrInfo = new TFile(Form("%s/Analysis.root", DirName.Data()), "READONLY");
    CorrInfo->GetObject("Header",header);
    header->SetBranchAddress("VersionBr",&version);
    header->SetBranchAddress("T1Br",&tx);       
    header->SetBranchAddress("T2Br",&ty);
    header->SetBranchAddress("tWindowSecBr",&dtWindow);       
    header->SetBranchAddress("LocalTimeBr" ,&ObsTime);
    header->SetBranchAddress("NumWindowsBr",&NumberOfFrames); // Number of Frames
    header->GetEntry(0);
    CorrInfo->GetObject("BaselineGraph", Baselines);
    CorrInfo->GetObject("JustShift", CFShifted);   
        
    //Read in ADCs and make sure nothing is empty 
    if(version > 0){
      CorrInfo->GetObject("ADC1", ADC1N); 
      CorrInfo->GetObject("ADC2", ADC2N);
      if(Baselines == 0 || CFShifted == 0 || ADC1N == 0 || ADC2N == 0){cout << "Oh no we have a problem!! " << DirName << endl;}
    }
    if(version < 0){

      CorrInfo->GetObject("ADC1", ADC1);
      CorrInfo->GetObject("ADC2", ADC2);
      if(Baselines == 0 || CFShifted == 0 || ADC1 == 0 || ADC2 == 0){cout << "Oh no we have a problem!! " << DirName << endl;};
    }
        
    //Get Relative Time Offset for this pair of telescopes
    double relTimePar = RelativeTimeOffset(tx, ty); 
        
    //If the run is longer than an hour split the run in half 
    if(dtWindow*NumberOfFrames > 3600 && splitNum == 1){splitNum = 2;}
    if(DirName == "UTC20211217/m12d16h23/T3T4" || DirName == "UTC20211217/m12d17h01/T3T4" || DirName == "UTC20220212/m02d11h20/T3T4" || DirName == "UTC20220214/m02d13h19/T3T4"){
        splitNum = 3; cout << "split = 3, dir: " << DirName << endl;}
    
    //Pull year, month, and day of the run      
    int year(ObsTime->GetYear()), month(ObsTime->GetMonth()), day(ObsTime->GetDay()), start(1), end(CFShifted->GetYaxis()->GetNbins()/splitNum); 
        
    //Loop over all splits
    for(int sectionloop = 0; sectionloop < splitNum; ++sectionloop){
      HumanFlag = 0;
      nOnCan++;
      ManyGraphs->cd(nOnCan); 

      //Fit HBT Peak
      TProfile * HBTProj = CFShifted->ProfileX(Form("HBTProj%d%d", idir, sectionloop), start, end);
      TF1 * PeakFit = HBTFit(HBTProj, relTimePar);

      // Calcuate Noise Level for a given runpair
      double SNRratio = CalcSNRRatio(HBTProj, PeakFit);
      double BaseErrorLow  = 1000;
      double BaseErrorHigh = -1;
      
      //Calculate Weighted Baselines
      double sumAveNum(0), sumAveDen(0), sumRMSTerm1Num(0);
      if(version > 0){ 
        for(int TimeSlice = start; TimeSlice <= end; TimeSlice++){
          if(Baselines->GetPointY(TimeSlice) < BaseErrorLow){BaseErrorLow = Baselines->GetPointY(TimeSlice);}
          if(Baselines->GetPointY(TimeSlice) > BaseErrorHigh){BaseErrorHigh = Baselines->GetPointY(TimeSlice);}
          sumAveNum += ADC1N->GetBinContent(TimeSlice)*ADC2N->GetBinContent(TimeSlice)*Baselines->GetPointY(TimeSlice);
          sumAveDen += ADC1N->GetBinContent(TimeSlice)*ADC2N->GetBinContent(TimeSlice);
          sumRMSTerm1Num += ADC1N->GetBinContent(TimeSlice)*ADC2N->GetBinContent(TimeSlice)*Baselines->GetPointY(TimeSlice)*Baselines->GetPointY(TimeSlice);
        }
      }
      if(version < 0){
        for(int TimeSlice = start; TimeSlice <= end; TimeSlice++){
          if(Baselines->GetPointY(TimeSlice) < BaseErrorLow){BaseErrorLow = Baselines->GetPointY(TimeSlice);}
          if(Baselines->GetPointY(TimeSlice) > BaseErrorHigh){BaseErrorHigh = Baselines->GetPointY(TimeSlice);}
          TH1D* ADC1y = ADC1->ProjectionX("temp1",TimeSlice,TimeSlice);
          TH1D* ADC2y = ADC2->ProjectionX("temp2",TimeSlice,TimeSlice);
          sumAveNum += ADC1y->GetMean()*ADC2y->GetMean()*Baselines->GetPointY(TimeSlice);
          sumAveDen += ADC1y->GetMean()*ADC2y->GetMean();
          sumRMSTerm1Num += ADC1y->GetMean()*ADC2y->GetMean()*Baselines->GetPointY(TimeSlice)*Baselines->GetPointY(TimeSlice);
        }
      }
            
      //Calculate the points and their errors
      double BasePoint  = abs(sumAveNum/sumAveDen);
      //double BaseError  = sqrt(sumRMSTerm1Num/sumAveDen - pow(BasePoint, 2));
      double PointVal   = PeakFit->GetParameter(0); 
      double PointError = PeakFit->GetParError(0);  
            
      //Declare Human Flags
      if(PeakFit->GetParameter(0) < 0){HumanFlag = 1;} //Peak went negative
      if(abs(PeakFit->GetParameter(2)- sigMax) < 0.01){HumanFlag = 2;} //Peak had too large of a width
      if(abs(PeakFit->GetParameter(2)- sigMin) < 0.01){HumanFlag = 3;} //Peak had too small of a width
      if((abs(PeakFit->GetParameter(1)-relTimePar-relTimeWindow)) < 0.01 || (abs(PeakFit->GetParameter(1) - relTimePar+relTimeWindow)) < 0.01){HumanFlag = 4;} //Peak hit limits on possible region
      if(tempImport != 0){HumanFlag = tempImport;}
      
      flags << tx << ty << "  " << HumanFlag << "  " << PeakFit->GetParameter(0) << "  " << PeakFit->GetParameter(1) << "  " << PeakFit->GetParameter(2) << endl;
      
      //Set Title and Offset
      HBTProj->SetTitle(Form("%d.%d) %s B =%6.1f Flag = %d ;Relative Time (ns); g^{2}(0)-1", idir, sectionloop, DirName.Data(), BasePoint, HumanFlag));
      HBTProj->GetYaxis()->SetTitleOffset(1.0); 
      
      //Skip saving points unless if they fail cuts            
      if(HumanFlag < 0 || HumanFlag > 1){goto cont;}
      
      //Save "good" full range HBT Peaks
      FullRelativeTimeGraphsGood->cd(); 
      HBTProj->Draw();
      FullRelativeTimeGraphsGood->Print("FullRelativeTimeHBTStamps.pdf"); 

      if(BasePoint > 90  && BasePoint <= 105){SumAveVisClump1 += PointVal; SumErrVisClump1 += pow(PointVal,2); SumAveBaseClump1 += BasePoint; SumErrBaseClump1 += pow(BasePoint,2); ++nptsClump1; ZeroDirectories << "/betUMa/" <<  DirName << endl;}
      if(BasePoint > 105 && BasePoint <= 150){SumAveVisClump2 += PointVal; SumErrVisClump2 += pow(PointVal,2); SumAveBaseClump2 += BasePoint; SumErrBaseClump2 += pow(BasePoint,2); ++nptsClump2;}
      if(BasePoint > 150)                    {SumAveVisClump3 += PointVal; SumErrVisClump3 += pow(PointVal,2); SumAveBaseClump3 += BasePoint; SumErrBaseClump3 += pow(BasePoint,2); ++nptsClump3;}
            
      //-----------------This Rationale in general will not work--------------------
      //----------------------Needs to be better-----------------------------------
      //Save all points until our uncertainty
      if(BasePoint <= ExpectPeaks){ //&& PointVal > RMSMuck  && PeakFit->GetParameter(0) > 2e-6
        gPad->SetFrameFillColor(33);
        tgGood->AddPoint(BasePoint, PointVal);
        tgGood->SetPointError(nptsGood, BaseErrorLow, BaseErrorHigh, PointVal-PointError, PointVal+PointError);
        ++nptsGood;
        GoodDirectories <<  DirName << "   " << splitNum << "   " << sectionloop << endl;
        GoodPointsSNR->Fill(SNRratio);
      }
      else{
        tgCut->AddPoint(BasePoint, PointVal);
        ++nptsCut;
        CutPointsSNR->Fill(SNRratio);
        if(BasePoint < 110){DataPointProjection->Fill(PointVal);}
      }
      
      if(SNRratio > 1.2){   
        tgSNR->AddPoint(BasePoint, PointVal);
        tgSNR->SetPointError(nptsSNR, BaseErrorLow, BaseErrorHigh, PointVal-PointError, PointVal+PointError);
        ++nptsSNR;
      }
      //---------------------------------------------------------------------------
      //---------------------------------------------------------------------------
            
      // Sort by Pair
      if((tx==1 && ty==2) || (tx==2 && ty==1)){tg12->AddPoint(BasePoint, PointVal); tg12->SetPointError(nptsT1T2, BaseErrorLow, BaseErrorHigh, PointVal-PointError, PointVal+PointError); ++nptsT1T2;}
      if((tx==1 && ty==3) || (tx==3 && ty==1)){tg13->AddPoint(BasePoint, PointVal); tg13->SetPointError(nptsT1T3, BaseErrorLow, BaseErrorHigh, PointVal-PointError, PointVal+PointError); ++nptsT1T3;}
      if((tx==1 && ty==4) || (tx==4 && ty==1)){tg14->AddPoint(BasePoint, PointVal); tg14->SetPointError(nptsT1T4, BaseErrorLow, BaseErrorHigh, PointVal-PointError, PointVal+PointError); ++nptsT1T4;}
      if((tx==2 && ty==3) || (tx==3 && ty==2)){tg23->AddPoint(BasePoint, PointVal); tg23->SetPointError(nptsT2T3, BaseErrorLow, BaseErrorHigh, PointVal-PointError, PointVal+PointError); ++nptsT2T3;}
      if((tx==2 && ty==4) || (tx==4 && ty==2)){tg24->AddPoint(BasePoint, PointVal); tg24->SetPointError(nptsT2T4, BaseErrorLow, BaseErrorHigh, PointVal-PointError, PointVal+PointError); ++nptsT2T4;}
      if((tx==3 && ty==4) || (tx==3 && ty==4)){tg34->AddPoint(BasePoint, PointVal); tg34->SetPointError(nptsT3T4, BaseErrorLow, BaseErrorHigh, PointVal-PointError, PointVal+PointError); ++nptsT3T4;}
      
      //Sort By Date
      if(year == 2021 && month == 12){tgmonth1->AddPoint(BasePoint, PointVal); tgmonth1->SetPointError(nptsmonth1, BaseErrorLow, BaseErrorHigh, PointVal-PointError, PointVal+PointError); ++nptsmonth1;}
      if(year == 2022 && month == 02){tgmonth2->AddPoint(BasePoint, PointVal); tgmonth2->SetPointError(nptsmonth2, BaseErrorLow, BaseErrorHigh, PointVal-PointError, PointVal+PointError); ++nptsmonth2;}
      if(year == 2022 && month == 03){tgmonth3->AddPoint(BasePoint, PointVal); tgmonth3->SetPointError(nptsmonth3, BaseErrorLow, BaseErrorHigh, PointVal-PointError, PointVal+PointError); ++nptsmonth3;}   
      if(year == 2022 && month == 05){tgmonth4->AddPoint(BasePoint, PointVal); tgmonth4->SetPointError(nptsmonth4, BaseErrorLow, BaseErrorHigh, PointVal-PointError, PointVal+PointError); ++nptsmonth4;}  

      //Print Publication Quality HBT Peaks 
      if((idir == 5 && sectionloop == 0) || (idir == 41 && sectionloop ==1)){
        gStyle->SetOptFit(0);
        TCanvas * GoodGraphs = new TCanvas("Good Graphs", "Postage Stamps", 1200, 900);
        HBTProj->SetTitle(";Relative Time (ns); g^{2}(0)-1");
        HBTProj->Draw(); //left->Draw(); right->Draw();
        TPaveText *statsbox=new TPaveText(0.63,0.7,0.9,0.9,"brNDC");
        statsbox->SetBorderSize(1);
        statsbox->SetTextAlign(12);
        statsbox->SetTextFont(42); statsbox->SetTextSize(0.033);
        statsbox->AddText(Form("#chi^{2}/ndf   %3.1lf / %d",PeakFit->GetChisquare(),PeakFit->GetNDF()));
        statsbox->AddText(Form("A   %3.1le #pm %3.1le",PeakFit->GetParameter(0), PeakFit->GetParError(0)));
        statsbox->AddText(Form("#tau_{0} (ns)   %3.1lf #pm %3.1lf",PeakFit->GetParameter(1), PeakFit->GetParError(1)));
        statsbox->AddText(Form("#sigma (ns)   %3.1lf #pm %3.1lf",PeakFit->GetParameter(2), PeakFit->GetParError(2)));
        statsbox->SetFillColor(0);
        statsbox->Draw();
        GoodGraphs->Print(Form("%d.%d.pdf", idir, sectionloop));
        gStyle->SetOptFit(1);
      }
      
      //Make the HBT Canvas image
      cont:
      ManyGraphs->cd(nOnCan); 
      HBTProj->SetMarkerColor(1);
      HBTProj->SetMaximum(1e-6); HBTProj->SetMinimum(-0.8e-6);  
      HBTProj->Rebin(4); HBTProj->GetXaxis()->SetRangeUser(-64,64); //This must come after the Rebin(4);
      HBTProj->GetYaxis()->SetTitleOffset(1.0); 
      HBTProj->Draw();
      
      TLine * left = new TLine(relTimePar-relTimeWindow,HBTProj->GetMinimum(),relTimePar-relTimeWindow,HBTProj->GetMaximum());  left->SetLineWidth(2);  left->Draw();
      TLine * right = new TLine(relTimePar+relTimeWindow,HBTProj->GetMinimum(),relTimePar+relTimeWindow,HBTProj->GetMaximum()); right->SetLineWidth(2); right->Draw();
      //TLine * Muck = new TLine(-256, AveBackgroundAmp, 256, AveBackgroundAmp);  Muck->SetLineWidth(2);  Muck->Draw();
      //TLine * Signal = new TLine(-256, PeakAmpReal, 256, PeakAmpReal); Signal->SetLineWidth(2); Signal->Draw();
      gPad->Update(); TPaveStats *st = (TPaveStats*)HBTProj->FindObject("stats"); st->Draw();

      //Save only the "Good" peaks to a seperate Canvas
      if(BasePoint <= ExpectPeaks && (HumanFlag == 0 || HumanFlag == 1)){
        nOnCanGood++;
        ManyGraphsGood->cd(nOnCanGood); 
        HBTProj->Draw();
        TLine * left = new TLine(relTimePar-relTimeWindow,HBTProj->GetMinimum(),relTimePar-relTimeWindow,HBTProj->GetMaximum());  left->SetLineWidth(2);  left->Draw();
        TLine * right = new TLine(relTimePar+relTimeWindow,HBTProj->GetMinimum(),relTimePar+relTimeWindow,HBTProj->GetMaximum()); right->SetLineWidth(2); right->Draw();
        gPad->Update(); TPaveStats *st = (TPaveStats*)HBTProj->FindObject("stats"); st->Draw();
      }
            
      //Draw Read line if point is removed 
      if(HumanFlag < 0 || HumanFlag > 1){
        TLine * BadRun = new TLine(HBTProj->GetXaxis()->GetXmin(), HBTProj->GetMinimum(), HBTProj->GetXaxis()->GetXmax(), HBTProj->GetMaximum());
        BadRun->SetLineColor(kRed); BadRun->Draw();
      }

      //Reset HBT Stamps if needed   
      if (nOnCanGood==nx*ny){ManyGraphsGood->Print("HBTStampsGood.pdf"); nOnCanGood=0; ManyGraphsGood->Clear(); ManyGraphsGood->Divide(nx,ny);}
      if (nOnCan==nx*ny){ManyGraphs->Print("HBTStamps.pdf"); nOnCan=0; ManyGraphs->Clear();    ManyGraphs->Divide(nx,ny);}
      start += CFShifted->GetYaxis()->GetNbins()/splitNum;
      end   += CFShifted->GetYaxis()->GetNbins()/splitNum;
    }   
  }  
    
  //Close Files and Stamps that Save Info about each runpair  
  GoodDirectories.close();
  ZeroDirectories.close();
  RelTauPars.close();
  flags.close();
  ManyGraphs->Print("HBTStamps.pdf");
  ManyGraphsGood->Print("HBTStampsGood.pdf)");
  FullRelativeTimeGraphsGood->Print("FullRelativeTimeHBTStamps.pdf)"); 
    
    SumAveVisClump1 /= nptsClump1;
    SumAveVisClump2 /= nptsClump2;
    SumAveVisClump3 /= nptsClump3;
    
    SumErrVisClump1 = sqrt(SumErrVisClump1/nptsClump1 - pow(SumAveVisClump1,2));
    SumErrVisClump2 = sqrt(SumErrVisClump2/nptsClump2 - pow(SumAveVisClump2,2));
    SumErrVisClump3 = sqrt(SumErrVisClump3/nptsClump3 - pow(SumAveVisClump3,2));
    
    SumErrVisClump1 /= sqrt(nptsClump1);
    SumErrVisClump2 /= sqrt(nptsClump2);
    SumErrVisClump3 /= sqrt(nptsClump3);
    
    
    SumAveBaseClump1 /= nptsClump1;
    SumAveBaseClump2 /= nptsClump2;
    SumAveBaseClump3 /= nptsClump3;
    
    SumErrBaseClump1 = sqrt(SumErrBaseClump1/nptsClump1 - pow(SumAveBaseClump1,2));
    SumErrBaseClump2 = sqrt(SumErrBaseClump2/nptsClump2 - pow(SumAveBaseClump2,2));
    SumErrBaseClump3 = sqrt(SumErrBaseClump3/nptsClump3 - pow(SumAveBaseClump3,2));
    
    SumErrBaseClump1 /= sqrt(nptsClump1);
    SumErrBaseClump2 /= sqrt(nptsClump2);
    SumErrBaseClump3 /= sqrt(nptsClump3);
    
    
    cout << "SumAveVisClump1:" << SumAveVisClump1 << "SumErrVisClump1:" << SumErrVisClump1 << "SumAveBaseClump1:" << SumAveBaseClump1 << "SumErrVisClump1:" << SumErrBaseClump1 << endl;
    cout << "SumAveVisClump2:" << SumAveVisClump2 << "SumErrVisClump2:" << SumErrVisClump2 << "SumAveBaseClump2:" << SumAveBaseClump2 << "SumErrVisClump2:" << SumErrBaseClump2 << endl;
    cout << "SumAveVisClump3:" << SumAveVisClump3 << "SumErrVisClump3:" << SumErrVisClump3 << "SumAveBaseClump3:" << SumAveBaseClump3 << "SumErrVisClump3:" << SumErrBaseClump3 << endl;
    
    
    tgGood->AddPoint(SumAveBaseClump1, 0); tgGood->SetPointError(nptsGood, SumErrBaseClump1, RMSMuck);
    tgGoodSingle->AddPoint(SumAveBaseClump1, 0); tgGoodSingle->SetPointError(0, SumErrBaseClump1, RMSMuck);
    
   // tgGood->AddPoint(SumAveBaseClump2, 0); tgGood->SetPointError(nptsGood, SumErrBaseClump2, RMSMuck); ++nptsGood;
    //tgGood->AddPoint(SumAveBaseClump3, 0); tgGood->SetPointError(nptsGood, SumErrBaseClump3, RMSMuck); ++nptsGood;
    

    //TF1* BessFit = new TF1("Bess","fabs([0])*pow(2.0*TMath::BesselJ1(TMath::Pi()*x*[1]/(416e-9*2.063e8))/(TMath::Pi()*x*[1]/(416e-9*2.063e8)),2)",0,200);
    TF1* BessFit = new TF1("Bess","fabs([0])*pow(2.0*TMath::BesselJ1(TMath::Pi()*x*[1]*1e-3/(4157e-10*206265))/(TMath::Pi()*x*[1]*1e-3/(4157e-10*206265)),2)",0,200);
    BessFit->SetParName(0,"Normalization"); //g^{2}(0)-1
    BessFit->SetParName(1,"#theta_{UD} (mas)");
    BessFit->SetParameter(0,1e-5); 
    BessFit->SetParameter(1,1);  
    tgGood->Fit(BessFit, "EX0");
    tgGood->RemovePoint(nptsGood);
 
    TH1D* UniformDiskModel = new TH1D("myFit","",100,-0.05,200);
    UniformDiskModel->SetMinimum(-4e-6);
    UniformDiskModel->SetMaximum(1.5e-5);
    UniformDiskModel->GetXaxis()->SetTitle("Projected Baseline (m)");
    UniformDiskModel->GetYaxis()->SetTitle("Area Under g^{2}(0)-1 Curve (ns)"); UniformDiskModel->GetYaxis()->SetTitleOffset(1.1); //g^{2}(0)-1
    
    TGaxis * NormalizedAxis = new TGaxis(200,-4e-6,200,1.5e-5,-4e-6/BessFit->GetParameter(0),1.5e-5/BessFit->GetParameter(0), 510,"+L");
    NormalizedAxis->SetTitle("Model Dependent Squared Visibilities"); 
    NormalizedAxis->SetLabelFont(42); NormalizedAxis->SetLabelSize(0.04); NormalizedAxis->SetLabelOffset(0.01);
    NormalizedAxis->SetTitleFont(42); NormalizedAxis->SetTitleSize(0.04); NormalizedAxis->SetTitleOffset(1.1);
    
    TPaveText *r=new TPaveText(0.50,0.7,0.9,0.9,"brNDC");
    r->SetBorderSize(1);
    r->SetTextFont(42); r->SetTextSize(0.035);
    r->AddText(Form("#chi^{2}/ndf   %3.1lf / %d",BessFit->GetChisquare(),BessFit->GetNDF()));
    r->AddText(Form("Normalization  %3.1le #pm %3.1le",BessFit->GetParameter(0), BessFit->GetParError(0)));
    r->AddText(Form("#theta_{UD} (mas)   %3.2lf #pm %3.2lf",BessFit->GetParameter(1), BessFit->GetParError(1)));
    r->SetFillColor(kWhite);
    gStyle->SetOptFit(1);
    
    ofstream out("PointsInVisibityCurve.txt");
    for(int i=0; i < tgGood->GetN(); ++i){out << tgGood->GetPointX(i) << "    " <<  tgGood->GetErrorX(i) << "    " << tgGood->GetPointY(i)  << "   " << tgGood->GetErrorY(i) << "   " << tgGood->GetPointY(i)/BessFit->GetParameter(0) << "   " <<  tgGood->GetErrorY(i)/BessFit->GetParameter(0) << endl;}

    //---------------- All Points ------------------
    TCanvas * AllPoints = new TCanvas("All Points","All Points",1200,900);
    UniformDiskModel->Draw();
    //MuckBounds->Draw("Same");
    //MuckBoundsOutter->Draw("Same");
    tgGood->Draw("P, same");
    tgCut->Draw("P, same");
    tgGoodSingle->Draw("P, same");
    BessFit->Draw("same");
    r->Draw();
    NormalizedAxis->Draw("same");
    AllPoints->Print("HBTStamps.pdf");

    //-------------- Pair Points ---------------
    TCanvas *TelescopePoints = new TCanvas("Telescope Pairs","Telescope Pairs",1200,900);
    
    TLegend * legend = new TLegend(0.65,0.65,0.9,0.9);
    legend->SetHeader("Telescope Pairs","C"); // option "C" allows to center the header
    legend->AddEntry(tg12,"T1T2","p"); legend->AddEntry(tg13,"T1T3","p");
    legend->AddEntry(tg14,"T1T4","p"); legend->AddEntry(tg23,"T2T3","p");
    legend->AddEntry(tg24,"T2T4","p"); legend->AddEntry(tg34,"T3T4","p");
    
    UniformDiskModel->Draw();
    UniformDiskModel->SetTitle("Bet UMa Points by Telescope Pair");
    tg12->Draw("P, same"); tg13->Draw("P, same");
    tg14->Draw("P, same"); tg23->Draw("P, same");
    tg24->Draw("P, same"); tg34->Draw("P, same");
    legend->Draw("same");  BessFit->Draw("same");
    TelescopePoints->Print("HBTStamps.pdf");
    
    // ---------- Date Points ---------------
    TCanvas *DatePoints = new TCanvas("Dates","Dates",1200,900);

    TLegend * legend2 = new TLegend(0.7,0.7,0.9,0.9);
    legend2->SetHeader("Telescope Pairs","C"); // option "C" allows to center the header
    legend2->AddEntry(tgmonth1,"Dec 2021","p"); legend2->AddEntry(tgmonth2,"Feb 2022","p");
    legend2->AddEntry(tgmonth3,"Mar 2022","p"); legend2->AddEntry(tgmonth4,"May 2022","p");
    
    UniformDiskModel->Draw();
    UniformDiskModel->SetTitle("Uniform Disk Visibility Curve for Bet UMa");
    tgmonth1->Draw("P, same"); tgmonth2->Draw("P, same");
    tgmonth3->Draw("P, same"); tgmonth4->Draw("P, same");
    legend2->Draw("same");     BessFit->Draw("same");    
    DatePoints->Print("HBTStamps.pdf");
    
    //---------- Individual Pairs ------------
    DrawGraph(tgSNR, BessFit, UniformDiskModel, "SNR Cuts");
    DrawGraph(tg12, BessFit, UniformDiskModel, "T1T2 Data Points Only");
    DrawGraph(tg13, BessFit, UniformDiskModel, "T1T3 Data Points Only"); 
    DrawGraph(tg14, BessFit, UniformDiskModel, "T1T4 Data Points Only");
    DrawGraph(tg23, BessFit, UniformDiskModel, "T2T3 Data Points Only");
    DrawGraph(tg24, BessFit, UniformDiskModel, "T2T4 Data Points Only");
    DrawGraph(tg34, BessFit, UniformDiskModel, "T3T4 Data Points Only");
    
    //--------- Individual Dates ------------
    DrawGraph(tgmonth1, BessFit, UniformDiskModel, "Dec 2021 Data Points Only");
    DrawGraph(tgmonth2, BessFit, UniformDiskModel, "Feb 2022 Data Points Only"); 
    DrawGraph(tgmonth3, BessFit, UniformDiskModel, "Mar 2022 Data Points Only");
    DrawGraph(tgmonth4, BessFit, UniformDiskModel, "May 2022 Data Points Only");
    
    TCanvas * BackGroundPoint = new TCanvas("All Background Points","All Background Points",1200,900);
    BackgroundPointProjectionAll->Draw();
    BackGroundPoint->Print("HBTStamps.pdf");
    
    TCanvas * BackGroundPointCut = new TCanvas("Cut Background Points","Cut Background Points",1200,900);
    BackgroundPointProjectionCut->Draw();
    BackGroundPointCut->Print("HBTStamps.pdf");
    
    TFile* outputFile = new TFile("DataPointProjection.root","RECREATE");
    DataPointProjection->Write();
    outputFile->Close();

    TCanvas * DataPointCan = new TCanvas("Long Baseline Point Projection","Long Baseline Point Projection",1200,900);
    DataPointProjection->Draw();
    DataPointCan->Print("HBTStamps.pdf");

    
    TCanvas * SNRDist = new TCanvas("SNRDist","SNRDist",1200,900);
    CutPointsSNR->SetLineColor(kRed);
    CutPointsSNR->Draw();
    GoodPointsSNR->Draw("SAME");
    SNRDist->Print("HBTStamps.pdf");

    DatePoints->Print("HBTStamps.pdf]");
    
    //Draw HeatMaps on Canvas
    int nOnCanHeat(0); int nOnCanFFT(0); int nOnCanDB(0);
    TProfile2D* HeatMap; TH1D* FFT;
    TGraph* pythonDelays; // josie 7 july 23
    TGraph* Delays;
    TGraph* pythonBaselines;
    TGraph* Baseline; 
    TCanvas * ManyHeatMaps = new TCanvas("Many HeatMaps", "Postage Maps", 1600, 1200);
    TCanvas * ManyFFTs = new TCanvas("Many FFTs", "Many FFTs", 1600, 1200);
    TCanvas * ManyDelays = new TCanvas("ManyDelays","opds",1600,1200);
    TCanvas * ManyBaselines = new TCanvas("ManyBaselines","baselines",1600,1200);
    ManyHeatMaps->Divide(nx,ny);
    ManyHeatMaps->Print("HBTMaps.pdf[");
    ManyFFTs->Divide(nx,ny);
    ManyFFTs->Print("FFTAll.pdf[");
    ManyDelays->Divide(nx,ny);
    ManyDelays->Print("delaysAll.pdf[");
    ManyBaselines->Divide(nx,ny);
    ManyBaselines->Print("baselinesAll.pdf[");
    
    //Save Heat Maps and FFTS to PDFs to be viewed
    for (int idir = 0; idir < pairNum; idir++){
      if (Setup(idir, DirName, HumanFlag, splitNum)==0) break;
      TFile * CorrInfo = new TFile(Form("%s/Analysis.root", DirName.Data()), "READONLY");
      CorrInfo->GetObject("HeatMap", HeatMap);
      CorrInfo->GetObject("FFTProjection2", FFT);
      
      CorrInfo->GetObject("BaselineGraph",Baseline); // josie added 7 july 23
      CorrInfo->GetObject("OPDGraph",Delays);
      //CorrInfo->GetObject("origBaselines",pythonBaselines);
      //CorrInfo->GetObject("origDelays",pythonDelays);
      //pythonDelays->SetMarkerStyle(9); pythonDelays->SetMarkerColor(4);   pythonBaselines->SetMarkerStyle(9);   pythonBaseli//nes->SetMarkerColor(4);
      Delays->SetMarkerStyle(9);    Delays->SetMarkerSize(2);    Baseline->SetMarkerStyle(9);   Baseline->SetMarkerSize(2);   
      
      nOnCanHeat++;
      ManyHeatMaps->cd(nOnCanHeat);
      HeatMap->SetTitle(Form("%d) %s",idir, DirName.Data()));
      HeatMap->Draw("COLZ");
      if (nOnCanHeat==nx*ny){
        ManyHeatMaps->Print("HBTMaps.pdf");  // note no parentheses nor bracket
        nOnCanHeat=0;
        ManyHeatMaps->Clear();    ManyHeatMaps->Divide(nx,ny);
      }
      
      nOnCanFFT++;
      ManyFFTs->cd(nOnCanFFT);
      FFT->SetTitle(Form("%d) %s",idir, DirName.Data()));
      FFT->Draw();
     if (nOnCanFFT==nx*ny){
        ManyFFTs->Print("FFTAll.pdf");  // note no parentheses nor bracket
        nOnCanFFT=0;
        ManyFFTs->Clear();    ManyFFTs->Divide(nx,ny);
      } 
      
      nOnCanDB++;
      ManyDelays->cd(nOnCanDB);
      Delays->SetTitle(Form("%d) %s",idir, DirName.Data()));
      Delays->Draw("AP");
      ManyBaselines->cd(nOnCanDB);
      Baseline->SetTitle(Form("%d) %s",idir, DirName.Data()));
      Baseline->Draw("AP");
      if (nOnCanDB==nx*ny){
          ManyDelays->Print("delaysAll.pdf");
          ManyBaselines->Print("baselinesAll.pdf");
          nOnCanDB=0;
          ManyDelays->Clear();  ManyDelays->Divide(nx,ny);
          ManyBaselines->Clear();   ManyBaselines->Divide(nx,ny);
      }
      
    }
    
    ManyHeatMaps->Print("HBTMaps.pdf)");   // note parenthesis
    ManyFFTs->Print("FFTAll.pdf)"); 
    ManyDelays->Print("delaysAll.pdf)");
    ManyBaselines->Print("baselinesAll.pdf)");

    //------------------------------------------------------------------------------------ 
    //-------------- Use this to determine if there are Type A or B glitches -------------
    //------------------------------------------------------------------------------------ 

    /*TCanvas * ManyNorms = new TCanvas("Many Norms", "Many Norms", 1600, 1200);
    ManyNorms->Print("TryNorm.pdf[");
    
    TProfile2D *Norm;
    for (int idir = 0; idir < pairNum; idir++){
      if (Setup(idir, DirName, HumanFlag, splitNum)==0) break;
      TFile * CorrInfo = new TFile(Form("%s/Analysis2.root", DirName.Data()), "READONLY");
      CorrInfo->GetObject("CFNormalized", Norm);

      ManyNorms->cd();
      Norm->Draw("LEGO");
      ManyNorms->Print("TryNorm.pdf");

      cout << idir << endl;
      CorrInfo->Close();
    }
    ManyNorms->Print("TryNorm.pdf)"); */
    
    //------------------------------------------------------------------------------------ 
    //-------------- Use this to determine if Normalization is Correct -------------------
    //------------------------------------------------------------------------------------
    
    /*TCanvas * CheckNorms = new TCanvas("Check Norms", "Check Norms", 1600, 1200);
    CheckNorms->Print("CheckNormAfterNoise.pdf[");
    
    TProfile2D * NormAfterNoise;
    for (int idir = 0; idir < pairNum; idir++){
      if (Setup(idir, DirName, HumanFlag, splitNum)==0) break;
      TFile * CorrInfo = new TFile(Form("%s/Analysis.root", DirName.Data()), "READONLY");
      CorrInfo->GetObject("CFMinusNoise1", NormAfterNoise);

      CheckNorms->cd();
      NormAfterNoise->Draw("LEGO");
      CheckNorms->Print("CheckNormAfterNoise.pdf");

      cout << idir << endl;
      CorrInfo->Close();
    }
    CheckNorms->Print("CheckNormAfterNoise.pdf)");*/

}    

TF1 * HBTFit(TProfile * HBT, double relTimePar){
    
    TF1* PeakFit = new TF1("HBTpeak","([0]*exp(-pow((x-[1])/[2],2)/2.0))/([2]*sqrt(2*pi))",-300,300); 
    PeakFit->SetParName(0,"Area");          PeakFit->SetParameter(0,0.0);     
    PeakFit->SetParName(1,"#tau_{0} (ns)"); PeakFit->SetParameter(1,relTimePar); PeakFit->SetParLimits(1,relTimePar-relTimeWindow, relTimePar+relTimeWindow); 
    PeakFit->SetParName(2,"#sigma (ns)");   PeakFit->SetParameter(2,4.5); PeakFit->SetParLimits(2, sigMin,  sigMax);
    TFitResultPtr fitResult = HBT->Fit(PeakFit, "SQ");
    return PeakFit;
}

TF1 * HBTFitAmp(TH1D * Noise, double relTimePar){
    
    TF1* PeakFit = new TF1("HBTpeak","([0]*exp(-pow((x-[1])/[2],2)/2.0))",-300,300); 
    PeakFit->SetParName(0,"Amp");           PeakFit->SetParameter(0,0.0);     
    PeakFit->SetParName(1,"#tau_{0} (ns)"); PeakFit->SetParameter(1,relTimePar); PeakFit->SetParLimits(1,relTimePar-relTimeWindow, relTimePar+relTimeWindow); 
    PeakFit->SetParName(2,"#sigma (ns)");   PeakFit->SetParameter(2,4.5); PeakFit->SetParLimits(2, sigMin,  sigMax);
    TFitResultPtr fitResult = Noise->Fit(PeakFit, "SQ");
    return PeakFit;
    
}

void DrawGraph(TGraph * tgtemp, TF1 *Bess, TH1D * UnifDisk, const char * title){
    TCanvas *TempCan = new TCanvas("temp","temp",1200,900);
    tgtemp->Fit(Bess, "EX0");
    UnifDisk->Draw();
    UnifDisk->SetTitle(title); 
    tgtemp->Draw("P, same");
    Bess->Draw("same");
    TempCan->Print("HBTStamps.pdf");
}


int RelativeTimeOffset(int tx, int ty){
    double relTimePar(-5);
    if((tx==1 && ty==2)){relTimePar = -20;} if((tx==2 && ty==1)){relTimePar =  20;}
    if((tx==1 && ty==3)){relTimePar = -10;} if((tx==3 && ty==1)){relTimePar =  10;}
    if((tx==1 && ty==4)){relTimePar =   2;} if((tx==4 && ty==1)){relTimePar = -2;}
    if((tx==2 && ty==3)){relTimePar =  10;} if((tx==3 && ty==2)){relTimePar = -10;}
    if((tx==2 && ty==4)){relTimePar =  18;} if((tx==4 && ty==2)){relTimePar = -18;}
    if((tx==3 && ty==4)){relTimePar =   8;} if((tx==4 && ty==3)){relTimePar = -8;}
    return relTimePar;
}


int newiteration(0);
double CalcSNRRatio(TProfile * HBT, TF1 *PeakFit){
  int NoiseDiff(20), NoiseWindow(40), nWindows((HBT->GetXaxis()->GetNbins()-160)/NoiseDiff);
  double AveBackgroundAmp(0);
  double PeakAmpReal  = PeakFit->GetParameter(0)/(sqrt(2*3.141592)*PeakFit->GetParameter(2));
  for(int NoiseValues = HBT->GetXaxis()->GetBinCenter(1)+30; NoiseValues < HBT->GetXaxis()->GetBinCenter(HBT->GetXaxis()->GetNbins())-30; NoiseValues += NoiseDiff ){   //ns 
    if(NoiseValues > -70 && NoiseValues < 70){continue;}  
    TH1D * NoiseChunk = new TH1D(Form("noisechunk%d%d",newiteration, NoiseValues), "Noise Hist", NoiseWindow/HBT->GetBinWidth(1), NoiseValues-0.5*NoiseWindow,  NoiseValues+0.5*NoiseWindow);
    for(int Noisebins = 1; Noisebins <= NoiseChunk->GetXaxis()->GetNbins(); ++Noisebins){
      NoiseChunk->SetBinContent(Noisebins, HBT->GetBinContent(HBT->FindBin(NoiseValues-0.5*NoiseWindow)+Noisebins-1));
      NoiseChunk->SetBinError(Noisebins,   HBT->GetBinError(HBT->FindBin(NoiseValues-0.5*NoiseWindow)+Noisebins-1));
    }
    TF1 * PeakFitAmp = HBTFitAmp(NoiseChunk, (double)NoiseValues);
    AveBackgroundAmp += pow(PeakFitAmp->GetParameter(0),2);
  }
  AveBackgroundAmp /= nWindows; 
  AveBackgroundAmp = sqrt(AveBackgroundAmp); 
  ++newiteration;
      
  //Signal to Noise Ratio for a given runpair
  return abs(PeakAmpReal/AveBackgroundAmp); 
}

int Setup(int irun, TString& DirName, int& HumanFlag, int& splitNum){
  splitNum = 1;  switch (irun) {
  case 0:
    DirName = "UTC20211217/m12d16h23/T1T2";
    HumanFlag = 0;
    break;
  case 1:
    DirName = "UTC20211217/m12d16h23/T1T3";
    HumanFlag = 0;
    break;
  case 2:
    DirName = "UTC20211217/m12d16h23/T1T4";
    HumanFlag = 0;
    break;
  case 3:
    DirName = "UTC20211217/m12d16h23/T2T3";
    HumanFlag = 0;
    break;
  case 4:
    DirName = "UTC20211217/m12d16h23/T2T4";
    HumanFlag = 0;
    break;
  case 5:
    DirName = "UTC20211217/m12d16h23/T3T4";
    HumanFlag = 0;
    break;
  case 6:
    DirName = "UTC20211217/m12d17h01/T1T2";
    HumanFlag = 0;
    break;
  case 7:
    DirName = "UTC20211217/m12d17h01/T1T3";
    HumanFlag = 0;
    break;
  case 8:
    DirName = "UTC20211217/m12d17h01/T1T4";
    HumanFlag = 0;
    break;
  case 9:
    DirName = "UTC20211217/m12d17h01/T2T3";
    HumanFlag = 0;
    break;
  case 10:
    DirName = "UTC20211217/m12d17h01/T2T4";
    HumanFlag = 0;
    break;
  case 11:
    DirName = "UTC20211217/m12d17h01/T3T4";
    HumanFlag = 0;
    break;
  case 12:
    DirName = "UTC20211217/m12d17h03/T1T2";
    HumanFlag = 0;
    break;
  case 13:
    DirName = "UTC20211217/m12d17h03/T1T3";
    HumanFlag = 0;
    break;
  case 14:
    DirName = "UTC20211217/m12d17h03/T1T4";
    HumanFlag = 0;
    break;
  case 15:
    DirName = "UTC20211217/m12d17h03/T2T3";
    HumanFlag = 0;
    break;
  case 16:
    DirName = "UTC20211217/m12d17h03/T2T4";
    HumanFlag = 0;
    break;
  case 17:
    DirName = "UTC20211217/m12d17h03/T3T4";
    HumanFlag = 0;
    break;
  case 18:
    DirName = "UTC20211217/m12d17h05/T1T2";
    HumanFlag = 0;
    break;
  case 19:
    DirName = "UTC20211217/m12d17h05/T1T3";
    HumanFlag = 0;
    break;
  case 20:
    DirName = "UTC20211217/m12d17h05/T1T4";
    HumanFlag = 0;
    break;
  case 21:
    DirName = "UTC20211217/m12d17h05/T2T3";
    HumanFlag = 0;
    break;
  case 22:
    DirName = "UTC20211217/m12d17h05/T2T4";
    HumanFlag = 0;
    break;
  case 23:
    DirName = "UTC20211217/m12d17h05/T3T4";
    HumanFlag = 0;
    break;
  case 24:
    DirName = "UTC20220212/m02d11h19/T1T2";
    HumanFlag = -1;
    break;
  case 25:
    DirName = "UTC20220212/m02d11h19/T1T3";
    HumanFlag = -1;
    break;
  case 26:
    DirName = "UTC20220212/m02d11h19/T1T4";
    HumanFlag = -1;
    break;
  case 27:
    DirName = "UTC20220212/m02d11h19/T2T3";
    HumanFlag = 0;
    break;
  case 28:
    DirName = "UTC20220212/m02d11h19/T2T4";
    HumanFlag = 0;
    break;
  case 29:
    DirName = "UTC20220212/m02d11h19/T3T4";
    HumanFlag = 0;
    break;
  case 30:
    DirName = "UTC20220212/m02d11h20/T1T2";
    HumanFlag = -1;
    break;
  case 31:
    DirName = "UTC20220212/m02d11h20/T1T3";
    HumanFlag = -1;
    break;
  case 32:
    DirName = "UTC20220212/m02d11h20/T1T4";
    HumanFlag = -1;
    break;
  case 33:
    DirName = "UTC20220212/m02d11h20/T2T3";
    HumanFlag = 0;
    break;
  case 34:
    DirName = "UTC20220212/m02d11h20/T2T4";
    HumanFlag = 0;
    break;
  case 35:
    DirName = "UTC20220212/m02d11h20/T3T4";
    HumanFlag = 0;
    break;
  case 36:
    DirName = "UTC20220212/m02d11h23/T1T2";
    HumanFlag = -2; // we decided this is an outlier and should not be included
    break;
  case 37:
    DirName = "UTC20220212/m02d11h23/T1T3";
    HumanFlag = 0;
    break;
  case 38:
    DirName = "UTC20220212/m02d11h23/T1T4";
    HumanFlag = 0;
    break;
  case 39:
    DirName = "UTC20220212/m02d11h23/T2T3";
    HumanFlag = 0;
    break;
  case 40:
    DirName = "UTC20220212/m02d11h23/T2T4";
    HumanFlag = 0;
    break;
  case 41:
    DirName = "UTC20220212/m02d11h23/T3T4";
    HumanFlag = 0;
    break;
  case 42:
    DirName = "UTC20220212/m02d12h01/T1T2";
    HumanFlag = 0;
    break;
  case 43:
    DirName = "UTC20220212/m02d12h01/T1T3";
    HumanFlag = 0;
    break;
  case 44:
    DirName = "UTC20220212/m02d12h01/T1T4";
    HumanFlag = 0;
    break;
  case 45:
    DirName = "UTC20220212/m02d12h01/T2T3";
    HumanFlag = 0;
    break;
  case 46:
    DirName = "UTC20220212/m02d12h01/T2T4";
    HumanFlag = 0;
    break;
  case 47:
    DirName = "UTC20220212/m02d12h01/T3T4";
    HumanFlag = 0;
    break;
  case 48:
    DirName = "UTC20220212/m02d12h03/T1T2";
    HumanFlag = 0;
    break;
  case 49:
    DirName = "UTC20220212/m02d12h03/T1T3";
    HumanFlag = 0;
    break;
  case 50:
    DirName = "UTC20220212/m02d12h03/T1T4";
    HumanFlag = 0;
    break;
  case 51:
    DirName = "UTC20220212/m02d12h03/T2T3";
    HumanFlag = 0;
    break;
  case 52:
    DirName = "UTC20220212/m02d12h03/T2T4";
    HumanFlag = 0;
    break;
  case 53:
    DirName = "UTC20220212/m02d12h03/T3T4";
    HumanFlag = 0;
    break;
  case 54:
    DirName = "UTC20220214/m02d13h19/T1T2";
    HumanFlag = 0;
    break;
  case 55:
    DirName = "UTC20220214/m02d13h19/T1T3";
    HumanFlag = 0;
    break;
  case 56:
    DirName = "UTC20220214/m02d13h19/T1T4";
    HumanFlag = 0;
    break;
  case 57:
    DirName = "UTC20220214/m02d13h19/T2T3";
    HumanFlag = 0;
    break;
  case 58:
    DirName = "UTC20220214/m02d13h19/T2T4";
    HumanFlag = 0;
    break;
  case 59:
    DirName = "UTC20220214/m02d13h19/T3T4";
    HumanFlag = 0;
    break;
  case 60:
    DirName = "UTC20220214/m02d13h21/T1T2";
    HumanFlag = 0;
    break;
  case 61:
    DirName = "UTC20220214/m02d13h21/T1T3";
    HumanFlag = 0;
    break;
  case 62:
    DirName = "UTC20220214/m02d13h21/T1T4";
    HumanFlag = 0;
    break;
  case 63:
    DirName = "UTC20220214/m02d13h21/T2T3";
    HumanFlag = 0;
    break;
  case 64:
    DirName = "UTC20220214/m02d13h21/T2T4";
    HumanFlag = 0;
    break;
  case 65:
    DirName = "UTC20220214/m02d13h21/T3T4";
    HumanFlag = 0;
    break;
  case 66:
    DirName = "UTC20220214/m02d13h23/T1T2";
    HumanFlag = 0;
    break;
  case 67:
    DirName = "UTC20220214/m02d13h23/T1T3";
    HumanFlag = 0;
    break;
  case 68:
    DirName = "UTC20220214/m02d13h23/T1T4";
    HumanFlag = 0;
    break;
  case 69:
    DirName = "UTC20220214/m02d13h23/T2T3";
    HumanFlag = 0;
    break;
  case 70:
    DirName = "UTC20220214/m02d13h23/T2T4";
    HumanFlag = 0;
    break;
  case 71:
    DirName = "UTC20220214/m02d13h23/T3T4";
    HumanFlag = 0;
    break;
  case 72:
    DirName = "UTC20220214/m02d14h01/T1T2";
    HumanFlag = 0;
    break;
  case 73:
    DirName = "UTC20220214/m02d14h01/T1T3";
    HumanFlag = 0;
    break;
  case 74:
    DirName = "UTC20220214/m02d14h01/T1T4";
    HumanFlag = 0;
    break;
  case 75:
    DirName = "UTC20220214/m02d14h01/T2T3";
    HumanFlag = 0;
    break;
  case 76:
    DirName = "UTC20220214/m02d14h01/T2T4";
    HumanFlag = 0;
    break;
  case 77:
    DirName = "UTC20220214/m02d14h01/T3T4";
    HumanFlag = 0;
    break;
  case 78:
    DirName = "UTC20220214/m02d14h03/T1T2";
    HumanFlag = 0;
    break;
  case 79:
    DirName = "UTC20220214/m02d14h03/T1T3";
    HumanFlag = 0;
    break;
  case 80:
    DirName = "UTC20220214/m02d14h03/T1T4";
    HumanFlag = 0;
    break;
  case 81:
    DirName = "UTC20220214/m02d14h03/T2T3";
    HumanFlag = 0;
    break;
  case 82:
    DirName = "UTC20220214/m02d14h03/T2T4";
    HumanFlag = 0;
    break;
  case 83:
    DirName = "UTC20220214/m02d14h03/T3T4";
    HumanFlag = 0;
    break;
  case 84:
    DirName = "UTC20220313/m03d12h22/T1T2";
    HumanFlag = 0;
    break;
  case 85:
    DirName = "UTC20220313/m03d12h22/T1T3";
    HumanFlag = 0;
    break;
  case 86:
    DirName = "UTC20220313/m03d12h22/T1T4";
    HumanFlag = 0;
    break;
  case 87:
    DirName = "UTC20220313/m03d12h22/T2T3";
    HumanFlag = 0;
    break;
  case 88:
    DirName = "UTC20220313/m03d12h22/T2T4";
    HumanFlag = 0;
    break;
  case 89:
    DirName = "UTC20220313/m03d12h22/T3T4";
    HumanFlag = 0;
    break;
  case 90:
    DirName = "UTC20220313/m03d13h01/T1T2";
    HumanFlag = 0;
    break;
  case 91:
    DirName = "UTC20220313/m03d13h01/T1T3";
    HumanFlag = 0;
    break;
  case 92:
    DirName = "UTC20220313/m03d13h01/T1T4";
    HumanFlag = 0;
    break;
  case 93:
    DirName = "UTC20220313/m03d13h01/T2T3";
    HumanFlag = 0;
    break;
  case 94:
    DirName = "UTC20220313/m03d13h01/T2T4";
    HumanFlag = 0;
    break;
  case 95:
    DirName = "UTC20220313/m03d13h01/T3T4";
    HumanFlag = 0;
    break;
  case 96:
    DirName = "UTC20220511/m05d10h21/T1T2";
    HumanFlag = 0;
    break;
  case 97:
    DirName = "UTC20220511/m05d10h21/T1T3";
    HumanFlag = 0;
    break;
  case 98:
    DirName = "UTC20220511/m05d10h21/T1T4";
    HumanFlag = 0;
    break;
  case 99:
    DirName = "UTC20220511/m05d10h21/T2T3";
    HumanFlag = 0;
    break;
   default:
     return 0;
     break;
   }
  return 1;
}


/*    SumAveVisClump1 /= nptsClump1;
    SumAveVisClump2 /= nptsClump2;
    SumAveVisClump3 /= nptsClump3;
    
    SumErrVisClump1 = sqrt(SumErrVisClump1/nptsClump1 - pow(SumAveVisClump1,2));
    SumErrVisClump2 = sqrt(SumErrVisClump2/nptsClump2 - pow(SumAveVisClump2,2));
    SumErrVisClump3 = sqrt(SumErrVisClump3/nptsClump3 - pow(SumAveVisClump3,2));
    
    SumErrVisClump1 /= sqrt(nptsClump1);
    SumErrVisClump2 /= sqrt(nptsClump2);
    SumErrVisClump3 /= sqrt(nptsClump3);
    
    
    SumAveBaseClump1 /= nptsClump1;
    SumAveBaseClump2 /= nptsClump2;
    SumAveBaseClump3 /= nptsClump3;
    
    SumErrBaseClump1 = sqrt(SumErrBaseClump1/nptsClump1 - pow(SumAveBaseClump1,2));
    SumErrBaseClump2 = sqrt(SumErrBaseClump2/nptsClump2 - pow(SumAveBaseClump2,2));
    SumErrBaseClump3 = sqrt(SumErrBaseClump3/nptsClump3 - pow(SumAveBaseClump3,2));
    
    SumErrBaseClump1 /= sqrt(nptsClump1);
    SumErrBaseClump2 /= sqrt(nptsClump2);
    SumErrBaseClump3 /= sqrt(nptsClump3);
    
    
    cout << "SumAveVisClump1:" << SumAveVisClump1 << "SumErrVisClump1:" << SumErrVisClump1 << "SumAveBaseClump1:" << SumAveBaseClump1 << "SumErrVisClump1:" << SumErrBaseClump1 << endl;
    cout << "SumAveVisClump2:" << SumAveVisClump2 << "SumErrVisClump2:" << SumErrVisClump2 << "SumAveBaseClump2:" << SumAveBaseClump2 << "SumErrVisClump2:" << SumErrBaseClump2 << endl;
    cout << "SumAveVisClump3:" << SumAveVisClump3 << "SumErrVisClump3:" << SumErrVisClump3 << "SumAveBaseClump3:" << SumAveBaseClump3 << "SumErrVisClump3:" << SumErrBaseClump3 << endl;
    
    
    //tgGood->AddPoint(SumAveBaseClump1, 0); tgGood->SetPointError(nptsGood, SumErrBaseClump1, RMSMuck);
    //tgGoodSingle->AddPoint(SumAveBaseClump1, 0); tgGoodSingle->SetPointError(0, SumErrBaseClump1, RMSMuck);
    
   // tgGood->AddPoint(SumAveBaseClump2, 0); tgGood->SetPointError(nptsGood, SumErrBaseClump2, RMSMuck); ++nptsGood;
    //tgGood->AddPoint(SumAveBaseClump3, 0); tgGood->SetPointError(nptsGood, SumErrBaseClump3, RMSMuck); ++nptsGood;*/






// oneClumpMax(85), oneClumpMin(105), twoClumpMax(105), twoClumpMin(150), threeClumpMax(150), threeClumpMin(300);

            //Create 3 large background clumps
           // if(BasePoint > 90  && BasePoint <= 105){if(BasePoint < oneClumpMin){oneClumpMin = BasePoint;} if(BasePoint > oneClumpMax){oneClumpMax = BasePoint;}}
            //if(BasePoint > 105 && BasePoint <= 150){if(BasePoint < twoClumpMin){twoClumpMin = BasePoint;} if(BasePoint > twoClumpMax){twoClumpMax = BasePoint;}}
            //if(BasePoint > 150){if(BasePoint < threeClumpMin){threeClumpMin = BasePoint;} if(BasePoint > threeClumpMax){threeClumpMax = BasePoint;}}


    //TLine * upperbound1 = new TLine(oneClumpMin,RMSMuck,oneClumpMax,RMSMuck);
    //TLine * upperbound2 = new TLine(twoClumpMin,RMSMuck,twoClumpMax,RMSMuck);
    //TLine * upperbound3 = new TLine(threeClumpMin,RMSMuck,threeClumpMax,RMSMuck);
    
    //TArrow * arrow1 = new TArrow((oneClumpMin+oneClumpMax)/2,RMSMuck,(oneClumpMin+oneClumpMax)/2,0, 0.02,"|>"); arrow1->SetAngle(30); 
    //TArrow * arrow2 = new TArrow((twoClumpMin+twoClumpMax)/2,RMSMuck,(twoClumpMin+twoClumpMax)/2,0,0.02, "|>"); arrow2->SetAngle(30);  
    //TArrow * arrow3 = new TArrow((threeClumpMin+threeClumpMax)/2,RMSMuck,(threeClumpMin+threeClumpMax)/2,0, 0.02, "|>"); arrow3->SetAngle(30);


    // ---------- Remove Points ----------------
    
   /* TCanvas *RemovePoints = new TCanvas("Remove","Remove",1200,900);
    UniformDiskModel->Draw();
    tgCutOff->Fit(BessFit);
    tgCutOff->Draw("P, same");
    upperbound1->Draw("Same");
    upperbound2->Draw("Same");
    upperbound3->Draw("Same");
    arrow1->Draw();
    arrow2->Draw();
    arrow3->Draw();
    BessFit->Draw("same"); 
    RemovePoints->Print("HBTStampsNew.pdf");*/
    
