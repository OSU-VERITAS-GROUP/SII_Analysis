// josie jan 29 2024
// 14 may 2024 - josie modifies "ugly" version to be 1) nicer and 2) add hour angle code from RoseAnalysis/hourAngles WITH data 
// 19 july 2024 - josie adding 10 mhz matching 
// 21 oct 2024 - josie adding a concrete threshold for throwing out low visibility points (new funct)

// fitting function using TMinuit
/*void uniformDiskFunc(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
        double norm = par[0];
        double angdiam = par[1];
        
        f = fabs(norm)*pow(2.0*TMath::BesselJ1(TMath::Pi()*x*angdiam*1e-3/(4157e-10*206265))/(TMath::Pi()*x*angdiam*1e-3/(4157e-10*206265)),2);
    }*/
    
// return visibility from the fitted curve
double calcVis(double baseline, double starsize, double cnorm){
    double vis = fabs(cnorm)*pow(2.0*TMath::BesselJ1(TMath::Pi()*baseline*starsize*1e-3/(4157e-10*206265))/(TMath::Pi()*baseline*starsize*1e-3/(4157e-10*206265)),2);
    return vis;
}

double calcThresh(TProfile* cf);

double calcExcessErr(TProfile* cf, double area, double sigma, int reps);

TH1D* simpeak(double area, double tau, double sigma, int windowsize);

void collectAllVisComplete_lowthresh(){
    
    TString source = "gam cas"; // only necessary for calculating hour angle with the python 
    
    int rebin(1);
    
    // in the nice future version there will be the subroutine to tell how many and what runs and which directories they live in - I am hard coding it
    const int nfiles(23);
    TString filename[nfiles] = {
        "2023DecDataVersii/gam cas_y2023m12d24h18m57s53_analysis.root",
        "2023DecDataVersii/gam cas_y2023m12d24h20m04s42_analysis.root",
        "2023DecDataVersii/gam cas_y2023m12d25h18m34s13_analysis.root",
        "2023DecDataVersii/gam cas_y2023m12d25h20m37s12_analysis.root",
        "2023DecDataVersii/gam cas_y2023m12d25h22m38s38_analysis.root",
        "2023DecDataVersii/gam cas_y2023m12d26h00m39s31_analysis.root",
        "2023DecDataVersii/gam Cas_y2023m12d26h18m43s56_analysis.root",
        "2023DecDataVersii/gam Cas_y2023m12d26h20m44s43_analysis.root",
        "2023DecDataVersii/gam Cas_y2023m12d26h22m47s50_analysis.root",
        "2023DataVersii/gam cas_y2023m01d07h22m23s24_analysis.root", 
        "2023DataVersii/gam Cas_y2023m02d01h20m10s16_analysis.root",
        "2023DataVersii/gam Cas_y2023m02d01h22m12s31_analysis.root",
        "2024FebDataVersii/gam cas_y2024m02d18h19m31s45_analysis.root",
        "2024FebDataVersii/gam cas_y2024m02d18h20m38s27_analysis.root",
        "2024FebDataVersii/gam cas_y2024m02d18h21m44s15_analysis.root",
        "2024FebDataVersii/gam cas_y2024m02d19h19m07s37_analysis.root",
        "2024FebDataVersii/gam cas_y2024m02d19h20m14s25_analysis.root",
        "2024FebDataVersii/gam cas_y2024m02d20h19m22s00_analysis.root",
        "2024FebDataVersii/gam cas_y2024m02d21h20m26s13_analysis.root",
        "2024FebDataVersii/gam cas_y2024m02d21h20m59s31_analysis.root",
        "2024FebDataVersii/gam cas_y2024m02d21h21m32s24_analysis.root",
        "2024MayDataVersii/gam cas_y2024m05d18h03m00s58_analysis.root",
        //"2024MayDataVersii/gam cas_y2024m05d18h04m03s24_analysis.root", // sunny run
        "2024MayDataVersii/gam cas_y2024m05d19h03m07s29_analysis.root"
    };
    
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    
    const double pi = TMath::Pi();
    
    int jcol[6] = {2, 4, 92, 3, 6, 1}; // jason's color palette
    
    // =============================================== set up canvases and plots ===========================================================
    
    // canvas and plots of visibility curve
    TCanvas* cvis = new TCanvas;
    TGraphAsymmErrors* allpts = new TGraphAsymmErrors();      allpts->SetMarkerStyle(20);     allpts->SetTitle("visibility curve;baseline (m);g^{2}(#tau)-1");
    TGraph* allptsLow = new TGraph;     allptsLow->SetMarkerStyle(20);      allptsLow->SetMarkerColor(16);
    TGraphAsymmErrors* viscurve1 = new TGraphAsymmErrors();   viscurve1->SetMarkerStyle(20);  //viscurve1->SetMarkerColor(
    TGraphAsymmErrors* viscurve2 = new TGraphAsymmErrors();   viscurve2->SetMarkerStyle(20);  viscurve2->SetMarkerColor(2);
    TGraphAsymmErrors* viscurve3 = new TGraphAsymmErrors();   viscurve3->SetMarkerStyle(20);  viscurve3->SetMarkerColor(3);
    TGraphAsymmErrors* viscurve4 = new TGraphAsymmErrors();   viscurve4->SetMarkerStyle(20);  viscurve4->SetMarkerColor(4);
    TGraphAsymmErrors* viscurveFeb2024 = new TGraphAsymmErrors();       viscurveFeb2024->SetMarkerStyle(20);        viscurveFeb2024->SetMarkerColor(2);     viscurveFeb2024->SetTitle("Feb 2024 data only;baseline (m);area from fit (ns)");
    TGraphAsymmErrors* viscurveDec2023 = new TGraphAsymmErrors();       viscurveDec2023->SetMarkerStyle(20);                                                viscurveDec2023->SetTitle("Dec 2023 data only;baseline (m);area from fit (ns)");
    TGraphAsymmErrors* viscurveJanFeb2023 = new TGraphAsymmErrors();    viscurveJanFeb2023->SetMarkerStyle(20);     viscurveJanFeb2023->SetMarkerColor(4);  viscurveJanFeb2023->SetTitle("Jan/Feb 2023 data only;baseline (m);area from fit (ns)");
    TGraphAsymmErrors* viscurveMay2024 = new TGraphAsymmErrors();       viscurveMay2024->SetMarkerStyle(20);        viscurveMay2024->SetMarkerColor(8);     viscurveMay2024->SetTitle("May 2024 data only;baseline (m);area from fit (ns)");
    TGraphAsymmErrors* viscurveFebMay2024 = new TGraphAsymmErrors();    viscurveFebMay2024->SetMarkerStyle(20);     viscurveFebMay2024->SetMarkerColor(6);  viscurveFebMay2024->SetTitle("Feb and May 2024 data;baseline (m);area from fit (ns)");
    TGraphAsymmErrors* viscurveNoJan = new TGraphAsymmErrors();         viscurveNoJan->SetMarkerStyle(20);          viscurveNoJan->SetMarkerColor(7);       viscurveNoJan->SetTitle("data since newest filters only;baseline (m);area from fit (ns)");
    TGraphAsymmErrors* viscurveNoFeb = new TGraphAsymmErrors();         viscurveNoFeb->SetMarkerStyle(20);          viscurveNoFeb->SetMarkerColor(94);      viscurveNoFeb->SetTitle("excluding Feb 2024;baseline (m);area from fit (ns)");
    TGraphAsymmErrors* viscurveDecMay = new TGraphAsymmErrors();        viscurveDecMay->SetMarkerStyle(20);         viscurveDecMay->SetMarkerColor(86);     viscurveDecMay->SetTitle("Dec 2023 and May 2024 only;baseline (m);area from fit (ns)");
    // these are by pair for each observation group (yes this is a ton of plots)
    TGraphAsymmErrors* vcjf23T1T2 = new TGraphAsymmErrors;        vcjf23T1T2->SetMarkerStyle(20);     vcjf23T1T2->SetMarkerColor(jcol[0]);  // 2, 8, 92, 6, 7, 4
    TGraphAsymmErrors* vcjf23T1T3 = new TGraphAsymmErrors;        vcjf23T1T3->SetMarkerStyle(20);     vcjf23T1T3->SetMarkerColor(jcol[1]);
    TGraphAsymmErrors* vcjf23T1T4 = new TGraphAsymmErrors;        vcjf23T1T4->SetMarkerStyle(20);     vcjf23T1T4->SetMarkerColor(jcol[2]);
    TGraphAsymmErrors* vcjf23T2T3 = new TGraphAsymmErrors;        vcjf23T2T3->SetMarkerStyle(20);     vcjf23T2T3->SetMarkerColor(jcol[3]);
    TGraphAsymmErrors* vcjf23T2T4 = new TGraphAsymmErrors;        vcjf23T2T4->SetMarkerStyle(20);     vcjf23T2T4->SetMarkerColor(jcol[4]);
    TGraphAsymmErrors* vcjf23T3T4 = new TGraphAsymmErrors;        vcjf23T3T4->SetMarkerStyle(20);     vcjf23T3T4->SetMarkerColor(jcol[5]);
    TGraphAsymmErrors* vcdec23T1T2 = new TGraphAsymmErrors;       vcdec23T1T2->SetMarkerStyle(20);    vcdec23T1T2->SetMarkerColor(2);
    TGraphAsymmErrors* vcdec23T1T3 = new TGraphAsymmErrors;       vcdec23T1T3->SetMarkerStyle(20);    vcdec23T1T3->SetMarkerColor(8);
    TGraphAsymmErrors* vcdec23T1T4 = new TGraphAsymmErrors;       vcdec23T1T4->SetMarkerStyle(20);    vcdec23T1T4->SetMarkerColor(92);
    TGraphAsymmErrors* vcdec23T2T3 = new TGraphAsymmErrors;       vcdec23T2T3->SetMarkerStyle(20);    vcdec23T2T3->SetMarkerColor(6);
    TGraphAsymmErrors* vcdec23T2T4 = new TGraphAsymmErrors;       vcdec23T2T4->SetMarkerStyle(20);    vcdec23T2T4->SetMarkerColor(7);
    TGraphAsymmErrors* vcdec23T3T4 = new TGraphAsymmErrors;       vcdec23T3T4->SetMarkerStyle(20);    vcdec23T3T4->SetMarkerColor(4);
    TGraphAsymmErrors* vcfeb24T1T2 = new TGraphAsymmErrors;       vcfeb24T1T2->SetMarkerStyle(20);    vcfeb24T1T2->SetMarkerColor(2);
    TGraphAsymmErrors* vcfeb24T1T3 = new TGraphAsymmErrors;       vcfeb24T1T3->SetMarkerStyle(20);    vcfeb24T1T3->SetMarkerColor(8);
    TGraphAsymmErrors* vcfeb24T1T4 = new TGraphAsymmErrors;       vcfeb24T1T4->SetMarkerStyle(20);    vcfeb24T1T4->SetMarkerColor(92);
    TGraphAsymmErrors* vcfeb24T2T3 = new TGraphAsymmErrors;       vcfeb24T2T3->SetMarkerStyle(20);    vcfeb24T2T3->SetMarkerColor(6);
    TGraphAsymmErrors* vcfeb24T2T4 = new TGraphAsymmErrors;       vcfeb24T2T4->SetMarkerStyle(20);    vcfeb24T2T4->SetMarkerColor(7);
    TGraphAsymmErrors* vcfeb24T3T4 = new TGraphAsymmErrors;       vcfeb24T3T4->SetMarkerStyle(20);    vcfeb24T3T4->SetMarkerColor(4);
    TGraphAsymmErrors* vcmay24T1T2 = new TGraphAsymmErrors;       vcmay24T1T2->SetMarkerStyle(20);    vcmay24T1T2->SetMarkerColor(2);
    TGraphAsymmErrors* vcmay24T1T3 = new TGraphAsymmErrors;       vcmay24T1T3->SetMarkerStyle(20);    vcmay24T1T3->SetMarkerColor(8);
    TGraphAsymmErrors* vcmay24T1T4 = new TGraphAsymmErrors;       vcmay24T1T4->SetMarkerStyle(20);    vcmay24T1T4->SetMarkerColor(92);
    TGraphAsymmErrors* vcmay24T2T3 = new TGraphAsymmErrors;       vcmay24T2T3->SetMarkerStyle(20);    vcmay24T2T3->SetMarkerColor(6);
    TGraphAsymmErrors* vcmay24T2T4 = new TGraphAsymmErrors;       vcmay24T2T4->SetMarkerStyle(20);    vcmay24T2T4->SetMarkerColor(7);
    TGraphAsymmErrors* vcmay24T3T4 = new TGraphAsymmErrors;       vcmay24T3T4->SetMarkerStyle(20);    vcmay24T3T4->SetMarkerColor(4);
    TGraph* axes = new TGraph();    axes->SetMarkerStyle(9);    axes->SetMarkerColor(0);
    axes->AddPoint(0,-2e-6);    axes->AddPoint(180,15e-6);
    int npts(0), npts1(0), npts2(0), npts3(0), npts4(0), npts2024(0), nptsd2023(0), nptsjf2023(0), nptsm2024(0), nptsfm2024(0), nptsnj(0), nptsnf(0), nptsdm(0);
    int npts12jf23(0), npts13jf23(0), npts14jf23(0), npts23jf23(0), npts24jf23(0), npts34jf23(0), npts12d23(0), npts13d23(0), npts14d23(0), npts23d23(0), npts24d23(0), npts34d23(0);
    int npts12f24(0), npts13f24(0), npts14f24(0), npts23f24(0), npts24f24(0), npts34f24(0), npts12m24(0), npts13m24(0), npts14m24(0), npts23m24(0), npts24m24(0), npts34m24(0);
    cvis->Print("viscurve.pdf[");
    
    // canvas and uv angle plots for each pair, each observation period
    TCanvas* chrangle = new TCanvas("","",1400,600);
    chrangle->Divide(3,2);
    TGraphAsymmErrors* hrangleT1T2 = new TGraphAsymmErrors();     hrangleT1T2->SetMarkerStyle(20);    hrangleT1T2->SetTitle("T1T2;hour angle;visibility");    hrangleT1T2->SetMarkerColor(2);  // 2,4,5,3,6,1
    TGraphAsymmErrors* hrangleT1T3 = new TGraphAsymmErrors();     hrangleT1T3->SetMarkerStyle(20);    hrangleT1T3->SetTitle("T1T3;hour angle;visibility");    hrangleT1T3->SetMarkerColor(3);
    TGraphAsymmErrors* hrangleT1T4 = new TGraphAsymmErrors();     hrangleT1T4->SetMarkerStyle(20);    hrangleT1T4->SetTitle("T1T4;hour angle;visibility");    hrangleT1T4->SetMarkerColor(92);
    TGraphAsymmErrors* hrangleT2T3 = new TGraphAsymmErrors();     hrangleT2T3->SetMarkerStyle(20);    hrangleT2T3->SetTitle("T2T3;hour angle;visibility");    hrangleT2T3->SetMarkerColor(7);
    TGraphAsymmErrors* hrangleT2T4 = new TGraphAsymmErrors();     hrangleT2T4->SetMarkerStyle(20);    hrangleT2T4->SetTitle("T2T4;hour angle;visibility");    hrangleT2T4->SetMarkerColor(6);
    TGraphAsymmErrors* hrangleT3T4 = new TGraphAsymmErrors();     hrangleT3T4->SetMarkerStyle(20);    hrangleT3T4->SetTitle("T3T4;hour angle;visibility");    hrangleT3T4->SetMarkerColor(4);
    int nptsha12(0), nptsha13(0), nptsha14(0), nptsha23(0), nptsha24(0), nptsha34(0);
    TGraphErrors* ha12jf23 = new TGraphErrors();        ha12jf23->SetMarkerStyle(20);       ha12jf23->SetTitle("T1T2;hour angle;visibility");       ha12jf23->SetMarkerColor(jcol[0]); // was all 4
    TGraphErrors* ha13jf23 = new TGraphErrors();        ha13jf23->SetMarkerStyle(20);       ha13jf23->SetTitle("T1T3;hour angle;visibility");       ha13jf23->SetMarkerColor(jcol[1]);
    TGraphErrors* ha14jf23 = new TGraphErrors();        ha14jf23->SetMarkerStyle(20);       ha14jf23->SetTitle("T1T4;hour angle;visibility");       ha14jf23->SetMarkerColor(jcol[2]);
    TGraphErrors* ha23jf23 = new TGraphErrors();        ha23jf23->SetMarkerStyle(20);       ha23jf23->SetTitle("T2T3;hour angle;visibility");       ha23jf23->SetMarkerColor(jcol[3]);
    TGraphErrors* ha24jf23 = new TGraphErrors();        ha24jf23->SetMarkerStyle(20);       ha24jf23->SetTitle("T2T4;hour angle;visibility");       ha24jf23->SetMarkerColor(jcol[4]);
    TGraphErrors* ha34jf23 = new TGraphErrors();        ha34jf23->SetMarkerStyle(20);       ha34jf23->SetTitle("T3T4;hour angle;visibility");       ha34jf23->SetMarkerColor(jcol[5]);
    int nha12jf23(0), nha13jf23(0), nha14jf23(0), nha23jf23(0), nha24jf23(0), nha34jf23(0);
    TGraphErrors* ha12dec23 = new TGraphErrors();       ha12dec23->SetMarkerStyle(20);      ha12dec23->SetTitle("T1T2;hour angle;visibility");      ha12dec23->SetMarkerColor(jcol[0]); // was 1,3,92,7,6,4
    TGraphErrors* ha13dec23 = new TGraphErrors();       ha13dec23->SetMarkerStyle(20);      ha13dec23->SetTitle("T1T3;hour angle;visibility");      ha13dec23->SetMarkerColor(jcol[1]);
    TGraphErrors* ha14dec23 = new TGraphErrors();       ha14dec23->SetMarkerStyle(20);      ha14dec23->SetTitle("T1T4;hour angle;visibility");      ha14dec23->SetMarkerColor(jcol[2]);
    TGraphErrors* ha23dec23 = new TGraphErrors();       ha23dec23->SetMarkerStyle(20);      ha23dec23->SetTitle("T2T3;hour angle;visibility");      ha23dec23->SetMarkerColor(jcol[3]);
    TGraphErrors* ha24dec23 = new TGraphErrors();       ha24dec23->SetMarkerStyle(20);      ha24dec23->SetTitle("T2T4;hour angle;visibility");      ha24dec23->SetMarkerColor(jcol[4]);
    TGraphErrors* ha34dec23 = new TGraphErrors();       ha34dec23->SetMarkerStyle(20);      ha34dec23->SetTitle("T3T4;hour angle;visibility");      ha34dec23->SetMarkerColor(jcol[5]);
    int nha12dec23(0), nha13dec23(0), nha14dec23(0), nha23dec23(0), nha24dec23(0), nha34dec23(0);
    TGraphErrors* ha12feb24 = new TGraphErrors();       ha12feb24->SetMarkerStyle(20);      ha12feb24->SetTitle("T1T2;hour angle;visibility");      ha12feb24->SetMarkerColor(jcol[0]);
    TGraphErrors* ha13feb24 = new TGraphErrors();       ha13feb24->SetMarkerStyle(20);      ha13feb24->SetTitle("T1T3;hour angle;visibility");      ha13feb24->SetMarkerColor(jcol[1]);
    TGraphErrors* ha14feb24 = new TGraphErrors();       ha14feb24->SetMarkerStyle(20);      ha14feb24->SetTitle("T1T4;hour angle;visibility");      ha14feb24->SetMarkerColor(jcol[2]);
    TGraphErrors* ha23feb24 = new TGraphErrors();       ha23feb24->SetMarkerStyle(20);      ha23feb24->SetTitle("T2T3;hour angle;visibility");      ha23feb24->SetMarkerColor(jcol[3]);
    TGraphErrors* ha24feb24 = new TGraphErrors();       ha24feb24->SetMarkerStyle(20);      ha24feb24->SetTitle("T2T4;hour angle;visibility");      ha24feb24->SetMarkerColor(jcol[4]);
    TGraphErrors* ha34feb24 = new TGraphErrors();       ha34feb24->SetMarkerStyle(20);      ha34feb24->SetTitle("T3T4;hour angle;visibility");      ha34feb24->SetMarkerColor(jcol[5]);
    int nha12feb24(0), nha13feb24(0), nha14feb24(0), nha23feb24(0), nha24feb24(0), nha34feb24(0);
    TGraphErrors* ha12may24 = new TGraphErrors();       ha12may24->SetMarkerStyle(20);      ha12may24->SetTitle("T1T2;hour angle;visibility");      ha12may24->SetMarkerColor(jcol[0]);
    TGraphErrors* ha13may24 = new TGraphErrors();       ha13may24->SetMarkerStyle(20);      ha13may24->SetTitle("T1T3;hour angle;visibility");      ha13may24->SetMarkerColor(jcol[1]);
    TGraphErrors* ha14may24 = new TGraphErrors();       ha14may24->SetMarkerStyle(20);      ha14may24->SetTitle("T1T4;hour angle;visibility");      ha14may24->SetMarkerColor(jcol[2]);
    TGraphErrors* ha23may24 = new TGraphErrors();       ha23may24->SetMarkerStyle(20);      ha23may24->SetTitle("T2T3;hour angle;visibility");      ha23may24->SetMarkerColor(jcol[3]);
    TGraphErrors* ha24may24 = new TGraphErrors();       ha24may24->SetMarkerStyle(20);      ha24may24->SetTitle("T2T4;hour angle;visibility");      ha24may24->SetMarkerColor(jcol[4]);
    TGraphErrors* ha34may24 = new TGraphErrors();       ha34may24->SetMarkerStyle(20);      ha34may24->SetTitle("T3T4;hour angle;visibility");      ha34may24->SetMarkerColor(jcol[5]);
    int nha12may24(0), nha13may24(0), nha14may24(0), nha23may24(0), nha24may24(0), nha34may24(0);
    TGraphAsymmErrors* haT1T2norm = new TGraphAsymmErrors;      haT1T2norm->SetMarkerStyle(20);    haT1T2norm->SetTitle("T1T2;hour angle;visibility");    haT1T2norm->SetMarkerColor(2);  // was my palette: 2,3,92,7,6,4
    TGraphAsymmErrors* haT1T3norm = new TGraphAsymmErrors;      haT1T3norm->SetMarkerStyle(20);    haT1T3norm->SetTitle("T1T3;hour angle;visibility");    haT1T3norm->SetMarkerColor(3);
    TGraphAsymmErrors* haT1T4norm = new TGraphAsymmErrors;      haT1T4norm->SetMarkerStyle(20);    haT1T4norm->SetTitle("T1T4;hour angle;visibility");    haT1T4norm->SetMarkerColor(92);
    TGraphAsymmErrors* haT2T3norm = new TGraphAsymmErrors;      haT2T3norm->SetMarkerStyle(20);    haT2T3norm->SetTitle("T2T3;hour angle;visibility");    haT2T3norm->SetMarkerColor(7);
    TGraphAsymmErrors* haT2T4norm = new TGraphAsymmErrors;      haT2T4norm->SetMarkerStyle(20);    haT2T4norm->SetTitle("T2T4;hour angle;visibility");    haT2T4norm->SetMarkerColor(6);
    TGraphAsymmErrors* haT3T4norm = new TGraphAsymmErrors;      haT3T4norm->SetMarkerStyle(20);    haT3T4norm->SetTitle("T3T4;hour angle;visibility");    haT3T4norm->SetMarkerColor(4);
    chrangle->Print("hourAngles.pdf[");
    
    // plots to hold predicted visibility vs hour angle, for given star model
    TGraph* visHA12all = new TGraph;    visHA12all->SetTitle("T1T2;hour angle (hrs);visibility");       visHA12all->SetMarkerStyle(20);     visHA12all->SetMarkerColor(2);      visHA12all->SetMarkerSize(0.5);     visHA12all->SetLineStyle(1);    visHA12all->SetLineColor(2);
    TGraph* visHA13all = new TGraph;    visHA13all->SetTitle("T1T3;hour angle (hrs);visibility");       visHA13all->SetMarkerStyle(20);     visHA13all->SetMarkerColor(3);      visHA13all->SetMarkerSize(0.5);     visHA13all->SetLineStyle(1);    visHA13all->SetLineColor(3);
    TGraph* visHA14all = new TGraph;    visHA14all->SetTitle("T1T4;hour angle (hrs);visibility");       visHA14all->SetMarkerStyle(20);     visHA14all->SetMarkerColor(92);     visHA14all->SetMarkerSize(0.5);     visHA14all->SetLineStyle(1);    visHA14all->SetLineColor(92);
    TGraph* visHA23all = new TGraph;    visHA23all->SetTitle("T2T3;hour angle (hrs);visibility");       visHA23all->SetMarkerStyle(20);     visHA23all->SetMarkerColor(7);      visHA23all->SetMarkerSize(0.5);     visHA23all->SetLineStyle(1);    visHA23all->SetLineColor(7);    
    TGraph* visHA24all = new TGraph;    visHA24all->SetTitle("T2T4;hour angle (hrs);visibility");       visHA24all->SetMarkerStyle(20);     visHA24all->SetMarkerColor(6);      visHA24all->SetMarkerSize(0.5);     visHA24all->SetLineStyle(1);    visHA24all->SetLineColor(6);
    TGraph* visHA34all = new TGraph;    visHA34all->SetTitle("T3T4;hour angle (hrs);visibility");       visHA34all->SetMarkerStyle(20);     visHA34all->SetMarkerColor(4);      visHA34all->SetMarkerSize(0.5);     visHA34all->SetLineStyle(1);    visHA34all->SetLineColor(4);
    
    TGraph* visHA12jf23 = new TGraph;   visHA12jf23->SetTitle("T1T2;hour angle (hrs);visibility");      visHA12jf23->SetMarkerStyle(20);    visHA12jf23->SetMarkerColor(15);     visHA12jf23->SetMarkerSize(0.5);    visHA12jf23->SetLineStyle(2);   visHA12jf23->SetLineColor(2);
    TGraph* visHA13jf23 = new TGraph;   visHA13jf23->SetTitle("T1T3;hour angle (hrs);visibility");      visHA13jf23->SetMarkerStyle(20);    visHA13jf23->SetMarkerColor(15);     visHA13jf23->SetMarkerSize(0.5);    visHA13jf23->SetLineStyle(2);   visHA13jf23->SetLineColor(3);
    TGraph* visHA14jf23 = new TGraph;   visHA14jf23->SetTitle("T1T4;hour angle (hrs);visibility");      visHA14jf23->SetMarkerStyle(20);    visHA14jf23->SetMarkerColor(15);     visHA14jf23->SetMarkerSize(0.5);    visHA14jf23->SetLineStyle(2);   visHA14jf23->SetLineColor(92);
    TGraph* visHA23jf23 = new TGraph;   visHA23jf23->SetTitle("T2T3;hour angle (hrs);visibility");      visHA23jf23->SetMarkerStyle(20);    visHA23jf23->SetMarkerColor(15);     visHA23jf23->SetMarkerSize(0.5);    visHA23jf23->SetLineStyle(2);   visHA23jf23->SetLineColor(7);
    TGraph* visHA24jf23 = new TGraph;   visHA24jf23->SetTitle("T2T4;hour angle (hrs);visibility");      visHA24jf23->SetMarkerStyle(20);    visHA24jf23->SetMarkerColor(15);     visHA24jf23->SetMarkerSize(0.5);    visHA24jf23->SetLineStyle(2);   visHA24jf23->SetLineColor(6);
    TGraph* visHA34jf23 = new TGraph;   visHA34jf23->SetTitle("T3T4;hour angle (hrs);visibility");      visHA34jf23->SetMarkerStyle(20);    visHA34jf23->SetMarkerColor(15);     visHA34jf23->SetMarkerSize(0.5);    visHA34jf23->SetLineStyle(2);   visHA34jf23->SetLineColor(4);
    
    TGraph* visHA12dec23 = new TGraph;  visHA12dec23->SetTitle("T1T2;hour angle (hrs);visibility");     visHA12dec23->SetMarkerStyle(20);   visHA12dec23->SetMarkerColor(15);    visHA12dec23->SetMarkerSize(0.5);   visHA12dec23->SetLineStyle(3);  visHA12dec23->SetLineColor(2);
    TGraph* visHA13dec23 = new TGraph;  visHA13dec23->SetTitle("T1T3;hour angle (hrs);visibility");     visHA13dec23->SetMarkerStyle(20);   visHA13dec23->SetMarkerColor(15);    visHA13dec23->SetMarkerSize(0.5);   visHA13dec23->SetLineStyle(3);  visHA13dec23->SetLineColor(3);
    TGraph* visHA14dec23 = new TGraph;  visHA14dec23->SetTitle("T1T4;hour angle (hrs);visibility");     visHA14dec23->SetMarkerStyle(20);   visHA14dec23->SetMarkerColor(15);    visHA14dec23->SetMarkerSize(0.5);   visHA14dec23->SetLineStyle(3);  visHA14dec23->SetLineColor(92);
    TGraph* visHA23dec23 = new TGraph;  visHA23dec23->SetTitle("T2T3;hour angle (hrs);visibility");     visHA23dec23->SetMarkerStyle(20);   visHA23dec23->SetMarkerColor(15);    visHA23dec23->SetMarkerSize(0.5);   visHA23dec23->SetLineStyle(3);  visHA23dec23->SetLineColor(7);
    TGraph* visHA24dec23 = new TGraph;  visHA24dec23->SetTitle("T2T4;hour angle (hrs);visibility");     visHA24dec23->SetMarkerStyle(20);   visHA24dec23->SetMarkerColor(15);    visHA24dec23->SetMarkerSize(0.5);   visHA24dec23->SetLineStyle(3);  visHA24dec23->SetLineColor(6);
    TGraph* visHA34dec23 = new TGraph;  visHA34dec23->SetTitle("T3T4;hour angle (hrs);visibility");     visHA34dec23->SetMarkerStyle(20);   visHA34dec23->SetMarkerColor(15);    visHA34dec23->SetMarkerSize(0.5);   visHA34dec23->SetLineStyle(3);  visHA34dec23->SetLineColor(4);
    
    TGraph* visHA12feb24 = new TGraph;  visHA12feb24->SetTitle("T1T2;hour angle (hrs);visibility");     visHA12feb24->SetMarkerStyle(20);   visHA12feb24->SetMarkerColor(15);    visHA12feb24->SetMarkerSize(0.5);   visHA12feb24->SetLineStyle(7);  visHA12feb24->SetLineColor(2);
    TGraph* visHA13feb24 = new TGraph;  visHA13feb24->SetTitle("T1T3;hour angle (hrs);visibility");     visHA13feb24->SetMarkerStyle(20);   visHA13feb24->SetMarkerColor(15);    visHA13feb24->SetMarkerSize(0.5);   visHA13feb24->SetLineStyle(7);  visHA13feb24->SetLineColor(3);
    TGraph* visHA14feb24 = new TGraph;  visHA14feb24->SetTitle("T1T4;hour angle (hrs);visibility");     visHA14feb24->SetMarkerStyle(20);   visHA14feb24->SetMarkerColor(15);   visHA14feb24->SetMarkerSize(0.5);   visHA14feb24->SetLineStyle(7);  visHA14feb24->SetLineColor(92);
    TGraph* visHA23feb24 = new TGraph;  visHA23feb24->SetTitle("T2T3;hour angle (hrs);visibility");     visHA23feb24->SetMarkerStyle(20);   visHA23feb24->SetMarkerColor(15);    visHA23feb24->SetMarkerSize(0.5);   visHA23feb24->SetLineStyle(7);  visHA23feb24->SetLineColor(7);
    TGraph* visHA24feb24 = new TGraph;  visHA24feb24->SetTitle("T2T4;hour angle (hrs);visibility");     visHA24feb24->SetMarkerStyle(20);   visHA24feb24->SetMarkerColor(15);    visHA24feb24->SetMarkerSize(0.5);   visHA24feb24->SetLineStyle(7);  visHA24feb24->SetLineColor(6);
    TGraph* visHA34feb24 = new TGraph;  visHA34feb24->SetTitle("T3T4;hour angle (hrs);visibility");     visHA34feb24->SetMarkerStyle(20);   visHA34feb24->SetMarkerColor(15);    visHA34feb24->SetMarkerSize(0.5);   visHA34feb24->SetLineStyle(7);  visHA34feb24->SetLineColor(4);
    
    TGraph* visHA12may24 = new TGraph;  visHA12may24->SetTitle("T1T2;hour angle (hrs);visibility");     visHA12may24->SetMarkerStyle(20);   visHA12may24->SetMarkerColor(15);    visHA12may24->SetMarkerSize(0.5);   visHA12may24->SetLineStyle(7);  visHA12may24->SetLineColor(2);
    TGraph* visHA13may24 = new TGraph;  visHA13may24->SetTitle("T1T3;hour angle (hrs);visibility");     visHA13may24->SetMarkerStyle(20);   visHA13may24->SetMarkerColor(15);    visHA13may24->SetMarkerSize(0.5);   visHA13may24->SetLineStyle(7);  visHA13may24->SetLineColor(3);
    TGraph* visHA14may24 = new TGraph;  visHA14may24->SetTitle("T1T4;hour angle (hrs);visibility");     visHA14may24->SetMarkerStyle(20);   visHA14may24->SetMarkerColor(15);   visHA14may24->SetMarkerSize(0.5);   visHA14may24->SetLineStyle(7);  visHA14may24->SetLineColor(92);
    TGraph* visHA23may24 = new TGraph;  visHA23may24->SetTitle("T2T3;hour angle (hrs);visibility");     visHA23may24->SetMarkerStyle(20);   visHA23may24->SetMarkerColor(15);    visHA23may24->SetMarkerSize(0.5);   visHA23may24->SetLineStyle(7);  visHA23may24->SetLineColor(7);
    TGraph* visHA24may24 = new TGraph;  visHA24may24->SetTitle("T2T4;hour angle (hrs);visibility");     visHA24may24->SetMarkerStyle(20);   visHA24may24->SetMarkerColor(15);    visHA24may24->SetMarkerSize(0.5);   visHA24may24->SetLineStyle(7);  visHA24may24->SetLineColor(6);
    TGraph* visHA34may24 = new TGraph;  visHA34may24->SetTitle("T3T4;hour angle (hrs);visibility");     visHA34may24->SetMarkerStyle(20);   visHA34may24->SetMarkerColor(15);    visHA34may24->SetMarkerSize(0.5);   visHA34may24->SetLineStyle(7);  visHA34may24->SetLineColor(4);
    
    // canvas and plots of uv angle vs vis in small baseline groups 
    TCanvas* cuvClumps = new TCanvas("","",1400,500);
    cuvClumps->Divide(4,2);
    TGraphErrors* uv55m = new TGraphErrors();   uv55m->SetMarkerStyle(20);      uv55m->SetTitle("baseline 50-60m;uv angle from E;visibility");
    TGraphErrors* uv65m = new TGraphErrors();   uv65m->SetMarkerStyle(20);      uv65m->SetTitle("baseline 60-68m;uv angle from E;visibility");
    TGraphErrors* uv75m = new TGraphErrors();   uv75m->SetMarkerStyle(20);      uv75m->SetTitle("baseline 68-78m;uv angle from E; visibility");
    TGraphErrors* uv80m = new TGraphErrors();   uv80m->SetMarkerStyle(20);      uv80m->SetTitle("baseline 78-83m;uv angle from E;visibility");
    TGraphErrors* uv85m = new TGraphErrors();   uv85m->SetMarkerStyle(20);      uv85m->SetTitle("baseline 83-88m;uv angle from E;visibility");
    TGraphErrors* uv90m = new TGraphErrors();   uv90m->SetMarkerStyle(20);      uv90m->SetTitle("baseline 88-95m;uv angle from E;visibility");
    TGraphErrors* uv100m = new TGraphErrors();  uv100m->SetMarkerStyle(20);     uv100m->SetTitle("baseline 95-103m;uv angle from E;visibility");
    TGraphErrors* uv110m = new TGraphErrors();  uv110m->SetMarkerStyle(20);     uv110m->SetTitle("baseline 103m+;uv angle from E; visibility");
    TGraph* uvaxes = new TGraph();  uvaxes->SetMarkerStyle(9);  uvaxes->SetMarkerColor(0);  uvaxes->SetTitle("visibility by uv angle;uv angle from E;visibility");
    
    TGraphErrors* new55m = new TGraphErrors();  new55m->SetMarkerStyle(20);     //new55m->SetMarkerColor(2);
    TGraphErrors* new65m = new TGraphErrors();  new65m->SetMarkerStyle(20);     //new65m->SetMarkerColor(2);
    TGraphErrors* new75m = new TGraphErrors();  new75m->SetMarkerStyle(20);     //new75m->SetMarkerColor(2);
    TGraphErrors* new80m = new TGraphErrors();  new80m->SetMarkerStyle(20);     //new80m->SetMarkerColor(2);
    TGraphErrors* new85m = new TGraphErrors();  new85m->SetMarkerStyle(20);     //new85m->SetMarkerColor(2);
    TGraphErrors* new90m = new TGraphErrors();  new90m->SetMarkerStyle(20);     //new90m->SetMarkerColor(2);
    TGraphErrors* new100m = new TGraphErrors(); new100m->SetMarkerStyle(20);    //new100m->SetMarkerColor(2);
    TGraphErrors* new110m = new TGraphErrors(); new110m->SetMarkerStyle(20);    //new110m->SetMarkerColor(2);
    int npts55m(0), npts65m(0), npts75m(0), npts80m(0), npts85m(0), npts90m(0), npts100m(0), npts110m(0);
    int npts55new(0), npts65new(0), npts75new(0), npts80new(0), npts85new(0), npts90new(0), npts100new(0), npts110new(0);
    cuvClumps->Print("uvClumps.pdf[");
    
    // canvas and plots of distributions of parameters - peak width from fit, baselines, and uv angles
    TCanvas* cdist = new TCanvas;
    TH1D* tau12 = new TH1D("tau12","T1T2",24,-20,4);    tau12->SetLineColor(2);
    TH1D* tau13 = new TH1D("tau13","T1T3",24,-12,12);   tau13->SetLineColor(2);
    TH1D* tau14 = new TH1D("tau14","T1T4",24,-10,14);   tau14->SetLineColor(2);
    TH1D* tau23 = new TH1D("tau23","T2T3",24,-3,21);    tau23->SetLineColor(2);
    TH1D* tau24 = new TH1D("tau24","T2T4",24,-3,21);    tau24->SetLineColor(2);
    TH1D* tau34 = new TH1D("tau34","T3T4",24,-14,10);   tau34->SetLineColor(2);
    
    TH1D* taug12 = new TH1D("taug12","T1T2",24,-20,4);      
    TH1D* taug13 = new TH1D("taug13","T1T3",24,-12,12);     
    TH1D* taug14 = new TH1D("taug14","T1T4",24,-10,14);
    TH1D* taug23 = new TH1D("taug23","T2T3",24,-3,21);
    TH1D* taug24 = new TH1D("taug24","T2T4",24,-3,21);
    TH1D* taug34 = new TH1D("taug34","T3T4",24,-14,10);
    
    TH1D* basedist = new TH1D("basedist","baselines;baseline;counts",45,0,180);
    TH1D* angledist = new TH1D("angledist","uv angles;angle from cel E;counts",20,-pi/2.0, pi/2.0);
    TGraph* opdwidth = new TGraph();    opdwidth->SetMarkerStyle(20);       opdwidth->SetTitle("avg change between OPD points vs peak width;#sigma from fit;avg #deltaOPD");
    TGraph* opdwidthfull = new TGraph();    opdwidthfull->SetMarkerStyle(20);   opdwidthfull->SetTitle("max difference of OPD vs peak width;#sigma from fit; #DeltaOPD");
    TGraph* startTimeWidth = new TGraph();  startTimeWidth->SetMarkerStyle(20); startTimeWidth->SetTitle("difference in 10 MHz time and OPD calc time;#sigma from fit (ns);10MHz start time - start time used for OPD (s)");
    TGraph* opdslopedist = new TGraph();    opdslopedist->SetMarkerStyle(20);   opdslopedist->SetTitle("abs difference between start slope and end slope of OPD;#sigma from fit (ns);|m_{start} - m_{end}|");
    TH1D* widthchangehist = new TH1D("widthchangehist","change in width over run (split runs only);#sigma_{2} - #sigma_{1} (ns);counts",50,-5,5);
    // noise distributions - peak mag vs size of the noise (following gain tests)
    TGraph* noiseAlldec23 = new TGraph;     noiseAlldec23->SetMarkerStyle(20);      noiseAlldec23->SetMarkerSize(0.5);  
    TGraph* noiseAlljf23 = new TGraph;      noiseAlljf23->SetMarkerStyle(20);       noiseAlljf23->SetMarkerSize(0.5);       noiseAlljf23->SetMarkerColor(4);
    TGraph* noiseAllfeb24 = new TGraph;     noiseAllfeb24->SetMarkerStyle(20);      noiseAllfeb24->SetMarkerSize(0.5);      noiseAllfeb24->SetMarkerColor(2);
    TGraph* noiseAllmay24 = new TGraph;     noiseAllmay24->SetMarkerStyle(20);      noiseAllmay24->SetMarkerSize(0.5);      noiseAllmay24->SetMarkerColor(8);
    TGraph* noise80dec23 = new TGraph;      noise80dec23->SetMarkerStyle(20);  
    TGraph* noise80jf23 = new TGraph;       noise80jf23->SetMarkerStyle(20);        noise80jf23->SetMarkerColor(4);     // bin 42, 56, 49, 52
    TGraph* noise80feb24 = new TGraph;      noise80feb24->SetMarkerStyle(20);       noise80feb24->SetMarkerColor(2);
    TGraph* noise80may24 = new TGraph;      noise80may24->SetMarkerStyle(20);       noise80may24->SetMarkerColor(8);
    TGraph* noise107dec23 = new TGraph;     noise107dec23->SetMarkerStyle(20);  
    TGraph* noise107jf23 = new TGraph;      noise107jf23->SetMarkerStyle(20);       noise107jf23->SetMarkerColor(4);
    TGraph* noise107feb24 = new TGraph;     noise107feb24->SetMarkerStyle(20);      noise107feb24->SetMarkerColor(2);
    TGraph* noise107may24 = new TGraph;     noise107may24->SetMarkerStyle(20);      noise107may24->SetMarkerColor(8);
    TGraph* noise99dec23 = new TGraph;      noise99dec23->SetMarkerStyle(20);  
    TGraph* noise99jf23 = new TGraph;       noise99jf23->SetMarkerStyle(20);        noise99jf23->SetMarkerColor(4);
    TGraph* noise99feb24 = new TGraph;      noise99feb24->SetMarkerStyle(20);       noise99feb24->SetMarkerColor(2);
    TGraph* noise99may24 = new TGraph;      noise99may24->SetMarkerStyle(20);       noise99may24->SetMarkerColor(8);
    TGraph* noise93dec23 = new TGraph;      noise93dec23->SetMarkerStyle(20);  
    TGraph* noise93jf23 = new TGraph;       noise93jf23->SetMarkerStyle(20);        noise93jf23->SetMarkerColor(4);
    TGraph* noise93feb24 = new TGraph;      noise93feb24->SetMarkerStyle(20);       noise93feb24->SetMarkerColor(2);
    TGraph* noise93may24 = new TGraph;      noise93may24->SetMarkerStyle(20);       noise93may24->SetMarkerColor(8);
    cdist->Print("dists.pdf[");
    
    // plots to look at off run vs peak area
    TGraph* offavgVsArea = new TGraph;  offavgVsArea->SetMarkerStyle(20);   offavgVsArea->SetTitle("avg of before/after background ADCs vs peak area;avg ADC (for both pairs, before/after);peak area from fit (ns)");
    TGraph* changeOffVsArea = new TGraph;   changeOffVsArea->SetMarkerStyle(20);    changeOffVsArea->SetTitle("change in background over run vs peak area;"); // this should be +/- so we can see patterns in rising/setting etc 
    // should this be avg of the two telescopes in the pair (for change) or plot 2 pts per pair?? -> prob plot 2 pts per pair?
    // should I also record the doff from the analysis? should try to interpolate to the start/end of split runs??
    
    // uv plots for all analysis data, with midpts
    TCanvas* cuvplot = new TCanvas("","",1000,900);
    TGraph* uv12 = new TGraph();    uv12->SetMarkerStyle(20);   uv12->SetMarkerColor(jcol[0]); // my palette: 2,3,4,92,6,7
    TGraph* uv13 = new TGraph();    uv13->SetMarkerStyle(20);   uv13->SetMarkerColor(jcol[1]);
    TGraph* uv14 = new TGraph();    uv14->SetMarkerStyle(20);   uv14->SetMarkerColor(jcol[2]);
    TGraph* uv23 = new TGraph();    uv23->SetMarkerStyle(20);   uv23->SetMarkerColor(jcol[3]);
    TGraph* uv24 = new TGraph();    uv24->SetMarkerStyle(20);   uv24->SetMarkerColor(jcol[4]);
    TGraph* uv34 = new TGraph();    uv34->SetMarkerStyle(20);   uv34->SetMarkerColor(jcol[5]);
    TGraph* uvmid12 = new TGraph();     uvmid12->SetMarkerStyle(8);     uvmid12->SetMarkerColor(jcol[0]);   uvmid12->SetMarkerSize(3);
    TGraph* uvmid13 = new TGraph();     uvmid13->SetMarkerStyle(21);    uvmid13->SetMarkerColor(jcol[1]);   uvmid13->SetMarkerSize(3);
    TGraph* uvmid14 = new TGraph();     uvmid14->SetMarkerStyle(33);    uvmid14->SetMarkerColor(jcol[2]);   uvmid14->SetMarkerSize(3);
    TGraph* uvmid23 = new TGraph();     uvmid23->SetMarkerStyle(22);    uvmid23->SetMarkerColor(jcol[3]);   uvmid23->SetMarkerSize(3);
    TGraph* uvmid24 = new TGraph();     uvmid24->SetMarkerStyle(34);    uvmid24->SetMarkerColor(jcol[4]);   uvmid24->SetMarkerSize(3);
    TGraph* uvmid34 = new TGraph();     uvmid34->SetMarkerStyle(29);    uvmid34->SetMarkerColor(jcol[5]);   uvmid34->SetMarkerSize(3);
    TGraph* uvallaxes = new TGraph();   uvallaxes->SetMarkerStyle(9);   uvallaxes->SetMarkerColor(0);   uvallaxes->SetTitle("uv coverage;u#lambda (m);v#lambda (m)");
    uvallaxes->AddPoint(-170,-150);     uvallaxes->AddPoint(170,150);
    uvallaxes->GetYaxis()->SetTitleOffset(1.1);
    cuvplot->Print("uvplots.pdf[");
    
    // canvas to print "hbt stamps"
    TCanvas* cpeaks = new TCanvas;  
    cpeaks->Divide(4,3);    cpeaks->cd(1);
    cpeaks->Print("allPeaksUgly.pdf[");
    int ncpeaks(1);
    
    // canvas to print all heatmaps
    TCanvas* cheatmaps = new TCanvas;
    cheatmaps->Divide(2,2);
    cheatmaps->Print("allheatmaps.pdf[");
    int nhmp(1);
    
    // printing each run's points w avg adcs
    TCanvas* cfileadcs = new TCanvas("","",1600,700);
    cfileadcs->Divide(2,1);
    cfileadcs->Print("ADCsByRun.pdf[");
    
    TCanvas* junkcan = new TCanvas;
    //junkcan->cd();
    
    // -------------------------------------- reading in some data files, declaring output data files -------------------------------------------
    
    ofstream angles;   angles.open("angles.txt");
    
    // read in tugdual's file of tau shifts from 10 MHz signal
    string line2, name10mhz;
    string files10mhz[28];
    double shift12(0.0), shift13(0.0), shift14(0.0), shift23(0.0), shift24(0.0), shift34(0.0);
    double shifts12[28], shifts13[28], shifts14[28], shifts23[28], shifts24[28], shifts34[28];
    int itug(0);
    
    ifstream tug10file;     tug10file.open("10MHzShifts.txt");   // this is the OLD 10 mhz file - now i am doing it differently, but I will leave this artifact temporarily to not mess anything up
    
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
            
        cout << name10mhz << endl;
    }
    tug10file.close();
    
    // reading NEW 10 mhz file
    ifstream data10mhz;     data10mhz.open("all10MHzInfo.txt");
    
    int date10, hr10, min10;
    string tel10;
    double val10;
    int match10dates[999];
    int match10times[999];
    string match10tels[999];
    double match10vals[999];
    
    int i10(0);
    while(data10mhz){
        data10mhz >> date10 >> hr10 >> min10 >> tel10 >> val10;
        match10dates[i10] = date10;
        match10times[i10] = (hr10*60) + min10;
        match10tels[i10] = tel10;
        match10vals[i10] = val10;
        i10++;
        //cout << date10 << "  " << hr10 << "  " << min10 << "  " << tel10 << "  " << val10 << endl;
    }
    
    // define absolute 10 mhz references, both tau from data and its matching 10 mhz refs, for T1T2, T1T3, T1T4, T2T3, T2T4, T3T4 in that order
    double refbin[6] = {-2.8, 8.0, 7.0, 10.8, 9.8, -1.0}; //{-2.8, -17.0, -18.0, -14.2, -15.2, -1.0}; // reference 10 mhz numbers (in BINS, must be multiplied by 4ns) corresp to ref taus run
    double reftau[6] = {-7.5, -4.5, 0, 4.0, 7.2, 2.0}; // reference taus from fit to peaks from "chosen" run - 2024m02d18h19

    TFile* outfile = new TFile("savecurve.root","RECREATE");
    
    ofstream peakinfo;  peakinfo.open("goodPeakInfo.txt");
    peakinfo << "#filename             pair   area          area err       baseline    u        v        hour angle" << endl;
    
    ofstream tempcheckbaseline;    tempcheckbaseline.open("tempcheckbaseline.txt");
    
    TFile* modeldatfile = new TFile("modelInfo.root","RECREATE");
    
    ofstream taus; taus.open("taus.txt");
    taus << "dir                   T1T2   err     T1T3   err     T2T3   err     T2T4   err     T3T4    err" << endl;
    
    // ============================================================== loop over data ====================================================================================
    
    // loop over all known root files
    for (int i=0; i<nfiles; i++){
        
        cout << endl << "starting file " << filename[i] << endl;
        
        // defining an array to keep the taus for this run
        double thesetaus[6];
        double thesetausErr[6];
        
        // josie fixed this 18 june 2024 - so it is no longer necessary to read the zipped frames file
        //TString zipfilename;   zipfilename = filename[i];   zipfilename = zipfilename.ReplaceAll("analysis","ZippedFrames");
        //TFile* zipfile = new TFile("/fs/ess/PAS1977/"+zipfilename,"READONLY");
        
        // open analysis data file
        TFile* analysisfile = new TFile("/fs/ess/PAS1977/"+filename[i],"READONLY");
        int split(1);
        
        // trim file string to just time of run for plot titles
        string namedate = filename[i].Data();
        namedate = namedate.substr(namedate.find("_y")+2, namedate.length()); // when I get edit permissions in other dir, can change to "/" char
        namedate = namedate.substr(0, namedate.find("_")); // for some reason, cannot combine these into one line - must trim front and back separately
        
        // print the avg adcs to the canvas
        cfileadcs->Clear();
        cfileadcs->Divide(2,1);
        cfileadcs->cd(2);
        gPad->Divide(2,2); 
        TH1D* tempADC1all = new TH1D;   analysisfile->GetDirectory("T1T2")->GetObject("avgADCT1whole",tempADC1all);      tempADC1all->SetLineColor(2);
        TH1D* tempADC1good = new TH1D;  analysisfile->GetDirectory("T1T2")->GetObject("avgADCT1rem",tempADC1good);
        TH1D* tempADC2all = new TH1D;   analysisfile->GetDirectory("T1T2")->GetObject("avgADCT2whole",tempADC2all);      tempADC2all->SetLineColor(2);
        TH1D* tempADC2good = new TH1D;  analysisfile->GetDirectory("T1T2")->GetObject("avgADCT2rem",tempADC2good);
        cfileadcs->cd(2);   gPad->cd(1);    tempADC1all->Draw();    tempADC1good->Draw("SAME");     
        cfileadcs->cd(2);   gPad->cd(2);    tempADC2all->Draw();    tempADC2good->Draw("SAME");
        TH1D* tempADC3all = new TH1D;   TH1D* tempADC3good = new TH1D;  TH1D* tempADC4all = new TH1D;   TH1D* tempADC4good = new TH1D;
        if(analysisfile->GetDirectory("T3T4")){ 
            analysisfile->GetDirectory("T3T4")->GetObject("avgADCT3whole",tempADC3all);      tempADC3all->SetLineColor(2);
            analysisfile->GetDirectory("T3T4")->GetObject("avgADCT3rem",tempADC3good);
            analysisfile->GetDirectory("T3T4")->GetObject("avgADCT4whole",tempADC4all);      tempADC4all->SetLineColor(2);
            analysisfile->GetDirectory("T3T4")->GetObject("avgADCT4rem",tempADC4good);
            cfileadcs->cd(2);   gPad->cd(3);    tempADC3all->Draw();    tempADC3good->Draw("SAME");
            cfileadcs->cd(2);   gPad->cd(4);    tempADC4all->Draw();    tempADC4good->Draw("SAME");
        }
        else{   
            analysisfile->GetDirectory("T1T4")->GetObject("avgADCT4whole",tempADC4all);      tempADC4all->SetLineColor(2);
            analysisfile->GetDirectory("T1T4")->GetObject("avgADCT4rem",tempADC4good);
            cfileadcs->cd(2);   gPad->cd(4);    tempADC4all->Draw();    tempADC4good->Draw("SAME");
        }
        //cfileadcs->Print("ADCsByRun.pdf");
        junkcan->cd();
        TGraphErrors* thisrunpts = new TGraphErrors;    thisrunpts->SetMarkerStyle(20); // graph for this run's data only
        thisrunpts->SetTitle(Form("%s;baseline;vis", namedate.c_str()));
        
        // looking at differences in start times of runs vs peak width - checking the OPD calculation
        /*TTree* fileheader = new TTree;
        analysisfile->GetObject("header",fileheader);
        TString* T1name = new TString;      fileheader->SetBranchAddress("FileListBr",&T1name);     fileheader->GetEvent(0);
        string t1name;      t1name = T1name->Data();
        t1name = t1name.substr(0,t1name.find("\n"));
        cout << "+++++++++       the T1 filename for this run is " << t1name << endl;*/
        
        // but wait we want the earliest tdatime which is what we passed to the python .... ok it is a tdatime object but we will convert
        TDatime* ltime = new TDatime;  ltime->Set(2020,1,1,23,59,59); // init to some time far in the future, so that all times in runs would be less
        cout << "check ltime " << ltime->AsSQLString() << endl;
        TDatime* timetemp = new TDatime;
        TList* keyListTemp = gDirectory->GetListOfKeys();
        for (int ik=0; ik<keyListTemp->GetSize(); ik++){
            TKey* key = (TKey*)keyListTemp->At(ik);
            TString keyName = key->GetName();
            TString className = key->GetClassName();
            if ((keyName.Data()[0]=='T')&&(className=="TDirectoryFile")){ // if we find directory that start w "T" - take it as a pair
            TTree* pairhead = new TTree;
            analysisfile->GetDirectory(keyName)->GetObject("pairheader",pairhead);
            pairhead->SetBranchAddress("LocalTimeBr",&timetemp);   pairhead->GetEvent(0);
            cout << "this time is " << timetemp->AsSQLString() << endl;
            if(timetemp->GetTime() < ltime->GetTime()){ ltime = timetemp; }  // ok I guessed the time calc needed fixed anyway but WTF!!!! it worked to compare to the initialized "ltime" 
            }                                                               // why is it just taking the last one?!?!? ... thought i fixed this :| - ANYWAY this is what it's taking for python so what we should compare
        }
        
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
        
        cout << pyyear << " " << pymonth << " " << pyday << " " << pyhour << " " << pymin << " " << pysec << endl;
        
        // match tugdual 10mhz files and analysis files for 2024 data (that it exists for)
        string yeartemp = namedate.substr(0,4); cout << yeartemp << endl;
        string monthtemp = namedate.substr(5,2); cout << "   month temp: " << monthtemp << endl;
        
        string tugyear, tugmonth, tugday, tughour, tugmin, tugsec;
        if(yeartemp == "2024" && monthtemp == "02"){
            
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
                
                if(pyday == tugday && pyhour == tughour){   // do not need to compare year, month here bc it's all feb 2024 - but in future maybe for the peak shift
                    // there is only one case where we have 2 runs in the same hour for gam cas - so we will treat it uniquely bc this doesn't need to be so generalized
                    if(pyday == "21" && pyhour == "20"){
                        if(stoi(pymin) < 40 && stoi(tugmin) < 40){ break; }
                        else if(stoi(pymin) > 40 && stoi(tugmin) > 40){ break; }
                    }
                    else{ break; }
                }
            }
            cout << tugyear << " " << tugmonth << " " << tugday << " " << tughour << " " << tugmin << " " << tugsec << endl;
            
            // calculate difference between 10 mhz timestamp and what im sending to python from the root file
            tugtimedif += ((stoi(tugmin) - stoi(pymin))*60.0);
            if(stoi(tugsec) > stoi(pysec)){ tugtimedif += (stoi(tugsec) - stoi(pysec)); }
            else{ tugtimedif += (((stoi(tugsec) + 60.0) - stoi(pysec)) - 60.0); }  // have to do this subtract 60 add 60 thing just in case there is more than a minute discrepancy, that it would be accounted for w mins
            
            cout << endl << "time dif !!!   " << tugtimedif << endl;
        }
        
        /*
        string tempfullstr = filestr->Data();
  const int nf = nfiles;
  string telnames[nf];
  for (int i=0; i<nf; i++){
    string file = tempfullstr.substr(0, tempfullstr.find("\n"));
    telnames[i] = GetTelescopeIdentifier(file);
    tempfullstr = tempfullstr.substr(tempfullstr.find("\n")+1, tempfullstr.length());
  }*/
        
        // loop over each tdirectory in root file to find each pair
        TList* keyList = gDirectory->GetListOfKeys();
        for (int ip=0; ip<keyList->GetSize(); ip++){
            TKey* key = (TKey*)keyList->At(ip);
            TString keyName = key->GetName();
            TString className = key->GetClassName();
            if ((keyName.Data()[0]=='T')&&(className=="TDirectoryFile")){
                TString tee1 = keyName[1];    int t1 = atoi(tee1.Data()); 
                TString tee2 = keyName[3];    int t2 = atoi(tee2.Data());
                TString pair = "T"+tee1+"T"+tee2;   cout << "pair = " << pair << endl; 
                split = 1;
                
                TTree* pairhead = new TTree; // maybe not necessary ...
                analysisfile->GetDirectory(keyName)->GetObject("pairheader",pairhead);
                
                TProfile2D* cf2D = new TProfile2D;      analysisfile->GetDirectory(keyName)->GetObject("shiftedCF", cf2D);
                double runlength;   runlength = cf2D->GetYaxis()->GetNbins()*cf2D->GetYaxis()->GetBinWidth(5);
                
                TH2D* adc12D = new TH2D;    analysisfile->GetDirectory("Singles")->GetObject(Form("ADCT%s",tee1.Data()), adc12D);     adc12D->RebinY(rebin);
                TH2D* adc22D = new TH2D;    analysisfile->GetDirectory("Singles")->GetObject(Form("ADCT%s",tee2.Data()), adc22D);     adc22D->RebinY(rebin);
                TH1D* adc1 = new TH1D;      analysisfile->GetDirectory(keyName)->GetObject(Form("avgADCT%srem",tee1.Data()), adc1);     
                TH1D* adc2 = new TH1D;      analysisfile->GetDirectory(keyName)->GetObject(Form("avgADCT%srem",tee2.Data()), adc2);
                
                TNtuple* geometry = new TNtuple;        analysisfile->GetDirectory(keyName)->GetObject("PairGeometryInfo",geometry);
                geometry->Draw("baseline:time");    TGraph* fullbaseline = new TGraph(geometry->GetSelectedRows(), geometry->GetV1(), geometry->GetV2());
                cout << keyName << "   mean baseline " << fullbaseline->GetMean() << endl;
                geometry->Draw("v:u");              TGraph* uv = new TGraph(geometry->GetSelectedRows(), geometry->GetV2(), geometry->GetV1());
                //geometry->Draw("baseline:u");       TGraph* uEval = new TGraph(geometry->GetSelectedRows(), geometry->GetV1(), geometry->GetV2());
                //geometry->Draw("baseline:v");       TGraph* vEval = new TGraph(geometry->GetSelectedRows(), geometry->GetV1(), geometry->GetV2());
                geometry->Draw("u:time");           TGraph* ufull = new TGraph(geometry->GetSelectedRows(), geometry->GetV1(), geometry->GetV2());
                geometry->Draw("v:time");           TGraph* vfull = new TGraph(geometry->GetSelectedRows(), geometry->GetV1(), geometry->GetV2());
                geometry->Draw("OPD:time");         TGraph* opdgraph = new TGraph(geometry->GetSelectedRows(), geometry->GetV1(), geometry->GetV2());
                geometry->Draw("hourAngle:time");   TGraph* hrangGraph = new TGraph(geometry->GetSelectedRows(), geometry->GetV1(), geometry->GetV2());
                
                // FFTs before/after some noise removal for distributions of noise (after gain tests discovery)
                TH1D* fft1 = new TH1D;      analysisfile->GetDirectory(keyName)->GetObject("FFTProjection1",fft1); // before 80 mhz removal
                TH1D* fft2 = new TH1D;      analysisfile->GetDirectory(keyName)->GetObject("FFTProjection2",fft2); // after 80 mhz removal, before 10/107 mhz removal (should get other freqs enough)
                
                TProfile2D* heatmap = new TProfile2D;
                analysisfile->GetDirectory(keyName)->GetObject("heatmap",heatmap);
                heatmap->SetTitle(Form("%s T%sT%s", namedate.c_str(), tee1.Data(), tee2.Data()));
                cheatmaps->cd(nhmp);
                heatmap->RebinY(2);
                heatmap->Draw("COLZ"); 
                if(nhmp == 4){ nhmp = 1;  cheatmaps->Print("allheatmaps.pdf"); }
                else if(nhmp < 4){ nhmp ++; }
                junkcan->cd();
                
                // write u, v, and adcs to root file to do the baseline (u,v) smearing model
                modeldatfile->cd();
                if(pair == "T1T2"){ modeldatfile->mkdir(Form("y%s", namedate.c_str())); }
                modeldatfile->cd(Form("y%s", namedate.c_str()));
                gDirectory->mkdir(pair);
                gDirectory->cd(pair);
                adc1->Write("avgADC1");//(Form("avgADCT%s",tee1.Data()));
                adc2->Write("avgADC2");//(Form("avgADC2%s",tee2.Data()));
                ufull->Write("ugraph");//(Form("ugraph%s",pair.Data()));
                vfull->Write("vgraph");//(Form("vgraph%s",pair.Data()));
                
                analysisfile->cd();
                
                // match this run to 10 mhz reference file
                double reft1(0.0), reft2(0.0);      //reftau1(0.0), reftau2(0.0);
                for(int i10=0; i10<size(match10vals); i10++){
                    if((abs((match10dates[i10] - stoi(pyyear+pymonth+pyday))) < 1) && (abs(match10times[i10] - ((stoi(pyhour)*60) + stoi(pymin))) < 3)){ 
                        if(match10tels[i10][1] == tee1){ reft1 = match10vals[i10]; }   //cout << endl << endl << "found match   " << match10dates[i10] << "  " << match10times[i10] << "    " << tee1 << " " << match10tels[i10] << "   " << match10vals[i10] << endl;}
                        if(match10tels[i10][1] == tee2){ reft2 = match10vals[i10]; }   //cout << endl << endl << "found match   " << match10dates[i10] << "  " << tee2 << " " << match10tels[i10] << "   " << match10vals[i10] << endl;}
                    }
                }
                double thisrefbin = reft2 - reft1;
                double thistau(0.0); //pairreftau(0.0), pairrefbin(0.0);
                int iref(0);
                if(t1 == 1 && t2 == 2){ iref = 0; }  //pairreftau = reftau[0];   pairrefbin = refbin[0];
                if(t1 == 1 && t2 == 3){ iref = 1; }
                if(t1 == 1 && t2 == 4){ iref = 2; }
                if(t1 == 2 && t2 == 3){ iref = 3; }
                if(t1 == 2 && t2 == 4){ iref = 4; }
                if(t1 == 3 && t2 == 4){ iref = 5; }
                cout << endl << "thisrefbin " << thisrefbin << "    abs ref bin " << refbin[iref] << endl;
                if(abs(thisrefbin) > 12.5 && thisrefbin < 0){ thisrefbin = thisrefbin + 25.0; } // THIS IS NOT JUSTIFIED NECESSARILY!!! ask tugdual!!
                if(abs(thisrefbin) > 12.5 && thisrefbin > 0){ thisrefbin = thisrefbin - 25.0; }
                if((yeartemp == "2024" && monthtemp == "02") && pair != "T1T4"){  thistau = reftau[iref] + (4.0*(thisrefbin - refbin[iref]));  }
                else { thistau = 0; }
                
                // split run in half (or more later) if longer than 1 hr and short baseline - SHOULD BE LONGER THAN HOUR AND A HALF??? (DISCUSS??)
                //if(runlength > 5400 && fullbaseline->GetMean() < 100){ split = 2; } // subject to change baseline - but I think it's justified NOT to split when we have long baselines
                int startbin(1);
                int endbin = cf2D->GetYaxis()->GetNbins()/split;
                
                double widthchange(0), width1(0), width2(0);
                
                if(split == 1){ widthchange = 10;}
                
                // loop over number of splits
                for (int is=0; is<split; is++){
                    
                    junkcan->cd();
                    
                    // project to get HBT peak
                    TProfile* cfProj = new TProfile;
                    cfProj = cf2D->ProfileX("",startbin, endbin);
                    int flag(0);
                    
                    // calc fit parameters
                    double rtWindow(12.0);
                    double sigmaMin(2.0), sigmaMax(7.0);
                    double reltimePar(-5);
                    if((t1==1 && t2==2)){reltimePar = -8;} 
                    if((t1==1 && t2==3)){reltimePar = 0;} 
                    if((t1==1 && t2==4)){reltimePar = 2;}  
                    if((t1==2 && t2==3)){reltimePar = 9;} //10
                    if((t1==2 && t2==4)){reltimePar = 9;} //5
                    if((t1==3 && t2==4)){reltimePar = -2;} 
                    
                    // fit function -- later acct for tugdual's thing where we get the 4ns unc
                    TF1* hbtFit = new TF1("hbtFit","([0]*exp(-pow((x-[1])/[2],2)/2.0))/([2]*sqrt(2*pi))",-300,300); 
                    hbtFit->SetParName(0,"A");       hbtFit->SetParameter(0,0.0);
                    if((yeartemp == "2024" && monthtemp == "02") && pair != "T1T4"){
                        hbtFit->SetParName(1,"#tau_{0}");   hbtFit->SetParameter(1, thistau);           hbtFit->SetParLimits(1,thistau-2.0, thistau+2.0); // starting w 4ns window - may change to edges of BIN (thistau % 4)
                        //hbtFit->SetParName(1,"#tau_{0}");       hbtFit->FixParameter(1, thistau);
                    }
                    else{ 
                        hbtFit->SetParName(1,"#tau_{0}");   hbtFit->SetParameter(1, reltimePar);        hbtFit->SetParLimits(1,reltimePar-rtWindow, reltimePar+rtWindow);
                    }
                    hbtFit->SetParName(2,"#sigma");     hbtFit->SetParameter(2, 4.5);               hbtFit->SetParLimits(2, sigmaMin, sigmaMax);
                    cfProj->Fit("hbtFit","Q"); // Q = quiet mode
                    
                    /*TF1* hbtFitAmp = new TF1("hbtFitAmp","([0]*exp(-pow((x-[1])/[2],2)/2.0))",-300,300);
                    hbtFitAmp->SetParName(0,"Amp");         hbtFitAmp->SetParameter(0,0.0);
                    hbtFitAmp->SetParName(1,"#tau_{0}");    hbtFitAmp->SetParameter(1, reltimePar);     hbtFitAmp->SetParLimits(1,reltimePar-rtWindow, reltimePar+rtWindow);
                    hbtFitAmp->SetParName(2,"#sigma");      hbtFitAmp->SetParameter(2, 4.5);            hbtFitAmp->SetParLimits(2, sigmaMin, sigmaMax);*/
                    //cfProj->Fit("hbtFitAmp");
                    
                    // weighted avg for baselines, by split
                    double avgBaseNum(0.0); //avgBaseDenom(0.0);
                    double avgHrAngNum(0.0); //avgHrAngDenom(0.0);
                    double avguNum(0.0), avgvNum(0.0); //avguDenom(0.0), , avgvDenom(0.0);  
                    double avgDenom(0.0);
                    
                    // graphs to hold u, v, for this split point only - for smeared baselines, u, v, etc model
                    //TGraph* uforModel = new TGraph;
                    //TGraph* vforModel = new TGraph;
                    
                    for(int ib = startbin; ib < endbin; ib++){
                        if(abs(adc1->GetBinContent(ib)) > 0.0001 && abs(adc2->GetBinContent(ib)) > 0.0001){
                            avgBaseNum += adc1->GetBinContent(ib)*adc2->GetBinContent(ib)*fullbaseline->GetPointX(ib);
                            avgHrAngNum += adc1->GetBinContent(ib)*adc2->GetBinContent(ib)*hrangGraph->GetPointX(ib);
                            avguNum += adc1->GetBinContent(ib)*adc2->GetBinContent(ib)*ufull->GetPointX(ib);
                            avgvNum += adc1->GetBinContent(ib)*adc2->GetBinContent(ib)*vfull->GetPointX(ib);
                            avgDenom += adc1->GetBinContent(ib)*adc2->GetBinContent(ib);
                        }
                    }
                        /*TH1D* adc1 = new TH1D;      adc1 = adc12D->ProjectionX("",ib,ib); // think if rebinning is needed (I think yes) - baselines already have the right number of y bins
                        TH1D* adc2 = new TH1D;      adc2 = adc22D->ProjectionX("",ib,ib);
                        ---avgBaseNum += adc1->GetMean()*adc2->GetMean()*fullbaseline->GetPointX(ib);
                        //tempcheckbaseline << ib << "   " << fullbaseline->GetPointX(ib) << endl;
                        ---avgBaseDenom += adc1->GetMean()*adc2->GetMean();
                        avgHrAngNum += adc1->GetMean()*adc2->GetMean()*hrangGraph->GetPointX(ib);
                        avgHrAngDenom += adc1->GetMean()*adc2->GetMean();
                        avguNum += adc1->GetMean()*adc2->GetMean()*ufull->GetPointX(ib);
                        avgvNum += adc1->GetMean()*adc2->GetMean()*vfull->GetPointX(ib);
                        avguDenom += adc1->GetMean()*adc2->GetMean();
                        avgvDenom += adc1->GetMean()*adc2->GetMean();*/
                    
                    double baseline = avgBaseNum/avgDenom;
                    basedist->Fill(baseline);
                    double hrang = avgHrAngNum/avgDenom;
                    double umid = avguNum/avgDenom;
                    double vmid = avgvNum/avgDenom;
                    
                    double opdChange(0), opdsum(0), opdct(0);
                    double opdchangeStart(0), opdchangeEnd(0);
                    for(int io = startbin; io < endbin-1; io++){
                        opdsum+= abs(opdgraph->GetPointX(io+1)-opdgraph->GetPointX(io)); // IS THIS ACTUALLY RIGHT OR SHOULD I DIVIDE BY BIN WIDTHS INSTEAD OF NUM PTS???
                        if(opdct < 5){ opdchangeStart += abs(opdgraph->GetPointX(io+1) - opdgraph->GetPointX(io)); /*cout << "start slope " << opdchangeStart << endl;*/}
                        if(endbin - io < 6){ opdchangeEnd += abs(opdgraph->GetPointX(io+1) - opdgraph->GetPointX(io));  /*cout << "end slope " << opdchangeEnd << endl;*/}
                        opdct++;
                    }
                    opdChange = opdsum/double(opdct);
                    opdchangeStart = opdchangeStart/5.0;
                    opdchangeEnd = opdchangeEnd/5.0;
                    double opdMaxDif = abs(opdgraph->GetXaxis()->GetXmax() - opdgraph->GetXaxis()->GetXmin());//opdgraph->GetMaximum() - opdgraph->GetMinimum());
                    //opdChange = abs(opdChange);
                    //cout << opdChange << " " << opdMaxDif << endl;
                    
                    // eval should find the u and v pts at the weighted base avg
                    //double umid = uEval->Eval(baseline);
                    //double vmid = vEval->Eval(baseline); // prev found uv at the TIME where they matched the weighted baseline - so this would corresp to an exact point on the uv curve
                    
                    double uvAngle = TMath::ATan(vmid/umid);   //ATan2(vmid, umid);  //ATan(vmid/umid);
                    //uvAngle = uvAngle + (pi/8.0);
                    angles << uvAngle << endl;
                    angledist->Fill(uvAngle);
                    
                    // calculate larger errorbar, set errorbars to either existing errorbar from fit or the new larger errorbar
                    double excessErr = calcExcessErr(cfProj,hbtFit->GetParameter(0), hbtFit->GetParameter(2), 100);
                    double errorbar = excessErr; // hbtFit->GetParError(0); // set the errorbar for every plot to either errors from fit or the errors from simulated peaks RMS - in one line
                    
                    cout << keyName << "  vis: " << hbtFit->GetParameter(0) << "  fit err: " << hbtFit->GetParError(0) << "  exc err: " << excessErr << "   baseline: " << baseline << "  u: " << umid << "  v: " << vmid << endl;
                    
                    // quality cuts
                    //if(abs(hbtFit->GetParameter(1)-reltimePar-rtWindow) < 0.01 || abs(hbtFit->GetParameter(1)-reltimePar+rtWindow) < 0.01){ flag = 2; }
                    if((yeartemp == "2024" && monthtemp == "02") && pair != "T1T4"){
                        if(abs(hbtFit->GetParameter(1)-thistau-2.0) < 0.01 || abs(hbtFit->GetParameter(1)-thistau+2.0) < 0.01){ flag = 2; }
                    }
                    else{
                        if(abs(hbtFit->GetParameter(1)-reltimePar-rtWindow) < 0.01 || abs(hbtFit->GetParameter(1)-reltimePar+rtWindow) < 0.01){ flag = 2; }
                    }
                    if(abs(hbtFit->GetParameter(2) - sigmaMin) < 0.01 || abs(hbtFit->GetParameter(2) - sigmaMax) < 0.01){ flag = 3; }
                    
                    //cout << endl << endl << "compare ref tau!!! " << hbtFit->GetParameter(1) << "     " << thistau << endl << endl << endl; 
                    
                    // calculate significance threshold for signal / noise ratio - determine if peak is in muck or not
                    double sigthresh(0);
                    if(hbtFit->GetParameter(0) < 5e-6){  sigthresh = calcThresh(cfProj); } // in reality - this should be done for ALL points. but it takes forever, so this is to save time
                    if((hbtFit->GetParameter(0) < sigthresh) && (flag < 2)){ 
                        flag = -2; 
                        cout << endl << "peak below significance threshold " << hbtFit->GetParameter(0) << "  thresh: " << sigthresh << endl << endl;
                        allptsLow->AddPoint(baseline, hbtFit->GetParameter(0));
                    }
                    
                    // draw each hbt peak, print to pdf
                    cpeaks->cd(ncpeaks);
                    cfProj->SetTitle(Form("%s T%sT%s B= %4.1f flag %d", namedate.c_str(), tee1.Data(), tee2.Data(), baseline, flag)); // why need c_str?
                    cfProj->Rebin(4);
                    cfProj->GetXaxis()->SetRangeUser(-64,64);
                    cpeaks->cd(ncpeaks);  
                    cfProj->Draw();
                    if(flag > 1){ TLine* bad = new TLine(-64,0,64,cfProj->GetMaximum());    bad->SetLineColor(2);   bad->Draw(); }
                    gPad->Update();
                    TLine* left;
                    TLine* right;
                    double stamptop = cfProj->GetMaximum()*1.2;
                    double stampbottom = cfProj->GetMinimum()*1.2;
                    TLine* threshline = new TLine(-64,sigthresh,64,sigthresh);  threshline->SetLineColor(8);
                    if((yeartemp == "2024" && monthtemp == "02") && pair != "T1T4"){
                        left = new TLine(thistau-4.0, stampbottom, thistau-4.0, stamptop); // 1e-6
                        right = new TLine(thistau+4.0, stampbottom, thistau+4.0, stamptop);
                    }
                    else{
                        left = new TLine(reltimePar-10, stampbottom, reltimePar-10, stamptop);
                        right = new TLine(reltimePar+10, stampbottom, reltimePar+10, stamptop);
                    }
                    //cpeaks->cd(ncpeaks);  
                    left->SetLineWidth(2);  right->SetLineWidth(2);
                    left->Draw("SAME");     right->Draw("SAME");
                    if(sigthresh > 1e-9){ threshline->Draw("SAME"); }
                    ncpeaks++;
                    if(ncpeaks == 13){ ncpeaks = 1;  cpeaks->Print("allPeaksUgly.pdf");     cpeaks->Clear();    cpeaks->Divide(4,3); cpeaks->cd(1);}
                    junkcan->cd();
                    
                    // make substring to check year, month and group by observation run
                    string year, month;
                    year = namedate.substr(0,4); // substr(start, #chars after start)
                    month = namedate.substr(5,2);
                    
                    // get off adcs saved in file - this is so clunky and ungraceful but better than matching the off times again 
                    TH1D* off1 = new TH1D;      analysisfile->GetDirectory(keyName)->GetObject("offADCs1",off1);
                    TH1D* off2 = new TH1D;      analysisfile->GetDirectory(keyName)->GetObject("offADCs2",off2);
                    double offbefore1(0), offafter1(0), offbefore2(0), offafter2(0), offavg(0);
                    offbefore1 = off1->GetBinContent(1);    offafter1 = off1->GetBinContent(2);
                    offbefore2 = off2->GetBinContent(1);    offafter2 = off2->GetBinContent(2);
                    if(offbefore1 != 0 && offafter1!= 0 && offbefore2 != 0 && offafter2 != 0){ offavg = (offbefore1 + offafter1 + offbefore2 + offafter2)/4.0; }
                    cout << "offavg " << offavg << endl;
                    
                    // fill noise plots (for both good AND BAD peaks, outside if "good")
                    if(year == "2024"){
                        if(month == "02"){
                            noise80feb24->AddPoint(hbtFit->GetParameter(2), fft1->GetBinContent(42));
                            noise107feb24->AddPoint(hbtFit->GetParameter(2), fft2->GetBinContent(56));
                            noise93feb24->AddPoint(hbtFit->GetParameter(2), fft2->GetBinContent(49));
                            noise99feb24->AddPoint(hbtFit->GetParameter(2), fft2->GetBinContent(52));
                            for(int in=1; in<=fft2->GetXaxis()->GetNbins(); in++){
                                noiseAllfeb24->AddPoint((fft2->GetBinCenter(in)*1e9)/(2.0*pi), fft2->GetBinContent(in));
                            }
                        }
                        else if(month == "05"){
                            noise80may24->AddPoint(hbtFit->GetParameter(2), fft1->GetBinContent(42));
                            noise107may24->AddPoint(hbtFit->GetParameter(2), fft2->GetBinContent(56));
                            noise93may24->AddPoint(hbtFit->GetParameter(2), fft2->GetBinContent(49));
                            noise99may24->AddPoint(hbtFit->GetParameter(2), fft2->GetBinContent(52));
                            for(int in=1; in<=fft2->GetXaxis()->GetNbins(); in++){
                                noiseAllmay24->AddPoint((fft2->GetBinCenter(in)*1e9)/(2.0*pi), fft2->GetBinContent(in));
                            }
                        }
                    }
                    else if(year == "2023"){
                        if(month == "01" || month == "02"){
                            noise80jf23->AddPoint(hbtFit->GetParameter(2), fft1->GetBinContent(42));
                            noise107jf23->AddPoint(hbtFit->GetParameter(2), fft2->GetBinContent(56));
                            noise93jf23->AddPoint(hbtFit->GetParameter(2), fft2->GetBinContent(49));
                            noise99jf23->AddPoint(hbtFit->GetParameter(2), fft2->GetBinContent(52));
                            for(int in=1; in<=fft2->GetXaxis()->GetNbins(); in++){
                                noiseAlljf23->AddPoint((fft2->GetBinCenter(in)*1e9)/(2.0*pi), fft2->GetBinContent(in));
                            }
                        }
                        else if(month == "12"){
                            noise80dec23->AddPoint(hbtFit->GetParameter(2), fft1->GetBinContent(42));
                            noise107dec23->AddPoint(hbtFit->GetParameter(2), fft2->GetBinContent(56));
                            noise93dec23->AddPoint(hbtFit->GetParameter(2), fft2->GetBinContent(49));
                            noise99dec23->AddPoint(hbtFit->GetParameter(2), fft2->GetBinContent(52));
                            for(int in=1; in<=fft2->GetXaxis()->GetNbins(); in++){
                                noiseAlldec23->AddPoint((fft2->GetBinCenter(in)*1e9)/(2.0*pi), fft2->GetBinContent(in));
                            }
                        }
                    }
                    
                    if(t1==1 && t2==2){ tau12->Fill(hbtFit->GetParameter(1));   thesetaus[0] = 0;   thesetausErr[0] = 0; }
                    if(t1==1 && t2==3){ tau13->Fill(hbtFit->GetParameter(1));   thesetaus[1] = 0;   thesetausErr[1] = 0; }
                    if(t1==1 && t2==4){ tau14->Fill(hbtFit->GetParameter(1));   thesetaus[2] = 0;   thesetausErr[2] = 0; }
                    if(t1==2 && t2==3){ tau23->Fill(hbtFit->GetParameter(1));   thesetaus[3] = 0;   thesetausErr[3] = 0; }
                    if(t1==2 && t2==4){ tau24->Fill(hbtFit->GetParameter(1));   thesetaus[4] = 0;   thesetausErr[4] = 0; }
                    if(t1==3 && t2==4){ tau34->Fill(hbtFit->GetParameter(1));   thesetaus[5] = 0;   thesetausErr[5] = 0; }
                    
                    // -------------------------------------------- fill plots with "good" points ---------------------------------------------------------
                    
                    if(hbtFit->GetParError(0) > 0.0001){ cout << endl << endl << endl << "HUGE ERRORBAR!! " << namedate << endl << endl << endl << endl;}
                    else{
                    // plot points
                    if(flag <=1){ 
                        allpts->AddPoint(baseline, hbtFit->GetParameter(0));     
                        allpts->SetPointError(npts, 0, 0, errorbar, errorbar);
                        npts++;
                        
                        peakinfo << namedate << "   T" << t1 << "T" << t2 << "   " << hbtFit->GetParameter(0) << "   " << errorbar << "   " << baseline << "   " << umid << "   " << vmid << "   " << hrang << endl;
                        thisrunpts->AddPoint(baseline, hbtFit->GetParameter(0));
                        thisrunpts->SetPointError(thisrunpts->GetN()-1, 0, hbtFit->GetParError(0));
                        
                        if(is == 0){ width1 = hbtFit->GetParameter(2); }
                        if(is == 1){ width2 = hbtFit->GetParameter(2); }
                        
                        if(t1==1 && t2==2){ taug12->Fill(hbtFit->GetParameter(1));  thesetaus[0] = hbtFit->GetParameter(1);     thesetausErr[0] = hbtFit->GetParError(1); }
                        if(t1==1 && t2==3){ taug13->Fill(hbtFit->GetParameter(1));  thesetaus[1] = hbtFit->GetParameter(1);     thesetausErr[1] = hbtFit->GetParError(1); }
                        if(t1==1 && t2==4){ taug14->Fill(hbtFit->GetParameter(1));  thesetaus[2] = hbtFit->GetParameter(1);     thesetausErr[2] = hbtFit->GetParError(1); }
                        if(t1==2 && t2==3){ taug23->Fill(hbtFit->GetParameter(1));  thesetaus[3] = hbtFit->GetParameter(1);     thesetausErr[3] = hbtFit->GetParError(1); }
                        if(t1==2 && t2==4){ taug24->Fill(hbtFit->GetParameter(1));  thesetaus[4] = hbtFit->GetParameter(1);     thesetausErr[4] = hbtFit->GetParError(1); }
                        if(t1==3 && t2==4){ taug34->Fill(hbtFit->GetParameter(1));  thesetaus[5] = hbtFit->GetParameter(1);     thesetausErr[5] = hbtFit->GetParError(1); }
                        
                        offavgVsArea->AddPoint(offavg, hbtFit->GetParameter(0));
                        
                         if(year == "2024" && month == "02"){
                            //if(baseline > 120 || hbtFit->GetParameter(0) > 0){
                                viscurveFeb2024->AddPoint(baseline, hbtFit->GetParameter(0));
                                viscurveFeb2024->SetPointError(npts2024, 0, 0, errorbar, errorbar);
                                npts2024++;
                                viscurveFebMay2024->AddPoint(baseline, hbtFit->GetParameter(0));
                                viscurveFebMay2024->SetPointError(nptsfm2024, 0, 0, errorbar, errorbar);
                                nptsfm2024++;
                                viscurveNoJan->AddPoint(baseline, hbtFit->GetParameter(0));
                                viscurveNoJan->SetPointError(nptsnj, 0, 0, errorbar, errorbar);
                                nptsnj++;
                                startTimeWidth->AddPoint(hbtFit->GetParameter(2), tugtimedif);
                                
                            //}
                        }
                        if(year == "2024" && month == "05"){
                            viscurveMay2024->AddPoint(baseline, hbtFit->GetParameter(0));
                            viscurveMay2024->SetPointError(nptsm2024, 0, 0, errorbar, errorbar);
                            nptsm2024++;
                            viscurveFebMay2024->AddPoint(baseline, hbtFit->GetParameter(0));
                            viscurveFebMay2024->SetPointError(nptsfm2024, 0, 0, errorbar, errorbar);
                            nptsfm2024++;
                            viscurveNoJan->AddPoint(baseline, hbtFit->GetParameter(0));
                            viscurveNoJan->SetPointError(nptsnj, 0, 0, errorbar, errorbar);
                            nptsnj++;
                            viscurveNoFeb->AddPoint(baseline, hbtFit->GetParameter(0));
                            viscurveNoFeb->SetPointError(nptsnf, 0, 0, errorbar, errorbar);
                            nptsnf++;
                            viscurveDecMay->AddPoint(baseline, hbtFit->GetParameter(0));
                            viscurveDecMay->SetPointError(nptsdm, 0, 0, errorbar, errorbar);
                            nptsdm++;
                        }
                        if((year == "2023" && month == "01") || (year == "2023" && month == "02")){
                            viscurveJanFeb2023->AddPoint(baseline, hbtFit->GetParameter(0));
                            viscurveJanFeb2023->SetPointError(nptsjf2023, 0, 0, errorbar, errorbar);
                            nptsjf2023++;
                            viscurveNoFeb->AddPoint(baseline, hbtFit->GetParameter(0));
                            viscurveNoFeb->SetPointError(nptsnf, 0, 0, errorbar, errorbar);
                            nptsnf++;
                        }
                        if(year == "2023" && month == "12"){
                            viscurveDec2023->AddPoint(baseline, hbtFit->GetParameter(0));
                            viscurveDec2023->SetPointError(nptsd2023, 0, 0, errorbar, errorbar);
                            nptsd2023++;
                            viscurveNoJan->AddPoint(baseline, hbtFit->GetParameter(0));
                            viscurveNoJan->SetPointError(nptsnj, 0, 0, errorbar, errorbar);
                            nptsnj++;
                            viscurveNoFeb->AddPoint(baseline, hbtFit->GetParameter(0));
                            viscurveNoFeb->SetPointError(nptsnf, 0, 0, errorbar, errorbar);
                            nptsnf++;
                            viscurveDecMay->AddPoint(baseline, hbtFit->GetParameter(0));
                            viscurveDecMay->SetPointError(nptsdm, 0, 0, errorbar, errorbar);
                            nptsdm++;
                        }
                        //viscurve2D->AddPoint(ucoord, vcoord, hbtFit->GetParameter(0));
                        // plot uv angle, split by baseline
                        if(baseline < 50){ cout << endl << endl << "BASELINE LESS THAN 50!" << endl << endl << endl;}
                        if(baseline > 50 && baseline <= 60){
                            uv55m->AddPoint(uvAngle, hbtFit->GetParameter(0));
                            uv55m->SetPointError(npts55m, 0, errorbar);
                            npts55m++;
                            if(year == "2024"){     
                                new55m->AddPoint(uvAngle, hbtFit->GetParameter(0));
                                new55m->SetPointError(npts55new, 0, errorbar);
                                npts55new++;
                            }
                        }
                        if(baseline > 60 && baseline <= 68){
                            uv65m->AddPoint(uvAngle, hbtFit->GetParameter(0));
                            uv65m->SetPointError(npts65m, 0, errorbar);
                            npts65m++;
                            if(year == "2024"){     
                                new65m->AddPoint(uvAngle, hbtFit->GetParameter(0));
                                new65m->SetPointError(npts65new, 0, errorbar);
                                npts65new++;
                            }
                        }
                        if(baseline > 68 && baseline <= 78){
                            uv75m->AddPoint(uvAngle, hbtFit->GetParameter(0));
                            uv75m->SetPointError(npts75m, 0, errorbar);
                            npts75m++;
                            if(year == "2024"){     
                                new75m->AddPoint(uvAngle, hbtFit->GetParameter(0));
                                new75m->SetPointError(npts75new, 0, errorbar);
                                npts75new++;
                            }
                        }
                        if(baseline > 78 && baseline <= 83){
                            cout << "special pt " << uvAngle << "  " << hbtFit->GetParameter(0) << endl << endl << endl;
                            uv80m->AddPoint(uvAngle, hbtFit->GetParameter(0));
                            uv80m->SetPointError(npts80m, 0, errorbar);
                            npts80m++;
                            if(year == "2024"){     
                                new80m->AddPoint(uvAngle, hbtFit->GetParameter(0));
                                new80m->SetPointError(npts80new, 0, errorbar);
                                npts80new++;
                            }
                        }
                        if(baseline > 83 && baseline <= 88){
                            uv85m->AddPoint(uvAngle, hbtFit->GetParameter(0));
                            uv85m->SetPointError(npts85m, 0, errorbar);
                            npts85m++;
                            if(year == "2024"){     
                                new85m->AddPoint(uvAngle, hbtFit->GetParameter(0));
                                new85m->SetPointError(npts85new, 0, errorbar);
                                npts85new++;
                            }
                        }
                        if(baseline > 88 && baseline <= 95){
                            uv90m->AddPoint(uvAngle, hbtFit->GetParameter(0));
                            uv90m->SetPointError(npts90m, 0, errorbar);
                            npts90m++;
                            if(year == "2024"){     
                                new90m->AddPoint(uvAngle, hbtFit->GetParameter(0));
                                new90m->SetPointError(npts90new, 0, errorbar);
                                npts90new++;
                            }
                        }
                        if(baseline > 95 && baseline <= 103){
                            uv100m->AddPoint(uvAngle, hbtFit->GetParameter(0));
                            uv100m->SetPointError(npts100m, 0, errorbar);
                            npts100m++;
                            if(year == "2024"){     
                                new100m->AddPoint(uvAngle, hbtFit->GetParameter(0));
                                new100m->SetPointError(npts100new, 0, errorbar);
                                npts100new++;
                            }
                        }
                        if(baseline > 103){
                            uv110m->AddPoint(uvAngle, hbtFit->GetParameter(0));
                            uv110m->SetPointError(npts110m, 0, errorbar);
                            npts110m++;
                            if(year == "2024"){     
                                new110m->AddPoint(uvAngle, hbtFit->GetParameter(0));
                                new110m->SetPointError(npts110new, 0, errorbar);
                                npts110new++;
                            }
                        }
                        /*
                        // plot vis curve, split by uv angle
                        if(uvAngle > 3.0*pi/8.0 || uvAngle < -3.0*pi/8.0){
                            viscurve1->AddPoint(baseline, hbtFit->GetParameter(0));
                            viscurve1->SetPointError(npts1, 0, hbtFit->GetParError(0));
                            npts1++;
                        }
                        if(uvAngle > -3.0*pi/8.0 && uvAngle < -pi/8.0){
                            viscurve2->AddPoint(baseline, hbtFit->GetParameter(0));
                            viscurve2->SetPointError(npts2, 0, hbtFit->GetParError(0));
                            npts2++;
                        }
                        if(uvAngle > -pi/8.0 && uvAngle < pi/8.0){
                            viscurve3->AddPoint(baseline, hbtFit->GetParameter(0));
                            viscurve3->SetPointError(npts3, 0, hbtFit->GetParError(0));
                            npts3++;
                        }
                        if(uvAngle > pi/8.0 && uvAngle < 3.0*pi/8.0){
                            viscurve4->AddPoint(baseline, hbtFit->GetParameter(0));
                            viscurve4->SetPointError(npts4, 0, hbtFit->GetParError(0));
                            npts4++;
                        }*/
                        
                        // plot uv points, split by pair
                        if(t1 == 1 && t2 == 2){ 
                            hrangleT1T2->AddPoint(hrang, hbtFit->GetParameter(0));
                            hrangleT1T2->SetPointError(nptsha12, 0, 0, errorbar, errorbar);
                            nptsha12++;
                            for(int ib=startbin; ib<=endbin; ib++){
                                uv12->AddPoint(uv->GetPointX(ib), uv->GetPointY(ib));
                                uv12->AddPoint(-uv->GetPointX(ib), -uv->GetPointY(ib));
                            }
                            uvmid12->AddPoint(umid, vmid);      uvmid12->AddPoint(-umid, -vmid);
                            if(year == "2024" && month == "02"){ 
                                ha12feb24->AddPoint(hrang, hbtFit->GetParameter(0));   ha12feb24->SetPointError(nha12feb24, 0, errorbar); nha12feb24++; 
                                vcfeb24T1T2->AddPoint(baseline, hbtFit->GetParameter(0));   vcfeb24T1T2->SetPointError(npts12f24, 0, 0, errorbar, errorbar);   npts12f24++;
                            }
                            if(year == "2024" && month == "05"){ 
                                ha12may24->AddPoint(hrang, hbtFit->GetParameter(0));   ha12may24->SetPointError(nha12may24, 0, errorbar); nha12may24++; 
                                vcmay24T1T2->AddPoint(baseline, hbtFit->GetParameter(0));   vcmay24T1T2->SetPointError(npts12m24, 0, 0, errorbar, errorbar);   npts12m24++;
                            }
                            if(year == "2023" && month == "12"){ 
                                ha12dec23->AddPoint(hrang, hbtFit->GetParameter(0));   ha12dec23->SetPointError(nha12dec23, 0, errorbar); nha12dec23++; 
                                vcdec23T1T2->AddPoint(baseline, hbtFit->GetParameter(0));   vcdec23T1T2->SetPointError(npts12d23, 0, 0, errorbar, errorbar);   npts12d23++;
                            }
                            if((year == "2023" && month == "01") || (year == "2023" && month == "02")){ 
                                ha12jf23->AddPoint(hrang, hbtFit->GetParameter(0));     ha12jf23->SetPointError(nha12jf23, 0, errorbar); nha12jf23++; 
                                vcjf23T1T2->AddPoint(baseline, hbtFit->GetParameter(0));   vcjf23T1T2->SetPointError(npts12jf23, 0, 0, errorbar, errorbar);   npts12jf23++;
                            }
                        }
                        if(t1 == 1 && t2 == 3){ 
                            hrangleT1T3->AddPoint(hrang, hbtFit->GetParameter(0));
                            hrangleT1T3->SetPointError(nptsha13, 0, 0, errorbar, errorbar);
                            nptsha13++;
                            for(int ib=startbin; ib<=endbin; ib++){
                                uv13->AddPoint(uv->GetPointX(ib), uv->GetPointY(ib));
                                uv13->AddPoint(-uv->GetPointX(ib), -uv->GetPointY(ib));
                            }
                            uvmid13->AddPoint(umid, vmid);      uvmid13->AddPoint(-umid, -vmid);
                            if(year == "2024" && month == "02"){ 
                                ha13feb24->AddPoint(hrang, hbtFit->GetParameter(0));   ha13feb24->SetPointError(nha13feb24, 0, errorbar); nha13feb24++; 
                                vcfeb24T1T3->AddPoint(baseline, hbtFit->GetParameter(0));   vcfeb24T1T3->SetPointError(npts13f24, 0, 0, errorbar, errorbar);   npts13f24++;
                            }
                            if(year == "2024" && month == "05"){ 
                                ha13may24->AddPoint(hrang, hbtFit->GetParameter(0));   ha13may24->SetPointError(nha13may24, 0, errorbar); nha13may24++; 
                                vcmay24T1T3->AddPoint(baseline, hbtFit->GetParameter(0));   vcmay24T1T3->SetPointError(npts13m24, 0, 0, errorbar, errorbar);   npts13m24++;
                            }
                            if(year == "2023" && month == "12"){ 
                                ha13dec23->AddPoint(hrang, hbtFit->GetParameter(0));   ha13dec23->SetPointError(nha13dec23, 0, errorbar); nha13dec23++; 
                                vcdec23T1T3->AddPoint(baseline, hbtFit->GetParameter(0));   vcdec23T1T3->SetPointError(npts13d23, 0, 0, errorbar, errorbar);   npts13d23++;
                            }
                            if((year == "2023" && month == "01") || (year == "2023" && month == "02")){ 
                                ha13jf23->AddPoint(hrang, hbtFit->GetParameter(0));     ha13jf23->SetPointError(nha13jf23, 0, errorbar); nha13jf23++; 
                                vcjf23T1T3->AddPoint(baseline, hbtFit->GetParameter(0));   vcjf23T1T3->SetPointError(npts13jf23, 0, 0, errorbar, errorbar);   npts13jf23++;
                            }
                        }
                        if(t1 == 1 && t2 == 4){ 
                            hrangleT1T4->AddPoint(hrang, hbtFit->GetParameter(0));
                            hrangleT1T4->SetPointError(nptsha14, 0, 0, errorbar, errorbar);
                            nptsha14++;
                            for(int ib=startbin; ib<=endbin; ib++){
                                uv14->AddPoint(uv->GetPointX(ib), uv->GetPointY(ib));
                                uv14->AddPoint(-uv->GetPointX(ib), -uv->GetPointY(ib));
                            }
                            uvmid14->AddPoint(umid, vmid);      uvmid14->AddPoint(-umid, -vmid);
                            if(year == "2024" && month == "02"){ 
                                ha14feb24->AddPoint(hrang, hbtFit->GetParameter(0));   ha14feb24->SetPointError(nha14feb24, 0, errorbar); nha14feb24++;
                                vcfeb24T1T4->AddPoint(baseline, hbtFit->GetParameter(0));   vcfeb24T1T4->SetPointError(npts14f24, 0, 0, errorbar, errorbar);   npts14f24++;
                            }
                            if(year == "2024" && month == "05"){ 
                                ha14may24->AddPoint(hrang, hbtFit->GetParameter(0));   ha14may24->SetPointError(nha14may24, 0, errorbar); nha14may24++;
                                vcmay24T1T4->AddPoint(baseline, hbtFit->GetParameter(0));   vcmay24T1T4->SetPointError(npts14m24, 0, 0, errorbar, errorbar);   npts14m24++;
                            }
                            if(year == "2023" && month == "12"){ 
                                ha14dec23->AddPoint(hrang, hbtFit->GetParameter(0));   ha14dec23->SetPointError(nha14dec23, 0, errorbar); nha14dec23++; 
                                vcdec23T1T4->AddPoint(baseline, hbtFit->GetParameter(0));   vcdec23T1T4->SetPointError(npts14d23, 0, 0, errorbar, errorbar);   npts14d23++;
                            }
                            if((year == "2023" && month == "01") || (year == "2023" && month == "02")){ 
                                ha14jf23->AddPoint(hrang, hbtFit->GetParameter(0));     ha14jf23->SetPointError(nha14jf23, 0, errorbar); nha14jf23++; 
                                vcjf23T1T4->AddPoint(baseline, hbtFit->GetParameter(0));   vcjf23T1T4->SetPointError(npts14jf23, 0, 0, errorbar, errorbar);   npts14jf23++;
                            }
                        }
                        if(t1 == 2 && t2 == 3){
                            hrangleT2T3->AddPoint(hrang, hbtFit->GetParameter(0));
                            hrangleT2T3->SetPointError(nptsha23, 0, 0, errorbar, errorbar);
                            nptsha23++;
                            for(int ib=startbin; ib<=endbin; ib++){
                                uv23->AddPoint(uv->GetPointX(ib), uv->GetPointY(ib));
                                uv23->AddPoint(-uv->GetPointX(ib), -uv->GetPointY(ib));
                            }
                            uvmid23->AddPoint(umid, vmid);      uvmid23->AddPoint(-umid, -vmid);
                            if(year == "2024" && month == "02"){ 
                                ha23feb24->AddPoint(hrang, hbtFit->GetParameter(0));   ha23feb24->SetPointError(nha23feb24, 0, errorbar); nha23feb24++; 
                                vcfeb24T2T3->AddPoint(baseline, hbtFit->GetParameter(0));   vcfeb24T2T3->SetPointError(npts23f24, 0, 0, errorbar, errorbar);   npts23f24++;
                            }
                            if(year == "2024" && month == "05"){ 
                                ha23may24->AddPoint(hrang, hbtFit->GetParameter(0));   ha23may24->SetPointError(nha23may24, 0, errorbar); nha23may24++;
                                vcmay24T2T3->AddPoint(baseline, hbtFit->GetParameter(0));   vcmay24T2T3->SetPointError(npts23m24, 0, 0, errorbar, errorbar);   npts23m24++;
                            }
                            if(year == "2023" && month == "12"){ 
                                ha23dec23->AddPoint(hrang, hbtFit->GetParameter(0));   ha23dec23->SetPointError(nha23dec23, 0, errorbar); nha23dec23++; 
                                vcdec23T2T3->AddPoint(baseline, hbtFit->GetParameter(0));   vcdec23T2T3->SetPointError(npts23d23, 0, 0, errorbar, errorbar);   npts23d23++;
                            }
                            if((year == "2023" && month == "01") || (year == "2023" && month == "02")){ 
                                ha23jf23->AddPoint(hrang, hbtFit->GetParameter(0));     ha23jf23->SetPointError(nha23jf23, 0, errorbar); nha23jf23++; 
                                vcjf23T2T3->AddPoint(baseline, hbtFit->GetParameter(0));   vcjf23T2T3->SetPointError(npts23jf23, 0, 0, errorbar, errorbar);   npts23jf23++;
                            }
                        }
                        if(t1 == 2 && t2 == 4){ 
                            hrangleT2T4->AddPoint(hrang, hbtFit->GetParameter(0));
                            hrangleT2T4->SetPointError(nptsha24, 0, 0, errorbar, errorbar);
                            nptsha24++;
                            for(int ib=startbin; ib<=endbin; ib++){
                                uv24->AddPoint(uv->GetPointX(ib), uv->GetPointY(ib));
                                uv24->AddPoint(-uv->GetPointX(ib), -uv->GetPointY(ib));
                            }
                            uvmid24->AddPoint(umid, vmid);      uvmid24->AddPoint(-umid, -vmid);
                            if(year == "2024" && month == "02"){ 
                                ha24feb24->AddPoint(hrang, hbtFit->GetParameter(0));   ha24feb24->SetPointError(nha24feb24, 0, errorbar); nha24feb24++; 
                                vcfeb24T2T4->AddPoint(baseline, hbtFit->GetParameter(0));   vcfeb24T2T4->SetPointError(npts24f24, 0, 0, errorbar, errorbar);   npts24f24++;
                            }
                            if(year == "2024" && month == "05"){ 
                                ha24may24->AddPoint(hrang, hbtFit->GetParameter(0));   ha24may24->SetPointError(nha24may24, 0, errorbar); nha24may24++; 
                                vcmay24T2T4->AddPoint(baseline, hbtFit->GetParameter(0));   vcmay24T2T4->SetPointError(npts24m24, 0, 0, errorbar, errorbar);   npts24m24++;
                            }
                            if(year == "2023" && month == "12"){ 
                                ha24dec23->AddPoint(hrang, hbtFit->GetParameter(0));   ha24dec23->SetPointError(nha24dec23, 0, errorbar); nha24dec23++; 
                                vcdec23T2T4->AddPoint(baseline, hbtFit->GetParameter(0));   vcdec23T2T4->SetPointError(npts24d23, 0, 0, errorbar, errorbar);   npts24d23++;
                            }
                            if((year == "2023" && month == "01") || (year == "2023" && month == "02")){ 
                                ha24jf23->AddPoint(hrang, hbtFit->GetParameter(0));     ha24jf23->SetPointError(nha24jf23, 0, errorbar); nha24jf23++; 
                                vcjf23T2T4->AddPoint(baseline, hbtFit->GetParameter(0));   vcjf23T2T4->SetPointError(npts24jf23, 0, 0, errorbar, errorbar);   npts24jf23++;
                            }
                        }
                        if(t1 == 3 && t2 == 4){ 
                            hrangleT3T4->AddPoint(hrang, hbtFit->GetParameter(0));
                            hrangleT3T4->SetPointError(nptsha34, 0, 0, errorbar, errorbar);
                            nptsha34++;
                            for(int ib=startbin; ib<=endbin; ib++){
                                uv34->AddPoint(uv->GetPointX(ib), uv->GetPointY(ib));
                                uv34->AddPoint(-uv->GetPointX(ib), -uv->GetPointY(ib));
                            }
                            uvmid34->AddPoint(umid, vmid);      uvmid34->AddPoint(-umid, -vmid);
                            if(year == "2024" && month == "02"){ 
                                ha34feb24->AddPoint(hrang, hbtFit->GetParameter(0));   ha34feb24->SetPointError(nha34feb24, 0, errorbar); nha34feb24++;
                                vcfeb24T3T4->AddPoint(baseline, hbtFit->GetParameter(0));   vcfeb24T3T4->SetPointError(npts34f24, 0, 0, errorbar, errorbar);   npts34f24++;}
                            if(year == "2024" && month == "05"){ 
                                ha34may24->AddPoint(hrang, hbtFit->GetParameter(0));   ha34may24->SetPointError(nha34may24, 0, errorbar); nha34may24++; 
                                vcmay24T3T4->AddPoint(baseline, hbtFit->GetParameter(0));   vcmay24T3T4->SetPointError(npts34m24, 0, 0, errorbar, errorbar);   npts34m24++;
                            }
                            if(year == "2023" && month == "12"){ 
                                ha34dec23->AddPoint(hrang, hbtFit->GetParameter(0));   ha34dec23->SetPointError(nha34dec23, 0, errorbar); nha34dec23++; 
                                vcdec23T3T4->AddPoint(baseline, hbtFit->GetParameter(0));   vcdec23T3T4->SetPointError(npts34d23, 0, 0, errorbar, errorbar);   npts34d23++;
                            }
                            if((year == "2023" && month == "01") || (year == "2023" && month == "02")){ 
                                ha34jf23->AddPoint(hrang, hbtFit->GetParameter(0));     ha34jf23->SetPointError(nha34jf23, 0, errorbar); nha34jf23++;
                                vcjf23T3T4->AddPoint(baseline, hbtFit->GetParameter(0));   vcjf23T3T4->SetPointError(npts34jf23, 0, 0, errorbar, errorbar);   npts34jf23++;}
                        }
                        
                        opdwidth->AddPoint(hbtFit->GetParameter(2), opdChange);
                        opdwidthfull->AddPoint(hbtFit->GetParameter(2), opdMaxDif);
                        opdslopedist->AddPoint(hbtFit->GetParameter(2), abs(opdchangeStart - opdchangeEnd));
                    } // end of if pass cuts
                    }
                    startbin += cf2D->GetYaxis()->GetNbins()/split;
                    endbin   += cf2D->GetYaxis()->GetNbins()/split;
                    
                } // end of split loop
                
                if(split == 2){ 
                    widthchange = width2-width1; 
                    widthchangehist->Fill(widthchange);
                }
                
            
            } // end of pair dir loop
        } // end of all objs in file loop
        
        cfileadcs->cd(1);
        thisrunpts->AddPoint(-1,-6e-6); thisrunpts->AddPoint(190,15e-6);
        thisrunpts->GetXaxis()->SetRangeUser(0,180);    thisrunpts->GetYaxis()->SetRangeUser(-5e-6,14e-6);
        thisrunpts->Draw("AP");
        TF1* refcurve = new TF1("refcurve","fabs([0])*pow(2.0*TMath::BesselJ1(TMath::Pi()*x*[1]*1e-3/(4157e-10*206265))/(TMath::Pi()*x*[1]*1e-3/(4157e-10*206265)),2)", 0, 180);
        refcurve->SetParName(0, "C_{UD}");                 refcurve->FixParameter(0, 8.07e-6);   //0.505877  8.06574e-06
        refcurve->SetParName(1, "#theta_{UD} (mas)");      refcurve->FixParameter(1, 0.506);
        refcurve->Draw("SAME");
        cfileadcs->Print("ADCsByRun.pdf");
        
        taus << namedate << "   " << thesetaus[0] << "  " << thesetausErr[0] << "   " << thesetaus[1] << "  " << thesetausErr[1] << "   " << thesetaus[3] << "  " << thesetausErr[3] << "   " << thesetaus[4] << "  " << thesetausErr[4] << "   " << thesetaus[5] << "  " << thesetausErr[5] << endl;
        
    } // end of root files loop
    // close pdf
    cpeaks->Print("allPeaksUgly.pdf");
    cpeaks->Print("allPeaksUgly.pdf]");
    
    cheatmaps->Print("allheatmaps.pdf");
    cheatmaps->Print("allheatmaps.pdf]");
    
    cfileadcs->Print("ADCsByRun.pdf]");
    
    // now read in pts from old data vis curve
    double oldbase(0.0), oldbaseerr(0.0), oldvis(0.0), oldviserr(0.0), junk1(0.0), junk2(0.0);
    string line;
    ifstream oldptsfile;    oldptsfile.open("PointsInVisibilityCurve.txt");
    
    // read in old fpga data (pre 2023) to add to plot -- eventually we will read FPGA data root files to do the full analysis
    TGraphErrors* oldpts = new TGraphErrors();  oldpts->SetMarkerColor(2);  oldpts->SetMarkerStyle(20);
    int nOld(0);
    while(getline(oldptsfile, line)){
        oldptsfile >> oldbase >> oldbaseerr >> oldvis >> oldviserr >> junk1 >> junk2;
        //cout << oldbase << " " << oldbaseerr << " " << oldvis << " " << oldviserr << " " << junk1 << " " << junk2 << endl;
        oldpts->AddPoint(oldbase, oldvis);      oldpts->SetPointError(nOld, oldbaseerr, oldviserr);
        nOld++;
    }
    
    angles.close();
    peakinfo.close();
    tempcheckbaseline.close();
    taus.close();
    
    // ======================================================= make plots =============================================================
    
    // draw uv plot with midpts
    cuvplot->cd();
    uvallaxes->Draw("APRX");
    uv12->Draw("PRXSAME");    uvmid12->Draw("PRXSAME");
    uv13->Draw("PRXSAME");    uvmid13->Draw("PRXSAME");
    uv14->Draw("PRXSAME");    uvmid14->Draw("PRXSAME");
    uv23->Draw("PRXSAME");    uvmid23->Draw("PRXSAME");
    uv24->Draw("PRXSAME");    uvmid24->Draw("PRXSAME");
    uv34->Draw("PRXSAME");    uvmid34->Draw("PRXSAME");
    TLegend* uvplotleg = new TLegend(0.75,0.6,0.9,0.9);
    uvplotleg->AddEntry(uvmid12,"T1T2","P");
    uvplotleg->AddEntry(uvmid13,"T1T3","P");
    uvplotleg->AddEntry(uvmid14,"T1T4","P");
    uvplotleg->AddEntry(uvmid23,"T2T3","P");
    uvplotleg->AddEntry(uvmid24,"T2T4","P");
    uvplotleg->AddEntry(uvmid34,"T3T4","P");
    uvplotleg->Draw("SAME");
    cuvplot->Print("uvplots.pdf");
    cuvplot->Print("uvplots.pdf]");
    
    //--------------------------------------------- "distribution"/comparison plots -------------------------------------------------------------
    
    // draw tau dists
    cdist->cd(1);
    cdist->Clear();
    cdist->Print("dists.pdf");
    cdist->Divide(3,2);
    cdist->cd(1);   tau12->Draw();  taug12->Draw("SAME");
    cdist->cd(2);   tau13->Draw();  taug13->Draw("SAME");
    cdist->cd(3);   tau14->Draw();  taug14->Draw("SAME");   
    cdist->cd(4);   tau23->Draw();  taug23->Draw("SAME");
    cdist->cd(5);   tau24->Draw();  taug24->Draw("SAME");
    cdist->cd(6);   tau34->Draw();  taug34->Draw("SAME");
    cout << "checking taus " << tau12->GetMean() << " " << tau13->GetMean() << " " << tau14->GetMean() << " " << tau23->GetMean() << " " << tau24->GetMean() << " " << tau34->GetMean() << endl;
    cout << "good only " << taug12->GetMean() << " " << taug13->GetMean() << " " << taug14->GetMean() << " " << taug23->GetMean() << " " << taug24->GetMean() << " " << taug34->GetMean() << endl;
    cdist->Print("dists.pdf");
    cdist->cd();
    // draw other dists
    TFile* distfile = new TFile("distfilerebin.root","RECREATE");
    basedist->Draw();
    cdist->Print("dists.pdf");
    angledist->Draw();
    cdist->Print("dists.pdf");
    opdwidth->Draw("AP");
    opdwidth->Write("delOPDwidth");
    cdist->Print("dists.pdf");
    opdwidthfull->Draw("AP");
    opdwidthfull->Write("maxdifOPD");
    cdist->Print("dists.pdf");
    startTimeWidth->Draw("AP");
    cdist->Print("dists.pdf");
    opdslopedist->Draw("AP");
    opdslopedist->Write("slopeOPD");
    cdist->Print("dists.pdf");
    widthchangehist->Draw();
    widthchangehist->Write("widthchange");
    cdist->Print("dists.pdf");
    offavgVsArea->Draw("AP");
    cdist->Print("dists.pdf");
    
    // draw noise plots
    TLegend* noisedates = new TLegend(0.8,0.7,0.95,0.95);
    noisedates->AddEntry(noiseAlljf23, "Jan/Feb 2023","P");
    noisedates->AddEntry(noiseAlldec23, "Dec 2023","P");
    noisedates->AddEntry(noiseAllfeb24, "Feb 2024","P");
    noisedates->AddEntry(noiseAllmay24, "May 2024","P");
    cout << "max of noise plots: feb24 " << noiseAllfeb24->GetMaximum() << "  may24 " << noiseAllmay24->GetMaximum() << "  jf23 " << noiseAlljf23->GetMaximum() << "  dec23 " << noiseAlldec23->GetMaximum() << endl;
    noiseAllfeb24->SetTitle("noise by observation (excl 80 MHz);frequency (Hz);magnitude of noise");
    noiseAllfeb24->Draw("AP");
    noiseAlldec23->Draw("PSAME");
    noiseAllmay24->Draw("PSAME");
    //noiseAlljf23->Draw("PSAME");
    noisedates->Draw("SAME");
    cdist->Print("dists.pdf");
    noise80feb24->SetTitle("mag of 80 MHz noise vs peak;peak width;mag of noise");
    noise80feb24->Draw("AP");
    noise80may24->Draw("PSAME");
    noise80jf23->Draw("PSAME");
    noise80dec23->Draw("PSAME");
    noisedates->Draw("SAME");
    cdist->Print("dists.pdf");
    noise107feb24->SetTitle("mag of 107 MHz noise vs peak;peak width;mag of noise");
    noise107feb24->Draw("AP");
    noise107may24->Draw("PSAME");
    noise107jf23->Draw("PSAME");
    noise107dec23->Draw("PSAME");
    noisedates->Draw("SAME");
    cdist->Print("dists.pdf");
    noise93feb24->SetTitle("mag of 93 MHz noise vs peak;peak width;mag of noise");
    noise93feb24->Draw("AP");
    noise93may24->Draw("PSAME");
    noise93jf23->Draw("PSAME");
    noise93dec23->Draw("PSAME");
    noisedates->Draw("SAME");
    cdist->Print("dists.pdf");
    noise99feb24->SetTitle("mag of 99 MHz noise vs peak;peak width;mag of noise");
    noise99feb24->Draw("AP");
    noise99may24->Draw("PSAME");
    noise99jf23->Draw("PSAME");
    noise99dec23->Draw("PSAME");
    noisedates->Draw("SAME");
    cdist->Print("dists.pdf");
    cdist->Print("dists.pdf]");
    
    // ------------------------------------------------- "vis curve" plots ---------------------------------------------------
    
    // draw vis curve
    allpts->AddPoint(0,0);
    cvis->cd();
    allpts->Draw("AP");
    allptsLow->Draw("PSAME");
    //oldpts->Draw("PSAME");
    TLine* zero = new TLine(0,0,180,0);
    zero->Draw("SAME");
    cvis->Print("viscurve.pdf");
    
    oldpts->Draw("PSAME");
    TLegend* corrleg = new TLegend(0.75,0.8,0.9,0.9);
    corrleg->AddEntry(allpts, "versii", "P");
    corrleg->AddEntry(oldpts, "FPGA", "P");
    corrleg->Draw("SAME");
    cvis->Print("viscurve.pdf");
    
    // fit function for uniform disk (not correct for gam cas but an ok first approx)
    TF1* udfit = new TF1("udfit","fabs([0])*pow(2.0*TMath::BesselJ1(TMath::Pi()*x*[1]*1e-3/(4157e-10*206265))/(TMath::Pi()*x*[1]*1e-3/(4157e-10*206265)),2)", 0, 180);
    udfit->SetParName(0, "C_{UD}");                 udfit->SetParameter(0, 1e-5);
    udfit->SetParName(1, "#theta_{UD} (mas)");      udfit->SetParameter(1, 1.0);
    
    // variables to keep the fit results for each set of data
    double diamAll, normAll, diamJF23, normJF23, diamDec23, normDec23, diamFeb24, normFeb24, diamMay24, normMay24;
    
    /*TMinuit* udminuit = new TMinuit();
    udminuit->SetFCN(uniformDiskFunc);
    Double_t arglist[10];       arglist[0] = 1;
    Int_t ierflg(0);
    udminuit->mnexcm("SET ERR", arglist, 1, ierflg);
    udminuit->mnparm(0,"normalization",1e-5,1e-8,0,0,ierflg);        //mnparm([par],"name",start,inc,min,max);
    udminuit->mnparm(1,"diameter",1.0,0.0001,0,0,ierflg);*/
    
    /*gMinuit->SetFCN(InterpolateJason);
    Double_t arglist[10];
    Int_t ierflg=0;
    arglist[0]=1;
    gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
    gMinuit->mnparm(0,"diameter",milliarcValues[lowestDir],0.003,0,0,ierflg);
    gMinuit->mnparm(1,"Scaling",ConstantGraph->GetPointY(lowestDir),0.1,0,0,ierflg);*/
    
    //SetFCN(void(*fcn)(Int_t &, Double_t *, Double_t &f, Double_t *, Int_t))
    
    /*void uniformDiskFunc(Int_t &npar, Double_t &f, Double_t *par, Int_t iflag){
        double norm = par[0];
        double angdiam = par[1];
        
        f = fabs(norm)*pow(2.0*TMath::BesselJ1(TMath::Pi()*x*angdiam*1e-3/(4157e-10*206265))/(TMath::Pi()*x*angdiam*1e-3/(4157e-10*206265)),2);
    }*/
    
    allpts->Draw("AP");
    //viscurveJanFeb2023->Draw("PSAME");
    //viscurveFeb2024->Draw("PSAME");
    //viscurveMay2024->Draw("PSAME");
    zero->Draw("SAME");
    TLegend* dateleg = new TLegend(0.7,0.6,0.9,0.75);
    dateleg->AddEntry(allpts, "Dec 2023", "P");
    dateleg->AddEntry(viscurveJanFeb2023, "Jan/Feb 2023", "P");
    dateleg->AddEntry(viscurveFeb2024, "Feb 2024", "P");
    dateleg->AddEntry(viscurveMay2024, "May 2024", "P");
    //dateleg->Draw("SAME");
    cvis->Print("viscurve.pdf");
    
    outfile->cd();
    allpts->Write("allpts");
    viscurveJanFeb2023->Write("jf2023");
    viscurveDec2023->Write("dec2023");
    viscurveFeb2024->Write("feb2024");
    viscurveMay2024->Write("may2024");
    //outfile->Close();
    
    allpts->Fit(udfit);
    normAll = udfit->GetParameter(0);
    diamAll = udfit->GetParameter(1);
    //udminuit->mnexcm("MIGRAD",arglist,2,ierflg);
    cvis->Print("viscurve.pdf");
    
    TH1D* temperr1sigAll = new TH1D("temperrs1sig","temp hist to hold confidence intervals from fit - 1 sigma", 100, 0, 10000);
    TH1D* temperr2sigAll = new TH1D("temperrs2sig","temp hist to hold confidence intervals from fit - 2 sigma", 100, 0, 10000);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(temperr1sigAll, 0.68); // this is 1 sigma - but with normal dist, not chi squared
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(temperr2sigAll, 0.95); // 2 sigma
    
    gMinuit->SetErrorDef(1);
    TGraph* errcont1sigAll = new TGraph;        errcont1sigAll = (TGraph*)gMinuit->Contour(80, 1, 0);       errcont1sigAll->SetLineWidth(3);    errcont1sigAll->SetLineColor(15);
    gMinuit->SetErrorDef(4);    
    TGraph* errcont2sigAll = new TGraph;        errcont2sigAll = (TGraph*)gMinuit->Contour(80, 1, 0);       errcont2sigAll->SetLineWidth(2);    errcont2sigAll->SetLineColor(15);   errcont2sigAll->SetLineStyle(2);
    
    // add baseline lines
 /*   TLine* b60 = new TLine(60,-2e-6,60,18e-6);      b60->Draw("SAME");
    TLine* b68 = new TLine(68,-2e-6,68,18e-6);      b68->Draw("SAME");
    TLine* b78 = new TLine(78,-2e-6,78,18e-6);      b78->Draw("SAME");
    TLine* b83 = new TLine(83,-2e-6,83,18e-6);      b83->Draw("SAME");
    TLine* b88 = new TLine(88,-2e-6,88,18e-6);      b88->Draw("SAME");
    TLine* b95 = new TLine(95,-2e-6,95,18e-6);      b95->Draw("SAME");
    TLine* b103 = new TLine(103,-2e-6,103,18e-6);   b103->Draw("SAME");
    cvis->Print("viscurve.pdf");*/
    
    viscurveJanFeb2023->AddPoint(0,0);
    viscurveJanFeb2023->Draw("AP");
    viscurveJanFeb2023->RemovePoint();
    udfit->SetParameter(0, 1e-5);       udfit->SetParameter(1, 1.0);
    viscurveJanFeb2023->Fit(udfit,"0");
    viscurveJanFeb2023->GetFunction("udfit")->Draw("SAME");
    normJF23 = udfit->GetParameter(0);
    diamJF23 = udfit->GetParameter(1);
    zero->Draw("SAME");
    cvis->Print("viscurve.pdf");
    
    TH1D* temperr1sigJF23 = new TH1D("temperrs1sigJF23","temp hist to hold confidence intervals from fit - 1 sigma", 100, 0, 10000);
    TH1D* temperr2sigJF23 = new TH1D("temperrs2sigJF23","temp hist to hold confidence intervals from fit - 2 sigma", 100, 0, 10000);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(temperr1sigJF23, 0.68); // this is 1 sigma - but with normal dist, not chi squared
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(temperr2sigJF23, 0.95); // 2 sigma
    
    gMinuit->SetErrorDef(1);
    TGraph* errcont1sigJF23 = new TGraph;       errcont1sigJF23 = (TGraph*)gMinuit->Contour(80, 1, 0);      errcont1sigJF23->SetLineWidth(3);       errcont1sigJF23->SetLineColor(4);
    gMinuit->SetErrorDef(4);    
    TGraph* errcont2sigJF23 = new TGraph;       errcont2sigJF23 = (TGraph*)gMinuit->Contour(80, 1, 0);      errcont2sigJF23->SetLineWidth(2);       errcont2sigJF23->SetLineColor(4);   errcont2sigJF23->SetLineStyle(2);
    
    viscurveDec2023->AddPoint(0,0);
    viscurveDec2023->Draw("AP");
    viscurveDec2023->RemovePoint();
    udfit->SetParameter(0, 1e-5);       udfit->SetParameter(1, 1.0);
    viscurveDec2023->Fit(udfit,"0");
    viscurveDec2023->GetFunction("udfit")->Draw("SAME");
    normDec23 = udfit->GetParameter(0);
    diamDec23 = udfit->GetParameter(1);
    zero->Draw("SAME");
    cvis->Print("viscurve.pdf");
    
    TH1D* temperr1sigDec23 = new TH1D("temperrs1sigDec23","temp hist to hold confidence intervals from fit - 1 sigma", 100, 0, 10000);
    TH1D* temperr2sigDec23 = new TH1D("temperrs2sigDec23","temp hist to hold confidence intervals from fit - 2 sigma", 100, 0, 10000);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(temperr1sigDec23, 0.68); // this is 1 sigma - but with normal dist, not chi squared
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(temperr2sigDec23, 0.95); // 2 sigma
    
    gMinuit->SetErrorDef(1);
    TGraph* errcont1sigDec23 = new TGraph;       errcont1sigDec23 = (TGraph*)gMinuit->Contour(80, 1, 0);      errcont1sigDec23->SetLineWidth(3);       errcont1sigDec23->SetLineColor(1);
    gMinuit->SetErrorDef(4);    
    TGraph* errcont2sigDec23 = new TGraph;       errcont2sigDec23 = (TGraph*)gMinuit->Contour(80, 1, 0);      errcont2sigDec23->SetLineWidth(2);       errcont2sigDec23->SetLineColor(1);   errcont2sigDec23->SetLineStyle(2);
    
    viscurveFeb2024->AddPoint(0,0);
    viscurveFeb2024->Draw("AP");
    viscurveFeb2024->RemovePoint();
    udfit->SetParameter(0, 1e-5);       udfit->SetParameter(1, 1.0);
    viscurveFeb2024->Fit(udfit,"0");
    viscurveFeb2024->GetFunction("udfit")->Draw("SAME");
    normFeb24 = udfit->GetParameter(0);
    diamFeb24 = udfit->GetParameter(1);
    zero->Draw("SAME");
    cvis->Print("viscurve.pdf");
    
    TH1D* temperr1sigFeb24 = new TH1D("temperrs1sigFeb24","temp hist to hold confidence intervals from fit - 1 sigma", 100, 0, 10000);
    TH1D* temperr2sigFeb24 = new TH1D("temperrs2sigFeb24","temp hist to hold confidence intervals from fit - 2 sigma", 100, 0, 10000);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(temperr1sigFeb24, 0.68); // this is 1 sigma - but with normal dist, not chi squared
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(temperr2sigFeb24, 0.95); // 2 sigma
    
    gMinuit->SetErrorDef(1);
    TGraph* errcont1sigFeb24 = new TGraph;       errcont1sigFeb24 = (TGraph*)gMinuit->Contour(80, 1, 0);      errcont1sigFeb24->SetLineWidth(2);       errcont1sigFeb24->SetLineColor(2);
    gMinuit->SetErrorDef(4);    
    TGraph* errcont2sigFeb24 = new TGraph;       errcont2sigFeb24 = (TGraph*)gMinuit->Contour(80, 1, 0);      errcont2sigFeb24->SetLineWidth(2);       errcont2sigFeb24->SetLineColor(2);   errcont2sigFeb24->SetLineStyle(2);
    
    viscurveMay2024->AddPoint(0,0);
    viscurveMay2024->Draw("AP");
    viscurveMay2024->RemovePoint();
    udfit->SetParameter(0, 1e-5);       udfit->SetParameter(1, 1.0);
    viscurveMay2024->Fit(udfit,"0");
    viscurveMay2024->GetFunction("udfit")->Draw("SAME");
    normMay24 = udfit->GetParameter(0);
    diamMay24 = udfit->GetParameter(1);
    zero->Draw("SAME");
    cvis->Print("viscurve.pdf");
    
    TH1D* temperr1sigMay24 = new TH1D("temperrs1sigMay24","temp hist to hold confidence intervals from fit - 1 sigma", 100, 0, 10000);
    TH1D* temperr2sigMay24 = new TH1D("temperrs2sigMay24","temp hist to hold confidence intervals from fit - 2 sigma", 100, 0, 10000);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(temperr1sigMay24, 0.68); // this is 1 sigma - but with normal dist, not chi squared
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(temperr2sigMay24, 0.95); // 2 sigma
    
    gMinuit->SetErrorDef(1);
    TGraph* errcont1sigMay24 = new TGraph;       errcont1sigMay24 = (TGraph*)gMinuit->Contour(80, 1, 0);      errcont1sigMay24->SetLineWidth(2);       errcont1sigMay24->SetLineColor(8);
    gMinuit->SetErrorDef(4);    
    TGraph* errcont2sigMay24 = new TGraph;       errcont2sigMay24 = (TGraph*)gMinuit->Contour(80, 1, 0);      errcont2sigMay24->SetLineWidth(2);       errcont2sigMay24->SetLineColor(8);   errcont2sigMay24->SetLineStyle(2);
    
    cvis->Clear();
    viscurveFebMay2024->AddPoint(0,0);
    viscurveFebMay2024->Draw("AP");
    viscurveFebMay2024->RemovePoint();
    udfit->SetParameter(0, 1e-5);       udfit->SetParameter(1, 1.0);
    viscurveFebMay2024->Fit(udfit,"0");
    viscurveFebMay2024->GetFunction("udfit")->Draw("SAME");
    double normFebMay24 = udfit->GetParameter(0);
    double diamFebMay24 = udfit->GetParameter(1);
    zero->Draw("SAME");
    viscurveFeb2024->Draw("PSAME"); 
    viscurveMay2024->Draw("PSAME"); 
    //TPaveStats* statboxFebMay = (TPaveStats*)viscurveFebMay2024->GetListOfFunctions()->FindObject("stats");
    //statboxFebMay->Draw("SAME");
    cvis->Print("viscurve.pdf");
    
    viscurveNoJan->AddPoint(0,0);
    viscurveNoJan->Draw("AP");
    viscurveNoJan->RemovePoint();
    udfit->SetParameter(0, 1e-5);       udfit->SetParameter(1, 1.0);
    viscurveNoJan->Fit(udfit,"0");
    viscurveNoJan->GetFunction("udfit")->Draw("SAME");
    double normNoJan = udfit->GetParameter(0);
    double diamNoJan = udfit->GetParameter(1);
    zero->Draw("SAME");
    viscurveDec2023->Draw("PSAME");     
    viscurveFeb2024->Draw("PSAME");     
    viscurveMay2024->Draw("PSAME");    
    //TPaveStats* statboxNoJan = (TPaveStats*)viscurveNoJan->GetListOfFunctions()->FindObject("stats");
    //statboxNoJan->Draw("SAME");
    cvis->Print("viscurve.pdf");
    
    viscurveNoFeb->AddPoint(0,0);
    viscurveNoFeb->Draw("AP");
    viscurveNoFeb->RemovePoint();
    udfit->SetParameter(0, 1e-5);       udfit->SetParameter(1, 1.0);
    viscurveNoFeb->Fit(udfit,"0");
    viscurveNoFeb->GetFunction("udfit")->Draw("SAME");
    double normNoFeb = udfit->GetParameter(0);
    double diamNoFeb = udfit->GetParameter(1);
    zero->Draw("SAME");
    viscurveJanFeb2023->Draw("PSAME");
    viscurveDec2023->Draw("PSAME");
    viscurveMay2024->Draw("PSAME");
    //TPaveStats* statboxNoFeb = (TPaveStats*)viscurveNoFeb->GetListOfFunctions()->FindObject("stats");
    //statboxNoFeb->Draw("SAME");
    cvis->Print("viscurve.pdf");
    
    viscurveDecMay->AddPoint(0,0);
    viscurveDecMay->Draw("AP");
    viscurveDecMay->RemovePoint();
    udfit->SetParameter(0, 1e-5);       udfit->SetParameter(1, 1.0);
    viscurveDecMay->Fit(udfit,"0");
    viscurveDecMay->GetFunction("udfit")->Draw("SAME");
    double normDecMay = udfit->GetParameter(0);
    double diamDecMay = udfit->GetParameter(1);
    zero->Draw("SAME");
    viscurveDec2023->Draw("PSAME");
    viscurveMay2024->Draw("PSAME");
    cvis->Print("viscurve.pdf");
    
    // plots by pair for each observation period
    TLegend* vcpairs = new TLegend(0.75,0.7,0.9,0.9);
    vcpairs->AddEntry(vcjf23T1T2, "T1T2", "P");
    vcpairs->AddEntry(vcjf23T1T3, "T1T3", "P");
    vcpairs->AddEntry(vcjf23T1T4, "T1T4", "P");
    vcpairs->AddEntry(vcjf23T2T3, "T2T3", "P");
    vcpairs->AddEntry(vcjf23T2T4, "T2T4", "P");
    vcpairs->AddEntry(vcjf23T3T4, "T3T4", "P");
    
    vcjf23T1T2->SetTitle("Jan/Feb 2023 by pairs;baseline (m);peak fit from area (ns)");
    viscurveJanFeb2023->Draw("AP");     viscurveJanFeb2023->GetFunction("udfit")->Draw("SAME");
    zero->Draw("SAME");
    vcjf23T1T2->Draw("PSAME");
    vcjf23T1T3->Draw("PSAME");
    vcjf23T1T4->Draw("PSAME");
    vcjf23T2T3->Draw("PSAME");
    vcjf23T2T4->Draw("PSAME");
    vcjf23T3T4->Draw("PSAME");
    vcpairs->Draw("SAME");
    cvis->Print("viscurve.pdf");
    
    vcdec23T1T2->SetTitle("Dec 2023 by pairs;baseline (m);peak fit from area (ns)");
    viscurveDec2023->Draw("AP");    viscurveDec2023->GetFunction("udfit")->Draw("SAME");
    zero->Draw("SAME");
    vcdec23T1T2->Draw("PSAME");
    vcdec23T1T3->Draw("PSAME");
    vcdec23T1T4->Draw("PSAME");
    vcdec23T2T3->Draw("PSAME");
    vcdec23T2T4->Draw("PSAME");
    vcdec23T3T4->Draw("PSAME");
    vcpairs->Draw("SAME");
    cvis->Print("viscurve.pdf");
    
    vcfeb24T1T2->SetTitle("Feb 2024 by pairs;baseline (m);peak fit from area (ns)");
    viscurveFeb2024->Draw("AP");    viscurveFeb2024->GetFunction("udfit")->Draw("SAME");
    zero->Draw("SAME");
    vcfeb24T1T2->Draw("PSAME");
    vcfeb24T1T3->Draw("PSAME");
    vcfeb24T1T4->Draw("PSAME");
    vcfeb24T2T3->Draw("PSAME");
    vcfeb24T2T4->Draw("PSAME");
    vcfeb24T3T4->Draw("PSAME");
    vcpairs->Draw("SAME");
    cvis->Print("viscurve.pdf");
    
    vcmay24T1T2->SetTitle("May 2024 by pairs;baseline (m);peak fit from area (ns)");
    viscurveMay2024->Draw("AP");    viscurveMay2024->GetFunction("udfit")->Draw("SAME");
    zero->Draw("SAME");
    vcmay24T1T2->Draw("PSAME");
    vcmay24T1T3->Draw("PSAME");
    vcmay24T1T4->Draw("PSAME");
    vcmay24T2T3->Draw("PSAME");
    vcmay24T2T4->Draw("PSAME");
    vcmay24T3T4->Draw("PSAME");
    vcpairs->Draw("SAME");
    cvis->Print("viscurve.pdf");
    
    TGraph* errcontaxes = new TGraph;   errcontaxes->SetMarkerStyle(9); errcontaxes->SetMarkerColor(0);     errcontaxes->AddPoint(0.4, 10e-6);   errcontaxes->AddPoint(0.68,25e-6);
    errcontaxes->SetTitle("error contours for uniform disk fit to visibility curve;#theta_{UD} (mas);normalization C_{UD}");
    errcontaxes->Draw("AP");
    errcont2sigAll->Draw("LSAME");
    errcont1sigAll->Draw("LSAME");
    errcont2sigJF23->Draw("LSAME");
    errcont1sigJF23->Draw("LSAME");
    errcont2sigDec23->Draw("LSAME");
    errcont1sigDec23->Draw("LSAME");
    errcont2sigFeb24->Draw("LSAME");
    errcont1sigFeb24->Draw("LSAME");
    errcont2sigMay24->Draw("LSAME");
    errcont1sigMay24->Draw("LSAME");
    
    TLegend* legcont = new TLegend(0.1,0.7,0.3,0.9);
    legcont->AddEntry(errcont1sigAll,"overall","L");
    legcont->AddEntry(errcont1sigJF23,"Jan/Feb 2023","L");
    legcont->AddEntry(errcont1sigDec23,"Dec 2023","L");
    legcont->AddEntry(errcont1sigFeb24,"Feb 2024","L");
    legcont->AddEntry(errcont1sigMay24,"May 2024","L");
    legcont->Draw("SAME");
    cvis->Print("viscurve.pdf");
    
    /*axes->Draw("AP");
    viscurve1->Draw("PSAME");
    viscurve2->Draw("PSAME");
    viscurve3->Draw("PSAME");
    viscurve4->Draw("PSAME");
    zero->Draw("SAME");
    TLegend* vcleg = new TLegend(0.7,0.7,0.9,0.9);
    vcleg->AddEntry(viscurve1, "-5#pi/8- -3#pi/8","P");
    vcleg->AddEntry(viscurve2, "-3#pi/8 - -#pi/8", "P");
    vcleg->AddEntry(viscurve3, "-#pi/8 - #pi/8", "P");
    vcleg->AddEntry(viscurve4, "#pi/8 - 3#pi/8","P");
    vcleg->Draw("SAME");
    cvis->Print("viscurve.pdf");*/
    cvis->Print("viscurve.pdf]");
    
    // ------------------------------------------------------- uv and hour angle plots ---------------------------------------------
    
    // normalize visibilities, using different month normalization constants
    int i12(0), i13(0), i14(0), i23(0), i24(0), i34(0);
    for(int ig=0; ig<ha12jf23->GetN(); ig++){   haT1T2norm->AddPoint(ha12jf23->GetPointX(ig), ha12jf23->GetPointY(ig)/normJF23);    haT1T2norm->SetPointError(i12, 0, 0, ha12jf23->GetErrorYlow(ig)/normJF23, ha12jf23->GetErrorYhigh(ig)/normJF23);    i12++; }
    for(int ig=0; ig<ha13jf23->GetN(); ig++){   haT1T3norm->AddPoint(ha13jf23->GetPointX(ig), ha13jf23->GetPointY(ig)/normJF23);    haT1T3norm->SetPointError(i13, 0, 0, ha13jf23->GetErrorYlow(ig)/normJF23, ha13jf23->GetErrorYhigh(ig)/normJF23);    i13++; }
    for(int ig=0; ig<ha14jf23->GetN(); ig++){   haT1T4norm->AddPoint(ha14jf23->GetPointX(ig), ha14jf23->GetPointY(ig)/normJF23);    haT1T4norm->SetPointError(i14, 0, 0, ha14jf23->GetErrorYlow(ig)/normJF23, ha14jf23->GetErrorYhigh(ig)/normJF23);    i14++; }
    for(int ig=0; ig<ha23jf23->GetN(); ig++){   haT2T3norm->AddPoint(ha23jf23->GetPointX(ig), ha23jf23->GetPointY(ig)/normJF23);    haT2T3norm->SetPointError(i23, 0, 0, ha23jf23->GetErrorYlow(ig)/normJF23, ha23jf23->GetErrorYhigh(ig)/normJF23);    i23++; }
    for(int ig=0; ig<ha24jf23->GetN(); ig++){   haT2T4norm->AddPoint(ha24jf23->GetPointX(ig), ha24jf23->GetPointY(ig)/normJF23);    haT2T4norm->SetPointError(i24, 0, 0, ha24jf23->GetErrorYlow(ig)/normJF23, ha24jf23->GetErrorYhigh(ig)/normJF23);    i24++; }
    for(int ig=0; ig<ha34jf23->GetN(); ig++){   haT3T4norm->AddPoint(ha34jf23->GetPointX(ig), ha34jf23->GetPointY(ig)/normJF23);    haT3T4norm->SetPointError(i34, 0, 0, ha34jf23->GetErrorYlow(ig)/normJF23, ha34jf23->GetErrorYhigh(ig)/normJF23);    i34++; }
    
    for(int ig=0; ig<ha12dec23->GetN(); ig++){   haT1T2norm->AddPoint(ha12dec23->GetPointX(ig), ha12dec23->GetPointY(ig)/normDec23);    haT1T2norm->SetPointError(i12, 0, 0, ha12dec23->GetErrorYlow(ig)/normDec23, ha12dec23->GetErrorYhigh(ig)/normDec23);    i12++; }
    for(int ig=0; ig<ha13dec23->GetN(); ig++){   haT1T3norm->AddPoint(ha13dec23->GetPointX(ig), ha13dec23->GetPointY(ig)/normDec23);    haT1T3norm->SetPointError(i13, 0, 0, ha13dec23->GetErrorYlow(ig)/normDec23, ha13dec23->GetErrorYhigh(ig)/normDec23);    i13++; }
    for(int ig=0; ig<ha14dec23->GetN(); ig++){   haT1T4norm->AddPoint(ha14dec23->GetPointX(ig), ha14dec23->GetPointY(ig)/normDec23);    haT1T4norm->SetPointError(i14, 0, 0, ha14dec23->GetErrorYlow(ig)/normDec23, ha14dec23->GetErrorYhigh(ig)/normDec23);    i14++; }
    for(int ig=0; ig<ha23dec23->GetN(); ig++){   haT2T3norm->AddPoint(ha23dec23->GetPointX(ig), ha23dec23->GetPointY(ig)/normDec23);    haT2T3norm->SetPointError(i23, 0, 0, ha23dec23->GetErrorYlow(ig)/normDec23, ha23dec23->GetErrorYhigh(ig)/normDec23);    i23++; }
    for(int ig=0; ig<ha24dec23->GetN(); ig++){   haT2T4norm->AddPoint(ha24dec23->GetPointX(ig), ha24dec23->GetPointY(ig)/normDec23);    haT2T4norm->SetPointError(i24, 0, 0, ha24dec23->GetErrorYlow(ig)/normDec23, ha24dec23->GetErrorYhigh(ig)/normDec23);    i24++; }
    for(int ig=0; ig<ha34dec23->GetN(); ig++){   haT3T4norm->AddPoint(ha34dec23->GetPointX(ig), ha34dec23->GetPointY(ig)/normDec23);    haT3T4norm->SetPointError(i34, 0, 0, ha34dec23->GetErrorYlow(ig)/normDec23, ha34dec23->GetErrorYhigh(ig)/normDec23);    i34++; }
    
    for(int ig=0; ig<ha12feb24->GetN(); ig++){   haT1T2norm->AddPoint(ha12feb24->GetPointX(ig), ha12feb24->GetPointY(ig)/normFeb24);    haT1T2norm->SetPointError(i12, 0, 0, ha12feb24->GetErrorYlow(ig)/normFeb24, ha12feb24->GetErrorYhigh(ig)/normFeb24);    i12++; }
    for(int ig=0; ig<ha13feb24->GetN(); ig++){   haT1T3norm->AddPoint(ha13feb24->GetPointX(ig), ha13feb24->GetPointY(ig)/normFeb24);    haT1T3norm->SetPointError(i13, 0, 0, ha13feb24->GetErrorYlow(ig)/normFeb24, ha13feb24->GetErrorYhigh(ig)/normFeb24);    i13++; }
    for(int ig=0; ig<ha14feb24->GetN(); ig++){   haT1T4norm->AddPoint(ha14feb24->GetPointX(ig), ha14feb24->GetPointY(ig)/normFeb24);    haT1T4norm->SetPointError(i14, 0, 0, ha14feb24->GetErrorYlow(ig)/normFeb24, ha14feb24->GetErrorYhigh(ig)/normFeb24);    i14++; }
    for(int ig=0; ig<ha23feb24->GetN(); ig++){   haT2T3norm->AddPoint(ha23feb24->GetPointX(ig), ha23feb24->GetPointY(ig)/normFeb24);    haT2T3norm->SetPointError(i23, 0, 0, ha23feb24->GetErrorYlow(ig)/normFeb24, ha23feb24->GetErrorYhigh(ig)/normFeb24);    i23++; }
    for(int ig=0; ig<ha24feb24->GetN(); ig++){   haT2T4norm->AddPoint(ha24feb24->GetPointX(ig), ha24feb24->GetPointY(ig)/normFeb24);    haT2T4norm->SetPointError(i24, 0, 0, ha24feb24->GetErrorYlow(ig)/normFeb24, ha24feb24->GetErrorYhigh(ig)/normFeb24);    i24++; }
    for(int ig=0; ig<ha34feb24->GetN(); ig++){   haT3T4norm->AddPoint(ha34feb24->GetPointX(ig), ha34feb24->GetPointY(ig)/normFeb24);    haT3T4norm->SetPointError(i34, 0, 0, ha34feb24->GetErrorYlow(ig)/normFeb24, ha34feb24->GetErrorYhigh(ig)/normFeb24);    i34++; }
    
    for(int ig=0; ig<ha12may24->GetN(); ig++){   haT1T2norm->AddPoint(ha12may24->GetPointX(ig), ha12may24->GetPointY(ig)/normMay24);    haT1T2norm->SetPointError(i12, 0, 0, ha12may24->GetErrorYlow(ig)/normMay24, ha12may24->GetErrorYhigh(ig)/normMay24);    i12++; }
    for(int ig=0; ig<ha13may24->GetN(); ig++){   haT1T3norm->AddPoint(ha13may24->GetPointX(ig), ha13may24->GetPointY(ig)/normMay24);    haT1T3norm->SetPointError(i13, 0, 0, ha13may24->GetErrorYlow(ig)/normMay24, ha13may24->GetErrorYhigh(ig)/normMay24);    i13++; }
    for(int ig=0; ig<ha14may24->GetN(); ig++){   haT1T4norm->AddPoint(ha14may24->GetPointX(ig), ha14may24->GetPointY(ig)/normMay24);    haT1T4norm->SetPointError(i14, 0, 0, ha14may24->GetErrorYlow(ig)/normMay24, ha14may24->GetErrorYhigh(ig)/normMay24);    i14++; }
    for(int ig=0; ig<ha23may24->GetN(); ig++){   haT2T3norm->AddPoint(ha23may24->GetPointX(ig), ha23may24->GetPointY(ig)/normMay24);    haT2T3norm->SetPointError(i23, 0, 0, ha23may24->GetErrorYlow(ig)/normMay24, ha23may24->GetErrorYhigh(ig)/normMay24);    i23++; }
    for(int ig=0; ig<ha24may24->GetN(); ig++){   haT2T4norm->AddPoint(ha24may24->GetPointX(ig), ha24may24->GetPointY(ig)/normMay24);    haT2T4norm->SetPointError(i24, 0, 0, ha24may24->GetErrorYlow(ig)/normMay24, ha24may24->GetErrorYhigh(ig)/normMay24);    i24++; }
    for(int ig=0; ig<ha34may24->GetN(); ig++){   haT3T4norm->AddPoint(ha34may24->GetPointX(ig), ha34may24->GetPointY(ig)/normMay24);    haT3T4norm->SetPointError(i34, 0, 0, ha34may24->GetErrorYlow(ig)/normMay24, ha34may24->GetErrorYhigh(ig)/normMay24);    i34++; }
    
    // save hour angle data to root file
    outfile->cd();  outfile->mkdir("hourAngles");   outfile->cd("hourAngles");
    hrangleT1T2->Write("haT1T2");
    hrangleT1T3->Write("haT1T3");
    hrangleT1T4->Write("haT1T4");
    hrangleT2T3->Write("haT2T3");
    hrangleT2T4->Write("haT2T4");
    hrangleT3T4->Write("haT3T4");
    outfile->cd();
    outfile->mkdir("normHourAngles");   outfile->cd("normHourAngles");
    haT1T2norm->Write("normhaT1T2");
    haT1T3norm->Write("normhaT1T3");
    haT1T4norm->Write("normhaT1T4");
    haT2T3norm->Write("normhaT2T3");
    haT2T4norm->Write("normhaT2T4");
    haT3T4norm->Write("normhaT3T4");
    outfile->cd();
    outfile->Close();
    
    
    // draw uv clumps all together
    TLegend* uvleg = new TLegend(0.85,0.7,0.9, 0.9);
    uvleg->AddEntry(uv55m, "50-60m","P");
    uvleg->AddEntry(uv65m, "60-68m","P");
    uvleg->AddEntry(uv75m, "68-78m","P");
    uvleg->AddEntry(uv80m, "78-83m","P");
    uvleg->AddEntry(uv85m, "83-88m","P");
    uvleg->AddEntry(uv90m, "88-95m","P");
    uvleg->AddEntry(uv100m, "95-103m","P");
    uvleg->AddEntry(uv110m, "103m+","P");
    
    cuvClumps->cd();    
    //TAxis* uvXaxis = uvaxes->GetXaxis();    uvXaxis->SetNdivisions(500);
    //uvXaxis->ChangeLabel(1,-1,-1,-1,-1,-1,"-#pi/4");     uvXaxis->ChangeLabel(-1,-1,-1,-1,-1,-1,"#pi/4");
    uvaxes->AddPoint(-pi/2.0, 16e-6);  uvaxes->AddPoint(pi/2.0, -2e-6);   uvaxes->Draw("AP");
    uv55m->SetMarkerColor(2);   uv55m->Draw("PSAME");   
    uv65m->SetMarkerColor(94);  uv65m->Draw("PSAME");
    uv75m->SetMarkerColor(5);   uv75m->Draw("PSAME");
    uv80m->SetMarkerColor(77);  uv80m->Draw("PSAME");
    uv85m->SetMarkerColor(7);   uv85m->Draw("PSAME");
    uv90m->SetMarkerColor(63);  uv90m->Draw("PSAME");
    uv100m->SetMarkerColor(51); uv100m->Draw("PSAME");
    uv110m->SetMarkerColor(6);  uv110m->Draw("PSAME");
    uvleg->Draw("SAME");
    cuvClumps->Print("uvClumps.pdf");
    
    // make axes all the same
    uv55m->AddPoint(-2,-6e-6);      uv55m->AddPoint(2,20e-6);   uv55m->GetXaxis()->SetRangeUser(-pi/2.0, pi/2.0);   uv55m->GetYaxis()->SetRangeUser(-4e-6,18e-6);
    uv65m->AddPoint(-2,-6e-6);      uv65m->AddPoint(2,20e-6);   uv65m->GetXaxis()->SetRangeUser(-pi/2.0, pi/2.0);   uv65m->GetYaxis()->SetRangeUser(-4e-6,18e-6);
    uv75m->AddPoint(-2,-6e-6);      uv75m->AddPoint(2,20e-6);   uv75m->GetXaxis()->SetRangeUser(-pi/2.0, pi/2.0);   uv75m->GetYaxis()->SetRangeUser(-4e-6,18e-6);
    uv80m->AddPoint(-2,-6e-6);      uv80m->AddPoint(2,20e-6);   uv80m->GetXaxis()->SetRangeUser(-pi/2.0, pi/2.0);   uv80m->GetYaxis()->SetRangeUser(-4e-6,18e-6);
    uv85m->AddPoint(-2,-6e-6);      uv85m->AddPoint(2,20e-6);   uv85m->GetXaxis()->SetRangeUser(-pi/2.0, pi/2.0);   uv85m->GetYaxis()->SetRangeUser(-4e-6,18e-6);
    uv90m->AddPoint(-2,-6e-6);      uv90m->AddPoint(2,20e-6);   uv90m->GetXaxis()->SetRangeUser(-pi/2.0, pi/2.0);   uv90m->GetYaxis()->SetRangeUser(-4e-6,18e-6);
    uv100m->AddPoint(-2,-6e-6);     uv100m->AddPoint(2,20e-6);  uv100m->GetXaxis()->SetRangeUser(-pi/2.0, pi/2.0);  uv100m->GetYaxis()->SetRangeUser(-4e-6,18e-6);
    uv110m->AddPoint(-2,-6e-6);     uv110m->AddPoint(2,20e-6);  uv110m->GetXaxis()->SetRangeUser(-pi/2.0, pi/2.0);  uv110m->GetYaxis()->SetRangeUser(-4e-6,18e-6);
    
    // draw uv vs vis plots grouped by baseline
    cuvClumps->Clear();     cuvClumps->Divide(4,2);
    cuvClumps->cd(1);   uv55m->Draw("AP");   new55m->Draw("PSAME");// set range user seems to only work if the range is SMALLER than the existing points
    cuvClumps->cd(2);   uv65m->Draw("AP");   new65m->Draw("PSAME");
    cuvClumps->cd(3);   uv75m->Draw("AP");   new75m->Draw("PSAME");
    cuvClumps->cd(4);   uv80m->Draw("AP");   new80m->Draw("PSAME");
    cuvClumps->cd(5);   uv85m->Draw("AP");   new85m->Draw("PSAME");
    cuvClumps->cd(6);   uv90m->Draw("AP");   new90m->Draw("PSAME");
    cuvClumps->cd(7);   uv100m->Draw("AP");  new100m->Draw("PSAME");
    cuvClumps->cd(8);   uv110m->Draw("AP");  new110m->Draw("PSAME");
    cuvClumps->Print("uvClumps.pdf");
    //cuvClumps->Print("uvClumps.pdf]");
    
    // --------------------------------- send info for hour angle calc to python ----------------------------------------------------- (just do it again bc I dont have a full 24 hrs for this data)
    
    TString datetime = "2024-02-20 12:00:00"; // start time - remember converted to UTC time in python - at local noon for now, so center = midnight (doesn't matter when)
    int length = 86400; // number of seconds in a day - length of time to calculate data for in seconds
    int inc = 60; // number of seconds per increment 
    const int nInc = length/inc; // number of data points = length/increment size
    
    // Calculate hour angle here first - I guess execute the python from here for some parameters?
    TString runpy;
    TString py = "python3 /users/PAS1977/jrose/macros/calcDelays/CalcMvt_hourangle.py";
    TString lengthStr, incStr, nIncStr;
    lengthStr = to_string(length);
    incStr = to_string(inc);
    nIncStr = to_string(nInc);
    
    runpy = py + " " + lengthStr + " " + nIncStr + " " + incStr + " " + source + " " + datetime;
    //86400 1440 60 gam cas 2024-02-20 12:00:00
    cout << "Running python script to calculate hour angles over all times:" << endl;
    gSystem->Exec(runpy);
    
    // read in hour angle info from python
    ifstream pyfile12, pyfile13, pyfile14, pyfile23, pyfile24, pyfile34;
    pyfile12.open("pyinfoT1T2.txt");
    pyfile13.open("pyinfoT1T3.txt");
    pyfile14.open("pyinfoT1T4.txt");
    pyfile23.open("pyinfoT2T3.txt");
    pyfile24.open("pyinfoT2T4.txt");
    pyfile34.open("pyinfoT3T4.txt");
    
    double frame[nInc], itime[nInc];
    double u12[nInc], u13[nInc], u14[nInc], u23[nInc], u24[nInc], u34[nInc];
    double v12[nInc], v13[nInc], v14[nInc], v23[nInc], v24[nInc], v34[nInc];
    double baseline12[nInc], baseline13[nInc], baseline14[nInc], baseline23[nInc], baseline24[nInc], baseline34[nInc];
    double opd12[nInc], opd13[nInc], opd14[nInc], opd23[nInc], opd24[nInc], opd34[nInc];
    double hrang12[nInc], hrang13[nInc], hrang14[nInc], hrang23[nInc], hrang24[nInc], hrang34[nInc];
    double junk;
    
    pyfile12 >> junk >> junk >> junk >> junk >> junk >> junk >> junk;
    pyfile13 >> junk >> junk >> junk >> junk >> junk >> junk >> junk;
    pyfile14 >> junk >> junk >> junk >> junk >> junk >> junk >> junk;
    pyfile23 >> junk >> junk >> junk >> junk >> junk >> junk >> junk;
    pyfile24 >> junk >> junk >> junk >> junk >> junk >> junk >> junk;
    pyfile34 >> junk >> junk >> junk >> junk >> junk >> junk >> junk;
    
    for (int i=0; i<nInc; i++){
        pyfile12 >> frame[i] >> itime[i] >> u12[i] >> v12[i] >> baseline12[i] >> opd12[i] >> hrang12[i];
        pyfile13 >> frame[i] >> itime[i] >> u13[i] >> v13[i] >> baseline13[i] >> opd13[i] >> hrang13[i];
        pyfile14 >> frame[i] >> itime[i] >> u14[i] >> v14[i] >> baseline14[i] >> opd14[i] >> hrang14[i];
        pyfile23 >> frame[i] >> itime[i] >> u23[i] >> v23[i] >> baseline23[i] >> opd23[i] >> hrang23[i];
        pyfile24 >> frame[i] >> itime[i] >> u24[i] >> v24[i] >> baseline24[i] >> opd24[i] >> hrang24[i];
        pyfile34 >> frame[i] >> itime[i] >> u34[i] >> v34[i] >> baseline34[i] >> opd34[i] >> hrang34[i];
    }
    
    // filling data for hour angle model plots based on fitted star size
    for(int i=0; i<nInc; i++){
        visHA12all->AddPoint(hrang12[i], calcVis(baseline12[i],diamAll,normAll));
        visHA13all->AddPoint(hrang13[i], calcVis(baseline13[i],diamAll,normAll));
        visHA14all->AddPoint(hrang14[i], calcVis(baseline14[i],diamAll,normAll));
        visHA23all->AddPoint(hrang23[i], calcVis(baseline23[i],diamAll,normAll));
        visHA24all->AddPoint(hrang24[i], calcVis(baseline24[i],diamAll,normAll));
        visHA34all->AddPoint(hrang34[i], calcVis(baseline34[i],diamAll,normAll));
        
        visHA12jf23->AddPoint(hrang12[i], calcVis(baseline12[i],diamJF23,normJF23));
        visHA13jf23->AddPoint(hrang13[i], calcVis(baseline13[i],diamJF23,normJF23));
        visHA14jf23->AddPoint(hrang14[i], calcVis(baseline14[i],diamJF23,normJF23));
        visHA23jf23->AddPoint(hrang23[i], calcVis(baseline23[i],diamJF23,normJF23));
        visHA24jf23->AddPoint(hrang24[i], calcVis(baseline24[i],diamJF23,normJF23));
        visHA34jf23->AddPoint(hrang34[i], calcVis(baseline34[i],diamJF23,normJF23));
        
        visHA12dec23->AddPoint(hrang12[i], calcVis(baseline12[i],diamDec23,normDec23));
        visHA13dec23->AddPoint(hrang13[i], calcVis(baseline13[i],diamDec23,normDec23));
        visHA14dec23->AddPoint(hrang14[i], calcVis(baseline14[i],diamDec23,normDec23));
        visHA23dec23->AddPoint(hrang23[i], calcVis(baseline23[i],diamDec23,normDec23));
        visHA24dec23->AddPoint(hrang24[i], calcVis(baseline24[i],diamDec23,normDec23));
        visHA34dec23->AddPoint(hrang34[i], calcVis(baseline34[i],diamDec23,normDec23));
        
        visHA12feb24->AddPoint(hrang12[i], calcVis(baseline12[i],diamFeb24,normFeb24));
        visHA13feb24->AddPoint(hrang13[i], calcVis(baseline13[i],diamFeb24,normFeb24));
        visHA14feb24->AddPoint(hrang14[i], calcVis(baseline14[i],diamFeb24,normFeb24));
        visHA23feb24->AddPoint(hrang23[i], calcVis(baseline23[i],diamFeb24,normFeb24));
        visHA24feb24->AddPoint(hrang24[i], calcVis(baseline24[i],diamFeb24,normFeb24));
        visHA34feb24->AddPoint(hrang34[i], calcVis(baseline34[i],diamFeb24,normFeb24));
        
        visHA12may24->AddPoint(hrang12[i], calcVis(baseline12[i],diamMay24,normMay24));
        visHA13may24->AddPoint(hrang13[i], calcVis(baseline13[i],diamMay24,normMay24));
        visHA14may24->AddPoint(hrang14[i], calcVis(baseline14[i],diamMay24,normMay24));
        visHA23may24->AddPoint(hrang23[i], calcVis(baseline23[i],diamMay24,normMay24));
        visHA24may24->AddPoint(hrang24[i], calcVis(baseline24[i],diamMay24,normMay24));
        visHA34may24->AddPoint(hrang34[i], calcVis(baseline34[i],diamMay24,normMay24));
    }
    
    // setting ranges of hour angle plots manually - adding extra points so we can use SetRangeUser to make axes the same for all
    //hrangleT1T2->AddPoint(-13,-6e-6);    hrangleT1T2->AddPoint(13,20e-6);    hrangleT1T2->GetXaxis()->SetRangeUser(-12.5,12.5);   hrangleT1T2->GetYaxis()->SetRangeUser(-4e-6,12e-6);
    //hrangleT1T3->AddPoint(-13,-6e-6);    hrangleT1T3->AddPoint(13,20e-6);    hrangleT1T3->GetXaxis()->SetRangeUser(-12.5,12.5);   hrangleT1T3->GetYaxis()->SetRangeUser(-2e-6,11e-6);
    //hrangleT1T4->AddPoint(-13,-6e-6);    hrangleT1T4->AddPoint(13,20e-6);    hrangleT1T4->GetXaxis()->SetRangeUser(-12.5,12.5);   hrangleT1T4->GetYaxis()->SetRangeUser(-5e-6,6e-6);
    //hrangleT2T3->AddPoint(-13,-6e-6);    hrangleT2T3->AddPoint(13,20e-6);    hrangleT2T3->GetXaxis()->SetRangeUser(-12.5,12.5);   hrangleT2T3->GetYaxis()->SetRangeUser(4e-6,14e-6);
    //hrangleT2T4->AddPoint(-13,-6e-6);    hrangleT2T4->AddPoint(13,20e-6);    hrangleT2T4->GetXaxis()->SetRangeUser(-12.5,12.5);   hrangleT2T4->GetYaxis()->SetRangeUser(3e-6,19e-6);
    //hrangleT3T4->AddPoint(-13,-6e-6);    hrangleT3T4->AddPoint(13,20e-6);    hrangleT3T4->GetXaxis()->SetRangeUser(-12.5,12.5);   hrangleT3T4->GetYaxis()->SetRangeUser(2e-6,18e-6);
    
    ha12jf23->AddPoint(-13,-5e-6);    ha12jf23->AddPoint(13,17e-6);    ha12jf23->GetXaxis()->SetRangeUser(-12.5,12.5);   ha12jf23->GetYaxis()->SetRangeUser(-4.5e-6,15.5e-6);
    ha13jf23->AddPoint(-13,-5e-6);    ha13jf23->AddPoint(13,17e-6);    ha13jf23->GetXaxis()->SetRangeUser(-12.5,12.5);   ha13jf23->GetYaxis()->SetRangeUser(-4.5e-6,15.5e-6);
    ha14jf23->AddPoint(-13,-5e-6);    ha14jf23->AddPoint(13,17e-6);    ha14jf23->GetXaxis()->SetRangeUser(-12.5,12.5);   ha14jf23->GetYaxis()->SetRangeUser(-4.5e-6,15.5e-6);
    ha23jf23->AddPoint(-13,-5e-6);    ha23jf23->AddPoint(13,17e-6);    ha23jf23->GetXaxis()->SetRangeUser(-12.5,12.5);   ha23jf23->GetYaxis()->SetRangeUser(-4.5e-6,15.5e-6);
    ha24jf23->AddPoint(-13,-5e-6);    ha24jf23->AddPoint(13,17e-6);    ha24jf23->GetXaxis()->SetRangeUser(-12.5,12.5);   ha24jf23->GetYaxis()->SetRangeUser(-4.5e-6,15.5e-6);
    ha34jf23->AddPoint(-13,-5e-6);    ha34jf23->AddPoint(13,17e-6);    ha34jf23->GetXaxis()->SetRangeUser(-12.5,12.5);   ha34jf23->GetYaxis()->SetRangeUser(-4.5e-6,15.5e-6);
    
    ha12dec23->AddPoint(-13,-5e-6);    ha12dec23->AddPoint(13,17e-6);    ha12dec23->GetXaxis()->SetRangeUser(-12.5,12.5);   ha12dec23->GetYaxis()->SetRangeUser(-4e-6,15e-6);
    ha13dec23->AddPoint(-13,-5e-6);    ha13dec23->AddPoint(13,17e-6);    ha13dec23->GetXaxis()->SetRangeUser(-12.5,12.5);   ha13dec23->GetYaxis()->SetRangeUser(-4e-6,15e-6);
    ha14dec23->AddPoint(-13,-5e-6);    ha14dec23->AddPoint(13,17e-6);    ha14dec23->GetXaxis()->SetRangeUser(-12.5,12.5);   ha14dec23->GetYaxis()->SetRangeUser(-4e-6,15e-6);
    ha23dec23->AddPoint(-13,-5e-6);    ha23dec23->AddPoint(13,17e-6);    ha23dec23->GetXaxis()->SetRangeUser(-12.5,12.5);   ha23dec23->GetYaxis()->SetRangeUser(-4e-6,15e-6);
    ha24dec23->AddPoint(-13,-5e-6);    ha24dec23->AddPoint(13,17e-6);    ha24dec23->GetXaxis()->SetRangeUser(-12.5,12.5);   ha24dec23->GetYaxis()->SetRangeUser(-4e-6,15e-6);
    ha34dec23->AddPoint(-13,-5e-6);    ha34dec23->AddPoint(13,17e-6);    ha34dec23->GetXaxis()->SetRangeUser(-12.5,12.5);   ha34dec23->GetYaxis()->SetRangeUser(-4e-6,15e-6);
    
    ha12feb24->AddPoint(-13,-6e-6);    ha12feb24->AddPoint(13,20e-6);    ha12feb24->GetXaxis()->SetRangeUser(-12.5,12.5);   ha12feb24->GetYaxis()->SetRangeUser(-5e-6,19e-6);
    ha13feb24->AddPoint(-13,-6e-6);    ha13feb24->AddPoint(13,20e-6);    ha13feb24->GetXaxis()->SetRangeUser(-12.5,12.5);   ha13feb24->GetYaxis()->SetRangeUser(-5e-6,19e-6);
    ha14feb24->AddPoint(-13,-6e-6);    ha14feb24->AddPoint(13,20e-6);    ha14feb24->GetXaxis()->SetRangeUser(-12.5,12.5);   ha14feb24->GetYaxis()->SetRangeUser(-5e-6,19e-6);
    ha23feb24->AddPoint(-13,-6e-6);    ha23feb24->AddPoint(13,20e-6);    ha23feb24->GetXaxis()->SetRangeUser(-12.5,12.5);   ha23feb24->GetYaxis()->SetRangeUser(-5e-6,19e-6);
    ha24feb24->AddPoint(-13,-6e-6);    ha24feb24->AddPoint(13,20e-6);    ha24feb24->GetXaxis()->SetRangeUser(-12.5,12.5);   ha24feb24->GetYaxis()->SetRangeUser(-5e-6,19e-6);
    ha34feb24->AddPoint(-13,-6e-6);    ha34feb24->AddPoint(13,20e-6);    ha34feb24->GetXaxis()->SetRangeUser(-12.5,12.5);   ha34feb24->GetYaxis()->SetRangeUser(-5e-6,19e-6);
    
    ha12may24->AddPoint(-13,-1e-6);    ha12may24->AddPoint(13,17e-6);    ha12may24->GetXaxis()->SetRangeUser(-12.5,12.5);   ha12may24->GetYaxis()->SetRangeUser(-0.5e-6,16e-6);
    ha13may24->AddPoint(-13,-1e-6);    ha13may24->AddPoint(13,17e-6);    ha13may24->GetXaxis()->SetRangeUser(-12.5,12.5);   ha13may24->GetYaxis()->SetRangeUser(-0.5e-6,16e-6);
    ha14may24->AddPoint(-13,-1e-6);    ha14may24->AddPoint(13,17e-6);    ha14may24->GetXaxis()->SetRangeUser(-12.5,12.5);   ha14may24->GetYaxis()->SetRangeUser(-0.5e-6,16e-6);
    ha23may24->AddPoint(-13,-1e-6);    ha23may24->AddPoint(13,17e-6);    ha23may24->GetXaxis()->SetRangeUser(-12.5,12.5);   ha23may24->GetYaxis()->SetRangeUser(-0.5e-6,16e-6);
    ha24may24->AddPoint(-13,-1e-6);    ha24may24->AddPoint(13,17e-6);    ha24may24->GetXaxis()->SetRangeUser(-12.5,12.5);   ha24may24->GetYaxis()->SetRangeUser(-0.5e-6,16e-6);
    ha34may24->AddPoint(-13,-1e-6);    ha34may24->AddPoint(13,17e-6);    ha34may24->GetXaxis()->SetRangeUser(-12.5,12.5);   ha34may24->GetYaxis()->SetRangeUser(-0.5e-6,16e-6);
    
    // read in data from fpga data file
    ifstream fpgapts;    fpgapts.open("fpgaHrAngles/matchedEverything.txt");
    const int nfpgapts(25);
    string fpgaDirs[nfpgapts];
    string fpgaPairs[nfpgapts];
    double fpgaVis[nfpgapts];
    double fpgaVisErr[nfpgapts];
    double fpgaBaseline[nfpgapts];
    double fpgaBaseErr[nfpgapts];
    double fpgaU[nfpgapts];
    double fpgaV[nfpgapts];
    double fpgaHA[nfpgapts];
    string thisdir;
    string thispair;
    double thisvis, thisviserr, thisbase, thisbaseerr, thisu, thisv, thisha;
    
    TGraphErrors* fpgaVisHA12 = new TGraphErrors;   fpgaVisHA12->SetMarkerStyle(47);    fpgaVisHA12->SetMarkerColor(2);     int nfpga12(0), nfpga13(0), nfpga14(0), nfpga23(0), nfpga24(0), nfpga34(0);
    TGraphErrors* fpgaVisHA13 = new TGraphErrors;   fpgaVisHA13->SetMarkerStyle(47);    fpgaVisHA13->SetMarkerColor(4);
    TGraphErrors* fpgaVisHA14 = new TGraphErrors;   fpgaVisHA14->SetMarkerStyle(47);    fpgaVisHA14->SetMarkerColor(92);
    TGraphErrors* fpgaVisHA23 = new TGraphErrors;   fpgaVisHA23->SetMarkerStyle(47);    fpgaVisHA23->SetMarkerColor(3);
    TGraphErrors* fpgaVisHA24 = new TGraphErrors;   fpgaVisHA24->SetMarkerStyle(47);    fpgaVisHA24->SetMarkerColor(6);
    TGraphErrors* fpgaVisHA34 = new TGraphErrors;   fpgaVisHA34->SetMarkerStyle(47);    fpgaVisHA34->SetMarkerColor(1);
    
    for(int ifpga=0; ifpga<nfpgapts; ifpga++){
        fpgapts >> thisdir >> thispair >> thisvis >> thisviserr >> thisbase >> thisbaseerr >> thisu >> thisv >> thisha;
        fpgaDirs[ifpga] = thisdir;
        fpgaPairs[ifpga] = thispair;
        fpgaVis[ifpga] = thisvis;
        fpgaVisErr[ifpga] = thisviserr;
        fpgaBaseline[ifpga] = thisbase;
        fpgaBaseErr[ifpga] = thisbaseerr;
        fpgaU[ifpga] = thisu;
        fpgaV[ifpga] = thisv;
        fpgaHA[ifpga] = thisha;
        if(thispair == "T1T2"){ fpgaVisHA12->AddPoint(thisha, thisvis);     fpgaVisHA12->SetPointError(nfpga12, 0, thisviserr);     nfpga12++; }
        if(thispair == "T1T3"){ fpgaVisHA13->AddPoint(thisha, thisvis);     fpgaVisHA13->SetPointError(nfpga13, 0, thisviserr);     nfpga13++; }
        if(thispair == "T1T4"){ fpgaVisHA14->AddPoint(thisha, thisvis);     fpgaVisHA14->SetPointError(nfpga14, 0, thisviserr);     nfpga14++; }
        if(thispair == "T2T3"){ fpgaVisHA23->AddPoint(thisha, thisvis);     fpgaVisHA23->SetPointError(nfpga23, 0, thisviserr);     nfpga23++; }
        if(thispair == "T2T4"){ fpgaVisHA24->AddPoint(thisha, thisvis);     fpgaVisHA24->SetPointError(nfpga24, 0, thisviserr);     nfpga24++; }
        if(thispair == "T3T4"){ fpgaVisHA34->AddPoint(thisha, thisvis);     fpgaVisHA34->SetPointError(nfpga34, 0, thisviserr);     nfpga34++; }
    }
    
    // draw hour angle plots
    cuvClumps->Clear();     cuvClumps->Divide(3,2);
    cuvClumps->cd(1);    hrangleT1T2->Draw("AP");       fpgaVisHA12->Draw("PSAME");
    cuvClumps->cd(2);    hrangleT1T3->Draw("AP");       fpgaVisHA13->Draw("PSAME");
    cuvClumps->cd(3);    hrangleT1T4->Draw("AP");       fpgaVisHA14->Draw("PSAME");
    cuvClumps->cd(4);    hrangleT2T3->Draw("AP");       fpgaVisHA23->Draw("PSAME");
    cuvClumps->cd(5);    hrangleT2T4->Draw("AP");       fpgaVisHA24->Draw("PSAME");
    cuvClumps->cd(6);    hrangleT3T4->Draw("AP");       fpgaVisHA34->Draw("PSAME");
    cuvClumps->Print("uvClumps.pdf");
    
    cuvClumps->Clear();     cuvClumps->Divide(3,2);
    cuvClumps->cd(1);   visHA12all->Draw("AP");
    cuvClumps->cd(2);   visHA13all->Draw("AP");
    cuvClumps->cd(3);   visHA14all->Draw("AP");
    cuvClumps->cd(4);   visHA23all->Draw("AP");
    cuvClumps->cd(5);   visHA24all->Draw("AP");
    cuvClumps->cd(6);   visHA34all->Draw("AP");
    cuvClumps->Print("uvClumps.pdf");
    
    // draw all predicted hour angles with visibilities for different fit results of star size
    cuvClumps->Clear();     cuvClumps->Divide(4,2);
    TPaveText* textHAlines = new TPaveText(0.1,0.1,0.9,0.9);
    textHAlines->AddText(source);
    textHAlines->AddText("hour angle vs visibility");
    textHAlines->AddText("for fits to different data groups");
    cuvClumps->cd(1);   textHAlines->Draw();
    
    TLegend* legHAlines = new TLegend(0.1,0.1,0.9,0.9);
    legHAlines->AddEntry(visHA12all,"all","L");
    legHAlines->AddEntry(visHA12jf23,"Jan/Feb 2023","L");
    legHAlines->AddEntry(visHA12dec23,"Dec 2023","L");
    legHAlines->AddEntry(visHA12feb24,"Feb 2024","L");
    legHAlines->AddEntry(visHA12may24,"May 2024","L");
    cuvClumps->cd(5);   legHAlines->Draw();
    
    cuvClumps->cd(2);   visHA12jf23->Draw("AL");    visHA12all->Draw("LSAME");      visHA12feb24->Draw("LSAME");     visHA12dec23->Draw("LSAME");   visHA12may24->Draw("LSAME");
    cuvClumps->cd(3);   visHA13jf23->Draw("AL");    visHA13all->Draw("LSAME");      visHA13feb24->Draw("LSAME");     visHA13dec23->Draw("LSAME");   visHA13may24->Draw("LSAME");
    cuvClumps->cd(4);   visHA14jf23->Draw("AL");    visHA14all->Draw("LSAME");      visHA14feb24->Draw("LSAME");     visHA14dec23->Draw("LSAME");   visHA14may24->Draw("LSAME");
    cuvClumps->cd(6);   visHA23jf23->Draw("AL");    visHA23all->Draw("LSAME");      visHA23feb24->Draw("LSAME");     visHA23dec23->Draw("LSAME");   visHA23may24->Draw("LSAME");
    cuvClumps->cd(7);   visHA24jf23->Draw("AL");    visHA24all->Draw("LSAME");      visHA24feb24->Draw("LSAME");     visHA24dec23->Draw("LSAME");   visHA24may24->Draw("LSAME");
    cuvClumps->cd(8);   visHA34jf23->Draw("AL");    visHA34all->Draw("LSAME");      visHA34feb24->Draw("LSAME");     visHA34dec23->Draw("LSAME");   visHA34may24->Draw("LSAME");
    cuvClumps->Print("uvClumps.pdf");
    
    // new sub-divided canvas plotting code to make very nice canvases (see testCanvasPads.C locally)
    // draw all predicted hour angles with visibilities for different fit results of star size
    cuvClumps->Clear();
    cuvClumps->cd(0);
    TPad *pTitle = new TPad("pTitle","pTitle",0.35,0.9,1.0,1.0);
    TPaveText* textalldata = new TPaveText(0.1,0.1,0.9,0.9);
    textalldata->SetTextFont(53); textalldata->SetTextSize(20);
    textalldata->AddText(source);
    textalldata->AddText("hour angle vs visibility all data, with expectation based on fit to all data");
    pTitle->Draw();
    pTitle->cd();
    textalldata->Draw("NB");
    cuvClumps->cd(0);
    TLegend* legpairs = new TLegend(0.21,0.05,0.3,0.32); //(0.65,0.01,0.99,0.99);
    legpairs->AddEntry(visHA12all,"T1T2","P");
    legpairs->AddEntry(visHA13all,"T1T3","P");
    legpairs->AddEntry(visHA14all,"T1T4","P");
    legpairs->AddEntry(visHA23all,"T2T3","P");
    legpairs->AddEntry(visHA24all,"T2T4","P");
    legpairs->AddEntry(visHA34all,"T3T4","P");
    legpairs->Draw();
    TLegend* legHApts = new TLegend(0.05,0.05,0.21,0.32); //(0.01,0.01,0.65,0.99);
    legHApts->AddEntry(ha12jf23,"Jan/Feb 2023","P");
    legHApts->AddEntry(ha12dec23,"Dec 2023","P");
    legHApts->AddEntry(ha12feb24,"Feb 2024","P");
    legHApts->AddEntry(ha12may24,"May 2024","P");
    legHApts->Draw();   
    cuvClumps->cd(0);
    TPad *pLvc = new TPad("pLvc","pLvc",0.0,0.35,0.35,1.0);
    pLvc->Draw();
    pLvc->cd();
    allpts->DrawClone("AP");
    zero->Draw("SAME");
    //viscurveJanFeb2023->Draw("PSAME");
    //viscurveFeb2024->Draw("PSAME");
    //viscurveMay2024->Draw("PSAME");
    cuvClumps->cd(0);
    TPad *pRight = new TPad("pRight","pRight",0.35,0.0,1.0,0.9);
    pRight->Divide(3,2);
    pRight->Draw();
    TPad *pR12 = (TPad*)pRight->cd(1);  pR12->Draw();   hrangleT1T2->Draw("AP");    visHA12all->Draw("PSAME");      hrangleT1T2->Draw("PSAME");
    TPad *pR13 = (TPad*)pRight->cd(2);  pR13->Draw();   hrangleT1T3->Draw("AP");    visHA13all->Draw("PSAME");      hrangleT1T3->Draw("PSAME");
    TPad *pR14 = (TPad*)pRight->cd(3);  pR14->Draw();   hrangleT1T4->Draw("AP");    visHA14all->Draw("PSAME");      hrangleT1T4->Draw("PSAME");
    TPad *pR23 = (TPad*)pRight->cd(4);  pR23->Draw();   hrangleT2T3->Draw("AP");    visHA23all->Draw("PSAME");      hrangleT2T3->Draw("PSAME");
    TPad *pR24 = (TPad*)pRight->cd(5);  pR24->Draw();   hrangleT2T4->Draw("AP");    visHA24all->Draw("PSAME");      hrangleT2T4->Draw("PSAME");
    TPad *pR34 = (TPad*)pRight->cd(6);  pR34->Draw();   hrangleT3T4->Draw("AP");    visHA34all->Draw("PSAME");      hrangleT3T4->Draw("PSAME");
    cuvClumps->Print("uvClumps.pdf");
    
    // draw all data sets on expected curve from overall fit
    cuvClumps->cd(0);
    TPaveText* texteachgroup = new TPaveText(0.1,0.1,0.9,0.9);
    texteachgroup->SetTextFont(53); texteachgroup->SetTextSize(10);
    texteachgroup->AddText(source);
    texteachgroup->AddText("curve from overall fit with data from each observation group");
    pTitle->Draw();     pTitle->cd();   textalldata->Draw("NB");
    //TPad *pLvc = new TPad("pLvc","pLvc",0.0,0.35,0.35,1.0);
    cuvClumps->cd(0);
    pLvc->Draw();   pLvc->cd();
    allpts->DrawClone("AP");
    zero->Draw("SAME");
    viscurveJanFeb2023->Draw("PSAME");
    viscurveDec2023->Draw("PSAME");
    viscurveFeb2024->Draw("PSAME");
    viscurveMay2024->Draw("PSAME");
    cuvClumps->cd(0);
    legpairs->Draw();   legHApts->Draw();
    cuvClumps->cd(0);
    pRight->Draw();
    (TPad*)pRight->cd(1);   pR12->Draw();   hrangleT1T2->Draw("AP");    visHA12all->Draw("PSAME");      ha12jf23->Draw("PSAME");    ha12dec23->Draw("PSAME");   ha12feb24->Draw("PSAME");      ha12may24->Draw("PSAME");
    (TPad*)pRight->cd(2);   pR13->Draw();   hrangleT1T3->Draw("AP");    visHA13all->Draw("PSAME");      ha13jf23->Draw("PSAME");    ha13dec23->Draw("PSAME");   ha13feb24->Draw("PSAME");      ha13may24->Draw("PSAME");
    (TPad*)pRight->cd(3);   pR14->Draw();   hrangleT1T4->Draw("AP");    visHA14all->Draw("PSAME");      ha14jf23->Draw("PSAME");    ha14dec23->Draw("PSAME");   ha14feb24->Draw("PSAME");      ha14may24->Draw("PSAME");
    (TPad*)pRight->cd(4);   pR23->Draw();   hrangleT2T3->Draw("AP");    visHA23all->Draw("PSAME");      ha23jf23->Draw("PSAME");    ha23dec23->Draw("PSAME");   ha23feb24->Draw("PSAME");      ha23may24->Draw("PSAME");
    (TPad*)pRight->cd(5);   pR24->Draw();   hrangleT2T4->Draw("AP");    visHA24all->Draw("PSAME");      ha24jf23->Draw("PSAME");    ha24dec23->Draw("PSAME");   ha24feb24->Draw("PSAME");      ha24may24->Draw("PSAME");
    (TPad*)pRight->cd(6);   pR34->Draw();   hrangleT3T4->Draw("AP");    visHA34all->Draw("PSAME");      ha34jf23->Draw("PSAME");    ha34dec23->Draw("PSAME");   ha34feb24->Draw("PSAME");      ha34may24->Draw("PSAME");
    cuvClumps->Print("uvClumps.pdf");
    
    // jan/feb 2023 data only
    cuvClumps->cd(0);
    TPaveText* textJF23 = new TPaveText(0.1,0.1,0.9,0.9);
    textJF23->SetTextFont(53);  textJF23->SetTextSize(20);
    textJF23->AddText("Jan/Feb 2023 data");     
    pTitle->Draw();     pTitle->cd();   textJF23->Draw("NB");
    cuvClumps->cd(0);
    pLvc->Draw();   pLvc->cd();
    viscurveJanFeb2023->Draw("AP");     viscurveJanFeb2023->GetFunction("udfit")->Draw("SAME");
    zero->Draw("SAME");
    cuvClumps->cd(0);
    legpairs->Draw();   legHApts->Draw();
    cuvClumps->cd(0);
    pRight->Draw();
    (TPad*)pRight->cd(1);   pR12->Draw();   ha12jf23->Draw("AP");   visHA12jf23->Draw("PSAME");    ha12jf23->Draw("PSAME");
    (TPad*)pRight->cd(2);   pR13->Draw();   ha13jf23->Draw("AP");   visHA13jf23->Draw("PSAME");    ha13jf23->Draw("PSAME");
    (TPad*)pRight->cd(3);   pR14->Draw();   ha14jf23->Draw("AP");   visHA14jf23->Draw("PSAME");    ha14jf23->Draw("PSAME");
    (TPad*)pRight->cd(4);   pR23->Draw();   ha23jf23->Draw("AP");   visHA23jf23->Draw("PSAME");    ha23jf23->Draw("PSAME");
    (TPad*)pRight->cd(5);   pR24->Draw();   ha24jf23->Draw("AP");   visHA24jf23->Draw("PSAME");    ha24jf23->Draw("PSAME");
    (TPad*)pRight->cd(6);   pR34->Draw();   ha34jf23->Draw("AP");   visHA34jf23->Draw("PSAME");    ha34jf23->Draw("PSAME");
    cuvClumps->Print("uvClumps.pdf");
    
    cuvClumps->cd(0);
    TPaveText* textDec23 = new TPaveText(0.1,0.1,0.9,0.9);
    textDec23->SetTextFont(53);  textDec23->SetTextSize(20);
    textDec23->AddText("Dec 2023 data");     
    pTitle->Draw();     pTitle->cd();   textDec23->Draw("NB");
    cuvClumps->cd(0);
    pLvc->Draw();   pLvc->cd();
    viscurveDec2023->Draw("AP");     viscurveDec2023->GetFunction("udfit")->Draw("SAME");
    zero->Draw("SAME");
    cuvClumps->cd(0);
    legpairs->Draw();   legHApts->Draw();
    cuvClumps->cd(0);
    pRight->Draw();
    (TPad*)pRight->cd(1);   pR12->Draw();   ha12dec23->Draw("AP");   visHA12dec23->Draw("PSAME");    ha12dec23->Draw("PSAME");
    (TPad*)pRight->cd(2);   pR13->Draw();   ha13dec23->Draw("AP");   visHA13dec23->Draw("PSAME");    ha13dec23->Draw("PSAME");
    (TPad*)pRight->cd(3);   pR14->Draw();   ha14dec23->Draw("AP");   visHA14dec23->Draw("PSAME");    ha14dec23->Draw("PSAME");
    (TPad*)pRight->cd(4);   pR23->Draw();   ha23dec23->Draw("AP");   visHA23dec23->Draw("PSAME");    ha23dec23->Draw("PSAME");
    (TPad*)pRight->cd(5);   pR24->Draw();   ha24dec23->Draw("AP");   visHA24dec23->Draw("PSAME");    ha24dec23->Draw("PSAME");
    (TPad*)pRight->cd(6);   pR34->Draw();   ha34dec23->Draw("AP");   visHA34dec23->Draw("PSAME");    ha34dec23->Draw("PSAME");
    cuvClumps->Print("uvClumps.pdf");
    
    cuvClumps->cd(0);
    TPaveText* textFeb24 = new TPaveText(0.1,0.1,0.9,0.9);
    textFeb24->SetTextFont(53);  textFeb24->SetTextSize(20);
    textFeb24->AddText("Feb 2024 data");     
    pTitle->Draw();     pTitle->cd();   textFeb24->Draw("NB");
    cuvClumps->cd(0);
    pLvc->Draw();   pLvc->cd();
    viscurveFeb2024->Draw("AP");     viscurveFeb2024->GetFunction("udfit")->Draw("SAME");
    zero->Draw("SAME");
    cuvClumps->cd(0);
    legpairs->Draw();   legHApts->Draw();
    cuvClumps->cd(0);
    pRight->Draw();
    (TPad*)pRight->cd(1);   pR12->Draw();   ha12feb24->Draw("AP");   visHA12feb24->Draw("PSAME");    ha12feb24->Draw("PSAME");
    (TPad*)pRight->cd(2);   pR13->Draw();   ha13feb24->Draw("AP");   visHA13feb24->Draw("PSAME");    ha13feb24->Draw("PSAME");
    (TPad*)pRight->cd(3);   pR14->Draw();   ha14feb24->Draw("AP");   visHA14feb24->Draw("PSAME");    ha14feb24->Draw("PSAME");
    (TPad*)pRight->cd(4);   pR23->Draw();   ha23feb24->Draw("AP");   visHA23feb24->Draw("PSAME");    ha23feb24->Draw("PSAME");
    (TPad*)pRight->cd(5);   pR24->Draw();   ha24feb24->Draw("AP");   visHA24feb24->Draw("PSAME");    ha24feb24->Draw("PSAME");
    (TPad*)pRight->cd(6);   pR34->Draw();   ha34feb24->Draw("AP");   visHA34feb24->Draw("PSAME");    ha34feb24->Draw("PSAME");
    cuvClumps->Print("uvClumps.pdf");
    
    cuvClumps->cd(0);
    TPaveText* textMay24 = new TPaveText(0.1,0.1,0.9,0.9);
    textMay24->SetTextFont(53);  textMay24->SetTextSize(20);
    textMay24->AddText("May 2024 data");     
    pTitle->Draw();     pTitle->cd();   textMay24->Draw("NB");
    cuvClumps->cd(0);
    pLvc->Draw();   pLvc->cd();
    viscurveMay2024->Draw("AP");     viscurveMay2024->GetFunction("udfit")->Draw("SAME");
    zero->Draw("SAME");
    cuvClumps->cd(0);
    legpairs->Draw();   legHApts->Draw();
    cuvClumps->cd(0);
    pRight->Draw();
    (TPad*)pRight->cd(1);   pR12->Draw();   ha12may24->Draw("AP");   visHA12may24->Draw("PSAME");    ha12may24->Draw("PSAME");
    (TPad*)pRight->cd(2);   pR13->Draw();   ha13may24->Draw("AP");   visHA13may24->Draw("PSAME");    ha13may24->Draw("PSAME");
    (TPad*)pRight->cd(3);   pR14->Draw();   ha14may24->Draw("AP");   visHA14may24->Draw("PSAME");    ha14may24->Draw("PSAME");
    (TPad*)pRight->cd(4);   pR23->Draw();   ha23may24->Draw("AP");   visHA23may24->Draw("PSAME");    ha23may24->Draw("PSAME");
    (TPad*)pRight->cd(5);   pR24->Draw();   ha24may24->Draw("AP");   visHA24may24->Draw("PSAME");    ha24may24->Draw("PSAME");
    (TPad*)pRight->cd(6);   pR34->Draw();   ha34may24->Draw("AP");   visHA34may24->Draw("PSAME");    ha34may24->Draw("PSAME");
    cuvClumps->Print("uvClumps.pdf");
    
    cuvClumps->cd(0);
    TPaveText* texthatogether = new TPaveText(0.1,0.1,0.9,0.9);
    texthatogether->SetTextFont(53);    texthatogether->SetTextSize(20);
    texthatogether->AddText("hour angles, all months together, colored by pair");
    texthatogether->AddText("note that points below 0 are cut off!");
    pTitle->Draw();     pTitle->cd();       texthatogether->Draw("NB");
    cuvClumps->cd(0);
    pLvc->Draw();   pLvc->cd();
    allpts->DrawClone("AP");     zero->Draw("SAME");
    vcjf23T1T2->Draw("PSAME");  vcjf23T1T3->Draw("PSAME");  vcjf23T1T4->Draw("PSAME");  vcjf23T2T3->Draw("PSAME");  vcjf23T2T4->Draw("PSAME");  vcjf23T3T4->Draw("PSAME");
    vcdec23T1T2->Draw("PSAME");  vcdec23T1T3->Draw("PSAME");  vcdec23T1T4->Draw("PSAME");  vcdec23T2T3->Draw("PSAME");  vcdec23T2T4->Draw("PSAME");  vcdec23T3T4->Draw("PSAME");
    vcfeb24T1T2->Draw("PSAME");  vcfeb24T1T3->Draw("PSAME");  vcfeb24T1T4->Draw("PSAME");  vcfeb24T2T3->Draw("PSAME");  vcfeb24T2T4->Draw("PSAME");  vcfeb24T3T4->Draw("PSAME");
    vcmay24T1T2->Draw("PSAME");  vcmay24T1T3->Draw("PSAME");  vcmay24T1T4->Draw("PSAME");  vcmay24T2T3->Draw("PSAME");  vcmay24T2T4->Draw("PSAME");  vcmay24T3T4->Draw("PSAME");
    cuvClumps->cd(0);
    pRight->Draw();     pRight->cd();
    //TPad *pRight2 = new TPad("pRight2","pRight2",0.35,0.0,1.0,1.0);     pRight2->Draw();    pRight2->cd();
    TGraph* haaxes = new TGraph;    haaxes->SetMarkerStyle(9);  haaxes->SetMarkerColor(0);  haaxes->SetTitle(";hour angle (hrs);area under peak (ns)");
    haaxes->AddPoint(-6, allpts->GetHistogram()->GetMaximum());     haaxes->AddPoint(6, allpts->GetHistogram()->GetMinimum());      haaxes->GetYaxis()->SetRangeUser(0,18e-6);    haaxes->Draw("AP");
    visHA12all->Draw("PSAME");      
    hrangleT1T2->Draw("PSAME");    
    visHA13all->Draw("PSAME");      
    hrangleT1T3->Draw("PSAME");    
    visHA14all->Draw("PSAME");      
    hrangleT1T4->Draw("PSAME");    
    visHA23all->Draw("PSAME");      
    hrangleT2T3->Draw("PSAME");    
    visHA24all->Draw("PSAME");      
    hrangleT2T4->Draw("PSAME");    
    visHA34all->Draw("PSAME");      
    hrangleT3T4->Draw("PSAME");    
    cuvClumps->Print("uvClumps.pdf");
    
    
    
    cuvClumps->Print("uvClumps.pdf]");
    //chrangle->Print("hourAngles.pdf");
    //chrangle->Print("hourAngles.pdf]");
    
    cout << "max all " << allpts->GetHistogram()->GetMaximum() << "   " << allpts->GetHistogram()->GetMinimum() << endl; 
    cout << "jf 23 " << viscurveJanFeb2023->GetHistogram()->GetMaximum() << "   " << viscurveJanFeb2023->GetHistogram()->GetMinimum() << endl;
    cout << "dec 23 " << viscurveDec2023->GetHistogram()->GetMaximum() << "   " << viscurveDec2023->GetHistogram()->GetMinimum() << endl;
    cout << "feb 24 " << viscurveFeb2024->GetHistogram()->GetMaximum() << "   " << viscurveFeb2024->GetHistogram()->GetMinimum() << endl;
    cout << "may 24 " << viscurveMay2024->GetHistogram()->GetMaximum() << "   " << viscurveMay2024->GetHistogram()->GetMinimum() << endl;
 
 cout << "       diam      norm" << endl;
 cout << "all    " << diamAll << "  " << normAll << endl;
 cout << "jf23   " << diamJF23 << "  " << normJF23 << endl;
 cout << "dec23  " << diamDec23 << "  " << normDec23 << endl;
 cout << "feb24  " << diamFeb24 << "  " << normFeb24 << endl;
 cout << "may24  " << diamMay24 << "  " << normMay24 << endl;
 cout << "feb/may24  " << diamFebMay24 << "   " << normFebMay24 << endl;
 cout << "no jan  " << diamNoJan << "   " << normNoJan << endl;
 cout << "no feb  " << diamNoFeb << "   " << normNoFeb << endl;
 cout << "dec/may  " << diamDecMay << "   " << normDecMay << endl;
    
} // end of macro


// functions
// =====================================================================================================================================================================

// calculate excess errorbars that come from RMS of fits to simulated peaks outside peak region
double calcExcessErr(TProfile* cf, double area, double sigma, int reps){
    
    //TH1::AddDirectory(kFALSE);
    
    //TCanvas* tempwindows = new TCanvas;     tempwindows->Divide(6,3);       //int iwindow(0);
    //tempwindows->Print("tempwindows.pdf["); 
    
    double simarea = area;//cf->GetFunction("hbtFit")->GetParameter(0);
    //double simsigma = sigma;
    
    // set up parameters for simulation peaks
    double simtau(0.0), simsigma(4.0);
    double taumin(-10), taumax(10), sigmamin(2.0), sigmamax(7.0);
    int nwindows(16), windowsize(40); // nwindows(18)
    
    // fit function
    TF1* areafit = new TF1("areafit","([0]*exp(-pow((x-[1])/[2],2)/2.0))/([2]*sqrt(2*pi))",-300,300);
    areafit->SetParName(0,"A");          areafit->SetParameter(0,0.0);
    areafit->SetParName(1,"#tau_{0}");   areafit->SetParameter(1, simtau);        areafit->SetParLimits(1, taumin, taumax);
    areafit->SetParName(2,"#sigma");     areafit->SetParameter(2, simsigma);      areafit->SetParLimits(2, sigmamin, sigmamax);
    
    TH1D* areadist = new TH1D("areadist",Form("difference between fit and real area, area = %f ns;input area - fitted area;counts", simarea), 50, -1e-5, 1e-5);
    
    // big loop to run "simulation peaks" more times and fill dist with more points --> IS THIS NEEDED?? you get about the same value only doing it once over each window, but it fluctuates
    for(int i=0; i<reps; i++){
        int startbin(21), endbin(40); // startbin(1)
        
        // loop over off peak region "windows"
        for(int iw=0; iw<nwindows; iw++){
            
            TH1D* fakepeak = simpeak(simarea, simtau, simsigma, windowsize);
            
            int flag(0);
            int thisbin = startbin;
            TH1D* thiswindow = new TH1D("thiswindow",Form("window %d",iw+1), windowsize, -(windowsize/2), (windowsize/2));
            
            // fill the window with off-peak noise plus fake peak
            for(int ix=1; ix<=windowsize; ix++){
                thiswindow->SetBinContent(ix, cf->GetBinContent(thisbin) + fakepeak->GetBinContent(ix));
                thiswindow->SetBinError(ix, cf->GetBinError(thisbin));
                thisbin++;
            }
            //tempwindows->cd(iw+1);
            thiswindow->Draw();
            
            //hbtFit->SetParameter(0,0.0);        hbtFit->SetParameter(1, simtau);        hbtFit->SetParameter(2, simsigma);
            thiswindow->Fit("areafit","Q"); //cout << areafit->GetParameter(0) << endl;
            
            // quality cuts as if this was a real peak in the analysis
            if(abs(areafit->GetParameter(1) - taumin) < 0.001 || abs(areafit->GetParameter(1) - taumax) < 0.001){ flag = 2; }
            if(abs(areafit->GetParameter(2) - sigmamin) < 0.001 || abs(areafit->GetParameter(2) - sigmamax) < 0.001){ flag = 3; }
            
            if(flag < 1){   areadist->Fill(simarea - areafit->GetParameter(0)); }
            
            if(i!=nwindows/2){ startbin = thisbin - 20; } 
            else{ startbin = thisbin + 90; }
            
            fakepeak->Delete();
            thiswindow->Delete();
            
        } // end of loop over windows
        //tempwindows->Print("tempwindows.pdf");
        
    } // end of loop over reps
    
    //tempwindows->Print("tempwindows.pdf]");
    
    double distRMS = areadist->GetRMS();
    
    return distRMS;
}

// -----------------------------------------------------------------------------------------------------------------------------------------------------------

// calculate signal to noise significance threshold, below which we discard all points, used to justify throwing out T1T4 and other low data in the model fitting 
double calcThresh(TProfile* cf){
    
    double threshold(0);
    int nsizes(80); // number of loops of size
    double startpeak(0.0);
    double incpeak(0.05e-6);
    int reps(10);
    
    TRandom3* tr = new TRandom3();
    tr->SetSeed(0);
    
    // set up parameters for peaks
    double simtau(0.0), simsigma(4.0);
    double taumin(-10), taumax(10), sigmamin(2.0), sigmamax(7.0);
    int nwindows(16), windowsize(40); // nwindows(18)
    double thisarea(0.0);
    
    // fit function
    TF1* areafit = new TF1("areafit","([0]*exp(-pow((x-[1])/[2],2)/2.0))/([2]*sqrt(2*pi))",-300,300);
    areafit->SetParName(0,"A");          areafit->SetParameter(0,0.0);
    areafit->SetParName(1,"#tau_{0}");   areafit->SetParameter(1, simtau);        areafit->SetParLimits(1, taumin, taumax);
    areafit->SetParName(2,"#sigma");     areafit->SetParameter(2, simsigma);      areafit->SetParLimits(2, sigmamin, sigmamax);
    
    TGraph* pcthresh = new TGraph; pcthresh->SetMarkerStyle(20);
    
    // big loop varying sizes of peak we insert to noise
    for(int i=0; i<nsizes; i++){
        int ngood(0); // count number of fits that find the correct peak
        int nbad(0);
        startpeak = startpeak+incpeak;
        thisarea = startpeak;
        
        TH1D* taudist = new TH1D("taudist","tau dist",30,-15,15);
        taudist->SetTitle(Form("area = %f e-6", thisarea/1.0e-6));
        
        // loop to fill tau dist many times, fit many peaks
        for(int i=0; i<reps; i++){
            int startbin(21), endbin(40); // startbin(1)
            
            // loop over off peak region "windows"
            for(int iw=0; iw<nwindows; iw++){
                
                double taushift = simtau - ((tr->Rndm()*10.0) - 5.0);
                TH1D* fakepeak = simpeak(thisarea, taushift, simsigma, windowsize);
                
                int flag(0);
                int thisbin = startbin;
                TH1D* thiswindow = new TH1D("thiswindow",Form("window %d",iw+1), windowsize, -(windowsize/2), (windowsize/2));
                
                // fill the window with off-peak noise plus fake peak
                for(int ix=1; ix<=windowsize; ix++){
                    thiswindow->SetBinContent(ix, cf->GetBinContent(thisbin) + fakepeak->GetBinContent(ix));
                    thiswindow->SetBinError(ix, cf->GetBinError(thisbin));
                    thisbin++;
                }
                
                areafit->SetParameter(0,0.0);       areafit->SetParameter(1, simtau);       areafit->SetParameter(2, simsigma);
                thiswindow->Fit("areafit","Q"); //cout << areafit->GetParameter(0) << endl;
                
                // quality cuts as if this was a real peak in the analysis
                if(abs(areafit->GetParameter(1) - taumin) < 0.001 || abs(areafit->GetParameter(1) - taumax) < 0.001){ flag = 2; nbad++;}
                else if(abs(areafit->GetParameter(2) - sigmamin) < 0.001 || abs(areafit->GetParameter(2) - sigmamax) < 0.001){ flag = 3; nbad++;}
                
                if(flag < 1){
                    // check if we found the "right" peak
                    taudist->Fill(areafit->GetParameter(1) - taushift);
                    if(abs(areafit->GetParameter(1) - taushift) > 3.0){ nbad++; }
                    else{ ngood++; }
                }
                
                if(i!=nwindows/2){ startbin = thisbin - 20; } 
                else{ startbin = thisbin + 90; }
                
                fakepeak->Delete();
                thiswindow->Delete();
                
            } // end of loop over windows
        } // end of loop over reps
        
        pcthresh->AddPoint(thisarea,double(ngood)/double(nbad+ngood));
        taudist->Delete();
        
    }// end of loop over sizes
    
    double thresh(0);
    for(int a=0; a<=pcthresh->GetN(); a++){
        if((pcthresh->GetPointY(a) - 0.5) >= 0){ thresh = pcthresh->GetPointX(a); break;}
    }
    
    cout << "threshold for +/- 3: " << thresh << endl;
    
    return thresh;
}

// -----------------------------------------------------------------------------------------------------------------------------------------------------------

// create simulated peaks for function to calculate excess errorbars
TH1D* simpeak(double area, double tau, double sigma, int windowsize){
    // get random tau shift for this fake peak
    TRandom3* tr = new TRandom3();          tr->SetSeed(0);
    double taushift = tau - ((tr->Rndm()*10.0) - 5.0);
    
    TH1D* thepeak = new TH1D("thepeak",Form("peak area %f",area), windowsize, -(windowsize/2), (windowsize/2));
    
    for(int i=1; i<=windowsize; i++){
        double x = i-1-(windowsize/2);
        double val = (area*exp(-pow((x-(taushift))/sigma,2)/2.0))/(sigma*sqrt(2*TMath::Pi()));
        thepeak->SetBinContent(i,val); // is there an extra assoc errorbar with the peak?
    }
    
    return thepeak;
}
