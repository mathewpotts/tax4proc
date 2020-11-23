const Int_t nphibins        = 20;
const Double_t phi_lo       = 0;
const Double_t phi_up       = 360;

const Int_t nthebins        = 40;
const Double_t the_lo       = 0;
const Double_t the_up       = 70;

const Int_t nrpbins         = 50;
const Double_t rp_lo        = 0;
const Double_t rp_up        = 50;

const Int_t npsibins         = 20;
const Double_t psi_lo        = 40;
const Double_t psi_up        = 130;

const Int_t nphpgtbins      = 50;
const Double_t phpgt_lo     = 150;
const Double_t phpgt_up     = 1400;

const Int_t nngmpebins      = 4;
const Double_t ngmpe_lo     = 0.5;
const Double_t ngmpe_up     = 4.5;

const Int_t ntrackbins      = 10;
const Double_t track_lo     = 8;
const Double_t track_up     = 28;

bool good_weather_cut(Int_t weatcode){
  Int_t n = weatcode / 1000000 % 10;
  Int_t e = weatcode / 100000 % 10;
  Int_t s = weatcode / 10000 % 10;
  Int_t w = weatcode / 1000 % 10;
  Int_t o = weatcode / 100 % 10;
  Int_t t = weatcode / 10 % 10;
  Int_t h = weatcode % 10;

  //printf("%d%d%d%d%d%d%d\n",n,e,s,w,o,t,h);
  if (o+t+n+e == 0){
          return 1;
    }
  return 0;
}

// Define Spectrum Constants
const Double_t log10en_min = 17.0;    // Minimum Thrown Energy
const Double_t log10en_Ank = 18.65;   // Ankle break point
const Double_t log10en_GZK = 19.75;   // GZK break point
const Double_t p1 = -3.3;             // Knee Spectral index
const Double_t p2 = -2.7;             // Ankle Spectral index
const Double_t p3 = -5.1;             // GZK Spectral index
const Double_t log10en_Ank_n = log10en_Ank - log10en_min;
const Double_t log10en_GZK_n = log10en_GZK - log10en_min;
const Double_t C = TMath::Power(10,(p1-p2)*log10en_Ank_n);
const Double_t K = C * TMath::Power(10.0,(p2-p3)*log10en_GZK_n);


Double_t spectrum_function(Double_t log10en){
  // This is the differenctial Spectrum Function
  Double_t log10en_n = log10en - log10en_min;

  if (log10en_n <= log10en_Ank_n)
    {    
      return TMath::Power(10.0,p1*log10en_n);
    }
  if (log10en_n > log10en_Ank_n && log10en_n <= log10en_GZK_n)
    {
      return C*TMath::Power(10.0,p2*log10en_n);
    }
  if (log10en_n > log10en_GZK_n)
    {
      return K * TMath::Power(10.0,p3*log10en_n);
    }
  return 0;
}

Double_t weight(Double_t log10en){
  return TMath::Power(10.0, 2*(log10en-log10en_min))*spectrum_function(log10en);
}

Double_t missing_E_corr(Double_t eng_EeV){
  Double_t eng_eV = 1e18 * eng_EeV;

  TF1 *JHK = new TF1("JHK","x / (-16.62 + 3.58 * log10(x) - 0.2786 * TMath::Power(log10(x),2) + 0.009775 * TMath::Power(log10(x),3) - 0.0001299 * TMath::Power(log10(x),4))",18.0,21.0);

  Double_t E = JHK->Eval(eng_eV);

  return E;  
}

void phi()
{
///// Data /////
  TFile *f = new TFile("/media/glados/Data_Storage_2/mdtax4/pass5/data/processing/root_all_files/20190625-20200215.mdtax4.ps5.ps2g.drmdpcgf.root","r");
  TTree *taTree = (TTree*)f->Get("taTree");

  // Make a Phi histogram
  TH1D *hPhidata = new TH1D("hPhidata","18.0 < log_{10}(E/eV);#phi (degrees);N_{EVENTS} / BIN",nphibins,phi_lo,phi_up);
  TH1D *hPhidata18 = new TH1D("hPhidata18","18.0 < log_{10}(E/eV) <= 18.5;#phi (degrees);N_{EVENTS} / BIN",nphibins,phi_lo,phi_up);
  TH1D *hPhidata185 = new TH1D("hPhidata185","18.5 < log_{10}(E/eV) <= 19.0;#phi (degrees);N_{EVENTS} / BIN",nphibins,phi_lo,phi_up);
  TH1D *hPhidata19 = new TH1D("hPhidata19","19.0 < log_{10}(E/eV) <= 19.5;#phi (degrees);N_{EVENTS} / BIN",nphibins,phi_lo,phi_up);
  TH1D *hPhidata195 = new TH1D("hPhidata195","19.5 < log_{10}(E/eV) <= 20.0;#phi (degrees);N_{EVENTS} / BIN",nphibins,phi_lo,phi_up);
  TH1D *hPhidata20 = new TH1D("hPhidata20","20.0 < log_{10}(E/eV) <= 20.5;#phi (degrees);N_{EVENTS} / BIN",nphibins,phi_lo,phi_up);

  taTree->Draw("57.296*etrack->phi>>hPhidata","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.0","goff");
  taTree->Draw("57.296*etrack->phi>>hPhidata18","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=18.5","goff");
  taTree->Draw("57.296*etrack->phi>>hPhidata185","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.0","goff");
  taTree->Draw("57.296*etrack->phi>>hPhidata19","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.5","goff");
  taTree->Draw("57.296*etrack->phi>>hPhidata195","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.0","goff");
  taTree->Draw("57.296*etrack->phi>>hPhidata20","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>20.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.5","goff");

///// MC /////
  TFile *f2 = new TFile("/home/glados/Data_Storage_2/mdtax4/mc/processing/root_all_files/pr.1e17.5-1e21.ps5.all.drmdpcgf.root");
  TTree *taTree2 = (TTree*)f2->Get("taTree");

  // Make a Phi histogram
  TH1D *hPhimc = new TH1D("hPhimc","18.0 < log_{10}(E/eV);#phi (degrees);N_{EVENTS} / BIN",nphibins,phi_lo,phi_up);
  TH1D *hPhimc18 = new TH1D("hPhimc18","18.0 < log_{10}(E/eV) <= 18.5;#phi (degrees);N_{EVENTS} / BIN",nphibins,phi_lo,phi_up);
  TH1D *hPhimc185 = new TH1D("hPhimc185","18.5 < log_{10}(E/eV) <= 19.0;#phi (degrees);N_{EVENTS} / BIN",nphibins,phi_lo,phi_up);
  TH1D *hPhimc19 = new TH1D("hPhimc19","19.0 < log_{10}(E/eV) <= 19.5;#phi (degrees);N_{EVENTS} / BIN",nphibins,phi_lo,phi_up);
  TH1D *hPhimc195= new TH1D("hPhimc195","19.5 < log_{10}(E/eV) <= 20.0;#phi (degrees);N_{EVENTS} / BIN",nphibins,phi_lo,phi_up);
  TH1D *hPhimc20 = new TH1D("hPhimc20","20.0 < log_{10}(E/eV) <= 20.5;#phi (degrees);N_{EVENTS} / BIN",nphibins,phi_lo,phi_up);

  taTree2->Draw("57.296*etrack->phi>>hPhimc","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.0)");
  const Double_t normal = (Double_t)hPhidata->Integral()/(Double_t)hPhimc->Integral();
  hPhimc->Scale(normal);
  gPad->Update();
  TPaveStats *sPhi = (TPaveStats*)hPhimc->FindObject("stats");
  taTree2->Draw("57.296*etrack->phi>>hPhimc18","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=18.5)");
  const Double_t normal18 = (Double_t)hPhidata18->Integral()/(Double_t)hPhimc18->Integral();
  hPhimc18->Scale(normal18);
  gPad->Update();
  TPaveStats *sPhi18 = (TPaveStats*)hPhimc18->FindObject("stats");
  taTree2->Draw("57.296*etrack->phi>>hPhimc185","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.0)");
  const Double_t normal185 = (Double_t)hPhidata185->Integral()/(Double_t)hPhimc185->Integral();
  hPhimc185->Scale(normal185);
  gPad->Update();
  TPaveStats *sPhi185 = (TPaveStats*)hPhimc185->FindObject("stats");
  taTree2->Draw("57.296*etrack->phi>>hPhimc19","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.5)");
  const Double_t normal19 = (Double_t)hPhidata19->Integral()/(Double_t)hPhimc19->Integral();
  hPhimc19->Scale(normal19);
  gPad->Update();
  TPaveStats *sPhi19 = (TPaveStats*)hPhimc19->FindObject("stats");
  taTree2->Draw("57.296*etrack->phi>>hPhimc195","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20)");
  const Double_t normal195 = (Double_t)hPhidata195->Integral()/(Double_t)hPhimc195->Integral();
  hPhimc195->Scale(normal195);
  gPad->Update();
  TPaveStats *sPhi195 = (TPaveStats*)hPhimc195->FindObject("stats");
  taTree2->Draw("57.296*etrack->phi>>hPhimc20","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>20.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.5)");
  const Double_t normal20 = (Double_t)hPhidata20->Integral()/(Double_t)hPhimc20->Integral();
  hPhimc20->Scale(normal20);
  gPad->Update();
  TPaveStats *sPhi20 = (TPaveStats*)hPhimc20->FindObject("stats");

//// Compare ////
  gStyle->SetOptFit(1);
  
  TLegend *legend = new TLegend(0.1,0.7,0.3,0.9);
  legend->SetHeader("Phi Data-MC Comparison","C");
  TLegendEntry *mcentry = legend->AddEntry("hPhimc","MC","l");
  mcentry->SetLineColor(kRed);
  TLegendEntry *dataentry = legend->AddEntry("hPhidata","Data","p");
  dataentry->SetMarkerStyle(2);
  dataentry->SetMarkerColor(kBlue+2);
  
  hPhimc->SetLineColor(kRed);
  TCanvas *c2Phi = new TCanvas("c2Phi","c2Phi",800,600);
  auto rpphi = new TRatioPlot(hPhidata,hPhimc,"divsym");
  rpphi->SetH1DrawOpt("E0");
  rpphi->SetH2DrawOpt("HIST");
  c2Phi->cd();
  rpphi->Draw();
  rpphi->GetUpperRefYaxis()->SetRangeUser(0,10);
  rpphi->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpphi->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpphi->GetLowerRefYaxis()->SetTitle("#phi_{data}/#phi_{mcrecon}");
  rpphi->GetUpperPad()->cd();
  sPhi->Draw();
  legend->Draw();
  rpphi->GetLowerRefGraph()->Fit("pol1","W");
  
  
  hPhimc18->SetLineColor(kRed);
  TCanvas *c2Phi18 = new TCanvas("c2Phi18","c2Phi18",800,600);
  auto rpphi18 = new TRatioPlot(hPhidata18,hPhimc18,"divsym");
  rpphi18->SetH1DrawOpt("E0");
  rpphi18->SetH2DrawOpt("HIST");
  c2Phi18->cd();
  rpphi18->Draw();
  rpphi18->GetUpperRefYaxis()->SetRangeUser(0,10);
  rpphi18->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpphi18->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpphi18->GetLowerRefYaxis()->SetTitle("#phi_{data}/#phi_{mcrecon}");
  rpphi18->GetUpperPad()->cd();
  sPhi18->Draw();
  legend->Draw();
  rpphi18->GetLowerRefGraph()->Fit("pol1","W");

  hPhimc185->SetLineColor(kRed);
  TCanvas *c2Phi185 = new TCanvas("c2Phi185","c2Phi185",800,600);
  auto rpphi185 = new TRatioPlot(hPhidata185,hPhimc185,"divsym");
  rpphi185->SetH1DrawOpt("E0");
  rpphi185->SetH2DrawOpt("HIST");
  c2Phi185->cd();
  rpphi185->Draw();
  rpphi185->GetUpperRefYaxis()->SetRangeUser(0,10);
  rpphi185->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpphi185->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpphi185->GetLowerRefYaxis()->SetTitle("#phi_{data}/#phi_{mcrecon}");
  rpphi185->GetUpperPad()->cd();
  sPhi185->Draw();
  legend->Draw();
  rpphi185->GetLowerRefGraph()->Fit("pol1","W");

  hPhimc19->SetLineColor(kRed);
  TCanvas *c2Phi19 = new TCanvas("c2Phi19","c2Phi19",800,600);
  auto rpphi19 = new TRatioPlot(hPhidata19,hPhimc19,"divsym");
  rpphi19->SetH1DrawOpt("E0");
  rpphi19->SetH2DrawOpt("HIST");
  c2Phi19->cd();
  rpphi19->Draw();
  rpphi19->GetUpperRefYaxis()->SetRangeUser(0,10);
  rpphi19->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpphi19->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpphi19->GetLowerRefYaxis()->SetTitle("#phi_{data}/#phi_{mcrecon}");
  rpphi19->GetUpperPad()->cd();
  sPhi19->Draw();
  rpphi19->GetLowerRefGraph()->Fit("pol1","W");

  hPhimc195->SetLineColor(kRed);
  TCanvas *c2Phi195 = new TCanvas("c2Phi195","c2Phi195",800,600);
  auto rpphi195 = new TRatioPlot(hPhidata195,hPhimc195,"divsym");
  rpphi195->SetH1DrawOpt("E0");
  rpphi195->SetH2DrawOpt("HIST");
  c2Phi195->cd();
  rpphi195->Draw();
  rpphi195->GetUpperRefYaxis()->SetRangeUser(0,10);
  rpphi195->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpphi195->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpphi195->GetLowerRefYaxis()->SetTitle("#phi_{data}/#phi_{mcrecon}");
  rpphi195->GetUpperPad()->cd();
  sPhi195->Draw();
  legend->Draw();
  rpphi195->GetLowerRefGraph()->Fit("pol1","W");

  hPhimc20->SetLineColor(kRed);
  TCanvas *c2Phi20 = new TCanvas("c2Phi20","c2Phi20",800,600);
  auto rpphi20 = new TRatioPlot(hPhidata20,hPhimc20,"divsym");
  rpphi20->SetH1DrawOpt("E0");
  rpphi20->SetH2DrawOpt("HIST");
  c2Phi20->cd();
  rpphi20->Draw();
  rpphi20->GetUpperRefYaxis()->SetRangeUser(0,10);
  rpphi20->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpphi20->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpphi20->GetLowerRefYaxis()->SetTitle("#phi_{data}/#phi_{mcrecon}");
  rpphi20->GetUpperPad()->cd();
  sPhi20->Draw();
  legend->Draw();
  rpphi20->GetLowerRefGraph()->Fit("pol1","W");
}

void theta()
{
///// Data /////
  TFile *f = new TFile("/media/glados/Data_Storage_2/mdtax4/pass5/data/processing/root_all_files/20190625-20200215.mdtax4.ps5.ps2g.drmdpcgf.root","r");
  TTree *taTree = (TTree*)f->Get("taTree");

   // Make a theta histogram
  TH1D *hThetadata = new TH1D("hThetadata","18.0 < log_{10}(E/eV);#theta (degrees);N_{EVENTS} / BIN",nthebins,the_lo,the_up);
  TH1D *hThetadata18 = new TH1D("hThetadata18","18.0 < log_{10}(E/eV) <= 18.5;#theta (degrees);N_{EVENTS} / BIN",nthebins,the_lo,the_up);
  TH1D *hThetadata185 = new TH1D("hThetadata185","18.5 < log_{10}(E/eV) <= 19.0;#theta (degrees);N_{EVENTS} / BIN",nthebins,the_lo,the_up);
  TH1D *hThetadata19 = new TH1D("hThetadata19","19.0 < log_{10}(E/eV) <= 19.5;#theta (degrees);N_{EVENTS} / BIN",nthebins,the_lo,the_up);
  TH1D *hThetadata195 = new TH1D("hThetadata195","19.5 < log_{10}(E/eV) <= 20.0;#theta (degrees);N_{EVENTS} / BIN",nthebins,the_lo,the_up);
  TH1D *hThetadata20 = new TH1D("hThetadata20","20.0 < log_{10}(E/eV) <= 20.5;#theta (degrees);N_{EVENTS} / BIN",nthebins,the_lo,the_up);
  
  taTree->Draw("57.296*hctim->mthe[int(etrack->udata[0])]>>hThetadata","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18","goff");
  taTree->Draw("57.296*hctim->mthe[int(etrack->udata[0])]>>hThetadata18","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&& log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.0 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=18.5","goff");
  taTree->Draw("57.296*hctim->mthe[int(etrack->udata[0])]>>hThetadata185","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&& log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.0","goff");  
  taTree->Draw("57.296*hctim->mthe[int(etrack->udata[0])]>>hThetadata19","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&& log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.0 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.5","goff");
  taTree->Draw("57.296*hctim->mthe[int(etrack->udata[0])]>>hThetadata195","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&& log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.5 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.0","goff");
  taTree->Draw("57.296*hctim->mthe[int(etrack->udata[0])]>>hThetadata20","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&& log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>20.0 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.5","goff");

///// MC /////
  TFile *f2 = new TFile("/home/glados/Data_Storage_2/mdtax4/mc/processing/root_all_files/pr.1e17.5-1e21.ps5.all.drmdpcgf.root");
  TTree *taTree2 = (TTree*)f2->Get("taTree");

  // Make a theta histogram
  TH1D *hThemc = new TH1D("hThemc","18.0 < log_{10}(E/eV);#theta (degrees);N_{EVENTS} / BIN",nthebins,the_lo,the_up);
  TH1D *hThemc18 = new TH1D("hThemc18","18.0 < log_{10}(E/eV) <= 18.5;#theta (degrees);N_{EVENTS} / BIN",nthebins,the_lo,the_up);
  TH1D *hThemc185 = new TH1D("hThemc185","18.5 < log_{10}(E/eV) <= 19.0;#theta (degrees);N_{EVENTS} / BIN",nthebins,the_lo,the_up);
  TH1D *hThemc19 = new TH1D("hThemc19","19.0 < log_{10}(E/eV) <= 19.5;#theta (degrees);N_{EVENTS} / BIN",nthebins,the_lo,the_up);
  TH1D *hThemc195 = new TH1D("hThemc195","19.5 < log_{10}(E/eV) <= 20.0;#theta (degrees);N_{EVENTS} / BIN",nthebins,the_lo,the_up);
  TH1D *hThemc20 = new TH1D("hThemc20","20.0 < log_{10}(E/eV) <= 20.5;#theta (degrees);N_{EVENTS} / BIN",nthebins,the_lo,the_up);
  
  taTree2->Draw("57.296*hctim->mthe[int(etrack->udata[0])]>>hThemc","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18)");
  const Double_t normal = (Double_t)hThetadata->Integral()/(Double_t)hThemc->Integral();
  hThemc->Scale(normal);
  gPad->Update();
  TPaveStats *sTheta = (TPaveStats*)hThemc->FindObject("stats");
  taTree2->Draw("57.296*hctim->mthe[int(etrack->udata[0])]>>hThemc18","weight(log10(mc04->energy))*(etrack->qualct==1&& log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.0 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=18.5)");
  const Double_t normal18 = (Double_t)hThetadata18->Integral()/(Double_t)hThemc18->Integral();
  hThemc18->Scale(normal18);
  gPad->Update();
  TPaveStats *sTheta18 = (TPaveStats*)hThemc18->FindObject("stats");
  taTree2->Draw("57.296*hctim->mthe[int(etrack->udata[0])]>>hThemc185","weight(log10(mc04->energy))*(etrack->qualct==1&& log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.0)");
  const Double_t normal185 = (Double_t)hThetadata185->Integral()/(Double_t)hThemc185->Integral();
  hThemc185->Scale(normal185);
  gPad->Update();
  TPaveStats *sTheta185 = (TPaveStats*)hThemc185->FindObject("stats");
  taTree2->Draw("57.296*hctim->mthe[int(etrack->udata[0])]>>hThemc19","weight(log10(mc04->energy))*(etrack->qualct==1&& log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.0 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.5)");
  const Double_t normal19 = (Double_t)hThetadata19->Integral()/(Double_t)hThemc19->Integral();
  hThemc19->Scale(normal19);
  gPad->Update();
  TPaveStats *sTheta19 = (TPaveStats*)hThemc19->FindObject("stats");
  taTree2->Draw("57.296*hctim->mthe[int(etrack->udata[0])]>>hThemc195","weight(log10(mc04->energy))*(etrack->qualct==1&& log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.5 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.0)");
  const Double_t normal195 = (Double_t)hThetadata195->Integral()/(Double_t)hThemc195->Integral();
  hThemc195->Scale(normal195);
  gPad->Update();
  TPaveStats *sTheta195 = (TPaveStats*)hThemc195->FindObject("stats");
  taTree2->Draw("57.296*hctim->mthe[int(etrack->udata[0])]>>hThemc20","weight(log10(mc04->energy))*(etrack->qualct==1&& log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>20.0 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.5)");
  const Double_t normal20 = (Double_t)hThetadata20->Integral()/(Double_t)hThemc20->Integral();
  hThemc20->Scale(normal20);
  gPad->Update();
  TPaveStats *sTheta20 = (TPaveStats*)hThemc20->FindObject("stats");


//// Compare ////
  gStyle->SetOptFit(1);
  
  TLegend *legend = new TLegend(0.1,0.7,0.3,0.9);
  legend->SetHeader("Theta Data-MC Comparison","C");
  TLegendEntry *mcentry = legend->AddEntry("hThemc","MC","l");
  mcentry->SetLineColor(kRed);
  TLegendEntry *dataentry = legend->AddEntry("hThetadata","Data","p");
  dataentry->SetMarkerStyle(2);
  dataentry->SetMarkerColor(kBlue+2);
  
  hThemc->SetLineColor(kRed);
  TCanvas *c2The = new TCanvas("c2The","c2The",800,600);
  auto rpthe = new TRatioPlot(hThetadata,hThemc,"divsym");
  rpthe->SetH1DrawOpt("E0");
  rpthe->SetH2DrawOpt("HIST");
  c2The->cd();
  rpthe->Draw();
  rpthe->GetUpperRefYaxis()->SetRangeUser(0,10);
  rpthe->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpthe->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpthe->GetLowerRefYaxis()->SetTitle("#theta_{data}/#theta_{mcrecon}");
  rpthe->GetUpperPad()->cd();
  sTheta->Draw();
  legend->Draw();
  rpthe->GetLowerRefGraph()->Fit("pol1","W","",9,65);

  hThemc18->SetLineColor(kRed);
  TCanvas *c2The18 = new TCanvas("c2The18","c2The18",800,600);
  auto rpthe18 = new TRatioPlot(hThetadata18,hThemc18,"divsym");
  rpthe18->SetH1DrawOpt("E0");
  rpthe18->SetH2DrawOpt("HIST");
  c2The18->cd();
  rpthe18->Draw();
  rpthe18->GetUpperRefYaxis()->SetRangeUser(0,10);
  rpthe18->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpthe18->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpthe18->GetLowerRefYaxis()->SetTitle("#theta_{data}/#theta_{mcrecon}");
  rpthe18->GetUpperPad()->cd();
  sTheta18->Draw();
  legend->Draw();
  rpthe18->GetLowerRefGraph()->Fit("pol1","W","",10,58);

  hThemc185->SetLineColor(kRed);
  TCanvas *c2The185 = new TCanvas("c2The185","c2The185",800,600);
  auto rpthe185 = new TRatioPlot(hThetadata185,hThemc185,"divsym");
  rpthe185->SetH1DrawOpt("E0");
  rpthe185->SetH2DrawOpt("HIST");
  c2The185->cd();
  rpthe185->Draw();
  rpthe185->GetUpperRefYaxis()->SetRangeUser(0,10);
  rpthe185->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpthe185->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpthe185->GetLowerRefYaxis()->SetTitle("#theta_{data}/#theta_{mcrecon}");
  rpthe185->GetUpperPad()->cd();
  sTheta185->Draw();
  legend->Draw();
  rpthe185->GetLowerRefGraph()->Fit("pol1","W","",10,60);

  hThemc19->SetLineColor(kRed);
  TCanvas *c2The19 = new TCanvas("c2The19","c2The19",800,600);
  auto rpthe19 = new TRatioPlot(hThetadata19,hThemc19,"divsym");
  rpthe19->SetH1DrawOpt("E0");
  rpthe19->SetH2DrawOpt("HIST");
  c2The19->cd();
  rpthe19->Draw();
  rpthe19->GetUpperRefYaxis()->SetRangeUser(0,10);
  rpthe19->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpthe19->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpthe19->GetLowerRefYaxis()->SetTitle("#theta_{data}/#theta_{mcrecon}");
  rpthe19->GetUpperPad()->cd();
  sTheta19->Draw();
  legend->Draw();
  rpthe19->GetLowerRefGraph()->Fit("pol1","W","",10,70);

  hThemc195->SetLineColor(kRed);
  TCanvas *c2The195 = new TCanvas("c2The195","c2The195",800,600);
  auto rpthe195 = new TRatioPlot(hThetadata195,hThemc195,"divsym");
  rpthe195->SetH1DrawOpt("E0");
  rpthe195->SetH2DrawOpt("HIST");
  c2The195->cd();
  rpthe195->Draw();
  rpthe195->GetUpperRefYaxis()->SetRangeUser(0,10);
  rpthe195->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpthe195->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpthe195->GetLowerRefYaxis()->SetTitle("#theta_{data}/#theta_{mcrecon}");
  rpthe195->GetUpperPad()->cd();
  sTheta195->Draw();
  legend->Draw();
  rpthe195->GetLowerRefGraph()->Fit("pol1","W","",10,70);

  hThemc20->SetLineColor(kRed);
  TCanvas *c2The20 = new TCanvas("c2The20","c2The20",800,600);
  auto rpthe20 = new TRatioPlot(hThetadata20,hThemc20,"divsym");
  rpthe20->SetH1DrawOpt("E0");
  rpthe20->SetH2DrawOpt("HIST");
  c2The20->cd();
  rpthe20->Draw();
  rpthe20->GetUpperRefYaxis()->SetRangeUser(0,10);
  rpthe20->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpthe20->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpthe20->GetLowerRefYaxis()->SetTitle("#theta_{data}/#theta_{mcrecon}");
  rpthe20->GetUpperPad()->cd();
  sTheta20->Draw();
  legend->Draw();
  rpthe20->GetLowerRefGraph()->Fit("pol1","W","",10,70);
}

void rp()
{
///// Data /////
  TFile *f = new TFile("/media/glados/Data_Storage_2/mdtax4/pass5/data/processing/root_all_files/20190625-20200215.mdtax4.ps5.ps2g.drmdpcgf.root","r");
  TTree *taTree = (TTree*)f->Get("taTree");

  // Make a Rp histogram
  TH1D *hRpdata = new TH1D("hRpdata","18.0 < log_{10}(E/eV);Rp (km);N_{EVENTS} / BIN",nrpbins,rp_lo,rp_up);
  TH1D *hRpdata18 = new TH1D("hRpdata18","18.0 < log_{10}(E/eV) <= 18.5;Rp (km);N_{EVENTS} / BIN",nrpbins,rp_lo,rp_up);
  TH1D *hRpdata185 = new TH1D("hRpdata185","18.5 < log_{10}(E/eV) <= 19.0;Rp (km);N_{EVENTS} / BIN",nrpbins,rp_lo,rp_up);
  TH1D *hRpdata19 = new TH1D("hRpdata19","19.0 < log_{10}(E/eV) <= 19.5;Rp (km);N_{EVENTS} / BIN",nrpbins,rp_lo,rp_up);
  TH1D *hRpdata195 = new TH1D("hRpdata195","19.5 < log_{10}(E/eV) <= 20.0;Rp (km);N_{EVENTS} / BIN",nrpbins,rp_lo,rp_up);
  TH1D *hRpdata20 = new TH1D("hRpdata20","20.0 < log_{10}(E/eV) <= 20.5;Rp (km);N_{EVENTS} / BIN",nrpbins,rp_lo,rp_up);

  taTree->Draw("1e-3*hctim->mrp[int(etrack->udata[0])]>>hRpdata","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18","goff");
  taTree->Draw("1e-3*hctim->mrp[int(etrack->udata[0])]>>hRpdata18","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=18.5","goff");
  taTree->Draw("1e-3*hctim->mrp[int(etrack->udata[0])]>>hRpdata185","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.0","goff");
  taTree->Draw("1e-3*hctim->mrp[int(etrack->udata[0])]>>hRpdata19","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.5","goff");
  taTree->Draw("1e-3*hctim->mrp[int(etrack->udata[0])]>>hRpdata195","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.0","goff");
  taTree->Draw("1e-3*hctim->mrp[int(etrack->udata[0])]>>hRpdata20","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>20.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.5","goff");


///// MC /////
  TFile *f2 = new TFile("/home/glados/Data_Storage_2/mdtax4/mc/processing/root_all_files/pr.1e17.5-1e21.ps5.all.drmdpcgf.root");
  TTree *taTree2 = (TTree*)f2->Get("taTree");

  // Make a Rp histogram
  TH1D *hRpmc = new TH1D("hRpmc","18.0 < log_{10}(E/eV);Rp (km);N_{EVENTS} / BIN",nrpbins,rp_lo,rp_up);
  TH1D *hRpmc18 = new TH1D("hRpmc18","18.0 < log_{10}(E/eV) <= 18.5;Rp (km);N_{EVENTS} / BIN",nrpbins,rp_lo,rp_up);
  TH1D *hRpmc185 = new TH1D("hRpmc185","18.5 < log_{10}(E/eV) <= 19.0;Rp (km);N_{EVENTS} / BIN",nrpbins,rp_lo,rp_up);
  TH1D *hRpmc19 = new TH1D("hRpmc19","19.0 < log_{10}(E/eV) <= 19.5;Rp (km);N_{EVENTS} / BIN",nrpbins,rp_lo,rp_up);
  TH1D *hRpmc195= new TH1D("hRpmc195","19.5 < log_{10}(E/eV) <= 20.0;Rp (km);N_{EVENTS} / BIN",nrpbins,rp_lo,rp_up);
  TH1D *hRpmc20 = new TH1D("hRpmc20","20.0 < log_{10}(E/eV) <= 20.5;Rp (km);N_{EVENTS} / BIN",nrpbins,rp_lo,rp_up);

  taTree2->Draw("1e-3*hctim->mrp[int(etrack->udata[0])]>>hRpmc","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18)");
  const Double_t normal = (Double_t)hRpdata->Integral()/(Double_t)hRpmc->Integral();
  hRpmc->Scale(normal);
  gPad->Update();
  TPaveStats *sRp = (TPaveStats*)hRpmc->FindObject("stats");
  taTree2->Draw("1e-3*hctim->mrp[int(etrack->udata[0])]>>hRpmc18","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=18.5)");
  const Double_t normal18 = (Double_t)hRpdata18->Integral()/(Double_t)hRpmc18->Integral();
  hRpmc18->Scale(normal18);
  gPad->Update();
  TPaveStats *sRp18 = (TPaveStats*)hRpmc18->FindObject("stats");
  taTree2->Draw("1e-3*hctim->mrp[int(etrack->udata[0])]>>hRpmc185","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.0)");
  const Double_t normal185 = (Double_t)hRpdata185->Integral()/(Double_t)hRpmc185->Integral();
  hRpmc185->Scale(normal185);
  gPad->Update();
  TPaveStats *sRp185 = (TPaveStats*)hRpmc185->FindObject("stats");
  taTree2->Draw("1e-3*hctim->mrp[int(etrack->udata[0])]>>hRpmc19","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.5)");
  const Double_t normal19 = (Double_t)hRpdata19->Integral()/(Double_t)hRpmc19->Integral();
  hRpmc19->Scale(normal19);
  gPad->Update();
  TPaveStats *sRp19 = (TPaveStats*)hRpmc19->FindObject("stats");
  taTree2->Draw("1e-3*hctim->mrp[int(etrack->udata[0])]>>hRpmc195","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20)");
  const Double_t normal195 = (Double_t)hRpdata195->Integral()/(Double_t)hRpmc195->Integral();
  hRpmc195->Scale(normal195);
  gPad->Update();
  TPaveStats *sRp195 = (TPaveStats*)hRpmc195->FindObject("stats");
  taTree2->Draw("1e-3*hctim->mrp[int(etrack->udata[0])]>>hRpmc20","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>20.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.5)");
  const Double_t normal20 = (Double_t)hRpdata20->Integral()/(Double_t)hRpmc20->Integral();
  hRpmc20->Scale(normal20);
  gPad->Update();
  TPaveStats *sRp20 = (TPaveStats*)hRpmc20->FindObject("stats");

//// Compare ////
  gStyle->SetOptFit(1);
   
  TLegend *legend = new TLegend(0.1,0.7,0.3,0.9);
  legend->SetHeader("Rp Data-MC Comparison","C");
  TLegendEntry *mcentry = legend->AddEntry("hRpmc","MC","l");
  mcentry->SetLineColor(kRed);
  TLegendEntry *dataentry = legend->AddEntry("hRpdata","Data","p");
  dataentry->SetMarkerStyle(2);
  dataentry->SetMarkerColor(kBlue+2);

  
  hRpmc->SetLineColor(kRed);
  TCanvas *c2Rp = new TCanvas("c2Rp","c2Rp",800,600);
  auto rprp = new TRatioPlot(hRpdata,hRpmc,"divsym");
  rprp->SetH1DrawOpt("E0");
  rprp->SetH2DrawOpt("HIST");
  c2Rp->cd();
  rprp->Draw();
  rprp->GetUpperRefYaxis()->SetRangeUser(0,15);
  rprp->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rprp->GetLowerRefYaxis()->SetLabelSize(0.025);
  rprp->GetLowerRefYaxis()->SetTitle("Rp_{data}/Rp_{mcrecon}");
  rprp->GetUpperPad()->cd();
  sRp->Draw();
  legend->Draw();
  rprp->GetLowerRefGraph()->Fit("pol1","W");
  
  hRpmc18->SetLineColor(kRed);
  TCanvas *c2Rp18 = new TCanvas("c2Rp18","c2Rp18",800,600);
  auto rprp18 = new TRatioPlot(hRpdata18,hRpmc18,"divsym");
  rprp18->SetH1DrawOpt("E0");
  rprp18->SetH2DrawOpt("HIST");
  c2Rp18->cd();
  rprp18->Draw();
  rprp18->GetUpperRefYaxis()->SetRangeUser(0,15);
  rprp18->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rprp18->GetLowerRefYaxis()->SetLabelSize(0.025);
  rprp18->GetLowerRefYaxis()->SetTitle("Rp_{data}/Rp_{mcrecon}");
  rprp18->GetUpperPad()->cd();
  sRp18->Draw();
  legend->Draw();
  rprp18->GetLowerRefGraph()->Fit("pol1","W");

  hRpmc185->SetLineColor(kRed);
  TCanvas *c2Rp185 = new TCanvas("c2Rp185","c2Rp185",800,600);
  auto rprp185 = new TRatioPlot(hRpdata185,hRpmc185,"divsym");
  rprp185->SetH1DrawOpt("E0");
  rprp185->SetH2DrawOpt("HIST");
  c2Rp185->cd();
  rprp185->Draw();
  rprp185->GetUpperRefYaxis()->SetRangeUser(0,15);
  rprp185->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rprp185->GetLowerRefYaxis()->SetLabelSize(0.025);
  rprp185->GetLowerRefYaxis()->SetTitle("Rp_{data}/Rp_{mcrecon}");
  rprp185->GetUpperPad()->cd();
  sRp185->Draw();
  legend->Draw();
  rprp185->GetLowerRefGraph()->Fit("pol1","W");

  hRpmc19->SetLineColor(kRed);
  TCanvas *c2Rp19 = new TCanvas("c2Rp19","c2Rp19",800,600);
  auto rprp19 = new TRatioPlot(hRpdata19,hRpmc19,"divsym");
  rprp19->SetH1DrawOpt("E0");
  rprp19->SetH2DrawOpt("HIST");
  c2Rp19->cd();
  rprp19->Draw();
  rprp19->GetUpperRefYaxis()->SetRangeUser(0,15);
  rprp19->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rprp19->GetLowerRefYaxis()->SetLabelSize(0.025);
  rprp19->GetLowerRefYaxis()->SetTitle("Rp_{data}/Rp_{mcrecon}");
  rprp19->GetUpperPad()->cd();
  sRp19->Draw();
  legend->Draw();
  rprp19->GetLowerRefGraph()->Fit("pol1","W");

  hRpmc195->SetLineColor(kRed);
  TCanvas *c2Rp195 = new TCanvas("c2Rp195","c2Rp195",800,600);
  auto rprp195 = new TRatioPlot(hRpdata195,hRpmc195,"divsym");
  rprp195->SetH1DrawOpt("E0");
  rprp195->SetH2DrawOpt("HIST");
  c2Rp195->cd();
  rprp195->Draw();
  rprp195->GetUpperRefYaxis()->SetRangeUser(0,15);
  rprp195->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rprp195->GetLowerRefYaxis()->SetLabelSize(0.025);
  rprp195->GetLowerRefYaxis()->SetTitle("Rp_{data}/Rp_{mcrecon}");
  rprp195->GetUpperPad()->cd();
  sRp195->Draw();
  legend->Draw();
  rprp195->GetLowerRefGraph()->Fit("pol1","W");

  hRpmc20->SetLineColor(kRed);
  TCanvas *c2Rp20 = new TCanvas("c2Rp20","c2Rp20",800,600);
  auto rprp20 = new TRatioPlot(hRpdata20,hRpmc20,"divsym");
  rprp20->SetH1DrawOpt("E0");
  rprp20->SetH2DrawOpt("HIST");
  c2Rp20->cd();
  rprp20->Draw();
  rprp20->GetUpperRefYaxis()->SetRangeUser(0,15);
  rprp20->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rprp20->GetLowerRefYaxis()->SetLabelSize(0.025);
  rprp20->GetLowerRefYaxis()->SetTitle("Rp_{data}/Rp_{mcrecon}");
  rprp20->GetUpperPad()->cd();
  sRp20->Draw();
  legend->Draw();
  rprp20->GetLowerRefGraph()->Fit("pol1","W");
}

void psi()
{
///// Data /////
  TFile *f = new TFile("/media/glados/Data_Storage_2/mdtax4/pass5/data/processing/root_all_files/20190625-20200215.mdtax4.ps5.ps2g.drmdpcgf.root","r");
  TTree *taTree = (TTree*)f->Get("taTree");

  // Make a Psi histogram
  TH1D *hPsidata = new TH1D("hPsidata","18.0 < log_{10}(E/eV);#psi (degrees);N_{EVENTS} / BIN",npsibins,psi_lo,psi_up);
  TH1D *hPsidata18 = new TH1D("hPsidata18","18.0 < log_{10}(E/eV) <= 18.5;#psi (degrees);N_{EVENTS} / BIN",npsibins,psi_lo,psi_up);
  TH1D *hPsidata185 = new TH1D("hPsidata185","18.5 < log_{10}(E/eV) <= 19.0;#psi (degrees);N_{EVENTS} / BIN",npsibins,psi_lo,psi_up);
  TH1D *hPsidata19 = new TH1D("hPsidata19","19.0 < log_{10}(E/eV) <= 19.5;#psi (degrees);N_{EVENTS} / BIN",npsibins,psi_lo,psi_up);
  TH1D *hPsidata195 = new TH1D("hPsidata195","19.5 < log_{10}(E/eV) <= 20.0;#psi (degrees);N_{EVENTS} / BIN",npsibins,psi_lo,psi_up);
  TH1D *hPsidata20 = new TH1D("hPsidata20","20.0 < log_{10}(E/eV) <= 20.5;#psi (degrees);N_{EVENTS} / BIN",npsibins,psi_lo,psi_up);

  taTree->Draw("57.296*hctim->mpsi[int(etrack->udata[0])]>>hPsidata","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18","goff");
  
  taTree->Draw("57.296*hctim->mpsi[int(etrack->udata[0])]>>hPsidata18","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=18.5","goff");
  taTree->Draw("57.296*hctim->mpsi[int(etrack->udata[0])]>>hPsidata185","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.0","goff");
  taTree->Draw("57.296*hctim->mpsi[int(etrack->udata[0])]>>hPsidata19","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.5","goff");
  taTree->Draw("57.296*hctim->mpsi[int(etrack->udata[0])]>>hPsidata195","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.0","goff");
  taTree->Draw("57.296*hctim->mpsi[int(etrack->udata[0])]>>hPsidata20","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>20.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.5","goff");

///// MC /////
  TFile *f2 = new TFile("/home/glados/Data_Storage_2/mdtax4/mc/processing/root_all_files/pr.1e17.5-1e21.ps5.all.drmdpcgf.root");
  TTree *taTree2 = (TTree*)f2->Get("taTree");

  // Make a Psi histogram
  TH1D *hPsimc = new TH1D("hPsimc","18.0 < log_{10}(E/eV);#psi (degrees);N_{EVENTS} / BIN",npsibins,psi_lo,psi_up);
  TH1D *hPsimc18 = new TH1D("hPsimc18","18.0 < log_{10}(E/eV) <= 18.5;#psi (degrees);N_{EVENTS} / BIN",npsibins,psi_lo,psi_up);
  TH1D *hPsimc185 = new TH1D("hPsimc185","18.5 < log_{10}(E/eV) <= 19.0;#psi (degrees);N_{EVENTS} / BIN",npsibins,psi_lo,psi_up);
  TH1D *hPsimc19 = new TH1D("hPsimc19","19.0 < log_{10}(E/eV) <= 19.5;#psi (degrees);N_{EVENTS} / BIN",npsibins,psi_lo,psi_up);
  TH1D *hPsimc195= new TH1D("hPsimc195","19.5 < log_{10}(E/eV) <= 20.0;#psi (degrees);N_{EVENTS} / BIN",npsibins,psi_lo,psi_up);
  TH1D *hPsimc20 = new TH1D("hPsimc20","20.0 < log_{10}(E/eV) <= 20.5;#psi (degrees);N_{EVENTS} / BIN",npsibins,psi_lo,psi_up);

  taTree2->Draw("57.296*hctim->mpsi[int(etrack->udata[0])]>>hPsimc","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18)");
  const Double_t normal = (Double_t)hPsidata->Integral()/(Double_t)hPsimc->Integral();
  hPsimc->Scale(normal);
  gPad->Update();
  TPaveStats *sPsi = (TPaveStats*)hPsimc->FindObject("stats");
  taTree2->Draw("57.296*hctim->mpsi[int(etrack->udata[0])]>>hPsimc18","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=18.5)");
  const Double_t normal18 = (Double_t)hPsidata18->Integral()/(Double_t)hPsimc18->Integral();
  hPsimc18->Scale(normal18);
  gPad->Update();
  TPaveStats *sPsi18 = (TPaveStats*)hPsimc18->FindObject("stats");
  taTree2->Draw("57.296*hctim->mpsi[int(etrack->udata[0])]>>hPsimc185","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.0)");
  const Double_t normal185 = (Double_t)hPsidata185->Integral()/(Double_t)hPsimc185->Integral();
  hPsimc185->Scale(normal185);
  gPad->Update();
  TPaveStats *sPsi185 = (TPaveStats*)hPsimc185->FindObject("stats");
  taTree2->Draw("57.296*hctim->mpsi[int(etrack->udata[0])]>>hPsimc19","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.5)");
  const Double_t normal19 = (Double_t)hPsidata19->Integral()/(Double_t)hPsimc19->Integral();
  hPsimc19->Scale(normal19);
  gPad->Update();
  TPaveStats *sPsi19 = (TPaveStats*)hPsimc19->FindObject("stats");
  taTree2->Draw("57.296*hctim->mpsi[int(etrack->udata[0])]>>hPsimc195","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20)");
  const Double_t normal195 = (Double_t)hPsidata195->Integral()/(Double_t)hPsimc195->Integral();
  hPsimc195->Scale(normal195);
  gPad->Update();
  TPaveStats *sPsi195 = (TPaveStats*)hPsimc195->FindObject("stats");
  taTree2->Draw("57.296*hctim->mpsi[int(etrack->udata[0])]>>hPsimc20","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>20.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.5)");
  const Double_t normal20 = (Double_t)hPsidata20->Integral()/(Double_t)hPsimc20->Integral();
  hPsimc20->Scale(normal20);
  gPad->Update();
  TPaveStats *sPsi20 = (TPaveStats*)hPsimc20->FindObject("stats");

//// Compare ////
  gStyle->SetOptFit(1);
   
  TLegend *legend = new TLegend(0.1,0.7,0.3,0.9);
  legend->SetHeader("Psi Data-MC Comparison","C");
  TLegendEntry *mcentry = legend->AddEntry("hPsimc","MC","l");
  mcentry->SetLineColor(kRed);
  TLegendEntry *dataentry = legend->AddEntry("hPsidata","Data","p");
  dataentry->SetMarkerStyle(2);
  dataentry->SetMarkerColor(kBlue+2);
  
  hPsimc->SetLineColor(kRed);
  TCanvas *c2Psi = new TCanvas("c2Psi","c2Psi",800,600);
  auto rppsi = new TRatioPlot(hPsidata,hPsimc,"divsym");
  rppsi->SetH1DrawOpt("E0");
  rppsi->SetH2DrawOpt("HIST");
  c2Psi->cd();
  rppsi->Draw();
  rppsi->GetUpperRefYaxis()->SetRangeUser(0,15);
  rppsi->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rppsi->GetLowerRefYaxis()->SetLabelSize(0.025);
  rppsi->GetLowerRefYaxis()->SetTitle("#psi_{data}/#psi_{mcrecon}");
  rppsi->GetUpperPad()->cd();
  sPsi->Draw();
  legend->Draw();
  rppsi->GetLowerRefGraph()->Fit("pol1","W");
  
  hPsimc18->SetLineColor(kRed);
  TCanvas *c2Psi18 = new TCanvas("c2Psi18","c2Psi18",800,600);
  auto rppsi18 = new TRatioPlot(hPsidata18,hPsimc18,"divsym");
  rppsi18->SetH1DrawOpt("E0");
  rppsi18->SetH2DrawOpt("HIST");
  c2Psi18->cd();
  rppsi18->Draw();
  rppsi18->GetUpperRefYaxis()->SetRangeUser(0,15);
  rppsi18->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rppsi18->GetLowerRefYaxis()->SetLabelSize(0.025);
  rppsi18->GetLowerRefYaxis()->SetTitle("#psi_{data}/#psi_{mcrecon}");
  rppsi18->GetUpperPad()->cd();
  sPsi18->Draw();
  legend->Draw();
  rppsi18->GetLowerRefGraph()->Fit("pol1","W");

  hPsimc185->SetLineColor(kRed);
  TCanvas *c2Psi185 = new TCanvas("c2Psi185","c2Psi185",800,600);
  auto rppsi185 = new TRatioPlot(hPsidata185,hPsimc185,"divsym");
  rppsi185->SetH1DrawOpt("E0");
  rppsi185->SetH2DrawOpt("HIST");
  c2Psi185->cd();
  rppsi185->Draw();
  rppsi185->GetUpperRefYaxis()->SetRangeUser(0,15);
  rppsi185->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rppsi185->GetLowerRefYaxis()->SetLabelSize(0.025);
  rppsi185->GetLowerRefYaxis()->SetTitle("#psi_{data}/#psi_{mcrecon}");
  rppsi185->GetUpperPad()->cd();
  sPsi185->Draw();
  legend->Draw();
  rppsi185->GetLowerRefGraph()->Fit("pol1","W");

  hPsimc19->SetLineColor(kRed);
  TCanvas *c2Psi19 = new TCanvas("c2Psi19","c2Psi19",800,600);
  auto rppsi19 = new TRatioPlot(hPsidata19,hPsimc19,"divsym");
  rppsi19->SetH1DrawOpt("E0");
  rppsi19->SetH2DrawOpt("HIST");
  c2Psi19->cd();
  rppsi19->Draw();
  rppsi19->GetUpperRefYaxis()->SetRangeUser(0,15);
  rppsi19->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rppsi19->GetLowerRefYaxis()->SetLabelSize(0.025);
  rppsi19->GetLowerRefYaxis()->SetTitle("#psi_{data}/#psi_{mcrecon}");
  rppsi19->GetUpperPad()->cd();
  sPsi19->Draw();
  legend->Draw();
  rppsi19->GetLowerRefGraph()->Fit("pol1","W");

  hPsimc195->SetLineColor(kRed);
  TCanvas *c2Psi195 = new TCanvas("c2Psi195","c2Psi195",800,600);
  auto rppsi195 = new TRatioPlot(hPsidata195,hPsimc195,"divsym");
  rppsi195->SetH1DrawOpt("E0");
  rppsi195->SetH2DrawOpt("HIST");
  c2Psi195->cd();
  rppsi195->Draw();
  rppsi195->GetUpperRefYaxis()->SetRangeUser(0,15);
  rppsi195->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rppsi195->GetLowerRefYaxis()->SetLabelSize(0.025);
  rppsi195->GetLowerRefYaxis()->SetTitle("#psi_{data}/#psi_{mcrecon}");
  rppsi195->GetUpperPad()->cd();
  sPsi195->Draw();
  legend->Draw();
  rppsi195->GetLowerRefGraph()->Fit("pol1","W");

  hPsimc20->SetLineColor(kRed);
  TCanvas *c2Psi20 = new TCanvas("c2Psi20","c2Psi20",800,600);
  auto rppsi20 = new TRatioPlot(hPsidata20,hPsimc20,"divsym");
  rppsi20->SetH1DrawOpt("E0");
  rppsi20->SetH2DrawOpt("HIST");
  c2Psi20->cd();
  rppsi20->Draw();
  rppsi20->GetUpperRefYaxis()->SetRangeUser(0,15);
  rppsi20->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rppsi20->GetLowerRefYaxis()->SetLabelSize(0.025);
  rppsi20->GetLowerRefYaxis()->SetTitle("#psi_{data}/#psi_{mcrecon}");
  rppsi20->GetUpperPad()->cd();
  sPsi20->Draw();
  legend->Draw();
  rppsi20->GetLowerRefGraph()->Fit("pol1","W");
}

void phpgt()
{
///// Data /////
  TFile *f = new TFile("/media/glados/Data_Storage_2/mdtax4/pass5/data/processing/root_all_files/20190625-20200215.mdtax4.ps5.ps2g.drmdpcgf.root","r");
  TTree *taTree = (TTree*)f->Get("taTree");

  // Make a Phpgt histogram
  TH1D *hPhpgtdata = new TH1D("hPhpgtdata","18.0 < log_{10}(E/eV);phpgt;N_{EVENTS} / BIN",nphpgtbins,phpgt_lo,phpgt_up);
  TH1D *hPhpgtdata18 = new TH1D("hPhpgtdata18","18.0 < log_{10}(E/eV) <= 18.5;phpgt;N_{EVENTS} / BIN",nphpgtbins,phpgt_lo,phpgt_up);
  TH1D *hPhpgtdata185 = new TH1D("hPhpgtdata185","18.5 < log_{10}(E/eV) <= 19.0;phpgt;N_{EVENTS} / BIN",nphpgtbins,phpgt_lo,phpgt_up);
  TH1D *hPhpgtdata19 = new TH1D("hPhpgtdata19","19.0 < log_{10}(E/eV) <= 19.5;phpgt;N_{EVENTS} / BIN",nphpgtbins,phpgt_lo,phpgt_up);
  TH1D *hPhpgtdata195 = new TH1D("hPhpgtdata195","19.5 < log_{10}(E/eV) <= 20.0;phpgt;N_{EVENTS} / BIN",nphpgtbins,phpgt_lo,phpgt_up);
  TH1D *hPhpgtdata20 = new TH1D("hPhpgtdata20","20.0 < log_{10}(E/eV) <= 20.5;phpgt;N_{EVENTS} / BIN",nphpgtbins,phpgt_lo,phpgt_up);

  taTree->Draw("stpln->ph_per_gtube[6]>>hPhpgtdata","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18","goff");
  taTree->Draw("stpln->ph_per_gtube[6]>>hPhpgtdata18","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=18.5","goff"); 
  taTree->Draw("stpln->ph_per_gtube[6]>>hPhpgtdata185","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.0","goff");
  taTree->Draw("stpln->ph_per_gtube[6]>>hPhpgtdata19","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.5","goff");
  taTree->Draw("stpln->ph_per_gtube[6]>>hPhpgtdata195","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.0","goff");
  taTree->Draw("stpln->ph_per_gtube[6]>>hPhpgtdata20","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>20.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.5","goff");

///// MC /////
  TFile *f2 = new TFile("/home/glados/Data_Storage_2/mdtax4/mc/processing/root_all_files/pr.1e17.5-1e21.ps5.all.drmdpcgf.root");
  TTree *taTree2 = (TTree*)f2->Get("taTree");

  // Make a Phpgt histogram
  TH1D *hPhpgtmc = new TH1D("hPhpgtmc","18.0 < log_{10}(E/eV);phpgt;N_{EVENTS} / BIN",nphpgtbins,phpgt_lo,phpgt_up);
  TH1D *hPhpgtmc18 = new TH1D("hPhpgtmc18","18.0 < log_{10}(E/eV) <= 18.5;phpgt;N_{EVENTS} / BIN",nphpgtbins,phpgt_lo,phpgt_up);
  TH1D *hPhpgtmc185 = new TH1D("hPhpgtmc185","18.5 < log_{10}(E/eV) <= 19.0;phpgt;N_{EVENTS} / BIN",nphpgtbins,phpgt_lo,phpgt_up);
  TH1D *hPhpgtmc19 = new TH1D("hPhpgtmc19","19.0 < log_{10}(E/eV) <= 19.5;phpgt;N_{EVENTS} / BIN",nphpgtbins,phpgt_lo,phpgt_up);
  TH1D *hPhpgtmc195= new TH1D("hPhpgtmc195","19.5 < log_{10}(E/eV) <= 20.0;phpgt;N_{EVENTS} / BIN",nphpgtbins,phpgt_lo,phpgt_up);
  TH1D *hPhpgtmc20 = new TH1D("hPhpgtmc20","20.0 < log_{10}(E/eV) <= 20.5;phpgt;N_{EVENTS} / BIN",nphpgtbins,phpgt_lo,phpgt_up);

  taTree2->Draw("stpln->ph_per_gtube[6]>>hPhpgtmc","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18)");
  const Double_t normal = (Double_t)hPhpgtdata->Integral()/(Double_t)hPhpgtmc->Integral();
  hPhpgtmc->Scale(normal);
  gPad->Update();
  TPaveStats *sPhpgt = (TPaveStats*)hPhpgtmc->FindObject("stats");
  taTree2->Draw("stpln->ph_per_gtube[6]>>hPhpgtmc18","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=18.5)");
  const Double_t normal18 = (Double_t)hPhpgtdata18->Integral()/(Double_t)hPhpgtmc18->Integral();
  hPhpgtmc18->Scale(normal18);
  gPad->Update();
  TPaveStats *sPhpgt18 = (TPaveStats*)hPhpgtmc18->FindObject("stats");
  taTree2->Draw("stpln->ph_per_gtube[6]>>hPhpgtmc185","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.0)");
  const Double_t normal185 = (Double_t)hPhpgtdata185->Integral()/(Double_t)hPhpgtmc185->Integral();
  hPhpgtmc185->Scale(normal185);
  gPad->Update();
  TPaveStats *sPhpgt185 = (TPaveStats*)hPhpgtmc185->FindObject("stats");
  taTree2->Draw("stpln->ph_per_gtube[6]>>hPhpgtmc19","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.5)");
  const Double_t normal19 = (Double_t)hPhpgtdata19->Integral()/(Double_t)hPhpgtmc19->Integral();
  hPhpgtmc19->Scale(normal19);
  gPad->Update();
  TPaveStats *sPhpgt19 = (TPaveStats*)hPhpgtmc19->FindObject("stats");
  taTree2->Draw("stpln->ph_per_gtube[6]>>hPhpgtmc195","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20)");
  const Double_t normal195 = (Double_t)hPhpgtdata195->Integral()/(Double_t)hPhpgtmc195->Integral();
  hPhpgtmc195->Scale(normal195);
  gPad->Update();
  TPaveStats *sPhpgt195 = (TPaveStats*)hPhpgtmc195->FindObject("stats");
  taTree2->Draw("stpln->ph_per_gtube[6]>>hPhpgtmc20","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>20.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.5)");
  const Double_t normal20 = (Double_t)hPhpgtdata20->Integral()/(Double_t)hPhpgtmc20->Integral();
  hPhpgtmc20->Scale(normal20);
  gPad->Update();
  TPaveStats *sPhpgt20 = (TPaveStats*)hPhpgtmc20->FindObject("stats");

//// Compare ////
  gStyle->SetOptFit(1);
   
  TLegend *legend = new TLegend(0.1,0.7,0.3,0.9);
  legend->SetHeader("PHPGT Data-MC Comparison","C");
  TLegendEntry *mcentry = legend->AddEntry("hPhpgtmc","MC","l");
  mcentry->SetLineColor(kRed);
  TLegendEntry *dataentry = legend->AddEntry("hPhpgtdata","Data","p");
  dataentry->SetMarkerStyle(2);
  dataentry->SetMarkerColor(kBlue+2);
  
  hPhpgtmc->SetLineColor(kRed);
  TCanvas *c2Phpgt = new TCanvas("c2Phpgt","c2Phpgt",800,600);
  auto rpphpgt = new TRatioPlot(hPhpgtdata,hPhpgtmc,"divsym");
  rpphpgt->SetH1DrawOpt("E0");
  rpphpgt->SetH2DrawOpt("HIST");
  c2Phpgt->cd();
  rpphpgt->Draw();
  rpphpgt->GetUpperRefYaxis()->SetRangeUser(0,20);
  rpphpgt->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpphpgt->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpphpgt->GetLowerRefYaxis()->SetTitle("phpgt_{data}/phpgt_{mcrecon}");
  rpphpgt->GetUpperPad()->cd();
  sPhpgt->Draw();
  legend->Draw();
  rpphpgt->GetLowerRefGraph()->Fit("pol1","W","",150,600);
  
  hPhpgtmc18->SetLineColor(kRed);
  TCanvas *c2Phpgt18 = new TCanvas("c2Phpgt18","c2Phpgt18",800,600);
  auto rpphpgt18 = new TRatioPlot(hPhpgtdata18,hPhpgtmc18,"divsym");
  rpphpgt18->SetH1DrawOpt("E0");
  rpphpgt18->SetH2DrawOpt("HIST");
  c2Phpgt18->cd();
  rpphpgt18->Draw();
  rpphpgt18->GetUpperRefYaxis()->SetRangeUser(0,15);
  rpphpgt18->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpphpgt18->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpphpgt18->GetLowerRefYaxis()->SetTitle("phpgt_{data}/phpgt_{mcrecon}");
  rpphpgt18->GetUpperPad()->cd();
  sPhpgt18->Draw();
  legend->Draw();
  rpphpgt18->GetLowerRefGraph()->Fit("pol1","W","",150,600);

  hPhpgtmc185->SetLineColor(kRed);
  TCanvas *c2Phpgt185 = new TCanvas("c2Phpgt185","c2Phpgt185",800,600);
  auto rpphpgt185 = new TRatioPlot(hPhpgtdata185,hPhpgtmc185,"divsym");
  rpphpgt185->SetH1DrawOpt("E0");
  rpphpgt185->SetH2DrawOpt("HIST");
  c2Phpgt185->cd();
  rpphpgt185->Draw();
  rpphpgt185->GetUpperRefYaxis()->SetRangeUser(0,15);
  rpphpgt185->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpphpgt185->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpphpgt185->GetLowerRefYaxis()->SetTitle("phpgt_{data}/phpgt_{mcrecon}");
  rpphpgt185->GetUpperPad()->cd();
  sPhpgt185->Draw();
  legend->Draw();
  rpphpgt185->GetLowerRefGraph()->Fit("pol1","W","",200,600);

  hPhpgtmc19->SetLineColor(kRed);
  TCanvas *c2Phpgt19 = new TCanvas("c2Phpgt19","c2Phpgt19",800,600);
  auto rpphpgt19 = new TRatioPlot(hPhpgtdata19,hPhpgtmc19,"divsym");
  rpphpgt19->SetH1DrawOpt("E0");
  rpphpgt19->SetH2DrawOpt("HIST");
  c2Phpgt19->cd();
  rpphpgt19->Draw();
  rpphpgt19->GetUpperRefYaxis()->SetRangeUser(0,15);
  rpphpgt19->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpphpgt19->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpphpgt19->GetLowerRefYaxis()->SetTitle("phpgt_{data}/phpgt_{mcrecon}");
  rpphpgt19->GetUpperPad()->cd();
  sPhpgt19->Draw();
  legend->Draw();
  rpphpgt19->GetLowerRefGraph()->Fit("pol1","W");

  hPhpgtmc195->SetLineColor(kRed);
  TCanvas *c2Phpgt195 = new TCanvas("c2Phpgt195","c2Phpgt195",800,600);
  auto rpphpgt195 = new TRatioPlot(hPhpgtdata195,hPhpgtmc195,"divsym");
  rpphpgt195->SetH1DrawOpt("E0");
  rpphpgt195->SetH2DrawOpt("HIST");
  c2Phpgt195->cd();
  rpphpgt195->Draw();
  rpphpgt195->GetUpperRefYaxis()->SetRangeUser(0,15);
  rpphpgt195->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpphpgt195->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpphpgt195->GetLowerRefYaxis()->SetTitle("phpgt_{data}/phpgt_{mcrecon}");
  rpphpgt195->GetUpperPad()->cd();
  sPhpgt195->Draw();
  legend->Draw();
  rpphpgt195->GetLowerRefGraph()->Fit("pol1","W");

  hPhpgtmc20->SetLineColor(kRed);
  TCanvas *c2Phpgt20 = new TCanvas("c2Phpgt20","c2Phpgt20",800,600);
  auto rpphpgt20 = new TRatioPlot(hPhpgtdata20,hPhpgtmc20,"divsym");
  rpphpgt20->SetH1DrawOpt("E0");
  rpphpgt20->SetH2DrawOpt("HIST");
  c2Phpgt20->cd();
  rpphpgt20->Draw();
  rpphpgt20->GetUpperRefYaxis()->SetRangeUser(0,15);
  rpphpgt20->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpphpgt20->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpphpgt20->GetLowerRefYaxis()->SetTitle("phpgt_{data}/phpgt_{mcrecon}");
  rpphpgt20->GetUpperPad()->cd();
  sPhpgt20->Draw();
  legend->Draw();
  rpphpgt20->GetLowerRefGraph()->Fit("pol1","W");

}

void track()
{
///// Data /////
  TFile *f = new TFile("/media/glados/Data_Storage_2/mdtax4/pass5/data/processing/root_all_files/20190625-20200215.mdtax4.ps5.ps2g.drmdpcgf.root","r");
  TTree *taTree = (TTree*)f->Get("taTree");

  // Make a Track histogram
  TH1D *hTrackdata = new TH1D("hTrackdata","18.0 < log_{10}(E/eV);track (degrees);N_{EVENTS} / BIN",ntrackbins,track_lo,track_up);
  TH1D *hTrackdata18 = new TH1D("hTrackdata18","18.0 < log_{10}(E/eV) <= 18.5;track (degrees);N_{EVENTS} / BIN",ntrackbins,track_lo,track_up);
  TH1D *hTrackdata185 = new TH1D("hTrackdata185","18.5 < log_{10}(E/eV) <= 19.0;track (degrees);N_{EVENTS} / BIN",ntrackbins,track_lo,track_up);
  TH1D *hTrackdata19 = new TH1D("hTrackdata19","19.0 < log_{10}(E/eV) <= 19.5;track (degrees);N_{EVENTS} / BIN",ntrackbins,track_lo,track_up);
  TH1D *hTrackdata195 = new TH1D("hTrackdata195","19.5 < log_{10}(E/eV) <= 20.0;track (degrees);N_{EVENTS} / BIN",ntrackbins,track_lo,track_up);
  TH1D *hTrackdata20 = new TH1D("hTrackdata20","20.0 < log_{10}(E/eV) <= 20.5;track (degrees);N_{EVENTS} / BIN",ntrackbins,track_lo,track_up);

  taTree->Draw("stpln->tracklength[6]>>hTrackdata","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18","goff");
  taTree->Draw("stpln->tracklength[6]>>hTrackdata18","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=18.5","goff");
  taTree->Draw("stpln->tracklength[6]>>hTrackdata185","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.0","goff");
  taTree->Draw("stpln->tracklength[6]>>hTrackdata19","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.5","goff");
  taTree->Draw("stpln->tracklength[6]>>hTrackdata195","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.0","goff");
  taTree->Draw("stpln->tracklength[6]>>hTrackdata20","good_weather_cut(tlweat->code)==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>20.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.5","goff");
  
///// MC /////
  TFile *f2 = new TFile("/home/glados/Data_Storage_2/mdtax4/mc/processing/root_all_files/pr.1e17.5-1e21.ps5.all.drmdpcgf.root");
  TTree *taTree2 = (TTree*)f2->Get("taTree");

  // Make a Track histogram
  TH1D *hTrackmc = new TH1D("hTrackmc","18.0 < log_{10}(E/eV);track (degrees);N_{EVENTS} / BIN",ntrackbins,track_lo,track_up);
  TH1D *hTrackmc18 = new TH1D("hTrackmc18","18.0 < log_{10}(E/eV) <= 18.5;track (degrees);N_{EVENTS} / BIN",ntrackbins,track_lo,track_up);
  TH1D *hTrackmc185 = new TH1D("hTrackmc185","18.5 < log_{10}(E/eV) <= 19.0;track (degrees);N_{EVENTS} / BIN",ntrackbins,track_lo,track_up);
  TH1D *hTrackmc19 = new TH1D("hTrackmc19","19.0 < log_{10}(E/eV) <= 19.5;track (degrees);N_{EVENTS} / BIN",ntrackbins,track_lo,track_up);
  TH1D *hTrackmc195= new TH1D("hTrackmc195","19.5 < log_{10}(E/eV) <= 20.0;track (degrees);N_{EVENTS} / BIN",ntrackbins,track_lo,track_up);
  TH1D *hTrackmc20 = new TH1D("hTrackmc20","20.0 < log_{10}(E/eV) <= 20.5;track (degrees);N_{EVENTS} / BIN",ntrackbins,track_lo,track_up);

  taTree2->Draw("stpln->tracklength[6]>>hTrackmc","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18)");
  const Double_t normal = (Double_t)hTrackdata->Integral()/(Double_t)hTrackmc->Integral();
  hTrackmc->Scale(normal);
  gPad->Update();
  TPaveStats *sTrack = (TPaveStats*)hTrackmc->FindObject("stats");
  taTree2->Draw("stpln->tracklength[6]>>hTrackmc18","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=18.5)");
  const Double_t normal18 = (Double_t)hTrackdata18->Integral()/(Double_t)hTrackmc18->Integral();
  hTrackmc18->Scale(normal18);
  gPad->Update();
  TPaveStats *sTrack18 = (TPaveStats*)hTrackmc18->FindObject("stats");
  taTree2->Draw("stpln->tracklength[6]>>hTrackmc185","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.0)");
  const Double_t normal185 = (Double_t)hTrackdata185->Integral()/(Double_t)hTrackmc185->Integral();
  hTrackmc185->Scale(normal185);
  gPad->Update();
  TPaveStats *sTrack185 = (TPaveStats*)hTrackmc185->FindObject("stats");
  taTree2->Draw("stpln->tracklength[6]>>hTrackmc19","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.5)");
  const Double_t normal19 = (Double_t)hTrackdata19->Integral()/(Double_t)hTrackmc19->Integral();
  hTrackmc19->Scale(normal19);
  gPad->Update();
  TPaveStats *sTrack19 = (TPaveStats*)hTrackmc19->FindObject("stats");
  taTree2->Draw("stpln->tracklength[6]>>hTrackmc195","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20)");
  const Double_t normal195 = (Double_t)hTrackdata195->Integral()/(Double_t)hTrackmc195->Integral();
  hTrackmc195->Scale(normal195);
  gPad->Update();
  TPaveStats *sTrack195 = (TPaveStats*)hTrackmc195->FindObject("stats");
  taTree2->Draw("stpln->tracklength[6]>>hTrackmc20","weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>20.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.5)");
  const Double_t normal20 = (Double_t)hTrackdata20->Integral()/(Double_t)hTrackmc20->Integral();
  hTrackmc20->Scale(normal20);
  gPad->Update();
  TPaveStats *sTrack20 = (TPaveStats*)hTrackmc20->FindObject("stats");

//// Compare ////
  gStyle->SetOptFit(1);
   
  TLegend *legend = new TLegend(0.1,0.7,0.3,0.9);
  legend->SetHeader("Track length Data-MC Comparison","C");
  TLegendEntry *mcentry = legend->AddEntry("hTrackmc","MC","l");
  mcentry->SetLineColor(kRed);
  TLegendEntry *dataentry = legend->AddEntry("hTrackdata","Data","p");
  dataentry->SetMarkerStyle(2);
  dataentry->SetMarkerColor(kBlue+2);
  
  hTrackmc->SetLineColor(kRed);
  TCanvas *c2Track = new TCanvas("c2Track","c2Track",800,600);
  auto rptrack = new TRatioPlot(hTrackdata,hTrackmc,"divsym");
  rptrack->SetH1DrawOpt("E0");
  rptrack->SetH2DrawOpt("HIST");
  c2Track->cd();
  rptrack->Draw();
  rptrack->GetUpperRefYaxis()->SetRangeUser(0,35);
  rptrack->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rptrack->GetLowerRefYaxis()->SetLabelSize(0.025);
  rptrack->GetLowerRefYaxis()->SetTitle("track_{data}/track_{mcrecon}");
  rptrack->GetUpperPad()->cd();
  sTrack->Draw();
  legend->Draw();
  rptrack->GetLowerRefGraph()->Fit("pol1","W");
  
  hTrackmc18->SetLineColor(kRed);
  TCanvas *c2Track18 = new TCanvas("c2Track18","c2Track18",800,600);
  auto rptrack18 = new TRatioPlot(hTrackdata18,hTrackmc18,"divsym");
  rptrack18->SetH1DrawOpt("E0");
  rptrack18->SetH2DrawOpt("HIST");
  c2Track18->cd();
  rptrack18->Draw();
  rptrack18->GetUpperRefYaxis()->SetRangeUser(0,35);
  rptrack18->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rptrack18->GetLowerRefYaxis()->SetLabelSize(0.025);
  rptrack18->GetLowerRefYaxis()->SetTitle("track_{data}/track_{mcrecon}");
  rptrack18->GetUpperPad()->cd();
  sTrack18->Draw();
  legend->Draw();
  rptrack18->GetLowerRefGraph()->Fit("pol1","W");

  hTrackmc185->SetLineColor(kRed);
  TCanvas *c2Track185 = new TCanvas("c2Track185","c2Track185",800,600);
  auto rptrack185 = new TRatioPlot(hTrackdata185,hTrackmc185,"divsym");
  rptrack185->SetH1DrawOpt("E0");
  rptrack185->SetH2DrawOpt("HIST");
  c2Track185->cd();
  rptrack185->Draw();
  rptrack185->GetUpperRefYaxis()->SetRangeUser(0,10);
  rptrack185->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rptrack185->GetLowerRefYaxis()->SetLabelSize(0.025);
  rptrack185->GetLowerRefYaxis()->SetTitle("track_{data}/track_{mcrecon}");
  rptrack185->GetUpperPad()->cd();
  sTrack185->Draw();
  legend->Draw();
  rptrack185->GetLowerRefGraph()->Fit("pol1","W");

  hTrackmc19->SetLineColor(kRed);
  TCanvas *c2Track19 = new TCanvas("c2Track19","c2Track19",800,600);
  auto rptrack19 = new TRatioPlot(hTrackdata19,hTrackmc19,"divsym");
  rptrack19->SetH1DrawOpt("E0");
  rptrack19->SetH2DrawOpt("HIST");
  c2Track19->cd();
  rptrack19->Draw();
  rptrack19->GetUpperRefYaxis()->SetRangeUser(0,10);
  rptrack19->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rptrack19->GetLowerRefYaxis()->SetLabelSize(0.025);
  rptrack19->GetLowerRefYaxis()->SetTitle("track_{data}/track_{mcrecon}");
  rptrack19->GetUpperPad()->cd();
  sTrack19->Draw();
  legend->Draw();
  rptrack19->GetLowerRefGraph()->Fit("pol1","W");

  hTrackmc195->SetLineColor(kRed);
  TCanvas *c2Track195 = new TCanvas("c2Track195","c2Track195",800,600);
  auto rptrack195 = new TRatioPlot(hTrackdata195,hTrackmc195,"divsym");
  rptrack195->SetH1DrawOpt("E0");
  rptrack195->SetH2DrawOpt("HIST");
  c2Track195->cd();
  rptrack195->Draw();
  rptrack195->GetUpperRefYaxis()->SetRangeUser(0,10);
  rptrack195->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rptrack195->GetLowerRefYaxis()->SetLabelSize(0.025);
  rptrack195->GetLowerRefYaxis()->SetTitle("track_{data}/track_{mcrecon}");
  rptrack195->GetUpperPad()->cd();
  sTrack195->Draw();
  legend->Draw();
  rptrack195->GetLowerRefGraph()->Fit("pol1","W");

  hTrackmc20->SetLineColor(kRed);
  TCanvas *c2Track20 = new TCanvas("c2Track20","c2Track20",800,600);
  auto rptrack20 = new TRatioPlot(hTrackdata20,hTrackmc20,"divsym");
  rptrack20->SetH1DrawOpt("E0");
  rptrack20->SetH2DrawOpt("HIST");
  c2Track20->cd();
  rptrack20->Draw();
  rptrack20->GetUpperRefYaxis()->SetRangeUser(0,10);
  rptrack20->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rptrack20->GetLowerRefYaxis()->SetLabelSize(0.025);
  rptrack20->GetLowerRefYaxis()->SetTitle("track_{data}/track_{mcrecon}");
  rptrack20->GetUpperPad()->cd();
  sTrack20->Draw();
  legend->Draw();
  rptrack20->GetLowerRefGraph()->Fit("pol1","W");
}

void ngmpe()
{
  TFile *f = new TFile("/media/glados/Data_Storage_2/mdtax4/pass5/data/processing/root_all_files/20190625-20200215.mdtax4.ps5.ps2g.drmdpcgf.root","r");
  TTree *taTree = (TTree*)f->Get("taTree");

  stpln_class *stplndata = new stpln_class;
  tlweat_class *tlweatdata = new tlweat_class;
  etrack_class *etrackdata = new etrack_class;
  prfc_class *prfcdata = new prfc_class;
  
  taTree->SetBranchAddress("stpln",&stplndata);
  taTree->SetBranchAddress("tlweat",&tlweatdata);
  taTree->SetBranchAddress("prfc",&prfcdata);
  taTree->SetBranchAddress("etrack",&etrackdata);

  const Int_t nentriesdata = taTree->GetEntries();

  TH1D *hngmpedata = new TH1D("hngmpedata","18.0 < log_{10}(E/eV);track (degrees);N_{EVENTS} / BIN",nngmpebins,ngmpe_lo,ngmpe_up);
  TH1D *hngmpedata18 = new TH1D("hngmpedata18","18.0 < log_{10}(E/eV) <= 18.5;track (degrees);N_{EVENTS} / BIN",nngmpebins,ngmpe_lo,ngmpe_up);
  TH1D *hngmpedata185 = new TH1D("hngmpedata185","18.5 < log_{10}(E/eV) <= 19.0;track (degrees);N_{EVENTS} / BIN",nngmpebins,ngmpe_lo,ngmpe_up);
  TH1D *hngmpedata19 = new TH1D("hngmpedata19","19.0 < log_{10}(E/eV) <= 19.5;track (degrees);N_{EVENTS} / BIN",nngmpebins,ngmpe_lo,ngmpe_up);
  TH1D *hngmpedata195 = new TH1D("hngmpedata195","19.5 < log_{10}(E/eV) <= 20.0;track (degrees);N_{EVENTS} / BIN",nngmpebins,ngmpe_lo,ngmpe_up);
  TH1D *hngmpedata20 = new TH1D("hngmpedata20","20.0 < log_{10}(E/eV) <= 20.5;track (degrees);N_{EVENTS} / BIN",nngmpebins,ngmpe_lo,ngmpe_up);
  
  for (Int_t i=0;i < nentriesdata; i++)
    {
      Int_t n_gmir = 0;
      Int_t n_gmir18 = 0;
      Int_t n_gmir185 = 0;
      Int_t n_gmir19 = 0;
      Int_t n_gmir195 = 0;
      Int_t n_gmir20 = 0;
      taTree->GetEntry(i);

      if (!taTree->Draw("etrack->energy","good_weather_cut(tlweat->code)==1 && etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18","goff",1,i)) continue;
      
      // Loop over all possible  mirror instances
      for (Int_t j=0; j<4;j++)
	{
	  if (stplndata->mir_ngtube[j] > 0 && stplndata->mir_ngtube[j] <= 256 && good_weather_cut(tlweatdata->code)==1 && etrackdata->qualct==1 && log10((missing_E_corr(prfcdata->eng[etrackdata->udata[0]])))>18)
	    {
	      n_gmir += 1;
	    }
	  if (stplndata->mir_ngtube[j] > 0 && stplndata->mir_ngtube[j] <= 256 && good_weather_cut(tlweatdata->code)==1 && etrackdata->qualct==1 && log10((missing_E_corr(prfcdata->eng[etrackdata->udata[0]])))<=18.5 && log10((missing_E_corr(prfcdata->eng[etrackdata->udata[0]])))>18)
	    {
	      n_gmir18 += 1;
	    }
	  if (stplndata->mir_ngtube[j] > 0 && stplndata->mir_ngtube[j] <= 256 && good_weather_cut(tlweatdata->code)==1 && etrackdata->qualct==1 && log10((missing_E_corr(prfcdata->eng[etrackdata->udata[0]])))<=19.0 && log10((missing_E_corr(prfcdata->eng[etrackdata->udata[0]])))>18.5)
	    {
	      n_gmir185 += 1;
	    }
	  if (stplndata->mir_ngtube[j] > 0 && stplndata->mir_ngtube[j] <= 256 && good_weather_cut(tlweatdata->code)==1 && etrackdata->qualct==1 && log10((missing_E_corr(prfcdata->eng[etrackdata->udata[0]])))<=19.5 && log10((missing_E_corr(prfcdata->eng[etrackdata->udata[0]])))>19)
	    {
	      n_gmir19 += 1;
	    }
	  if (stplndata->mir_ngtube[j] > 0 && stplndata->mir_ngtube[j] <= 256 && good_weather_cut(tlweatdata->code)==1 && etrackdata->qualct==1 && log10((missing_E_corr(prfcdata->eng[etrackdata->udata[0]])))<=20.0 && log10((missing_E_corr(prfcdata->eng[etrackdata->udata[0]])))>19.5)
	    {
	      n_gmir195 += 1;
	    }
	  if (stplndata->mir_ngtube[j] > 0 && stplndata->mir_ngtube[j] <= 256 && good_weather_cut(tlweatdata->code)==1 && etrackdata->qualct==1 && log10((missing_E_corr(prfcdata->eng[etrackdata->udata[0]])))<=20.5 && log10((missing_E_corr(prfcdata->eng[etrackdata->udata[0]])))>20.0)
	    {
	      n_gmir20 += 1;
	    }
	}
      hngmpedata->Fill(n_gmir);
      hngmpedata18->Fill(n_gmir18);
      hngmpedata185->Fill(n_gmir185);
      hngmpedata19->Fill(n_gmir19);
      hngmpedata195->Fill(n_gmir195);
      hngmpedata20->Fill(n_gmir20);
    }
  hngmpedata->Draw("goff");

  // MC
  TFile *f2 = new TFile("/home/glados/Data_Storage_2/mdtax4/mc/processing/root_all_files/pr.1e17.5-1e21.ps5.all.drmdpcgf.root");
  TTree *taTree2 = (TTree*)f2->Get("taTree");

  stpln_class *stplnmc   = new stpln_class;
  mc04_class *mc04       = new mc04_class;
  etrack_class *etrackmc = new etrack_class;
  prfc_class *prfcmc = new prfc_class;

  taTree2->SetBranchAddress("stpln",&stplnmc);
  taTree2->SetBranchAddress("mc04",&mc04);
  taTree2->SetBranchAddress("prfc",&prfcmc);
  taTree2->SetBranchAddress("etrack",&etrackmc);

  const Int_t nentriesmc = taTree2->GetEntries();

  TH1D *hngmpemc = new TH1D("hngmpemc","18.0 < log_{10}(E/eV);NGMPE;N_{EVENTS} / BIN",nngmpebins,ngmpe_lo,ngmpe_up);
  TH1D *hngmpemc18 = new TH1D("hngmpemc18","18.0 < log_{10}(E/eV);NGMPE;N_{EVENTS} / BIN",nngmpebins,ngmpe_lo,ngmpe_up);
  TH1D *hngmpemc185 = new TH1D("hngmpemc185","18.0 < log_{10}(E/eV);NGMPE;N_{EVENTS} / BIN",nngmpebins,ngmpe_lo,ngmpe_up);
  TH1D *hngmpemc19 = new TH1D("hngmpemc19","18.0 < log_{10}(E/eV);NGMPE;N_{EVENTS} / BIN",nngmpebins,ngmpe_lo,ngmpe_up);
  TH1D *hngmpemc195 = new TH1D("hngmpemc195","18.0 < log_{10}(E/eV);NGMPE;N_{EVENTS} / BIN",nngmpebins,ngmpe_lo,ngmpe_up);
  TH1D *hngmpemc20 = new TH1D("hngmpemc20","18.0 < log_{10}(E/eV);NGMPE;N_{EVENTS} / BIN",nngmpebins,ngmpe_lo,ngmpe_up);
  
  for (Int_t i=0;i < nentriesmc; i++)
    {
      Int_t n_gmir = 0;
      Int_t n_gmir18 = 0;
      Int_t n_gmir185 = 0;
      Int_t n_gmir19 = 0;
      Int_t n_gmir195 = 0;
      Int_t n_gmir20 = 0;
      taTree2->GetEntry(i);

      if (!taTree2->Draw("etrack->energy","etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18","goff",1,i)) continue;
      
      // Loop over all possible  mirror instances
      for (Int_t j=0; j<4;j++)
  	{
  	  if (stplnmc->mir_ngtube[j] > 0 && stplnmc->mir_ngtube[j] <= 256 && etrackmc->qualct==1 && log10((missing_E_corr(prfcmc->eng[etrackmc->udata[0]])))>18)
  	    {
  	      n_gmir += 1;
  	    }
	  if (stplnmc->mir_ngtube[j] > 0 && stplnmc->mir_ngtube[j] <= 256 && etrackmc->qualct==1 && log10((missing_E_corr(prfcmc->eng[etrackmc->udata[0]])))<=18.5 && log10((missing_E_corr(prfcmc->eng[etrackmc->udata[0]])))>18)
	    {
	      n_gmir18 += 1;
	    }
	  if (stplnmc->mir_ngtube[j] > 0 && stplnmc->mir_ngtube[j] <= 256 && etrackmc->qualct==1 && log10((missing_E_corr(prfcmc->eng[etrackmc->udata[0]])))<=19.0 && log10((missing_E_corr(prfcmc->eng[etrackmc->udata[0]])))>18.5)
	    {
	      n_gmir185 += 1;
	    }
	  if (stplnmc->mir_ngtube[j] > 0 && stplnmc->mir_ngtube[j] <= 256 &&  etrackmc->qualct==1 && log10((missing_E_corr(prfcmc->eng[etrackmc->udata[0]])))<=19.5 && log10((missing_E_corr(prfcmc->eng[etrackmc->udata[0]])))>19)
	    {
	      n_gmir19 += 1;
	    }
	  if (stplnmc->mir_ngtube[j] > 0 && stplnmc->mir_ngtube[j] <= 256 &&  etrackmc->qualct==1 && log10((missing_E_corr(prfcmc->eng[etrackmc->udata[0]])))<=20.0 && log10((missing_E_corr(prfcmc->eng[etrackmc->udata[0]])))>19.5)
	    {
	      n_gmir195 += 1;
	    }
	  if (stplndata->mir_ngtube[j] > 0 && stplndata->mir_ngtube[j] <= 256 && etrackmc->qualct==1 && log10((missing_E_corr(prfcmc->eng[etrackmc->udata[0]])))<=20.5 && log10((missing_E_corr(prfcmc->eng[etrackmc->udata[0]])))>20.0)
	    {
	      n_gmir20 += 1;
	    }
  	}
      hngmpemc->Fill(n_gmir,weight(log10(mc04->energy)));
      hngmpemc18->Fill(n_gmir18,weight(log10(mc04->energy)));
      hngmpemc185->Fill(n_gmir185,weight(log10(mc04->energy)));
      hngmpemc19->Fill(n_gmir19,weight(log10(mc04->energy)));
      hngmpemc195->Fill(n_gmir195,weight(log10(mc04->energy)));
      hngmpemc20->Fill(n_gmir20,weight(log10(mc04->energy)));
    }
  
  Double_t normal =  hngmpedata->Integral()/hngmpemc->Integral();
  Double_t normal18 =  hngmpedata18->Integral()/hngmpemc18->Integral();
  Double_t normal185 =  hngmpedata185->Integral()/hngmpemc185->Integral();
  Double_t normal19 =  hngmpedata19->Integral()/hngmpemc19->Integral();
  Double_t normal195 =  hngmpedata195->Integral()/hngmpemc195->Integral();
  Double_t normal20 =  hngmpedata20->Integral()/hngmpemc20->Integral();
  
  hngmpemc->Scale(normal);
  hngmpemc->Draw();
  gPad->Update();
  TPaveStats *sNgmpe = (TPaveStats*)hngmpemc->FindObject("stats");
  hngmpemc18->Scale(normal18);
  hngmpemc18->Draw("");
  gPad->Update();
  TPaveStats *sNgmpe18 = (TPaveStats*)hngmpemc18->FindObject("stats");
  hngmpemc185->Scale(normal185);
  hngmpemc185->Draw("");
  gPad->Update();
  TPaveStats *sNgmpe185 = (TPaveStats*)hngmpemc185->FindObject("stats");
  hngmpemc19->Scale(normal19);
  hngmpemc19->Draw("");
  gPad->Update();
  TPaveStats *sNgmpe19 = (TPaveStats*)hngmpemc19->FindObject("stats");
  hngmpemc195->Scale(normal);
  hngmpemc195->Draw("");
  gPad->Update();
  TPaveStats *sNgmpe195 = (TPaveStats*)hngmpemc195->FindObject("stats");
  hngmpemc20->Scale(normal);
  hngmpemc20->Draw("");
  gPad->Update();
  TPaveStats *sNgmpe20 = (TPaveStats*)hngmpemc20->FindObject("stats");

  //// Compare ////
  gStyle->SetOptFit(1);
   
  TLegend *legend = new TLegend(0.1,0.7,0.3,0.9);
  legend->SetHeader("NGMPE Data-MC Comparison","C");
  TLegendEntry *mcentry = legend->AddEntry("hngmpemc","MC","l");
  mcentry->SetLineColor(kRed);
  TLegendEntry *dataentry = legend->AddEntry("hngmpedata","Data","p");
  dataentry->SetMarkerStyle(2);
  dataentry->SetMarkerColor(kBlue+2);
  
  hngmpemc->SetLineColor(kRed);
  TCanvas *c2ngmpe = new TCanvas("c2ngmpe","c2ngmpe",800,600);
  auto rpngmpe = new TRatioPlot(hngmpedata,hngmpemc,"divsym");
  rpngmpe->SetH1DrawOpt("E0");
  rpngmpe->SetH2DrawOpt("HIST");
  c2ngmpe->cd();
  rpngmpe->Draw();
  rpngmpe->GetUpperRefYaxis()->SetRangeUser(0,80);
  rpngmpe->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpngmpe->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpngmpe->GetLowerRefYaxis()->SetTitle("ngmpe_{data}/ngmpe_{mcrecon}");
  rpngmpe->GetUpperPad()->cd();
  sNgmpe->Draw();
  legend->Draw();
  rpngmpe->GetLowerRefGraph()->Fit("pol1","W","",0,2.5);

  hngmpemc18->SetLineColor(kRed);
  TCanvas *c2ngmpe18 = new TCanvas("c2ngmpe18","c2ngmpe18",800,600);
  auto rpngmpe18 = new TRatioPlot(hngmpedata18,hngmpemc18,"divsym");
  rpngmpe18->SetH1DrawOpt("E0");
  rpngmpe18->SetH2DrawOpt("HIST");
  c2ngmpe18->cd();
  rpngmpe18->Draw();
  rpngmpe18->GetUpperRefYaxis()->SetRangeUser(0,80);
  rpngmpe18->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpngmpe18->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpngmpe18->GetLowerRefYaxis()->SetTitle("ngmpe_{data}/ngmpe_{mcrecon}");
  rpngmpe18->GetUpperPad()->cd();
  sNgmpe18->Draw();
  legend->Draw();
  rpngmpe18->GetLowerRefGraph()->Fit("pol1","W","",0,2.5);

  hngmpemc185->SetLineColor(kRed);
  TCanvas *c2ngmpe185 = new TCanvas("c2ngmpe185","c2ngmpe185",800,600);
  auto rpngmpe185 = new TRatioPlot(hngmpedata185,hngmpemc185,"divsym");
  rpngmpe185->SetH1DrawOpt("E0");
  rpngmpe185->SetH2DrawOpt("HIST");
  c2ngmpe185->cd();
  rpngmpe185->Draw();
  rpngmpe185->GetUpperRefYaxis()->SetRangeUser(0,80);
  rpngmpe185->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpngmpe185->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpngmpe185->GetLowerRefYaxis()->SetTitle("ngmpe_{data}/ngmpe_{mcrecon}");
  rpngmpe185->GetUpperPad()->cd();
  sNgmpe185->Draw();
  legend->Draw();
  rpngmpe185->GetLowerRefGraph()->Fit("pol1","W","",0,2.5);

  hngmpemc19->SetLineColor(kRed);
  TCanvas *c2ngmpe19 = new TCanvas("c2ngmpe19","c2ngmpe19",800,600);
  auto rpngmpe19 = new TRatioPlot(hngmpedata19,hngmpemc19,"divsym");
  rpngmpe19->SetH1DrawOpt("E0");
  rpngmpe19->SetH2DrawOpt("HIST");
  c2ngmpe19->cd();
  rpngmpe19->Draw();
  rpngmpe19->GetUpperRefYaxis()->SetRangeUser(0,80);
  rpngmpe19->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpngmpe19->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpngmpe19->GetLowerRefYaxis()->SetTitle("ngmpe_{data}/ngmpe_{mcrecon}");
  rpngmpe19->GetUpperPad()->cd();
  sNgmpe19->Draw();
  legend->Draw();
  rpngmpe19->GetLowerRefGraph()->Fit("pol1","W","",0,2.5);

  hngmpemc195->SetLineColor(kRed);
  TCanvas *c2ngmpe195 = new TCanvas("c2ngmpe195","c2ngmpe195",800,600);
  auto rpngmpe195 = new TRatioPlot(hngmpedata195,hngmpemc195,"divsym");
  rpngmpe195->SetH1DrawOpt("E0");
  rpngmpe195->SetH2DrawOpt("HIST");
  c2ngmpe195->cd();
  rpngmpe195->Draw();
  rpngmpe195->GetUpperRefYaxis()->SetRangeUser(0,80);
  rpngmpe195->GetLowerRefYaxis()->SetRangeUser(0.1,3);
  rpngmpe195->GetLowerRefYaxis()->SetLabelSize(0.025);
  rpngmpe195->GetLowerRefYaxis()->SetTitle("ngmpe_{data}/ngmpe_{mcrecon}");
  rpngmpe195->GetUpperPad()->cd();
  sNgmpe195->Draw();
  legend->Draw();
  rpngmpe195->GetLowerRefGraph()->Fit("pol1","W","",0,2.5);
}

void data_mc_compare()
{
///// Data /////


  //phi();
  
  //theta();

  //rp();

  //psi();

  //phpgt();

  //track();

  ngmpe();

}
