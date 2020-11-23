Double_t missing_E_corr(Double_t eng_EeV)
{
  Double_t eng_eV = 1e18 * eng_EeV;

  TF1 *JHK = new TF1("JHK","x / (-16.62 + 3.58 * log10(x) - 0.2786 * TMath::Power(log10(x),2) + 0.009775 * TMath::Power(log10(x),3) - 0.0001299 * TMath::Power(log10(x),4))",18.0,21.0);

  Double_t E = JHK->Eval(eng_eV);

  return E; 
}


void energy(TTree *taTree)
{
  /////// Energy Resolution Plots ///////
  TCanvas *ce18 = new TCanvas("ce18","ce18",600,600);
  TCanvas *ce185 = new TCanvas("ce185","ce185",600,600);
  TCanvas *ce19 = new TCanvas("ce19","ce19",600,600);
  TCanvas *ce195 = new TCanvas("ce195","ce195",600,600);
  TCanvas *ce20 = new TCanvas("ce20","ce20",600,600);
  
  TH1F *hE18 = new TH1F("hE18","18.0 < log_{10}(E_{RECON}/eV) <= 18.5;ln(E_{RECON}/E_{THROWN});N_{EVENTS}/bin",100,-1,1);
  TH1F *hE185 = new TH1F("hE185","18.5 < log_{10}(E_{RECON}/eV) <= 19.0;ln(E_{RECON}/E_{THROWN});N_{EVENTS}/bin",100,-1,1);
  TH1F *hE19 = new TH1F("hE19","19.0 < log_{10}(E_{RECON}/eV) <= 19.5;ln(E_{RECON}/E_{THROWN});N_{EVENTS}/bin",100,-1,1);
  TH1F *hE195 = new TH1F("hE195","19.5 < log_{10}(E_{RECON}/eV) <= 20.0;ln(E_{RECON}/E_{THROWN});N_{EVENTS}/bin",100,-1,1);
  TH1F *hE20 = new TH1F("hE20","20.0 < log_{10}(E_{RECON}/eV) <= 20.5;ln(E_{RECON}/E_{THROWN});N_{EVENTS}/bin",100,-1,1);
  
  // hE18->SetLineColor(kBlack);hE18->SetFillColorAlpha(42,0.35);
  // hE185->SetLineColor(kBlack);hE185->SetFillColorAlpha(46,0.35);
  // hE19->SetLineColor(kBlack);hE19->SetFillColorAlpha(30,0.35);
  // hE195->SetLineColor(kBlack);hE195->SetFillColorAlpha(44,0.35);
  // hE20->SetLineColor(kBlack);hE20->SetFillColorAlpha(32,0.35);

  gStyle->SetOptFit(1);

  ce18->cd();
  taTree->Draw("(log((missing_E_corr(prfc->eng[etrack->udata[0]]))/mc04->energy))>>hE18","etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.0 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=18.5");
  hE18->Fit("gaus");hE18->GetFunction("gaus")->SetLineColor(kRed);hE18->GetFunction("gaus")->SetLineWidth(3);
  ce18->Update();
  ce18->SaveAs("/media/glados/Data_Storage_2/plots/resolution/energy_resolution_hist_18.png");
  
  ce185->cd();
  taTree->Draw("(log((missing_E_corr(prfc->eng[etrack->udata[0]]))/mc04->energy))>>hE185","etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.0");
  hE185->Fit("gaus");hE185->GetFunction("gaus")->SetLineColor(kRed);hE185->GetFunction("gaus")->SetLineWidth(3);
  ce185->Update();
  ce185->SaveAs("/media/glados/Data_Storage_2/plots/resolution/energy_resolution_hist_185.png");

  ce19->cd();
  taTree->Draw("(log((missing_E_corr(prfc->eng[etrack->udata[0]]))/mc04->energy))>>hE19","etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.5");
  hE19->Fit("gaus");hE19->GetFunction("gaus")->SetLineColor(kRed);hE19->GetFunction("gaus")->SetLineWidth(3);
  ce19->Update();
  ce19->SaveAs("/media/glados/Data_Storage_2/plots/resolution/energy_resolution_hist_19.png");

  ce195->cd();
  taTree->Draw("(log((missing_E_corr(prfc->eng[etrack->udata[0]]))/mc04->energy))>>hE195","etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.5 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20");
  hE195->Fit("gaus");hE195->GetFunction("gaus")->SetLineColor(kRed);hE195->GetFunction("gaus")->SetLineWidth(3);
  ce195->Update();
  ce195->SaveAs("/media/glados/Data_Storage_2/plots/resolution/energy_resolution_hist_195.png");

  ce20->cd();
  taTree->Draw("(log((missing_E_corr(prfc->eng[etrack->udata[0]]))/mc04->energy))>>hE20","etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>20 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.5");
  hE20->Fit("gaus");hE20->GetFunction("gaus")->SetLineColor(kRed);hE20->GetFunction("gaus")->SetLineWidth(3);
  ce20->Update();
  ce20->SaveAs("/media/glados/Data_Storage_2/plots/resolution/energy_resolution_hist_20.png");


  // TPaveStats *se18 = (TPaveStats*)hE18->FindObject("stats");
  // se18->SetLineColor(42);se18->SetTextColor(42);se18->Draw();
  // TPaveStats *sE185 = (TPaveStats*)hE185->FindObject("stats");
  // sE185->SetLineColor(46);sE185->SetTextColor(46); sE185->Draw();
  // TPaveStats *sE19 = (TPaveStats*)hE19->FindObject("stats");
  // sE19->SetLineColor(30);sE19->SetTextColor(30); sE19->Draw();
  // TPaveStats *sE195 = (TPaveStats*)hE195->FindObject("stats");
  // sE195->SetLineColor(44);sE195->SetTextColor(44); sE195->Draw();
  // TPaveStats *sE20 = (TPaveStats*)hE20->FindObject("stats");
  // sE20->SetLineColor(32);sE20->SetTextColor(32); sE20->Draw();

  

  TCanvas *ceprof = new TCanvas("ceprof","ceprof",600,600);
  TProfile *profe = new TProfile("profe",";log_{10}(E_THROWN}/eV);ln(E_{RECON}/E_{THROWN})",30,18.0,20.5);
  TH2F *hebox = new TH2F("hebox",";log_{10}(E_{thrown}/eV);ln(E_{RECON}/E_{THROWN})",50,17,21,100,-1,1);
  ceprof->cd();
  taTree->Draw("(log((missing_E_corr(prfc->eng[etrack->udata[0]]))/mc04->energy)):log10(mc04->energy)>>hebox","etrack->qualct==1","colz");
  taTree->Draw("(log((missing_E_corr(prfc->eng[etrack->udata[0]]))/mc04->energy)):log10(mc04->energy)>>profe","etrack->qualct==1","prof same");

  profe->SetLineColor(kRed);
  profe->SetLineWidth(2);
  ceprof->Update();
  
  return hE18;
}

void psi(TTree *taTree)
{
  /////// Psi Resolution Plots ///////
  TCanvas *cpsi18 = new TCanvas("cpsi18","cpsi18",600,600);
  TCanvas *cpsi185 = new TCanvas("cpsi185","cpsi185",600,600);
  TCanvas *cpsi19 = new TCanvas("cpsi19","cpsi19",600,600);
  TCanvas *cpsi195 = new TCanvas("cpsi195","cpsi195",600,600);
  TCanvas *cpsi20 = new TCanvas("cpsi20","cpsi20",600,600);
  
  TH1F *hpsi18 = new TH1F("hpsi18","18.0 < log_{10}(E_{RECON}/eV) <= 18.5;#psi_{RECON} - #psi_{THROWN};N_{EVENTS}/bin",100,-40,40);
  TH1F *hpsi185 = new TH1F("hpsi185","18.5 < log_{10}(E_{RECON}/eV) <= 19.0;#psi_{RECON} - #psi_{THROWN};N_{EVENTS}/bin",100,-40,40);
  TH1F *hpsi19 = new TH1F("hpsi19","19.0 < log_{10}(E_{RECON}/eV) <= 19.5;#psi_{RECON} - #psi_{THROWN};N_{EVENTS}/bin",100,-40,40);
  TH1F *hpsi195 = new TH1F("hpsi195","19.5 < log_{10}(E_{RECON}/eV) <= 20.0;#psi_{RECON} - #psi_{THROWN};N_{EVENTS}/bin",100,-40,40);
  TH1F *hpsi20 = new TH1F("hpsi20","20.0 < log_{10}(E_{RECON}/eV) <= 20.5;#psi_{RECON} - #psi_{THROWN};N_{EVENTS}/bin",100,-40,40);

  // hpsi18->SetLineColor(kBlack);hpsi18->SetFillColorAlpha(42,0.35);
  // hpsi185->SetLineColor(kBlack);hpsi185->SetFillColorAlpha(46,0.35);
  // hpsi19->SetLineColor(kBlack);hpsi19->SetFillColorAlpha(30,0.35);
  // hpsi195->SetLineColor(kBlack);hpsi195->SetFillColorAlpha(44,0.35);
  // hpsi20->SetLineColor(kBlack);hpsi20->SetFillColorAlpha(32,0.35);

  gStyle->SetOptFit(1);
    
  cpsi18->cd();
  taTree->Draw("57.296*hctim->mpsi[int(etrack->udata[0])]-mc04->psi[6]*57.296>>hpsi18","etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.0 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=18.5");
  hpsi18->Fit("gaus");hpsi18->GetFunction("gaus")->SetLineColor(kRed);hpsi18->GetFunction("gaus")->SetLineWidth(3);
  cpsi18->Update();
  cpsi18->SaveAs("/media/glados/Data_Storage_2/plots/resolution/psi_resolution_hist_18.png");

  cpsi185->cd();
  taTree->Draw("57.296*hctim->mpsi[int(etrack->udata[0])]-mc04->psi[6]*57.296>>hpsi185","etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.0");
  hpsi185->Fit("gaus");hpsi185->GetFunction("gaus")->SetLineColor(kRed);hpsi185->GetFunction("gaus")->SetLineWidth(3);
  cpsi185->Update();
  cpsi185->SaveAs("/media/glados/Data_Storage_2/plots/resolution/psi_resolution_hist_185.png");

  cpsi19->cd();
  taTree->Draw("57.296*hctim->mpsi[int(etrack->udata[0])]-mc04->psi[6]*57.296>>hpsi19","etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.5");
  hpsi19->Fit("gaus");hpsi19->GetFunction("gaus")->SetLineColor(kRed);hpsi19->GetFunction("gaus")->SetLineWidth(3);
  cpsi19->Update();
  cpsi19->SaveAs("/media/glados/Data_Storage_2/plots/resolution/psi_resolution_hist_19.png");

  cpsi195->cd();
  taTree->Draw("57.296*hctim->mpsi[int(etrack->udata[0])]-mc04->psi[6]*57.296>>hpsi195","etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.5 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20");
  hpsi195->Fit("gaus");hpsi195->GetFunction("gaus")->SetLineColor(kRed);hpsi195->GetFunction("gaus")->SetLineWidth(3);
  cpsi195->Update();
  cpsi195->SaveAs("/media/glados/Data_Storage_2/plots/resolution/psi_resolution_hist_195.png");
  
  cpsi20->cd();
  taTree->Draw("57.296*hctim->mpsi[int(etrack->udata[0])]-mc04->psi[6]*57.296>>hpsi20","etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>20 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.5");
  hpsi20->Fit("gaus");hpsi20->GetFunction("gaus")->SetLineColor(kRed);hpsi20->GetFunction("gaus")->SetLineWidth(3);
  cpsi20->Update();
  cpsi20->SaveAs("/media/glados/Data_Storage_2/plots/resolution/psi_resolution_hist_20.png");
  

  
  // TPaveStats *spsi18 = (TPaveStats*)hpsi18->FindObject("stats");
  // spsi18->SetLineColor(42);spsi18->SetTextColor(42); spsi18->Draw();
  // TPaveStats *spsi185 = (TPaveStats*)hpsi185->FindObject("stats");
  // spsi185->SetLineColor(46);spsi185->SetTextColor(46); spsi185->Draw();
  // TPaveStats *spsi19 = (TPaveStats*)hpsi19->FindObject("stats");
  // spsi19->SetLineColor(30);spsi19->SetTextColor(30); spsi19->Draw();
  // TPaveStats *spsi195 = (TPaveStats*)hpsi195->FindObject("stats");
  // spsi195->SetLineColor(44);spsi195->SetTextColor(44); spsi195->Draw();
  // TPaveStats *spsi20 = (TPaveStats*)hpsi20->FindObject("stats");
  // spsi20->SetLineColor(32);spsi20->SetTextColor(32); spsi20->Draw();

  TCanvas *cpsiprof = new TCanvas("cpsiprof","cpsiprof",600,600);
  TProfile *profpsi = new TProfile("profpsi",";log_{10}(E_{RECON}/eV);#psi_{RECON} - #psi_{THROWN}",40,18.0,20.5);
  TH2F *hpsibox = new TH2F("hpsibox",";log_{10}(E_{RECON}/eV);#psi_{RECON} - #psi_{THROWN}",50,17,21,100,-40,40);
  cpsiprof->cd();
  taTree->Draw("57.296*hctim->mpsi[int(etrack->udata[0])]-mc04->psi[6]*57.296:log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>>hpsibox","etrack->qualct==1","colz");
  taTree->Draw("57.296*hctim->mpsi[int(etrack->udata[0])]-mc04->psi[6]*57.296:log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>>profpsi","etrack->qualct==1","profs sames");

  profpsi->SetLineColor(kRed);
  profpsi->SetLineWidth(2);
  cpsiprof->Update();
}

void theta(TTree *taTree)
{
 /////// Theta Resolution Plots ///////
  TCanvas *cthe18 = new TCanvas("cthe18","cthe18",600,600);
  TCanvas *cthe185 = new TCanvas("cthe185","cthe185",600,600);
  TCanvas *cthe19 = new TCanvas("cthe19","cthe19",600,600);
  TCanvas *cthe195 = new TCanvas("cthe195","cthe195",600,600);
  TCanvas *cthe20 = new TCanvas("cthe20","cthe20",600,600);
  
  TH1F *hThe18 = new TH1F("hThe18","18.0 < log_{10}(E_{RECON}/eV) <= 18.5;#theta_{RECON}-#theta_{THROWN};N_{EVENTS}/bin",100,-40,40);
  TH1F *hThe185 = new TH1F("hThe185","18.5 < log_{10}(E_{RECON}/eV) <= 19.0;#theta_{RECON}-#theta_{THROWN};N_{EVENTS}/bin",100,-40,40);
  TH1F *hThe19 = new TH1F("hThe19","19.0 < log_{10}(E_{RECON}/eV) <= 19.5;#theta_{RECON}-#theta_{THROWN};N_{EVENTS}/bin",100,-40,40);
  TH1F *hThe195 = new TH1F("hThe195","19.5 < log_{10}(E_{RECON}/eV) <= 20.0;#theta_{RECON}-#theta_{THROWN};N_{EVENTS}/bin",100,-40,40);
  TH1F *hThe20 = new TH1F("hThe20","20.0 < log_{10}(E_{RECON}/eV) <= 20.5;#theta_{RECON}-#theta_{THROWN};N_{EVENTS}/bin",100,-40,40);

  // hThe18->SetLineColor(kBlack);hThe18->SetFillColorAlpha(42,0.35);
  // hThe185->SetLineColor(kBlack);hThe185->SetFillColorAlpha(46,0.35);
  // hThe19->SetLineColor(kBlack);hThe19->SetFillColorAlpha(30,0.35);
  // hThe195->SetLineColor(kBlack);hThe195->SetFillColorAlpha(44,0.35);
  // hThe20->SetLineColor(kBlack);hThe20->SetFillColorAlpha(32,0.35);

  gStyle->SetOptFit(1);
  
  cthe18->cd();
  taTree->Draw("57.296*hctim->mthe[int(etrack->udata[0])]-mc04->theta*57.296>>hThe18","etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=18.5");
  hThe18->Fit("gaus");hThe18->GetFunction("gaus")->SetLineColor(kRed);hThe18->GetFunction("gaus")->SetLineWidth(3);
  cthe18->Update();
  cthe18->SaveAs("/media/glados/Data_Storage_2/plots/resolution/theta_resolution_hist_18.png");

  cthe185->cd();
  taTree->Draw("57.296*hctim->mthe[int(etrack->udata[0])]-mc04->theta*57.296>>hThe185","etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.0");
  hThe185->Fit("gaus");hThe185->GetFunction("gaus")->SetLineColor(kRed);hThe185->GetFunction("gaus")->SetLineWidth(3);
  cthe185->Update();
  cthe185->SaveAs("/media/glados/Data_Storage_2/plots/resolution/theta_resolution_hist_185.png");
  
  cthe19->cd();
  taTree->Draw("57.296*hctim->mthe[int(etrack->udata[0])]-mc04->theta*57.296>>hThe19","etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.5");
  hThe19->Fit("gaus");hThe19->GetFunction("gaus")->SetLineColor(kRed);hThe19->GetFunction("gaus")->SetLineWidth(3);
  cthe19->Update();
  cthe19->SaveAs("/media/glados/Data_Storage_2/plots/resolution/theta_resolution_hist_19.png");

  cthe195->cd();
  taTree->Draw("57.296*hctim->mthe[int(etrack->udata[0])]-mc04->theta*57.296>>hThe195","etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.5 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20");
  hThe195->Fit("gaus");hThe195->GetFunction("gaus")->SetLineColor(kRed);hThe195->GetFunction("gaus")->SetLineWidth(3);
  cthe195->Update();
  cthe195->SaveAs("/media/glados/Data_Storage_2/plots/resolution/theta_resolution_hist_195.png");

  cthe20->cd();
  taTree->Draw("57.296*hctim->mthe[int(etrack->udata[0])]-mc04->theta*57.296>>hThe20","etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>20 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.5");
  hThe20->Fit("gaus");hThe20->GetFunction("gaus")->SetLineColor(kRed);hThe20->GetFunction("gaus")->SetLineWidth(3);
  cthe20->Update();
  cthe20->SaveAs("/media/glados/Data_Storage_2/plots/resolution/theta_resolution_hist_20.png");
  
  // TPaveStats *sThe18 = (TPaveStats*)hThe18->FindObject("stats");
  // sThe18->SetLineColor(42);sThe18->SetTextColor(42); sThe18->Draw();
  // TPaveStats *sThe185 = (TPaveStats*)hThe185->FindObject("stats");
  // sThe185->SetLineColor(46);sThe185->SetTextColor(46); sThe185->Draw();
  // TPaveStats *sThe19 = (TPaveStats*)hThe19->FindObject("stats");
  // sThe19->SetLineColor(30);sThe19->SetTextColor(30); sThe19->Draw();
  // TPaveStats *sThe195 = (TPaveStats*)hThe195->FindObject("stats");
  // sThe195->SetLineColor(44);sThe195->SetTextColor(44); sThe195->Draw();
  // TPaveStats *sThe20 = (TPaveStats*)hThe20->FindObject("stats");
  // sThe20->SetLineColor(32);sThe20->SetTextColor(32); sThe20->Draw();

  TCanvas *ctheprof = new TCanvas("ctheprof","ctheprof",600,600);
  TProfile *profThe = new TProfile("profThe",";log_{10}(E_{THROWN}/eV);#theta_{RECON} - #theta_{THROWN}",40,18.0,20.5);
  TH2F *hThebox = new TH2F("hThebox",";log_{10}(E_{THROWN}/eV);#theta_{RECON} - #theta_{THROWN}",50,17,21,100,-40,40);
  ctheprof->cd();
  taTree->Draw("57.296*hctim->mthe[int(etrack->udata[0])]-mc04->theta*57.296:log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>>hThebox","etrack->qualct==1","colz");
  taTree->Draw("57.296*hctim->mthe[int(etrack->udata[0])]-mc04->theta*57.296:log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>>profThe","etrack->qualct==1","profs same");

  profThe->SetLineColor(kRed);
  profThe->SetLineWidth(2);
  ctheprof->Update();
}

void rp(TTree *taTree)
{
   /////// Rp resolution plots ///////
  TCanvas *crp18 = new TCanvas("crp18","crp18",600,600);
  TCanvas *crp185 = new TCanvas("crp185","crp185",600,600);
  TCanvas *crp19 = new TCanvas("crp19","crp19",600,600);
  TCanvas *crp195 = new TCanvas("crp195","crp195",600,600);
  TCanvas *crp20 = new TCanvas("crp20","crp20",600,600);
  
  TH1F *hRp18 = new TH1F("hRp18","18.0 < log_{10}(E_{RECON}/eV) <= 18.5;ln(Rp_{RECON}/Rp_{THROWN})",100,-1,1);
  TH1F *hRp185 = new TH1F("hRp185","18.5 < log_{10}(E_{RECON}/eV) <= 19.0;ln(Rp_{RECON}/Rp_{THROWN})",100,-1,1);
  TH1F *hRp19 = new TH1F("hRp19","19.0 < log_{10}(E_{RECON}/eV) <= 19.5;ln(Rp_{RECON}/Rp_{THROWN})",100,-1,1);
  TH1F *hRp195 = new TH1F("hRp195","19.5 < log_{10}(E_{RECON}/eV) <= 20.0;ln(Rp_{RECON}/Rp_{THROWN})",100,-1,1);
  TH1F *hRp20 = new TH1F("hRp20","20.0 < log_{10}(E_{RECON}/eV) <= 20.5;ln(Rp_{RECON}/Rp_{THROWN})",100,-1,1);

  // hRp18->SetLineColor(kBlack);hRp18->SetFillColorAlpha(42,0.35);
  // hRp185->SetLineColor(kBlack);hRp185->SetFillColorAlpha(46,0.35);
  // hRp19->SetLineColor(kBlack);hRp19->SetFillColorAlpha(30,0.35);
  // hRp195->SetLineColor(kBlack);hRp195->SetFillColorAlpha(44,0.35);
  // hRp20->SetLineColor(kBlack);hRp20->SetFillColorAlpha(32,0.35);

  gStyle->SetOptFit(1);
  
  crp18->cd();
  taTree->Draw("log(hctim->mrp[int(etrack->udata[0])]/mc04->rp[6])>>hRp18","etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=18.5");
  hRp18->Fit("gaus");hRp18->GetFunction("gaus")->SetLineColor(kRed);hRp18->GetFunction("gaus")->SetLineWidth(3);
  crp18->Update();
  crp18->SaveAs("/media/glados/Data_Storage_2/plots/resolution/rp_resolution_hist_18.png");
 
  crp185->cd();
  taTree->Draw("log(hctim->mrp[int(etrack->udata[0])]/mc04->rp[6])>>hRp185","etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.0");
  hRp185->Fit("gaus");hRp185->GetFunction("gaus")->SetLineColor(kRed);hRp185->GetFunction("gaus")->SetLineWidth(3);
  crp185->Update();
  crp185->SaveAs("/media/glados/Data_Storage_2/plots/resolution/rp_resolution_hist_185.png");

  crp19->cd();
  taTree->Draw("log(hctim->mrp[int(etrack->udata[0])]/mc04->rp[6])>>hRp19","etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.5");
  hRp19->Fit("gaus");hRp19->GetFunction("gaus")->SetLineColor(kRed);hRp19->GetFunction("gaus")->SetLineWidth(3);
  crp19->Update();
  crp19->SaveAs("/media/glados/Data_Storage_2/plots/resolution/rp_resolution_hist_19.png");

  crp195->cd();
  taTree->Draw("log(hctim->mrp[int(etrack->udata[0])]/mc04->rp[6])>>hRp195","etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.5 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20");
  hRp195->Fit("gaus");hRp195->GetFunction("gaus")->SetLineColor(kRed);hRp195->GetFunction("gaus")->SetLineWidth(3);
  crp195->Update();
  crp195->SaveAs("/media/glados/Data_Storage_2/plots/resolution/rp_resolution_hist_195.png");

  crp20->cd();
  taTree->Draw("log(hctim->mrp[int(etrack->udata[0])]/mc04->rp[6])>>hRp20","etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>20 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.5");
  hRp20->Fit("gaus");hRp20->GetFunction("gaus")->SetLineColor(kRed);hRp20->GetFunction("gaus")->SetLineWidth(3);
  crp20->Update();
  crp20->SaveAs("/media/glados/Data_Storage_2/plots/resolution/rp_resolution_hist_20.png");
  
  // TPaveStats *sRp18 = (TPaveStats*)hRp18->FindObject("stats");
  // sRp18->SetLineColor(42);sRp18->SetTextColor(42); sRp18->Draw();
  // TPaveStats *sRp185 = (TPaveStats*)hRp185->FindObject("stats");
  // sRp185->SetLineColor(46);sRp185->SetTextColor(46); sRp185->Draw();
  // TPaveStats *sRp19 = (TPaveStats*)hRp19->FindObject("stats");
  // sRp19->SetLineColor(30);sRp19->SetTextColor(30); sRp19->Draw();
  // TPaveStats *sRp195 = (TPaveStats*)hRp195->FindObject("stats");
  // sRp195->SetLineColor(44);sRp195->SetTextColor(44); sRp195->Draw();
  // TPaveStats *sRp20 = (TPaveStats*)hRp20->FindObject("stats");
  // sRp20->SetLineColor(32);sRp20->SetTextColor(32); sRp20->Draw();

  TCanvas *crpprof = new TCanvas("crpprof","crpprof",600,600);
  TProfile *profRp = new TProfile("profRp",";log_{10}(E_{RECON}/eV);ln(Rp_{RECON} - Rp_{THROWN})",40,18.0,20.5);
  TH2F *hRpbox = new TH2F("hRpbox",";log_{10}(E_{RECON}/eV);ln(Rp_{RECON} / Rp_{THROWN})",50,17,21,40,-1,1);
  crpprof->cd();
  taTree->Draw("log(hctim->mrp[int(etrack->udata[0])]/mc04->rp[6]):log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>>hRpbox","etrack->qualct==1","colz");
  taTree->Draw("log(hctim->mrp[int(etrack->udata[0])]/mc04->rp[6]):log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>>profRp","etrack->qualct==1","prof same");

  profRp->SetLineColor(kRed);
  profRp->SetLineWidth(2);
  crpprof->Update();
}

void resolution_plots()
{
  // Load TDSTio libs
  gSystem->Load("libTDSTio.so");

  // Open MC file
  TFile *f = TFile::Open("/home/glados/Data_Storage_2/mdtax4/mc/processing/root_all_files/pr.1e17.5-1e21.ps5.all.drmdpcgf.root");
  TTree *taTree = (TTree*)f->Get("taTree");

  //energy(taTree);

  //psi(taTree);

  //theta(taTree);

  rp(taTree);
  
}
