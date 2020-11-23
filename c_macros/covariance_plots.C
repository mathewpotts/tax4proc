Double_t missing_E_corr(Double_t eng_EeV)
{
  Double_t eng_eV = 1e18 * eng_EeV;

  TF1 *JHK = new TF1("JHK","x / (-16.62 + 3.58 * log10(x) - 0.2786 * TMath::Power(log10(x),2) + 0.009775 * TMath::Power(log10(x),3) - 0.0001299 * TMath::Power(log10(x),4))",18.0,21.0);

  Double_t E = JHK->Eval(eng_eV);

  return E; 
}

void covariance_plots()
{
  // Load TDSTio libs
  gSystem->Load("libTDSTio.so");

  // Open MC file
  TFile *f = TFile::Open("/home/glados/Data_Storage_2/mdtax4/mc/processing/root_all_files/pr.1e17.5-1e21.ps5.all.drmdpcgf.root");
  TTree *taTree = (TTree*)f->Get("taTree");

  // // Energy vs. Rp Covariance
  // TCanvas *cerpcoprof = new TCanvas("cerpcoprof","cerpcoprof",600,600);
  // TProfile *proferpco = new TProfile("proferpco",";ln(Rp_{RECON}/Rp_{THROWN});ln(E_{RECON}/E_{THROWN})",20,-0.5,0.5);
  // TH2F *herpcobox = new TH2F("herpcobox",";ln(Rp_{RECON}/Rp_{THROWN});ln(E_{RECON}/E_{THROWN})",50,-1,1,100,-1,1);
  // cerpcoprof->cd();
  // taTree->Draw("(log((1e18*etrack->energy)/mc04->energy)):log(hctim->mrp[int(etrack->udata[0])]/mc04->rp[6])>>herpcobox","etrack->qualct==1","colz");
  // taTree->Draw("(log((1e18*etrack->energy)/mc04->energy)):log(hctim->mrp[int(etrack->udata[0])]/mc04->rp[6])>>proferpco","etrack->qualct==1","prof sames");
  // proferpco->SetLineColor(kRed);
  // proferpco->SetLineWidth(2);
  // cerpcoprof->Update();

  // TCanvas *cerpco2prof = new TCanvas("cerpco2prof","cerpco2prof",600,600);
  // TProfile *proferpco2 = new TProfile("proferpco2",";Rp_{RECON}-Rp_{THROWN};ln(E_{RECON}/E_{THROWN})",20,-5,5);
  // TH2F *herpco2box = new TH2F("herpco2box",";Rp_{RECON}-Rp_{THROWN};ln(E_{RECON}/E_{THROWN})",50,-5,5,100,-1,1);
  // cerpco2prof->cd();
  // taTree->Draw("(log((1e18*etrack->energy)/mc04->energy)):(1e-3*hctim->mrp[int(etrack->udata[0])]-1e-3*mc04->rp[6])>>herpco2box","etrack->qualct==1","colz");
  // taTree->Draw("(log((1e18*etrack->energy)/mc04->energy)):(1e-3*hctim->mrp[int(etrack->udata[0])]-1e-3*mc04->rp[6])>>proferpco2","etrack->qualct==1","prof sames");
  // proferpco2->SetLineColor(kRed);
  // proferpco2->SetLineWidth(2);
  // cerpco2prof->Update();

  // Energy vs. Psi Covariance
  TCanvas *cepsico2prof = new TCanvas("cepsico2prof","cepsico2prof",600,600);
  TProfile *profepsico2 = new TProfile("profepsico2",";#psi_{RECON};ln(E_{RECON}/E_{THROWN})",30,0,120);
  TH2F *hepsico2box = new TH2F("hepsico2box",";#psi_{RECON};ln(E_{RECON}/E_{THROWN})",100,0,120,100,-1,1);
  cepsico2prof->cd();
  taTree->Draw("(log((1e18*etrack->energy)/mc04->energy)):(57.296*hctim->mpsi[int(etrack->udata[0])])>>hepsico2box","","colz");
  taTree->Draw("(log((1e18*etrack->energy)/mc04->energy)):(57.296*hctim->mpsi[int(etrack->udata[0])])>>profepsico2","","prof sames");
  profepsico2->SetLineColor(kRed);
  profepsico2->SetLineWidth(2);
  cepsico2prof->Update();

  TCanvas *cepsico3prof = new TCanvas("cepsico3prof","cepsico3prof",600,600);
  TProfile *profepsico3 = new TProfile("profepsico3",";#psi_{THROWN};ln(E_{RECON}/E_{THROWN})",30,0,120);
  TH2F *hepsico3box = new TH2F("hepsico3box",";#psi_{THROWN};ln(E_{THROWN}/E_{THROWN})",100,0,120,100,-1,1);
  cepsico3prof->cd();
  taTree->Draw("(log((1e18*etrack->energy)/mc04->energy)):(mc04->psi[6]*57.296)>>hepsico3box","","colz");
  taTree->Draw("(log((1e18*etrack->energy)/mc04->energy)):(mc04->psi[6]*57.296)>>profepsico3","","prof sames");
  profepsico3->SetLineColor(kRed);
  profepsico3->SetLineWidth(2);
  cepsico3prof->Update();


  //Energy Vs. XM




}
