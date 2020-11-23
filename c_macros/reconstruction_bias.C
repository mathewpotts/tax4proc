Double_t missing_E_corr(Double_t eng_EeV)
{
  Double_t eng_eV = 1e18 * eng_EeV;

  TF1 *JHK = new TF1("JHK","x / (-16.62 + 3.58 * log10(x) - 0.2786 * TMath::Power(log10(x),2) + 0.009775 * TMath::Power(log10(x),3) - 0.0001299 * TMath::Power(log10(x),4))",18.0,21.0);

  Double_t E = JHK->Eval(eng_eV);

  return E;   
}

Double_t bias_corr(Double_t eng_eV)
{
  TF1 *bias = new TF1("bias","x/(1.01 - 16.87 + 1.751 * log10(x) - 0.04523 * TMath::Power(log10(x),2))",18.0,21.0);

  Double_t E = bias->Eval(eng_eV);

  return E;
}

void reconstruction_bias()
{
  TFile *f = new TFile("/home/glados/Data_Storage_2/mdtax4/mc/processing/root_all_files/pr.1e17.5-1e21.ps5.all.drmdpcgf.root","r");
  TTree *taTree = (TTree*)f->Get("taTree");
   
  TCanvas *Ebias = new TCanvas("Ebias","Ebias",800,600);
  Ebias->cd();
		 
  TProfile *pEbias = new TProfile("pEbias",";log_{10}(E_{RECON}/EeV);ln(E_{RECON}/E_{THROWN})",10,18,20.5);
  taTree->Draw("log((missing_E_corr(prfc->eng[etrack->udata[0]]))/mc04->energy):log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>>pEbias","(etrack->qualct==1)","prof");
  pEbias->Fit("pol0","","",18,20.5);
  //Ebias->SetLogy();

  // TCanvas *Rpbias = new TCanvas("Rpbias","Rpbias",800,600);
  // Rpbias->cd();
  // TProfile *pRpbias = new TProfile("pRpbias",";log_{10}(E_{RECON}/EeV);Rp_{THROWN}/RP_{RECON}",30,-0.5,2.5);
  // taTree->Draw("(mc04->rp[6]/hctim->mrp[int(etrack->udata[0])]):log10(prfc->eng[etrack->udata[0]])>>pRpbias","(etrack->qualct==1)","prof");
  // pRpbias->Fit("pol1","","",0,2.5);

}
