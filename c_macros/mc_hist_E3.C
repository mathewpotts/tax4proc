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

Double_t missing_E_corr(Double_t eng_EeV)
{
  Double_t eng_eV = 1e18 * eng_EeV;

  TF1 *JHK = new TF1("JHK","x / (-16.62 + 3.58 * log10(x) - 0.2786 * TMath::Power(log10(x),2) + 0.009775 * TMath::Power(log10(x),3) - 0.0001299 * TMath::Power(log10(x),4))",18.0,21.0);

  Double_t E = JHK->Eval(eng_eV);

  return E;  
}

Double_t spectrum_function(Double_t log10en)
{
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

Double_t weight(Double_t log10en)
{
  return TMath::Power(10.0, 2*(log10en-log10en_min))*spectrum_function(log10en);
}

void mc_hist_E3()
{
  // Concatinate files
  system("dstcat -o /home/glados/Data_Storage_2/mdtax4/mc/processing/concat_dst/pr.1e17-1e21.ps5.all.dst.gz /home/glados/Data_Storage_2/mdtax4/mc/*/pass5/*.dst");
  system("dstcat -o /home/glados/Data_Storage_2/mdtax4/mc/processing/concat_dst/pr.1e17-1e21.ps0.all.dst.gz /home/glados/Data_Storage_2/mdtax4/mc/*/pass0/*.dst");

  // Apply DCR cuts to ps5 file
  system("dcrmdpcgf.run -f -fd TAX4_MIDDLE_DRUM /home/glados/Data_Storage_2/mdtax4/mc/processing/concat_dst/pr.1e17-1e21.ps5.all.dst.gz -o /home/glados/Data_Storage_2/mdtax4/mc/processing/concat_drmdpcgf_dst/");

  // Convert the DST files into ROOT trees
  system("dst2rt_ta.run -f /home/glados/Data_Storage_2/mdtax4/mc/processing/concat_drmdpcgf_dst/pr.1e17-1e21.ps5.all.drmdpcgf.dst.gz -o /home/glados/Data_Storage_2/mdtax4/mc/processing/root_all_files/");
  system("dst2rt_ta.run -f /home/glados/Data_Storage_2/mdtax4/mc/processing/concat_dst/pr.1e17-1e21.ps0.all.dst.gz -o /home/glados/Data_Storage_2/mdtax4/mc/processing/root_all_files/");

  // Load Dmitri's TDSTio libs
  gSystem->Load("libTDSTio.so");

  // // Get Reconstructed histogram and save it to root file
  TFile *f = new TFile("/home/glados/Data_Storage_2/mdtax4/mc/processing/root_all_files/pr.1e17-1e21.ps5.all.drmdpcgf.root","r");
  TTree *taTree1 = (TTree*)f->Get("taTree");
  TH1D *hErecon = new TH1D("hErecon",";log_{10}(E/eV);N_{EVENTS} / BIN",40,17,21);
  //TH1D *hErecon_E2 = new TH1D("hErecon_E2",";log_{10}(E/eV);N_{EVENTS} / BIN",40,17,21);
  taTree1->Draw("log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>>hErecon","(etrack->qualct==1) * weight(log10(mc04->energy))");
  //taTree1->Draw("log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>>hErecon_E2","(etrack->qualct==1)","sames");
  hErecon->SaveAs("/home/glados/Data_Storage_2/mdtax4/mc/processing/hErecon_E3.root"); 
  
  // Get Thrown histogram and save it to root file
  TFile *g = new TFile("/home/glados/Data_Storage_2/mdtax4/mc/processing/root_all_files/pr.1e17-1e21.ps0.all.root","r");
  TTree *taTree2 = (TTree*)g->Get("taTree");
  TH1D *hEthrown = new TH1D("hEthrown",";log_{10}(E/eV);N_{EVENTS} / BIN",40,17,21);
  //TH1D *hEthrown_E2 = new TH1D("hEthrown_E2",";log_{10}(E/eV);N_{EVENTS} / BIN",40,17,21);
  taTree2->Draw("log10(mc04->energy)>>hEthrown", "weight(log10(mc04->energy))","sames");
  //taTree2->Draw("log10(mc04->energy)>>hEthrown_E2", "","sames");
  hEthrown->SaveAs("/home/glados/Data_Storage_2/mdtax4/mc/processing/hEthrown_E3.root");

}
