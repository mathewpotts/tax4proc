const Int_t nxbins        = 20;
const Double_t log10en_lo = 18.0;
const Double_t log10en_up = 20.0;

bool good_weather_cut(Int_t weatcode,Int_t mir)
{
  Int_t n = weatcode / 1000000 % 10;
  Int_t e = weatcode / 100000 % 10;
  Int_t s = weatcode / 10000 % 10;
  Int_t w = weatcode / 1000 % 10;
  Int_t o = weatcode / 100 % 10;
  Int_t t = weatcode / 10 % 10;
  Int_t h = weatcode % 10;
  //printf("%d%d%d%d%d%d%d\n",n,e,s,w,o,t,h);
  if (mir >= 25 && mir <= 28){
    if (o+t+n+e == 0){ // MDTAx4 weather codes
      return 1;
    }
  }
  if (mir >= 29 && mir <= 36){ 
    if (o+t+s+w == 0){ // BRTAx4 weather codes
      return 1;
    }
  }
  return 0;
}

Double_t missing_E_corr(Double_t eng_EeV)
{
  Double_t eng_eV = 1e18 * eng_EeV;

  TF1 *JHK = new TF1("JHK","x / (-16.62 + 3.58 * log10(x) - 0.2786 * TMath::Power(log10(x),2) + 0.009775 * TMath::Power(log10(x),3) - 0.0001299 * TMath::Power(log10(x),4))",18.0,21.0);

  Double_t E = JHK->Eval(eng_eV);

  return E; 
}

void data_hist()
{
  // // Concatinate files
  //system("dstcat -o /home/glados/Data_Storage_2/mdtax4/pass5/data/processing/concat_dst/20190625-20200215.mdtax4.ps5.ps2g.dst.gz /home/glados/Data_Storage_2/mdtax4/pass5/data/*/*ps2g.dst.gz");
  //system("dstcat -o /home/glados/Data_Storage_2/mdtax4/pass5/data/processing/concat_dst/20190625-20200215.mdtax4.ps0.dst.gz /home/glados/Data_Storage_2/mdtax4/pass0/data/*/*fraw1.dst.gz");
  
  // // Apply DCR to ps5 file
  //system("dcrmdpcgf.run -f -fd TAX4_MIDDLE_DRUM /home/glados/Data_Storage_2/mdtax4/pass5/data/processing/concat_dst/20190625-20200215.mdtax4.ps5.ps2g.dst.gz -o /home/glados/Data_Storage_2/mdtax4/pass5/data/processing/concat_drmdpcgf_dst/");

  // // Convert the DST files into ROOT trees
  //system("dst2rt_ta.run -f /home/glados/Data_Storage_2/mdtax4/pass5/data/processing/concat_drmdpcgf_dst/20190625-20200215.mdtax4.ps5.ps2g.drmdpcgf.dst.gz -o /home/glados/Data_Storage_2/mdtax4/pass5/data/processing/root_all_files/");
  // system("dst2rt_ta.run -f /home/glados/Data_Storage_2/mdtax4/pass5/data/processing/concat_dst/20190625-20200215.mdtax4.ps0.dst.gz -o /home/glados/Data_Storage_2/mdtax4/pass5/data/processing/root_all_files/");
  
  // Load Dmitri's TDSTio libs
  gSystem->Load("libTDSTio.so");

  // Get Reconstructed histogram and save it to root file
  TFile *f = new TFile("/media/glados/Data_Storage_2/mdtax4/pass5/data/processing/root_all_files/20190625-20200215.mdtax4.ps5.ps2g.drmdpcgf.root","r");
  TTree *taTree = (TTree*)f->Get("taTree");
  TH1D *hEdata = new TH1D("hEdata",";log_{10}(E/eV);N_{EVENTS} / BIN",nxbins,log10en_lo,log10en_up);
  taTree->Draw("log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>>hEdata","good_weather_cut(tlweat->code,fscn1->mir[etrack->udata[0]])==1 && etrack.qualct==1 && log10(missing_E_corr(prfc->eng[etrack->udata[0]]))>=18.5");
  hEdata->SaveAs("/media/glados/Data_Storage_2/mdtax4/pass5/data/processing/hEdata.root");
}
