void mc_hist()
{
  // Concatinate files
  // system("dstcat -o /home/glados/Data_Storage_2/mdtax4/mc/processing/concat_dst/pr.1e17.5-1e21.ps5.all.dst.gz /home/glados/Data_Storage_2/mdtax4/mc/*/pass5/*.dst");
  // system("dstcat -o /home/glados/Data_Storage_2/mdtax4/mc/processing/concat_dst/pr.1e17.5-1e21.ps0.all.dst.gz /home/glados/Data_Storage_2/mdtax4/mc/*/pass0/*.dst");

  // // Apply DCR cuts to ps5 file
  // system("dcrmdpcgf.run -fd TAX4_MIDDLE_DRUM /home/glados/Data_Storage_2/mdtax4/mc/processing/concat_dst/pr.1e17.5-1e21.ps5.all.dst.gz -o /home/glados/Data_Storage_2/mdtax4/mc/processing/concat_drmdpcgf_dst/");

  // // Convert the DST files into ROOT trees
  // system("dst2rt_ta.run /home/glados/Data_Storage_2/mdtax4/mc/processing/concat_drmdpcgf_dst/pr.1e17.5-1e21.ps5.all.drmdpcgf.dst.gz -o /home/glados/Data_Storage_2/mdtax4/mc/processing/root_all_files/");
  // system("dst2rt_ta.run /home/glados/Data_Storage_2/mdtax4/mc/processing/concat_dst/pr.1e17.5-1e21.ps0.all.dst.gz -o /home/glados/Data_Storage_2/mdtax4/mc/processing/root_all_files/");

  // Load Dmitri's TDSTio libs
  gSystem->Load("libTDSTio.so");

  // Get Reconstructed histogram and save it to root file
  TFile *f = new TFile("/home/glados/Data_Storage_2/mdtax4/mc/processing/root_all_files/pr.1e17.5-1e21.ps5.all.drmdpcgf.root","r");
  TTree *taTree1 = (TTree*)f->Get("taTree");
  TH1D *hErecon = new TH1D("",";log_{10}(E/eV);N_{EVENTS} / BIN",40,17,21);
  taTree1->Draw("18+log10(etrack.energy)>>hErecon","etrack.qualct==1");
  hErecon->SaveAs("/home/glados/Data_Storage_2/mdtax4/mc/processing/hErecon.root");
  hErecon->SaveAs("/home/glados/Data_Storage_2/mdtax4/mc/processing/hErecon.C");  
  
  
  // Get Thrown histogram and save it to root file
  TFile *g = new TFile("/home/glados/Data_Storage_2/mdtax4/mc/processing/root_all_files/pr.1e17.5-1e21.ps0.all.root","r");
  TTree *taTree2 = (TTree*)g->Get("taTree");
  TH1D *hEthrown = new TH1D("hEthrown",";log_{10}(E/eV);N_{EVENTS} / BIN",40,17,21);
  taTree2->Draw("log10(mc04.energy)>>hEthrown");
  hEthrown->SaveAs("/home/glados/Data_Storage_2/mdtax4/mc/processing/hEthrown.root");
  hErecon->SaveAs("/home/glados/Data_Storage_2/mdtax4/mc/processing/hEthrown.C"); 
}
