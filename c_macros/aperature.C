void aperature()
{
  
  TFile *g=new TFile("/home/glados/Data_Storage_2/mdtax4/mc/processing/hErecon_E3.root","r");
  TH1D *hErecon = (TH1D *)g->Get("hErecon");
  hErecon->Sumw2();

  TFile *f = new TFile("/home/glados/Data_Storage_2/mdtax4/mc/processing/hEthrown_E3.root","r");
  TH1D *hEthrown = (TH1D *)f->Get("hEthrown");
  hEthrown->Sumw2();
  
  TH1D *hEratio=(TH1D*)hErecon->Clone("hEratio");
  hEratio->Sumw2();
  hEratio->Divide(hEthrown);

  const Double_t rpmin = 100.0;
  const Double_t rpmax = 50.0e3;
  const Double_t theta_max = 70.0;

  Double_t  area = TMath::Pi() * (rpmax*rpmax - rpmin*rpmin);
  Double_t omega = TMath::TwoPi() * (1 - TMath::Cos(TMath::DegToRad()*theta_max));
  Double_t area_times_omega = area * omega;

  // fprintf(stderr,"CHECK: %.3e %.3e\n", x, area_times_omega);
  
  hEratio->Scale(area_times_omega);
  hEratio->Draw();
  hEratio->GetYaxis()->SetTitle("Reconstructed Aperture (m^{2}sr)");

  TF1 *func=new TF1("fit_func","TMath::Power(10,[0]*(1-TMath::Exp(([1]-x)/[2])))",18.,20.5);
  func->SetParameters(9,16,1);
  hEratio->Fit("fit_func","","",18.5,20.5);

  TH1D *hAp=(TH1D*)hEratio->Clone("hAp");
  
  hAp->Draw();
  gPad->SetLogy();
  hAp->SaveAs("/home/glados/Data_Storage_2/mdtax4/mc/processing/hAp.root");
  //hAp->SaveAs("/home/glados/Data_Storage_2/mdtax4/mc/processing/hAp.C");

  // const Int_t n = 26;
  // Float_t x[n]={18.05,18.15,18.25,18.35,18.45,18.55,18.65,18.75,18.85,18.95,19.05,19.15,19.25,19.35,19.45,19.55,19.65,19.75,19.85,19.95,20.05,20.15,20.25,20.35,20.45};
  // Float_t y[n]{68415665.7611779,110848062.989398,191775201.735584,286192474.850177,413658939.349657,515935262.625619,689092271.961884,887683982.339774,998963984.494787,992783535.159917,1170167839.75128,1481980190.67113,1610491876.43851,1711108125.41266,2034043345.60373,1909335223.56464,1939320172.31519,2455500971.79196,2833579977.58436,3048749509.97625,5199346495.05097,4730340259.35345,4059836884.62511,5285597053.29796,13011888720.8389};
  // Float_t ey[n]={4581451.84340906,6941574.96636977,10771168.4325757,14979998.7801881,20657141.6586599,26260413.4400607,34031841.5919631,43522287.9601371,50519780.1862879,54404202.5979076,66568522.5765147,84719168.6364415,98011414.3782437,114328220.661057,137447911.790919,150946203.081276,165687303.998207,217890477.119719,286234802.454903,367026310.603312,655056085.998104,767362519.304377,667433072.457305,1126893078.13023,2774143996.42233};
  // TGraphErrors *gr = new TGraphErrors(n,x,y,0,ey);
  // gr->Draw("same");
}
