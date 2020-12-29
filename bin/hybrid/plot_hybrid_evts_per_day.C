void plot_hybrid_evts_per_day(){

  TFile *f2 = new TFile("/home/potts/data/hyb_mdtax4/pass5/hybrid_190625_201207.tlhbgeomnp.root","r");
  TTree *taTree = (TTree*)f2->Get("taTree");

  talex00_class *talex00 = new talex00_class;
  taTree->SetBranchAddress("talex00",&talex00);

  TDatime X0(2019,6,25,0,0,0);
  auto T0 = X0.Convert();
  gStyle->SetTimeOffset(T0);
  
  TDatime X1(2019,6,25,0,0,0);
  auto T1 = X1.Convert()-T0;
  
  TDatime X2(2020,12,30,0,0,0);
  auto T2 = X2.Convert(1)-T0;

  auto h1 = new TH1F("h1","test",17,T1,T2);

  Int_t nentries=taTree->GetEntries();
  Int_t date;
  Int_t yr;
  Int_t mn;
  Int_t d;
  
  for (Int_t i=0;i<nentries;i++){
    taTree->GetEntry(i);

    date = talex00->yymmdd;
    yr   = 2000 + date / 10000;
    mn   = date / 100 % 100;
    d    = date % 4;
    //cout << yr << " " << " " << mn << " " << d << endl;
    
    TDatime x(yr,mn,d,0,0,0);
    auto t = x.Convert()-T0;
    h1->Fill(t);
    
  }
  
  h1->GetXaxis()->SetTimeDisplay(1);
  h1->GetXaxis()->SetTimeFormat("%Y/%m/%d");
  h1->GetXaxis()->SetLabelSize(0.03);
  h1->GetXaxis()->SetLabelOffset(0.01);
  h1->SetTitle("# of Hybrid Events per Day");
  h1->Draw();

}
