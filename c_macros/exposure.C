const Double_t gwontime = 294.1825;
//const Double_t gwontime = 222.6961 * 2;
//const Double_t gwontime = 750;

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
  // else if (o < 2)
  //   {
  //     if (n + e < 2)
  // 	{
  // 	  return 1;
  // 	}
  //     else if (o + n + e < 3)
  // 	{
  // 	  return 1;
  // 	}
  //   }
  
  return 0;
}

// double ontime_calc()
// {
//   TFile *f = new TFile("/media/glados/Data_Storage_2/mdtax4/pass5/data/processing/root_all_files/20190625-20200215.mdtax4.ps0.root","r");
//   TTree *taTree = (TTree*)f->Get("taTree");
  
//   tlweat_class *tlweat = new tlweat_class;
//   //etrack_class *etrack = new etrack_class;
  
//   taTree->SetBranchAddress("tlweat",&tlweat);
//   //taTree->SetBranchAddress("etrack",&etrack);
  
//   const Int_t nentries = taTree->GetEntries();
//   Double_t ontime_gw =0;

//   // Set up maps
//   map<Long64_t,Int_t> eventmap;
//   map<Long64_t,Int_t> partdur;
//   map<Long64_t,Int_t> eventmap_gw;
  
//   for (Int_t i=0;i < nentries; i++)
//     {
//       taTree->GetEntry(i);
//       if (tlweat->yymmdd <= 191125 )
// 	{
// 	  Long64_t date_part = ((Long64_t)tlweat->yymmdd)*(Long64_t)1000+((Long64_t)tlweat->part_num);
      
// 	  partdur[date_part] = tlweat->dt;
// 	  eventmap[date_part] += 1;

	  
// 	  if(good_weather_cut(tlweat->code))
// 	    {
// 	      eventmap_gw[date_part] += 1;
// 	    }
// 	}
//     }
//   for (auto i0=eventmap_gw.begin(); i0!=eventmap_gw.end(); ++i0)
//     {
//       Long64_t date         = i0->first;
//       Int_t nevents_part    = eventmap[date];    
//       Int_t nevents_gw_part = i0->second;
//       Double_t dt           = partdur[date];
//       Double_t f  = (Double_t)nevents_gw_part / (Double_t)nevents_part;	
//       ontime_gw += f*dt;

//       std::cout<<f<<"*"<<dt<<"="<< ontime_gw<<endl;
//     }
//   Double_t totalontime;
//   for (auto i1=eventmap.begin();i1!=eventmap.end();++i1)
//     {
//       Long64_t date         = i1->first;
//       Double_t ontime       = partdur[date];
      
//       totalontime += ontime;
//       std::cout<<totalontime<<endl;
//     }
//   return ontime_gw;
// }

void exposure()
{
  // Grabing the Aperture histogram
  TFile *f = new TFile("/media/glados/Data_Storage_2/mdtax4/mc/processing/hAp.root","r");
  TH1D *hAp = (TH1D*)f->Get("hAp");
  
  
  // Scale histogram
  TH1D *hExposure = (TH1D*)hAp->Clone("hExposure");
  
  //hExposure->Scale(ontime_calc());
  hExposure->Scale(gwontime*3600);
  
  TF1 *func=new TF1("fit_func","TMath::Power(10,[0]*(1-TMath::Exp(([1]-x)/[2])))",18.,20.5);
  func->SetParameters(9,16,1);
  hExposure->Fit("fit_func","","",18.0,20.5);
  hExposure->GetYaxis()->SetTitle("Exposure (m^{2}sr s)");
  hExposure->Draw();
  hExposure->SaveAs("/media/glados/Data_Storage_2/mdtax4/mc/processing/hExposure.root");

  // const Int_t n = 15;
  // Double_t x[n] = {18.55,18.65,18.75,18.85,18.95,19.05,19.15,19.25,19.35,19.45,19.55,19.65,19.75,19.85,19.95};
  // Double_t y[n] = {679327316881789,907321592661765,1168805510390939,1.31532688770843E+015,1.30718914569325E+015,1.54074946308995E+015,1.95130997923496E+015,2120520159279694,2253000668757720,2678206566766324,2.51400450485055E+015,2.55348541700586E+015,3233136030049725,3.73094925426552E+015,4.01426104104015E+015};
  // Double_t ey[n] = {34576850032806.9,44809419537563.4,57305404856520.9,66518939893263.6,71633523922841.6,87650174559793.8,111548966871085,129050747209104,150534939190428,180976428423896,198749507081289,218158981988704,286894430209240,376882788279149,483260239934586};
  // TGraphErrors *gr = new TGraphErrors(n,x,y,0,ey);
  // gr->Draw("same");
}

  
