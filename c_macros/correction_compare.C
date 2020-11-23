void correction_compare()
{
  TFile *f = new TFile("/media/glados/Data_Storage_2/mdtax4/pass5/data/processing/root_all_files/20190625-20200215.mdtax4.ps5.ps2g.drmdpcgf.root","r");
  TTree *taTree = (TTree*)f->Get("taTree");

  // DCR
  TF1 *fdcr = new TF1("fdcr","log(((TMath::Power(10,x))*(1.086-0.01896*log10(TMath::Power(10,x))))/TMath::Power(10,x))*100",-0.5,2.5); 
  TCanvas *dcr = new TCanvas("dcr","dcr",800,600);
  dcr->cd();  
  taTree->Draw("log(etrack->energy/prfc->eng[etrack->udata[0]])>>DCR","etrack->qualct==1");

  // Tom's
  TF1 *ftom = new TF1("ftom","log(((TMath::Power(10,x))/(-.5717 + .1416 * log10(TMath::Power(10,x)*1e18) +-.003328 * TMath::Power(log10(TMath::Power(10,x)*1e18),2)))/TMath::Power(10,x))*100",-0.5,2.5);
  TCanvas *tom = new TCanvas("tom","tom",800,600);
  tom->cd();
  taTree->Draw("log(((1e18*prfc->eng[etrack->udata[0]])/(-.5717 + .1416 * log10(1e18*prfc->eng[etrack->udata[0]]) -.003328 * TMath::Power(log10(1e18*prfc->eng[etrack->udata[0]]),2)))/(1e18*prfc->eng[etrack->udata[0]]))>>TOM","etrack->qualct==1");

  // JiHee
  TF1 *fjhk = new TF1("fjhk","log(((TMath::Power(10,x))/(-16.62 + 3.58 * log10(1e18*TMath::Power(10,x)) + -.2786  * TMath::Power(log10(1e18*TMath::Power(10,x)),2) + .009775 * TMath::Power(log10(1e18*TMath::Power(10,x)),3) - .0001299 * TMath::Power(log10(1e18*TMath::Power(10,x)),4)))/TMath::Power(10,x))*100",-0.5,2.5);
  TCanvas *jhk = new TCanvas("jhk","jhk",800,600);
  jhk->cd();
  taTree->Draw("log(((1e18*prfc->eng[etrack->udata[0]])/(-16.62 + 3.58 * log10(1e18*prfc->eng[etrack->udata[0]]) + -.2786  * TMath::Power(log10(1e18*prfc->eng[etrack->udata[0]]),2) + .009775 * TMath::Power(log10(1e18*prfc->eng[etrack->udata[0]]),3) - .0001299 * TMath::Power(log10(1e18*prfc->eng[etrack->udata[0]]),4)))/(1e18*prfc->eng[etrack->udata[0]]))>>JHK","etrack->qualct==1");

  TCanvas *eq_comp = new TCanvas("eq_comp","eq_comp",800,600);
  eq_comp->cd();
  fdcr->SetLineColor(kGreen+2);
  fdcr->Draw();
  fdcr->GetYaxis()->SetTitle("ln(E_{CORRECTED} / E_{CAL})");
  fdcr->GetXaxis()->SetTitle("log_{10}(E_{CAL}/EeV)");
  fdcr->SetTitle("");
  ftom->SetTitle("");
  ftom->Draw("same");
  fjhk->SetLineColor(kBlue);
  fjhk->SetTitle("");
  fjhk->Draw("same");
  eq_comp->BuildLegend();
}
