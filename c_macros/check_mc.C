const Double_t en_min = TMath::Power(10.0,17.0);    // Minimum Thrown Energy
const Double_t log10en_max = 21.0;

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

Double_t e3j_n(Double_t log10en)
{
return TMath::Power(10.0,3*(log10en-log10en_min)) * spectrum_function(log10en);
}

void check_mc()
{
  TFile *f = new TFile("/home/glados/Data_Storage_2/mdtax4/mc/processing/hEthrown_E3.root","r");
  TH1D *hEthrown = (TH1D*)f->Get("hEthrown");
  
  Int_t nxbins = hEthrown->GetNbinsX();

  TGraphAsymmErrors *gE3thrown = new TGraphAsymmErrors(0);
  
  Int_t npts = 0;
  for (int i=1 ; i <= nxbins ; i++)
    {
      // Defining variables for each bin
      Double_t log_bin_center = hEthrown->GetBinCenter(i);
      Double_t log_bin_width = hEthrown->GetBinWidth(i);
      Double_t log_bin_lower_edge = (log_bin_center - log_bin_width / 2.0);
      Double_t log_bin_upper_edge = (log_bin_center + log_bin_width / 2.0);
      Double_t bin_center = TMath::Power(10,log_bin_center);
      Double_t bin_lower_edge = TMath::Power(10,log_bin_lower_edge);
      Double_t bin_upper_edge = TMath::Power(10,log_bin_upper_edge);

      
      // Each piece of the spectrum calc for each bin
      Double_t nevent = (hEthrown->GetBinContent(i));
      Double_t Ec3 = (TMath::Power(bin_center/en_min,3));
      Double_t dE = (bin_upper_edge - bin_lower_edge);

      gE3thrown->SetPoint(npts, log_bin_center, Ec3 * nevent / dE );
     
      npts ++;
    }
  gE3thrown->SetMarkerStyle(21);
  gE3thrown->Draw("a,e1p");

  
  TF1 *fE3J = new TF1("fE3J", "[0] * e3j_n(x)",log10en_min,log10en_max);
  fE3J->SetParameter(0,2.1e-10);
  fE3J->Draw("same");
}
