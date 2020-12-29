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

Double_t weight(Double_t log10en)
{
  return TMath::Power(10.0, 2*(log10en-log10en_min))*spectrum_function(log10en);
}

Int_t bracketing_cut(Double_t dep_first, Double_t dep_last, Double_t xmax, Double_t dx1, Double_t dx2)
{
  Double_t lower_bracket = dep_first + dx1;
  Double_t upper_bracket = dep_last - dx2;
  if (lower_bracket < xmax && upper_bracket > xmax){
    return 1;
  }
  return 0;
}

Double_t get_border_distance(Int_t fdsiteid,
			     Int_t sd_det_id,
			     Double_t nx, Double_t ny, Double_t nz,
			     Double_t rp, Double_t psi)
{
  // Convert FD core coordinates to CLF coordinates
  Double_t sdp_n[3] = {nx, ny, nz};
  tlevent_track t_fd(sdp_n,rp,psi,0);
  tlevent_track t_clf = TLUTI_REC::fd2clf(fdsiteid,t_fd);

  Double_t xcore = t_clf.xyz[0];
  Double_t ycore = t_clf.xyz[1];

  Double_t bdist = TLUTI_DETECTORS::get_sd_bdist(TASD_MAIN, xcore, ycore);
  
  return bdist;
}

// Check that the calculated core positions between the FD and SD are within a certain tolerance
Int_t is_core_position_diff(Int_t fdsiteid,
			     Double_t xcore_SD, Double_t ycore_SD,
			     Double_t nx, Double_t ny, Double_t nz,
			     Double_t rp, Double_t psi,
			     Int_t tolerance)
{
  // Convert FD core coordinates to CLF coordinates 
  Double_t sdp_n[3] = {nx, ny, nz};
  tlevent_track t_fd(sdp_n,rp,psi,0);
  tlevent_track t_clf = TLUTI_REC::fd2clf(fdsiteid,t_fd);

  Double_t xcore = t_clf.xyz[0];
  Double_t ycore = t_clf.xyz[1];
  
  Double_t xcore_diff = abs(xcore-xcore_SD);
  Double_t ycore_diff = abs(ycore-ycore_SD);

  if (xcore_diff < tolerance && ycore_diff < tolerance)
    {
      return 1;
    }
  
  return 0;
}

// Cut out events that have a certain percentage of saturated tubes
Int_t event_saturation_cut(vector<int> tubeArray,Int_t saturation_limit){
  Double_t nentries = tubeArray.size();
  Double_t saturation_count = 0.0;

  for(Int_t i = 0;i < nentries;i++){
      if (tubeArray[i] == -2){
  	  saturation_count += 1.0;
  	}
    }
  
  Double_t saturation_percent = saturation_count / nentries;
  Double_t saturation_plimit  = saturation_limit / 100.;

  if (saturation_percent <= saturation_plimit){
      return 1;
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

bool good_weather_cut(Int_t weatcode,Int_t mir)
{
  Int_t n = weatcode / 1000000 % 10;
  Int_t e = weatcode / 100000 % 10;
  Int_t s = weatcode / 10000 % 10;
  Int_t w = weatcode / 1000 % 10;
  Int_t o = weatcode / 100 % 10;
  Int_t t = weatcode / 10 % 10;
  Int_t h = weatcode % 10;

  if (mir >= 25 && mir <= 28){
    if (o+t+n+e == 0){ // MDTAx4 weather codes
      return 1;
    }
  }
  if (mir >= 29 && mir <= 36){ 
    if (o+t+s+e == 0){ // BRTAx4 weather codes
      return 1;
    }
  }
  return 0;
}
