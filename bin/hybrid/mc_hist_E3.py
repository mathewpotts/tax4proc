#!/usr/bin/env python2

# Import Libs
import ROOT
import os
import argparse

def read_in_args():
    parser = argparse.ArgumentParser(description = 'Concatinate data files (optional) and plot the mc energy histogram.')
    parser.add_argument('-det',metavar='detector',action='store',help='Detector that you want to produce a data histogram for. (brtax4/mdtax4)',required=True)
    parser.add_argument('-concat',action='store_true',help='Use this flag if you want to concatinate all the data files into a single file.',default=False)
    args = parser.parse_args()
    concat=args.concat
    det=args.det
    return concat,det

def concat_dst(det,date):
    if det == 'mdtax4':
        fd = 'TAX4_MIDDLE_DRUM'
    if det == 'brtax4':
        fd = 'TAX4_BLACK_ROCK'

    # Generate directories
    os.system("mkdir -p /home/glados/Data_Storage_2/{0}/mc/processing/concat_dst/".format(det))
    os.system("mkdir -p /home/glados/Data_Storage_2/{0}/mc/processing/concat_drmdpcgf_dst/".format(det))
    os.system("mkdir -p /home/glados/Data_Storage_2/{0}/mc/processing/root_all_files/".format(det))
        
    # Concatinate files
    os.system("dstcat -o /home/glados/Data_Storage_2/{0}/mc/processing/concat_dst/pr.1e17-1e21.20190625-{1}.{0}.ps5.dst.gz /home/glados/Data_Storage_2/{0}/mc/*/pass5/*.dst".format(det,date))
    os.system("dstcat -o /home/glados/Data_Storage_2/{0}/mc/processing/concat_dst/pr.1e17-1e21.20190625-{1}.{0}.ps0.dst.gz /home/glados/Data_Storage_2/{0}/mc/*/pass0/*.dst".format(det,date))
    # Apply DCR to ps5 file
    os.system("dcrmdpcgf.run -f -fd {0} /home/glados/Data_Storage_2/{1}/mc/processing/concat_dst/pr.1e17-1e21.20190625-{2}.{1}.ps5.dst.gz -o /home/glados/Data_Storage_2/{1}/mc/processing/concat_drmdpcgf_dst/".format(fd,det,date))
    # Convert the DST files into ROOT trees
    os.system("dst2rt_ta.run -f /home/glados/Data_Storage_2/{0}/mc/processing/concat_drmdpcgf_dst/pr.1e17-1e21.20190625-{1}.{0}.ps5.drmdpcgf.dst.gz -o /home/glados/Data_Storage_2/{0}/mc/processing/root_all_files/".format(det,date))
    os.system("dst2rt_ta.run -f /home/glados/Data_Storage_2/{0}/mc/processing/concat_dst/pr.1e17-1e21.20190625-{1}.{0}.ps0.dst.gz -o /home/glados/Data_Storage_2/{0}/mc/processing/root_all_files/".format(det,date))

def find_last_date(det):
    file = os.popen('ls -drt /home/glados/Data_Storage_2/{0}/mc/[0-9][0-9]* | sort -n | tail -n1'.format(det)).read()
    date = file[-9:-1]
    return date

missing_E_corr_def = '''
Double_t missing_E_corr(Double_t eng_EeV)
{
  Double_t eng_eV = 1e18 * eng_EeV;

  TF1 *JHK = new TF1("JHK","x / (-16.62 + 3.58 * log10(x) - 0.2786 * TMath::Power(log10(x),2) + 0.009775 * TMath::Power(log10(x),3) - 0.0001299 * TMath::Power(log10(x),4))",18.0,21.0);

  Double_t E = JHK->Eval(eng_eV);

  return E; 
}
'''
ROOT.gInterpreter.Declare(missing_E_corr_def)

spectrum_fun_def = '''
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
'''
ROOT.gInterpreter.Declare(spectrum_fun_def)

weight_def = '''
Double_t weight(Double_t log10en)
{
  return TMath::Power(10.0, 2*(log10en-log10en_min))*spectrum_function(log10en);
}
'''
ROOT.gInterpreter.Declare(weight_def)

if __name__=='__main__':
    # Define variables
    concat,det = read_in_args()
    last_date = find_last_date(det)
    
    # Concatinate files?
    if concat is True:
        concat_dst(det,last_date)

    # Load Dmitri's DST libs
    ROOT.gSystem.Load("libTDSTio.so")

    # Create Canvas
    canvas = ROOT.TCanvas("c1")
    canvas.cd()
    
    # Load mc pass5 root tree
    filepath = "/media/glados/Data_Storage_2/{0}/mc/processing/root_all_files/pr.1e17-1e21.20190625-{1}.{0}.ps5.drmdpcgf.root".format(det,last_date)
    file = ROOT.TFile(filepath,"r")
    taTree = file.Get("taTree")
    
    # Make TH1
    hErecon  = ROOT.TH1D("hErecon",";log_{10}(E/eV);N_{EVENTS} / BIN",40,17,21)
    
    
    # Draw Recon and save energy histogram
    taTree.Draw("log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>>hErecon","(etrack->qualct==1) * weight(log10(mc04->energy))","sames")
    hErecon.SaveAs("/home/glados/Data_Storage_2/{0}/mc/processing/hErecon_E3.root".format(det))

    # Prompt exit
    #print('\nPress enter to continue.\n')
    #os.system('read gonext')
    #print('Plotting thrown...')

    # Load mc pass0 root tree
    filepath = "/media/glados/Data_Storage_2/{0}/mc/processing/root_all_files/pr.1e17-1e21.20190625-{1}.{0}.ps0.root".format(det,last_date)
    file = ROOT.TFile(filepath,"r")
    taTree2 = file.Get("taTree")

    # Draw thrown and save energy histogram
    hEthrown = ROOT.TH1D("hEthrown",";log_{10}(E/eV);N_{EVENTS} / BIN",40,17,21)
    taTree2.Draw("log10(mc04->energy)>>hEthrown", "weight(log10(mc04->energy))")
    hEthrown.SaveAs("/home/glados/Data_Storage_2/{0}/mc/processing/hEthrown_E3.root".format(det))

    # Prompt exit
    print('\nPress enter to continue.\n')
    os.system('read gonext')
