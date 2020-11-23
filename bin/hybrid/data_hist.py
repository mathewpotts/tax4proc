#!/usr/bin/env python2

# Import Libs
import ROOT
import os
import argparse

def read_in_args():
    parser = argparse.ArgumentParser(description = 'Concatinate data files (optional) and plot the data energy histogram.')
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
    os.system("mkdir -p /home/glados/Data_Storage_2/{0}/pass5/data/processing/concat_dst/".format(det))
    os.system("mkdir -p /home/glados/Data_Storage_2/{0}/pass5/data/processing/concat_drmdpcgf_dst/".format(det))
    os.system("mkdir -p /home/glados/Data_Storage_2/{0}/pass5/data/processing/root_all_files/".format(det))
    
    # Concatinate files
    os.system("dstcat -o /home/glados/Data_Storage_2/{0}/pass5/data/processing/concat_dst/20190625-{1}.{0}.ps5.ps2g.dst.gz /home/glados/Data_Storage_2/{0}/pass5/data/*/*ps2g.dst.gz".format(det,date))
    #os.system("dstcat -o /home/glados/Data_Storage_2/{0}/pass5/data/processing/concat_dst/20190625-{1}.{0}.ps0.dst.gz /home/glados/Data_Storage_2/{0}/pass0/data/*/*fraw1.dst.gz".format(det,date))
    # Apply DCR to ps5 file
    os.system("dcrmdpcgf.run -f -fd {0} /home/glados/Data_Storage_2/{1}/pass5/data/processing/concat_dst/20190625-{2}.{1}.ps5.ps2g.dst.gz -o /home/glados/Data_Storage_2/{1}/pass5/data/processing/concat_drmdpcgf_dst/".format(fd,det,date))
    # Convert the DST files into ROOT trees
    os.system("dst2rt_ta.run -f /home/glados/Data_Storage_2/{0}/pass5/data/processing/concat_drmdpcgf_dst/20190625-{1}.{0}.ps5.ps2g.drmdpcgf.dst.gz -o /home/glados/Data_Storage_2/{0}/pass5/data/processing/root_all_files/".format(det,date))
    #os.system("dst2rt_ta.run -f /home/glados/Data_Storage_2/{0}/pass5/data/processing/concat_dst/20190625-{1}.{0}.ps0.dst.gz -o /home/glados/Data_Storage_2/{0}/pass5/data/processing/root_all_files/".format(det,date))

def find_last_date(det):
    file = os.popen('ls -drt /home/glados/Data_Storage_2/{0}/pass5/data/[0-9][0-9]* | sort -n | tail -n1'.format(det)).read()
    date = file[-9:-1]
    return date
    
gw_cut_def = '''
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
    if (o+t+s+w == 0){ // BRTAx4 weather codes
      return 1;
    }
  }
  return 0;
}
'''
ROOT.gInterpreter.Declare(gw_cut_def)

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
    
if __name__=='__main__':
    # Define variables
    nxbins        = 20
    log10en_lo    = 18.0
    log10en_up    = 20.0
    concat,det = read_in_args()
    last_date = find_last_date(det)
    
    # Concatinate files?
    if concat is True:
        concat_dst(det,last_date)
    
    # Load Dmitri's DST libs
    ROOT.gSystem.Load("libTDSTio.so")
    
    # Load pass5 root tree
    filepath = "/media/glados/Data_Storage_2/{0}/pass5/data/processing/root_all_files/20190625-{1}.{0}.ps5.ps2g.drmdpcgf.root".format(det,last_date)
    file = ROOT.TFile(filepath,"r")
    taTree = file.Get("taTree")
    
    # Make TH1
    hEdata = ROOT.TH1D("hEdata",";log_{10}(E/eV);N_{EVENTS} / BIN",nxbins,log10en_lo,log10en_up);
    
    # Draw and save energy histogram
    taTree.Draw("log10((missing_E_corr(prfc.eng[etrack.udata[0]])))>>hEdata","good_weather_cut(tlweat.code,fscn1->mir[1])==1 && etrack.qualct==1 && log10(missing_E_corr(prfc.eng[etrack.udata[0]]))>=18.5")
    hEdata.SaveAs("/media/glados/Data_Storage_2/{0}/pass5/data/processing/hEdata.root".format(det))
    
    # Prompt exit
    print('\nPress enter to continue.\n')
    os.system('read gonext')
