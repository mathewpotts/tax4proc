#!/usr/bin/env python3

# Import Libs
import ROOT
import TDSTio
import TLPY
import os
import argparse

def read_in_args():
    parser = argparse.ArgumentParser(description = 'Concatinate data files (optional) and plot the data energy histogram.')
    parser.add_argument('-det',metavar='detector',action='store',help='Detector that you want to produce a data histogram for. (brtax4/mdtax4)',required=True)
    parser.add_argument('-plot',metavar='plot',action='store',help='Choose which resolution plots you want to plot. (all,energy,psi,theta,rp)',default='all')
    args = parser.parse_args()
    det=args.det
    if det == 'mdtax4':
        det_id  = 6
        fd_path = os.environ['MDTAX4_DATA_ROOT']
    if det == 'brtax4':
        det_id  = 7
        fd_path = os.environ['BRTAX4_DATA_ROOT']
    os.system('mkdir -p {0}/plots/resolution'.format(fd_path)
    plot=args.plot
    print('Plotting ' + plot + ' resolution plot(s). This will take a couple minutes.')
    return det,plot,det_id,fd_path

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

def plot_resolution(taTree,plot,fd_path,det_id):
    kRed = 2
    if plot == 'energy' or plot == 'rp':
        hist_range     = 1
        neg_hist_range = -1
        if plot == 'energy':
            hist_title  = 'ln(E_{RECON}/E_{THROWN});N_{EVENTS}/bin'
            prof_titles = ';log_{10}(E_THROWN}/eV);ln(E_{RECON}/E_{THROWN})'
            hist_draw   = '(log((missing_E_corr(prfc->eng[etrack->udata[0]]))/mc04->energy))'
        if plot == 'rp':
            hist_title  = 'ln(Rp_{RECON}/Rp_{THROWN})'
            prof_titles = ';log_{10}(E_{RECON}/eV);ln(Rp_{RECON} - Rp_{THROWN})'
            hist_draw   = 'log(hctim->mrp[int(etrack->udata[0])]/mc04->rp[{0}])'.format(det_id)
    if plot == 'psi' or plot == 'theta':
        hist_range     = 40
        neg_hist_range = -40
        if plot == 'psi':
            hist_title  = '#psi_{RECON} - #psi_{THROWN};N_{EVENTS}/bin'
            prof_titles = ';log_{10}(E_{RECON}/eV);#psi_{RECON} - #psi_{THROWN}'
            hist_draw   = '57.296*hctim->mpsi[int(etrack->udata[0])]-mc04->psi[{0}]*57.296'.format(det_id)
        if plot == 'theta':
            hist_title  = '#theta_{RECON}-#theta_{THROWN};N_{EVENTS}/bin'
            prof_titles = ';log_{10}(E_{THROWN}/eV);#theta_{RECON} - #theta_{THROWN}'
            hist_draw   = '57.296*hctim->mthe[int(etrack->udata[0])]-mc04->theta*57.296'
        
    c18  = ROOT.TCanvas("c{0}18".format(plot),"c{0}18".format(plot),600,600)
    c185 = ROOT.TCanvas("c{0}185".format(plot),"c{0}185".format(plot),600,600)
    c19  = ROOT.TCanvas("c{0}19".format(plot),"c{0}19".format(plot),600,600)
    c195 = ROOT.TCanvas("c{0}195".format(plot),"c{0}195".format(plot),600,600)
    c20  = ROOT.TCanvas("c{0}20".format(plot),"c{0}20".format(plot),600,600)

    h18  = ROOT.TH1F("h{0}18".format(plot),"18.0 < log_{10}(E_{RECON}/eV) <= 18.5;%s}"%(hist_title),100,neg_hist_range,hist_range)
    h185 = ROOT.TH1F("h{0}185".format(plot),"18.5 < log_{10}(E_{RECON}/eV) <= 19.0;%s"%(hist_title),100,neg_hist_range,hist_range)
    h19  = ROOT.TH1F("h{0}19".format(plot),"19.0 < log_{10}(E_{RECON}/eV) <= 19.5;%s"%(hist_title),100,neg_hist_range,hist_range)
    h195 = ROOT.TH1F("h{0}195".format(plot),"19.5 < log_{10}(E_{RECON}/eV) <= 20.0;%s"%(hist_title),100,neg_hist_range,hist_range)
    h20  = ROOT.TH1F("h{0}20".format(plot),"20.0 < log_{10}(E_{RECON}/eV) <= 20.5;%s"%(hist_title),100,neg_hist_range,hist_range)

    ROOT.gStyle.SetOptFit(1)
    
    c18.cd()
    taTree.Draw("{0}>>h{1}18".format(hist_draw,plot),"etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.0 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=18.5")
    h18.Fit("gaus");h18.GetFunction("gaus").SetLineColor(kRed);h18.GetFunction("gaus").SetLineWidth(3)
    c18.Update()
    c18.SaveAs("{0}/plots/resolution/{1}_resolution_hist_18.png".format(fd_path,plot))

    c185.cd()
    taTree.Draw("{0}>>h{1}185".format(hist_draw,plot),"etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.0")
    h185.Fit("gaus");h185.GetFunction("gaus").SetLineColor(kRed);h185.GetFunction("gaus").SetLineWidth(3)
    c185.Update()
    c185.SaveAs("{0}/plots/resolution/{1}_resolution_hist_185.png".format(fd_path,plot))

    c19.cd()
    taTree.Draw("{0}>>h{1}19".format(hist_draw,plot),"etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.5")
    h19.Fit("gaus");h19.GetFunction("gaus").SetLineColor(kRed);h19.GetFunction("gaus").SetLineWidth(3)
    c19.Update()
    c19.SaveAs("{0}/plots/resolution/{1}_resolution_hist_19.png".format(fd_path,plot))

    c195.cd()
    taTree.Draw("{0}>>h{1}195".format(hist_draw,plot),"etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.5 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20")
    h195.Fit("gaus");h195.GetFunction("gaus").SetLineColor(kRed);h195.GetFunction("gaus").SetLineWidth(3)
    c195.Update()
    c195.SaveAs("{0}/plots/resolution/{1}_resolution_hist_195.png".format(fd_path,plot))

    c20.cd()
    taTree.Draw("{0}>>h{1}20".format(hist_draw,plot),"etrack->qualct==1 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>20 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.5")
    h20.Fit("gaus");h20.GetFunction("gaus").SetLineColor(kRed);h20.GetFunction("gaus").SetLineWidth(3)
    c20.Update()
    c20.SaveAs("{0}/plots/resolution/{1}_resolution_hist_20.png".format(fd_path,plot))

def prof_declare(det_id):
    energy_def = '''
    void energy(TTree *taTree)
    {
    TCanvas *ceprof = new TCanvas("ceprof","ceprof",600,600);
    TProfile *profe = new TProfile("profe",";log_{10}(E_THROWN}/eV);ln(E_{RECON}/E_{THROWN})",30,18.0,20.5);
    TH2F *hebox = new TH2F("hebox",";log_{10}(E_{thrown}/eV);ln(E_{RECON}/E_{THROWN})",50,17,21,100,-1,1);
    ceprof->cd();
    taTree->Draw("(log((missing_E_corr(prfc->eng[etrack->udata[0]]))/mc04->energy)):log10(mc04->energy)>>hebox","etrack->qualct==1","colz");
    taTree->Draw("(log((missing_E_corr(prfc->eng[etrack->udata[0]]))/mc04->energy)):log10(mc04->energy)>>profe","etrack->qualct==1","prof same");
    
    profe->SetLineColor(kRed);
    profe->SetLineWidth(2);
    ceprof->Update();
    }
    '''
    ROOT.gInterpreter.Declare(energy_def)
    
    psi_def = '''
    void psi(TTree *taTree)
    {
    TCanvas *cpsiprof = new TCanvas("cpsiprof","cpsiprof",600,600);
    TProfile *profpsi = new TProfile("profpsi",";log_{10}(E_{RECON}/eV);#psi_{RECON} - #psi_{THROWN}",40,18.0,20.5);
    TH2F *hpsibox = new TH2F("hpsibox",";log_{10}(E_{RECON}/eV);#psi_{RECON} - #psi_{THROWN}",50,17,21,100,-40,40);
    cpsiprof->cd();
    taTree->Draw("57.296*hctim->mpsi[int(etrack->udata[0])]-mc04->psi[%d]*57.296:log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>>hpsibox","etrack->qualct==1","colz");
    taTree->Draw("57.296*hctim->mpsi[int(etrack->udata[0])]-mc04->psi[%d]*57.296:log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>>profpsi","etrack->qualct==1","profs sames");
    
    profpsi->SetLineColor(kRed);
    profpsi->SetLineWidth(2);
    cpsiprof->Update();
    }
    ''' % (det_id,det_id) 
    ROOT.gInterpreter.Declare(psi_def)
    
    theta_def = '''
    void theta(TTree *taTree)
    {
    TCanvas *ctheprof = new TCanvas("ctheprof","ctheprof",600,600);
    TProfile *profThe = new TProfile("profThe",";log_{10}(E_{THROWN}/eV);#theta_{RECON} - #theta_{THROWN}",40,18.0,20.5);
    TH2F *hThebox = new TH2F("hThebox",";log_{10}(E_{THROWN}/eV);#theta_{RECON} - #theta_{THROWN}",50,17,21,100,-40,40);
    ctheprof->cd();
    taTree->Draw("57.296*hctim->mthe[int(etrack->udata[0])]-mc04->theta*57.296:log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>>hThebox","etrack->qualct==1","colz");
    taTree->Draw("57.296*hctim->mthe[int(etrack->udata[0])]-mc04->theta*57.296:log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>>profThe","etrack->qualct==1","profs same");
    
    profThe->SetLineColor(kRed);
    profThe->SetLineWidth(2);
    ctheprof->Update();
    }
    '''
    ROOT.gInterpreter.Declare(theta_def)

    rp_df = '''
    void rp(TTree *taTree)
    {
    TCanvas *crpprof = new TCanvas("crpprof","crpprof",600,600);
    TProfile *profRp = new TProfile("profRp",";log_{10}(E_{RECON}/eV);ln(Rp_{RECON} - Rp_{THROWN})",40,18.0,20.5);
    TH2F *hRpbox = new TH2F("hRpbox",";log_{10}(E_{RECON}/eV);ln(Rp_{RECON} / Rp_{THROWN})",50,17,21,40,-1,1);
    crpprof->cd();
    taTree->Draw("log(hctim->mrp[int(etrack->udata[0])]/mc04->rp[%d]):log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>>hRpbox","etrack->qualct==1","colz");
    taTree->Draw("log(hctim->mrp[int(etrack->udata[0])]/mc04->rp[%d]):log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>>profRp","etrack->qualct==1","prof same");
    
    profRp->SetLineColor(kRed);
    profRp->SetLineWidth(2);
    crpprof->Update();
    }
    ''' % (det_id,det_id)
    ROOT.gInterpreter.Declare(rp_df)

def find_last_date(fd_path):
    file = os.popen('ls -drt {0}/mc/[0-9][0-9]* | sort -n | tail -n1'.format(fd_path)).read()
    date = file[-9:-1]
    return date


if __name__=='__main__':
    # Define variables
    det,plot,det_id,fd_path = read_in_args()
    last_date = find_last_date(fd_path)
    filename = '{0}/mc/processing/root_all_files/pr.1e17-1e21.20190625-{1}.{0}.ps5.drmdpcgf.root'.format(fd_path,last_date)

    # Declare C prof plot functions
    prof_declare(det_id)

    # Open MC file
    file = ROOT.TFile(filename,"r")
    taTree = file.Get("taTree")

    # Choose which resolution plots to plot
    if plot == 'all':
        plot_resolution(taTree,'energy',fd_path,det_id)
        ROOT.energy(taTree)
        print('\nPress enter to continue.\n')
        os.system('read gonext')
        plot_resolution(taTree,'psi',fd_path,det_id)
        ROOT.psi(taTree)
        print('\nPress enter to continue.\n')
        os.system('read gonext')
        plot_resolution(taTree,'theta',fd_path,det_id)
        ROOT.theta(taTree)
        print('\nPress enter to continue.\n')
        os.system('read gonext')
        plot_resolution(taTree,'rp',fd_path,det_id)
        ROOT.rp(taTree)
        print('\nPress enter to continue.\n')
        os.system('read gonext')
    if plot == 'energy':
        plot_resolution(taTree,plot,fd_path,det_id)
        ROOT.energy(taTree)
        print('\nPress enter to continue.\n')
        os.system('read gonext')
    if plot == 'psi':
        plot_resolution(taTree,plot,fd_path,det_id)
        ROOT.psi(taTree)
        print('\nPress enter to continue.\n')
        os.system('read gonext')
    if plot == 'theta':
        plot_resolution(taTree,plot,fd_path,det_id)
        ROOT.theta(taTree)
        print('\nPress enter to continue.\n')
        os.system('read gonext')
    if plot == 'rp':
        plot_resolution(taTree,plot,fd_path,det_id)
        ROOT.rp(taTree)
        print('\nPress enter to continue.\n')
        os.system('read gonext')
    
