#!/usr/bin/env python3

# Import Libs
import ROOT
import TLPY
import TDSTio
import os
import argparse
from sys import exit

################################### Quality Cuts ###################################
# Call in external C script with quality cut functions
ROOT.gROOT.LoadMacro("qualct.C")
        
#upward_events         = "stps2.upward[6] == 0"
#rayleigh_filter       = "stps2.plog[6] > 2.0"
#theta_cut             = "hctim.mthe[6] <= (65 * (TMath::Pi() / 180.0))"
#rms_dev_offpln        = "stpln.rmsdevpln[6] < 1.2"
#tracklength_cut       = "stpln.tracklength[6] > 8.0"
#inverse_angular_speed = "(stpln.crossingtime[6]/stpln.tracklength[6]) > 0.05"
#ph_per_gtube_cut      = "stpln.ph_per_gtube[6] > 50"
#good_tubes_per_degree = "(stpln.eye_ngtube[6]/stpln.tracklength[6]) < 5 && (stpln.eye_ngtube[6]/stpln.tracklength[6]) > 0.75"
#psi_cut               = "hctim.mpsi[6] < (120 * (TMath::Pi() / 180.0))"
#first_interaction_cut = "prfc.dep[6][0] < 1000"
#geom_chi2_dof_cut     = "(hctim.mchi2[6]/(stpln.eye_ngtube[6]-3)) <=5"
#prof_chi2_dof_cut     = "(prfc.chi2[6]/prfc.ndf[6]) <= 5"

qtcuts = "qualct==1"
###################################################################################

def read_in_args():
    parser = argparse.ArgumentParser(description = 'Concatinate data files (optional) and plot the data energy histogram.')
    parser.add_argument('-det',metavar='detector',action='store',help='Detector that you want to produce a data histogram for. (brtax4/mdtax4)',required=True)
    parser.add_argument('-plot',metavar='plot',action='store',help='Choose which resolution plots you want to plot. (all,energy,psi,theta,rp)',default='all')
    args = parser.parse_args()
    det=args.det
    if det == 'mdtax4':
        det_id = 6
    if det == 'brtax4':
        det_id = 7
    plot=args.plot
    print('Plotting ' + plot + ' resolution plot(s). This will take a couple minutes.')
    return det,plot,det_id

def plot_resolution(taTree,plot,det,det_id):
    kRed = 2
    if plot == 'energy' or plot == 'rp':
        hist_range     = 1
        neg_hist_range = -1
        if plot == 'energy':
            hist_title  = 'ln(E_{RECON}/E_{THROWN});N_{EVENTS}/bin'
            prof_titles = ';log_{10}(E_THROWN}/eV);ln(E_{RECON}/E_{THROWN})'
            hist_draw   = '(log((missing_E_corr(prfc->eng[{0}]))/mc04->energy))'.format(det_id)
        if plot == 'rp':
            hist_title  = 'ln(Rp_{RECON}/Rp_{THROWN})'
            prof_titles = ';log_{10}(E_{RECON}/eV);ln(Rp_{RECON} - Rp_{THROWN})'
            hist_draw   = 'log(hctim->mrp[{0}]/mc04->rp[{0}])'.format(det_id)
    if plot == 'psi' or plot == 'theta':
        hist_range     = 40
        neg_hist_range = -40
        if plot == 'psi':
            hist_title  = '#psi_{RECON} - #psi_{THROWN};N_{EVENTS}/bin'
            prof_titles = ';log_{10}(E_{RECON}/eV);#psi_{RECON} - #psi_{THROWN}'
            hist_draw   = '57.296*hctim->mpsi[{0}]-mc04->psi[{0}]*57.296'.format(det_id)
        if plot == 'theta':
            hist_title  = '#theta_{RECON}-#theta_{THROWN};N_{EVENTS}/bin'
            prof_titles = ';log_{10}(E_{THROWN}/eV);#theta_{RECON} - #theta_{THROWN}'
            hist_draw   = '57.296*hctim->mthe[{0}]-mc04->theta*57.296'.format(det_id)
        
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
    taTree.Draw("{0}>>h{1}18".format(hist_draw,plot),"{0} && log10((missing_E_corr(prfc->eng[{1}])))>18.0 && log10((missing_E_corr(prfc->eng[{1}])))<=18.5".format(qtcuts,det_id))
    h18.Fit("gaus");h18.GetFunction("gaus").SetLineColor(kRed);h18.GetFunction("gaus").SetLineWidth(3)
    c18.Update()
    c18.SaveAs("/media/glados/Data_Storage_2/plots/resolution/hybrid/{0}/{1}_resolution_hist_18.png".format(det,plot))

    c185.cd()
    taTree.Draw("{0}>>h{1}185".format(hist_draw,plot),"{0} && log10((missing_E_corr(prfc->eng[{1}])))>18.5 && log10((missing_E_corr(prfc->eng[{1}])))<=19.0".format(qtcuts,det_id))
    h185.Fit("gaus");h185.GetFunction("gaus").SetLineColor(kRed);h185.GetFunction("gaus").SetLineWidth(3)
    c185.Update()
    c185.SaveAs("/media/glados/Data_Storage_2/plots/resolution/hybrid/{0}/{1}_resolution_hist_185.png".format(det,plot))

    c19.cd()
    taTree.Draw("{0}>>h{1}19".format(hist_draw,plot),"{0} && log10((missing_E_corr(prfc->eng[{1}])))>19 && log10((missing_E_corr(prfc->eng[{1}])))<=19.5".format(qtcuts,det_id))
    h19.Fit("gaus");h19.GetFunction("gaus").SetLineColor(kRed);h19.GetFunction("gaus").SetLineWidth(3)
    c19.Update()
    c19.SaveAs("/media/glados/Data_Storage_2/plots/resolution/hybrid/{0}/{1}_resolution_hist_19.png".format(det,plot))

    c195.cd()
    taTree.Draw("{0}>>h{1}195".format(hist_draw,plot),"{0} && log10((missing_E_corr(prfc->eng[{1}])))>19.5 && log10((missing_E_corr(prfc->eng[{1}])))<=20".format(qtcuts,det_id))
    h195.Fit("gaus");h195.GetFunction("gaus").SetLineColor(kRed);h195.GetFunction("gaus").SetLineWidth(3)
    c195.Update()
    c195.SaveAs("/media/glados/Data_Storage_2/plots/resolution/hybrid/{0}/{1}_resolution_hist_195.png".format(det,plot))

    c20.cd()
    taTree.Draw("{0}>>h{1}20".format(hist_draw,plot),"{0} && log10((missing_E_corr(prfc->eng[{1}])))>20 && log10((missing_E_corr(prfc->eng[{1}])))<=20.5".format(qtcuts,det_id))
    h20.Fit("gaus");h20.GetFunction("gaus").SetLineColor(kRed);h20.GetFunction("gaus").SetLineWidth(3)
    c20.Update()
    c20.SaveAs("/media/glados/Data_Storage_2/plots/resolution/hybrid/{0}/{1}_resolution_hist_20.png".format(det,plot))

def prof_declare(det_id):
    energy_def = '''
    void energy(TTree *taTree)
    {
    TCanvas *ceprof = new TCanvas("ceprof","ceprof",600,600);
    TProfile *profe = new TProfile("profe",";log_{10}(E_THROWN}/eV);ln(E_{RECON}/E_{THROWN})",30,18.0,20.5);
    TH2F *hebox = new TH2F("hebox",";log_{10}(E_{thrown}/eV);ln(E_{RECON}/E_{THROWN})",50,17,21,100,-1,1);
    ceprof->cd();
    taTree->Draw("(log((missing_E_corr(prfc->eng[6]))/mc04->energy)):log10(mc04->energy)>>hebox","%s","colz");
    taTree->Draw("(log((missing_E_corr(prfc->eng[6]))/mc04->energy)):log10(mc04->energy)>>profe","%s","prof same");
    
    profe->SetLineColor(kRed);
    profe->SetLineWidth(2);
    ceprof->Update();
    }
    ''' % (qtcuts,qtcuts)
    ROOT.gInterpreter.Declare(energy_def)
    
    psi_def = '''
    void psi(TTree *taTree)
    {
    TCanvas *cpsiprof = new TCanvas("cpsiprof","cpsiprof",600,600);
    TProfile *profpsi = new TProfile("profpsi",";log_{10}(E_{RECON}/eV);#psi_{RECON} - #psi_{THROWN}",40,18.0,20.5);
    TH2F *hpsibox = new TH2F("hpsibox",";log_{10}(E_{RECON}/eV);#psi_{RECON} - #psi_{THROWN}",50,17,21,100,-40,40);
    cpsiprof->cd();
    taTree->Draw("57.296*hctim->mpsi[6]-mc04->psi[%d]*57.296:log10((missing_E_corr(prfc->eng[6])))>>hpsibox","%s","colz");
    taTree->Draw("57.296*hctim->mpsi[6]-mc04->psi[%d]*57.296:log10((missing_E_corr(prfc->eng[6])))>>profpsi","%s","prof sames");
    
    profpsi->SetLineColor(kRed);
    profpsi->SetLineWidth(2);
    cpsiprof->Update();
    }
    ''' % (det_id,qtcuts,det_id,qtcuts) 
    ROOT.gInterpreter.Declare(psi_def)
    
    theta_def = '''
    void theta(TTree *taTree)
    {
    TCanvas *ctheprof = new TCanvas("ctheprof","ctheprof",600,600);
    TProfile *profThe = new TProfile("profThe",";log_{10}(E_{THROWN}/eV);#theta_{RECON} - #theta_{THROWN}",40,18.0,20.5);
    TH2F *hThebox = new TH2F("hThebox",";log_{10}(E_{THROWN}/eV);#theta_{RECON} - #theta_{THROWN}",50,17,21,100,-40,40);
    ctheprof->cd();
    taTree->Draw("57.296*hctim->mthe[6]-mc04->theta*57.296:log10((missing_E_corr(prfc->eng[6])))>>hThebox","%s","colz");
    taTree->Draw("57.296*hctim->mthe[6]-mc04->theta*57.296:log10((missing_E_corr(prfc->eng[6])))>>profThe","%s","prof same");
    
    profThe->SetLineColor(kRed);
    profThe->SetLineWidth(2);
    ctheprof->Update();
    }
    ''' % (qtcuts,qtcuts)
    ROOT.gInterpreter.Declare(theta_def)

    rp_df = '''
    void rp(TTree *taTree)
    {
    TCanvas *crpprof = new TCanvas("crpprof","crpprof",600,600);
    TProfile *profRp = new TProfile("profRp",";log_{10}(E_{RECON}/eV);ln(Rp_{RECON} - Rp_{THROWN})",40,18.0,20.5);
    TH2F *hRpbox = new TH2F("hRpbox",";log_{10}(E_{RECON}/eV);ln(Rp_{RECON} / Rp_{THROWN})",50,17,21,40,-1,1);
    crpprof->cd();
    taTree->Draw("log(hctim->mrp[6]/mc04->rp[%d]):log10((missing_E_corr(prfc->eng[6])))>>hRpbox","%s","colz");
    taTree->Draw("log(hctim->mrp[6]/mc04->rp[%d]):log10((missing_E_corr(prfc->eng[6])))>>profRp","%s","prof same");
    
    profRp->SetLineColor(kRed);
    profRp->SetLineWidth(2);
    crpprof->Update();
    }
    ''' % (det_id,qtcuts,det_id,qtcuts)
    ROOT.gInterpreter.Declare(rp_df)

def find_last_date(det):
    file = os.popen('ls -drt {0}/[0-9][0-9]* | sort -n | tail -n1'.format(os.environ['MDTAX4_MC_ROOT'])).read()
    date = file[-9:-1]
    return date


if __name__=='__main__':
    # Define variables
    det,plot,det_id = read_in_args()
    last_date = find_last_date(det)
    if det == 'mdtax4':
        filename = '{0}/mc/hybrid_mc.all.cut.root'.format(os.environ['MDTAX4_HYBRID_ROOT'])
    elif det == 'brtax4':
        filename = '{0}/mc/hybrid_mc.all.cut.root'.format(os.environ['BRTAX4_HYBRID_ROOT'])
    else:
        exit(1)
        
    # Declare C prof plot functions
    prof_declare(det_id)

    # Open MC file
    file = ROOT.TFile(filename,"r")
    taTree = file.Get("taTree")

    # Choose which resolution plots to plot
    if plot == 'all':
        plot_resolution(taTree,'energy',det,det_id)
        ROOT.energy(taTree)
        print('\nPress enter to continue.\n')
        os.system('read gonext')
        plot_resolution(taTree,'psi',det,det_id)
        ROOT.psi(taTree)
        print('\nPress enter to continue.\n')
        os.system('read gonext')
        plot_resolution(taTree,'theta',det,det_id)
        ROOT.theta(taTree)
        print('\nPress enter to continue.\n')
        os.system('read gonext')
        plot_resolution(taTree,'rp',det,det_id)
        ROOT.rp(taTree)
        print('\nPress enter to continue.\n')
        os.system('read gonext')
    if plot == 'energy':
        plot_resolution(taTree,plot,det,det_id)
        ROOT.energy(taTree)
        print('\nPress enter to continue.\n')
        os.system('read gonext')
    if plot == 'psi':
        plot_resolution(taTree,plot,det,det_id)
        ROOT.psi(taTree)
        print('\nPress enter to continue.\n')
        os.system('read gonext')
    if plot == 'theta':
        plot_resolution(taTree,plot,det,det_id)
        ROOT.theta(taTree)
        print('\nPress enter to continue.\n')
        os.system('read gonext')
    if plot == 'rp':
        plot_resolution(taTree,plot,det,det_id)
        ROOT.rp(taTree)
        print('\nPress enter to continue.\n')
        os.system('read gonext')
    
