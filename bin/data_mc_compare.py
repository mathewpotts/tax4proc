#!/usr/bin/env python3

# Import Libs
import ROOT
import TDSTio
import TLPY
import os
import argparse
import numpy as np

# Call in external C script with quality cut functions
ROOT.gROOT.LoadMacro("hybrid/qualct.C")

def read_in_args():
    parser = argparse.ArgumentParser(description = 'Concatinate data files (optional) and plot the data energy histogram.')
    parser.add_argument('-det',metavar='detector',action='store',help='Detector that you want to produce a data histogram for. (brtax4/mdtax4)',required=True)
    parser.add_argument('-plot',metavar='plot',action='store',help='Choose which data-mc comparison plots you want to plot. (phi,theta,rp,psi,phpgt,track,ngmpe)',default='all')
    args = parser.parse_args()
    det=args.det
    if det == 'mdtax4':
        fd_path = os.environ['MDTAX4_DATA_ROOT']
    if det == 'brtax4':
        fd_path = os.environ['BRTAX4_DATA_ROOT']
    plot=args.plot
    print('Plotting ' + plot + ' data-mc comparison plot(s). This will take a couple minutes.')
    return det,plot,fd_path

def ngmpe(taTree1,taTree2,det):
    # Define constants
    nentriesdata = taTree1.GetEntries()
    nentriesmc   = taTree2.GetEntries()
    if det == 'mdtax4':
        nngmpebins = 4
        ngmpe_lo = 0.5
        ngmpe_up = 4
    elif det == 'brtax4':
        nngmpebins = 8
        ngmpe_lo = 0.5
        ngmpe_up = 8

    # Make data histogram
    hngmpedata    = ROOT.TH1D("hngmpedata","18.5 < log_{10}(E/eV);track (degrees);N_{EVENTS} / BIN",nngmpebins,ngmpe_lo,ngmpe_up);

    for i in range(0,nentriesdata):
        n_gmir = 0;
        taTree1.GetEntry(i);

        # Check if event exist
        try:
            taTree1.Draw("etrack->energy","good_weather_cut(tlweat.code,fscn1->mir[1])==1 && etrack->qualct==1 && (prfc->chi2[etrack->udata[0]]/prfc->ndf[etrack->udata[0]])<20 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5","goff",1,i)
        except:
            continue
      
        # Loop over all possible  mirror instances
        for j in taTree1.stpln.mir_ngtube:
            if j > 0 and j <= 256 and ROOT.good_weather_cut(taTree1.tlweat.code,taTree1.fscn1.mir[1])==1 and taTree1.etrack.qualct==1 and np.log10((ROOT.missing_E_corr(taTree1.prfc.eng[int(taTree1.etrack.udata[0])])))>18.5:
                n_gmir += 1

        hngmpedata.Fill(n_gmir)
        print("Data Progress: ",round(i/nentriesdata * 100,4),"%")

    #### MC ####
    # Make mc histogram
    hngmpemc    = ROOT.TH1D("hngmpemc","18.5 < log_{10}(E/eV);track (degrees);N_{EVENTS} / BIN",nngmpebins,ngmpe_lo,ngmpe_up);


    for i in range(0,nentriesmc):
        n_gmir = 0;
        taTree2.GetEntry(i);

        # Check if event exist
        try:
            taTree2.Draw("etrack->energy","etrack->qualct==1 && (prfc->chi2[etrack->udata[0]]/prfc->ndf[etrack->udata[0]])<20 && log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5","goff",1,i)
        except:
            continue
      
        # Loop over all possible  mirror instances
        for j in taTree2.stpln.mir_ngtube:
            if j > 0 and j <= 256 and taTree2.etrack.qualct==1 and np.log10((ROOT.missing_E_corr(taTree2.prfc.eng[int(taTree2.etrack.udata[0])])))>18.5:
                n_gmir += 1

        hngmpemc.Fill(n_gmir,ROOT.weight(np.log10(taTree2.mc04.energy)))
        print("MC Progress: ",round(i/nentriesmc * 100,4),"%")

    normal    =  hngmpedata.Integral()/hngmpemc.Integral()
    hngmpemc.Scale(normal)
    hngmpemc.Draw()
    ROOT.gPad.Update()
    sNgmpe = hngmpemc.FindObject("stats")

    #### Compare ####
    ROOT.gStyle.SetOptFit(1)
    kRed = 2
    kBlue = 4
    legend = ROOT.TLegend(0.1,0.7,0.3,0.9)
    legend.SetHeader("NGMPE Data-MC Comparison","C")
    mcentry = legend.AddEntry("hngmpemc","MC","l")
    mcentry.SetLineColor(kRed)
    dataentry = legend.AddEntry("hngmpedata","Data","p")
    dataentry.SetMarkerStyle(2)
    dataentry.SetMarkerColor(kBlue)
    ##########
    hngmpemc.SetLineColor(kRed)
    c2ngmpe = ROOT.TCanvas("c2ngmpe","c2ngmpe",800,600);
    rpngmpe = ROOT.TRatioPlot(hngmpedata,hngmpemc,"divsym")
    rpngmpe.SetH1DrawOpt("E0")
    rpngmpe.SetH2DrawOpt("HIST")
    c2ngmpe.cd()
    rpngmpe.Draw()
    rpngmpe.GetUpperRefYaxis().SetRangeUser(0,80)
    rpngmpe.GetLowerRefYaxis().SetRangeUser(0.1,3)
    rpngmpe.GetLowerRefYaxis().SetLabelSize(0.025)
    rpngmpe.GetLowerRefYaxis().SetTitle("ngmpe_{data}/ngmpe_{mcrecon}")
    rpngmpe.GetUpperPad().cd()
    sNgmpe.Draw()
    legend.Draw()
    rpngmpe.GetLowerRefGraph().Fit("pol1","W","",0,2.5)
    ROOT.gPad.Update()
    rpngmpe.SaveAs('{0}/{1}/plots/data-mc-compare/ngmpe.png'.format(os.environ['DATA_ROOT_DIR'],det))
    print('\nPress enter to continue.\n')
    os.system('read gonext')

def find_last_date(fd_path):
    file = os.popen('ls -drt {0}/pass5/data/[0-9][0-9]* | sort -n | tail -n1'.format(fd_path)).read()
    date_data = file[-9:-1]
    file = os.popen('ls -drt {0}/mc/[0-9][0-9]* | sort -n | tail -n1'.format(fd_path)).read()
    date_mc = file[-9:-1]
    return date_data,date_mc

def data_mc_compare(taTree1,taTree2,plot):
    # Define variables based on function input
    if plot == 'phi':
        nbins        = 20
        bin_lo       = 0
        bin_up       = 360
        hist_title   = '#phi (degrees);N_{EVENTS} / BIN'
        hist_draw    = '57.296*etrack->phi'
        compare_sym  = '#phi'
    if plot == 'theta':
        nbins        = 40
        bin_lo       = 0
        bin_up       = 70
        hist_title   = '#theta (degrees);N_{EVENTS} / BIN'
        hist_draw    = '57.296*hctim->mthe[int(etrack->udata[0])]'
        compare_sym  = '#theta'
    if plot == 'rp':
        nbins        = 50
        bin_lo       = 0
        bin_up       = 50
        hist_title   = 'Rp (km);N_{EVENTS} / BIN'
        hist_draw    = '1e-3*hctim->mrp[int(etrack->udata[0])]'
        compare_sym  = 'Rp'
    if plot == 'psi':
        nbins        = 20
        bin_lo       = 40
        bin_up       = 130
        hist_title   = '#psi (degrees);N_{EVENTS} / BIN'
        hist_draw    = '57.296*hctim->mpsi[int(etrack->udata[0])]'
        compare_sym  = '#psi'
    if plot == 'phpgt':
        nbins        = 50
        bin_lo       = 150
        bin_up       = 1400
        hist_title   = 'phpgt;N_{EVENTS} / BIN'
        hist_draw    = 'stpln->ph_per_gtube[6]'
        compare_sym  = 'phpgt'
    if plot == 'track':
        nbins        = 10
        bin_lo       = 8
        bin_up       = 28
        hist_title   = 'track (degrees);N_{EVENTS} / BIN'
        hist_draw    = 'stpln->tracklength[6]'
        compare_sym  = 'track'
    
    # Make data histograms
    hdata    = ROOT.TH1D("h{0}data".format(plot),"18.5 < log_{10}(E/eV);%s"%(hist_title),nbins,bin_lo,bin_up)
    #hdata18  = ROOT.TH1D("h{0}data18".format(plot),"18.0 < log_{10}(E/eV) <= 18.5;%s"%(hist_title),nbins,bin_lo,bin_up)
    hdata185 = ROOT.TH1D("h{0}data185".format(plot),"18.5 < log_{10}(E/eV) <= 19.0;%s"%(hist_title),nbins,bin_lo,bin_up)
    hdata19  = ROOT.TH1D("h{0}data19".format(plot),"19.0 < log_{10}(E/eV) <= 19.5;%s"%(hist_title),nbins,bin_lo,bin_up)
    hdata195 = ROOT.TH1D("h{0}data195".format(plot),"19.5 < log_{10}(E/eV) <= 20.0;%s"%(hist_title),nbins,bin_lo,bin_up)
    hdata20  = ROOT.TH1D("h{0}data20".format(plot),"20.0 < log_{10}(E/eV) <= 20.5;%s"%(hist_title),nbins,bin_lo,bin_up)

    # Draw data histograms
    taTree1.Draw("{0}>>h{1}data".format(hist_draw,plot),"good_weather_cut(tlweat.code,fscn1->mir[1])==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5","goff")
    #taTree1.Draw("{0}>>h{1}data18".format(hist_draw,plot),"good_weather_cut(tlweat.code,fscn1->mir[1])==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=18.5 && (prfc->chi2[etrack->udata[0]]/prfc->ndf[etrack->udata[0]])<20","goff")
    taTree1.Draw("{0}>>h{1}data185".format(hist_draw,plot),"good_weather_cut(tlweat.code,fscn1->mir[1])==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.0 && (prfc->chi2[etrack->udata[0]]/prfc->ndf[etrack->udata[0]])<20","goff")
    taTree1.Draw("{0}>>h{1}data19".format(hist_draw,plot),"good_weather_cut(tlweat.code,fscn1->mir[1])==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.5 && (prfc->chi2[etrack->udata[0]]/prfc->ndf[etrack->udata[0]])<20","goff")
    taTree1.Draw("{0}>>h{1}data195".format(hist_draw,plot),"good_weather_cut(tlweat.code,fscn1->mir[1])==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.0 && (prfc->chi2[etrack->udata[0]]/prfc->ndf[etrack->udata[0]])<20","goff")
    taTree1.Draw("{0}>>h{1}data20".format(hist_draw,plot),"good_weather_cut(tlweat.code,fscn1->mir[1])==1&&etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>20.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.5 && (prfc->chi2[etrack->udata[0]]/prfc->ndf[etrack->udata[0]])<20","goff")

    # Make MC histograms
    hmc    = ROOT.TH1D("h{0}mc".format(plot),"18.5 < log_{10}(E/eV);%s"%(hist_title),nbins,bin_lo,bin_up)
    #hmc18  = ROOT.TH1D("h{0}mc18".format(plot),"18.0 < log_{10}(E/eV) <= 18.5;%s"%(hist_title),nbins,bin_lo,bin_up)
    hmc185 = ROOT.TH1D("h{0}mc185".format(plot),"18.5 < log_{10}(E/eV) <= 19.0;%s"%(hist_title),nbins,bin_lo,bin_up)
    hmc19  = ROOT.TH1D("h{0}mc19".format(plot),"19.0 < log_{10}(E/eV) <= 19.5;%s"%(hist_title),nbins,bin_lo,bin_up)
    hmc195 = ROOT.TH1D("h{0}mc195".format(plot),"19.5 < log_{10}(E/eV) <= 20.0;%s"%(hist_title),nbins,bin_lo,bin_up)
    hmc20  = ROOT.TH1D("h{0}mc20".format(plot),"20.0 < log_{10}(E/eV) <= 20.5;%s"%(hist_title),nbins,bin_lo,bin_up)

    # Draw and scale MC histograms
    taTree2.Draw("{0}>>h{1}mc".format(hist_draw,plot),"weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5 && (prfc->chi2[etrack->udata[0]]/prfc->ndf[etrack->udata[0]])<20)")
    normal = hdata.Integral()/hmc.Integral()
    hmc.Scale(normal)
    ROOT.gPad.Update()
    s = hmc.FindObject("stats")
    ###########################
    #taTree2.Draw("{0}>>h{1}mc18".format(hist_draw,plot),"weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=18.5 && (prfc->chi2[etrack->udata[0]]/prfc->ndf[etrack->udata[0]])<20)")
    #normal = hdata18.Integral()/hmc18.Integral()
    #hmc18.Scale(normal)
    #ROOT.gPad.Update()
    #s18 = hmc18.FindObject("stats")
    ###########################
    taTree2.Draw("{0}>>h{1}mc185".format(hist_draw,plot),"weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>18.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.0 && (prfc->chi2[etrack->udata[0]]/prfc->ndf[etrack->udata[0]])<20)")
    normal = hdata185.Integral()/hmc185.Integral()
    hmc185.Scale(normal)
    ROOT.gPad.Update()
    s185 = hmc185.FindObject("stats")
    ###########################
    taTree2.Draw("{0}>>h{1}mc19".format(hist_draw,plot),"weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=19.5 && (prfc->chi2[etrack->udata[0]]/prfc->ndf[etrack->udata[0]])<20)")
    normal = hdata19.Integral()/hmc19.Integral()
    hmc19.Scale(normal)
    ROOT.gPad.Update()
    s19 = hmc19.FindObject("stats")
    ###########################
    taTree2.Draw("{0}>>h{1}mc195".format(hist_draw,plot),"weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>19.5&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20 && (prfc->chi2[etrack->udata[0]]/prfc->ndf[etrack->udata[0]])<20)")
    normal = hdata195.Integral()/hmc195.Integral()
    hmc195.Scale(normal)
    ROOT.gPad.Update()
    s195 = hmc195.FindObject("stats")
    ###########################
    taTree2.Draw("{0}>>h{1}mc20".format(hist_draw,plot),"weight(log10(mc04->energy))*(etrack->qualct==1&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>20.0&&log10((missing_E_corr(prfc->eng[etrack->udata[0]])))<=20.5 && (prfc->chi2[etrack->udata[0]]/prfc->ndf[etrack->udata[0]])<20)")
    normal = hdata20.Integral()/hmc20.Integral()
    hmc20.Scale(normal)
    ROOT.gPad.Update()
    s20 = hmc20.FindObject("stats")
    ###########################

    ### Compare ###
    ROOT.gStyle.SetOptFit(1)
    kRed = 2
    kBlue = 4
    
    # Set up legend
    legend = ROOT.TLegend(0.1,0.7,0.3,0.9)
    legend.SetHeader("{0} Data-MC Comparison".format(plot),"C")
    mcentry = legend.AddEntry("h{0}mc".format(plot),"MC","l")
    mcentry.SetLineColor(kRed)
    dataentry = legend.AddEntry("h{0}data".format(plot),"Data","p")
    dataentry.SetMarkerStyle(2)
    dataentry.SetMarkerColor(kBlue)

    hmc.SetLineColor(kRed)
    c2 = ROOT.TCanvas("c2{0}".format(plot),"c2{0}".format(plot),800,600)
    rp = ROOT.TRatioPlot(hdata,hmc,"divsym")
    rp.SetH1DrawOpt("E0")
    rp.SetH2DrawOpt("HIST")
    c2.cd()
    rp.Draw()
    rp.GetUpperRefYaxis().SetRangeUser(0,10)
    rp.GetLowerRefYaxis().SetRangeUser(0.1,3)
    rp.GetLowerRefYaxis().SetLabelSize(0.025)
    rp.GetLowerRefYaxis().SetTitle("%s_{data}/%s_{mcrecon}"%(compare_sym,compare_sym))
    rp.GetUpperPad().cd()
    s.Draw()
    legend.Draw()
    rp.GetLowerRefGraph().Fit("pol1","W")
    rp.SaveAs('{0}/{1}/data-mc-compare/{2}.png'.format(os.environ['DATA_ROOT_DIR'],det,plot))
    print('\nPress enter to continue.\n')
    os.system('read gonext')
    ###########################
    #hmc18.SetLineColor(kRed)
    #c218 = ROOT.TCanvas("c2{0}18".format(plot),"c2{0}18".format(plot),800,600)
    #rp18 = ROOT.TRatioPlot(hdata18,hmc18,"divsym")
    #rp18.SetH1DrawOpt("E0")
    #rp18.SetH2DrawOpt("HIST")
    #c218.cd()
    #rp18.Draw()
    #rp18.GetUpperRefYaxis().SetRangeUser(0,10)
    #rp18.GetLowerRefYaxis().SetRangeUser(0.1,3)
    #rp18.GetLowerRefYaxis().SetLabelSize(0.025)
    #rp18.GetLowerRefYaxis().SetTitle("%s_{data}/%s_{mcrecon}"%(compare_sym,compare_sym))
    #rp18.GetUpperPad().cd()
    #s18.Draw()
    #legend.Draw()
    #rp18.GetLowerRefGraph().Fit("pol1","W")
    #rp18.SaveAs('{0}/{1}/data-mc-compare/{2}_18-185.png'.format(os.environ['DATA_ROOT_DIR'],det,plot))
    #print('\nPress enter to continue.\n')
    #os.system('read gonext')
    ###########################
    hmc185.SetLineColor(kRed)
    c2185 = ROOT.TCanvas("c2{0}185".format(plot),"c2{0}185".format(plot),800,600)
    rp185 = ROOT.TRatioPlot(hdata185,hmc185,"divsym")
    rp185.SetH1DrawOpt("E0")
    rp185.SetH2DrawOpt("HIST")
    c2185.cd()
    rp185.Draw()
    rp185.GetUpperRefYaxis().SetRangeUser(0,10)
    rp185.GetLowerRefYaxis().SetRangeUser(0.1,3)
    rp185.GetLowerRefYaxis().SetLabelSize(0.025)
    rp185.GetLowerRefYaxis().SetTitle("%s_{data}/%s_{mcrecon}"%(compare_sym,compare_sym))
    rp185.GetUpperPad().cd()
    s185.Draw()
    legend.Draw()
    rp185.GetLowerRefGraph().Fit("pol1","W")
    rp185.SaveAs('{0}/{1}/data-mc-compare/{2}_185-19.png'.format(os.environ['DATA_ROOT_DIR'],det,plot))
    print('\nPress enter to continue.\n')
    os.system('read gonext')
    ###########################
    hmc19.SetLineColor(kRed)
    c219 = ROOT.TCanvas("c2{0}19".format(plot),"c2{0}19".format(plot),800,600)
    rp19 = ROOT.TRatioPlot(hdata19,hmc19,"divsym")
    rp19.SetH1DrawOpt("E0")
    rp19.SetH2DrawOpt("HIST")
    c219.cd()
    rp19.Draw()
    rp19.GetUpperRefYaxis().SetRangeUser(0,10)
    rp19.GetLowerRefYaxis().SetRangeUser(0.1,3)
    rp19.GetLowerRefYaxis().SetLabelSize(0.025)
    rp19.GetLowerRefYaxis().SetTitle("%s_{data}/%s_{mcrecon}"%(compare_sym,compare_sym))
    rp19.GetUpperPad().cd()
    s19.Draw()
    legend.Draw()
    rp19.GetLowerRefGraph().Fit("pol1","W")
    rp19.SaveAs('{0}/{1}/data-mc-compare/{2}_19-195.png'.format(os.environ['DATA_ROOT_DIR'],det,plot))
    print('\nPress enter to continue.\n')
    os.system('read gonext')
    ###########################
    hmc195.SetLineColor(kRed)
    c2195 = ROOT.TCanvas("c2{0}195".format(plot),"c2{0}195".format(plot),800,600)
    rp195 = ROOT.TRatioPlot(hdata195,hmc195,"divsym")
    rp195.SetH1DrawOpt("E0")
    rp195.SetH2DrawOpt("HIST")
    c2195.cd()
    rp195.Draw()
    rp195.GetUpperRefYaxis().SetRangeUser(0,10)
    rp195.GetLowerRefYaxis().SetRangeUser(0.1,3)
    rp195.GetLowerRefYaxis().SetLabelSize(0.025)
    rp195.GetLowerRefYaxis().SetTitle("%s_{data}/%s_{mcrecon}"%(compare_sym,compare_sym))
    rp195.GetUpperPad().cd()
    s195.Draw()
    legend.Draw()
    rp195.GetLowerRefGraph().Fit("pol1","W")
    rp195.SaveAs('{0}/{1}/data-mc-compare/{2}_195-20.png'.format(os.environ['DATA_ROOT_DIR'],det,plot))
    print('\nPress enter to continue.\n')
    os.system('read gonext')
    ###########################
    hmc20.SetLineColor(kRed)
    c220 = ROOT.TCanvas("c2{0}20".format(plot),"c2{0}20".format(plot),800,600)
    rp20 = ROOT.TRatioPlot(hdata20,hmc20,"divsym")
    rp20.SetH1DrawOpt("E0")
    rp20.SetH2DrawOpt("HIST")
    c220.cd()
    rp20.Draw()
    rp20.GetUpperRefYaxis().SetRangeUser(0,10)
    rp20.GetLowerRefYaxis().SetRangeUser(0.1,3)
    rp20.GetLowerRefYaxis().SetLabelSize(0.025)
    rp20.GetLowerRefYaxis().SetTitle("%s_{data}/%s_{mcrecon}"%(compare_sym,compare_sym))
    rp20.GetUpperPad().cd()
    s20.Draw()
    legend.Draw()
    rp20.GetLowerRefGraph().Fit("pol1","W")
    rp20.SaveAs('{0}/{1}/data-mc-compare/{2}_20-205.png'.format(os.environ['DATA_ROOT_DIR'],det,plot))
    print('\nPress enter to continue.\n')
    os.system('read gonext')
    ###########################

if __name__=='__main__':
    # Define variables
    det,plot,fd_path  = read_in_args()
    last_date_data,last_date_mc = find_last_date(fd_path)
    filename1 = '{2}/pass5/data/processing/root_all_files/20190625-{1}.{0}.ps5.ps2g.drmdpcgf.root'.format(det,last_date_data,fd_path)
    filename2 = '{2}/mc/processing/root_all_files/pr.1e17-1e21.20190625-{1}.{0}.ps5.drmdpcgf.root'.format(det,last_date_mc,fd_path)

    # Open MC file
    file1   = ROOT.TFile(filename1,"r")
    taTree1 = file1.Get("taTree")
    file2   = ROOT.TFile(filename2,"r")
    taTree2 = file2.Get("taTree")

    if plot == 'all':
        data_mc_compare(taTree1,taTree2,'phi')
        data_mc_compare(taTree1,taTree2,'theta')
        data_mc_compare(taTree1,taTree2,'psi')
        data_mc_compare(taTree1,taTree2,'phpgt')
        data_mc_compare(taTree1,taTree2,'track')
        ngmpe(taTree1,taTree2,det)
    elif plot != 'ngmpe':
        data_mc_compare(taTree1,taTree2,plot)
    else:
        ngmpe(taTree1,taTree2,det)
