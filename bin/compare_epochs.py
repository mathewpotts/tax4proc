#!/usr/bin/env python3

# Import Libs
import ROOT
import TDSTio
import TLPY
import os
import argparse

# PASS5 FILE PATHS FOR EACH EPOCH
epoch_path = os.environ['MDTAX4_DATA_PASS5'] + '/processing/root_all_files/'
epoch0     = epoch_path + 'epoch0_20180608-20190409.mdtax4ps5.ps2g.drmdpcgf.root'
epoch1     = epoch_path + 'epoch1_20190424-20190607.mdtax4ps5.ps2g.drmdpcgf.root'
epoch1_2   = epoch_path + 'epoch0-1.mdtax4.ps5.root'
epoch2     = epoch_path + '20190625-20201206.mdtax4.ps5.ps2g.drmdpcgf.root'

# Drawing in-plane angle (n_z)
files = [epoch2,epoch1_2,epoch0,epoch1]
for i,f in enumerate(files):
    print(i,f)
    file = ROOT.TFile('{0}'.format(f),'r')
    taTree = file.Get('taTree')
    #if i == 0:
        #taTree.Draw('stpln.n_ampwt[6][2]')
    #else:
        #taTree.Draw('stpln.n_ampwt[6][2]','','same')
    #ROOT.gPad.Update()

# Prompt exit
print('\nPress enter to continue.\n')
os.system('read gonext')



# 40 bins e0 vs e2
file = ROOT.TFile(epoch0,'r')
taTree_e0 = file.Get('taTree')
nentries_e0 = taTree_e0.GetEntries()
hepoch0 = ROOT.TH1D('hepoch0','Epoch 0 vs Epoch 2 cos^{-1}(n_{z}) comparison;cos^{-1}(n_{z}) degrees;',30,-90,90)
taTree_e0.Draw('(180/TMath::Pi())*(TMath::ACos(stpln.n_ampwt[6][2]))- 90>>hepoch0')

file1 = ROOT.TFile(epoch2,'r')
taTree_e2 = file1.Get('taTree')
nentries_e2 = taTree_e2.GetEntries()
norm = nentries_e0/nentries_e2
hepoch2 = ROOT.TH1D('hepoch2','hepoch2',30,-90,90)
taTree_e2.Draw('(180/TMath::Pi())*(TMath::ACos(stpln.n_ampwt[6][2]))- 90>>hepoch2','{0}'.format(norm),'p sames')
ROOT.gPad.BuildLegend()
ROOT.gPad.Update()

# Prompt exit
print('\nPress enter to continue.\n')
os.system('read gonext')

# 40 bins e0 vs e1
file = ROOT.TFile(epoch0,'r')
taTree_e0 = file.Get('taTree')
nentries_e0 = taTree_e0.GetEntries()
hepoch0 = ROOT.TH1D('hepoch0','Epoch 0 vs Epoch 1 cos^{-1}(n_{z}) comparison;cos^{-1}(n_{z}) degree;',30,-90,90)
taTree_e0.Draw('(180/TMath::Pi())*(TMath::ACos(stpln.n_ampwt[6][2])) - 90>>hepoch0')

file1 = ROOT.TFile(epoch1,'r')
taTree_e1 = file1.Get('taTree')
nentries_e1 = taTree_e1.GetEntries()
norm = nentries_e0/nentries_e1
hepoch1 = ROOT.TH1D('hepoch1','hepoch1',30,-90,90)
taTree_e1.Draw('(180/TMath::Pi())*(TMath::ACos(stpln.n_ampwt[6][2])) - 90>>hepoch1','{0}'.format(norm),'p sames')
ROOT.gPad.BuildLegend()
ROOT.gPad.Update()

# Prompt exit
print('\nPress enter to continue.\n')
os.system('read gonext')



    
