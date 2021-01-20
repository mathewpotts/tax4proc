#!/bin/env python3

import ROOT
import TLPY
from TLPY import *
import TDSTio
import array
import numpy as np
from os import system
from sys import exit
import numpy as np
import argparse

def read_in_args():
    parser = argparse.ArgumentParser(description = 'Add the branch to taTree that flags event as quality.')
    parser.add_argument('-det',metavar='detector',action='store',help='Detector that you want that the root file is for. (brtax4/mdtax4)',required=True)
    parser.add_argument('-i',metavar='infile',action='store',help='Input ROOT file that needs quality events flagged.',required=True)
    args     = parser.parse_args()
    det      = args.det
    infile   = args.i
    #datatype = args.t
    if det == 'mdtax4':
        det_id = 6
        sd_id  = 1 # TAX4SD_208_N
    if det == 'brtax4':
        det_id = 7
        sd_id  = 2 # TAX4SD_208_S
    return det,det_id,sd_id,infile

def get_border_distance(fdsiteid,sd_det_id,nx,ny,nz,rp,psi):
    # Convert FD core coordinates to CLF core coordinates
    sdp_n = np.array([nx,ny,nz])
    t_fd  = tlevent_track(sdp_n,rp,psi,0)
    t_clf = tlevent_track(TLUTI_REC.fd2clf(fdsiteid,t_fd))
    xcore = t_clf.xyz[0]
    ycore = t_clf.xyz[1]
    bdist = TLUTI_DETECTORS.get_sd_bdist(sd_det_id, xcore, ycore)
    return bdist

det,det_id,sd_id,infile = read_in_args()

# Load ROOT file
filename = infile
file     = ROOT.TFile(filename,"update")
taTree   = file.Get("taTree")

# Define starting variables
events                    = taTree.GetEntries()
bracketing_cut_loss       = 0
border_cut_loss           = 0
core_position_cut_loss    = 0
event_saturation_cut_loss = 0
tracklength_cut_loss      = 0
profile_chi2dof_cut_loss  = 0
theta_cut_loss            = 0
psi_cut_loss              = 0
dep_first_cut_loss        = 0

# Call in external C script with quality cut functions
ROOT.gROOT.LoadMacro("qualct.C")

# Created an array to input values into a new branch
leafValue = array.array("f",[0.,0.])
qualct = taTree.Branch("qualct",leafValue,"qualct\F")

for i in range(events):
    taTree.GetEntry(i)

    # Getting Values from taTree
    dep_first = taTree.prfc.dep[det_id][0]
    dep_last  = taTree.prfc.dep[det_id][taTree.prfc.nbin[det_id]]
    xmax      = taTree.prfc.xm[det_id]
    nx        = taTree.hctim.mshwn[det_id][0]
    ny        = taTree.hctim.mshwn[det_id][1]
    nz        = taTree.hctim.mshwn[det_id][2]
    rp        = taTree.hctim.mrp[det_id]
    psi       = taTree.hctim.mpsi[det_id]
    ig        = taTree.stpln.ig
    xcore_SD  = taTree.tlfptn.tyro_xymoments[0][0]
    ycore_SD  = taTree.tlfptn.tyro_xymoments[0][1]
    track     = taTree.stpln.tracklength[det_id]
    prfc_chi2 = taTree.prfc.chi2[det_id]
    prfc_dof  = taTree.prfc.ndf[det_id]
    theta     = taTree.hctim.mthe[det_id]

    # Applying functions that part of C macro
    bracketing_cut        = ROOT.bracketing_cut(dep_first,dep_last,xmax,10,10)
    bdist                 = get_border_distance(det_id,sd_id,nx,ny,nz,rp,psi)
    if bdist > 100:
        border_cut = 1.0
    else:
        border_cut = 0.0
    core_position_cut = ROOT.is_core_position_diff(det_id,xcore_SD,ycore_SD,nx,ny,nz,rp,psi,2080)
    event_saturation_cut  = ROOT.event_saturation_cut(ig,20)

    # Applying misc cuts
    if track > 8.0:
        tracklength_cut = 1.0
    else:
        tracklength_cut = 0.0
    if (prfc_chi2/prfc_dof) <= 20:
        profile_chi2dof_cut = 1.0
    else:
        profile_chi2dof_cut = 0.0
    if theta <= np.deg2rad(65):
        theta_cut = 1.0
    else:
        theta_cut = 0.0
    if psi < np.deg2rad(120):
        psi_cut = 1.0
    else:
        psi_cut = 0.0
    if dep_first < 1000:
        dep_first_cut = 1.0
    else:
        dep_first_cut = 0.0


    # Checking if cuts drops event
    if bracketing_cut == 0:
        bracketing_cut_loss += 1
    if border_cut == 0:
        border_cut_loss += 1
    if core_position_cut == 0:
        core_position_cut_loss += 1
    if event_saturation_cut == 0:
        event_saturation_cut_loss += 1
    if tracklength_cut == 0:
        tracklength_cut_loss += 1
    if profile_chi2dof_cut == 0:
        profile_chi2dof_cut_loss += 1
    if theta_cut == 0:
        theta_cut_loss += 1
    if psi_cut == 0:
        psi_cut_loss += 1
    if dep_first_cut == 0:
        dep_first_cut_loss += 1
        
    leafValue[0] = bracketing_cut * border_cut * core_position_cut * event_saturation_cut * tracklength_cut * profile_chi2dof_cut * theta_cut * psi_cut * dep_first_cut
    qualct.Fill()
    

taTree.Write()

print('''
******** Number of Events Cut ********
  Total Number of Events: {0}
  Bracketing Cut        : {1}
  Border Cut            : {2}
  Core Position Cut     : {3}
  Event Saturation Cut  : {4}
  Tracklength Cut       : {5}
  Profile chi2/dof Cut  : {6}
  Theta Cut             : {7}
  Psi Cut               : {8}
  First interaction cut : {9}
**************************************
'''.format(events,bracketing_cut_loss,border_cut_loss,core_position_cut_loss,event_saturation_cut_loss,tracklength_cut_loss,profile_chi2dof_cut_loss,theta_cut_loss,psi_cut_loss,dep_first_cut_loss))

#taTree.Draw("(log((missing_E_corr(prfc->eng[6]))/mc04->energy)):log10(mc04->energy)","qualct ==1","prof")
system('read gonext')
