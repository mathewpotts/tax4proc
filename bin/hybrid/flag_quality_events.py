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

#Creating the tax4 branch class
ROOT.gROOT.ProcessLine('''
struct tax4 { 
    Int_t    fit_num;        // best fit number
    Int_t    qualct;         // 1 means passed choosen quality cuts
    Float_t  fd_clf_xcore;   // FD xcore position in CLF corrdinates (meters)
    Float_t  fd_clf_ycore;   // FD ycore position in CLF corrdinates (meters)
    Float_t  bdist;          // SD border distance (meters) 
};
''')
tax4 = ROOT.tax4()

def is_core_pos_diff(fdsiteid,xcore_SD,ycore_SD,nx,ny,nz,rp,psi,limit):
    sdp_n            = np.array([nx,ny,nz])
    t_fd             = tlevent_track(sdp_n,rp,psi,0)
    t_clf            = tlevent_track(TLUTI_REC.fd2clf(fdsiteid,t_fd))
    xcore            = t_clf.xyz[0]
    ycore            = t_clf.xyz[0]
    xcore_diff       = abs(xcore - xcore_SD)
    ycore_diff       = abs(ycore - ycore_SD)
    is_core_pos_diff = 1 if xcore_diff > limit and ycore_diff > limit else 0
    return is_core_pos_diff,xcore,ycore

def get_border_distance(fdsiteid,sd_det_id,nx,ny,nz,rp,psi):
    # Convert FD core coordinates to CLF core coordinates
    sdp_n = np.array([nx,ny,nz])
    t_fd  = tlevent_track(sdp_n,rp,psi,0)
    t_clf = tlevent_track(TLUTI_REC.fd2clf(fdsiteid,t_fd))
    xcore = t_clf.xyz[0]
    ycore = t_clf.xyz[1]
    bdist = TLUTI_DETECTORS.get_sd_bdist(sd_det_id, xcore, ycore)
    return bdist

def is_xm_FOV(dep_first,dep_last,xmax,dx1 = 0,dx2 = 0):
    lower_bracket = dep_first + dx1
    upper_bracket = dep_last - dx2
    return True if lower_bracket < xmax and upper_bracket > xmax else False

def is_event_saturated(tubeArray,limit = 20):
    nentries             = len(tubeArray)
    saturation_count     = 0.0
    for i in range(nentries):
        saturation_count += 1.0 if tubeArray[i] == -2 else 0.0
    saturation_percent   = saturation_count / nentries
    saturation_plimit    = limit / 100.
    return True if saturation_percent > saturation_plimit else False

def get_fit_num(taTree):
    FIT_START    = 10               # used to determine imin[] values
    HCTIM_MAXFIT = 16               
    ndfp         = 0.0              
    ndft         = 0.0              
    avgcfc       = [0]*HCTIM_MAXFIT # determination of a good fit (avg. correction fact)
    ndft_fit     = [0]*HCTIM_MAXFIT # # degrees of freedom in timing fit
    imin         = [0]*3            # the minimum chi-squares of combined, profile, and timing
    igfit        = [0]*HCTIM_MAXFIT # determination of a good fit (fit good)
    c2p          = 1.0e20           
    c2t          = 1.0e20           
    c2c          = 1.0e20           
    ################ IMINC DETERMINATION (for event) ################
    # number of good tubes as determined by plane-fit
    ng = 0
    for itube in range(0,taTree.stpln.ntube):
        ng += 1 if taTree.stpln.ig[itube] == 1 else 0
    for ifit in range(HCTIM_MAXFIT):
        # determine if particular fit was good
        igfit[ifit] = 1 if taTree.hctim.timinfo[ifit] == 1 and taTree.hctim.failmode[ifit] == 0 else 0
        # determine average correction factor for particular fit
        nbins,sums = 0,0
        for i in range(taTree.hcbin.nbin[ifit]):
            if taTree.hcbin.ig[ifit][i] != 1:
                continue
            sums   += taTree.hcbin.cfc[ifit][i]
            nbins += 1
        if nbins > 0:
            sums /= nbins
        avgcfc[ifit] = sums
        # number of degrees of freedom for the timing fit
        ndft_fit[ifit] = ng - 3
    # determine imin[] values: combined, profile, timing   
    for ifit in range(FIT_START,HCTIM_MAXFIT):
        if igfit[ifit] != 1:
            continue
        if taTree.prfc.chi2[ifit] < c2p:
            c2p = taTree.prfc.chi2[ifit]
            ndfp = taTree.prfc.ndf[ifit]
            imin[1] = ifit
        if taTree.hctim.mchi2[ifit] < c2t:
            c2t = taTree.hctim.mchi2[ifit]
            ndft = ndft_fit[ifit]
            imin[2] = ifit
    normp = 1.0 if c2p == 0.0 else ndfp/c2p
    normt = 1.0 if c2t == 0.0 else ndft/c2t
    for ifit in range(FIT_START,HCTIM_MAXFIT):
        if igfit[ifit] != 1:
            continue
        c2p = normp * taTree.prfc.chi2[ifit]
        c2t = normt * taTree.hctim.mchi2[ifit]
        if (c2p + c2t) < c2c:
            c2c = c2p + c2t
            imin[0] = ifit
    return imin[0],avgcfc

def dcr_qualcut(infile,det_id,sd_id):
    ##########################################################################
    # DCR Quality cuts:
    # Pass DCR During Reconstruction
    #  *Remove upward events
    #  *Track length >= 4.5⁰
    #  *0.75/⁰ =< Good tubes per degree =< 5.0/⁰ 
    #  *Average photons per good tube >= 50
    #  *RMS Deviation in the off-plane angle =< 1.2⁰
    #  *Inverse angular speed >= 0.05 µs/◦
    # Pass DCR Post Reconstruction 
    #  *Angular tracks are ≥ 7.9⁰
    #  *The depth of the first interaction in the atmosphere happens before 1000 g/cm<sup>2</sup>
    #  *𝜓 is < 120⁰
    #  *The average correction factor (amount of fitted track away from the mirror edges) is > 0.9 m^2
    #  *prfc chi2/dof <= 20
    ##########################################################################
    filename = infile
    file     = ROOT.TFile(filename,"update")                                              # Load ROOT file
    taTree   = file.Get("taTree")
    nevents   = taTree.GetEntries()                                                       # Define starting variables
    taTree.Branch("tax4",tax4,"fit_num/I:qualct/I:fd_clf_xcore/F:fd_clf_ycore/F:bdist/F") # Created a new branch
    branches = [branch.GetName() for branch in taTree.GetListOfBranches()]
    for i in range(nevents):
        taTree.GetEntry(i)
        if 'tlfptn' not in branches:
            fit_num,avgcfc = get_fit_num(taTree)
            tax4.fit_num   = fit_num
        else:
            print("fit_num set to detector id..")
            fit_num,avgcfc = det_id,1.0
            tax4.fit_num   = fit_num
        # Getting Values from taTree
        upward                  = taTree.stps2.upward[det_id]
        plog                    = taTree.stps2.plog[det_id]
        dep_first               = taTree.prfc.dep[fit_num][0]
        nx                      = taTree.hctim.mshwn[fit_num][0]
        ny                      = taTree.hctim.mshwn[fit_num][1]
        nz                      = taTree.hctim.mshwn[fit_num][2]
        rp                      = taTree.hctim.mrp[fit_num]
        psi                     = taTree.hctim.mpsi[fit_num]
        track                   = taTree.stpln.tracklength[det_id]
        rmsdvepln               = taTree.stpln.rmsdevpln[det_id]
        ph_per_gtube            = taTree.stpln.ph_per_gtube[det_id]
        prfc_chi2               = taTree.prfc.chi2[fit_num]
        prfc_dof                = taTree.prfc.ndf[fit_num]
        ig                      = taTree.stpln.ig
        stpln_ngtube            = ig.count(1)
        stpln_ngtube_per_degree = stpln_ngtube / track
        try:
            xcore_SD,ycore_SD  = taTree.tlfptn.tyro_xymoments[0][0],taTree.tlfptn.tyro_xymoments[0][1]
        except:
            xcore_SD,ycore_SD = 0,0 
        tax4.bdist                = get_border_distance(det_id,sd_id,nx,ny,nz,rp,psi)
        tmp,tax4.xcore,tax4.ycore = is_core_pos_diff(det_id,xcore_SD,ycore_SD,nx,ny,nz,rp,psi,2080)
        
        # Applying misc qual
        downward_qual          = True if upward == 1 else False
        plog_qual              = True if plog >= 2 else False
        rmsdevpln_qual         = True if rmsdevpln <= 1.2 else False
        tracklength_qual       = True if track > 7.9 else False
        ias_qual               = True if (crosstim/track) >= 0.05 else False
        ph_per_gtube_qual      = True if ph_per_gtube >= 50.0 else False
        profile_chi2dof_qual   = True if (prfc_chi2/prfc_dof) <= 20 else False
        psi_qual               = True if psi < np.deg2rad(120) else False
        ngtube_per_degree_qual = True if stpln_ngtube_per_degree >= 0.75 and stpln_ngtube_per_degree <= 5.0 else False
        avgcfc_qual            = True if avgcfc[fit_num] > 0.9 else False
        dep_first_qual         = True if dep_first < 1000 else False
        
        tax4.qualct = downward_qual*plog_qual*rmsdevpln_qual*tracklength_qual*ias_qual*ph_per_gtube_qual*profile_chi2dof_qual*psi_qual*ngtube_per_degree_qual*avgcfc_qual*dep_first_qual
        taTree.Fill()
    file.Write()
    file.Close()

def gregs_tax4_qualcut(infile,det_id,sd_id):
    ##########################################################################
    # Greg's TAx4 MD FD Quality cuts (based on HiRes-1):
    # Rayleigh Filter    P_{log10} >=2
    # Brightness Cut     sum(phtons)/N_{gt} >= 200   
    # Track Length       stpln.tracklength[det_id] > 7.9 
    # Track Width RMS    <= 1 degree
    # Angular Speed      <= 5.73 degree/microsec 
    # Profile Fit        prfc.chi2[fit_num]/prfc.dof[fit_num] < 14
    # Cherenkov Fraction crnk_sum / sig_sum < 0.20
    # First Interaction  prfc.dep[fit_num][0] <= 1200
    ##########################################################################
    filename = infile
    file     = ROOT.TFile(filename,"update") # Load ROOT file
    taTree   = file.Get("taTree")
    nevents  = taTree.GetEntries()
    # Define starting variables
    branch     = taTree.Branch("tax4",tax4,"fit_num/I:qualct/I:fd_clf_xcore/F:fd_clf_ycore/F:bdist/F") # Created a new branch
    branches = [branch.GetName() for branch in taTree.GetListOfBranches()]
    hyb      = 'tlfptn' in branches
    if hyb:
        print("ROOT tree contains hybrid information. Setting hybrid mode...")
    for i in range(0,nevents):
        taTree.GetEntry(i)
        if not hyb:
            fit_num,avgcfc = get_fit_num(taTree)
            tax4.fit_num   = fit_num
        else:
            print("fit_num set to detector id..")
            fit_num,avgcfc = det_id,1.0
            tax4.fit_num   = fit_num
        # Getting Values from taTree
        dep_first    = taTree.prfc.dep[fit_num][0]
        track        = taTree.stpln.tracklength[det_id]
        crosstime    = taTree.stpln.crossingtime[det_id]
        prfc_chi2    = taTree.prfc.chi2[fit_num]
        prfc_dof     = taTree.prfc.ndf[fit_num]
        theta        = taTree.hctim.mthe[fit_num]
        rmsdevpln    = 57.296 * taTree.stpln.rmsdevpln[det_id]
        plog         = taTree.stps2.plog[det_id]
        ph_per_gtube = taTree.stpln.ph_per_gtube[det_id]
        nx           = taTree.hctim.mshwn[fit_num][0]
        ny           = taTree.hctim.mshwn[fit_num][1]
        nz           = taTree.hctim.mshwn[fit_num][2]
        rp           = taTree.hctim.mrp[fit_num]
        psi          = taTree.hctim.mpsi[fit_num]
        if hyb:
            xcore_SD,ycore_SD = taTree.tlfptn.tyro_xymoments[0][0],taTree.tlfptn.tyro_xymoments[0][1]
        else:
            xcore_SD,ycore_SD = 0,0
        # Reset Sum values
        sig_sum   = 0
        scin_sum  = 0
        crnk_sum  = 0
        aero_sum  = 0
        rayl_sum  = 0
        for tube in range(0,len(taTree.prfc.sig[fit_num])):
            sig_sum  += taTree.prfc.sig[fit_num][tube]
            scin_sum += taTree.prfc.scin[fit_num][tube]
            crnk_sum += taTree.prfc.crnk[fit_num][tube]
            aero_sum += taTree.prfc.aero[fit_num][tube]
            rayl_sum += taTree.prfc.rayl[fit_num][tube]
        try:
            cv_frac = (crnk_sum) / sig_sum
        except:
            print("Signal sum is zero!! Setting Cherenkov Fraction to 1")
            print("Event:",i,"Signal Sum:",sig_sum)
            cv_frac = 1
        # Flagging Quality Events
        plog_qual                 = True if plog >= 2.0 else False
        ph_per_gtube_qual         = True if ph_per_gtube >= 200 else False
        track_qual                = True if track > 7.9 else False
        as_qual                   = True if track/crosstime <= 5.73 else False
        rmsdevpln_qual            = True if rmsdevpln < 1.0 else False
        try:
            prof_chi2dof_qual     = True if (prfc_chi2/prfc_dof) < 14 else False
        except:
            print("prfc DOF equals 0!! Setting prof_chi2dof_qual to FALSE")
            print("Event:",i,"prfc fit DOF:",prfc_dof)
            prof_chi2dof_qual     = False
        cv_frac_qual              = True if cv_frac < 0.20 else False
        dep_first_qual            = True if dep_first <= 1200 else False
        tax4.bdist                = get_border_distance(det_id,sd_id,nx,ny,nz,rp,psi)
        try:
            tmp,tax4.xcore,tax4.ycore = is_core_pos_diff(det_id,xcore_SD,ycore_SD,nx,ny,nz,rp,psi,2080)
        except:
            print("Error with core position function!!")
            tax4.xcore,tax4.ycore     = 0,0
        tax4.qualct               = plog_qual * ph_per_gtube_qual * track_qual * as_qual * rmsdevpln_qual * prof_chi2dof_qual * cv_frac_qual * dep_first_qual
        branch.Fill()
        pass
    taTree.Write()
    file.Close()

def current_hybrid_qualcut(infile,det_id,sd_id):
    ##########################################################################
    # My TAx4 Hybrid Quality cuts (based on Monica and Elliot theses):
    # Bracketing         XM is in FOV
    # Border Distance    Core of Shower is in 100m of array border
    # Core Position      Are the SD and FD calculated different by more than 2080m
    # Event Saturation   Are more than 20% of the tubes labeled as saturated
    # Track Length       stpln.tracklength[det_id] > 8 
    # Psi                psi < 120 
    # Profile Fit        prfc.chi2[fit_num]/prfc.dof[fit_num] < 14
    # Theta              theta <= 65 degrees
    # First Interaction  prfc.dep[det_id][0] <= 1000
    ##########################################################################
    filename = infile
    file     = ROOT.TFile(filename,"update")                                              # Load ROOT file
    taTree   = file.Get("taTree")
    nevents   = taTree.GetEntries()                                                       # Define starting variables
    taTree.Branch("tax4",tax4,"fit_num/I:qualct/I:fd_clf_xcore/F:fd_clf_ycore/F:bdist/F") # Created a new branch
    branches = [branch.GetName() for branch in taTree.GetListOfBranches()]
    for i in range(nevents):
        taTree.GetEntry(i)
        if 'tlfptn' not in branches:
            tax4.fit_num,avgcfc = get_fit_num(taTree)
        else:
            tax4.fit_num,avgcfc = det_id,1.0
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
        track     = taTree.stpln.tracklength[det_id]
        prfc_chi2 = taTree.prfc.chi2[det_id]
        prfc_dof  = taTree.prfc.ndf[det_id]
        theta     = taTree.hctim.mthe[det_id]
        try:
            xcore_SD,ycore_SD = taTree.tlfptn.tyro_xymoments[0][0],taTree.tlfptn.tyro_xymoments[0][1]
        except:
            xcore_SD,ycore_SD = 0,0
        # Flag quality events
        bracketing_qual                          = is_xm_FOV(dep_first,dep_last,xmax,10,10)
        tax4.bdist                               = get_border_distance(det_id,sd_id,nx,ny,nz,rp,psi)
        border_qual                              = True if tax4.bdist > 100 else False
        core_position_qual                       = 0 if core_position_bool else 1
        event_saturation_qual                    = not is_event_saturated(ig,20)
        tracklength_qual                         = True if track > 8.0 else False
        profile_chi2dof_qual                     = True if (prfc_chi2/prfc_dof) <= 20 else False
        theta_qual                               = True if theta <= np.deg2rad(65) else False
        psi_qual                                 = True if psi < np.deg2rad(120) else False
        dep_first_qual                           = True if dep_first < 1000 else False
        try:
            core_position_bool,tax4.xcore,tax4.ycore = is_core_pos_diff(det_id,xcore_SD,ycore_SD,nx,ny,nz,rp,psi,2080)
        except:
            core_position_bool,tax4.xcore,tax4.ycore = 0,0,0
        tax4.qualct = bracketing_qual * border_qual * core_position_qual * event_saturation_qual * tracklength_qual * profile_chi2dof_qual * theta_qual * psi_qual * dep_first_qual
        taTree.Fill()
        pass
    file.Write()
    file.Close()

if __name__ == '__main__':
    def read_in_args():
        parser = argparse.ArgumentParser(description = 'Add the branch to taTree that flags event as quality.')
        parser.add_argument('-det',metavar='detector',action='store',help='Detector that you want that the root file is for. (brtax4/mdtax4)',required=True)
        parser.add_argument('-q',metavar='Quality Cuts',action='store',help='Choose the quality cuts you want to use. (mine/gregs)',required=True)
        parser.add_argument('-i',metavar='infile',action='store',help='Input ROOT file that needs quality events flagged.',required=True)
        args     = parser.parse_args()
        det      = args.det
        infile   = args.i
        qualct   = args.q
        if det == 'mdtax4':
            det_id = 6
            sd_id  = 1 # TAX4SD_208_N
        if det == 'brtax4':
            det_id = 7
            sd_id  = 2 # TAX4SD_208_S
        return det,det_id,sd_id,infile,qualct

    det,det_id,sd_id,infile,qualct = read_in_args()
    print('''
    RUNNING flag_quality_events.py
    FILE   : {0}
    DET ID : {1}
    SD ID  : {2}
    QUALCT : {3}
    '''.format(infile,det_id,sd_id,qualct))

    # Call in external C script with quality cut functions
    ROOT.gROOT.LoadMacro("qualct.C")

    #Fill qualct tree
    if qualct == 'mine':
        current_hybrid_qualcut(infile,det_id,sd_id)
    elif qualct == 'gregs':
        gregs_tax4_qualcut(infile,det_id,sd_id)
    elif qualct == 'dcr':
        dcr_qualcut(infile,det_id,sd_id)
