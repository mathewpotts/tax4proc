#!/usr/bin/env python3

# Read in a list of ASCII ontime files for MD/BR TAx4 and output total ontime for each  mirror and total ontime for the detector. 

# Import libs
import ROOT
import os
import glob
import sys
import numpy as np
import argparse

def read_in_ascii_list():
    parser = argparse.ArgumentParser(description = 'Read in a list of ASCII ontime files for MDTAx4 and output total ontime for each  mirror and total ontime for the detector.')
    parser.add_argument('-i',metavar='infile',action='store',help='A list of the ontime ASCII files that are to be used to calculate ontime.',required=True)
    parser.add_argument('-det',metavar='detector',action='store',help='Detector that you want to calculate the ontime for. (brtax4/mdtax4)',default='mdtax4')
    args = parser.parse_args()
    infile = args.i
    print('Reading %s...' % infile)
    det=args.det
    filelist=[line.rstrip('\n') for line in open(infile)]
    return infile,filelist,det

def calc_ontime(file,det):
    if det == 'mdtax4':
        max = 7
        mirror_ontime=[0,0,0,0]
    if det == 'brtax4':
        max = 11
        mirror_ontime=[0,0,0,0,0,0,0,0]
    line_list=[line for line in open(file)][0].split()
    start=int(line_list[0])
    stop=int(line_list[1])
    det_ontime=stop-start
    for idx in range(3,max):
        mirror_ontime[idx-3]=(stop-start)*float(line_list[idx])
    return mirror_ontime,det_ontime

def is_gw(file,det):
    date = file[-32:-28]+file[-27:-25]+file[-24:-22]
    part = int(file[-20:-18])
    if det == 'mdtax4':
        filename = file[-33:-20] + str(part).rjust(2,'0') + '.weather.mdtax4.txt'
        filename_prev = file[-33:-20] + str(part - 1).rjust(2,'0') + '.weather.mdtax4.txt'
        filepath = os.environ['MDTAX4_DATA_ROOT']+'/weather/'+date+'/'+filename
        filepath_prev = os.environ['MDTAX4_DATA_ROOT']+'/weather/'+date+'/'+ filename_prev
    if det == 'brtax4':
        filename = file[-33:-20] + str(part).rjust(2,'0') + '.weather.brtax4.txt'
        filename_prev = file[-33:-20] + str(part - 1).rjust(2,'0') + '.weather.brtax4.txt'
        filepath = os.environ['BRTAX4_DATA_ROOT']+'/weather/'+date+'/'+filename
        filepath_prev = os.environ['BRTAX4_DATA_ROOT']+'/weather/'+date+'/'+ filename_prev
    # If weat file doesn't exist for current part use weat code from last part
    try:
        line_list=[line for line in open(filepath)][0].split()
    except:
        line_list=[line for line in open(filepath_prev)][0].split() 
    weather_code=line_list[2]
    n= int(weather_code[0]) # North
    e= int(weather_code[1]) # East
    s= int(weather_code[2]) # South
    w= int(weather_code[3]) # West
    o= int(weather_code[4]) # Overcast
    t= int(weather_code[5]) # Thickness of clouds
    h= int(weather_code[6]) # Haze
    if det == 'mdtax4':
        if o+t+n+e == 0:
            return 1
        else:
            return 0
    if det == 'brtax4':
        if o+t+s+e == 0:
            return 1
        else:
            return 0

def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False
    
if __name__=='__main__':
    infile,filelist,det=read_in_ascii_list()
    total_mirror_ontime=[0,0,0,0,0,0,0,0]
    total_ontime=0
    gw_ontime=0
    tmp_mirror_ontime=[0,0,0,0,0,0,0,0]
    tmp_ontime=0
    if det == 'mdtax4':
        max = 4
        tmp_mirror_ontime=[0,0,0,0]
        total_mirror_ontime=[0,0,0,0]
    if det == 'brtax4':
        max = 8
        tmp_mirror_ontime=[0,0,0,0,0,0,0,0]
        total_mirror_ontime=[0,0,0,0,0,0,0,0]
    # Remove gw time from end of filelist
    filelist = [file for file in filelist if not isfloat(file)]

    # List varibles for the TGraph to be filled in for loop
    date_ontime_list   = []
    gw_ontime_list     = []
    by_day_ontime_list = []

    # Open each file and calculate the ontime for each mirror 
    # Make sure to delete epoch 0,1 ontime files or there will be error with MDTAx4
    for file in filelist:
        date = file[-32:-28]+file[-27:-25]+file[-24:-22]
        part = int(file[-20:-18])
        tmp_mirror_ontime,tmp_ontime=calc_ontime(file,det)
        total_ontime+=tmp_ontime
        gw_ontime+=is_gw(file,det)*tmp_ontime

        # Form a list of date and ontime for that day
        if (date == '20190625' or date == '20191022') and part == 1:
            date_ontime    = total_ontime
            gw_date_ontime = gw_ontime
        elif (date != '20190625' or date != '20191022') and part == 1:
            date_ontime    = total_ontime
            gw_date_ontime = gw_ontime
            date_ontime_list.append([date_prev,(total_ontime - tmp_ontime)/3600])
            gw_ontime_list.append([date_prev,(gw_ontime - tmp_ontime)/3600])
            if len(date_ontime_list) > 1:
                by_day_ontime_list.append([date_prev,float(date_ontime_list[-1][1]) - float(date_ontime_list[-2][1]),float(gw_ontime_list[-1][1]) - float(gw_ontime_list[-2][1])])
            else:
                by_day_ontime_list.append([date_prev,float(date_ontime_list[-1][1]),float(gw_ontime_list[-1][1])])
        elif date == '20201222' and part == 2:
            date_ontime_list.append([int(date),total_ontime/3600])
            gw_ontime_list.append([int(date),gw_ontime/3600])
            by_day_ontime_list.append([int(date),float(date_ontime_list[-1][1]) - float(date_ontime_list[-2][1]),float(gw_ontime_list[-1][1]) - float(gw_ontime_list[-2][1])])
        elif date != '20201222' and part != 1:
            date_ontime    += tmp_ontime
            gw_date_ontime += is_gw(file,det)*tmp_ontime
        date_prev = int(date)
        
        for mirror in range(0,max):
            total_mirror_ontime[mirror]+=tmp_mirror_ontime[mirror]

    # Drawing TGraph of ontime per day
    x = np.array(date_ontime_list).T[0]
    tmp = []
    for i in x:
        yr = int(i / 10000)
        mn = int(i / 100 % 100)
        d  = int(i % 100)
        tmp.append(ROOT.TDatime(yr,mn,d,0,0,0).Convert() - ROOT.TDatime(1995,1,1,0,0,0).Convert())
    floatx = np.array(tmp,dtype=np.float32) 
    y = np.array(date_ontime_list).T[1]
    y2= np.array(gw_ontime_list).T[1]
    y3= np.array(by_day_ontime_list).T[1]
    y4= np.array(by_day_ontime_list).T[2]
    floaty = np.array(y,dtype=np.float32)
    floaty2= np.array(y2,dtype=np.float32)
    floaty3= np.array(y3,dtype=np.float32)
    floaty4= np.array(y4,dtype=np.float32)
    ontime_per_day   = ROOT.TGraph(int(len(date_ontime_list)),floatx,floaty)
    gwtime_per_day   = ROOT.TGraph(int(len(gw_ontime_list)),floatx,floaty2)
    if det == 'mdtax4':
        daily_bin = 516
    if det == 'brtax4':
        daily_bin = 397
    by_day_ontime    = ROOT.TH1F('by_day_ontime','by_day_ontime',daily_bin,floatx[0],floatx[-1])
    by_gw_day_ontime = ROOT.TH1F('by_gw_day_ontime','by_gw_day_ontime',daily_bin,floatx[0],floatx[-1])
    for i,date in enumerate(floatx):
        by_day_ontime.Fill(date,floaty3[i])
        by_gw_day_ontime.Fill(date,floaty4[i])
    #by_day_ontime    = ROOT.TGraph(int(len(by_day_ontime_list)),floatx,floaty3)
    #by_gw_day_ontime = ROOT.TGraph(int(len(by_day_ontime_list)),floatx,floaty4)
    ontime_per_day.GetXaxis().SetTimeDisplay(1)
    ontime_per_day.GetXaxis().SetTimeFormat("%Y/%m/%d")
    ontime_per_day.GetXaxis().SetLabelSize(0.022)
    ontime_per_day.GetYaxis().SetTitle('Hours')
    ontime_per_day.SetMarkerStyle(20)
    ontime_per_day.SetMarkerSize(0.5)
    ontime_per_day.SetTitle('{0} Ontime'.format(det.upper()))
    ontime_per_day.Draw("apl")
    gwtime_per_day.SetMarkerStyle(20)
    gwtime_per_day.SetMarkerSize(0.5)
    gwtime_per_day.SetMarkerColor(2)
    gwtime_per_day.SetLineColor(2)
    gwtime_per_day.Draw("pl")
    by_day_ontime.SetFillColor(38)
    #by_day_ontime.Draw('same,hist')
    by_gw_day_ontime.SetFillColor(2)
    #by_gw_day_ontime.Draw('same,hist')
    #ROOT.gPad.SetLogy()
    ROOT.gPad.SaveAs("{0}/plots/ontime_per_day.png".format(os.environ['{0}_DATA_ROOT'.format(det.upper())]))
    

    # Prompt exit
    print('\nPress enter to continue.\n')
    os.system('read gonext')
    
    if det == 'mdtax4':
        mirror_dict={0:'m25',1:'m26',2:'m27',3:'m28'}
        for mirror in range(0,len(total_mirror_ontime)):
            print('The total mirror ontime for %s   :  %10.4f' % (mirror_dict[mirror],total_mirror_ontime[mirror]/3600),'hrs')
    if det == 'brtax4':
        mirror_dict_brtax4={0:'m29',1:'m30',2:'m31',3:'m32',4:'m33',5:'m34',6:'m35',7:'m36'}
        for mirror in range(0,len(total_mirror_ontime)):
            print('The total mirror ontime for %s   :  %10.4f' % (mirror_dict_brtax4[mirror],total_mirror_ontime[mirror]/3600),'hrs')
    print('The total ontime of the detector  :  %10.4f' % (total_ontime/3600),'hrs')
    print('The GW ontime of the detector : %10.4f' % (gw_ontime/3600),'hrs')

    # Print GW time to the end of the ontime file
    with open(infile,'r') as f:
        lines = f.readlines()
    with open(infile,'w') as f:
        for line in lines:
            if not isfloat(line):
                f.write(line)
    print(gw_ontime/3600,file=open(infile,'a'))
