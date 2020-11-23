#!/usr/bin/env python3

# Read in a list of ASCII ontime files for MD/BR TAx4 and output total ontime for each  mirror and total ontime for the detector. 

# Import libs
import os
import glob
import sys
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
        if o+t+s+w == 0:
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

    # Open each file and calculate the ontime for each mirror 
    for file in filelist:
        #print(file)
        tmp_mirror_ontime,tmp_ontime=calc_ontime(file,det)
        total_ontime+=tmp_ontime
        gw_ontime+=is_gw(file,det)*tmp_ontime
        for mirror in range(0,max):
            total_mirror_ontime[mirror]+=tmp_mirror_ontime[mirror]
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
