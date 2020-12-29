def concat_dst(fd,fd_path,date,det):
    # Generate directories
    os.system("mkdir -p {0}/processing/concat_dst/".format(fd_path))
    os.system("mkdir -p {0}/processing/concat_drmdpcgf_dst/".format(fd_path))
    os.system("mkdir -p {0}/processing/root_all_files/".format(fd_path))
    
    # Concatinate files
    os.system("dstcat.run -o {0}/processing/concat_dst/20190625-{1}.{2}.ps5.ps2g.dst.gz {0}/[0-9]*[0-9]/*ps2g.dst.gz".format(fd_path,date,det))

    # Apply DCR to ps5 file
    os.system("dcrmdpcgf.run -f -fd {0} {1}/processing/concat_dst/20190625-{2}.{3}.ps5.ps2g.dst.gz -o {1}/processing/concat_drmdpcgf_dst/".format(fd,fd_path,date,det))
    
    # Convert the DST files into ROOT trees
    os.system("dst2rt_ta.run -f {0}/processing/concat_drmdpcgf_dst/20190625-{1}.{2}.ps5.ps2g.drmdpcgf.dst.gz -o {0}/processing/root_all_files/".format(fd_path,date,det))

def find_last_date(fd_path):
    file = os.popen('ls -drt {0}/[0-9][0-9]* | sort -n | tail -n1'.format(fd_path)).read()
    date = file[-9:-1]
    return date
