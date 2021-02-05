# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 14:27:28 2021

@author: Prayat
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

EPOCH = 2450000
Nv = 12 # total data of filter v
Ni = 14 # total data of filter i
dvmag = [0.00,-0.04,0.05,0.04,0.06,0.04,-0.02,-0.09,0.01,-0.01,0.00,0.02]
dimag = [-0.49,-0.54,-0.44,-0.44,-0.44,-0.44,-0.47,-0.51,-0.43,-0.43,-0.43,-0.42,-0.43,-0.42]
flag1,flag2 = 0,1

def read_file(path):
    data = open(path,'r')
    return data    

def main():
    #check file
    file = "Data/data_for_rrfit/fitlc_star_name_txt.txt"
    if os.path.isfile(file):
        os.remove(file)
    
    new_file = open(file,"w+")
    new_file.write("star_filenames\n")
    new_file.close()
    # 1. read file
    montage_dir = 'Data/data_for_fitlc/montage.raw'
    vhjd_dir = 'Data/data_for_fitlc/vhjd.dat'
    ihjd_dir = 'Data/data_for_fitlc/ihjd.dat'
    
        # bring hjd both filter - EPOCH
    vhjd_data = read_file(vhjd_dir)
    vhjd = [float(x.strip()) - EPOCH for x in vhjd_data.readlines()]
    
    ihjd_data = read_file(ihjd_dir)
    ihjd = [float(x.strip()) - EPOCH for x in ihjd_data.readlines()]
    
    montage_data = read_file(montage_dir)
    vmag,imag,verr,ierr,filter_i,filter_v = [],[],[],[],[],[]
    filter_v_check = True
    

    
    #2. access data from montage.raw
    for i,line in enumerate(montage_data.readlines()):
        if i < 3:
            continue
        else:
            n_line = [float(x) for x in line.split()]
            if len(n_line) == 15: #first line 
                id_ = int(n_line[0])
                
                for index in range(3,Nv+3,2):
                    vmag.append(n_line[index])
                    verr.append(n_line[index+1])
                    filter_v.append('v')
            elif len(n_line) == 12:
                if filter_v_check:
                    for index in range(0,Nv,2):
                        vmag.append(n_line[index])
                        verr.append(n_line[index+1])
                        filter_v.append('v')
                    filter_v_check = False
                elif not(filter_v_check):
                    for index in range(0,Ni-2,2):
                        imag.append(n_line[index])
                        ierr.append(n_line[index+1])
                        filter_i.append('i')
            elif len(n_line) == 6:
                filter_v_check = True
                for index in range(0,4,2):
                    imag.append(n_line[index])
                    ierr.append(n_line[index+1])
                    filter_i.append('i')
                
                
                
                vmag = [vmag[x] + dvmag[x] for x in range(len(dvmag))]
                imag = [imag[x] + dimag[x] for x in range(len(dimag))]
                
          
                        
                vmag_avg = np.mean(vmag)
                imag_avg = np.mean(imag)
                
                # verr_sig = np.mean(verr)
                # ierr_sig = np.mean(ierr)
                
                vicol = vmag_avg - imag_avg
                
                if (13<vmag_avg<18) and (-1<vicol<0.1): 

                    vdevsum = 0
                    for j in range(Nv):
                        if vmag[j] != 99.999:
                            # adjust to z value 
                            vdev = (vmag[j] - vmag_avg)/verr[j]
                            vdevsum += vdev**2

                   
                    idevsum = 0
                    for j in range(Ni):
                        if imag[j] != 99.999:
                            idev = (imag[j] - imag_avg)/ierr[j]
                            idevsum += idev**2

                    
                    # find boundary of magnitude by mean and std
                    vmin = vmag_avg - (3 * np.std(vmag))
                    vmax = vmag_avg + (3 * np.std(vmag))
                    
                    imin = imag_avg - (2 * np.std(imag))
                    imax = imag_avg + (2 * np.std(imag))

                
                
                    vindex = (vdevsum+idevsum)/(Nv+Ni)
                    if vindex > 3:
                        if id_ <= 9 :
                            file_name = "Data/data_from_fitlc/txt/V00000"+str(id_)+".fitlc"
                        elif id_ >= 10 and id_ <= 99:
                            file_name = "Data/data_from_fitlc/txt/V0000"+str(id_)+".fitlc"
                        elif id_ >= 100 and id_ <= 999:
                            file_name = "Data/data_from_fitlc/txt/V000"+str(id_)+".fitlc"
                        elif id_ >= 1000 and id_ <= 9999:
                            file_name = "Data/data_from_fitlc/txt/V00"+str(id_)+".fitlc"
                        elif id_ >= 10000 and id_ <= 99999:
                            file_name = "Data/data_from_fitlc/txt/V0"+str(id_)+".fitlc"
                        elif id_ >= 100000 and id_ <= 999999:
                            file_name = "Data/data_from_fitlc/txt/V"+str(id_)+".fitlc"

                        
                        star_name_file = open(file,"a")
                        star_name_file.write(file_name.lstrip('Data/data_from_fitlc/txt/')+"\n")
                        
                        
                        star_detail_file = open(file_name,'w')
                        
                        for j in range(Nv):
                            if vmag[j] != 99.999:
                                if vmin < vmag[j] < vmax:
                                    star_detail_file.write("{0} {1:.5f} {2:.3f} {3:.3f}\n".format(flag1,vhjd[j],vmag[j],verr[j]))
                        
                        for p in range(Ni):
                                if imag[p] != 99.999:
                                    if imin < imag[p] < imax:
                                        star_detail_file.write("{0} {1:.5f} {2:.3f} {3:.3f}\n".format(flag2,ihjd[p],imag[p],ierr[p]))
               
                vmag,imag,verr,ierr,filter_v,filter_i = [],[],[],[],[],[]
                
            

    
if __name__ == "__main__":
    main()
