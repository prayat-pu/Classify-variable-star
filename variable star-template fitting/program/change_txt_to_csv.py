# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 14:43:23 2021

@author: Prayat
"""


import glob
# import numpy as np
import pandas as pd

def handle_file(path):
    star_data = open(path,'r')
    col = ['filter','hjd','mag','err']
    df1 = pd.DataFrame(columns=col)
    df2 = pd.DataFrame(columns=col)

    flag1,vhjd,vmag,verr = [],[],[],[]
    flag2,ihjd,imag,ierr = [],[],[],[]

    for i,data in enumerate(star_data.readlines()):
        line_data_ = data.split()
        line_data = [float(info) for info in line_data_]
        flag = int(line_data[0])
        hjd = line_data[1]
        mag = line_data[2]
        err = line_data[3]
        if flag == 0:
            flag1.append(flag)
            vhjd.append(hjd)
            vmag.append(mag)
            verr.append(err)
        elif flag == 1:
            flag2.append(flag)
            ihjd.append(hjd)
            imag.append(mag)
            ierr.append(err)


    df1['filter'] = flag1
    df1['hjd'] = vhjd
    df1['mag'] = vmag
    df1['err'] = verr

    df2['filter'] = flag2
    df2['hjd'] = ihjd
    df2['mag'] = imag
    df2['err'] = ierr

    df = pd.concat([df1,df2])
    df.reset_index(inplace=True,drop=True)

    # df.replace(0,'v',inplace=True)
    # df.replace(1,'i',inplace=True)
    df.loc[df['filter'] == 0,'filter'] = 'v'
    df.loc[df['filter'] == 1, 'filter'] = 'i'

    # star_name_file = open("Data/data_for_rrfit/fitlc_star_name_csv.txt","a+")
    # star_name_file.write(path.split('.')[0][-7:]+'.csv\n')


    df.to_csv('Data/data_from_fitlc/csv/'+path.split('.')[0][-7:]+'.csv',index=False)



def main():
    # col1  col2    col3    col4
    # flag  hjd     mag     err
    # flag 0 = v, 1 = i
    path = 'Data/data_from_fitlc/txt/*'
    for name in glob.glob(path):
        handle_file(name)

if __name__ == "__main__":
    main()