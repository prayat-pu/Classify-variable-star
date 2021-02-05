# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 14:16:50 2021

@author: Prayat
"""
# import sys
# sys.path.append(r'C:\Users\Prayat\Desktop\variable star-template fitting\program\lib')

import pandas as pd
from tqdm import tqdm
from lib.help_function_fitting import *
import warnings
warnings.filterwarnings("ignore")
import os


# Parameters 
#------------------------------------------------------------------------------#
RUNNING_FILE_NAME = "fitlc_star_name_txt.txt"
ERAT = 1.0
ERRMAX = 0.5
EMIN = 0.05
PERIOD1 = 0.1
PERIOD2 = 1.0
DPER = 0.1
#------------------------------------------------------------------------------#
def main():
    
    output_filename = 'all_result.txt'
    if os.path.isfile(output_filename):
        os.remove(output_filename)
    else:
        file = open(output_filename,'w+')
        file.close()
    
    all_output_result_file = open(output_filename,"a+")
    all_output_result_file.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n"%("Star_id","template_number",'period','v_amp','i_amp','v_med','i_med','score','type',))
        
    #=====================================================================#
    #constant variable
    running_file_name = RUNNING_FILE_NAME
    
    ALL_STAR_DATA_NAME = 'Data/data_for_rrfit/'+running_file_name
    STAR_DATA_FOLDER = 'Data/data_from_fitlc/csv/'
    TEMPLATE_PATH = 'Data/data_for_rrfit/Template/template.csv'
    
    template_df = pd.read_csv(TEMPLATE_PATH,usecols=[0,1,2,3,4,5,6,7,8])
    all_filename_list = pd.read_csv(ALL_STAR_DATA_NAME).star_filenames.to_list()

    #Parameter-------------------------------------#
    TEMPLATE_NUMBER = len(template_df.columns) - 1
    #----------------------------------------------#
    erat = ERAT  # default 1.0
    errmax = ERRMAX
    emin = EMIN
    
    filter_i = 'i'
    filter_v = 'v'

    template_number = TEMPLATE_NUMBER
    #=====================================================================#
    
    #=====================================================================#
    #template fitting each star
    print('process progress...')
    for i in tqdm(range(len(all_filename_list))):
        filename = all_filename_list[i]
         # --> set data
        STAR_PATH = STAR_DATA_FOLDER + filename.replace('fitlc','csv')
        star_df = pd.read_csv(STAR_PATH)
        
        # step1. check error of data
        star_df = star_df.loc[(star_df['err']<1)&(star_df['err']<errmax)]
        star_df.loc[star_df['err']<emin,'err'] = emin
        
        # add adjust erat to error for handle bad seeing
        star_df['err'] = star_df['err'] * erat
        
        # select only filter that you want with variable that create above
        star_df_v = star_df.loc[star_df['filter'] == filter_v]
        star_df_i = star_df.loc[star_df['filter'] == filter_i].reset_index(drop=True)
        
        
        vmin = star_df_v['mag'].min()
        v_t0 = star_df_v.loc[star_df_v['mag']==vmin,'hjd'].values[0]
        
        imin = star_df_i['mag'].min()
        i_t0 = star_df_i.loc[star_df_i['mag']==imin,'hjd'].values[0]
        
        template_dict_inform = dict()
        #----> set each period that choose boundary by period1,period2
        for round_ in range(template_number):
            period1 = PERIOD1
            period2 = PERIOD2
            dper = DPER
            current_template = round_ + 1
            
            period_out_dict = dict()
            
            periods = []
            periods.append(period1)
            period = period1
            while(period < period2):
                period = (1+dper)*period1
                if period > period2:
                    continue
                periods.append(period)
                period1 = period
            
            periods = [round(x,3)for x in periods] # list of all periods
        
            # find period of data
            for index in range(len(periods)):
                period = periods[index]
                #filter v
                v_cyc = (star_df_v[['hjd']] - v_t0)/ period
                v_cyc['hjd'] = v_cyc['hjd'] - v_cyc['hjd'].astype(int)
                
                v_cyc.loc[v_cyc['hjd']<0,'hjd'] = v_cyc.loc[v_cyc['hjd']<0,'hjd'] + 1
                v_cyc.loc[v_cyc['hjd']>1,'hjd'] = v_cyc.loc[v_cyc['hjd']>1,'hjd'] - 1
                
                star_df_v['phase'] = v_cyc
                
                #filter i
                i_cyc = (star_df_i[['hjd']] - i_t0)/ period
                i_cyc['hjd'] = i_cyc['hjd'] - i_cyc['hjd'].astype(int)
                
                i_cyc.loc[i_cyc['hjd']<0,'hjd'] = i_cyc.loc[i_cyc['hjd']<0,'hjd'] + 1
                i_cyc.loc[i_cyc['hjd']>1,'hjd'] = i_cyc.loc[i_cyc['hjd']>1,'hjd'] - 1
                
                star_df_i['phase'] = i_cyc
                
                #---->fitting data
                x_optimize,score = fitlc(star_df_v,star_df_i,template_df.copy(), current_template)
                
                # save fitting results
                period_out_dict[period] = (x_optimize,score)
            min_dict = find_min_score(period_out_dict)
            template_dict_inform[current_template] = min_dict
            
        best_template, best_period, best_x , best_score = find_min_score_all_template(template_dict_inform)
        
        #save result to output folder
        output_path = "output/"+filename.split('.')[0]
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        
        if not os.path.exists('output/img'):
            os.makedirs('output/img')
            
        save_star_data(output_path,filename,best_template,best_period,template_dict_inform.copy(),star_df_v.copy(),star_df_i.copy())
        plot_best_fitting(output_path,best_period,filename,best_template,template_dict_inform.copy(),template_df.copy(),star_df_v.copy(),star_df_i.copy())
        data_result = save_fitting_result(output_path,template_dict_inform.copy(),template_df.columns)
        
        save_all_result(best_template,output_filename,all_output_result_file,data_result,filename.split('.')[0])
        
    all_output_result_file.close()
    print('-----------complete fitting all star-----------------')
    #=====================================================================#
    

    
    
if __name__ == "__main__":
    main()
