# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 15:28:48 2020

@author: Prayat

"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import minimize


def create_plot(row,col,size,title,df):
    fig,ax = plt.subplots(nrows=row,ncols=col,figsize=size)
    i = 0
    for r in range(len(ax)):
        for c in range(len(ax[r])):
            ax[r][c].plot(df[df.columns[i]])
            ax[r][c].set_title(df.columns[i])
            ax[r][c].invert_yaxis()
            i += 1
def find_min_chi(dict_):
    
    # dict format 
    #  period : ((x_optimize1,x_optimize2), score)
    return min(dict_.items(),key= lambda x:x[1])
            
def find_min_score(dict_):
    
    # dict format 
    #  period : ((x_optimize1,x_optimize2), score)
    return min(dict_.items(),key= lambda x:x[1][1])

def find_max_score(dict_):
    return max(dict_.items(),key= lambda x:x[1][1])


def find_min_score_all_template(dict_):
    # dict format
    # template_number: (period,((x_optimize1,x_optimize2),score))
    dict_min = min(dict_.items(), key=lambda x:x[1][1][1])
    best_template = dict_min[0]
    best_period = dict_min[1][0]
    best_x = dict_min[1][1][0]
    best_score = dict_min[1][1][1]
    return best_template, best_period, best_x, best_score

def vintrp(phase,sph,spl,svh,svl):
    dx = sph - phase
    xxx = sph - spl
    yyy = svh - svl
    vint = svh - dx*yyy/xxx
    return vint

def fittest_lightcurve(output_path,x,template_df,template_number,flag):
    template_columns_list = template_df.columns.to_list()
    template_df.set_index(['phase'],inplace=True)
    
    # expand window by adjust phase  > 0.5 - 1 ==>  range -0.5 to 0 
    A = pd.DataFrame()
    A['phase'] = template_df.index.to_list()
    # A.set_index(['phase'],inplace=True)
    A['mag'] = (x[0] * template_df[template_columns_list[template_number]].values) + x[1]
    A.loc[A['phase']>=0.5,['phase']] = A.loc[A['phase']>=0.5,['phase']] - 1 
    
    #   range 0 to 1
    B = pd.DataFrame()
    B['phase'] = template_df.index.to_list()
    # B.set_index(['phase'],inplace=True)
    B['mag'] = (x[0] * template_df[template_columns_list[template_number]].values) + x[1]
    
    # expand window by adjust phase  < 0.5 + 1  ==>  range 1 to 0.5 
    C = pd.DataFrame()
    C['phase'] = template_df.index.to_list()
    # C.set_index(['phase'],inplace=True)
    C['mag'] = (x[0] * template_df[template_columns_list[template_number]].values) + x[1]
    C.loc[C['phase']<=0.5,['phase']] = C.loc[C['phase']<=0.5,['phase']] + 1
    
    fitted_df = pd.concat([A,B,C]).sort_values(by=['phase'])
    fitted_df1 = fitted_df.copy()
    fitted_df1.set_index(['phase'],inplace=True)
    fitted_df1.to_csv(output_path+'/'+flag+'-lightcurve-fitting-result.csv')
    
    return fitted_df

def fitlc(star_df_v,star_df_i,template_df,template_number):
    mag_min_i = star_df_i['mag'].min()
    mag_min_v = star_df_v['mag'].min()
    
    mag_max_i = star_df_i['mag'].max()
    mag_max_v = star_df_v['mag'].max()
    
    median_guess_i = (mag_max_v-mag_min_v)/2
    amplitude_guess_i = (mag_max_i - mag_min_i)
    
    median_guess_v = (mag_max_v-mag_min_v)/2
    amplitude_guess_v = mag_max_v - mag_min_v
    
    x0_v = np.array([amplitude_guess_v,median_guess_v])
    
    std_v = np.std(star_df_v['mag'])
    std_i = np.std(star_df_i['mag'])

    sol_v = minimize(objective,x0_v,args=(star_df_v,template_df,template_number),method='SLSQP',bounds=((0.1,None),(0,None)))

    x0_i = np.array([sol_v.x[0]*1.6,sol_v.x[1]])
    sol_i = minimize(objective,x0_i,args=(star_df_i,template_df,template_number),method='SLSQP',bounds=((0.1,None),(0,None)))

    x0_v = sol_v.x
    x0_i = sol_i.x
    
    score = round(objective(x0_v, star_df_v, template_df, template_number) + objective(x0_i, star_df_i, template_df, template_number),4)
   
    x_optimize = [x0_v[0],x0_v[1],x0_i[0],x0_i[1]]    
    return x_optimize,score

def objective(x,star_df,template_df,template_number):
    
    col_tem_list = template_df.columns.to_list()
    squared_error = 0
    
    for i in range(1,len(star_df)):
        p_v = star_df['phase'].iloc[i]
        mag_v = star_df['mag'].iloc[i]
        if(template_df['phase']>=p_v).any() and (template_df['phase']<p_v).any():
            sp_vin = template_df.loc[template_df['phase']>=p_v].iloc[0].values[0]
            sp_vin_minus1 = template_df.loc[template_df['phase']<p_v].iloc[-1].values[0]
            
            mag_vin = template_df.loc[template_df['phase']==sp_vin,col_tem_list[template_number]].values[0]
            mag_vin_minus1 = template_df.loc[template_df['phase']==sp_vin_minus1,col_tem_list[template_number]].values[0]
            
            svh_v = (x[0]*mag_vin) + x[1]
            svl_v = (x[0]*mag_vin_minus1) + x[1]
            
            vint_v = vintrp(p_v,sp_vin,sp_vin_minus1,svh_v,svl_v)
            squared_error += round((star_df['mag'].iloc[i] - vint_v)**2,8)
            
    return squared_error

def save_star_data(output_path,filename,best_template,best_period,template_dict_inform,star_df_v,star_df_i):
    vmin = star_df_v['mag'].min()
    v_t0 = star_df_v.loc[star_df_v['mag']==vmin,'hjd'].values[0]
    imin = star_df_i['mag'].min()
    i_t0 = star_df_i.loc[star_df_i['mag']==imin,'hjd'].values[0]
    
    x_optimize_i = template_dict_inform[best_template][1][0][2],template_dict_inform[best_template][1][0][3]
    x_optimize_v = template_dict_inform[best_template][1][0][0],template_dict_inform[best_template][1][0][1]
    
    #plot fitting results with grapth to output folder
    v_cyc = (star_df_v[['hjd']] - v_t0)/best_period
    v_cyc['hjd'] = v_cyc['hjd'] - v_cyc['hjd'].astype(int)
    
    v_cyc.loc[v_cyc['hjd']<0,'hjd'] = v_cyc.loc[v_cyc['hjd']<0,'hjd'] + 1
    v_cyc.loc[v_cyc['hjd']>1,'hjd'] = v_cyc.loc[v_cyc['hjd']>1,'hjd'] - 1
    
    i_cyc = (star_df_i[['hjd']] - i_t0)/best_period
    i_cyc['hjd'] = i_cyc['hjd'] - i_cyc['hjd'].astype(int)
    
    i_cyc.loc[i_cyc['hjd']<0,'hjd'] = i_cyc.loc[i_cyc['hjd']<0,'hjd'] + 1
    i_cyc.loc[i_cyc['hjd']>1,'hjd'] = i_cyc.loc[i_cyc['hjd']>1,'hjd'] - 1
    
    star_df_v['phase'] = v_cyc
    star_df_i['phase'] = i_cyc
    
    df = pd.concat([star_df_v,star_df_i])
    df.set_index(['phase'],inplace=True)
    df.to_csv(output_path+'/'+filename.split('.')[0]+'.csv')
    
def plot_best_fitting(output_path,best_period,filename,best_template,template_dict_inform,template_df,star_df_v,star_df_i):
    #find min of magnitude
    vmin = star_df_v['mag'].min()
    v_t0 = star_df_v.loc[star_df_v['mag']==vmin,'hjd'].values[0]
    imin = star_df_i['mag'].min()
    i_t0 = star_df_i.loc[star_df_i['mag']==imin,'hjd'].values[0]
    
    
    x_optimize_i = template_dict_inform[best_template][1][0][2],template_dict_inform[best_template][1][0][3]
    x_optimize_v = template_dict_inform[best_template][1][0][0],template_dict_inform[best_template][1][0][1]
    
    #plot fitting results with grapth to output folder
    v_cyc = (star_df_v[['hjd']] - v_t0)/best_period
    v_cyc['hjd'] = v_cyc['hjd'] - v_cyc['hjd'].astype(int)
    
    v_cyc.loc[v_cyc['hjd']<0,'hjd'] = v_cyc.loc[v_cyc['hjd']<0,'hjd'] + 1
    v_cyc.loc[v_cyc['hjd']>1,'hjd'] = v_cyc.loc[v_cyc['hjd']>1,'hjd'] - 1
    
    i_cyc = (star_df_i[['hjd']] - i_t0)/best_period
    i_cyc['hjd'] = i_cyc['hjd'] - i_cyc['hjd'].astype(int)
    
    i_cyc.loc[i_cyc['hjd']<0,'hjd'] = i_cyc.loc[i_cyc['hjd']<0,'hjd'] + 1
    i_cyc.loc[i_cyc['hjd']>1,'hjd'] = i_cyc.loc[i_cyc['hjd']>1,'hjd'] - 1
    
    title = 'filename: %s, best_template: %d'%(filename,best_template)
    
    star_df_v['phase'] = v_cyc
    star_df_i['phase'] = i_cyc
    

    plt.figure()
    fit_lightcurve_df_v = fittest_lightcurve(output_path,x_optimize_v, template_df.copy(),best_template,'v')
    fit_lightcurve_df_v.set_index(['phase'],inplace=True)

    fit_lightcurve_df_i = fittest_lightcurve(output_path,x_optimize_i, template_df,best_template,'i')
    fit_lightcurve_df_i.set_index(['phase'],inplace=True)

    plt.errorbar(star_df_v['phase'],
                      star_df_v['mag'],
                      yerr=star_df_v['err'],
                      xerr=0.03,
                      fmt='o',
                      color='black',
                      capsize=5,label='data filter v')
        
    fit_lightcurve_df_v['mag'].plot(x='phase',
                                  y='mag',
                                  color='r',
                                  title=title,
                                  lw=2,label='lightcurve filter v')
    
    plt.errorbar(star_df_i['phase'],
                      star_df_i['mag'],
                      yerr=star_df_i['err'],
                      xerr=0.03,
                      fmt='o',
                      color='blue',
                      capsize=5,label='data filter i')
        
    fit_lightcurve_df_i['mag'].plot(x='phase',
                                  y='mag',
                                  color='m',
                                  lw=2,label='lightcurve filter i')
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.gca().invert_yaxis()
    plt.savefig(output_path+"/result-fitting= "+filename.split('.')[0],bbox_inches='tight',dpi=100)
    plt.savefig("output/img/"+"/result-fitting= "+filename.split('.')[0],bbox_inches='tight',dpi=100)
    
def save_fitting_result(output_path,template_dict_inform,template_df_columns):
    data_result = pd.DataFrame(columns=['template','name','period','v_amp','v_median','i_amp','i_median','score'])
    
    templates = []
    templates_name = template_df_columns[1:]
    periods = []
    v_amplitude = []
    i_amplitude = []
    v_median = []
    i_median = []
    scores = []
    for template, data in template_dict_inform.items():
        template_num = template
        period_ = data[0]
        v_amp = data[1][0][0]
        i_amp = data[1][0][2]
        v_med = data[1][0][1]
        i_med = data[1][0][3]
        score = data[1][1]
        
        templates.append(template_num)
        periods.append(period_)
        v_amplitude.append(v_amp)
        i_amplitude.append(i_amp)
        v_median.append(v_med)
        i_median.append(i_med)
        scores.append(score)
    data_result['template'] = templates
    data_result['name'] = templates_name
    data_result['period'] = periods
    data_result['v_amp'] = v_amplitude
    data_result['i_amp'] = i_amplitude
    data_result['v_median'] = v_median
    data_result['i_median'] = i_median
    data_result['score'] = scores
    
    # data_result.set_index(['template'],inplace=True)
    data_result.to_csv(output_path+'/fitting_result.csv',index=False)
    return data_result

def save_all_result(best_template,output_filename,all_output_result_file,data_result,filename):
    # filename period v_amplitude iamplitude error
    df = data_result.loc[data_result['template']==best_template]
    all_output_result_file.write("%s,%d,%f,%f,%f,%f,%f,%f,%s\n"%(filename,df['template'],df['period'],df['v_amp'],df['i_amp'],df['v_median'],df['i_median'],df['score'],df['name'].values[0]))