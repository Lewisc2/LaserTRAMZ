#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  6 07:45:41 2025

@author: ctlewis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


file_path='/Users/ctlewis/Documents/Projects/Chaxas/zircon_data/zircon_05242023/Zircon_05242023zircon05242023_LT_ready.xlsx'
df = pd.read_excel(file_path,sheet_name='Buffer')
df['Time_s'] = df['Time']/1000


bstart = 10
bend=30
tstart = 33
tend=61
arrayofdwelltimes=[0.01,0.01,0.043,0.045,0.001,0.001,0.015,0.01]

data_toapprove = df[(df['Time_s'] <= tend) & (df['Time_s'] >= bstart)  & (df['SampleLabel'] == 'Temora_-_2')].reset_index(drop=True)

def threesigoutlierremoval(data,variable,intervaltype):
    threesig = 3*np.std(data[variable])
    mean = np.mean(data[variable])
    trigger = True
    while trigger is True:
        mask = (data.loc[:,variable] < mean - threesig) | (data.loc[:,variable] > mean + threesig)
        data.loc[mask,variable] = np.nan
        data.loc[:,variable] = data.loc[:,variable].infer_objects(copy=False).interpolate(method='linear')
        threesig = 3*np.std(data[variable])
        mean = np.mean(data[variable])
        if any(mask):
            pass
        else:
            break
               
    return data[variable]

def backgroundsubtract_convert_lod(data,bstart,bend,tstart,tend,arrayofdwelltimes):
    """
    Function for getting background subtracted counts data. 
    Takes raw intensities in cps, background subtracts the intensities, then converts to counts based on the integration time
    
    Parameters
    ----------
    data : pandas dataframe
        dataframe of time resolved intensities
    bstart : float
        user selected value for the background interval start
    bend : float
        user selected value for the background interval end
    tstart : float
        user selected value for the ablation interval start
    tend : float
        user selected value for the ablation interval end
    integrationtime : float
        integration time - stripped from the input file

    Returns
    -------
    ablation_backsub_counts : pandas dataframe
        pandas dataframe of background subtracted data converted into counts

    """
    data_counts = data.loc[:,'202Hg':'238U'] * arrayofdwelltimes
    data_counts.insert(0,'Time_s',data['Time_s'])
    backgrounds = data_counts[(data_counts['Time_s']>=bstart) & (data_counts['Time_s']<=bend)]
    ablation = data_counts[(data_counts['Time_s']>=tstart) & (data_counts['Time_s']<=tend)]
    
    for analyte in backgrounds.loc[:,'202Hg':'238U'].columns:
        outlier_removed_ablation = threesigoutlierremoval(ablation,analyte,'Ablation')
        backgrounds.loc[:,analyte] = backgrounds
        ablation.loc[:,analyte] = outlier_removed_ablation
        
    meanbackgrounds = backgrounds.loc[:,'202Hg':'238U'].mean()
    sebackgrounds = backgrounds.loc[:,'202Hg':'238U'].std()/np.sqrt(len(backgrounds))
    lods = 3*backgrounds.loc[:,'202Hg':'238U'].std()
    ablation_backsub = ablation.loc[:,'202Hg':'238U'].sub(meanbackgrounds,axis='columns').clip(lower=0)
    ablation_backsub.insert(0,'Time_s',data['Time_s'])
    ablation_backsub = ablation_backsub.reset_index(drop=True)
        
    return ablation_backsub,meanbackgrounds,sebackgrounds,lods


counts_data,backgrounds,sebackgrounds,lods = backgroundsubtract_convert_lod(data_toapprove,bstart,bend,tstart,tend,arrayofdwelltimes)

isotopes_means = pd.DataFrame([counts_data.loc[:,'202Hg':'238U'].mean()]) # put means of isotopic values in a dataframe
isotopes_SE = pd.DataFrame([counts_data.loc[:,'202Hg':'238U'].sem()])
isotopes_full_SE = ((sebackgrounds)**2 + (isotopes_SE)**2)**(1/2)
isotopes_full_SE = isotopes_full_SE.add_suffix('_1SE')
isotopes_full_SE.loc[:,'204Pb_1SE'] = ((isotopes_full_SE.loc[:,'204Pb_1SE'])**2 + (isotopes_full_SE.loc[:,'202Hg_1SE'])**2)**(1/2)



dist1 = np.random.normal(100,100*0.05,size=200)
dist1_mean = np.mean(dist1)
dist2 = np.random.normal(200,100*0.05,size=200)
dist2_mean = np.mean(dist2)

covar_sum_term = []
for i in range(0,200):
    dist1_i = dist1[i]-dist1_mean
    dist2_i = dist2[i]-dist2_mean
    productterm = dist1_i*dist2_i
    covar_sum_term.append(productterm)
    
sigmaxy2 = 1/(len(dist1-1)) * np.sum(covar_sum_term)
covar = 2*sigmaxy2/(dist1_mean * dist2_mean)

(print('LOOPED COVAR'))
print(covar)


covarinline = 2 * (1/len(dist1-1) * np.sum([(dist1[i] - dist1_mean)*(dist2[i]-dist2_mean) if dist1[i]>0 and dist2[i]>0 and 1<0 else 0.0 for i in range(0,200)]) ) / (dist1_mean*dist2_mean)

(print('IN LINE COVAR'))
print(covarinline)


mu_207Pb206Pb = 0.0577339520515288
data = counts_data


summationlist = []
for i in range(0,len(data)):
    pb207term = data.loc[i,'207Pb']-float(isotopes_means['207Pb']) 
    pb206term = data.loc[i,'206Pb']-float(isotopes_means['206Pb'])
    summationterm76 = pb207term*pb206term
    summationlist.append(summationterm76)
summationtermfull = np.sum(summationlist)

covmatrix = np.cov(data['207Pb'],data['206Pb'])
pb76covar = covmatrix[0,1]
pb76covartest = np.cov(data['207Pb'],data['206Pb'])[0,1]
pb76covar==pb76covartest

pb207contribution = (isotopes_full_SE['207Pb_1SE']*np.sqrt(len(data))/isotopes_means['207Pb'])**2
pb206contribution = (isotopes_full_SE['206Pb_1SE']*np.sqrt(len(data))/isotopes_means['206Pb'])**2
unbiasedestvalue = 2*1/(len(data)-1)*summationtermfull/(isotopes_means['207Pb']*isotopes_means['206Pb'])
covarvalue = 2*pb76covar/(isotopes_means['207Pb']*isotopes_means['206Pb'])
print(pb207contribution)
print(pb206contribution)
print(unbiasedestvalue)
print(covarvalue)

print(pb207contribution+pb206contribution-unbiasedestvalue)
print(pb207contribution+pb206contribution-covarvalue)

print(mu_207Pb206Pb*np.sqrt(pb207contribution+pb206contribution-unbiasedestvalue)/np.sqrt(len(data)))
print(pb207contribution+pb206contribution-covarvalue)


summationlist = []
for i in range(0,len(data)):
    pb207term = data.loc[i,'238U']-float(isotopes_means['238U']) 
    pb206term = data.loc[i,'235U']-float(isotopes_means['235U'])
    summationterm76 = pb207term*pb206term
    summationlist.append(summationterm76)
summationtermfull = np.sum(summationlist)

covmatrix = np.cov(data['238U'],data['235U'])
pb76covar = covmatrix[0,1]
pb76covartest = np.cov(data['238U'],data['235U'])[0,1]
pb76covar==pb76covartest

pb207contribution = (isotopes_full_SE['238U_1SE']*np.sqrt(len(data))/isotopes_means['238U'])**2
pb206contribution = (isotopes_full_SE['235U_1SE']*np.sqrt(len(data))/isotopes_means['235U'])**2
unbiasedestvalue = 2*1/(len(data)-1)*summationtermfull/(isotopes_means['238U']*isotopes_means['235U'])
covarvalue = 2*pb76covar/(isotopes_means['238U']*isotopes_means['235U'])
print(pb207contribution)
print(pb206contribution)
print(unbiasedestvalue)
print(covarvalue)

print(pb207contribution+pb206contribution-unbiasedestvalue)
print(pb207contribution+pb206contribution-covarvalue)

0.000986/mu_207Pb206Pb*100

isotopes_full_SE*np.sqrt(len(data))

import numpy as np
pb_bias_dict = {'NIST-610':{'206Pb/204Pb': [17.047,0.0018/2],'207Pb/204Pb': [15.509,0.001/2],'208Pb/204Pb': [36.975,0.0026/2],'207Pb/206Pb': 0.9098},
                'NIST-612':{'206Pb/204Pb': [17.094,0.0026],'207Pb/204Pb': [15.510,0.0036],'208Pb/204Pb': [37.000,0.0094],'207Pb/206Pb': 0.9073},
                'NIST-614':{'206Pb/204Pb': [17.833,0.0134],'207Pb/204Pb': [15.533,0.0066],'208Pb/204Pb': [37.472,0.0214],'207Pb/206Pb': 0.8710}
                }


std = 'NIST-614'
pb74 = pb_bias_dict[std]['207Pb/204Pb'][0]
pb64 = pb_bias_dict[std]['206Pb/204Pb'][0]
s_pb74 = pb_bias_dict[std]['207Pb/204Pb'][1]
s_pb64 = pb_bias_dict[std]['206Pb/204Pb'][1]
pb76 = pb_bias_dict[std]['207Pb/206Pb']

pb76 * np.sqrt((s_pb74/pb74)**2 + (s_pb64/pb64)**2)


u_bias_dict = {'NIST-610':{'238U/235U': 419.4992},
               'NIST-612':{'238U/235U': 418.2650},
               'NIST-614':{'238U/235U': 374.4964}
               }



std = 'NIST-614'
ratio = u_bias_dict[std]['238U/235U']
u38 = 99.7284
u35 = 0.2663
s_u38 = 0.0003/2
s_u35 = 0.0003/2

ratio * np.sqrt((s_u38/u38)**2 + (s_u35/u35)**2)




f = 0.5
value = 15
othervalue = 100

value**f*othervalue == (value**f)*othervalue





df['1S 206Pb/238U Bias Corrected Age'] = np.sqrt((1/(lambda_238*(df['206Pb/238U Bias Corrected']+1)))**2*(df['1S 206Pb/238U Bias Corrected'])**2 + 
                                                 (-(np.log(df['206Pb/238U Bias Corrected']+1))/(lambda_238**2))**2*(lambda_238_2sig_percent/2/100)**2)
                
                
                
# define constants for reducing U-Pb data
lambda_238 = 1.55125e-10 # Jaffey et al. 1971
lambda_235 = 9.8485e-10 # Jaffey et al. 1971
lambda_232 = 4.9475e-11 # Le Roux and Glendenin 1963 / Steiger and Jager 1977
lambda_230 = 9.158e-6 # Cheng et al. 2000
# Errors on uraniumn decay constants from Mattionson(1987)
lambda_238_2sig_percent = 0.16 
lambda_235_2sig_percent = 0.21
lambda_230_2sig_percent = 0.30

SK74_2sig = 0.3
SK64_2sig = 1 
                
df['1S Concordant Age'] = np.sqrt((1/(lambda_238*(df['Concordant 206Pb/238U']+1)))**2*(df['1S Concordant 206Pb/238U'])**2 + 
                                  (-np.log(df['Concordant 206Pb/238U']+1)/(lambda_238**2))**2*(lambda_238_2sig_percent/2/100*lambda_238)**2)


common76 = 0.847332
scommon76 = common76*np.sqrt(((SK74_2sig/2)/15.616423)**2 + ((SK64_2sig/2)/18.430111)**2)
r638 = 0.007993
s638 = 0.000424
r76 = 0.03582
s76 = 0.034867
timsage = 416780000
stimsage = 330000
r638concordant = 0.008105


s638concordant = r638*np.sqrt((s638/r638)**2 + (scommon76/common76)**2)
s638concordant

np.sqrt((1/(lambda_238*(r638concordant+1)))**2*(s638concordant)**2)










