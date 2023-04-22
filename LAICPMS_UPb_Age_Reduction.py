#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 15:15:39 2023

@author: ctlewis
"""

# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('agg')
from matplotlib.figure import Figure
from matplotlib.patches import Ellipse
import panel as pn
import statistics
# import math
import param
# import io
# import os
import sys
import statsmodels.api as sm
# from statsmodels.graphics.api import abline_plot
from patsy import dmatrices
# from panel.viewable import Viewer
import scipy
from scipy import stats
# import base64

#%%
# define constants for reducing U-Pb data
lambda_238 = 1.55125e-10 #
lambda_235 = 9.8485e-10 #
lambda_232 = 4.9475e-11 # 
lambda_230 = 9.158e-6

#%%
class calc_fncs:
    """ Class that holds all of the functions for reducing the reduced LAICPMS data"""
    def __init__(self,*args):
        for a in args:
            self.__setattr__(str(a), args[0])

    # this needs to be a function that receives the selected sample
    def get_data_TW_regressions(df):
        """
        """
        pts_x = np.array([df['U238_Pb206'],np.zeros_like(df['SK207_206'])]).T
        pts_y = np.array([df['207Pb/206Pb'],df['SK207_206']]).T
        discordia_t = np.zeros((len(df),2))
    
        for i in range(0,len(df)):
            discordia_t[i] = np.poly1d(np.polyfit(pts_x[i],pts_y[i],1))
    
        return discordia_t
    
    
    
    def get_TW_concordia():
        """
        """
        # calculate T-W discordia - put in plot function
        t = np.linspace(1,4.6e9,100000)
        u238_pb206 = np.zeros(len(t))
        pb207_206r = np.zeros(len(t))
        for i in range(0,len(t)):
            u238_pb206[i] = 1/(np.exp(lambda_238*t[i])-1)
            pb207_206r[i] = (1/137.88) * ((np.exp(lambda_235*t[i])-1) / (np.exp(lambda_238*t[i])-1))
            
        return u238_pb206,pb207_206r
    
    
    
    
    def get_projections(df,ellipse_mode_selector,power):
        """
        """
        # get the TW concordia values
        x_TW,y_TW = calc_fncs.get_TW_concordia()
        discorida_regressions = calc_fncs.get_data_TW_regressions(df) # get regressions
        # array of xvalues to project over
        x_vals = np.linspace(min(x_TW),max(x_TW),100000)
        # pts_pb_r = [] # set up array to be filled with calculated radiogenic lead component
        pts_pb_r = np.zeros(len(discorida_regressions))
        concordia_238_206 = np.zeros(len(discorida_regressions))
    
        for i in range(0,len(discorida_regressions)):
    
            discordia_207206 = discorida_regressions[i][0]*x_vals+discorida_regressions[i][1]
            discordia_238206 = (discordia_207206-discorida_regressions[i][1])/discorida_regressions[i][0]
    
            delta_y = (discorida_regressions[i][1] + x_TW * discorida_regressions[i][0]) - y_TW # distance of y value from line
            indx = np.where(delta_y[1:]*delta_y[:-1]<0)[0] # index the regression for where the curve cross the regression
            d_ratio = delta_y[indx] / (delta_y[indx] - delta_y[indx + 1]) # similar triangles geometry gets points
            points = np.zeros((len(indx), 2)) # empty array for crossing points
            points[:,0] = x_TW[indx] + d_ratio * (x_TW[indx+1] - x_TW[indx]) # x crossings
            points[:,1] = y_TW[indx] + d_ratio * (y_TW[indx+1] - y_TW[indx]) # y crossings
            y_point = y_TW[indx] + d_ratio * (y_TW[indx+1] - y_TW[indx]) # y crossings
            x_point = x_TW[indx] + d_ratio * (x_TW[indx+1] - x_TW[indx]) # x crossings
            
            if len(y_point) >= 1:
                pts_pb_r[i] = min(y_point)
            elif len(y_point) < 0:
                pts_pb_r[i] = 0
                
            
            if len(x_point) > 1:
                concordia_238_206[i] = min(x_point)
            elif len(x_point) == 1:
                concordia_238_206[i] = x_point
        
        # if ellipse_mode_selector == 'Point Estimates':
        # points = points
        # concordia_238_206 = concordia_238_206
        # pts_pb_r = pts_pb_r
        # discordia_207206 = discordia_207206
        # discordia_238206 = discordia_238206
        return points,concordia_238_206,pts_pb_r,discordia_207206,discordia_238206
        # elif ellipse_mode_selector == 'Ellipses':
        #     x2 = 1/df['U238_Pb206']
        #     y2 = df['Pb207_Pb206']
            
        #     cov2 = np.cov(x2,y2)
        #     eigval2,eigvec2 = np.linalg.eig(cov2)
        #     order2 = eigval2.argsort()[::-1]
        #     eigvals_order2 = eigval2[order2]
        #     eigvecs_order2 = eigvec2[:,order2]
            
        #     c2 = (np.mean(x2),np.mean(y2))
        #     wid2 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power),df=2)*eigvals_order2[0])
        #     hgt2 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power),df=2)*eigvals_order2[1])
        #     theta2 = np.degrees(np.arctan2(*eigvecs_order2[:,0][::-1]))
            
        #     ell2_params = [c2,wid2,hgt2,theta2]
            
        #     return ell2_params
    
    def get_ellipse(df,power):
        callingmethod = sys._getframe().f_back.f_code.co_name
        if callingmethod == 'correct_sample_ages':
            df.replace([np.inf,-np.inf],np.nan,inplace=True)
            df.dropna(inplace=True)
            x2 = 1/df['206Pb/238U Corrected']
        else:
            x2 = df['U238_Pb206']
        y2 = df['Pb207_Pb206']
        
        cov2 = np.cov(x2,y2)
        eigval2,eigvec2 = np.linalg.eig(cov2)
        order2 = eigval2.argsort()[::-1]
        eigvals_order2 = eigval2[order2]
        eigvecs_order2 = eigvec2[:,order2]
        
        c2 = (np.mean(x2),np.mean(y2))
        wid2 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power),df=2)*eigvals_order2[0])
        hgt2 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power),df=2)*eigvals_order2[1])
        theta2 = np.degrees(np.arctan2(*eigvecs_order2[:,0][::-1]))
        
        ell2_params = [c2,wid2,hgt2,theta2]
        
        return ell2_params
    
    def plot_TW(df,xlim_start,xlim_stop,ylim_start,ylim_stop,label_toggle,ellipse_mode_selector,power):
        """
        """
        plt.style.use('seaborn-colorblind')
        fig = Figure(figsize=(10,14))
        # fig = Figure()
        ax = fig.add_subplot()
        
        x_TW,y_TW = calc_fncs.get_TW_concordia()
        
        disc_regressions = calc_fncs.get_data_TW_regressions(df)
        
        x_vals = np.linspace(min(x_TW),max(x_TW),100000)
        pts_pb_r = np.zeros(len(disc_regressions))
        
        for i in range(0,len(disc_regressions)):
            discordia_207206 = disc_regressions[i][0]*x_vals+disc_regressions[i][1]
            discordia_238206 = (discordia_207206-disc_regressions[i][1])/disc_regressions[i][0]
    
            delta_y = (disc_regressions[i][1] + x_TW * disc_regressions[i][0]) - y_TW # distance of y value from line
            indx = np.where(delta_y[1:]*delta_y[:-1]<0)[0] # index the regression for where the curve cross the regression
            d_ratio = delta_y[indx] / (delta_y[indx] - delta_y[indx + 1]) # similar triangles geometry gets points
            points = np.zeros((len(indx), 2)) # empty array for crossing points
            points[:,0] = x_TW[indx] + d_ratio * (x_TW[indx+1] - x_TW[indx]) # x crossings
            points[:,1] = y_TW[indx] + d_ratio * (y_TW[indx+1] - y_TW[indx]) # y crossings
            y_point = y_TW[indx] + d_ratio * (y_TW[indx+1] - y_TW[indx]) # y crossings
            if len(y_point) >= 1:
                pts_pb_r[i] = min(y_point)
            elif len(y_point) < 0:
                pts_pb_r[i] = 0
            
            if ellipse_mode_selector == 'Point Estimates':
                ax.plot(discordia_238206,discordia_207206,'-k',lw=0.5)
                ax.plot(points[:,0],points[:,1],'og',lw=1,mec='k')
            else:
                pass
        if ellipse_mode_selector == 'Point Estimates':
            ax.plot(df.U238_Pb206,df.Pb207_Pb206,'kd')
            ax.errorbar(df.U238_Pb206,df.Pb207_Pb206,xerr=df['206/238 Reg. err'],yerr=df['SE 207/206'],fmt='none',ecolor='k',elinewidth=0.5)
            ax.plot(x_TW,y_TW,'k',lw=1)
        elif ellipse_mode_selector == 'Ellipses':
            for i in df.SampleLabel.unique():
                conf_ellipse = calc_fncs.get_ellipse(df[df['SampleLabel']==i],power)
                ell = Ellipse(conf_ellipse[0],conf_ellipse[1],conf_ellipse[2],conf_ellipse[3],color='darkslategray',alpha=0.3,ec='k')
                ax.add_artist(ell)
            ax.plot(x_TW,y_TW,'k',lw=1)
        ax.set_xlabel('$^{238}$U / $^{206}$Pb$^{*}$',fontsize=12)
        ax.set_ylabel('[$^{207}$Pb / $^{206}$Pb]$^{*}$',fontsize=12)
        ax.tick_params(axis='both',labelsize=10)
        ax.set_ylim(ylim_start,ylim_stop)
        ax.set_xlim(xlim_start,xlim_stop)
        
        # zip joins x and y coordinates in pairs
        if ellipse_mode_selector == 'Point Estimates':
            if 'Concordia' in label_toggle:
                for x,y,t in zip(df.U238_Pb206,df.Pb207_Pb206,df.SampleLabel):
        
                    label = t
        
                    ax.annotate(label, # this is the text
                                 (x,y), # these are the coordinates to position the label
                                 textcoords="offset points", # how to position the text
                                 xytext=(0,10), # distance from text to points (x,y)
                                 ha='center',
                                fontsize=8) # horizontal alignment can be left, right or center
        else:
            pass
    
        return fig
    
    
    
    
    def plot_boxplot(ages,analysis_ID,label_toggle,ellipse_mode_selector):
        """
        """
        if ellipse_mode_selector == 'Point Estimates':
            plt.style.use('seaborn-colorblind')
            fig = Figure(figsize=(2,4))
            ax = fig.add_subplot()
        
            bp = ax.boxplot(ages, patch_artist=True, boxprops=dict(facecolor='slategray',color='k'), 
                       medianprops=dict(color='limegreen'), meanprops=dict(marker='d',mfc='limegreen',mec='k',markersize=4),
                            flierprops=dict(marker='o',mfc='None',mec='k',markersize=4),
                       showmeans=True)
        
        
            ax.text(0.05,0.8,'Mean ='+str(round(ages.mean(),2)),fontsize=5,transform=ax.transAxes)
            ax.text(0.05,0.7,'Med ='+str(round(ages.median(),2)),fontsize=5,transform=ax.transAxes)
            ax.text(0.05,0.6,'Min ='+str(round(ages.min(),2)),fontsize=5,transform=ax.transAxes)
            ax.text(0.05,0.5,'Max ='+str(round(ages.max(),2)),fontsize=5,transform=ax.transAxes)
            ax.text(0.05,0.4,'n = '+str(len(ages)),fontsize=5,transform=ax.transAxes)
            
            ax.set_ylabel('Age (Ma)',fontsize=8)
            ax.set_xlabel(' ',fontsize=1)
            ax.tick_params(axis='both',labelsize=5)
            
            if 'Box + Whisker' in label_toggle:
                fliers = [item.get_ydata() for item in bp['fliers']]
                for x,t in zip(ages,analysis_ID):
                    
                    if t in fliers:
                        label = t
                        age = x
                        ax.annotate(label,age,textcoords='offset points',ha='center',fontsize=5)
                
            return fig
        else:
            return print('In Ellipse Mode')
    
    
    
    
    def correct_standard_ages(df,std_txt,regression_selector,ellipse_mode_selector,power,Pb_Th_std_crct_selector):
        """
        function to calculate Pb and fractionation factor corrected ages
        """
        # create a dictionary that holds known or estimated U/Th ratios of zircon and associated magma for standards, as well as common Pb ratios
        stds_dict = {'Temora': [2.4,0.79200,18.0528,15.5941], # Black et al., 2004. Ambiguous Th correction > assumed D = 0.33
                     'FishCanyon': [1.496558,0.454545,18.4275,15.5425], # Schmitz and Bowring 2001; only one with measured common Pb so far
                     '94-35': [1,0.33,18.6191,15.626],
                     'Plesovice': [10.7,0.25,18.1804,15.6022], # Slama et al 2008
                     'R33': [1.4,0.46200,18.0487,15.5939], # Black et al., 2004. Ambiguous Th correction > assumed D = 0.33
                     '91500': [1,0.33,16.9583,15.4995], # Wiedenbeck et al 1995
                     'FC1': [1.7,0.56100,16.892,15.492], # Paces and Miller 1993. Ambiguous Th correction > assumed D = 0.33
                     'Oracle': [2.2,0.725999,16.2726,15.4099], # unpublished; Bowring > assumed D = 0.33
                     'Tan-Bra': [1.2,0.39600,14.0716,14.8653], # Pecha unpublished > assumed D = 0.33
                     'OG1': [1.3,0.42900,11.8337,13.6071]} # Stern et al 2009. Ambiguous Th correction > assumed D = 0.33
        # get points needed for correction
        points,concordia_238_206,pts_pb_r,na,naII = calc_fncs.get_projections(df,ellipse_mode_selector,power)
        
        # set up empty arrays to be filled for common lead corrections
        common_filter = []
        
        pb_m = df['207Pb/206Pb'] # measured 207/206
        if Pb_Th_std_crct_selector == 'Common Pb':
            common = df['SK207_206']
        elif Pb_Th_std_crct_selector == 'Common Pb + Th Disequil.':
            df['SK207_206'] = stds_dict.get(std_txt)[3]/stds_dict.get(std_txt)[2] # calculated Stacey-Kramers 207/206 overwriting the input 207/206 should manual values be requested
            common = df['SK207_206']
            
        
        for i in common:
            if i <= 0:
                common_filter.append(0)
            else:
                common_filter.append(i)
        
        # calculate fraction of common Pb
        f_ = (pb_m - pts_pb_r) / (common - pts_pb_r)
        # set up array to set f = 0 if point lies on or below Concordia (i.e., no common Pb present)
        f = []
        
        for k,j in zip(common_filter,f_):
            if k <= 0:
                f.append(0)
            elif j < 0:
                f.append(0)
            else:
                f.append(j)
                
        # append the calculations to the sample dataframe
        df['207Pb/206Pbr'] = pts_pb_r
        df['f'] = f
        
        df['counts_pb206r'] = df['206Pb'] * (1-df['f'])
        df['206Pb/238Upbc_numerical'] = 1/df['U238_Pb206']-(1/df['U238_Pb206']*f)
        df['206Pb/238Upbc'] = 1/concordia_238_206
        df['206Pb/238Uc_age'] = np.log(df['206Pb/238Upbc'] + 1) / lambda_238
        UThstd,UThstd_rx  = stds_dict.get(std_txt)[0],stds_dict.get(std_txt)[1]
        DThU = (1/UThstd)/(1/UThstd_rx)
        df['206Pb/238UPbThc'] = df['206Pb/238Upbc'] - (lambda_238/lambda_230*(DThU-1))
        df['206Pb/238UPbTh_age'] = np.log(df['206Pb/238UPbThc'] + 1) / lambda_238
        UTh_std = stds_dict.get(std_txt)[0]
        UTh_std_m = df['238U'].mean()/df['232Th'].mean()
        
        
        avg_std_age = df['206Pb/238Uc_age'].mean()
        avg_std_age_Thcrct = df['206Pb/238UPbTh_age'].mean()
        
        if ellipse_mode_selector == 'Point Estiamtes':
            if regression_selector == '1st Order':
                avg_reg_err = df['SE 206/238 1st Order'].mean()
            elif regression_selector == '2nd Order':
                avg_reg_err = df['SE 206/238 2nd Order'].mean()
            elif regression_selector == 'Exp.':
                avg_reg_err = df['SE 206/238 Exp'].mean()
        else:
            avg_reg_err = 1
        
        return avg_std_age,avg_std_age_Thcrct,avg_reg_err,UTh_std,UTh_std_m
        
        
    
    def correct_sample_ages(df,std,std_txt,ThU_zrn,ThU_magma,Pb_Th_std_crct_selector,regression_selector,ellipse_mode_selector,power,UTh_std_norm):
        """
        function to calculate Pb and fractionation factor corrected ages
        """
        
        if ellipse_mode_selector == 'Point Estimates':
            # df = df.reset_index()
            # get points needed for correction
            points,concordia_238_206,pts_pb_r,n,a = calc_fncs.get_projections(df,ellipse_mode_selector,power)
            
            # set up empty arrays to be filled for common lead corrections
            common_filter = []
            
            pb_m = df['207Pb/206Pb'] # measured 207/206
            common = df['SK207_206']
            
            for i in common:
                if i <= 0:
                    common_filter.append(0)
                else:
                    common_filter.append(i)
            
            # calculate fraction of common Pb
            f_ = (pb_m - pts_pb_r) / (common - pts_pb_r)
            # set up array to set f = 0 if point lies on or below Concordia (i.e., no common Pb present)
            f = []
            
            for k,j in zip(common_filter,f_):
                if k <= 0:
                    f.append(0)
                elif j < 0:
                    f.append(0)
                else:
                    f.append(j)
                    
            # append the calculations to the sample dataframe
            df['207Pb/206Pbr'] = pts_pb_r
            # df['207Pb/206Pbcommon'] = common_filter
            df['f'] = f
            
            df['counts_pb206r'] = df['206Pb'] * (1-df['f'])
            df['206Pb/238Upbc_numerical'] = 1/df['U238_Pb206']-(1/df['U238_Pb206']*f)
            df['206Pb/238Upbc'] = 1/concordia_238_206
        
            
            # will need input for which standard to use now! standard will need to be calculated by a separate function and called/input here.
            frac_factor,tims_age,tims_error,std_avg,avg_std_age_Thcrct,std_avgerr,UTh_std,UTh_std_m = calc_fncs.get_standard_fracfctr(std,std_txt,Pb_Th_std_crct_selector,regression_selector,
                                                                                                                    ellipse_mode_selector,power)
            
            if Pb_Th_std_crct_selector == 'Common Pb':
                df['206Pb/238Uc_age'] = np.log(df['206Pb/238Upbc'] + 1) / lambda_238
                df['206Pb/238U_correctedage'] = df['206Pb/238Uc_age']*frac_factor
                #propagate errors
                #error on fractionation factor. includes error from ID-TIMS and ICPMS
                dfrac = np.abs(frac_factor)*((1/tims_age)**2*(tims_error/2)**2 + (1/std_avg)**2*std_avgerr**2)**(1/2)
                #error on age equation. error on decay constant includes 1.5* counting stats (Mattinson 1987)
                dage = np.abs(df['206Pb/238Uc_age'])*((1/(1/df['U238_Pb206']))**2*df['206/238 Reg. err']**2 + (0.16/2/100)**2)**(1/2)
                #error on estimation for common lead using 207 method. Uses conservaitve estimates of 1.0 for 206/204 and 0.3 for 207/204 (Mattionson, 1987)
                dinit = np.abs(common)*((1/df['SK207_204'])**2*(0.3/2)**2 + (1/df['SK206_204'])**2*(1/2)**2)**(1/2)
                # total propagated error
                # dagetot = (dage**2 +  df['SE 206/204']**2 + dinit**2)**(1/2)
                dagetot = np.abs(df['206Pb/238Uc_age'])*((1/tims_age)**2*(tims_error/2)**2 + (1/std_avg)**2*std_avgerr**2 + 
                                                         (1/(1/df['U238_Pb206']))**2*df['206/238 Reg. err']**2 + (0.16/2/100)**2 +
                                                         (1/df['SK207_204'])**2*(0.3/2)**2 + (1/df['SK206_204'])**2*(1/2)**2)**(1/2)
                
                df['∆206/238 age (meas.)'] = dage
                df['∆206/238 age (tot.)'] = dagetot
            elif Pb_Th_std_crct_selector == 'Common Pb + Th Disequil.':
                if UTh_std_norm == 'Off':
                    df['206Pb/238UPbThc'] = df['206Pb/238Upbc'] - (lambda_238/lambda_230*((ThU_zrn/ThU_magma)-1))
                    df['206Pb/238Uc_age'] = (np.log(df['206Pb/238UPbThc'] + 1) / lambda_238)
                    df['206Pb/238U_correctedage'] = df['206Pb/238Uc_age']*frac_factor
                    #propagate errors
                    #error on fractionation factor. includes error from ID-TIMS and ICPMS
                    dfrac = np.abs(frac_factor)*((1/tims_age)**2*(tims_error/2)**2 + (1/avg_std_age_Thcrct)**2*std_avgerr**2)**(1/2)
                    #error on age equation. error on decay constant includes 1.5* counting stats (Mattinson 1987)
                    dage = np.abs(df['206Pb/238U_correctedage'])*((1/(1/df['U238_Pb206']))**2*df['206/238 Reg. err']**2 + (0.16/2/100)**2)**(1/2)
                    #error on estimation for common lead using 207 method. Uses conservaitve estimates of 1.0 for 206/204 and 0.3 for 207/204 (Mattionson, 1987)
                    dinit = np.abs(common)*((1/df['SK207_204'])**2*(0.3/2)**2 + (1/df['SK206_204'])**2*(1/2)**2)**(1/2)
                    # UTh errors - only using errors from measured zircon here, as this will by and large be the largest error contribution
                    # should probably add error for U / Th in icpms glass analyses, though rock measurements will undoubtedly be incorrectly used by users of the
                    # program (in absence of other data) and so added error would overall be fairly misleading anyway
                    # errors on absolute concentrations from ID-TIMS are overall negligible and often not reported either. 
                    # need to put in possibility of putting in errors on Th/U measurements
                    # total propagated error
                    # dagetot = (dage**2 +  df['SE 206/204']**2 + dinit**2)**(1/2)
                    dagetot = np.abs(df['206Pb/238Uc_age'])*((1/tims_age)**2*(tims_error/2)**2 + (1/std_avg)**2*std_avgerr**2 + 
                                                             (1/(1/df['U238_Pb206']))**2*df['206/238 Reg. err']**2 + (0.16/2/100)**2 +
                                                             (1/df['SK207_204'])**2*(0.3/2)**2 + (1/df['SK206_204'])**2*(1/2)**2)**(1/2)
                    df['∆206/238 age (meas.)'] = dage
                    df['∆206/238 age (tot.)'] = dagetot
                elif UTh_std_norm == 'Calc U/Th from Std.':
                    I232_std = std['232Th'].mean()
                    I238_std = std['238U'].mean()
                    std232err = std['232Th_1SE'].mean()
                    std238err = std['238U_1SE'].mean()
                    ThUzrn_calc = (((1/UTh_std)/(1/UTh_std_m))/(I232_std/I238_std)) * (df['232Th']/df['238U'])
                    df['206Pb/238UPbThc'] = df['206Pb/238Upbc'] - (lambda_238/lambda_230*((ThUzrn_calc/ThU_magma)-1))
                    df['206Pb/238Uc_age'] = (np.log(df['206Pb/238UPbThc'] + 1) / lambda_238)
                    df['206Pb/238U_correctedage'] = df['206Pb/238Uc_age']*frac_factor
                    #propagate errors
                    #error on fractionation factor. includes error from ID-TIMS and ICPMS
                    dfrac = np.abs(frac_factor)*((1/tims_age)**2*(tims_error/2)**2 + (1/avg_std_age_Thcrct)**2*std_avgerr**2)**(1/2)
                    #error on age equation. error on decay constant includes 1.5* counting stats (Mattinson 1987)
                    dage = np.abs(df['206Pb/238U_correctedage'])*((1/(1/df['U238_Pb206']))**2*df['206/238 Reg. err']**2 + (0.16/2/100)**2)**(1/2)
                    #error on estimation for common lead using 207 method. Uses conservaitve estimates of 1.0 for 206/204 and 0.3 for 207/204 (Mattionson, 1987)
                    dinit = np.abs(common)*((1/df['SK207_204'])**2*(0.3/2)**2 + (1/df['SK206_204'])**2*(1/2)**2)**(1/2)
                    # UTh errors
                    dThU = np.abs(ThUzrn_calc) * ((df['232Th_1SE']/df['232Th'])**2 + (df['238U_1SE']/df['238U'])**2 + (std232err/I232_std)**2 + (std238err/I238_std)**2)**(1/2)
                    # total propagated error
                    # dagetot = (dage**2 +  df['SE 206/204']**2 + dinit**2 + dThU**2)**(1/2)
                    dagetot = np.abs(df['206Pb/238Uc_age'])*((1/tims_age)**2*(tims_error/2)**2 + (1/std_avg)**2*std_avgerr**2 + 
                                                             (1/(1/df['U238_Pb206']))**2*df['206/238 Reg. err']**2 + (0.16/2/100)**2 +
                                                             (1/df['SK207_204'])**2*(0.3/2)**2 + (1/df['SK206_204'])**2*(1/2)**2 +
                                                             (df['232Th_1SE']/df['232Th'])**2 + (df['238U_1SE']/df['238U'])**2 + (std232err/I232_std)**2 + (std238err/I238_std)**2)**(1/2)
                    df['∆206/238 age (meas.)'] = dage
                    df['∆206/238 age (tot.)'] = dagetot
            
            return df
                
        elif ellipse_mode_selector == 'Ellipses':
            # ellipse_params = pd.DataFrame([np.zeros(5)],columns=['206Pb/238U Center','207Pb/206Pb Center','Ell. Width','Ell. Height','Ell. Rotation'])
            frac_factor,tims_age,tims_error,std_avg,avg_std_age_Thcrct,std_avgerr,UTh_std,UTh_std_m = calc_fncs.get_standard_fracfctr(std,std_txt,Pb_Th_std_crct_selector,regression_selector,
                                                                                                                    ellipse_mode_selector,power)
        # for i in df.SampleLabel.unique():
            # conf_ellipse = calc_fncs.get_ellipse(df[df['SampleLabel']==i],power)
            
            
            # ell = Ellipse(conf_ellipse[0],conf_ellipse[1],conf_ellipse[2],conf_ellipse[3],color='darkslategray',alpha=0.3,ec='k')
            # ax.add_artist(ell)
            
            # c2 = (np.mean(x2),np.mean(y2))
            # wid2 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power),df=2)*eigvals_order2[0])
            # hgt2 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power),df=2)*eigvals_order2[1])
            # theta2 = np.degrees(np.arctan2(*eigvecs_order2[:,0][::-1]))
            # ell2_params = [c2,wid2,hgt2,theta2]
            
            points,concordia_238_206,pts_pb_r,n,a = calc_fncs.get_projections(df,ellipse_mode_selector,power)
            
            # set up empty arrays to be filled for common lead corrections
            common_filter = []
            
            pb_m = df['207Pb/206Pb'] # measured 207/206
            common = df['SK207_206']
            
            for i in common:
                if i <= 0:
                    common_filter.append(0)
                else:
                    common_filter.append(i)
            
            # calculate fraction of common Pb
            f_ = (pb_m - pts_pb_r) / (common - pts_pb_r)
            # set up array to set f = 0 if point lies on or below Concordia (i.e., no common Pb present)
            f = []
            
            for k,j in zip(common_filter,f_):
                if k <= 0:
                    f.append(0)
                elif j < 0:
                    f.append(0)
                else:
                    f.append(j)
                    
            df['207Pb/206Pbr'] = pts_pb_r
            # df['207Pb/206Pbcommon'] = common_filter
            df['f'] = f
            
            df['counts_pb206r'] = df['206Pb'] * (1-df['f'])
            df['206Pb/238Upbc_numerical'] = 1/df['U238_Pb206']-(1/df['U238_Pb206']*f)
            df['206Pb/238Upbc'] = 1/concordia_238_206
        
            
            # will need input for which standard to use now! standard will need to be calculated by a separate function and called/input here
            
            if Pb_Th_std_crct_selector == 'Common Pb':
                df['206Pb/238Uc_age'] = np.log(df['206Pb/238Upbc'] + 1) / lambda_238
                df['206Pb/238U_correctedage'] = df['206Pb/238Uc_age']*frac_factor
                
                df['206Pb/238U Corrected'] = np.exp(df['206Pb/238U_correctedage']*lambda_238) - 1
                
                ellparams_i = calc_fncs.get_ellipse(df, power)
                # ellipse_params = ellipse_params.append(pd.DataFrame([ellparams_i],columns=['Ellipse Center','Ell. Width','Ell. Height','Ell. Rotation']),ignore_index=True)
                ellparams_i = pd.DataFrame([ellparams_i],columns=['Ellipse Center','Ell. Width','Ell. Height','Ell. Rotation'])
            
            elif Pb_Th_std_crct_selector == 'Common Pb + Th Disequil.':
                if UTh_std_norm == 'Off':
                    df['206Pb/238UPbThc'] = df['206Pb/238Upbc'] - (lambda_238/lambda_230*((ThU_zrn/ThU_magma)-1))
                    df['206Pb/238Uc_age'] = (np.log(df['206Pb/238UPbThc'] + 1) / lambda_238)
                    df['206Pb/238U_correctedage'] = df['206Pb/238Uc_age']*frac_factor
                    
                    df['206Pb/238U Corrected'] = np.exp(df['206Pb/238U_correctedage']*lambda_238) - 1
                    
                    ellparams_i = calc_fncs.get_ellipse(df, power)
                    # ellipse_params = ellipse_params.append(pd.DataFrame([ellparams_i],columns=['Ellipse Center','Ell. Width','Ell. Height','Ell. Rotation']),ignore_index=True)
                    ellparams_i = pd.DataFrame([ellparams_i],columns=['Ellipse Center','Ell. Width','Ell. Height','Ell. Rotation'])
                    
                elif UTh_std_norm == 'Calc U/Th from Std.':
                    I232_std = std['232Th'].mean()
                    I238_std = std['238U'].mean()
                    std232err = std['232Th_1SE'].mean()
                    std238err = std['238U_1SE'].mean()
                    ThUzrn_calc = (((1/UTh_std)/(1/UTh_std_m))/(I232_std/I238_std)) * (df['232Th']/df['238U'])
                    df['206Pb/238UPbThc'] = df['206Pb/238Upbc'] - (lambda_238/lambda_230*((ThUzrn_calc/ThU_magma)-1))
                    df['206Pb/238Uc_age'] = (np.log(df['206Pb/238UPbThc'] + 1) / lambda_238)
                    df['206Pb/238U_correctedage'] = df['206Pb/238Uc_age']*frac_factor
                    
                    df['206Pb/238U Corrected'] = np.exp(df['206Pb/238U_correctedage']*lambda_238) - 1
                    
                    ellparams_i = calc_fncs.get_ellipse(df, power)
                    ellparams_i = pd.DataFrame([ellparams_i],columns=['Ellipse Center','Ell. Width','Ell. Height','Ell. Rotation'])
        return ellparams_i
                
                
                    
    
    
    
    
    def get_standard_fracfctr(std,std_txt,Pb_Th_std_crct_selector,regression_selector,ellipse_mode_selector,power):
        """
        """
        # selected_std = samples.get_group(std)
        
        accepted_ages = {
            'Temora':416780000,
            'FishCanyon': 28478000,
            '94-35': 55500000,
            'Plesovice': 337100000,
            'R33': 419300000,
            '91500': 1062400000,
            'FC1': 1099500000,
            'Oracle': 1436200000,
            'Tan-Bra': 2507800000,
            'OG1': 3440700000
        }
        
        TIMS_errors = {
            'Temora':330000,
            'FishCanyon': 24000,
            '94-35': 80000,
            'Plesovice': 200000,
            'R33': 400000,
            '91500': 1900000,
            'FC1': 330000,
            'Oracle': 1300000,
            'Tan-Bra': 1500000,
            'OG1': 3200000
        }
        
        if Pb_Th_std_crct_selector == 'Common Pb':
            std_avg,avg_std_age_Thcrct,std_avgerr,UTh_std,UTh_std_m = calc_fncs.correct_standard_ages(std,std_txt,regression_selector,ellipse_mode_selector,power,Pb_Th_std_crct_selector)
            frac_factor = accepted_ages.get(std_txt)/std_avg
            tims_age = accepted_ages.get(std_txt)
            tims_error = TIMS_errors.get(std_txt)
        elif Pb_Th_std_crct_selector == 'Common Pb + Th Disequil.':
            std_avg,avg_std_age_Thcrct,std_avgerr,UTh_std,UTh_std_m = calc_fncs.correct_standard_ages(std,std_txt,regression_selector,ellipse_mode_selector,power,Pb_Th_std_crct_selector)
            frac_factor = accepted_ages.get(std_txt)/avg_std_age_Thcrct
            tims_age = accepted_ages.get(std_txt)
            tims_error = TIMS_errors.get(std_txt)

        return frac_factor,tims_age,tims_error,std_avg,avg_std_age_Thcrct,std_avgerr,UTh_std,UTh_std_m
            
    
#%%
class finalize_ages(param.Parameterized):
    
    """ class that parameterizes inputs and sends them to the above functions to be rendered in a GUI"""
    file_path = param.String(default='Insert File Path')
    input_data = param.DataFrame(precedence=-1)
    output_data = param.DataFrame(precedence=-1)
    file_path_ellipse = param.String(default='Insert Ellipse File Path')
    input_data_ellipse = param.DataFrame(precedence=-1)
    output_data_ellipse = param.DataFrame(precedence=-1)
    regression_selector = param.Selector(objects=['1st Order','2nd Order','Exp.'])
    ellipse_mode_selector = param.Selector(default='Point Estimates',objects=['Point Estimates','Ellipses'])
    
    x_axis_TW_min = param.Number(default=0.5)
    x_axis_TW_max = param.Number(default=25)
    y_axis_TW = param.Range(default=(0.0,0.5),bounds=(0.0,1.0)) # for y-axis
    text_sample_selector = param.String(default='Input Sample ID')# for getting samples. need to split to get all of given name... or have them do text input
    text_standard_selector = param.String(default='Input Standard ID')# for getting samples. need to split to get all of given name... or have them do text input
    common_206204_input = param.Number()
    common_207204_input = param.Number()
    ThU_zrn_input = param.Number()
    ThU_magma_input = param.Number()
    UTh_std_norm = param.Selector(default='Off',objects=['Calc U/Th from Std.','Off'])
    Pb_Th_std_crct_selector = param.Selector(objects=['Common Pb','Common Pb + Th Disequil.'])
    power = param.Number(default=0.05)
    
    update_output_button = param.Action(lambda x: x.add_output_data(),label='Approve Data')
    export_data_button = param.Action(lambda x: x.export_data(),label='DDDT!')
    export_TWplot_button = param.Action(lambda x: x.export_plot(),label='Save Plot')
    label_toggle = param.ListSelector(default=['Concordia'], objects=['Concordia','Box + Whisker'])
    
    def __init__(self,**params):
        super().__init__(**params)
        self.input_data_widget = pn.Param(self.param.input_data),
        self.output_data_widget = pn.Param(self.param.output_data),
        self.input_data_ellipse_widget = pn.Param(self.param.input_data_ellipse)
        self.output_data_ellipse_widget = pn.Param(self.param.output_data_ellipse)
        self.widgets = pn.Param(self,parameters=['text_standard_selector','label_toggle','regression_selector','ellipse_mode_selector',
                                                 'update_output_button','export_data_button','export_TWplot_button',
                                                 'y_axis_TW',
                                                 'x_axis_TW_min','x_axis_TW_max',
                                                 'common_206204_input','common_207204_input','ThU_zrn_input','ThU_magma_input','power'])
    

    @pn.depends('file_path','regression_selector','ellipse_mode_selector',watch=True)
    def _uploadfile(self):
        if self.ellipse_mode_selector == 'Point Estimates':
            if self.file_path != 'Insert File Path':
                df = pd.read_excel(self.file_path,sheet_name='Sheet1')
                # df = pd.read_excel(self.file_path,sheet_name='5mJ_7Hz')
                self.input_data = df
                self.input_data['Pb207_Pb206'] = self.input_data['207Pb/206Pb']
                if '206/238 1st Order' in self.input_data.columns and self.regression_selector == '1st Order':
                    self.input_data['U238_Pb206'] = 1/self.input_data['206/238 1st Order']
                    self.input_data['206/238 Reg. err'] = self.input_data['SE 206/238 1st Order']
                    self.input_data['206/238U_age_init'] = np.log(self.input_data['206/238 1st Order'] + 1) / lambda_238
                elif '206/238 2nd Order' in self.input_data.columns and self.regression_selector == '2nd Order':
                    self.input_data['U238_Pb206'] = 1/self.input_data['206/238 2nd Order']
                    self.input_data['206/238 Reg. err'] = self.input_data['SE 206/238 2nd Order']
                    self.input_data['206/238U_age_init'] = np.log(self.input_data['206/238 2nd Order'] + 1) / lambda_238
                elif '206/238 Exp.' in self.input_data.columns and self.regression_selector == 'Exp.':
                    self.input_data['U238_Pb206'] = 1/self.input_data['206/238 Exp.']
                    self.input_data['206/238 Reg. err'] = self.input_data['SE 206/238 Exp']
                    self.input_data['206/238U_age_init'] = np.log(self.input_data['206/238 Exp.'] + 1) / lambda_238
                else:
                    pass
                self.output_data = pd.DataFrame([np.zeros(len(self.input_data.columns))],columns=list(self.input_data.columns))
            else:
                return print('No Input Data Available')
        else:
            pass
        
    @pn.depends('file_path_ellipse','ellipse_mode_selector',watch=True)
    def _uploadelllipsefile(self):
        if self.ellipse_mode_selector == 'Ellipses':
            if self.file_path_ellipse != 'Insert Ellipse File Path':
                df = pd.read_excel(self.file_path_ellipse,sheet_name='Sheet1')
                col_bdl_condn = df[(df['206Pb/238U']=='bdl') | (df['207Pb/206Pb']=='bdl') | (df['207Pb/235U'] == 'bdl')].index
                df.drop(col_bdl_condn,inplace=True) # drop rows that have 'bdl' in the 206/238 or 207/206 columns
                df = df.reset_index(drop=True)
                df['206Pb/238U'] = pd.to_numeric(df['206Pb/238U']) # having 'bdl's in columns makes them not numeric and thus they need to be changed to be used in calculations
                df['207Pb/206Pb'] = pd.to_numeric(df['207Pb/206Pb'])
                df['207Pb/235U'] = pd.to_numeric(df['207Pb/235U'])
                self.input_data_ellipse = df
                self.input_data_ellipse['Pb207_Pb206'] = self.input_data_ellipse['207Pb/206Pb']
                self.input_data_ellipse['U238_Pb206'] = 1/self.input_data_ellipse['206Pb/238U']
                self.input_data_ellipse['206/238U_age_init'] = np.log(self.input_data_ellipse['206Pb/238U'] + 1) / lambda_238
            self.output_data_ellipse = pd.DataFrame([np.zeros(4)],columns=['Ellipse Center','Ell. Width','Ell. Height','Ell. Rotation'])
        else:
            pass
    
    @pn.depends('input_data','input_data_ellipse','common_206204_input','common_207204_input','ellipse_mode_selector','text_sample_selector',watch=True)
    def _updateCommonPb(self):
        if self.text_sample_selector != 'Input Sample ID':
            if self.ellipse_mode_selector == 'Point Estimates':
                if self.common_206204_input != 0:
                    self.input_data['SK206_204'] = self.common_206204_input
                elif self.common_206204_input == 0:
                    self.input_data['SK206_204'] = 11.152 + 9.74*(np.exp(lambda_238*3.7e9)-np.exp(lambda_238*self.input_data['206/238U_age_init']))
                    # self.input_data['SK206_204'] = 11.152 + 10*(np.exp(lambda_238*3.7e9)-np.exp(lambda_238*self.input_data['206/238U_age_init']))
                if self.common_207204_input != 0:
                    self.input_data['SK207_204'] = self.common_207204_input
                elif self.common_207204_input == 0:
                    self.input_data['SK207_204'] = 12.998 + 9.74/137.82*(np.exp(lambda_235*3.7e9)-np.exp(lambda_235*self.input_data['206/238U_age_init']))
                    # self.input_data['SK207_204'] = 12.998 + 10/137.82*(np.exp(lambda_235*3.7e9)-np.exp(lambda_235*self.input_data['206/238U_age_init']))
                self.input_data['SK207_206'] = self.input_data['SK207_204'] / self.input_data['SK206_204']
            elif self.ellipse_mode_selector == 'Ellipses':
                if self.common_206204_input != 0:
                    self.input_data_ellipse['SK206_204'] = self.common_206204_input
                elif self.common_206204_input == 0:
                    self.input_data_ellipse['SK206_204'] = 11.152 + 9.74*(np.exp(lambda_238*3.7e9)-np.exp(lambda_238*self.input_data_ellipse['206/238U_age_init']))
                if self.common_207204_input != 0:
                    self.input_data_ellipse['SK207_204'] = self.common_207204_input
                elif self.common_207204_input == 0:
                    self.input_data_ellipse['SK207_204'] = 12.998 + 9.74/137.82*(np.exp(lambda_235*3.7e9)-np.exp(lambda_235*self.input_data_ellipse['206/238U_age_init']))
                self.input_data_ellipse['SK207_206'] = self.input_data_ellipse['SK207_204'] / self.input_data_ellipse['SK206_204']
        
    
    @pn.depends('input_data','input_data_ellipse','text_sample_selector','y_axis_TW','label_toggle',
                'x_axis_TW_min','x_axis_TW_max','ellipse_mode_selector','power')
    def call_TW(self):
        if self.ellipse_mode_selector == 'Point Estimates':
            if self.text_sample_selector != 'Input Sample ID':
                data_toplot = self.input_data[self.input_data['SampleLabel'].str.contains(self.text_sample_selector)]
                return calc_fncs.plot_TW(data_toplot,
                               self.x_axis_TW_min,self.x_axis_TW_max,
                               self.y_axis_TW[0],self.y_axis_TW[1],
                               self.label_toggle,self.ellipse_mode_selector,self.power)
        elif self.ellipse_mode_selector == 'Ellipses':
            if self.text_sample_selector != 'Input Sample ID':
                data_toplot = self.input_data_ellipse[self.input_data_ellipse['SampleLabel'].str.contains(self.text_sample_selector)]
                return calc_fncs.plot_TW(data_toplot,
                               self.x_axis_TW_min,self.x_axis_TW_max,
                               self.y_axis_TW[0],self.y_axis_TW[1],
                               self.label_toggle,self.ellipse_mode_selector,self.power)

    @pn.depends('input_data','text_sample_selector','text_standard_selector','label_toggle','ThU_zrn_input','ThU_magma_input','Pb_Th_std_crct_selector','regression_selector',
                'ellipse_mode_selector','power','UTh_std_norm')
    def call_boxplot(self):
        if self.ellipse_mode_selector == 'Point Estimates':
            if self.text_sample_selector != 'Input Sample ID':
                data_toplot = self.input_data[self.input_data['SampleLabel'].str.contains(self.text_sample_selector)]
                chosen_std = self.input_data[self.input_data['SampleLabel'].str.contains(self.text_standard_selector)]
                if self.text_sample_selector == 'Input Sample ID':
                    return 'Placeholder'
                else:
                    ages = calc_fncs.correct_sample_ages(data_toplot,chosen_std,self.text_standard_selector,self.ThU_zrn_input,self.ThU_magma_input,self.Pb_Th_std_crct_selector,self.regression_selector,
                                                         self.ellipse_mode_selector,self.power,self.UTh_std_norm)
                    return calc_fncs.plot_boxplot(ages['206Pb/238U_correctedage']/(1e6),ages['SampleLabel'],self.label_toggle,self.ellipse_mode_selector)
        else:
            pass
    
    
    
    
    def add_output_data(self,event=None):
        if self.ellipse_mode_selector == 'Point Estimates':
            data_to_update = self.input_data[self.input_data['SampleLabel'].str.contains(self.text_sample_selector)]
            chosen_std = self.input_data[self.input_data['SampleLabel'].str.contains(self.text_standard_selector)]
            ages = calc_fncs.correct_sample_ages(data_to_update,chosen_std,self.text_standard_selector,self.ThU_zrn_input,self.ThU_magma_input,self.Pb_Th_std_crct_selector,self.regression_selector,
                                                 self.ellipse_mode_selector,self.power,self.UTh_std_norm)
            if self.output_data is None:
                self.output_data = ages
            else:
                self.output_data = self.output_data.append(ages,ignore_index=True)
        elif self.ellipse_mode_selector == 'Ellipses':
            chosen_std = self.input_data_ellipse[self.input_data_ellipse['SampleLabel'].str.contains(self.text_standard_selector)]
            data_to_update = self.input_data_ellipse[self.input_data_ellipse['SampleLabel'].str.contains(self.text_sample_selector)]
            for i in data_to_update.SampleLabel.unique():
                data = calc_fncs.correct_sample_ages(data_to_update[data_to_update['SampleLabel']==i],chosen_std,self.text_standard_selector,self.ThU_zrn_input,self.ThU_magma_input,self.Pb_Th_std_crct_selector,self.regression_selector,
                                                     self.ellipse_mode_selector,self.power,self.UTh_std_norm)
                if self.output_data_ellipse is None:
                    self.output_data_ellipse = data
                else:
                    self.output_data_ellipse = self.output_data_ellipse.append(data,ignore_index=True)
        
    
    @pn.depends('output_data',watch=True)
    def _update_data_widget(self):
        if self.output_data is not None:
            self.output_data_widget = self.output_data
            self.output_data_widget.height = 40
            self.output_data_widget.heightpolicy = 'Fixed'
            return pn.widgets.Tabulator(self.output_data_widget,width=800)
    
    
    @pn.depends('output_data')
    def export_data(self,event=None):
        if self.ellipse_mode_selector == 'Point Estimates':
            self.output_data.to_excel('output_lasertramZ_ages.xlsx')
        elif self.ellipse_mode_selector == 'Ellipses':
            self.output_data_ellipse.to_excel('output_lasertramZ_ellipses.xlsx')
        
    
    def export_plot(self,event=None):
        if self.ellipse_mode_selector == 'Point Estimates':
            data_toplot = self.input_data[self.input_data['SampleLabel'].str.contains(self.text_sample_selector)]
            plot = calc_fncs.plot_TW(data_toplot,
                           self.x_axis_TW_min,self.x_axis_TW_max,
                           self.y_axis_TW[0],self.y_axis_TW[1],
                           self.label_toggle,self.ellipse_mode_selector,self.power)
            plot.savefig('LaserTRAMZ_TW.pdf',format='pdf',dpi=250)
        elif self.ellipse_mode_selector == 'Ellipses':
            data_toplot = self.input_data_ellipse[self.input_data_ellipse['SampleLabel'].str.contains(self.text_sample_selector)]
            plot = calc_fncs.plot_TW(data_toplot,
                           self.x_axis_TW_min,self.x_axis_TW_max,
                           self.y_axis_TW[0],self.y_axis_TW[1],
                           self.label_toggle,self.ellipse_mode_selector,self.power)
            plot.savefig('LaserTRAMZ_TW.pdf',format='pdf',dpi=250)

reduce_ages = finalize_ages(name='Reduce Ages')

#%%
# need to put regression selector widget in here
pn.extension('tabulator','mathjax')

width_ratios=[10,5]
grid_layout = pn.GridSpec(sizing_mode='scale_both')
# grid_layout = pn.GridSpec()

grid_layout[0,0] = pn.Column(reduce_ages.call_TW)
grid_layout[1,0] = pn.Row(reduce_ages._update_data_widget)

grid_layout[0,1] = pn.Row(pn.WidgetBox(pn.Param(reduce_ages.param,
                                         widgets={'label_toggle': pn.widgets.CheckBoxGroup,
                                                  'export_data_button': pn.widgets.Button(name='DDDT!',button_type='success'),
                                                  'regression_selector': pn.widgets.RadioButtonGroup,
                                                  'ellipse_mode_selector':pn.widgets.RadioButtonGroup,
                                                  'Pb_Th_std_crct_selector':pn.widgets.RadioButtonGroup,
                                                  'UTh_std_norm':pn.widgets.RadioBoxGroup}
                                               )
                                      )
                         )
grid_layout[1,1] = pn.Column(reduce_ages.call_boxplot)

grid_layout.show()
