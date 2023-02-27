#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 08:24:38 2023

@author: ctlewis
"""

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



class calc_fncs:
    """ Class that holds all of the functions for reducing the time resolved data"""
    def __init__(self,*args):
        for a in args:
            self.__setattr__(str(a), args[0])
    
    def get_ratios(data):
        """
        function that calculates relevant isotopic ratios for the U-Pb decay system 235/238 is returned strictly to make 
        plotting easy.
        
        Parameters
        ----------
        data : pandas dataframe
            pandas dataframe holding the observed time resolved LAICPMS measurements
    
        Returns
        -------
        data_ratio : pandas dataframe
            pandas dataframe with calculated isotopic ratios of interest. Note these are note the Hg corrected ratios for 204.
            These are the ratios that are used specifically for plotting and visualizing the measured ratios.
        """
        
        og_len = len(data.columns) # get the length of columns
        data = data.reset_index(drop=True) # reset index to always set the data into a format that allows data to be manipulated
        data_ratio = data.copy() # create copy of df so that nothing is overwritten
        # create empty lists to be filled with calculated ratios from the data
        pb206_pb204 = []
        pb206_u238 = []
        pb207_u235 = []
        pb207_pb206 = []
        u235_u238 = []
        u238_u235 = []
        
        data_ratio['Time_s2'] = data_ratio['Time_s']**2 # create a squared time term for second order regressions
        
        # create loops that append 1) calculated ratio for each observation if denominator > 0 or 
        # 2) zero in the case that denominator = 0 (avoids dividing by zero)
        
        
        # need to test for both cases being bdl separately. use or statement?
        for i in range(0,len(data_ratio)): # loop through range of the data
            if data_ratio['238U'][i] != 'bdl' and data_ratio['206Pb'][i] != 'bdl':
                if data_ratio['238U'][i] > 0 and data_ratio['206Pb'][i] > 0: # test for denominator > 0
                    pb206_u238.append((data_ratio['206Pb'][i]/data_ratio['238U'][i])) # append ratio if condition satisfied
                else: # append zero otherwise
                    pb206_u238.append(0)
            else:
                pb206_u238.append(0)
    
        for i in range(0,len(data_ratio)):
            if data_ratio['204Pb'][i] != 'bdl' and data_ratio['206Pb'][i] != 'bdl':
                if data_ratio['204Pb'][i] > 0 and data_ratio['206Pb'][i] > 0:
                    pb206_pb204.append((data_ratio['206Pb'][i]/data_ratio['204Pb'][i]))
                else:
                    pb206_pb204.append(0)
            else:
                pb206_pb204.append(0)
            
        for i in range(0,len(data_ratio)):
            if data_ratio['235U'][i] != 'bdl' and data_ratio['207Pb'][i] != 'bdl':
                if data_ratio['235U'][i] > 0 and data_ratio['207Pb'][i] > 0:
                    pb207_u235.append((data_ratio['207Pb'][i]/data_ratio['235U'][i]))
                else:
                    pb207_u235.append(0)
            else:
                pb207_u235.append(0)
            
        for i in range(0,len(data_ratio)):
            if data_ratio['206Pb'][i] != 'bdl' and data_ratio['207Pb'][i] != 'bdl':
                if data_ratio['206Pb'][i] > 0 and data_ratio['207Pb'][i] > 0:
                    pb207_pb206.append((data_ratio['207Pb'][i]/data_ratio['206Pb'][i]))
                else:
                    pb207_pb206.append(0)
            else:
                pb207_pb206.append(0)
            
        for i in range(0,len(data_ratio)):
            if data_ratio['238U'][i] != 'bdl' and data_ratio['235U'][i] != 'bdl':
                if data_ratio['238U'][i] > 0 and data_ratio['235U'][i] > 0:
                    u235_u238.append((data_ratio['235U'][i]/data_ratio['238U'][i]))
                else:
                    u235_u238.append(0)
            else:
                u235_u238.append(0)
                    
                
        for i in range(0,len(data_ratio)):
            if data_ratio['235U'][i] != 'bdl' and data_ratio['238U'][i] != 'bdl':
                if data_ratio['235U'][i] > 0 and data_ratio['238U'][i] > 0:
                    u238_u235.append((data_ratio['238U'][i]/data_ratio['235U'][i]))
                else:
                    u238_u235.append(0)
            else:
                u238_u235.append(0)
        
        # insert the lists into the copied dataframe
        data_ratio['Pb206_Pb204'] = pb206_pb204
        data_ratio['Pb206_U238'] = pb206_u238
        data_ratio['Pb207_U235'] = pb207_u235
        data_ratio['Pb207_Pb206'] = pb207_pb206
        data_ratio['U235_U238'] = u235_u238
        data_ratio['U238_U235'] = u238_u235
        
        data_ratio = data_ratio.iloc[:,(og_len-1):] # insert the calculated ratios onto the end of the copied dataframe
        
        return data_ratio
    
    
    def get_regressions(data,regression_buttons,ablation_start_true):
        """
        get_regressions is a function that handles all of the 206/238 regression data to deal with downhole fractionation.
        Various class methods can call this function. The return values depend on the calling method.
        
        Parameters
        ----------
        data : Pandas dataframe
            A dataframe that contains the calculated isotopic ratios (from the get_ratios function above)
        regression_buttons : list
            list of strings that indicate which regressions are requested by the user. These are selected by the panel 
            buttons in the code line param.regression_buttons = param.ListSelector
        ablation_start_true : Integer
            Integer value that that the user inputs to represent the start of the ablation time. Regression is projected back
            to this value
    
        Returns
        -------
        Return depends on which function called the regressions class method.
        
        If get_approved: 
        ratios_to_return - float: 
            list of values that are inserted into the exportable dataframe. These are means of the background subtracted data
        stats_to_return - float:
            list of statistics including various R-squared values, standard errors, and % standard errors as
            
        If ratio_plot or get_regression_stats:
        regressions_to_return - list of floats:
            list of float values that are put into the plot to visualize
        stats_to_report - list of floats &/or strings:
            list of floats and strings of various regression statistics. Updated in realtime
        
        If residuals_plot:
            predicted_to_return - list of floats:
                predicted values from the regression
            resid_to_return - list of floats:
                list of float values hosting the regression residuals
            names - list of strings:
                list of strings that indicates which order / type of regression belongs to which set of residuals.
        """
        if '1st Order' in regression_buttons:
            y1, X1 = dmatrices('Pb206_U238 ~ Time_s', data=data, return_type='dataframe') # get y and x regression data
            mod1 = sm.OLS(y1, X1) # fit a linear model on the data
            fit1 = mod1.fit() # get the list of fit parameters
            predicted1 = fit1.params[0] + fit1.params[1]*data.Time_s # get the predicted y values for the given x values
            rsquared1 = fit1.rsquared # get the R2 value of the regression
            predicted_b01 = fit1.params[0] + fit1.params[1]*ablation_start_true # get the predicted value at the ablation start that is input by the user
            sigma1 = np.sqrt(fit1.ssr/fit1.df_resid) # get the 1SD (Sum Squared residuals / residual degrees of freedom)^(1/2)
            # get standard error for a single point estiamted by a regression model
            # = sigma_regression * sqrt[1/n + (x_predicted - x_mean)^2 / ((df)*V(x))]
            # where sigma_regression is calculated above, x_predicted is the value at which you want to predict, x_mean is hte mean of all x
            # df is the degrees of freedom (n-#params), and V(x) is the variance of the predicted values
            SE_b01 = sigma1*np.sqrt(1/len(data)+(ablation_start_true-data['Time_s'].mean())**2/((len(data)-1)*statistics.variance(data['Time_s'])))
            SE_b01_percent = SE_b01/predicted_b01*100 # get the % 1SE
            resid1 = fit1.resid # get the residuals of the regression
        else:
            # fill the above with blank values so that there is not an error in the output if '1st Order' is not wanted by the user
            predicted1 = np.zeros_like(data['Time_s'])
            rsquared1 = 0
            predicted_b01 = 'NA'
            SE_b01 = 0
            SE_b01_percent = 0
            resid1 = np.zeros_like(data['Time_s'])
            
        # run an if loop as above. Virtually the same except we are using a 2nd order regression
        if '2nd Order' in regression_buttons:
            y2, X2 = dmatrices('Pb206_U238 ~ Time_s + Time_s2', data=data, return_type='dataframe')
            mod2 = sm.OLS(y2, X2)
            fit2 = mod2.fit()
            predicted2 = fit2.params[0] + fit2.params[1]*data.Time_s + fit2.params[2]*data.Time_s2
            rsquared2 = fit2.rsquared
            rsquared2_adj = fit2.rsquared_adj
            predicted_b02 = fit2.params[0] + fit2.params[1]*ablation_start_true + fit2.params[2]*ablation_start_true
            sigma2 = np.sqrt(fit2.ssr/fit2.df_resid)
            SE_b02 = sigma2*np.sqrt(1/len(data)+(ablation_start_true-data['Time_s'].mean())**2/((len(data)-1)*statistics.variance(data['Time_s'])))
            SE_b02_percent = SE_b02/predicted_b02*100
            resid2 = fit2.resid
        else:
            predicted2 = np.zeros_like(data['Time_s'])
            rsquared2 = 0
            rsquared2_adj = 0
            predicted_b02 = 'NA'
            SE_b02 = 0
            SE_b02_percent = 0
            resid2 = np.zeros_like(data['Time_s'])
    
        # get the method that called up regressions. f_back gets the function that called. Removing this gives current method
        callingmethod = sys._getframe().f_back.f_code.co_name
        # set a series of if statements that causes the appropriate return depending on the function that called up regresssions
        if callingmethod == 'get_approved':
            # get means of the calculated ratios
            pb206_204 = data['Pb206_Pb204'].mean()
            pb207_pb206 = data['Pb207_Pb206'].mean()
            u238_u235 = data['U238_U235'].mean()
            
            # get SEs and %SEs of the ratios
            pb206_204SE = data['Pb206_Pb204'].sem()
            pb206_204SE_percent = data['Pb206_Pb204'].sem()/data['Pb206_Pb204'].mean()*100
            pb207_206SE = data['Pb207_Pb206'].sem()
            pb207_206SE_percent = data['Pb207_Pb206'].sem()/data['Pb207_Pb206'].mean()*100
            U238_U235SE = data['U238_U235'].sem()
            U238_U235SE_percent = data['U238_U235'].sem()/data['U238_U235'].mean()*100
            
            # put the calculated values and statistics in lists to be returned
            ratios_to_return = [predicted_b01,predicted_b02,pb206_204,pb207_pb206,u238_u235]
            stats_to_return = [rsquared1,rsquared2,rsquared2_adj,SE_b01,SE_b01_percent,SE_b02,SE_b02_percent,pb206_204SE,pb206_204SE_percent,pb207_206SE,pb207_206SE_percent,U238_U235SE,U238_U235SE_percent]
            
            return ratios_to_return,stats_to_return
        
        elif callingmethod == 'ratio_plot':
            regressions_to_return = [predicted1,predicted2]
            stats_to_report = [rsquared1,rsquared2,rsquared2_adj,SE_b01_percent,SE_b02_percent]
            
            return regressions_to_return,stats_to_report
        
        elif callingmethod == 'get_regression_stats':
            # ""
            regressions_to_return = [predicted1,predicted2]
            stats_to_report = [rsquared1,rsquared2,rsquared2_adj,SE_b01_percent,SE_b02_percent]
            
            return regressions_to_return,stats_to_report
        
        elif callingmethod == 'residuals_plot':
            predicted_to_return = [predicted1,predicted2]
            resid_to_return = [resid1,resid2]
            
            return predicted_to_return,resid_to_return
        
        else:
            pass
        
        
        
    def get_approved(data,bckgrnd_start_input,bckgrnd_stop_input,
                     ablation_start_input,ablation_stop_input,ablation_start_true,
                    regression_buttons):
        """
        Function that finalizes regression interecepts / stats and calculated ratios, then puts them into the output dataframe
        
        Parameters
        ----------
        data : pandas dataframe
            pandas dataframe that has the data and calculated ratios
        bckgrnd_start_input : integer
            Integer value input by the user in the background_slider param. This is the lower value on the slider
        bckgrund_stop_input : integer
            Integer value input by the user in the background_slider param. This is the higher value on the slider
        ablation_start_input : integer
            Integer value input by the user in the ablation_slider param. This is the lower value on the slider
        ablation_stop_input : integer
            Integer value input by the user in the ablation_slider param. This is the higher value on the slider
        ablation_start_true : integer
            integer value input by the user in the ablation_start_true param. Regression is projected back to this value.
        regression_buttons : list of strings
            This is a list of strings that holds the type / order of regressions requested by the user through the 
            regression_bottoms param.
    
        Returns
        -------
        data_approved : pandas object
            A single observation with all observed and calculated data and regression statistics. Appended into the output data
    
        """
        data_toapprove = data.reset_index(drop=True) # reset the index of the data so that it can be altered/indexed appropriately
        og_col_length = len(data_toapprove.columns[2:]) # get original column length
        analyte_cols = data_toapprove.columns[2:og_col_length] # get the column names of the analytes that will go into ratios
        # need to send background subtracted data to regression, so do the following:
        background = data_toapprove[(data_toapprove.Time_s >= bckgrnd_start_input) & (data_toapprove.Time_s <= bckgrnd_stop_input)] # get the measured background across the selected interval
        lod = [] # create a list for detection limits to be filled
        # create a for loop to calculate the detection limits: LOD = (3SD * 2^(1/2)) / n^(1/2) (Longerich et al., 1996)
        for i in background.columns[2:-1]:
            limit = 3*background[i].std()/np.sqrt(len(background[i]))*np.sqrt(2)
            lod.append(limit)
        background = background.mean()[1:-1] # calculate the mean background
        # subtract the mean background from all data. Any values < 0 are assigned a value of zero to avoid errors when passing dataframe through functions.
        # these are assigned as 'bdl' later
        background_subtracted_data = data_toapprove.iloc[:,2:-1].sub(background,axis='columns').clip(lower=0)
        data_toapprove.iloc[:,2:-1] = background_subtracted_data # insert the background subtracted data into the copied dataframe
        

        data_approved_ratio = data_toapprove.copy() # copy the dataframe as to not overwrite the input data
        data_approved_ratio = calc_fncs.get_ratios(data_approved_ratio) # get calculated ratios of the background subtracted data. This is why we need 0 and not 'bdl'
        # index the background subtracted data for the input ablation periods.
        data_approved_ratio = data_approved_ratio[(data_approved_ratio.Time_s >= ablation_start_input) & (data_approved_ratio.Time_s <= ablation_stop_input)]
        # ellipse_data_approved = data_approved_ratio
        data_ratios,data_stats = calc_fncs.get_regressions(data_approved_ratio,regression_buttons,ablation_start_true) # get estimated intercepts, background subtracted ratios, regression statistics, and SE's of the ratios
        data_ratios = pd.DataFrame([data_ratios],columns=['206/238 1st Order','206/238 2nd Order','206Pb/204Pb','207Pb/206Pb','238U/235U']) # put the relevant calculations in df with appropriate column headers
        # put the relevant calculations in df with appropriate column headers
        data_stats = pd.DataFrame([data_stats],columns=['R2 1st Order','R2 2nd Order','R2 Adj 2nd Order','SE 206/238 1st Order','SE% 206/238 1st Order','SE 206/238 2nd Order','SE% 206/238 2nd Order',
                                                        'SE 206/204','SE% 206/204','SE 207/206','SE% 207/206','SE 238/235','SE% 238/235'])
        
        # get the ablation period indexed, background subtracted, indvidual isotope measurements
        data_toapprove = data_toapprove[(data_toapprove.Time_s >= ablation_start_input) & (data_toapprove.Time_s <= ablation_stop_input)]
        ellipse_data_toapprove = data_toapprove.copy()
        ellipse_data_toapprove = ellipse_data_toapprove.reset_index(drop=True)
        data_toapprove_SE = pd.DataFrame([data_toapprove.iloc[:,2:-1].sem()]).add_suffix('_1SE') # get SE's of the individual isotopes
        data_toapprove = pd.DataFrame([data_toapprove.iloc[:,2:-1].mean()]) # put means of isotopic values in a dataframe
        
        # 202Hg: 29.86%, 204Hg: 6.87% https://www.nndc.bnl.gov/nudat3/
        Hgratio = 6.87/100 / 29.86/100
        
        for i,k in zip(ellipse_data_toapprove.loc[:,list(analyte_cols)],lod):
            for j in range(0,len(ellipse_data_toapprove.iloc[:,0])):
                if ellipse_data_toapprove.loc[j,i]<k:
                    ellipse_data_toapprove.loc[j,i]='bdl'
                elif ellipse_data_toapprove.loc[j,i]>k:
                    ellipse_data_toapprove.loc[j,i]=ellipse_data_toapprove.loc[j,i]
        
        for i in range(0,len(ellipse_data_toapprove.iloc[:,0])):
            loc202hg = ellipse_data_toapprove.columns.get_loc('202Hg')
            loc204pb = ellipse_data_toapprove.columns.get_loc('204Pb')

            if ellipse_data_toapprove.iloc[i,loc202hg] != 'bdl':
                ellipse_data_toapprove['204Hg'] = np.zeros_like(ellipse_data_toapprove['202Hg'])
                loc204hg = ellipse_data_toapprove.columns.get_loc('204Hg')
                ellipse_data_toapprove.iloc[i,loc204hg] = ellipse_data_toapprove.iloc[i,loc202hg]*Hgratio

                if ellipse_data_toapprove.iloc[i,loc204pb] != 'bdl':
                    ellipse_data_toapprove.iloc[i,loc204pb] = ellipse_data_toapprove.iloc[i,loc204pb] - ellipse_data_toapprove.iloc[i,loc204hg]
                    if ellipse_data_toapprove.iloc[i,loc204pb] <= lod[loc204pb]:
                        ellipse_data_toapprove.iloc[i,loc204pb] = 'bdl'
                    elif ellipse_data_toapprove.iloc[i,loc204pb] > lod[loc204pb]:
                        ellipse_data_toapprove.iloc[i,loc204pb] = ellipse_data_toapprove.iloc[i,loc204pb]      
            else:
                pass
                        
        # take the means of the individual isotopes and run them through a for loop with the detection limits. Assign a value of 'bdl' if below detection limit. Otherwise, leave the 
        # value unchanged
        for i,k in zip(range(0,len(data_toapprove.iloc[0])),lod):
            if data_toapprove.iloc[0,i]<k:
                data_toapprove.iloc[0,i]='bdl'
            elif data_toapprove.iloc[0,i]>k:
                data_toapprove.iloc[0,i]=data_toapprove.iloc[0,i]
                
        # subtract the isobaric interference of 204Hg on 204Pb using the measured 202Hg, so long as 202 > 'bdl'
        if data_toapprove.loc[0,'202Hg'] != 'bdl':
            # 202Hg: 29.86%, 204Hg: 6.87% https://www.nndc.bnl.gov/nudat3/
            data_toapprove['204Hg'] = data_toapprove['202Hg']*Hgratio # calculate 204Hg from measured 202Hg based on isotopic abundance
            # subtract the 204Hg from the 204 signal, so long as the 204 signal > 'bdl'
            if data_toapprove.loc[0,'204Pb'] != 'bdl':
                data_toapprove['204Pb'] = data_toapprove['204Pb'] - data_toapprove['204Hg']
                # Recheck to make sure newly calculated if the newly calculated 204 signal is greater or less than 'bdl'. Assign 'bdl' or leave unchanged appropriately.
                loc204 = data_toapprove.columns.get_loc('204Pb')
                if data_toapprove.iloc[0,loc204] <= lod[loc204]:
                    data_toapprove['204Pb'] = 'bdl'
                elif data_toapprove.iloc[0,loc204] > lod[loc204]:
                    data_toapprove['204Pb'] = data_toapprove['204Pb']
        else:
            pass
        
        data_toapprove.insert(0,'SampleLabel',data.iloc[0,0]) # reinsert the sample label into the calculations
        # stitch the individual isotopic measurements, their errors, ratios, ratio errors, and regression results / regression statistics into a df together.
        # these are then appeneded into the output df
        data_approved = data_toapprove.join([data_toapprove_SE,data_ratios,data_stats])
        
        ratio_cols = ['206Pb/204Pb','207Pb/206Pb','238U/235U']
        for i in analyte_cols:
            for k in ratio_cols:
                if k.__contains__(str(i)) and data_approved[i].item() == 'bdl':
                    data_approved[k] = 'bdl'
                else:
                    pass
                    
        ellipse_data_toapprove_ratio = calc_fncs.get_ratios(ellipse_data_toapprove).drop(['Time_s2','204Hg'],axis=1)
        ellipse_data_toapprove_ratio.rename({'Pb206_Pb204':'206Pb/204Pb' , 'Pb206_U238':'206Pb/238U', 'Pb207_U235':'207Pb/235U', 'Pb207_Pb206':'207Pb/206Pb', 'U235_U238':'235U/238U', 'U238_U235':'238U/235U'}, axis=1, inplace=True)
        ellipse_data_toapprove = pd.concat([ellipse_data_toapprove,ellipse_data_toapprove_ratio],axis=1)
        
        ratio_colsII = ['206Pb/204Pb','207Pb/206Pb','207Pb/235U','235U/238U','238U/235U']
        for i in analyte_cols:
            for k in ratio_colsII:
                if k.__contains__(str(i)):
                    for j in range(0,len(ellipse_data_toapprove)):
                        if ellipse_data_toapprove.loc[j,i] == 'bdl':
                            ellipse_data_toapprove.loc[j,k] = 'bdl'
                    
        
        return data_approved,ellipse_data_toapprove
    
    def get_ellipse(data,power):
        
        # data = data.reset_index(drop=True)
        # data = calc_fncs.get_ratios(data)
        x1 = data['Pb207_U235']
        y1 = data['Pb206_U238']
        x2 = 1/data['Pb206_U238']
        y2 = data['Pb207_Pb206']
        x3 = data['Pb206_Pb204']
        y3 = data['Pb207_Pb206']
        
        cov1 = np.cov(x1,y1)
        cov2 = np.cov(x2,y2)
        cov3 = np.cov(x3,y3)
        eigval1,eigvec1 = np.linalg.eig(cov1)
        eigval2,eigvec2 = np.linalg.eig(cov2)
        eigval3,eigvec3 = np.linalg.eig(cov3)
        order1 = eigval1.argsort()[::-1]
        order2 = eigval2.argsort()[::-1]
        order3 = eigval3.argsort()[::-1]
        eigvals_order1 = eigval1[order1]
        eigvals_order2 = eigval2[order2]
        eigvals_order3 = eigval3[order3]
        eigvecs_order1 = eigvec1[:,order1]
        eigvecs_order2 = eigvec2[:,order2]
        eigvecs_order3 = eigvec3[:,order3]
        
        c1 = (np.mean(x1),np.mean(y1))
        c2 = (np.mean(x2),np.mean(y2))
        c3 = (np.mean(x3),np.mean(y3))
        wid1 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power),df=2)*eigvals_order1[0])
        hgt1 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power),df=2)*eigvals_order1[1])
        wid2 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power),df=2)*eigvals_order2[0])
        hgt2 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power),df=2)*eigvals_order2[1])
        wid3 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power),df=2)*eigvals_order3[0])
        hgt3 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power),df=2)*eigvals_order3[1])
        theta1 = np.degrees(np.arctan2(*eigvecs_order1[:,0][::-1]))
        theta2 = np.degrees(np.arctan2(*eigvecs_order2[:,0][::-1]))
        theta3 = np.degrees(np.arctan2(*eigvecs_order3[:,0][::-1]))
        
        ell1_params = [c1,wid1,hgt1,theta1]
        ell2_params = [c2,wid2,hgt2,theta2]
        ell3_params = [c3,wid3,hgt3,theta3]
        
        return ell1_params,ell2_params,ell3_params
        
        
class plots(calc_fncs):
    """ Class that holds all of the functions for reducing the time resolved data"""
    def __init__(self,*args):
        super().__init__(*args)
        for a in args:
            self.__setattr__(str(a), args[0])
            
    
    def ratio_plot(data,ratio_buttons,regression_buttons,
                   start_bckgrnd,stop_bckgrnd,
                   start_ablation,stop_ablation,ablation_start_true,
                   ratio_plot_ylim_slider_min,ratio_plot_ylim_slider_max):
        """
        Function for displaying the plot of relavent isotopic ratios.
        
        Parameters
        ----------
        data : pandas dataframe
            pandas dataframe of observed data
        ratio_buttons : list of strings
            list of strings that hosts the ratios the user requests to plot
        regression_buttons : list of strings
            list of strings that hoststhe regressions the user requests
        start_bckgrnd : integer
            Integer value input by the user in the background_slider param. This is the lower value on the slider 
        stop_bckgrnd : integer
            Integer value input by the user in the background_slider param. This is the higher value on the slider
        start_ablation : integer
            Integer value input by the user in the ablation_slider param. This is the lower value on the slider
        stop_ablation : integer
            Integer value input by the user in the ablation_slider param. This is the higher value on the slider
        ablation_start_true : integer
            integer value input by the user in the ablation_start_true param. Regression is projected back to this value
        ratio_plot_ylim_slider_min : float
            float values input by the user in the ratio_plot_ylim_min param. Sets min ylim on ratio plot
        ratio_plot_ylim_slider_max : float
            float values input by the user in the ratio_plot_ylim_max param. Sets max ylim on ratio plot
    
        Returns
        -------
        fig : matplotlib figure
            plot that shows the relevant isotopic ratios for the entire measurement (background+ablation+washout)
    
        """
        
        plt.style.use('seaborn-colorblind') # use colorscheme that is viewable by colorblind folks
        fig = Figure(figsize=(12,8)) # create a figure that is 12Wx8H
        ax = fig.add_subplot() # add a subplot. This initates the first set of axes onto the renderable figure
        data = calc_fncs.get_ratios(data) # get calculated ratios from the data
        data_to_regress = data[(data.Time_s >= start_ablation) & (data.Time_s <= stop_ablation)] # get the selected ablation period
        var_cols = ratio_buttons # get the ratios selected by the user. These are used to get the columns from the calculated ratios
        regressions_to_plot = regression_buttons# get the regression type requested by the user
        # get regression parameters and stats for the requested regressions across the specified ablation period
        regressions,stats = calc_fncs.get_regressions(data_to_regress,regression_buttons,ablation_start_true)
        
        # plot a line for each selected ratio.
        for i in var_cols:
            ax.plot(data.Time_s,data[i],'-',lw=1,label='{}'.format(i))
        
        # plot a line for each selected regression
        if regressions_to_plot is not None:
            for i in regressions:
                ax.plot(data_to_regress.Time_s,i)
        
        # plot vertical lines for selected background period, ablation period, and ablation start values
        ax.plot([start_bckgrnd,start_bckgrnd],[0,max(data['Pb206_U238'])],'--k',lw=0.75)
        ax.plot([stop_bckgrnd,stop_bckgrnd],[0,max(data['Pb206_U238'])],'--k',lw=0.75)
        ax.plot([start_ablation,start_ablation],[0,max(data['Pb206_U238'])],'-k',lw=0.75)
        ax.plot([stop_ablation,stop_ablation],[0,max(data['Pb206_U238'])],'-k',lw=0.75)
        ax.plot([ablation_start_true,ablation_start_true],[0,1e8],':',c='k',lw=0.75)
        ax.legend(loc='upper right',fontsize=18) # get a legened of ratios. These are assinged their uniqeu value with label='{}'.format(i)
        ax.set_ylim(ratio_plot_ylim_slider_min,ratio_plot_ylim_slider_max) # set the ylim to the requested value
        ax.set_xlabel('Time (s)',fontsize=18) # put xlabel on plot
        ax.set_ylabel('Selected Ratios',fontsize=18) # put ylabel on plot
        ax.tick_params(axis='both',labelsize=14) # change some plot aesthetics
        fig.tight_layout() # set the figure to be 'tight', as opposed to lots of blank white space
        return fig
    
    
    def ablation_plot(data_ablation_,bckgrnd_start_input,bckgrund_stop_input,ablation_start_input,ablation_stop_input,
                      ablation_start_true,ablation_plot_ylim_slider):  
        """
        Function to plot the measured isotopic data across the entire analysis (background+ablation+washout)
    
        Parameters
        ----------
        data_ablation_ : pandas dataframe
            dataframe of the measured values.
        bckgrnd_start_input : integer
            lower integer value of the background_slider param.
        bckgrund_stop_input : integer
            upper integer value of the background_slider param.
        ablation_start_input : integer
            lower value of the abaltion_slider param.
        ablation_stop_input : integer
            upper value of the ablation_slider param.
        ablation_start_true : integer
            integer value that indicates where the regression will be projected back towards (i.e., start of ablation).
        ablation_plot_ylim_slider : integer
            integer value that sets the max y-lim of the plot.
    
        Returns
        -------
        fig : matplotlib fig
            Figure showing all of the ablation data.
    
        """
        plt.style.use('seaborn-colorblind') # use colorscheme that is viewable by colorblind folks
        fig = Figure(figsize=(12,8)) # create a figure that is 12Wx8H
        ax = fig.add_subplot() # add a subplot. This initates the first set of axes onto the renderable figure
        var_cols = data_ablation_.columns[3:-1] # get all of the measured isotopes, starting at col 3. This is one of hte moany reasons the input format is specific
        # plot a line for each of the measured isotopes
        for i in var_cols:
            ax.plot(data_ablation_.Time_s,data_ablation_[i],'-',lw=1,label='{}'.format(i))
        # plot a vertical line with various parameters for background period, ablation period, and ablation start value
        ax.plot([bckgrnd_start_input,bckgrnd_start_input],[0,1e8],'--k',lw=0.75)
        ax.plot([bckgrund_stop_input,bckgrund_stop_input],[0,1e8],'--k',lw=0.75)
        ax.plot([ablation_start_input,ablation_start_input],[0,1e8],'-k',lw=0.75)
        ax.plot([ablation_stop_input,ablation_stop_input],[0,1e8],'-k',lw=0.75)
        ax.plot([ablation_start_true,ablation_start_true],[0,1e8],':',c='k',lw=0.75)
        ax.legend(loc='upper right',fontsize=18) # plot a legened
        ax.set_ylim(0,ablation_plot_ylim_slider) # set the ylim from 0 to max, where max is the requested value
        ax.set_xlabel('Time (s)',fontsize=18) # xlabel
        ax.set_ylabel('Counts',fontsize=18) # ylabel
        ax.tick_params(axis='both',labelsize=14) # plot aesthetics
        fig.tight_layout() # make the fig layout tight, as opposed to a bunch of white space
        
        return fig
        
        
    def residuals_plot(data,regression_buttons,start_ablation,stop_ablation,ablation_start_true):
        plt.style.use('seaborn-colorblind') # use colorscheme that is viewable by colorblind folks
        fig = Figure(figsize=(6,3)) # create a figure that is 6Wx3H
        ax = fig.add_subplot() # add a subplot. This initates the first set of axes onto the renderable figure
        data = calc_fncs.get_ratios(data) # calculate relevant isotopic ratios from the data
        data_to_regress = data[(data.Time_s >= start_ablation) & (data.Time_s <= stop_ablation)] # get data across the requested ablation interval
        fitted,residuals = calc_fncs.get_regressions(data_to_regress,regression_buttons,ablation_start_true) # get fitted regression values, regression names, and residuals
        
        # plot the residuals for the regressions and the fitted value
        for i,k,j in zip(fitted,residuals,regression_buttons):
            if i is not 0:
                ax.plot(i,k,'o',mec='k',label='{}'.format(j))
                ax.plot([min(i),max(i)],[0,0],'-k',lw=0.5)
        ax.legend(loc='upper right',fontsize=8) # get the legend
        ax.set_xlabel('Fitted Value',fontsize=10) # xlabel
        ax.set_ylabel('Residuals',fontsize=10) # ylabel
        ax.tick_params(axis='both',labelsize=8) # plot aesthetics
        fig.tight_layout() # make the figure tight, as opposed to lots of white background
        
        return fig
    
    def ellipse_plot(data,power,start_ablation,stop_ablation):
        data = calc_fncs.get_ratios(data)
        data = data[(data.Time_s >= start_ablation) & (data.Time_s <= stop_ablation)]
        ell1p,ell2p,ell3p = calc_fncs.get_ellipse(data, power)
        ell1 = Ellipse(ell1p[0],ell1p[1],ell1p[2],ell1p[3],color='darkslategray',alpha=0.5)
        ell2 = Ellipse(ell2p[0],ell2p[1],ell2p[2],ell2p[3],color='darkslategray',alpha=0.5)
        ell3 = Ellipse(ell3p[0],ell3p[1],ell3p[2],ell3p[3],color='darkslategray',alpha=0.5)
        plt.style.use('seaborn-colorblind') # use colorscheme that is viewable by colorblind folks
        fig1 = Figure(figsize=(8,8)) # create a figure that is 8Wx8H
        fig2 = Figure(figsize=(8,8))
        fig3 = Figure(figsize=(8,8))
        ax1 = fig1.add_subplot()
        ax2 = fig2.add_subplot()
        ax3 = fig3.add_subplot()
        
        ax1.add_artist(ell1)
        ax1.plot(data['Pb207_U235'],data['Pb206_U238'],'.k')
        ax2.add_artist(ell2)
        ax2.plot(1/data['Pb206_U238'],data['Pb207_Pb206'],'.k')
        ax3.add_artist(ell3)
        ax3.plot(data['Pb206_Pb204'],data['Pb207_Pb206'],'.k')
        
        ax1.set_xlabel('207/235',fontsize=10)
        ax1.set_ylabel('206/238',fontsize=10)
        ax2.set_xlabel('238/206',fontsize=10)
        ax2.set_ylabel('207/206',fontsize=10)
        ax3.set_xlabel('206/204',fontsize=10)
        ax3.set_ylabel('207/206',fontsize=10)
        ax1.tick_params(axis='both',labelsize=8)
        ax2.tick_params(axis='both',labelsize=8)
        ax3.tick_params(axis='both',labelsize=8)
        ax1.set_xlim(ell1p[0][0]-ell1p[1]/1.5,ell1p[0][0]+ell1p[1]/1.5)
        ax1.set_ylim(ell1p[0][1]-ell1p[2]/1.5,ell1p[0][1]+ell1p[2]/1.5)
        ax2.set_xlim(ell2p[0][0]-ell2p[1]/1.5,ell2p[0][0]+ell2p[1]/1.5)
        ax2.set_ylim(ell2p[0][1]-ell2p[2]/1.5,ell2p[0][1]+ell2p[2]/1.5)
        ax3.set_xlim(ell3p[0][0]-ell3p[1]/1.5,ell3p[0][0]+ell3p[1]/1.5)
        ax3.set_ylim(ell3p[0][1]-ell3p[2]/1.5,ell3p[0][1]+ell3p[2]/1.5)
        
        fig1.tight_layout()
        fig2.tight_layout()
        fig3.tight_layout()
        
        return fig1,fig2,fig3

        
class make_plots(param.Parameterized):
    """ class that parameterizes inputs and sends them to the above functions to be rendered in a GUI"""
    input_data = param.DataFrame(precedence=-1)
    file_path = param.String(default='Insert File Path')
    output_data = param.DataFrame(precedence=-1)
    output_data_ellipse = param.DataFrame(precedence=-1)
    sample_subset = param.Selector(objects=[])
    
    update_output_button = param.Action(lambda x: x.add_output_data(),label='Approve Interval')
    export_data_button = param.Action(lambda x: x.export_data(),label='DDDT!')
    export_plots_button = param.Action(lambda x: x.export_plots(),label ='Export All Plots')
    
    ablation_start_true = param.Integer(20)
    ablation_slider = param.Range(default=(23,33),bounds=(0,100))
    background_slider = param.Range(default=(3,19),bounds=(0,100))
    ablation_plot_ylim_slider = param.Integer(default=1000000,bounds=(0,10000000))
    ratio_plot_ylim_min = param.Number(default=0.01)
    ratio_plot_ylim_max = param.Number(default=0.5)
    ratio_buttons = param.ListSelector(default=['Pb206_U238'], objects=['Pb206_U238','Pb206_Pb204','Pb207_U235','Pb207_Pb206','U235_U238'])
    regression_buttons = param.ListSelector(default=['1st Order'], objects=['1st Order','2nd Order'])
    power = param.Number(default=0.05)
        
        
    def __init__(self,**params):
        super().__init__(**params)
        self.file_input_widget = pn.Param(self.param.input_data)
        self.output_data_widget = pn.Param(self.param.output_data)
        self.output_data_ellipse_widget = pn.Param(self.param.output_data_ellipse)
        self.widgets = pn.Param(self,parameters=['file_path','sample_subset','update_output_button','export_data_button',
                                                 'ablation_start_true','ablation_slider','background_slider','ablation_plot_ylim_slider',
                                                 'ratio_plot_ylim_slider','ratio_plot_ylim_max','ratio_buttons','regression_buttons','power'])
    
    @pn.depends('file_path',watch=True)
    def _uploadfile(self):
        if self.file_path != 'Insert File Path':
            df = pd.read_excel(self.file_path,sheet_name='Sheet1')
            self.input_data = df
            self.input_data['Time_s'] = self.input_data['Time']/1000
            options = list(self.input_data.SampleLabel.unique())
            self.param.sample_subset.objects = options
            self.output_data = pd.DataFrame([np.zeros(len(self.input_data.columns))],columns=list(self.input_data.columns))
            self.output_data_ellipse = pd.DataFrame([np.zeros(len(self.input_data.columns))],columns=list(self.input_data.columns))
        else:
            return print('No Input Data Available')
    
        
    @pn.depends('input_data','background_slider','ablation_slider','ablation_start_true','sample_subset','ablation_plot_ylim_slider')
    def call_ablation_plot(self):
        if self.sample_subset is not None:
            data_toplot = self.input_data[self.input_data['SampleLabel'] == self.sample_subset]
            return plots.ablation_plot(data_toplot,self.background_slider[0],self.background_slider[1],
                                  self.ablation_slider[0],self.ablation_slider[1],self.ablation_start_true,self.ablation_plot_ylim_slider)
        else:
            pass
        
    @pn.depends('input_data','ratio_buttons','regression_buttons',
                'background_slider','ablation_slider','ablation_start_true',
                'sample_subset',
                'ratio_plot_ylim_min','ratio_plot_ylim_max')
    def call_ratio_plot(self):
        if self.sample_subset is not None:
            data_toplot = self.input_data[self.input_data['SampleLabel'] == self.sample_subset]
            return plots.ratio_plot(data_toplot,self.ratio_buttons,self.regression_buttons,self.background_slider[0],self.background_slider[1],
                              self.ablation_slider[0],self.ablation_slider[1],self.ablation_start_true,
                              self.ratio_plot_ylim_min,self.ratio_plot_ylim_max)
    
    @pn.depends('input_data','regression_buttons','ablation_slider','sample_subset','ablation_start_true')
    def call_residuals_plot(self):
        if self.sample_subset is not None:
            data_toplot = self.input_data[self.input_data['SampleLabel'] == self.sample_subset]
            return plots.residuals_plot(data_toplot,self.regression_buttons,self.ablation_slider[0],self.ablation_slider[1],self.ablation_start_true)
        
    @pn.depends('input_data','power','sample_subset','ablation_slider')
    def call_ellipse_plot(self):
        if self.sample_subset is not None:
            data_toplot = self.input_data[self.input_data['SampleLabel'] == self.sample_subset]
            Weth_ell,TW_ell,Pb_ell = plots.ellipse_plot(data_toplot, self.power,self.ablation_slider[0],self.ablation_slider[1])
            tabs = pn.Tabs(('TW',TW_ell),('Weth.',Weth_ell),('PbPb',Pb_ell),dynamic=True)
            return tabs
            
        
    def export_plots(self,event=None):
        data_toplot = self.input_data[self.input_data['SampleLabel'] == self.sample_subset]
        ablationplot = plots.ablation_plot(data_toplot,self.background_slider[0],self.background_slider[1],self.ablation_slider[0],self.ablation_slider[1],self.ablation_start_true,self.ablation_plot_ylim_slider)
        ratioplot = plots.ratio_plot(data_toplot,self.ratio_buttons,self.regression_buttons,self.background_slider[0],self.background_slider[1],self.ablation_slider[0],self.ablation_slider[1],self.ablation_start_true,
                          self.ratio_plot_ylim_min,self.ratio_plot_ylim_max)
        residualsplot = plots.residuals_plot(data_toplot,self.regression_buttons,self.ablation_slider[0],self.ablation_slider[1],self.ablation_start_true)
        ablationplot.savefig('Ablation plot.pdf',format='pdf',dpi=250)
        ratioplot.savefig('Ratio plot.pdf',format='pdf',dpi=250)
        residualsplot.savefig('Residuals plot.pdf',format='pdf',dpi=250)
    
    @pn.depends('input_data','regression_buttons','ablation_slider','sample_subset','ablation_start_true')
    def get_regression_stats(self):
        if self.sample_subset is not None:
            data_toanalyze = self.input_data[self.input_data['SampleLabel'] == self.sample_subset]
            data_toanalyze = plots.get_ratios(data_toanalyze)
            data_toanalyze = data_toanalyze[(data_toanalyze.Time_s >= self.ablation_slider[0]) & (data_toanalyze.Time_s <= self.ablation_slider[1])]
            regressions,stats = calc_fncs.get_regressions(data_toanalyze,self.regression_buttons,self.ablation_start_true)
            r2_1 = pn.pane.Markdown('$$R^{2}$$ = '+str(round(stats[0],3)))
            r2_2 = pn.pane.Markdown('$$R^{2} 2^{nd} order$$ = '+str(round(stats[1],3)))
            r2adj_2 = pn.pane.Markdown('$$R^{2}_{adj}$$ = '+str(round(stats[2],3)))
            SE_1stper = pn.pane.Markdown('$$1^{st} Order SE \%$$ = '+str(round(stats[3],3)))
            SE_2ndper = pn.pane.Markdown('$$2^{nd} Order SE \%$$ = '+str(round(stats[4],3)))
            stats_markdown = pn.Column(r2_1,r2_2,r2adj_2,SE_1stper,SE_2ndper)
            # stats_markdown = pn.Column(trythis)
            return stats_markdown
        

    @pn.depends('output_data',watch=True)
    def _update_output_widget(self):
        if self.output_data is not None:
            self.output_data_widget = self.output_data
            self.output_data_widget.height = 400
            self.output_data_widget.heightpolicy = 'Fixed'
            return pn.widgets.Tabulator(self.output_data_widget,width=600) # use 600 for large screen, 100-150 for small screen 
    
    def add_output_data(self, event=None):
        data_tosend = self.input_data[self.input_data['SampleLabel'] == self.sample_subset]
        data_approved,ellipse_datapproved = calc_fncs.get_approved(data_tosend,self.background_slider[0],self.background_slider[1],self.ablation_slider[0],self.ablation_slider[1],
                                     self.ablation_start_true,self.regression_buttons)
        if self.output_data is None:
            self.output_data = data_approved
        else:
            self.output_data = self.output_data.append(data_approved,ignore_index=True)
            
        if self.output_data_ellipse is None:
            self.output_data_ellipse = ellipse_datapproved
        else:
            self.output_data_ellipse = self.output_data_ellipse.append(ellipse_datapproved,ignore_index=True)
            
    @pn.depends('output_data','output_data_ellipse')
    def export_data(self,event=None):
        self.output_data.to_excel('output_lasertramZ.xlsx')
        self.output_data_ellipse.to_excel('output_CEllipse_lasertramZ.xlsx')
        

callapp = make_plots(name='Reduce Ablation Data')

pn.extension('tabulator','mathjax')

grid_layout = pn.GridSpec(sizing_mode='scale_both')
grid_layout[:,2] = pn.Column(pn.WidgetBox(pn.Param(callapp.param,
                                                   widgets = {'ratio_buttons': pn.widgets.CheckBoxGroup,
                                                              'regression_buttons': pn.widgets.CheckBoxGroup,
                                                              'export_data_button': pn.widgets.Button(name='DDDT!',button_type='success')}))
                             )
grid_layout[0,0] = pn.Column(callapp.call_ablation_plot)
grid_layout[0,1] = pn.Column(callapp.call_ratio_plot)
grid_layout[1,0] = pn.Column(callapp.call_residuals_plot)
# grid_layout[1,1] = pn.Column(callapp.get_regression_stats)
grid_layout[2,0] = pn.Column(callapp._update_output_widget)
grid_layout[1,1] = pn.Column(callapp.call_ellipse_plot)


grid_layout.show()







    