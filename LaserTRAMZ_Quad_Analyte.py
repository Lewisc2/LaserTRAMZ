#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 07:20:05 2023

@author: Chuck Lewis, Oregon State University
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 08:24:38 2023

@author: ctlewis
"""

import pandas as pd
import numpy as np
import bokeh
from bokeh.plotting import figure
from bokeh.layouts import column, row, gridplot
from bokeh.models import Ellipse, Plot
import holoviews as hv
import panel as pn
import statistics
import param
import sys
import statsmodels.api as sm
from patsy import dmatrices
import scipy
from scipy import stats
from scipy.optimize import curve_fit
from itertools import cycle

color_palette = bokeh.palettes.Accent8
color_palette_regressions = bokeh.palettes.Dark2_3
hv.extension('bokeh')

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
            if data_ratio['235U'][i] != 'bdl' and data_ratio['238U'][i] != 'bdl':
                if data_ratio['235U'][i] > 0 and data_ratio['238U'][i] > 0:
                    u238_u235.append((data_ratio['238U'][i]/data_ratio['235U'][i]))
                else:
                    u238_u235.append(0)
            else:
                u238_u235.append(0)
                
        
        
        # insert the lists into the copied dataframe
        data_ratio['206Pb/204Pb'] = pb206_pb204
        data_ratio['206Pb/238U'] = pb206_u238
        data_ratio['207Pb/235U'] = pb207_u235
        data_ratio['207Pb/206Pb'] = pb207_pb206
        data_ratio['238U/235U'] = u238_u235
        
        data_ratio = data_ratio.iloc[:,(og_len-1):] # insert the calculated ratios onto the end of the copied dataframe
        
        return data_ratio
    
    
    def get_counts(data,counts_mode,arrayofdwelltimes,ablation_start_true,bckgrnd_start_input,bckgrnd_stop_input,ablation_start_input,ablation_stop_input):
        """
        Function that gets the ratio / counts data for the relevant ratios and analytes.

        Parameters
        ----------
        data : dataframe
            pandas dataframe holding the data from the current analysis.
        counts_mode : string
            string denoting which method the using wants to reduce data with.
        arrayofdwelltimes : array
            Array of dwell times in s.
        ablation_start_true : float
            float value of the projected ablation start time.
        bckgrnd_start_input : integer
            integer from the slider defining when to start taking the gas blank.
        bckgrnd_stop_input : integer
            integer from the slider defining when to stop taking the gas blank.
        ablation_start_input : integer
            integer from the slider defining when to start the regression / counts.
        ablation_stop_input : integer
            integer form the slider defining when to stop the regression / counts.

        Returns
        -------
        ratios_to_return : array
            array of values including relevant ratios and analyte intensities when applicable.
        stats_to_return : array
            array of values including relevant errors and regression statistics.

        """
        
        if counts_mode == 'Total Counts':
            data_totalcounts = np.sum(data.loc[:,'202Hg':'238U'] * arrayofdwelltimes)
            data_ratios_ = calc_fncs.get_ratios(data)
            pb206_204 = data_totalcounts['206Pb']/data_totalcounts['204Pb']
            pb207_206 = data_totalcounts['207Pb']/data_totalcounts['206Pb']
            u238_235 = data_totalcounts['238U']/data_totalcounts['235U']
            pb207_u235 = data_totalcounts['207Pb']/data_totalcounts['235U']
            pb206_u238 = data_totalcounts['206Pb']/data_totalcounts['238U']
            
            pb206_204SE = data_ratios_['206Pb/204Pb'].sem()
            pb206_204SE_percent = data_ratios_['206Pb/204Pb'].sem()/pb206_204*100
            pb207_206SE = data_ratios_['207Pb/206Pb'].sem()
            pb207_206SE_percent = data_ratios_['207Pb/206Pb'].sem()/pb207_206*100
            u238_235SE = data_ratios_['238U/235U'].sem()
            u238_235SE_percent = data_ratios_['238U/235U'].sem()/u238_235*100
            
            pb207_u235SE = np.sqrt(pb207_u235)/np.sqrt(len(data_ratios_))
            pb207_u235SE_percent = np.sqrt(pb207_u235)/np.sqrt(len(data_ratios_))/pb207_u235*100
            
            pb206_u238SE = np.sqrt(pb206_u238)/np.sqrt(len(data_ratios_))
            pb206_u238SE_percent = np.sqrt(pb206_u238)/np.sqrt(len(data_ratios_))/pb206_u238*100
            
            ratios_to_return = [pb206_u238,pb207_u235,pb206_204,pb207_206,u238_235]
            stats_to_return = [pb206_u238SE,pb206_u238SE_percent,pb207_u235SE,pb207_u235SE_percent,pb206_204SE,pb206_204SE_percent,pb207_206SE,pb207_206SE_percent,u238_235SE,u238_235SE_percent]
            
            return ratios_to_return,stats_to_return
            
        elif counts_mode == 'Poisson':    
            data_ratios_ = calc_fncs.get_ratios(data)
            x = data['Time_s']
            pb207poisition = data.columns.get_loc('207Pb')
            pb206poisition = data.columns.get_loc('206Pb')
            pb204poisition = data.columns.get_loc('204Pb')
            pb207 = data['207Pb']*arrayofdwelltimes[pb207poisition-2]
            pb207 = np.ceil(pb207)
            pb206 = data['206Pb']*arrayofdwelltimes[pb206poisition-2]
            pb206 = np.ceil(pb206)
            pb204 = data['204Pb']*arrayofdwelltimes[pb204poisition-2]
            pb204 = np.ceil(pb204)
            
            y207,x207 = dmatrices('pb207 ~ x',return_type='dataframe')
            y206,x206 = dmatrices('pb206 ~ x',return_type='dataframe')
            y204,x204 = dmatrices('pb204 ~ x',return_type='dataframe')
            
            mod207zero = sm.ZeroInflatedPoisson(y207,x207)
            mod206pois = sm.Poisson(y206,x206)
            mod204zero = sm.ZeroInflatedPoisson(y204,x204)
            fit207zero = mod207zero.fit(method='basinhopping')
            fit206pois = mod206pois.fit(method='basinhopping')
            fit204zero = mod204zero.fit(method='basinhopping')
            pb207zero = np.mean(fit207zero.predict())
            pb206pois = np.mean(fit206pois.predict())
            pb204zero = np.mean(fit204zero.predict())
            variance207zero = fit207zero.predict(which='var')
            variance206pois = fit206pois.predict(which='var')
            variance204zero = fit204zero.predict(which='var')
            
            pb207SEzero = np.sqrt(np.mean(variance207zero))/np.sqrt(len(variance207zero)-1)
            pb207SEpercentzero = pb207SEzero / np.mean(fit207zero.predict()) * 100
            pb206SEpois = np.sqrt(np.mean(variance206pois)) / np.sqrt(len(variance206pois)-1)
            pb206SEpercentpois = pb206SEpois / np.mean(fit206pois.predict()) * 100
            pb204SEzero = np.sqrt(np.mean(variance204zero))/np.sqrt(len(variance204zero)-1)
            pb204SEpercentzero = pb204SEzero / np.mean(fit204zero.predict()) * 100

            pb207_206 = pb207zero/pb206pois
            pb206_204 = pb206pois/pb204zero
            u238_235 = data_ratios_['238U/235U'].mean()
            
            u238_235SE = data_ratios_['238U/235U'].sem()
            u238_235SE_percent = data_ratios_['238U/235U'].sem()/u238_235*100
            pb207_206SE = pb207_206*((pb207SEpercentzero/100)**2+(pb206SEpercentpois/100)**2)**(1/2)
            pb207_206SE_percent = pb207_206SE/pb207_206*100
            pb206_204SE = pb206_204*((pb206SEpercentpois/100)**2+(pb204SEpercentzero/100)**2)**(1/2)
            pb206_204SE_percent = pb206_204SE/pb206_204*100
            
            ratios_to_return = [pb206_204,pb207_206,u238_235,pb207zero,pb206pois,pb204zero]
            stats_to_return = [pb207SEzero,pb207SEpercentzero,pb206SEpois,pb206SEpercentpois,pb204SEzero,pb204SEpercentzero,u238_235SE,u238_235SE_percent,
                               pb207_206SE,pb207_206SE_percent,pb206_204SE,pb206_204SE_percent]

            return ratios_to_return,stats_to_return
            
        elif counts_mode == 'Means & Regression':
            data_ratios_ = calc_fncs.get_ratios(data)
            pb206_204 = data_ratios_['206Pb/204Pb'].mean()
            pb207_206 = data_ratios_['207Pb/206Pb'].mean()
            u238_235 = data_ratios_['238U/235U'].mean()
            
            pb206_204SE = data_ratios_['206Pb/204Pb'].sem()
            pb206_204SE_percent = data_ratios_['206Pb/204Pb'].sem()/pb206_204*100
            pb207_206SE = data_ratios_['207Pb/206Pb'].sem()
            pb207_206SE_percent = data_ratios_['207Pb/206Pb'].sem()/pb207_206*100
            u238_235SE = data_ratios_['238U/235U'].sem()
            u238_235SE_percent = data_ratios_['238U/235U'].sem()/u238_235*100
            
            ratios_to_return = [pb206_204,pb207_206,u238_235]
            stats_to_return = [pb206_204SE,pb206_204SE_percent,pb207_206SE,pb207_206SE_percent,u238_235SE,u238_235SE_percent]
            
            return ratios_to_return,stats_to_return
        
        
    
    
    def get_regressions(data,regression_buttons,ablation_start_true,force_exp_params,joggle_params,aguess,bguess,aguess_207,bguess_207):
        """
        function that gets the 206Pb/238U regression

        Parameters
        ----------
        data : dataframe
            pandas dataframe hosting the measured data.
        regression_buttons : string
            string object defining which type of regression to use.
        ablation_start_true : float
            float of where to project the ratio back to.
        force_exp_params : Boolean
            Allows user to provide the default initial guess for the exponential regression.
        joggle_params : Boolean
            Allows user to override the default initial guess for the exponential regression and define their own.
        aguess : float
            input initial a parameter to optimize from. This is the intercept
        bguess : float
            input initial b parameter to optimize from. This is effectively the slope
        cguess : float
            input initial c parameter to optimize from. This is an additive term

        Returns
        -------
        Returned objects depend on the function that calls get_regressions
            get_approved: returns regressed 206/238 and error
            ratio_plot: returns regression to plot
            get_regression stats: returns regression statistics to display them to user
            residuals_plot: returns regression residuals to plot

        """
        y = data['206Pb/238U']
        y207 = data['207Pb/235U']
        if '1st Order' in regression_buttons:
            y1, X1 = dmatrices('y ~ Time_s', data=data, return_type='dataframe') # get y and x regression data
            mod1 = sm.OLS(y1, X1) # fit a linear model on the data
            fit1 = mod1.fit() # get the list of fit parameters
            predicted1 = fit1.params[0] + fit1.params[1]*data.Time_s # get the predicted y values for the given x values
            rsquared1 = fit1.rsquared # get the R2 value of the regression
            predicted_b01 = fit1.params[0] + fit1.params[1]*ablation_start_true # get the predicted value at the ablation start that is input by the user
            sigma1 = np.sqrt(fit1.ssr/fit1.df_resid) # get the 1SD (Sum Squared residuals / residual degrees of freedom)^(1/2)
            # get standard error for a single point estiamted by a regression model
            SE_b01 = sigma1*np.sqrt(1/len(data)+(ablation_start_true-data['Time_s'].mean())**2/((len(data)-1)*statistics.variance(data['Time_s'])))
            SE_b01_percent = SE_b01/predicted_b01*100 # get the % 1SE
            resid1 = fit1.resid # get the residuals of the regression
            
            y1_207,X1_207 = dmatrices('y207 ~ Time_s', data=data, return_type='dataframe') # get y and x regression data
            fit1_207 = sm.OLS(y1_207,X1_207).fit()
            predicted1_207 = fit1_207.params[0] + fit1_207.params[1]*data.Time_s # get the predicted y values for the given x values
            predicted_b01_207 = fit1_207.params[0] + fit1_207.params[1]*ablation_start_true # get the predicted value at the ablation start that is input by the user
            sigma1_207 = np.sqrt(fit1_207.ssr/fit1_207.df_resid) # get the 1SD (Sum Squared residuals / residual degrees of freedom)^(1/2)
            SE_b01_207 = sigma1_207*np.sqrt(1/len(data)+(ablation_start_true-data['Time_s'].mean())**2/((len(data)-1)*statistics.variance(data['Time_s'])))
            SE_b01_percent_207 = SE_b01_207/predicted_b01_207*100 # get the % 1SE
            resid1_207 = fit1_207.resid # get the residuals of the regression
        else:
            # fill the above with blank values so that there is not an error in the output if '1st Order' is not wanted by the user
            predicted1 = np.zeros_like(data['Time_s'])
            rsquared1 = 0
            predicted_b01 = 0
            SE_b01 = 0
            SE_b01_percent = 0
            resid1 = np.zeros_like(data['Time_s'])
            
            predicted1_207 = np.zeros_like(data['Time_s'])
            predicted_b01_207 = 0
            SE_b01_207 = 0
            SE_b01_percent_207 = 0
            resid1_207 = np.zeros_like(data['Time_s'])
            
        # run an if loop as above. Virtually the same except we are using a 2nd order regression
        if '2nd Order' in regression_buttons:
            y2, X2 = dmatrices('y ~ Time_s + Time_s2', data=data, return_type='dataframe')
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
            
            y2_207, X2_207 = dmatrices('y207 ~ Time_s + Time_s2', data=data, return_type='dataframe')
            fit2_207 = sm.OLS(y2_207,X2_207).fit()
            predicted2_207 = fit2_207.params[0] + fit2_207.params[1]*data.Time_s + fit2_207.params[2]*data.Time_s2
            predicted_b02_207 = fit2_207.params[0] + fit2_207.params[1]*ablation_start_true + fit2_207.params[2]*ablation_start_true
            sigma2_207 = np.sqrt(fit2_207.ssr/fit2_207.df_resid)
            SE_b02_207 = sigma2_207*np.sqrt(1/len(data)+(ablation_start_true-data['Time_s'].mean())**2/((len(data)-1)*statistics.variance(data['Time_s'])))
            SE_b02_percent_207 = SE_b02_207/predicted_b02_207*100
            resid2_207 = fit2_207.resid
        else:
            predicted2 = np.zeros_like(data['Time_s'])
            rsquared2 = 0
            rsquared2_adj = 0
            predicted_b02 = 0
            SE_b02 = 0
            SE_b02_percent = 0
            resid2 = np.zeros_like(data['Time_s'])
            
            predicted2_207 = np.zeros_like(data['Time_s'])
            predicted_b02_207 = 0
            SE_b02_207 = 0
            SE_b02_percent_207 = 0
            resid2_207 = np.zeros_like(data['Time_s'])
        
        if 'Exp. Regression' in regression_buttons:
            def exp_func(x,a,b):
                return a*np.exp(-b*x)
            if force_exp_params is True:
                if joggle_params is True:
                    popt,pcov = curve_fit(exp_func,data['Time_s'],data['206Pb/238U'],p0=(aguess,bguess))
                    popt_207,pcov_207 = curve_fit(exp_func,data['Time_s'],data['207Pb/235U'],p0=(aguess_207,bguess_207))
                elif joggle_params is False:
                    popt,pcov = curve_fit(exp_func,data['Time_s'],data['206Pb/238U'],p0=(0.1,0.05))
                    popt_207,pcov_207 = curve_fit(exp_func,data['Time_s'],data['207Pb/235U'],p0=(0.01,0.05))
            else:
                popt,pcov = curve_fit(lambda t, a, b: exp_func, data['Time_s'],data['206Pb/238U'])
                popt_207,pcov_207 = curve_fit(lambda t, a, b: exp_func, data['Time_s'],data['207Pb/235U'])
            predictedexp = exp_func(data['Time_s'],*popt)
            predictedexp_207 = exp_func(data['Time_s'],*popt_207)
            predicted_b0exp = popt[0]*np.exp(-popt[1]*ablation_start_true)
            predicted_b0exp_207 = popt_207[0]*np.exp(-popt_207[1]*ablation_start_true)
            resid = []
            sq_resid = []
            resid_207 = []
            sq_resid_207 = []
            for i,k,m,l in zip(predictedexp,data['206Pb/238U'],predictedexp_207,data['207Pb/235U']):
                resid.append((k - i))
                sq_resid.append(((k - i)**2))
                resid_207.append((l-m))
                sq_resid_207.append((l-m)**2)
            sum_sq_resid = np.sum(sq_resid)
            sum_sq_resid_207 = np.sum(sq_resid_207)
            sigmaexp = np.sqrt(sum_sq_resid/(len(data['206Pb/238U'])-2)) # denominator = d.f. = n-#params
            sigmaexp_207 = np.sqrt(sum_sq_resid_207/(len(data['207Pb/235U'])-2)) # denominator = d.f. = n-#params
            SE_b0exp = sigmaexp*np.sqrt(1/len(data)+(ablation_start_true-data['Time_s'].mean())**2/((len(data)-1)*statistics.variance(data['Time_s'])))
            SE_b0exp_207 = sigmaexp_207*np.sqrt(1/len(data)+(ablation_start_true-data['Time_s'].mean())**2/((len(data)-1)*statistics.variance(data['Time_s'])))
            SE_b0exp_percent = SE_b0exp/predicted_b0exp
            SE_b0exp_percent_207 = SE_b0exp_207/predicted_b0exp_207
            residexp = resid
            residexp_207 = resid_207
            tss = ((data['206Pb/238U'] - np.mean(data['206Pb/238U']))**2).sum()
            rsquared_exp = 1 - (sum_sq_resid/tss)
            rsquaredexp_adj = 0
        else:
            predictedexp = np.zeros_like(data['Time_s'])
            rsquared_exp = 0
            rsquaredexp_adj = 0
            predicted_b0exp = 0
            SE_b0exp = 0
            SE_b0exp_percent = 0
            residexp = np.zeros_like(data['Time_s'])
            predictedexp_207 = np.zeros_like(data['Time_s'])
            predicted_b0exp_207 = 0
            SE_b0exp_207 = 0
            SE_b0exp_percent_207 = 0
            residexp_207 = np.zeros_like(data['Time_s'])
    
        # get the method that called up regressions. f_back gets the function that called. Removing this gives current method
        callingmethod = sys._getframe().f_back.f_code.co_name
        # set a series of if statements that causes the appropriate return depending on the function that called up regresssions
        if callingmethod == 'get_approved':
            
            # put the calculated values and statistics in lists to be returned
            ratios_to_return = [predicted_b01,predicted_b02,predicted_b0exp,predicted_b01_207,predicted_b02_207,predicted_b0exp_207]
            stats_to_return = [rsquared1,rsquared2,rsquared_exp,rsquared2_adj,rsquaredexp_adj,
                               SE_b01,SE_b01_percent,SE_b02,SE_b02_percent,SE_b0exp,SE_b0exp_percent,
                               SE_b01_207,SE_b01_percent_207,SE_b02_207,SE_b02_percent_207,SE_b0exp_207,SE_b0exp_percent_207]
            
            return ratios_to_return,stats_to_return
        
        elif callingmethod == 'ratio_plot':
            regressions_to_return = [predicted1,predicted2,predictedexp,predicted1_207,predicted2_207,predictedexp_207]
            stats_to_report = [rsquared1,rsquared2,rsquared2_adj,rsquaredexp_adj,SE_b01_percent,SE_b02_percent,SE_b0exp_percent]
            
            return regressions_to_return,stats_to_report
        
        elif callingmethod == 'get_regression_stats':
            # ""
            regressions_to_return = [predicted1,predicted2,predictedexp]
            stats_to_report = [rsquared1,rsquared2,rsquared_exp,SE_b01_percent,SE_b02_percent,SE_b0exp_percent,SE_b01_percent_207,SE_b02_percent_207,SE_b0exp_percent_207]
            
            return regressions_to_return,stats_to_report
        
        elif callingmethod == 'residuals_plot':
            predicted_to_return = [predicted1,predicted2,predictedexp]
            predicted_to_return_207 = [predicted1_207,predicted2_207,predictedexp_207]
            resid_to_return = [resid1,resid2,residexp]
            resid_to_return_207 = [resid1_207,resid2_207,residexp_207]
            
            return predicted_to_return,predicted_to_return_207,resid_to_return,resid_to_return_207
        
        else:
            pass
        
        
        
    def get_approved(data,bckgrnd_start_input,bckgrnd_stop_input,
                     ablation_start_input,ablation_stop_input,ablation_start_true,
                    regression_buttons,ellipsemode_selector,force_exp_params,joggle_params,aguess,bguess,aguess_207,bguess_207,
                    counts_mode,arrayofdwelltimes,sample_names):
        """

        Parameters
        ----------
        data : dataframe
            pandas dataframe hosting the measured data.
        bckgrnd_start_input : integer
            integer from the slider defining when to start taking the gas blank.
        bckgrnd_stop_input : integer
            integer from the slider defining when to stop taking the gas blank.
        ablation_start_input : integer
            integer from the slider defining when to start the regression / counts.
        ablation_stop_input : integer
            integer form the slider defining when to stop the regression / counts.
        ablation_start_true : float
            float value of the projected ablation start time.
        regression_buttons : string
            string object defining which type of regression to use.
        ellipsemode_selector : Boolean
            Object defining weather or not to include confidence ellipses in output.
        force_exp_params : Boolean
            Allows user to provide the default initial guess for the exponential regression.
        joggle_params : Boolean
            Allows user to override the default initial guess for the exponential regression and define their own.
        aguess : float
            input initial a parameter to optimize from. This is the intercept
        bguess : float
            input initial b parameter to optimize from. This is effectively the slope
        cguess : TYPE
            input initial c parameter to optimize from. This is an additive term
        counts_mode : string
            string denoting which method the using wants to reduce data with.
        arrayofdwelltimes : array
            Array of dwell times in s.
            
        Returns
        -------
        data_approved : array
            arary of all reduced ratios and analytes.
        ellipse_data_toapprove : dataframe
            dataframe holding all the passes of the detector across the selected ablation period.

        """
        data_toapprove = data.reset_index(drop=True) # reset the index of the data so that it can be altered/indexed appropriately
        og_col_length = len(data_toapprove.columns[2:])+1 # get original column length
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
        data_toapprove = data_toapprove[(data_toapprove.Time_s >= ablation_start_input) & (data_toapprove.Time_s <= ablation_stop_input)]
        
        
        if counts_mode == 'Total Counts':
            data_ratios_all,data_stats_all = calc_fncs.get_counts(data_toapprove, counts_mode, arrayofdwelltimes, ablation_start_true, bckgrnd_start_input, bckgrnd_stop_input, ablation_start_input, ablation_stop_input)
            data_ratios = pd.DataFrame([data_ratios_all],columns=['206/238 Exp.','207Pb/235U','206Pb/204Pb','207Pb/206Pb','238U/235U'])
            data_stats = pd.DataFrame([data_stats_all],columns=['SE 206/238 Exp','SE% 206/238 Exp','SE 207/235','SE% 207/235','SE 206/204','SE% 206/204','SE 207/206','SE% 207/206',
                                                                'SE 238/235','SE% 238/235'])
            
        elif counts_mode == 'Poisson':
            data_approved_ratio = data_toapprove.copy() # copy the dataframe as to not overwrite the input data
            data_approved_ratio = calc_fncs.get_ratios(data_approved_ratio) # get calculated ratios of the background subtracted data. This is why we need 0 and not 'bdl'
            # index the background subtracted data for the input ablation periods.
            data_approved_ratio = data_approved_ratio[(data_approved_ratio.Time_s >= ablation_start_input) & (data_approved_ratio.Time_s <= ablation_stop_input)]
            # ellipse_data_approved = data_approved_ratio
            data_ratios_reg,data_stats_reg = calc_fncs.get_regressions(data_approved_ratio,regression_buttons,ablation_start_true,force_exp_params,joggle_params,aguess,bguess,aguess_207,bguess_207) # get estimated regression intercepts, background subtracted ratios, regression statistics, and SE's of the ratios3
            data_ratios_counts,data_stats_counts = calc_fncs.get_counts(data_toapprove, counts_mode, arrayofdwelltimes, ablation_start_true, bckgrnd_start_input, bckgrnd_stop_input, ablation_start_input, ablation_stop_input)
            data_ratios = data_ratios_reg + data_ratios_counts
            data_stats = data_stats_reg + data_stats_counts            
            data_ratios = pd.DataFrame([data_ratios],columns=['206/238 1st Order','206/238 2nd Order','206/238 Exp.','207Pb/235U 1st Order','207Pb/235U 2nd Order','207Pb/235U Exp.','206Pb/204Pb','207Pb/206Pb','238U/235U','207Pb counts','206Pb counts','204Pb counts']) # put the relevant calculations in df with appropriate column headers
            # put the relevant calculations in df with appropriate column headers
            data_stats = pd.DataFrame([data_stats],columns=['R2 1st Order','R2 2nd Order','R2 Exp','R2 Adj 2nd Order','R2 Adj Exp','SE 206/238 1st Order','SE% 206/238 1st Order','SE 206/238 2nd Order','SE% 206/238 2nd Order',
                                                            'SE 206/238 Exp','SE% 206/238 Exp','SE 207/235 1st Order','SE% 207/235 1st Order','SE 207/235 2nd Order','SE% 207/235 2nd Order','SE 207/235 Exp','SE% 207/235 Exp',
                                                            'SE 207Pb','SE% 207Pb','SE 206Pb','SE% 206Pb','SE 204Pb','SE% 204Pb','SE 238/235','SE% 238/235','SE 207/206','SE% 207/206','SE 206/204','SE% 206/204'])
            
        elif counts_mode == 'Means & Regression':
            data_approved_ratio = data_toapprove.copy() # copy the dataframe as to not overwrite the input data
            data_approved_ratio = calc_fncs.get_ratios(data_approved_ratio) # get calculated ratios of the background subtracted data. This is why we need 0 and not 'bdl'
            # index the background subtracted data for the input ablation periods.
            data_approved_ratio = data_approved_ratio[(data_approved_ratio.Time_s >= ablation_start_input) & (data_approved_ratio.Time_s <= ablation_stop_input)]
            # ellipse_data_approved = data_approved_ratio
            data_ratios_reg,data_stats_reg = calc_fncs.get_regressions(data_approved_ratio,regression_buttons,ablation_start_true,force_exp_params,joggle_params,aguess,bguess,aguess_207,bguess_207) # get estimated regression intercepts, background subtracted ratios, regression statistics, and SE's of the ratios3
            data_ratios_counts,data_stats_counts = calc_fncs.get_counts(data_toapprove, counts_mode, arrayofdwelltimes, ablation_start_true, bckgrnd_start_input, bckgrnd_stop_input, ablation_start_input, ablation_stop_input)
            data_ratios = data_ratios_reg + data_ratios_counts
            data_stats = data_stats_reg + data_stats_counts
            data_ratios = pd.DataFrame([data_ratios],columns=['206/238 1st Order','206/238 2nd Order','206/238 Exp.','207Pb/235U 1st Order','207Pb/235U 2nd Order','207Pb/235U Exp.','206Pb/204Pb','207Pb/206Pb','238U/235U']) # put the relevant calculations in df with appropriate column headers
            # put the relevant calculations in df with appropriate column headers
            data_stats = pd.DataFrame([data_stats],columns=['R2 1st Order','R2 2nd Order','R2 Exp','R2 Adj 2nd Order','R2 Adj Exp','SE 206/238 1st Order','SE% 206/238 1st Order','SE 206/238 2nd Order','SE% 206/238 2nd Order',
                                                            'SE 206/238 Exp','SE% 206/238 Exp','SE 207/235 1st Order','SE% 207/235 1st Order','SE 207/235 2nd Order','SE% 207/235 2nd Order','SE 207/235 Exp','SE% 207/235 Exp',
                                                            'SE 206/204','SE% 206/204','SE 207/206','SE% 207/206','SE 238/235','SE% 238/235'])
            
        # get the ablation period indexed, background subtracted, indvidual isotope measurements
        data_toapprove = data_toapprove[(data_toapprove.Time_s >= ablation_start_input) & (data_toapprove.Time_s <= ablation_stop_input)]
        if ellipsemode_selector is True:
            ellipse_data_toapprove = data_toapprove.copy()
            ellipse_data_toapprove = ellipse_data_toapprove.reset_index(drop=True)
        data_toapprove_SE = pd.DataFrame([data_toapprove.iloc[:,2:-1].sem()]).add_suffix('_1SE') # get SE's of the individual isotopes
        data_toapprove = pd.DataFrame([data_toapprove.iloc[:,2:-1].mean()]) # put means of isotopic values in a dataframe
        
        # 202Hg: 29.86%, 204Hg: 6.87% https://www.nndc.bnl.gov/nudat3/
        Hgratio = (6.87/100) / (29.86/100)
        
        if ellipsemode_selector is True:
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
            ellipse_data_toapprove.insert(0,'measurementindex',sample_names.index(data.iloc[0,0]))
                        
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
            data_toapprove['204Hg'] = 'bdl'
        
        data_toapprove.insert(0,'measurementindex',sample_names.index(data.iloc[0,0]))
        data_toapprove.insert(1,'SampleLabel',data.iloc[0,0]) # reinsert the sample label into the calculation df
        data_toapprove.insert(2,'t start',[ablation_start_input]) # insert ablation start time into df
        data_toapprove.insert(3,'t end',[ablation_stop_input]) # insert ablation stop time into df
        data_toapprove.insert(4,'t project',[ablation_start_true]) # insert projected regression start time into df
        # stitch the individual isotopic measurements, their errors, ratios, ratio errors, and regression results / regression statistics into a df together.
        # these are then appeneded into the output df
        data_approved = data_toapprove.join([data_toapprove_SE,data_ratios,data_stats])
        
        # data_approved = data_toapprove.join([data_toapprove_SE,data_ratios_reg,data_ratios_counts,data_stats_reg,data_stats_counts])
        
        ratio_cols = ['206Pb/204Pb','207Pb/206Pb','238U/235U']
        for i in analyte_cols:
            for k in ratio_cols:
                if k.__contains__(str(i)) and data_approved[i].item() == 'bdl':
                    data_approved[k] = 'bdl'
                else:
                    pass
                
        if ellipsemode_selector is True:
            ellipse_data_toapprove_ratio = calc_fncs.get_ratios(ellipse_data_toapprove).drop(['Time_s2','204Hg'],axis=1)
            # ellipse_data_toapprove_ratio.rename({'Pb206_Pb204':'206Pb/204Pb' , 'Pb206_U238':'206Pb/238U', 'Pb207_U235':'207Pb/235U', 'Pb207_Pb206':'207Pb/206Pb', 'U235_U238':'235U/238U', 'U238_U235':'238U/235U'}, axis=1, inplace=True)
            ellipse_data_toapprove = pd.concat([ellipse_data_toapprove,ellipse_data_toapprove_ratio],axis=1)
            
            ratio_colsII = ['206Pb/204Pb','207Pb/206Pb','238U/235U','207Pb/235U 1st Order','207Pb/235U 2nd Order','207Pb/235U Exp.']
            for i in analyte_cols:
                for k in ratio_colsII:
                    if k.__contains__(str(i)):
                        for j in range(0,len(ellipse_data_toapprove)):
                            if ellipse_data_toapprove.loc[j,i] == 'bdl':
                                ellipse_data_toapprove.loc[j,k] = 'bdl'
        else:
            ellipse_data_toapprove = 0
                    
        
        return data_approved,ellipse_data_toapprove
    
    def get_ellipse(data,power):
        """
        Function that gets the confidence ellipses

        Parameters
        ----------
        data : dataframe
            pandas dataframe holding the observed data.
        power : float
            float object defining the power of the confidence. Recommended to keep this at 0.05

        Returns
        -------
        ell1_params : array
            array of float objects defining the dimensions of the confidence ellipse for Wetherhill concordia.
        ell2_params : array
            array of float objects defining the dimensions of the confidence ellipse for Tera-Wasserburg concordia.

        """
        data = data.dropna()
        drop_condn = data[(data['206Pb/238U'] == 0) | (data['207Pb/235U'] == 0) | (data['207Pb/206Pb'] == 0)].index
        data.drop(drop_condn,inplace=True)
        x1 = data['207Pb/235U']
        y1 = data['206Pb/238U']
        x2 = 1/data['206Pb/238U']
        y2 = data['207Pb/206Pb']
        
        cov1 = np.cov(x1,y1)
        cov2 = np.cov(x2,y2)
        eigval1,eigvec1 = np.linalg.eig(cov1)
        eigval2,eigvec2 = np.linalg.eig(cov2)
        order1 = eigval1.argsort()[::-1]
        order2 = eigval2.argsort()[::-1]
        eigvals_order1 = eigval1[order1]
        eigvals_order2 = eigval2[order2]
        eigvecs_order1 = eigvec1[:,order1]
        eigvecs_order2 = eigvec2[:,order2]
        
        c1 = (np.mean(x1),np.mean(y1))
        c2 = (np.mean(x2),np.mean(y2))
        wid1 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power),df=2)*eigvals_order1[0])
        hgt1 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power),df=2)*eigvals_order1[1])
        wid2 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power),df=2)*eigvals_order2[0])
        hgt2 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power),df=2)*eigvals_order2[1])
        theta1 = np.degrees(np.arctan2(*eigvecs_order1[:,0][::-1]))
        theta2 = np.degrees(np.arctan2(*eigvecs_order2[:,0][::-1]))
        
        ell1_params = [c1,wid1,hgt1,theta1]
        ell2_params = [c2,wid2,hgt2,theta2]
        
        return ell1_params,ell2_params
        
        
class plots(calc_fncs):
    """ Class that holds all of the functions for reducing the time resolved data"""
    def __init__(self,*args):
        super().__init__(*args)
        for a in args:
            self.__setattr__(str(a), args[0])
            
    
    def ratio_plot(data,ratio_buttons,regression_buttons,background_slider,ablation_slider,ablation_start_true,
                   ratio_plot_ylim_slider,
                   force_exp_params,joggle_params,aguess,bguess,aguess_207,bguess_207,
                   ):
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
        fig = figure(height=500,width=800,title='Ratios',tools='pan,reset,save,wheel_zoom,xwheel_zoom,ywheel_zoom',toolbar_location='above',
                     x_axis_label='Time (s)',y_axis_label='Isotope Ratios',
                     y_range=[ratio_plot_ylim_slider[0],ratio_plot_ylim_slider[1]],x_range=[min(data.Time_s),max(data.Time_s)])
        fig.xgrid.grid_line_color = None
        fig.ygrid.grid_line_color = None
        data = calc_fncs.get_ratios(data) # get calculated ratios from the data
        data_to_regress = data[(data.Time_s >= ablation_slider[0]) & (data.Time_s <= ablation_slider[1])] # get the selected ablation period
        var_cols = ratio_buttons # get the ratios selected by the user. These are used to get the columns from the calculated ratios
        regressions_to_plot = regression_buttons# get the regression type requested by the user
        # get regression parameters and stats for the requested regressions across the specified ablation period
        regressions,stats = calc_fncs.get_regressions(data_to_regress,regression_buttons,ablation_start_true,force_exp_params,
                                                      joggle_params,aguess,bguess,aguess_207,bguess_207)
        # plot a line for each selected ratio
        for i,c in zip(var_cols,cycle(color_palette)):
            fig.line(data.Time_s,data[i],line_width=1,legend_label='{}'.format(i),color=c)
            # print(data[i])
        # plot a line for each selected regression
        if regressions_to_plot is not None:
            for i,c in zip(regressions,cycle(color_palette_regressions)):
                fig.line(data_to_regress.Time_s,i,line_width=0.5,color=c)
        fig.line([background_slider[0],background_slider[0]],[0,1],line_width=0.4,line_dash='dashed',color='black')
        fig.line([background_slider[1],background_slider[1]],[0,1],line_width=0.4,line_dash='dashed',color='black')
        fig.line([ablation_slider[0],ablation_slider[0]],[0,1],line_width=0.4,line_dash='solid',color='black')
        fig.line([ablation_slider[1],ablation_slider[1]],[0,1],line_width=0.4,line_dash='solid',color='black')
        fig.line([ablation_start_true,ablation_start_true],[0,1],line_width=0.4,line_dash='dotted',color='black')
        
        return fig
    
    
    def ablation_plot(data_ablation_,background_slider,ablation_slider,ablation_start_true,ablation_plot_ylim_slider,logdata,analytes_):  
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
        # # ax.tick_params(axis='both',labelsize=14) # plot aesthetics
        if logdata == True:
            y_type = 'log'
        else:
            y_type = 'auto'
        fig = figure(height=400,width=1000,title='All Time Resolved Data',tools='pan,reset,save,wheel_zoom,xwheel_zoom,ywheel_zoom',toolbar_location='left',
                      y_axis_type=y_type,x_axis_label = 'Time (s)', y_axis_label = 'Intensities (cps)',
                      y_range=[ablation_plot_ylim_slider[0],ablation_plot_ylim_slider[1]],
                      x_range=[min(data_ablation_.Time_s),max(data_ablation_.Time_s)]
                      )
        
        var_cols=analytes_
        for i,c in zip(var_cols,cycle(color_palette)):
            fig.line(data_ablation_.Time_s,data_ablation_[i],line_width=0.7,legend_label='{}'.format(i),color=c)
        fig.line([background_slider[0],background_slider[0]],[0,1e8],line_width=0.4,line_dash='dashed',color='black')
        fig.line([background_slider[1],background_slider[1]],[0,1e8],line_width=0.4,line_dash='dashed',color='black')
        fig.line([ablation_slider[0],ablation_slider[0]],[0,1e8],line_width=0.4,line_dash='solid',color='black')
        fig.line([ablation_slider[1],ablation_slider[1]],[0,1e8],line_width=0.4,line_dash='solid',color='black')
        fig.line([ablation_start_true,ablation_start_true],[0,1e8],line_width=0.4,line_dash='dotted',color='black')
        
        return fig
        
        
    def residuals_plot(data,regression_buttons,start_ablation,stop_ablation,ablation_start_true,force_exp_params,joggle_params,aguess,bguess,aguess_207,bguess_207):
        fig206 = figure(height=250,width=250,title='206Pb/238U Residuals',tools='pan,reset,save,wheel_zoom,xwheel_zoom,ywheel_zoom',toolbar_location='above',
                        x_axis_label='Fitted Value',y_axis_label='Residuals')
        fig206.legend.label_text_font_size='8pt'
        fig207 = figure(height=250,width=250,title='207Pb/235U Residuals',tools='pan,reset,save,wheel_zoom,xwheel_zoom,ywheel_zoom',toolbar_location='above',
                        x_axis_label='Fitted Value',y_axis_label='Residuals')
        fig207.legend.label_text_font_size='8pt'
        data = calc_fncs.get_ratios(data) # calculate relevant isotopic ratios from the data
        data_to_regress = data[(data.Time_s >= start_ablation) & (data.Time_s <= stop_ablation)] # get data across the requested ablation interval
        fitted,fitted207,residuals,residuals207 = calc_fncs.get_regressions(data_to_regress,regression_buttons,ablation_start_true,force_exp_params,joggle_params,
                                                                            aguess,bguess,aguess_207,bguess_207) # get fitted regression values, regression names, and residuals
        # plot the residuals for the regressions and the fitted value
        for j,c in zip(range(0,len(regression_buttons)),color_palette_regressions):
            fig206.circle(fitted[j],residuals[j],color=c,legend_label='{}'.format(regression_buttons[j]))
            fig207.circle(fitted207[j],residuals207[j],color=c,legend_label='{}'.format(regression_buttons[j]))
        fig206.line([min(fitted[j]),max(fitted[j])],[0,0],color='black',line_width=0.4)
        fig207.line([min(fitted207[j]),max(fitted207[j])],[0,0],color='black',line_width=0.4)
        fig = gridplot([[fig206],[fig207]],width=500,height=250)
        
        return fig
    
    def ellipse_plot(data,power,start_ablation,stop_ablation):
        data = calc_fncs.get_ratios(data)
        data = data[(data.Time_s >= start_ablation) & (data.Time_s <= stop_ablation)]
        data = data.dropna()
        drop_condn = data[(data['206Pb/238U'] == 0) | (data['207Pb/235U'] == 0) | (data['207Pb/206Pb'] == 0)].index
        data.drop(drop_condn,inplace=True)
        # ell1p,ell2p,ell3p = calc_fncs.get_ellipse(data, power)
        ell1p,ell2p = calc_fncs.get_ellipse(data, power)
        figell1 = figure(height=400,width=400,title='Confidence Ellipse',toolbar_location=None,
                         y_axis_label='206Pb/238U',x_axis_label='207Pb/235U')
        figell2 = figure(height=400,width=400,title='Confidence Ellipse',toolbar_location=None,
                         y_axis_label='207Pb/206Pb',x_axis_label='238U/206Pb')
        figell1.dot(data['207Pb/235U'],data['206Pb/238U'],color='black',size=8)
        figell2.dot(1/data['206Pb/238U'],data['207Pb/206Pb'],color='black',size=8)
        figell1.ellipse(x=ell1p[0][0],y=ell1p[0][1],width=ell1p[1],height=ell1p[2],angle=ell1p[3],fill_color=color_palette_regressions[0],fill_alpha=0.5)
        figell2.ellipse(x=ell2p[0][0],y=ell2p[0][1],width=ell2p[1],height=ell2p[2],angle=ell2p[3],fill_color=color_palette_regressions[0],fill_alpha=0.5)
        return figell1,figell2
    

        
class make_plots(param.Parameterized):
    """ class that parameterizes inputs and sends them to the above functions to be rendered in a GUI"""
    sample_subset = param.Selector(objects=[])
    
    update_output_button = param.Action(lambda x: x.add_output_data(),label='Approve Interval')
    export_data_button = param.Action(lambda x: x.export_data(),label='DDDT!')
    export_plots_button = param.Action(lambda x: x.export_plots(),label ='Export All Plots')

    ablation_start_true = param.Number(31.5)
    ablation_slider = param.Range(default=(32,60),bounds=(0,100))
    background_slider = param.Range(default=(5,29),bounds=(0,100))
    
    ratio_buttons = param.ListSelector(default=['206Pb/238U'], objects=['206Pb/238U','206Pb/204Pb','207Pb/235U','207Pb/206Pb','238U/235U'])
    
    regression_buttons = param.ListSelector(default=['1st Order'], objects=['1st Order','2nd Order','Exp. Regression'])
    force_exp_params = param.Boolean(False,label='Force Exp. Params')
    joggle_params = param.Boolean(False,label='Joggle Exp. Params')
    aguess = param.Number(default=0.001)
    bguess = param.Number(default=-0.08)
    aguess_207 = param.Number(default=0.0001)
    bguess_207 = param.Number(default=-0.01)
    
    ablation_plot_ylim_slider = param.Range(default=(0,1000000),bounds=(0,10000000))
    ratio_plot_ylim_slider = param.Range(default=(0.0,0.3),bounds=(-0.05,3))
    # ratio_plot_xlim_slider = param.Range(default=(0,60),bounds=(0,100))
    
    logcountsdata = param.Boolean(False,label='Log Intensities')
    analytes_ = param.ListSelector(default=['202Hg','204Pb','206Pb','207Pb','208Pb','232Th','235U','238U'],objects=['202Hg','204Pb','206Pb','207Pb','208Pb','232Th','235U','238U'])
    
    counts_mode = param.Selector(default='Means & Regression',objects=['Total Counts','Poisson','Means & Regression'])
    
    ellipsemode_selector = param.Boolean(True,label='Generate Ellipse')
    power = param.Number(default=0.05)
    
    arrayofdwelltimes = param.Array()
    
    input_data = param.DataFrame(precedence=-1)
    file_path = param.String(default='Insert File Path')
    output_data = param.DataFrame(precedence=-1)
    output_data_ellipse = param.DataFrame(precedence=-1)
    
        
    def __init__(self,**params):
        super().__init__(**params)
        self.file_input_widget = pn.Param(self.param.input_data)
        self.output_data_widget = pn.Param(self.param.output_data)
        self.output_data_ellipse_widget = pn.Param(self.param.output_data_ellipse)
        self.widgets = pn.Param(self,parameters=['sample_subset','update_output_button','export_data_button',
                                                 'ablation_start_true','ablation_slider','background_slider','ablation_plot_ylim_slider','logcountsdata','analytes_',
                                                 'ratio_plot_ylim_slider','ratio_plot_ylim_max','ratio_buttons','regression_buttons','counts_mode',
                                                 'ellipsemode_selector','power',
                                                 'arrayofdwelltimes','file_path'])
    
    @pn.depends('file_path',watch=True)
    def _uploadfile(self):
        if self.file_path != 'Insert File Path':
            df = pd.read_excel(self.file_path,sheet_name='Buffer')
            self.input_data = df
            self.input_data['Time_s'] = self.input_data['Time']/1000
            options = list(self.input_data.SampleLabel.unique())
            self.param.sample_subset.objects = options
            self.output_data = pd.DataFrame([np.zeros(len(self.input_data.columns))],columns=list(self.input_data.columns))
            self.output_data.insert(0,'measurementindex',0)
            self.output_data.insert(2,'t start',0)
            self.output_data.insert(3,'t end',0)
            self.output_data.insert(4,'t project',0)
            self.output_data.drop('Time',axis=1)
            # self.output_data.iloc[0,1:3] = pd.DataFrame(np.zeros(3),columns=['t start','t end','t project'])
            if self.ellipsemode_selector is True:
                self.output_data_ellipse = pd.DataFrame([np.zeros(len(self.input_data.columns))],columns=list(self.input_data.columns))
                self.output_data_ellipse.insert(0,'measurementindex',0)
            else:
                pass
            analytelength = len(self.input_data.columns)-3
            dwellarray = np.full_like(np.arange(analytelength,dtype=float),0.01)
            self.arrayofdwelltimes = dwellarray
        else:
            return print('No Input Data Available')
    
        
    @pn.depends('input_data','background_slider','ablation_slider','ablation_start_true','sample_subset','ablation_plot_ylim_slider')
    def call_ablation_plot(self):
        if self.sample_subset is not None:
            data_toplot = self.input_data[self.input_data['SampleLabel'] == self.sample_subset]
            return plots.ablation_plot(data_toplot,self.background_slider,self.ablation_slider,self.ablation_start_true,
                                       self.ablation_plot_ylim_slider,self.logcountsdata,self.analytes_)
        else:
            pass
        
    @pn.depends('input_data','ratio_buttons','regression_buttons',
                'background_slider','ablation_slider','ablation_start_true','ratio_plot_ylim_slider',
                'sample_subset',
                'force_exp_params','joggle_params','aguess','bguess','aguess_207','bguess_207')
    def call_ratio_plot(self):
        if self.sample_subset is not None:
            data_toplot = self.input_data[self.input_data['SampleLabel'] == self.sample_subset]
            return plots.ratio_plot(data_toplot,self.ratio_buttons,self.regression_buttons,self.background_slider,
                              self.ablation_slider,self.ablation_start_true,self.ratio_plot_ylim_slider,
                              self.force_exp_params,self.joggle_params,self.aguess,self.bguess,self.aguess_207,self.bguess_207)
    
    @pn.depends('input_data','regression_buttons','ablation_slider','sample_subset','ablation_start_true','force_exp_params','joggle_params','aguess','bguess','aguess_207','bguess_207')
    def call_residuals_plot(self):
        if self.sample_subset is not None:
            data_toplot = self.input_data[self.input_data['SampleLabel'] == self.sample_subset]
            return plots.residuals_plot(data_toplot,self.regression_buttons,self.ablation_slider[0],self.ablation_slider[1],self.ablation_start_true,self.force_exp_params,self.joggle_params,
                                        self.aguess,self.bguess,self.aguess_207,self.bguess_207)
        
    @pn.depends('input_data','power','sample_subset','ablation_slider')
    def call_ellipse_plot(self):
        if self.sample_subset is not None:
            if self.ellipsemode_selector is True:
                data_toplot = self.input_data[self.input_data['SampleLabel'] == self.sample_subset]
                Weth_ell,TW_ell = plots.ellipse_plot(data_toplot, self.power,self.ablation_slider[0],self.ablation_slider[1])
                tabs = pn.Tabs(('TW',TW_ell),('Weth.',Weth_ell),dynamic=True)
                return tabs
            else:
                pass
            
        
    def export_plots(self,event=None):
        data_toplot = self.input_data[self.input_data['SampleLabel'] == self.sample_subset]
        ablationplot = plots.ablation_plot(data_toplot,self.background_slider[0],self.background_slider[1],
                              self.ablation_slider[0],self.ablation_slider[1],self.ablation_start_true,self.ablation_plot_ylim_slider,self.logcountsdata,self.analytes_)
        ratioplot = plots.ratio_plot(data_toplot,self.ratio_buttons,self.regression_buttons,self.background_slider[0],self.background_slider[1],
                          self.ablation_slider[0],self.ablation_slider[1],self.ablation_start_true,
                          self.force_exp_params,self.joggle_params,self.aguess,self.bguess,self.aguess_207,self.bguess_207)
        residualsplot = plots.residuals_plot(data_toplot,self.regression_buttons,self.ablation_slider[0],self.ablation_slider[1],self.ablation_start_true,self.force_exp_params,self.joggle_params,
                                             self.aguess,self.bguess,self.aguess_207,self.bguess_207)
        ablationplot.savefig('Ablation plot.pdf',format='pdf',dpi=250)
        ratioplot.savefig('Ratio plot.pdf',format='pdf',dpi=250)
        residualsplot.savefig('Residuals plot.pdf',format='pdf',dpi=250)
        if self.ellipsemode_selector is True:
            Weth_ell,TW_ell = plots.ellipse_plot(data_toplot, self.power,self.ablation_slider[0],self.ablation_slider[1])
            Weth_ell.savefig('Wetherhill Ellipse.pdf',format='pdf',dpi=250)
            TW_ell.savefig('Tara-Wasserburg Ellipse.pdf',format='pdf',dpi=250)
        else:
            pass
        
    
    @pn.depends('input_data','regression_buttons','ablation_slider','sample_subset','ablation_start_true','force_exp_params','joggle_params')
    def get_regression_stats(self):
        if self.sample_subset is not None:
            data_toanalyze = self.input_data[self.input_data['SampleLabel'] == self.sample_subset]
            data_toanalyze = plots.get_ratios(data_toanalyze)
            data_toanalyze = data_toanalyze[(data_toanalyze.Time_s >= self.ablation_slider[0]) & (data_toanalyze.Time_s <= self.ablation_slider[1])]
            regressions,stats = calc_fncs.get_regressions(data_toanalyze,self.regression_buttons,self.ablation_start_true,self.force_exp_params,self.joggle_params,
                                                          self.aguess,self.bguess,self.aguess_207,self.bguess_207)
            
            divider1 = pn.pane.Markdown('206/238 Reg. Stats:')
            r2_1 = pn.pane.Markdown('$$R^{2}$$ = '+str(round(stats[0],3)))
            r2_2 = pn.pane.Markdown('$$R^{2} 2^{nd} order$$ = '+str(round(stats[1],3)))
            r2_exp = pn.pane.Markdown('$$R^{2}_{exp}$$ = '+str(round(stats[2],3)))
            SE_1stper = pn.pane.Markdown('$$1^{st} Order SE \%$$ = '+str(round(stats[3],3)))
            SE_2ndper = pn.pane.Markdown('$$2^{nd} Order SE \%$$ = '+str(round(stats[4],3)))
            SE_expper = pn.pane.Markdown('$$Exp SE \%$$ = '+str(round(stats[5],3)))
            divider2 = pn.pane.Markdown('207/235 Reg. Stats:')
            SE_1stper_207 = pn.pane.Markdown('$$1^{st} Order SE \%$$ = '+str(round(stats[6],3)))
            SE_2ndper_207 = pn.pane.Markdown('$$2^{nd} Order SE \%$$ = '+str(round(stats[7],3)))
            SE_expper_207 = pn.pane.Markdown('$$Exp SE %$$ = '+str(round(stats[8],3)))
            
            stats_markdown = pn.Row(divider1,r2_1,r2_2,r2_exp,SE_1stper,SE_2ndper,SE_expper,divider2,SE_1stper_207,SE_2ndper_207,SE_expper_207)
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
                                                                   self.ablation_start_true,self.regression_buttons,self.ellipsemode_selector,self.force_exp_params,self.joggle_params,
                                                                   self.aguess,self.bguess,self.aguess_207,self.bguess_207,
                                                                   self.counts_mode,self.arrayofdwelltimes,self.param.sample_subset.objects)
        if self.output_data is None:
            self.output_data = data_approved
        else:
            self.output_data = self.output_data.append(data_approved,ignore_index=True)
        
        if self.ellipsemode_selector is True:
            if self.output_data_ellipse is None:
                self.output_data_ellipse = ellipse_datapproved
            else:
                self.output_data_ellipse = self.output_data_ellipse.append(ellipse_datapproved,ignore_index=True)
        else:
            pass
            
    @pn.depends('output_data','output_data_ellipse')
    def export_data(self,event=None):
        self.output_data = self.output_data.drop(columns=['Time','Time_s'])
        self.output_data.to_excel('output_lasertramZ.xlsx')
        if self.ellipsemode_selector is True and self.output_data_ellipse is not None:
            self.output_data_ellipse.to_excel('output_CEllipse_lasertramZ.xlsx')
        else:
            pass
        

callapp = make_plots(name='Reduce Ablation Data')

pn.extension('tabulator','mathjax')

# buttons_=pn.WidgetBox(pn.Param(callapp.param.accept_array_button,
#                                widgets={'accept_array_button': pn.widgets.Button(name='Accept Detector Array',button_type='success')}))
# buttons_sample=pn.WidgetBox(pn.Param(callapp.param.accept_samplename_button,
#                                widgets={'accept_samplename_button': pn.widgets.Button(name='Accept Sample Name',button_type='success')}))

widgets = {'ratio_buttons': pn.widgets.CheckBoxGroup,
            'regression_buttons': pn.widgets.CheckBoxGroup,
           'counts_mode': pn.widgets.RadioButtonGroup,
           'export_data_button': pn.widgets.Button(name='DDDT!',button_type='success'),
           'analytes_':pn.widgets.CheckBoxGroup,
           'ablation_plot_ylim_slider':pn.widgets.EditableRangeSlider(name='Ablation Plot Y-lim',start=0,end=1e6,value=(0,1e6),step=0.01),
           'ratio_plot_ylim_slider':pn.widgets.EditableRangeSlider(name='Ratio Plot Y-lim',start=0,end=2,value=(0,0.3),step=0.1)
           }

fastgrid_layout = pn.template.VanillaTemplate(title='LaserTRAMZ: LA-ICP-Q-MS',
                                                sidebar=pn.Column(pn.WidgetBox(pn.Param(callapp.param,widgets=widgets))),sidebar_width=380)


fastgrid_layout.main.append(pn.Column(pn.Row(callapp.call_ratio_plot,callapp.call_residuals_plot))) # for vanilla
# fastgrid_layout.main.append(pn.Column(hv.DynamicMap(callapp.call_ratio_plot),callapp.call_residuals_plot)) # for vanilla
fastgrid_layout.main.append(pn.Row(callapp.get_regression_stats))
fastgrid_layout.main.append(pn.Column(callapp.call_ablation_plot)) # for vanilla
fastgrid_layout.main.append(pn.Column(pn.Row(callapp.call_ellipse_plot))) # for vanilla
fastgrid_layout.main.append(pn.Column(callapp._update_output_widget))
fastgrid_layout.show();









    
