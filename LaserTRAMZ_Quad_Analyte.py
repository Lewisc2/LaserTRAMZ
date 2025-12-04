#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Chuck Lewis
"""

# %% Import Dependables
import pandas as pd
import numpy as np
import bokeh
from bokeh.plotting import figure
from bokeh.layouts import column, row, gridplot
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
import matplotlib as mpl
mpl.use('agg')
from matplotlib.figure import Figure
from matplotlib.patches import Ellipse
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

color_palette = bokeh.palettes.Muted9
color_palette_regressions = bokeh.palettes.Light3
hv.extension('bokeh')

# %% Calculation Functions
class calc_fncs:
    """ Class that holds all of the functions for reducing the time resolved data"""
    def __init__(self,*args):
        for a in args:
            self.__setattr__(str(a), args[0])
    
    
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
            outlier_removed_ablation = calc_fncs.threesigoutlierremoval(ablation,analyte,'Ablation')
            backgrounds.loc[:,analyte] = backgrounds
            ablation.loc[:,analyte] = outlier_removed_ablation
            
        meanbackgrounds = backgrounds.loc[:,'202Hg':'238U'].mean()
        lods = 3*backgrounds.loc[:,'202Hg':'238U'].std()
        ablation_backsub = ablation.loc[:,'202Hg':'238U'].sub(meanbackgrounds,axis='columns').clip(lower=0)
        ablation_backsub.insert(0,'Time_s',data['Time_s'])
        ablation_backsub = ablation_backsub.reset_index(drop=True)
        
        for analyte in ablation_backsub.loc[:,'202Hg':'238U'].columns:
            if ablation_backsub[analyte].mean() < lods[analyte] and analyte in ['238U', '206Pb']:
                pn.state.notifications.error('Warning! '+str(analyte)+' is b.d.l. - Can Not Reduce',duration=3000)
            else:
                pass
            
        return ablation_backsub,meanbackgrounds,lods
    
    
    def get_tresolved_ratios(data):
        """
        function that returns time resolved ratios for the entire analysis from background start to ablation end
        used for visualization and thre sig outlier removal
        
        Parameters
        ----------
        data : pandas dataframe
            dataframe of time resolved intensities
        bstart : float
            user selected value for the background interval start
        tend : float
            user selected value for the ablation interval end
    
        Returns
        -------
        tresolvedr : pandas dataframe
            dataframe with time resovled ratios
    
        """
        r206238 = np.divide(data['206Pb'].astype(float),data['238U'].astype(float),out=np.zeros_like(data['238U'].astype(float)),where=data['238U'].astype(float)>0) if np.mean(data['238U']>0) else np.zeros_like(data.iloc[:,1])
        r238206 = np.divide(data['238U'].astype(float),data['206Pb'].astype(float),out=np.zeros_like(data['206Pb'].astype(float)),where=data['206Pb'].astype(float)>0) if np.mean(data['206Pb']>0) else np.zeros_like(data.iloc[:,1])
        r207235 = np.divide(data['207Pb'].astype(float),data['235U'].astype(float),out=np.zeros_like(data['235U'].astype(float)),where=data['235U'].astype(float)>0) if np.mean(data['235U']>0) else np.divide(data['207Pb'].astype(float),(data['238U']/137.818).astype(float),out=np.zeros_like(data['238U'].astype(float)),where=data['238U'].astype(float)!=0) if np.mean(data['238U']>0) else np.zeros_like(data.iloc[:,1])
        r208232 = np.divide(data['208Pb'].astype(float),data['232Th'].astype(float),out=np.zeros_like(data['232Th'].astype(float)),where=data['232Th'].astype(float)>0) if np.mean(data['232Th']>0) else np.zeros_like(data.iloc[:,1])
        r207206 = np.divide(data['207Pb'].astype(float),data['206Pb'].astype(float),out=np.zeros_like(data['206Pb'].astype(float)),where=data['206Pb'].astype(float)>0) if np.mean(data['206Pb']>0) else np.zeros_like(data.iloc[:,1])
        r207204 = np.divide(data['207Pb'].astype(float),data['204Pb'].astype(float),out=np.zeros_like(data['204Pb'].astype(float)),where=data['204Pb'].astype(float)>0) if np.mean(data['204Pb']>0) else np.zeros_like(data.iloc[:,1])
        r206204 = np.divide(data['206Pb'].astype(float),data['204Pb'].astype(float),out=np.zeros_like(data['204Pb'].astype(float)),where=data['204Pb'].astype(float)>0) if np.mean(data['204Pb']>0) else np.zeros_like(data.iloc[:,1])
        r238232 = np.divide(data['238U'].astype(float),data['232Th'].astype(float),out=np.zeros_like(data['232Th'].astype(float)),where=data['232Th'].astype(float)>0) if np.mean(data['232Th']>0) else np.zeros_like(data.iloc[:,1])
        r238235 = np.divide(data['238U'].astype(float),data['235U'].astype(float),out=np.zeros_like(data['235U'].astype(float)),where=data['235U'].astype(float)>0) if np.mean(data['235U']>0) else np.full_like(data.iloc[:,1], 137.818)
     
        tresolved_ratio_list = [r206238,r238206,r207235,r208232,r207206,r207204,r206204,r238232,r238235]
        ratiolist = ['206Pb/238U','238U/206Pb','207Pb/235U','208Pb/232Th','207Pb/206Pb','207Pb/204Pb','206Pb/204Pb','238U/232Th','238U/235U']
        stacked_tresolvedr = np.stack(tresolved_ratio_list,axis=-1)
        tresolvedr = pd.DataFrame(stacked_tresolvedr,columns=ratiolist)
        tresolvedata = pd.concat([data.reset_index(drop=True),tresolvedr.reset_index(drop=True)],axis=1).reset_index(drop=True)
        
        return tresolvedata
        
        
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
        if len(data[variable]) <= 0 and intervaltype == 'Ablation':
            pn.state.notifications.warning('Sporadic '+str(data[variable].name)+' Signal - Outlier Removal Failed',duration=5000)
                   
        return data[variable]
    
    
    
    def get_mean_ratios(data,ratiotype):
        """
        Function used to get reduced isotope ratios of data. Values are done with either ratio of means method or geometric mean
        Uncertainties are calculated according to selection  - standard error of time resolved ratio or geometric standard error of time resolved ratio

        Parameters
        ----------
        data : pandas dataframe
            dataframe containing values to reduce. Must be counts data as output from the function backgroundsubtract_convert()
        ratiotype : string
            string value specifying which method will be used to reduce ratios.
            Must be either 'Ratio of Means' or 'Geometric'

        Returns
        -------
        reduced_ratios : pandas dataframe
            dataframe with reduced ratios and uncertainties in 1SE%

        """
        if ratiotype == 'Ratio of Means':
            mu_206Pb238U = np.mean(data['206Pb'])/np.mean(data['238U']) if np.mean(data['238U'])>0 else 0
            se_206Pb238U = 0 if mu_206Pb238U==0 else np.std(np.divide(data['206Pb'].astype(float),data['238U'].astype(float),out=np.zeros_like(data['238U'].astype(float)),where=data['238U'].astype(float)!=0),ddof=1)/np.sqrt(len(data))
            mu_238U206Pb = np.mean(data['238U'])/np.mean(data['206Pb']) if np.mean(data['206Pb'])>0 else 0
            se_238U206Pb = 0 if mu_238U206Pb==0 else np.std(np.divide(data['238U'].astype(float),data['206Pb'].astype(float),out=np.zeros_like(data['238U'].astype(float)),where=data['206Pb'].astype(float)!=0),ddof=1)/np.sqrt(len(data))
            try:
                mu_207Pb235U = np.mean(data['207Pb'])/np.mean(data['235U']) if np.mean(data['235U'])>0 else np.mean(data['207Pb'])/np.mean(data['238U']/137.818)
                se_207Pb235U = np.std(np.divide(data['207Pb'].astype(float),data['238U'].to_numpy()/137.818,out=np.zeros_like(data['238U'].astype(float)),where=data['238U'].astype(float)!=0),ddof=1)/np.sqrt(len(data)) if np.mean(data['235U'])<=0 else np.std(np.divide(data['207Pb'].astype(float),data['235U'].astype(float),out=np.zeros_like(data['235U'].astype(float)),where=data['235U'].astype(float)!=0),ddof=1)/np.sqrt(len(data))
            except:
                mu_207Pb235U = 0
                se_207Pb235U = 0
            mu_208Pb232Th = np.mean(data['208Pb'])/np.mean(data['232Th']) if np.mean(data['232Th'])>0 else 0
            se_208Pb232Th = 0 if mu_208Pb232Th==0 else np.std(np.divide(data['208Pb'].astype(float),data['232Th'].astype(float),out=np.zeros_like(data['232Th'].astype(float)),where=data['232Th'].astype(float)!=0),ddof=1)/np.sqrt(len(data))
            mu_207Pb206Pb = np.mean(data['207Pb'])/np.mean(data['206Pb']) if np.mean(data['206Pb'])>0 else 0
            se_207Pb206Pb = 0 if mu_207Pb206Pb==0 else np.std(np.divide(data['207Pb'].astype(float),data['206Pb'].astype(float),out=np.zeros_like(data['206Pb'].astype(float)),where=data['206Pb'].astype(float)!=0),ddof=1)/np.sqrt(len(data))
            mu_207Pb204Pb = np.mean(data['207Pb'])/np.mean(data['204Pb']) if np.mean(data['204Pb'])>0 else 0
            se_207Pb204Pb = 0 if mu_207Pb204Pb==0 else np.std(np.divide(data['207Pb'].astype(float),data['204Pb'].astype(float),out=np.zeros_like(data['204Pb'].astype(float)),where=data['204Pb'].astype(float)!=0))/np.sqrt(len(data))
            mu_206Pb204Pb = np.mean(data['206Pb'])/np.mean(data['204Pb']) if np.mean(data['204Pb'])>0 else 0
            se_206Pb204Pb = 0 if mu_206Pb204Pb==0 else np.std(np.divide(data['206Pb'].astype(float),data['204Pb'].astype(float),out=np.zeros_like(data['204Pb'].astype(float)),where=data['204Pb'].astype(float)!=0),ddof=1)/np.sqrt(len(data))
            mu_238U232Th = np.mean(data['238U'])/np.mean(data['232Th']) if np.mean(data['232Th'])>0 else 0
            se_238U232Th = 0 if mu_238U232Th==0 else np.std(np.divide(data['238U'].astype(float),data['232Th'].astype(float),out=np.zeros_like(data['232Th'].astype(float)),where=data['232Th'].astype(float)!=0),ddof=1)/np.sqrt(len(data))
            mu_238U235U = np.mean(data['238U'])/np.mean(data['235U']) if np.mean(data['235U'])>0 else 137.818
            se_238U235U = 137.818 if mu_238U235U==0 else np.std(np.divide(data['238U'].astype(float),data['235U'].astype(float),out=np.zeros_like(data['235U'].astype(float)),where=data['235U'].astype(float)!=0),ddof=1)/np.sqrt(len(data))
            
        elif ratiotype == 'Geometric':
            mu_206Pb238U = stats.gmean(np.divide(data['206Pb'].astype(float),data['238U'].astype(float),out=np.ones_like(data['238U'].astype(float)),where=(data['238U'].astype(float)>0) & (data['206Pb'].astype(float)>0))) if np.mean(data['238U'])>0 else 0
            se_206Pb238U = 0 if mu_206Pb238U==0 else stats.gstd(np.divide(data['206Pb'].astype(float),data['238U'].astype(float),out=np.ones_like(data['238U'].astype(float)),where=(data['238U'].astype(float)>0) & (data['206Pb'].astype(float)>0)))/np.sqrt(len(data))
            mu_238U206Pb = 1/mu_206Pb238U
            se_238U206Pb = se_206Pb238U
            try:
                mu_207Pb235U = stats.gmean(np.divide(data['207Pb'].astype(float),data['235U'].astype(float),out=np.ones_like(data['235U'].astype(float)),where=(data['235U'].astype(float)>0) & (data['207Pb'].astype(float)>0))) if np.mean(data['235U'].astype(float))>0 else stats.gmean(np.divide(data['207Pb'].astype(float),data['238U'].to_numpy()/137.818,out=np.ones_like(data['238U'].astype(float)),where=(data['238U'].astype(float)>0)) & (data['207Pb'].astype(float)>0))
                se_207Pb235U = stats.gstd(np.divide(data['207Pb'].astype(float),data['238U'].to_numpy()/137.818,out=np.ones_like(data['238U'].astype(float)),where=(data['235U'].astype(float)>0) & (data['207Pb'].astype(float)>0))) if np.mean(data['235U'].astype(float))<=0 else stats.gstd(np.divide(data['207Pb'].astype(float),data['235U'].astype(float),out=np.ones_like(data['235U'].astype(float)),where=(data['235U'].astype(float)>0)) & (data['207Pb'].astype(float)>0))/np.sqrt(len(data))
            except:
                mu_207Pb235U = 0
                se_207Pb235U = 0
            mu_208Pb232Th = stats.gmean(np.divide(data['208Pb'].astype(float),data['232Th'].astype(float),out=np.ones_like(data['232Th'].astype(float)),where=(data['232Th'].astype(float)>0) & (data['208Pb'].astype(float)>0))) if np.mean(data['232Th'])>0 else 0
            se_208Pb232Th = 0 if mu_208Pb232Th==0 else stats.gstd(np.divide(data['208Pb'].astype(float),data['232Th'].astype(float),out=np.ones_like(data['232Th'].astype(float)),where=(data['232Th'].astype(float)>0) & (data['208Pb'].astype(float)>0)))/np.sqrt(len(data))
            mu_207Pb206Pb = stats.gmean(np.divide(data['207Pb'].astype(float),data['206Pb'].astype(float),out=np.ones_like(data['206Pb'].astype(float)),where=(data['207Pb'].astype(float)>0) & (data['206Pb'].astype(float)>0))) if np.mean(data['206Pb'])>0 else 0
            se_207Pb206Pb = 0 if mu_207Pb206Pb==0 else stats.gstd(np.divide(data['207Pb'].astype(float),data['206Pb'].astype(float),out=np.ones_like(data['206Pb'].astype(float)),where=(data['207Pb'].astype(float)>0) & (data['206Pb'].astype(float)>0)))/np.sqrt(len(data))
            mu_207Pb204Pb = stats.gmean(np.divide(data['207Pb'].astype(float),data['204Pb'].astype(float),out=np.ones_like(data['204Pb'].astype(float)),where=(data['207Pb'].astype(float)>0) & (data['204Pb'].astype(float)>0))) if np.mean(data['204Pb'])>0 else 0
            se_207Pb204Pb = 0 if mu_207Pb204Pb==0 else stats.gstd(np.divide(data['207Pb'].astype(float),data['204Pb'].astype(float),out=np.ones_like(data['204Pb'].astype(float)),where=(data['207Pb'].astype(float)>0) & (data['204Pb'].astype(float)>0)))/np.sqrt(len(data))
            mu_206Pb204Pb = stats.gmean(np.divide(data['206Pb'].astype(float),data['204Pb'].astype(float),out=np.ones_like(data['204Pb'].astype(float)),where=(data['206Pb'].astype(float)>0) & (data['204Pb'].astype(float)>0))) if np.mean(data['204Pb'])>0 else 0
            se_206Pb204Pb = 0 if mu_206Pb204Pb==0 else stats.gstd(np.divide(data['206Pb'].astype(float),data['204Pb'].astype(float),out=np.ones_like(data['204Pb'].astype(float)),where=(data['206Pb'].astype(float)>0) & (data['204Pb'].astype(float)>0)))/np.sqrt(len(data))
            mu_238U232Th = stats.gmean(np.divide(data['238U'].astype(float),data['232Th'].astype(float),out=np.ones_like(data['232Th'].astype(float)),where=(data['238U'].astype(float)>0) & (data['232Th'].astype(float)>0))) if np.mean(data['232Th'])>0 else 0
            se_238U232Th = 0 if mu_238U232Th==0 else stats.gstd(np.divide(data['238U'].astype(float),data['232Th'].astype(float),out=np.ones_like(data['232Th'].astype(float)),where=(data['238U'].astype(float)>0) & (data['232Th'].astype(float)>0)))/np.sqrt(len(data))
            mu_238U235U = stats.gmean(np.divide(data['238U'].astype(float),data['235U'].astype(float),out=np.ones_like(data['235U'].astype(float)),where=(data['238U'].astype(float)>0) & (data['235U'].astype(float)>0))) if np.mean(data['235U'])>0 else 137.818
            se_238U235U = 137.818 if mu_238U235U==0 else stats.gstd(np.divide(data['238U'].astype(float),data['235U'].astype(float),out=np.ones_like(data['235U'].astype(float)),where=(data['238U'].astype(float)>0) & (data['235U'].astype(float)>0)))/np.sqrt(len(data))
            
        means_array = np.array([mu_206Pb238U,mu_238U206Pb,mu_207Pb235U,mu_208Pb232Th,mu_207Pb206Pb,mu_207Pb204Pb,mu_206Pb204Pb,mu_238U232Th,mu_238U235U])
        uncertainties_array = np.array([se_206Pb238U,se_238U206Pb,se_207Pb235U,se_208Pb232Th,se_207Pb206Pb,se_207Pb204Pb,se_206Pb204Pb,se_238U232Th,se_238U235U])
        uncertainties_array = uncertainties_array/means_array*100
        full_array = np.concatenate((means_array,uncertainties_array))
        ratio_uncertainties_list = ['206Pb/238U','238U/206Pb','207Pb/235U','208Pb/232Th','207Pb/206Pb','207Pb/204Pb','206Pb/204Pb','238U/232Th','238U/235U',
                                    'SE% 206Pb/238U','SE% 238U/206Pb','SE% 207Pb/235U','SE% 208Pb/232Th','SE% 207Pb/206Pb','SE% 207Pb/204Pb','SE% 206Pb/204Pb','SE% 238U/232Th','SE% 238U/235U'
                                    ]
        reduced_ratios = pd.DataFrame([full_array],columns=ratio_uncertainties_list)
        
        return reduced_ratios
    

    # Note: Input data must be counts data only - no ratios output from the function used to get ratios
    def get_regressions(data,regression_buttons,ablation_start_true):
        """
        Function used to get regressions of time-resolved Pb/U data. 
        Returns either 1st order regression or exponential regresssion
        Stats and regression parameters returned depends on calling function

        Parameters
        ----------
        data : pandas dataframe
            pandas dataframe of coutns data. Note that this cannot be data that already has ratios in the dataframe as these are calculated here
        regression_buttons : list
            list of arguments specifying if the user wants 1st order or exp. regression. Currently get 1st or both - need to implement way to get only one
        ablation_start_true : float
            float value indicating where to set time zero intercept in ablation. by default set to tstart in program

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        data = data.reset_index(drop=True)
        t0 = data.loc[0,'Time_s']
        data['Time_s'] = data['Time_s']-t0
        ablation_start_true = ablation_start_true-t0
        y = data['206Pb/238U'].to_numpy(dtype=float)
        y207 = data['207Pb/235U'].to_numpy(dtype=float)
        x = data['Time_s'].to_numpy(dtype=float)
        X = sm.add_constant(x)
        if '1st Order' in regression_buttons:
            # 206Pb/238U 1st order regression
            linmod1 = sm.OLS(y, X).fit() # fit a linear model on the data
            predicted1 = linmod1.predict() # predicted values from y=mx+b
            rsquared1 = linmod1.rsquared # get the R2 value of the regression
            predicted_b01 = linmod1.params[0] + linmod1.params[1]*ablation_start_true # get the predicted value at the ablation start that is input by the user - ntoe this allows projection if desired - otherwise could use fit.params[0] or fit.params['const']
            sigma1 = np.sqrt(linmod1.ssr/linmod1.df_resid) # get the 1SD (Sum Squared residuals / residual degrees of freedom)^(1/2)
            SE_b01 = sigma1*np.sqrt(1/len(data)+(data['Time_s'].mean())**2/((len(data)-1)*np.var(data['Time_s'],ddof=2))) # get standard error for a single point estiamted by a regression model
            SE_b01_percent = SE_b01/predicted_b01*100 # get the % 1SE
            resid1 = linmod1.resid # get the residuals of the regression
            # 207Pb/235U 1st order regression - equations and functions as above
            linmod1_207 = sm.OLS(y207,X).fit()
            predicted1_207 = linmod1_207.predict()
            predicted_b01_207 = linmod1_207.params[0] + linmod1_207.params[1]*ablation_start_true
            sigma1_207 = np.sqrt(linmod1_207.ssr/linmod1_207.df_resid) 
            SE_b01_207 = sigma1_207*np.sqrt(1/len(data)+(data['Time_s'].mean())**2/((len(data)-1)*np.var(data['Time_s'],ddof=2)))
            SE_b01_percent_207 = SE_b01_207/predicted_b01_207*100 
            resid1_207 = linmod1_207.resid
            
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
            
        if 'Exp. Regression' in regression_buttons:
            # define the exponential function adn a simplified exponential function if iterations fail
            def exp_func(x,a,b,c):
                return a*np.exp(-b*x)+c
            def simple_exp_func(x,a,b):
                return a*np.exp(-b*x)
            # define a string variable to recognize if runtime error was excepted or not
            curve638 = 'Three Variables'
            curve735 = 'Three Variables'
            
            default_206238_initparams = [0.1,0.05,0.02]
            default_207235_initparams = [1,0.05,0.02]
            
            try:
                popt,pcov = curve_fit(exp_func,data['Time_s'].to_numpy(dtype=float),data['206Pb/238U'].to_numpy(dtype=float),p0=default_206238_initparams) # fit data to an exponential curve
            except RuntimeError:
                print('Runtime Error: simplifying exponential function 206/238')
                newparams_638 = [0.1,0.05]
                popt,pcov = curve_fit(simple_exp_func,data['Time_s'].to_numpy(dtype=float),data['206Pb/238U'].to_numpy(dtype=float),p0=newparams_638) # fit data to two parameter exponential curve in case runtime error occurred
                curve638 = 'Two Variables'
            except Exception as e:
                pn.state.notifications.error('206Pb/238U Exp. regression failing due to: '+str(e),duration=10000)
                
            try:
                popt_207,pcov_207 = curve_fit(exp_func,data['Time_s'].to_numpy(dtype=float),data['207Pb/235U'].to_numpy(dtype=float),p0=default_207235_initparams)
            except RuntimeError:
                print('Runtime Error: simplifying exponential function 207/235')
                newparams_735 = [0.9,0.05]
                popt_207,pcov_207 = curve_fit(simple_exp_func,data['Time_s'].to_numpy(dtype=float),data['207Pb/235U'].to_numpy(dtype=float),p0=newparams_735)
                curve735 = 'Two Variables'
            except Exception as e:
                pn.state.notifications.error('207Pb/235U Exp. regression failing due to: '+str(e),duration=10000)
                
            if curve638 == 'Three Variables':
                predictedexp = exp_func(data['Time_s'].to_numpy(dtype=float),*popt) # get predicted values for exponential curve
                predicted_b0exp = popt[0]*np.exp(-popt[1]*ablation_start_true)+popt[2] # get zero-intercept for exponential curve
            elif curve638 == 'Two Variables':
                predictedexp = simple_exp_func(data['Time_s'].to_numpy(dtype=float),*popt)
                predicted_b0exp = popt[0]*np.exp(-popt[1]*ablation_start_true)
            if curve735 == 'Three Variables':
                predictedexp_207 = exp_func(data['Time_s'].to_numpy(dtype=float),*popt_207)
                predicted_b0exp_207 = popt_207[0]*np.exp(-popt_207[1]*ablation_start_true)+popt_207[2]
            elif curve735 == 'Two Variables':
                predictedexp_207 = simple_exp_func(data['Time_s'].to_numpy(dtype=float),*popt_207)
                predicted_b0exp_207 = popt_207[0]*np.exp(-popt_207[1]*ablation_start_true)
            
            
            # initialize arrays to be filled with residuals and squared residuals
            resid = np.zeros(len(predictedexp))
            sq_resid = np.zeros(len(predictedexp))
            resid_207 = np.zeros(len(predictedexp))
            sq_resid_207 = np.zeros(len(predictedexp))
            # for i,k,m,l in zip(range(0,len(predictedexp)),range(0,len(data['206Pb/238U'])),range(0,len(predictedexp_207)),range(0,len(data['207Pb/235U']))):
            for i in range(0,len(predictedexp)):
                resid[i] = ((data['206Pb/238U'][i] - predictedexp[i])) # get residuals
                sq_resid[i] = (((data['206Pb/238U'][i] - predictedexp[i])**2)) # square them
                resid_207[i] = ((data['207Pb/235U'][i]-predictedexp_207[i]))
                sq_resid_207[i] = ((data['207Pb/235U'][i]-predictedexp_207[i])**2)
            sum_sq_resid = np.sum(sq_resid) # sum the squared residuals
            sum_sq_resid_207 = np.sum(sq_resid_207)
            sigmaexp = np.sqrt(sum_sq_resid/(len(data['206Pb/238U'])-2)) # denominator = d.f. = n-#params
            sigmaexp_207 = np.sqrt(sum_sq_resid_207/(len(data['207Pb/235U'])-2)) # denominator = d.f. = n-#params
            SE_b0exp = sigmaexp*np.sqrt(1/len(data['206Pb/238U'])+(np.mean(data['Time_s']))**2/((len(data)-1)*np.var(data['Time_s'],ddof=2)))
            SE_b0exp_207 = sigmaexp_207*np.sqrt(1/len(data['207Pb/235U'])+(np.mean(data['Time_s']))**2/((len(data)-1)*np.var(data['Time_s'],ddof=2)))
            SE_b0exp_percent = SE_b0exp/predicted_b0exp*100
            SE_b0exp_percent_207 = SE_b0exp_207/predicted_b0exp_207*100
            residexp = resid # reassign variables
            residexp_207 = resid_207
            tss = ((data['206Pb/238U'] - np.mean(data['206Pb/238U']))**2).sum() # get the total sum of squared residuals
            rsquared_exp = 1 - (sum_sq_resid/tss) # get the r-squared for the exponential regression
            rsquaredexp_adj = 1 - ( (1-rsquared_exp)*(len(data['206Pb/238U'])-1) / (len(data['206Pb/238U'])-2-1) ) # get the adjusted r squared for the exponential regression
        else:
            # set everything to zeros if exponenetial regression not chosen
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
            ratios_to_return = [predicted_b01,predicted_b0exp,predicted_b01_207,predicted_b0exp_207]
            stats_to_return = [SE_b01_percent,SE_b0exp_percent,
                               SE_b01_percent_207,SE_b0exp_percent_207]
            
            return ratios_to_return,stats_to_return
        
        elif callingmethod == 'ratio_plot':
            regressions_to_return = [predicted1,predictedexp,predicted1_207,predictedexp_207]
            stats_to_report = [rsquared1,rsquaredexp_adj,SE_b01_percent,SE_b0exp_percent]
            
            return regressions_to_return,stats_to_report
        
        elif callingmethod == 'evaluate_output_data':
            # ""
            regressions_to_return = [predicted1,predictedexp]
            stats_to_report = [rsquared1,rsquared_exp,SE_b01_percent,SE_b0exp_percent,SE_b01_percent_207,SE_b0exp_percent_207]
            
            return regressions_to_return,stats_to_report
        
        elif callingmethod == 'residuals_plot':
            predicted_to_return = [predicted1,predictedexp]
            predicted_to_return_207 = [predicted1_207,predictedexp_207]
            resid_to_return = [resid1,residexp]
            resid_to_return_207 = [resid1_207,residexp_207]
            
            return predicted_to_return,predicted_to_return_207,resid_to_return,resid_to_return_207
        
        elif callingmethod == 'get_ellipse':
            predicted1_ellipse = resid1 + linmod1.params[0]
            predicted1_207_ellipse = resid1_207 + linmod1_207.params[0]
            if 'Exp. Regression' in regression_buttons:
                if curve638 == 'Three Variables':
                    predictedexp_ellipse = residexp + popt[0] + popt[2]
                else:
                    predictedexp_ellipse = residexp + popt[0]
                if curve735 == 'Three Variables':
                    predictedexp_207_ellipse = residexp_207 + popt_207[0] + popt_207[2]
                else:
                    predictedexp_207_ellipse = residexp_207 + popt_207[0]
            else:
                predictedexp_ellipse = np.zeros_like(data['Time_s'])
                predictedexp_207_ellipse = np.zeros_like(data['Time_s'])

            predicted_to_return = [predicted1_ellipse,predictedexp_ellipse]
            predicted_to_return_207 = [predicted1_207_ellipse,predictedexp_207_ellipse]
            
            return predicted_to_return,predicted_to_return_207
        
        else:
            pass
        
        
        
    def get_approved(data,bstart,bend,tstart,tend,tproject,regression_buttons,
                    ratio_type,counts_mode,arrayofdwelltimes,sample_names,power):
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
        counts_data,backgrounds,lods = calc_fncs.backgroundsubtract_convert_lod(data_toapprove,bstart,bend,tstart,tend,arrayofdwelltimes) # background subtracted counts, mean background intensities
        counts_data = counts_data.fillna(0) # fill any na values with zero in case 'ends' of data has nans after masking for outlier removal
        mean_ratios = calc_fncs.get_mean_ratios(counts_data,ratio_type) # mean ratios - either ratio of means or geometric depending on user input
        t_resolved_ratios = calc_fncs.get_tresolved_ratios(counts_data) # get time resovled ratios of the ablation interval
        
        if counts_mode == 'Total Counts':
            mean_ratios = mean_ratios
            
        elif counts_mode == 'Means & LIEF':
            data_ratios_reg,data_stats_reg = calc_fncs.get_regressions(t_resolved_ratios,regression_buttons,tproject) # get estimated regression intercepts, background subtracted ratios, regression statistics, and SE's of the ratios3
            if ('1st Order' in regression_buttons) and ('Exp. Regression' not in regression_buttons):
                mean_ratios['206Pb/238U'] = data_ratios_reg[0]
                mean_ratios['238U/206Pb'] = 1/data_ratios_reg[0]
                mean_ratios['207Pb/235U'] = data_ratios_reg[2]
                mean_ratios['SE% 206Pb/238U'] = data_stats_reg[0]
                mean_ratios['SE% 238U/206Pb'] = data_stats_reg[0]
                mean_ratios['SE% 207Pb/235U'] = data_stats_reg[2]
            else:
                mean_ratios['206Pb/238U'] = data_ratios_reg[1]
                mean_ratios['238U/206Pb'] = 1/data_ratios_reg[1]
                mean_ratios['207Pb/235U'] = data_ratios_reg[3]
                mean_ratios['SE% 206Pb/238U'] = data_stats_reg[1]
                mean_ratios['SE% 238U/206Pb'] = data_stats_reg[1]
                mean_ratios['SE% 207Pb/235U'] = data_stats_reg[3]


        isotopes_SE = pd.DataFrame([counts_data.loc[:,'202Hg':'238U'].sem()]).add_suffix('_1SE') # get SE's of the individual isotopes
        isotopes_means = pd.DataFrame([counts_data.loc[:,'202Hg':'238U'].mean()]) # put means of isotopic values in a dataframe
    
        Weth_ellparams,TW_ellparams,x1,y1,y2 = calc_fncs.get_ellipse(t_resolved_ratios,power,tproject,regression_buttons) # get confidence ellipse paramters
        # turn lists into dataframe to get joined into one large dataframe that gets sent to the output data
        Weth_ellparams = pd.DataFrame([Weth_ellparams],columns=['Weth C','Weth Wid1','Weth Wid2','Weth rho'])
        TW_ellparams = pd.DataFrame([TW_ellparams],columns=['TW C','TW Wid1','TW Wid2','TW rho'])
        
        # take the means of the individual isotopes and run them through a for loop with the detection limits. Assign a value of 'bdl' if below detection limit. Otherwise, leave unchanged
        for analyte in isotopes_means.loc[:,'202Hg':'238U'].columns:
            if isotopes_means.loc[0,analyte]<=lods[analyte]:
                isotopes_means.loc[0,analyte]='bdl'
                pn.state.notifications.warning(str(analyte)+' is b.d.l.',duration=2000)
            elif isotopes_means.loc[0,analyte]>lods[analyte]:
                isotopes_means.loc[0,analyte]=isotopes_means.loc[0,analyte]
             
        # 202Hg: 29.86%, 204Hg: 6.87% https://www.nndc.bnl.gov/nudat3/
        Hgratio = (6.87/100) / (29.86/100)
        # subtract the isobaric interference of 204Hg on 204Pb using the measured 202Hg, so long as 202 > 'bdl'
        if isotopes_means.loc[0,'202Hg'] != 'bdl':
            # 202Hg: 29.86%, 204Hg: 6.87% https://www.nndc.bnl.gov/nudat3/
            isotopes_means['204Hg'] = isotopes_means['202Hg']*Hgratio # calculate 204Hg from measured 202Hg based on isotopic abundance
            # subtract the 204Hg from the 204 signal, so long as the 204 signal > 'bdl'
            if isotopes_means['204Pb'] != 'bdl':
                isotopes_means['204Pb'] = isotopes_means['204Pb'] - isotopes_means['204Hg']
                # Recheck to make sure newly calculated if the newly calculated 204 signal is greater or less than 'bdl'. Assign 'bdl' or leave unchanged appropriately.
                if isotopes_means['204Pb'] <= lods['204Pb']:
                    isotopes_means['204Pb'] = 'bdl'
                elif isotopes_means['204Pb'] > lods['204Pb']:
                    isotopes_means['204Pb'] = isotopes_means['204Pb']
        else:
            isotopes_means['204Hg'] = 'bdl'
        
        isotopes_means.insert(0,'measurementindex',sample_names.index(data.iloc[0,0]))
        isotopes_means.insert(1,'SampleLabel',data.iloc[0,0]) # reinsert the sample label into the calculation df
        isotopes_means.insert(2,'t start',[tstart]) # insert ablation start time into df
        isotopes_means.insert(3,'t end',[tend]) # insert ablation stop time into df
        isotopes_means.insert(4,'t project',[tproject]) # insert projected regression start time into df
        isotopes_means.insert(5,'b start',[bstart])
        isotopes_means.insert(6,'b end',[bend])
        # stitch the individual isotopic ratios, their errors, and ellipsoid information into a dataframe
        # these are then appeneded into the output df
        data_approved = isotopes_means.join([isotopes_SE,mean_ratios,Weth_ellparams,TW_ellparams])
        
        
        ratio_cols = ['206Pb/238U','238U/206Pb','207Pb/235U','208Pb/232Th','207Pb/206Pb','207Pb/204Pb','206Pb/204Pb','238U/232Th','238U/235U']
        for i in isotopes_means.columns:
            for k in ratio_cols:
                if k.__contains__(str(i)) and data_approved[i].item() == 'bdl':
                    data_approved[k] = 'bdl'
                else:
                    pass
                    
        
        return data_approved
    
    def get_ellipse(data,power,ablation_start_true,regression_buttons):
        data = data.fillna(0)
        drop_condn = data[(data['206Pb/238U'] == 0) | (data['207Pb/235U'] == 0) | (data['207Pb/206Pb'] == 0)].index
        data.drop(drop_condn,inplace=True)
        data = data.reset_index(drop=True)
        
        try:
            predicted,predicted_207 = calc_fncs.get_regressions(data,regression_buttons,ablation_start_true)
        except KeyError:
            meansdict = {'206Pb/238U':np.mean(data['206Pb/238U']),'207Pb/235U':np.mean(data['207Pb/235U']),'207Pb/206Pb':np.mean(data['207Pb/206Pb'])}
            for key,value in meansdict.items():
                if value is np.nan or value <= 0:
                    pn.state.notifications.warning(str(key)+' has an isotope b.d.l. - assigning standard normal distribution to pass ellipse calculations',duration=10000)
                    data[key] = np.random.standard_normal(size=len(data))
                    predicted,predicted_207 = calc_fncs.get_regressions(data,regression_buttons,ablation_start_true)
                    
            
        if ('1st Order' in regression_buttons) and ('Exp. Regression' not in regression_buttons):
            data['207Pb/235U'] = predicted_207[0]
            data['206Pb/238U'] = predicted[0]
        else:
            data['207Pb/235U'] = predicted_207[1]
            data['206Pb/238U'] = predicted[1]
        
        x1 = data['207Pb/235U']
        y1 = data['206Pb/238U']
        x2 = 1/data['206Pb/238U']
        y2 = data['207Pb/206Pb']
        
        cov1 = np.cov(x1,y1)
        cov2 = np.cov(x2,y2)
        if np.isnan(np.min(cov1)):
            cov1 = np.nan_to_num(cov1)
        if np.isnan(np.min(cov2)):
            cov2 = np.nan_to_num(cov2)
        
        
        eigval1,eigvec1 = np.linalg.eig(cov1)
        eigval2,eigvec2 = np.linalg.eig(cov2)
        order1 = eigval1.argsort()[::-1]
        order2 = eigval2.argsort()[::-1]
        eigvals_order1 = eigval1[order1]
        eigvals_order2 = eigval2[order2]
        eigvecs_order1 = eigvec1[:,order1]
        eigvecs_order2 = eigvec2[:,order2]
        
        c1x = np.mean(x1)
        c1y = np.mean(y1)
        c2x = np.mean(x2)
        c2y = np.mean(y2)
        c1 = (c1x,c1y)
        c2 = (c2x,c2y)
        wid1 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power),df=2)*eigvals_order1[0])
        hgt1 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power),df=2)*eigvals_order1[1])
        wid2 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power),df=2)*eigvals_order2[0])
        hgt2 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power),df=2)*eigvals_order2[1])
        theta1 = np.degrees(np.arctan2(*eigvecs_order1[:,0][::-1]))
        theta2 = np.degrees(np.arctan2(*eigvecs_order2[:,0][::-1]))
        
        ell1_params = [c1,wid1,hgt1,theta1]
        ell2_params = [c2,wid2,hgt2,theta2]
        
        return ell1_params,ell2_params,x1,y1,y2
        
# %% Plots functions  
class plots(calc_fncs):
    """ Class that holds all of the functions for reducing the time resolved data"""
    def __init__(self,*args):
        super().__init__(*args)
        for a in args:
            self.__setattr__(str(a), args[0])
            
    
    def ratio_plot(data,ratio_buttons,regression_buttons,background_slider,ablation_slider,ablation_start_true,
                   ratio_plot_ylim_slider,displayframe
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
        max_list = []
        min_list = []
        for r in ratio_buttons:
            max_list.append(max(data[r]))
            min_list.append(min(data[r]))
        
        min_val = min(min_list)
        max_val = max(max_list)
        if displayframe == 'main':
            height,width=350,600
            fig = figure(height=height,width=width,title='Isotope Ratios',tools='pan,reset,save,wheel_zoom,xwheel_zoom,ywheel_zoom',toolbar_location='above',
                         x_axis_label='Time (s)',
                         y_range=[ratio_plot_ylim_slider[0],ratio_plot_ylim_slider[1]],
                         x_range=[min(data.Time_s)-3,max(data.Time_s)+3])
            
        else:
            height,width=500,500
            fig = figure(height=height,width=width,title='Isotope Ratios',tools='pan,reset,save,wheel_zoom,xwheel_zoom,ywheel_zoom',toolbar_location='above',
                         x_axis_label='Time (s)',
                         y_range=[min_val-min_val*0.1,max_val+max_val*0.1],
                         x_range=[min(data.Time_s)-3,max(data.Time_s)+3])
            
        fig.xgrid.grid_line_color = None
        fig.ygrid.grid_line_color = None
        var_cols = ratio_buttons # get the ratios selected by the user. These are used to get the columns from the calculated ratios
        # plot a line for each selected ratio
        for i,c in zip(var_cols,cycle(color_palette_regressions)):
            fig.line(data.Time_s,data[i],line_width=1,legend_label='{}'.format(i),color=c)
        
        if displayframe == 'main':
            pass
        else:
            pass
            regressions_to_plot = regression_buttons# get the regression type requested by the user
            # get regression parameters and stats for the requested regressions across the specified ablation period
            regressions,stats = calc_fncs.get_regressions(data,regression_buttons,ablation_start_true)
            # plot a line for each selected regression
            if regressions_to_plot is not None:
                for i,c in zip(regressions,cycle(color_palette_regressions)):
                    fig.line(data.Time_s,i,line_width=0.5,color=c)
                    
        fig.line([background_slider[0],background_slider[0]],[0,1],line_width=0.4,line_dash='dashed',color='black')
        fig.line([background_slider[1],background_slider[1]],[0,1],line_width=0.4,line_dash='dashed',color='black')
        fig.line([ablation_slider[0],ablation_slider[0]],[0,1],line_width=0.4,line_dash='solid',color='black')
        fig.line([ablation_slider[1],ablation_slider[1]],[0,1],line_width=0.4,line_dash='solid',color='black')
        fig.line([ablation_start_true,ablation_start_true],[0,1],line_width=0.4,line_dash='dotted',color='black')
        
        return fig
    
    
    def ratio_plot_76(data,background_slider,ablation_slider,ablation_start_true
                      ):
        """
        Function for displaying the plot of relavent isotopic ratios.
        
        Parameters
        ----------
        data : pandas dataframe
            pandas dataframe of observed data
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
    
        Returns
        -------
        fig : matplotlib figure
            plot that shows the relevant isotopic ratios for the entire measurement (background+ablation+washout)
    
        """
        
        fig = figure(height=350,width=600,title='207Pb/206Pb',tools='pan,reset,save,wheel_zoom,xwheel_zoom,ywheel_zoom',toolbar_location='above',
                     x_axis_label='Time (s)',
                     y_range=[min(data['207Pb/206Pb'])-0.04,max(data['207Pb/206Pb'])+0.04],
                     x_range=[min(data.Time_s)-3,max(data.Time_s)+3])
        fig.xgrid.grid_line_color = None
        fig.ygrid.grid_line_color = None
        # plot the 207Pb/206Pb ratio
        fig.line(data.Time_s,data['207Pb/206Pb'],line_width=1,legend_label='207Pb/206Pb',color='teal')
                    
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
        fig : bokeh fig
            Figure showing all of the ablation data.
    
        """
        # # ax.tick_params(axis='both',labelsize=14) # plot aesthetics
        if logdata == True:
            y_type = 'log'
        else:
            y_type = 'auto'
        max_list = []
        min_list = []
        for a in analytes_:
            max_list.append(max(data_ablation_[a]))
            min_list.append(min(data_ablation_[a]))
        
        min_val = min(min_list)
        max_val = max(max_list)
        
        fig = figure(height=350,width=1200,title='All Time Resolved Data',tools='pan,reset,save,wheel_zoom,xwheel_zoom,ywheel_zoom',toolbar_location='left',active_drag=None,
                      y_axis_type=y_type,x_axis_label = 'Time (s)', y_axis_label = 'Intensities (cps)',
                      y_range=[min_val-min_val*0.05,max_val+max_val*0.01],
                      x_range=[min(data_ablation_.Time_s),max(data_ablation_.Time_s)]
                      )
        
        var_cols=analytes_
        for i,c in zip(var_cols,cycle(color_palette)):
            fig.line(data_ablation_.Time_s,data_ablation_[i],line_width=0.7,legend_label='{}'.format(i),color=c)
        fig.line([background_slider[0],background_slider[0]],[1,1e8],line_width=0.4,line_dash='dashed',color='black')
        fig.line([background_slider[1],background_slider[1]],[1,1e8],line_width=0.4,line_dash='dashed',color='black')
        fig.line([ablation_slider[0],ablation_slider[0]],[1,1e8],line_width=0.4,line_dash='solid',color='black')
        fig.line([ablation_slider[1],ablation_slider[1]],[1,1e8],line_width=0.4,line_dash='solid',color='black')
        fig.line([ablation_start_true,ablation_start_true],[1,1e8],line_width=0.4,line_dash='dotted',color='black')
        
        return fig
        
        
    def residuals_plot(data,regression_buttons,start_ablation,stop_ablation,ablation_start_true):
        fig206 = figure(height=250,width=250,title='206Pb/238U Residuals',tools='pan,reset,save,wheel_zoom,xwheel_zoom,ywheel_zoom',toolbar_location='above',
                        x_axis_label='Fitted Value',y_axis_label='Residuals')
        fig207 = figure(height=250,width=250,title='207Pb/235U Residuals',
                        x_axis_label='Fitted Value',y_axis_label='Residuals')
        
        
        fitted,fitted207,residuals,residuals207 = calc_fncs.get_regressions(data,regression_buttons,ablation_start_true) # get fitted regression values, regression names, and residuals
        # plot the residuals for the regressions and the fitted value
        for j,c in zip(range(0,len(regression_buttons)),color_palette_regressions):
            fig206.scatter(fitted[j],residuals[j],marker='circle',size=8,color=c,legend_label='{}'.format(regression_buttons[j]))
            fig207.scatter(fitted207[j],residuals207[j],marker='circle',size=8,color=c,legend_label='{}'.format(regression_buttons[j]))
        fig206.line([min(fitted[j]),max(fitted[j])],[0,0],color='black',line_width=0.4)
        fig207.line([min(fitted207[j]),max(fitted207[j])],[0,0],color='black',line_width=0.4)
        fig206.legend.label_text_font_size='8pt'
        fig207.legend.label_text_font_size='8pt'
        fig = gridplot([[fig206],[fig207]],width=500,height=250)
        
        return fig
    
    def ellipse_plot(data,power,start_ablation,stop_ablation,ablation_start_true,lock_ablation_start_true,regression_buttons,counts_mode):
        ell1p,ell2p,x1,y1,y2 = calc_fncs.get_ellipse(data, power, ablation_start_true, regression_buttons)
        
        ell1 = Ellipse(xy=ell1p[0],width=ell1p[1],height=ell1p[2],angle=ell1p[3],color='steelblue',ec='k',alpha=0.5) # set the parameters into a plotable 'patch'
        ell2 = Ellipse(xy=ell2p[0],width=ell2p[1],height=ell2p[2],angle=ell2p[3],color='steelblue',ec='k',alpha=0.5)
        fig1 = Figure(figsize=(3,3)) # create a figure that is 8Wx8H
        fig2 = Figure(figsize=(3,3))
        ax1 = fig1.add_subplot() # add an axis to the mpl figure
        ax2 = fig2.add_subplot()
        
        ax1.add_artist(ell1) # adde the ellipsoid patch to the axis
        ax1.plot(x1,y1,'.k',markersize=1) # plot individual observations as dots
        ax2.add_artist(ell2)
        ax2.plot(1/y1,y2,'.k',markersize=1)
        
        ax1.set_xlabel('$^{207}$Pb/$^{235}$U',fontsize=6) # set xlabel
        ax1.set_ylabel('$^{206}$Pb/$^{238}$U',fontsize=6) # set ylabel
        ax2.set_xlabel('$^{238}$U/$^{206}$Pb',fontsize=6)
        ax2.set_ylabel('$^{207}$Pb/$^{206}$Pb',fontsize=6)
        ax1.tick_params(axis='both',labelsize=5) # set tick parameters on the axes
        ax2.tick_params(axis='both',labelsize=5)
        
        ax1.set_xlim(min(x1)-min(x1)*0.1,max(x1)+max(x1)*0.1)
        ax1.set_ylim(min(y1)-min(y1)*0.1,max(y1)+max(y1)*0.1)
        ax2.set_xlim(min(1/y1)-min(1/y1)*0.1,max(1/y1)+max(1/y1)*0.1)
        ax2.set_ylim(min(y2)-min(y2)*0.1,max(y2)+max(y2)*0.1)
        
        fig1.tight_layout() 
        fig2.tight_layout()
        
        return fig1,fig2
    

# %% Display Plots and Widgets
class make_plots(param.Parameterized):
    """ class that parameterizes inputs and sends them to the above functions to be rendered in a GUI"""
    sample_subset = param.Selector(objects=[])
    
    update_output_button = param.Action(lambda x: x.evaluate_output_data(),label='Evaluate Interval')
    accept_interval_button = param.Action(lambda x: x.send_reduction(),label='Accept Interval') # button that triggers sample name to be accepted and ablation to be reduced
    export_data_button = param.Action(lambda x: x.export_data(),label='Export')
    
    lock_ablation_start_true = param.Boolean(True,label='Lock Back Projection')
    ablation_start_true = param.Number(24)
    ablation_slider = param.Range(default=(24,37),bounds=(0,100))
    background_slider = param.Range(default=(3,19),bounds=(0,100))
    
    ratio_buttons = param.ListSelector(default=['206Pb/238U'], objects=['206Pb/238U','206Pb/204Pb','207Pb/235U','207Pb/206Pb','238U/235U'])
    
    regression_buttons = param.ListSelector(default=['1st Order'], objects=['1st Order','Exp. Regression'])
    
    ablation_plot_ylim_slider = param.Range(default=(0,1000000),bounds=(0,10000000))
    ratio_plot_ylim_slider = param.Range(default=(0.0,0.3),bounds=(-0.05,3))
    # ratio_plot_xlim_slider = param.Range(default=(0,60),bounds=(0,100))
    
    logcountsdata = param.Boolean(False,label='Log Intensities')
    analytes_ = param.ListSelector(default=['202Hg','204Pb','206Pb','207Pb','208Pb','232Th','235U','238U'],objects=['202Hg','204Pb','206Pb','207Pb','208Pb','232Th','235U','238U'])
    
    counts_mode = param.Selector(default='Means & LIEF',objects=['Total Counts','Means & LIEF'])
    ratio_type = param.Selector(default='Ratio of Means',objects=['Ratio of Means','Geometric'])
    
    ellipsemode_selector = param.Boolean(True,label='Generate Ellipse')
    power = param.Number(default=0.05)
    
    arrayofdwelltimes = param.Array()
    
    input_data = param.DataFrame(precedence=-1)
    file_path = param.String(default='Insert File Path')
    output_data = param.DataFrame(precedence=-1)
    
        
    def __init__(self,**params):
        super().__init__(**params)
        self.file_input_widget = pn.Param(self.param.input_data)
        self.output_data_widget = pn.Param(self.param.output_data)
        self.widgets = pn.Param(self,parameters=['sample_subset','update_output_button','export_data_button',
                                                 'ablation_start_true','lock_ablation_start_true','ablation_slider','background_slider','ablation_plot_ylim_slider','logcountsdata','analytes_',
                                                 'ratio_plot_ylim_slider','ratio_plot_ylim_max','ratio_buttons','regression_buttons','counts_mode','ratio_type',
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
            self.output_data.insert(5,'b start',0)
            self.output_data.insert(6,'b end',0)
            self.output_data.drop('Time',axis=1)
            analytelength = len(self.input_data.columns)-3
            dwellarray = np.full_like(np.arange(analytelength,dtype=float),0.01)
            self.arrayofdwelltimes = dwellarray
        else:
            return print('No Input Data Available')
    
        
    @pn.depends('input_data','background_slider','ablation_slider','ablation_start_true','lock_ablation_start_true','sample_subset','ablation_plot_ylim_slider')
    def call_ablation_plot(self):
        if self.sample_subset is not None:
            if self.lock_ablation_start_true == True:
                self.ablation_start_true = self.ablation_slider[0]
            data_toplot = self.input_data[self.input_data['SampleLabel'] == self.sample_subset]
            return plots.ablation_plot(data_toplot,self.background_slider,self.ablation_slider,self.ablation_start_true,
                                       self.ablation_plot_ylim_slider,self.logcountsdata,self.analytes_)
        else:
            pass
        
    @pn.depends('input_data','ratio_buttons','regression_buttons',
                'background_slider','ablation_slider','ablation_start_true','lock_ablation_start_true','ratio_plot_ylim_slider',
                'sample_subset')
    def call_ratio_plot(self):
        if self.sample_subset is not None:
            if self.lock_ablation_start_true == True:
                self.ablation_start_true = self.ablation_slider[0]
            data_toplot = self.input_data[self.input_data['SampleLabel'] == self.sample_subset]
            t_resolved_ratio_data = calc_fncs.get_tresolved_ratios(data_toplot)
            return plots.ratio_plot(t_resolved_ratio_data,self.ratio_buttons,self.regression_buttons,self.background_slider,
                                    self.ablation_slider,self.ablation_start_true,self.ratio_plot_ylim_slider,'main')
    
                  
    @pn.depends('input_data','background_slider','ablation_slider','ablation_start_true','lock_ablation_start_true','sample_subset')
    def call_ratio_plot76(self):
        if self.sample_subset is not None:
            if self.lock_ablation_start_true == True:
                self.ablation_start_true = self.ablation_slider[0]
            data_toplot = self.input_data[self.input_data['SampleLabel'] == self.sample_subset]
            t_resolved_ratio_data = calc_fncs.get_tresolved_ratios(data_toplot)
            return plots.ratio_plot_76(t_resolved_ratio_data,self.background_slider,self.ablation_slider,self.ablation_start_true)
    
    
    
    @pn.depends('input_data','ratio_buttons','regression_buttons','ratio_plot_ylim_slider',
                'background_slider','ablation_slider','ablation_start_true','lock_ablation_start_true',
                'power','counts_mode','arrayofdwelltimes')
    def evaluate_output_data(self, event=None):
        """
        Function that clears any residuals input from uploading data or previously approved analyses, then generates fresh ones in a modal

        Parameters
        ----------
        event : open panel modal

        """
        if self.ablation_start_true == True:
            self.ablation_start_true = self.ablation_slider[0]
            
        data_toevaluate = self.input_data[self.input_data['SampleLabel'] == self.sample_subset].reset_index(drop=True)
        counts_data,backgrounds,lods = calc_fncs.backgroundsubtract_convert_lod(data_toevaluate,self.background_slider[0],self.background_slider[1],self.ablation_slider[0],self.ablation_slider[1],self.arrayofdwelltimes)
        counts_data = counts_data.fillna(0)
        t_resolved_ratios = calc_fncs.get_tresolved_ratios(counts_data) # get time resovled ratios of the ablation interval
        
        ratioplot_pane = pn.pane.Bokeh(row(plots.ratio_plot(t_resolved_ratios,self.ratio_buttons,self.regression_buttons,
                                                            self.background_slider,self.ablation_slider,self.ablation_start_true,
                                                            self.ratio_plot_ylim_slider,'modal')))

        residuals_plot_pane = pn.pane.Bokeh(row(plots.residuals_plot(t_resolved_ratios,self.regression_buttons,self.ablation_slider[0],self.ablation_slider[1],
                                                                     self.ablation_start_true)))
        
        Weth_ell,TW_ell = plots.ellipse_plot(t_resolved_ratios, self.power, self.ablation_slider[0], self.ablation_slider[1], self.ablation_start_true, 
                                             self.lock_ablation_start_true, self.regression_buttons, self.counts_mode)
        ellipse_tabs = pn.Tabs(('TW',TW_ell),('Weth.',Weth_ell),dynamic=True)
        
        
        regressions,stats = calc_fncs.get_regressions(t_resolved_ratios,self.regression_buttons,self.ablation_start_true)
        divider1 = pn.pane.Markdown('206/238 Reg. Stats:')
        r2_1 = pn.pane.Markdown('$$R^{2}$$ = '+str(round(stats[0],3)))
        r2_exp = pn.pane.Markdown('$$R^{2}_{exp}$$ = '+str(round(stats[1],3)))
        SE_1stper = pn.pane.Markdown('$$1^{st} Order SE %% $$ = '+str(round(stats[2],3)))
        SE_expper = pn.pane.Markdown('$$Exp SE %% $$ = '+str(round(stats[3],3)))
        divider2 = pn.pane.Markdown('207/235 Reg. Stats:')
        SE_1stper_207 = pn.pane.Markdown('$$1^{st} Order SE %% $$ = '+str(round(stats[4],3)))
        SE_expper_207 = pn.pane.Markdown('$$Exp SE %% $$ = '+str(round(stats[5],3)))
         
        stats_markdown = pn.Row(pn.Column(divider1,r2_1,r2_exp,SE_1stper,SE_expper),pn.Column(divider2,SE_1stper_207,SE_expper_207))
        
        
        fastgrid_layout.modal[0].clear()
        fastgrid_layout.modal[1].clear()
        fastgrid_layout.modal[2].clear()
        fastgrid_layout.modal[3].clear()
        fastgrid_layout.modal[4].clear()
        
        fastgrid_layout.modal[0].append(buttons_sample)
        fastgrid_layout.modal[1].append(pn.Row(ratioplot_pane,pn.Spacer(width=50, height=500),residuals_plot_pane))
        fastgrid_layout.modal[2].append(pn.Row(ellipse_tabs,pn.Spacer(width=50, height=500),stats_markdown))
        
        
        fastgrid_layout.open_modal()
        
    
    def send_reduction(self,event=None):
        """
        Function that gets the fully reduced data and sends it to the output data file that will be exported. Closes modal.

        """
        if self.lock_ablation_start_true == True:
            self.ablation_start_true = self.ablation_slider[0]
        # get the approved data by calling functions
        data_tosend = self.input_data[self.input_data['SampleLabel'] == self.sample_subset]
        data_approved = calc_fncs.get_approved(data_tosend,
                                               self.background_slider[0],self.background_slider[1],self.ablation_slider[0],self.ablation_slider[1],self.ablation_start_true,
                                               self.regression_buttons,self.ratio_type,self.counts_mode,self.arrayofdwelltimes,self.param.sample_subset.objects,self.power
                                               )
        # if there is no data in the first column (assuming you can see 238U) assign the first sample measurement as -1 and make this global measurement (measurmentindex) 1
        # otherwise, add 1 to the last global measurement number and to the last sample measurement number
        if self.output_data is None:
            self.output_data = data_approved
        else:
            self.output_data = pd.concat([self.output_data,data_approved],ignore_index=True)
    
        fastgrid_layout.close_modal()
        pn.state.notifications.success('Interval Reduced',duration=2000)
    
    
    @pn.depends('output_data',watch=True)
    def _update_output_widget(self):
        if self.output_data is not None:
            self.output_data_widget = self.output_data
            self.output_data_widget.height = 400
            self.output_data_widget.heightpolicy = 'Fixed'
            return pn.widgets.Tabulator(self.output_data_widget,width=600) # use 600 for large screen, 100-150 for small screen 
        
            
    @pn.depends('output_data')
    def export_data(self,event=None):
        self.output_data = self.output_data.drop(columns=['Time','Time_s'])
        self.output_data = self.output_data.drop(0,axis=0)
        self.output_data.to_excel('output_lasertramZ.xlsx',index=False)
        

# %% Initialize and call app

callapp = make_plots(name='Reduce Ablation Data')

pn.extension('tabulator','mathjax',notifications=True)
pn.state.notifications.position = 'bottom-right'

buttons_sample=pn.WidgetBox(pn.Param(callapp.param.accept_interval_button,
                               widgets={'accept_interval_button': pn.widgets.Button(name='Accept Interval',button_type='success')}))

widgets = {'ratio_buttons': pn.widgets.CheckBoxGroup,
            'regression_buttons': pn.widgets.CheckBoxGroup,
           'counts_mode': pn.widgets.RadioButtonGroup,
           'ratio_type': pn.widgets.RadioButtonGroup,
           'export_data_button': pn.widgets.Button(name='DDDT!',button_type='success'),
           'analytes_':pn.widgets.CheckBoxGroup,
           'ablation_plot_ylim_slider':pn.widgets.EditableRangeSlider(name='Ablation Plot Y-lim',start=0,end=1e6,value=(0,1e6),step=0.01),
           'ratio_plot_ylim_slider':pn.widgets.EditableRangeSlider(name='Ratio Plot Y-lim',start=0,end=2,value=(0,0.3),step=0.1)
           }

fastgrid_layout = pn.template.VanillaTemplate(title='LaserTRAMZ: LA-Q-ICP-MS',
                                                sidebar=pn.Column(pn.WidgetBox(pn.Param(callapp.param,widgets=widgets))),sidebar_width=380)

fastgrid_layout.modal.append(pn.Row())
fastgrid_layout.modal.append(pn.Row())
fastgrid_layout.modal.append(pn.Row())
fastgrid_layout.modal.append(pn.Row())
fastgrid_layout.modal.append(pn.Row())
fastgrid_layout.modal.append(pn.Row())

fastgrid_layout.main.append(pn.Column(callapp.call_ablation_plot,pn.Row(callapp.call_ratio_plot,callapp.call_ratio_plot76))) # for vanilla
fastgrid_layout.main.append(pn.Column(callapp._update_output_widget))
fastgrid_layout.show();









    
