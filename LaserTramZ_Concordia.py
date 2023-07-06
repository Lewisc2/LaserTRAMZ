#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 15:15:39 2023

@author: Chuck Lewis, Oregon State University
"""

# %%
from scipy import stats
import scipy
from patsy import dmatrices
import statsmodels.api as sm
import sys
import param
import statistics
import panel as pn
from matplotlib.patches import Ellipse
from matplotlib.figure import Figure
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('agg')

# %%
# define constants for reducing U-Pb data
lambda_238 = 1.55125e-10
lambda_235 = 9.8485e-10
lambda_232 = 4.9475e-11
lambda_230 = 9.158e-6

# %%


class calc_fncs:
    """ Class that holds all of the functions for reducing the reduced LAICPMS data"""

    def __init__(self, *args):
        for a in args:
            self.__setattr__(str(a), args[0])

    def get_driftstats(df, analytes):
        df = df.replace('bdl', 0)
        slope_array = []
        intercept_array = []
        f_p_array = []
        p_slope_array = []
        for i in analytes:
            y = pd.to_numeric(df[i])
            X = pd.to_numeric(df['measurementindex'])
            X1 = sm.add_constant(X)
            linmod = sm.OLS(y, X1).fit()
            f_p = linmod.f_pvalue
            p_slope = linmod.pvalues[1]
            if f_p < 0.01 and p_slope < 0.01:
                slope_array.append(linmod.params[1])
                intercept_array.append(linmod.params[0])
                f_p_array.append(f_p)
                p_slope_array.append(p_slope)
            else:
                slope_array.append(0)
                intercept_array.append(0)
                f_p_array.append('-')
                p_slope_array.append('-')
        return slope_array, intercept_array, f_p_array, p_slope_array
    
    def calc_pbpb_bias(df,selected_bias_sample):
        if selected_bias_sample == 'NIST614':
            biasdf = df[df['SampleLabel'].str.contains(selected_bias_sample)]
            biasdf = biasdf.reset_index(drop=True)
            biasdf = biasdf.replace('bdl',0)
            nist614_206204 = 17.84
            nist614_207206 = 0.8707
            bias_206204 = nist614_206204/np.mean(biasdf['206Pb/204Pb'])
            bias_207206 = nist614_207206/np.mean(biasdf['207Pb/206Pb'])
        else:
            bias_206204 = 1
            bias_207206 = 1
        return bias_206204,bias_207206
        

    def get_data_TW_regressions(df):
        """
        """
        pts_x = np.array(
            [df['238U/206Pb'], np.zeros_like(df['SK 207Pb/206Pb'])]).T
        pts_y = np.array([df['207Pb/206Pb'], df['SK 207Pb/206Pb']]).T
        discordia_t = np.zeros((len(df), 2))

        for i in range(0, len(df)):
            discordia_t[i] = np.poly1d(np.polyfit(pts_x[i], pts_y[i], 1))

        return discordia_t

    def get_TW_concordia():
        """
        """
        # calculate T-W discordia - put in plot function
        t = np.linspace(1, 4.6e9, 100000)
        u238_pb206 = np.zeros(len(t))
        pb207_206r = np.zeros(len(t))
        for i in range(0, len(t)):
            u238_pb206[i] = 1/(np.exp(lambda_238*t[i])-1)
            pb207_206r[i] = (
                1/137.88) * ((np.exp(lambda_235*t[i])-1) / (np.exp(lambda_238*t[i])-1))

        return u238_pb206, pb207_206r

    def get_projections(df, ellipse_mode_selector, power):
        """
        """
        # get the TW concordia values
        x_TW, y_TW = calc_fncs.get_TW_concordia()
        discorida_regressions = calc_fncs.get_data_TW_regressions(
            df)  # get regressions
        # array of xvalues to project over
        x_vals = np.linspace(min(x_TW), max(x_TW), 100000)
        # set up array to be filled with calculated radiogenic lead component
        pts_pb_r = np.zeros(len(discorida_regressions))
        concordia_238_206 = np.zeros(len(discorida_regressions))

        for i in range(0, len(discorida_regressions)):

            discordia_207206 = discorida_regressions[i][0] * \
                x_vals+discorida_regressions[i][1]
            discordia_238206 = (
                discordia_207206-discorida_regressions[i][1])/discorida_regressions[i][0]

            # distance of y value from line
            delta_y = (discorida_regressions[i][1] +
                       x_TW * discorida_regressions[i][0]) - y_TW
            # index the regression for where the curve cross the regression
            indx = np.where(delta_y[1:]*delta_y[:-1] < 0)[0]
            # similar triangles geometry gets points
            d_ratio = delta_y[indx] / (delta_y[indx] - delta_y[indx + 1])
            # empty array for crossing points
            points = np.zeros((len(indx), 2))
            points[:, 0] = x_TW[indx] + d_ratio * \
                (x_TW[indx+1] - x_TW[indx])  # x crossings
            points[:, 1] = y_TW[indx] + d_ratio * \
                (y_TW[indx+1] - y_TW[indx])  # y crossings
            y_point = y_TW[indx] + d_ratio * \
                (y_TW[indx+1] - y_TW[indx])  # y crossings
            x_point = x_TW[indx] + d_ratio * \
                (x_TW[indx+1] - x_TW[indx])  # x crossings

            if len(y_point) >= 1:
                pts_pb_r[i] = min(y_point)
            elif len(y_point) < 0:
                pts_pb_r[i] = 0

            if len(x_point) > 1:
                concordia_238_206[i] = min(x_point)
            elif len(x_point) == 1:
                concordia_238_206[i] = x_point

        return points, concordia_238_206, pts_pb_r, discordia_207206, discordia_238206

    def get_ellipse(df, power):
        callingmethod = sys._getframe().f_back.f_code.co_name
        if callingmethod == 'correct_sample_ages':
            df.replace([np.inf, -np.inf], np.nan, inplace=True)
            df.dropna(inplace=True)
            x2 = 1/df['206Pb/238U Corrected']
        else:
            x2 = df['238U/206Pb']
        y2 = df['207Pb/206Pb']

        cov2 = np.cov(x2, y2)
        eigval2, eigvec2 = np.linalg.eig(cov2)
        order2 = eigval2.argsort()[::-1]
        eigvals_order2 = eigval2[order2]
        eigvecs_order2 = eigvec2[:, order2]

        c2 = (np.mean(x2), np.mean(y2))
        wid2 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power), df=2)
                         * eigvals_order2[0])
        hgt2 = 2*np.sqrt(scipy.stats.chi2.ppf((1-power), df=2)
                         * eigvals_order2[1])
        theta2 = np.degrees(np.arctan2(*eigvecs_order2[:, 0][::-1]))

        ell2_params = [c2, wid2, hgt2, theta2]

        return ell2_params

    def plot_TW(df, xlim_start, xlim_stop, ylim_start, ylim_stop, label_toggle, ellipse_mode_selector, power):
        """
        """
        plt.style.use('seaborn-colorblind')
        fig = Figure(figsize=(10, 14))
        ax = fig.add_subplot()

        x_TW, y_TW = calc_fncs.get_TW_concordia()

        disc_regressions = calc_fncs.get_data_TW_regressions(df)

        x_vals = np.linspace(min(x_TW), max(x_TW), 100000)
        pts_pb_r = np.zeros(len(disc_regressions))

        for i in range(0, len(disc_regressions)):
            discordia_207206 = disc_regressions[i][0] * \
                x_vals+disc_regressions[i][1]
            discordia_238206 = (
                discordia_207206-disc_regressions[i][1])/disc_regressions[i][0]

            # distance of y value from line
            delta_y = (disc_regressions[i][1] +
                       x_TW * disc_regressions[i][0]) - y_TW
            # index the regression for where the curve cross the regression
            indx = np.where(delta_y[1:]*delta_y[:-1] < 0)[0]
            # similar triangles geometry gets points
            d_ratio = delta_y[indx] / (delta_y[indx] - delta_y[indx + 1])
            # empty array for crossing points
            points = np.zeros((len(indx), 2))
            points[:, 0] = x_TW[indx] + d_ratio * \
                (x_TW[indx+1] - x_TW[indx])  # x crossings
            points[:, 1] = y_TW[indx] + d_ratio * \
                (y_TW[indx+1] - y_TW[indx])  # y crossings
            y_point = y_TW[indx] + d_ratio * \
                (y_TW[indx+1] - y_TW[indx])  # y crossings
            if len(y_point) >= 1:
                pts_pb_r[i] = min(y_point)
            elif len(y_point) < 0:
                pts_pb_r[i] = 0

            if ellipse_mode_selector == 'Point Estimates':
                ax.plot(discordia_238206, discordia_207206, '-k', lw=0.5)
                ax.plot(points[:, 0], points[:, 1], 'og', lw=1, mec='k')
            else:
                pass
        if ellipse_mode_selector == 'Point Estimates':
            ax.plot(df['238U/206Pb'], df['207Pb/206Pb'], 'kd')
            ax.errorbar(df['238U/206Pb'], df['207Pb/206Pb'], xerr=df['206/238 Reg. err'],
                        yerr=df['SE 207/206'], fmt='none', ecolor='k', elinewidth=0.5)
            ax.plot(x_TW, y_TW, 'k', lw=1)
        elif ellipse_mode_selector == 'Ellipses':
            for i in df.SampleLabel.unique():
                conf_ellipse = calc_fncs.get_ellipse(
                    df[df['SampleLabel'] == i], power)
                ell = Ellipse(conf_ellipse[0], conf_ellipse[1], conf_ellipse[2],
                              conf_ellipse[3], color='darkslategray', alpha=0.3, ec='k')
                ax.add_artist(ell)
            ax.plot(x_TW, y_TW, 'k', lw=1)
        ax.set_xlabel('$^{238}$U / $^{206}$Pb$^{*}$', fontsize=12)
        ax.set_ylabel('[$^{207}$Pb / $^{206}$Pb]$^{*}$', fontsize=12)
        ax.tick_params(axis='both', labelsize=10)
        ax.set_ylim(ylim_start, ylim_stop)
        ax.set_xlim(xlim_start, xlim_stop)

        # zip joins x and y coordinates in pairs
        if ellipse_mode_selector == 'Point Estimates':
            if 'Concordia' in label_toggle:
                for x, y, t in zip(df['238U/206Pb'], df['207Pb/206Pb'], df['SampleLabel']):

                    label = t

                    ax.annotate(label,  # this is the text
                                # these are the coordinates to position the label
                                (x, y),
                                textcoords="offset points",  # how to position the text
                                # distance from text to points (x,y)
                                xytext=(0, 10),
                                ha='center',
                                fontsize=8)  # horizontal alignment can be left, right or center
        else:
            pass

        return fig

    def plot_boxplot(ages, analysis_ID, label_toggle, ellipse_mode_selector):
        """
        """
        if ellipse_mode_selector == 'Point Estimates':
            plt.style.use('seaborn-colorblind')
            fig = Figure(figsize=(4, 10))
            ax = fig.add_subplot()

            bp = ax.boxplot(ages, patch_artist=True, boxprops=dict(facecolor='slategray', color='k'),
                            medianprops=dict(color='limegreen'), meanprops=dict(marker='d', mfc='limegreen', mec='k', markersize=4),
                            flierprops=dict(
                                marker='o', mfc='None', mec='k', markersize=4),
                            showmeans=True)

            ax.text(0.05, 0.8, 'Mean ='+str(round(ages.mean(), 2)),
                    fontsize=12, transform=ax.transAxes)
            ax.text(0.05, 0.7, 'Med ='+str(round(ages.median(), 2)),
                    fontsize=12, transform=ax.transAxes)
            ax.text(0.05, 0.6, 'Min ='+str(round(ages.min(), 2)),
                    fontsize=12, transform=ax.transAxes)
            ax.text(0.05, 0.5, 'Max ='+str(round(ages.max(), 2)),
                    fontsize=12, transform=ax.transAxes)
            ax.text(0.05, 0.4, 'n = '+str(len(ages)),
                    fontsize=12, transform=ax.transAxes)

            ax.set_ylabel('Age (Ma)', fontsize=12)
            ax.set_xlabel(' ', fontsize=1)
            ax.tick_params(axis='both', labelsize=8)

            if 'Box + Whisker' in label_toggle:
                fliers = [item.get_ydata() for item in bp['fliers']]
                for x, t in zip(ages, analysis_ID):

                    if t in fliers:
                        label = t
                        age = x
                        ax.annotate(
                            label, age, textcoords='offset points', ha='center', fontsize=5)

            return fig
        else:
            return print('In Ellipse Mode')

    def plot_drift(df, analyte, drift_slope, drift_intercept, drift_fp, drift_pslope):
        plt.style.use('seaborn-colorblind')
        fig = Figure(figsize=(12, 4))
        ax = fig.add_subplot()

        yvals = df[analyte].replace('bdl', 0)
        ax.plot(df['measurementindex'], yvals, 'd',
                mfc='forestgreen', mec='k', ms=10,label='Observed')
        regression_y = df[analyte] - df['measurementindex']*drift_slope
        regression_x = df['measurementindex']
        if drift_slope != 0:
            ax.plot(regression_x, regression_y, 's',mfc='steelblue',mec='k', lw=0.6,label='Drift corrected',ms=10)
            ax.text(min(df['measurementindex']), min(df[analyte]),
                    'Slope: '+str(drift_slope)+'\n' + 'Slope pval: '+str(drift_pslope)+'\n'+'F pval:'+str(drift_fp)+'\n'+'Intercept'+str(drift_intercept))
        else:
            pass
        ax.legend(loc='best')
        ax.set_xlabel('Measurement #')
        ax.set_ylabel(str(analyte))

        return fig

    def correct_standard_ages(df, std_txt, Pb_Th_std_crct_selector, regression_selector, ellipse_mode_selector, power):
        """
        function to calculate Pb and fractionation factor corrected ages
        """
        df = df.reset_index(drop=True)
        # create a dictionary that holds known or estimated U/Th ratios of zircon and associated magma for standards, as well as common Pb ratios
        stds_dict = {'Temora': [2.4, 0.79200, 18.0528, 15.5941],  # Black et al., 2004. Ambiguous Th correction > assumed D = 0.33
                     # Schmitz and Bowring 2001; only one with measured common Pb so far
                     'FishCanyon': [1.496558, 0.454545, 18.4275, 15.5425],
                     '94-35': [1, 0.33, 18.6191, 15.626],
                     # Slama et al 2008
                     'Plesovice': [10.7, 0.25, 18.1804, 15.6022],
                     # Black et al., 2004. Ambiguous Th correction > assumed D = 0.33
                     'R33': [1.4, 0.46200, 18.0487, 15.5939],
                     # Wiedenbeck et al 1995
                     '91500': [1, 0.33, 16.9583, 15.4995],
                     # Paces and Miller 1993. Ambiguous Th correction > assumed D = 0.33
                     'FC1': [1.7, 0.56100, 16.892, 15.492],
                     # unpublished; Bowring > assumed D = 0.33
                     'Oracle': [2.2, 0.725999, 16.2726, 15.4099],
                     # Pecha unpublished > assumed D = 0.33
                     'Tan-Bra': [1.2, 0.39600, 14.0716, 14.8653],
                     'OG1': [1.3, 0.42900, 11.8337, 13.6071]}  # Stern et al 2009. Ambiguous Th correction > assumed D = 0.33
        # get points needed for correction
        points, concordia_238_206, pts_pb_r, na, naII = calc_fncs.get_projections(
            df, ellipse_mode_selector, power)

        # set up empty arrays to be filled for common lead corrections
        common_filter = []

        pb_m = df['207Pb/206Pb']  # measured 207/206
        if Pb_Th_std_crct_selector == 'Common Pb':
            common = df['SK 207Pb/206Pb']
        elif Pb_Th_std_crct_selector == 'Common Pb + Th Disequil.':
            # calculated Stacey-Kramers 207/206 overwriting the input 207/206 should manual values be requested
            df['SK 207Pb/206Pb'] = stds_dict.get(std_txt)[3] / \
                stds_dict.get(std_txt)[2]
            common = df['SK 207Pb/206Pb']

        for i in common:
            if i <= 0:
                common_filter.append(0)
            else:
                common_filter.append(i)

        # calculate fraction of common Pb
        f_ = (pb_m - pts_pb_r) / (common - pts_pb_r)
        # set up array to set f = 0 if point lies on or below Concordia (i.e., no common Pb present)
        f = []

        for k, j in zip(common_filter, f_):
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
        df['206Pb/238Upbc_numerical'] = 1 / \
            df['238U/206Pb']-(1/df['238U/206Pb']*f)
        df['206Pb/238Upbc'] = 1/concordia_238_206
        df['206Pb/238Uc_age'] = np.log(df['206Pb/238Upbc'] + 1) / lambda_238
        UThstd, UThstd_rx = stds_dict.get(
            std_txt)[0], stds_dict.get(std_txt)[1]
        DThU = (1/UThstd)/(1/UThstd_rx)
        df['206Pb/238UPbThc'] = df['206Pb/238Upbc'] - \
            (lambda_238/lambda_230*(DThU-1))
        df['206Pb/238UPbTh_age'] = np.log(
            df['206Pb/238UPbThc'] + 1) / lambda_238
        UTh_std = stds_dict.get(std_txt)[0]
        UTh_std_m = df['238U'].mean()/df['232Th'].mean()

        avg_std_age = df['206Pb/238Uc_age'].mean()
        avg_std_age_Thcrct = df['206Pb/238UPbTh_age'].mean()

        avg_std_ratio = 1/df['238U/206Pb'].mean()
        avg_std_ratio_Thcrct = df['206Pb/238UPbThc'].mean()

        wts = []
        if ellipse_mode_selector == 'Point Estimates':
            for i in range(0, len(df)):
                if regression_selector == '1st Order':
                    wt_se_i = 1/df['SE 206/238 1st Order'][i]
                    wts.append(wt_se_i)
                elif regression_selector == '2nd Order':
                    wt_se_i = 1/df['SE 206/238 2nd Order'][i]
                    wts.append(wt_se_i)
                elif regression_selector == 'Exp.':
                    wt_se_i = 1/df['SE 206/238 Exp'][i]
                    wts.append(wt_se_i)

                # if regression_selector == '1st Order':
                #     avg_reg_err = df['SE 206/238 1st Order'].mean()
                # elif regression_selector == '2nd Order':
                #     avg_reg_err = df['SE 206/238 2nd Order'].mean()
                # elif regression_selector == 'Exp.':
                #     avg_reg_err = df['SE 206/238 Exp'].mean()
        else:
            pass

        avg_reg_err = 1/np.sum(wts)

        return avg_std_age, avg_std_age_Thcrct, avg_std_ratio, avg_std_ratio_Thcrct, avg_reg_err, UTh_std, UTh_std_m

    def correct_sample_ages(df, std, std_txt, ThU_zrn, ThU_magma, Pb_Th_std_crct_selector, regression_selector, ellipse_mode_selector, power, UTh_std_norm):
        """
        function to calculate Pb and fractionation factor corrected ages
        """

        if ellipse_mode_selector == 'Point Estimates':
            # get points needed for correction
            points, concordia_238_206, pts_pb_r, n, a = calc_fncs.get_projections(
                df, ellipse_mode_selector, power)

            # set up empty arrays to be filled for common lead corrections
            common_filter = []

            pb_m = df['207Pb/206Pb']  # measured 207/206
            common = df['SK 207Pb/206Pb']

            for i in common:
                if i <= 0:
                    common_filter.append(0)
                else:
                    common_filter.append(i)

            # calculate fraction of common Pb
            f_ = (pb_m - pts_pb_r) / (common - pts_pb_r)
            # set up array to set f = 0 if point lies on or below Concordia (i.e., no common Pb present)
            f = []

            for k, j in zip(common_filter, f_):
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
            df['206Pb/238Upbc_numerical'] = 1 / \
                df['238U/206Pb']-(1/df['238U/206Pb']*f)
            df['206Pb/238Upbc'] = 1/concordia_238_206

            frac_factor, tims_age, tims_error, std_avg, avg_std_age_Thcrct, avg_std_ratio, avg_std_ratio_Thcrct, std_avgerr, UTh_std, UTh_std_m = calc_fncs.get_standard_fracfctr(
                std, std_txt, Pb_Th_std_crct_selector, regression_selector, ellipse_mode_selector, power)

            if Pb_Th_std_crct_selector == 'Common Pb':
                df['206Pb/238Uc_age'] = np.log(
                    df['206Pb/238Upbc'] + 1) / lambda_238
                df['206Pb/238U_correctedage'] = df['206Pb/238Uc_age']*frac_factor
                # propagate errors
                # error on fractionation factor. includes error from ID-TIMS and ICPMS
                dfrac = np.abs(frac_factor)*(((tims_error/2)/tims_age)
                                             ** 2 + (std_avgerr/avg_std_ratio)**2)**(1/2)
                dage = np.abs(df['206Pb/238U_correctedage'])*((1/(1/df['238U/206Pb']))
                                                              ** 2*df['206/238 Reg. err']**2 + (0.16/2/100)**2)**(1/2)
                # error on age equation. error on decay constant includes 1.5* counting stats (Mattinson 1987)
                # error on estimation for common lead using 207 method. Uses conservaitve estimates of 1.0 for 206/204 and 0.3 for 207/204 (Mattionson, 1987)
                # total propagated error
                dagetot = np.abs(df['206Pb/238U_correctedage'])*(dfrac**2 +
                                                                 (df['206/238 Reg. err']/(1/df['238U/206Pb']))**2 + (0.16/2/100)**2 +
                                                                 ((0.3/2)/df['SK 207Pb/204Pb'])**2 + ((1/2)/df['SK 206Pb/204Pb'])**2 +
                                                                 (df['SE 207/206']/df['207Pb/206Pb'])**2)**(1/2)

                df['∆206/238 age (meas.)'] = dage
                df['∆206/238 age (tot.)'] = dagetot
            elif Pb_Th_std_crct_selector == 'Common Pb + Th Disequil.':
                if UTh_std_norm == 'Off':
                    df['206Pb/238UPbThc'] = df['206Pb/238Upbc'] - \
                        (lambda_238/lambda_230*((ThU_zrn/ThU_magma)-1))
                    df['206Pb/238Uc_age'] = (
                        np.log(df['206Pb/238UPbThc'] + 1) / lambda_238)
                    df['206Pb/238U_correctedage'] = df['206Pb/238Uc_age']*frac_factor
                    # propagate errors
                    # error on fractionation factor. includes error from ID-TIMS and ICPMS
                    dfrac = np.abs(frac_factor)*(((tims_error/2)/tims_age)
                                                 ** 2 + (std_avgerr/avg_std_ratio_Thcrct)**2)**(1/2)
                    dage = np.abs(df['206Pb/238U_correctedage'])*((1/(1/df['238U/206Pb']))
                                                                  ** 2*df['206/238 Reg. err']**2 + (0.16/2/100)**2)**(1/2)
                    # error on age equation. error on decay constant includes 1.5* counting stats (Mattinson 1987)
                    # error on estimation for common lead using 207 method. Uses conservaitve estimates of 1.0 for 206/204 and 0.3 for 207/204 (Mattionson, 1987)
                    # UTh errors - only using errors from measured zircon here, as this will by and large be the largest error contribution
                    # should probably add error for U / Th in icpms glass analyses, though rock measurements will undoubtedly be incorrectly used by users of the
                    # program (in absence of other data) and so added error would overall be fairly misleading anyway
                    # errors on absolute concentrations from ID-TIMS are overall negligible and often not reported either.
                    # need to put in possibility of putting in errors on Th/U measurements
                    # total propagated error
                    dagetot = np.abs(df['206Pb/238U_correctedage'])*(dfrac**2 +
                                                                     (df['206/238 Reg. err']/(1/df['238U/206Pb']))**2 + (0.16/2/100)**2 +
                                                                     ((0.3/2)/df['SK 207Pb/204Pb'])**2 + ((1/2)/df['SK 206Pb/204Pb'])**2 +
                                                                     (df['SE 207/206']/df['207Pb/206Pb'])**2)**(1/2)

                    df['∆206/238 age (meas.)'] = dage
                    df['∆206/238 age (tot.)'] = dagetot
                elif UTh_std_norm == 'Calc U/Th from Std.':
                    I232_std = std['232Th'].mean()
                    I238_std = std['238U'].mean()
                    std232err = std['232Th_1SE'].mean()
                    std238err = std['238U_1SE'].mean()
                    ThUzrn_calc = (((1/UTh_std)/(1/UTh_std_m)) /
                                   (I232_std/I238_std)) * (df['232Th']/df['238U'])
                    df['206Pb/238UPbThc'] = df['206Pb/238Upbc'] - \
                        (lambda_238/lambda_230*((ThUzrn_calc/ThU_magma)-1))
                    df['206Pb/238Uc_age'] = (
                        np.log(df['206Pb/238UPbThc'] + 1) / lambda_238)
                    df['206Pb/238U_correctedage'] = df['206Pb/238Uc_age']*frac_factor
                    # propagate errors
                    # error on fractionation factor. includes error from ID-TIMS and ICPMS
                    dfrac = np.abs(frac_factor)*(((tims_error/2)/tims_age)
                                                 ** 2 + (std_avgerr/avg_std_ratio_Thcrct)**2)**(1/2)
                    # error on age equation. error on decay constant includes 1.5* counting stats (Mattinson 1987)
                    dage = np.abs(df['206Pb/238U_correctedage'])*((1/(1/df['238U/206Pb']))
                                                                  ** 2*df['206/238 Reg. err']**2 + (0.16/2/100)**2)**(1/2)
                    # error on estimation for common lead using 207 method. Uses conservaitve estimates of 1.0 for 206/204 and 0.3 for 207/204 (Mattionson, 1987)
                    # UTh errors
                    # total propagated error
                    dagetot = np.abs(df['206Pb/238U_correctedage'])*(dfrac**2 +
                                                                     (df['206/238 Reg. err']/(1/df['238U/206Pb']))**2 + (0.16/2/100)**2 +
                                                                     ((0.3/2)/df['SK 207Pb/204Pb'])**2 + ((1/2)/df['SK 206Pb/204Pb'])**2 +
                                                                     (df['232Th_1SE']/df['232Th'])**2 + (df['238U_1SE']/df['238U'] + (df['SE 207/206']/df['207Pb/206Pb']))**2 +
                                                                     (std232err/I232_std)**2 + (std238err/I238_std)**2)**(1/2)
                    df['∆206/238 age (meas.)'] = dage
                    df['∆206/238 age (tot.)'] = dagetot

            return df

        elif ellipse_mode_selector == 'Ellipses':
            frac_factor, tims_age, tims_error, std_avg, avg_std_age_Thcrct, std_avgerr, UTh_std, UTh_std_m = calc_fncs.get_standard_fracfctr(std, std_txt, Pb_Th_std_crct_selector, regression_selector,
                                                                                                                                             ellipse_mode_selector, power)

            points, concordia_238_206, pts_pb_r, n, a = calc_fncs.get_projections(
                df, ellipse_mode_selector, power)

            # set up empty arrays to be filled for common lead corrections
            common_filter = []

            pb_m = df['207Pb/206Pb']  # measured 207/206
            common = df['SK 207Pb/206Pb']

            for i in common:
                if i <= 0:
                    common_filter.append(0)
                else:
                    common_filter.append(i)

            # calculate fraction of common Pb
            f_ = (pb_m - pts_pb_r) / (common - pts_pb_r)
            # set up array to set f = 0 if point lies on or below Concordia (i.e., no common Pb present)
            f = []

            for k, j in zip(common_filter, f_):
                if k <= 0:
                    f.append(0)
                elif j < 0:
                    f.append(0)
                else:
                    f.append(j)

            df['207Pb/206Pbr'] = pts_pb_r
            df['f'] = f

            df['counts_pb206r'] = df['206Pb'] * (1-df['f'])
            df['206Pb/238Upbc_numerical'] = 1 / \
                df['238U/206Pb']-(1/df['238U/206Pb']*f)
            df['206Pb/238Upbc'] = 1/concordia_238_206

            # will need input for which standard to use now! standard will need to be calculated by a separate function and called/input here

            if Pb_Th_std_crct_selector == 'Common Pb':
                df['206Pb/238Uc_age'] = np.log(
                    df['206Pb/238Upbc'] + 1) / lambda_238
                df['206Pb/238U_correctedage'] = df['206Pb/238Uc_age']*frac_factor

                df['206Pb/238U Corrected'] = np.exp(
                    df['206Pb/238U_correctedage']*lambda_238) - 1

                ellparams_i = calc_fncs.get_ellipse(df, power)
                ellparams_i = pd.DataFrame([ellparams_i], columns=[
                                           'Ellipse Center', 'Ell. Width', 'Ell. Height', 'Ell. Rotation'])

            elif Pb_Th_std_crct_selector == 'Common Pb + Th Disequil.':
                if UTh_std_norm == 'Off':
                    df['206Pb/238UPbThc'] = df['206Pb/238Upbc'] - \
                        (lambda_238/lambda_230*((ThU_zrn/ThU_magma)-1))
                    df['206Pb/238Uc_age'] = (
                        np.log(df['206Pb/238UPbThc'] + 1) / lambda_238)
                    df['206Pb/238U_correctedage'] = df['206Pb/238Uc_age']*frac_factor

                    df['206Pb/238U Corrected'] = np.exp(
                        df['206Pb/238U_correctedage']*lambda_238) - 1

                    ellparams_i = calc_fncs.get_ellipse(df, power)
                    ellparams_i = pd.DataFrame([ellparams_i], columns=[
                                               'Ellipse Center', 'Ell. Width', 'Ell. Height', 'Ell. Rotation'])

                elif UTh_std_norm == 'Calc U/Th from Std.':
                    I232_std = std['232Th'].mean()
                    I238_std = std['238U'].mean()
                    std232err = std['232Th_1SE'].mean()
                    std238err = std['238U_1SE'].mean()
                    ThUzrn_calc = (((1/UTh_std)/(1/UTh_std_m)) /
                                   (I232_std/I238_std)) * (df['232Th']/df['238U'])
                    df['206Pb/238UPbThc'] = df['206Pb/238Upbc'] - \
                        (lambda_238/lambda_230*((ThUzrn_calc/ThU_magma)-1))
                    df['206Pb/238Uc_age'] = (
                        np.log(df['206Pb/238UPbThc'] + 1) / lambda_238)
                    df['206Pb/238U_correctedage'] = df['206Pb/238Uc_age']*frac_factor

                    df['206Pb/238U Corrected'] = np.exp(
                        df['206Pb/238U_correctedage']*lambda_238) - 1

                    ellparams_i = calc_fncs.get_ellipse(df, power)
                    ellparams_i = pd.DataFrame([ellparams_i], columns=[
                                               'Ellipse Center', 'Ell. Width', 'Ell. Height', 'Ell. Rotation'])
        return ellparams_i

    def get_standard_fracfctr(std, std_txt, Pb_Th_std_crct_selector, regression_selector, ellipse_mode_selector, power):
        """
        """

        accepted_ages = {
            'Temora': 416780000,
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
            'Temora': 330000,
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
            avg_std_age, avg_std_age_Thcrct, avg_std_ratio, avg_std_ratio_Thcrct, avg_reg_err, UTh_std, UTh_std_m = calc_fncs.correct_standard_ages(
                std, std_txt, Pb_Th_std_crct_selector, regression_selector, ellipse_mode_selector, power)
            frac_factor = accepted_ages.get(std_txt)/avg_std_age
            tims_age = accepted_ages.get(std_txt)
            tims_error = TIMS_errors.get(std_txt)
        elif Pb_Th_std_crct_selector == 'Common Pb + Th Disequil.':
            avg_std_age, avg_std_age_Thcrct, avg_std_ratio, avg_std_ratio_Thcrct, avg_reg_err, UTh_std, UTh_std_m = calc_fncs.correct_standard_ages(
                std, std_txt, Pb_Th_std_crct_selector, regression_selector, ellipse_mode_selector, power)
            frac_factor = accepted_ages.get(std_txt)/avg_std_age_Thcrct
            tims_age = accepted_ages.get(std_txt)
            tims_error = TIMS_errors.get(std_txt)

        return frac_factor, tims_age, tims_error, avg_std_age, avg_std_age_Thcrct, avg_std_ratio, avg_std_ratio_Thcrct, avg_reg_err, UTh_std, UTh_std_m


# %%
class finalize_ages(param.Parameterized):

    """ class that parameterizes inputs and sends them to the above functions to be rendered in a GUI"""
    file_path = param.String(default='Insert File Path')
    input_data = param.DataFrame(precedence=-1)
    output_data = param.DataFrame(precedence=-1)
    drift_data = param.DataFrame(precedence=-1)
    file_path_ellipse = param.String(default='Insert Ellipse File Path')
    input_data_ellipse = param.DataFrame(precedence=-1)
    output_data_ellipse = param.DataFrame(precedence=-1)
    regression_selector = param.Selector(
        objects=['1st Order', '2nd Order', 'Exp.'])
    ellipse_mode_selector = param.Selector(default='Point Estimates', objects=[
                                           'Point Estimates', 'Ellipses'])
    drift_analyte_dropdown = param.Selector(objects=[])

    x_axis_TW_min = param.Number(default=0.5)
    x_axis_TW_max = param.Number(default=25)
    y_axis_TW = param.Range(
        default=(0.0, 0.5), bounds=(0.0, 1.0))  # for y-axis
    text_sample_selector = param.String(default='Input Sample ID')
    text_standard_selector = param.String(default='Input Standard ID')
    text_drift_selector = param.String(default='Input Drift ID')
    common_206204_input = param.Number()
    common_207204_input = param.Number()
    ThU_zrn_input = param.Number()
    ThU_magma_input = param.Number()
    UTh_std_norm = param.Selector(default='Off', objects=[
                                  'Calc U/Th from Std.', 'Off'])
    Pb_Th_std_crct_selector = param.Selector(
        objects=['Common Pb', 'Common Pb + Th Disequil.'])
    mass_bias_nist_selector = param.Selector(objects=['No Pb-Pb Mass Bias','Pb-Pb Mass Bias'])
    power = param.Number(default=0.05)

    update_output_button = param.Action(
        lambda x: x.add_output_data(), label='Approve Data')
    export_data_button = param.Action(lambda x: x.export_data(), label='DDDT!')
    export_TWplot_button = param.Action(
        lambda x: x.export_plot(), label='Save Plot')
    label_toggle = param.ListSelector(default=['Concordia'], objects=[
                                      'Concordia', 'Box + Whisker'])

    def __init__(self, **params):
        super().__init__(**params)
        self.input_data_widget = pn.Param(self.param.input_data),
        self.output_data_widget = pn.Param(self.param.output_data),
        self.input_data_ellipse_widget = pn.Param(
            self.param.input_data_ellipse)
        self.output_data_ellipse_widget = pn.Param(
            self.param.output_data_ellipse)
        self.widgets = pn.Param(self, parameters=['text_standard_selector', 'label_toggle', 'regression_selector', 'ellipse_mode_selector',
                                                  'update_output_button', 'export_data_button', 'export_TWplot_button',
                                                  'y_axis_TW',
                                                  'x_axis_TW_min', 'x_axis_TW_max',
                                                  'common_206204_input', 'common_207204_input', 'ThU_zrn_input', 'ThU_magma_input', 'power'])

    @pn.depends('file_path', 'regression_selector', 'ellipse_mode_selector', 'text_drift_selector','mass_bias_nist_selector', watch=True)
    def _uploadfile(self):
        if self.ellipse_mode_selector == 'Point Estimates':
            if self.file_path != 'Insert File Path':
                df = pd.read_excel(self.file_path, sheet_name='Sheet1')
                self.input_data = df
                self.input_data = self.input_data.replace('bdl',0)
                if '206/238 1st Order' in self.input_data.columns and self.regression_selector == '1st Order':
                    self.input_data['238U/206Pb'] = 1 / self.input_data['206/238 1st Order']
                    self.input_data['206/238 Reg. err'] = self.input_data['SE 206/238 1st Order']
                elif '206/238 2nd Order' in self.input_data.columns and self.regression_selector == '2nd Order':
                    self.input_data['238U/206Pb'] = 1 / self.input_data['206/238 2nd Order']
                    self.input_data['206/238 Reg. err'] = self.input_data['SE 206/238 2nd Order']
                elif '206/238 Exp.' in self.input_data.columns and self.regression_selector == 'Exp.':
                    self.input_data['238U/206Pb'] = 1 / self.input_data['206/238 Exp.']
                    self.input_data['206/238 Reg. err'] = self.input_data['SE 206/238 Exp']
                else:
                    pass

            if self.text_drift_selector != 'Input Drift ID':
                chosen_drift_std = self.input_data[self.input_data['SampleLabel'].str.contains(self.text_drift_selector)]
                drift_analytes_df = pd.concat([chosen_drift_std.loc[:, 'measurementindex'], chosen_drift_std.loc[:, '238U/206Pb'], chosen_drift_std.loc[:, '206Pb/204Pb'],
                                               chosen_drift_std.loc[:, '207Pb/206Pb'], chosen_drift_std.loc[:, '238U/235U'], chosen_drift_std.loc[:, '202Hg':'238U']], axis=1)
                drift_analytes = list(drift_analytes_df.columns[1:])
                drift_slope, intercepts, f_parray, p_slopearray = calc_fncs.get_driftstats(drift_analytes_df, drift_analytes)
                self.drift_data = chosen_drift_std
                self.drift_data = self.drift_data.reset_index(drop=True)

                for i, j in zip(drift_slope, drift_analytes):
                    if i > 0:
                        self.input_data[j] = self.input_data[j] - self.input_data['measurementindex']*i
                    else:
                        pass
                self.param.drift_analyte_dropdown.objects = drift_analytes
                
            if self.mass_bias_nist_selector == 'Pb-Pb Mass Bias':
                pb206204bias,pb207206bias = calc_fncs.calc_pbpb_bias(self.input_data,self.text_drift_selector)
                self.input_data['206Pb/204Pb'] = self.input_data['206Pb/204Pb']*pb206204bias
                self.input_data['207Pb/206Pb'] = self.input_data['207Pb/206Pb']*pb207206bias

            self.input_data['206/238U_age_init'] = np.log((1/self.input_data['238U/206Pb']) + 1) / lambda_238
            self.output_data = pd.DataFrame([np.zeros(len(self.input_data.columns))], columns=list(self.input_data.columns))

        else:
            pass

    @pn.depends('drift_data', 'text_drift_selector', 'drift_analyte_dropdown', 'regression_selector', 'ellipse_mode_selector')
    def call_drift_plot(self):
        if self.ellipse_mode_selector == 'Point Estimates':
            if self.drift_analyte_dropdown is not None:
                chosen_drift_std = self.drift_data[self.drift_data['SampleLabel'].str.contains(self.text_drift_selector)]
                drift_analytes_df = pd.concat([chosen_drift_std.loc[:, 'measurementindex'], chosen_drift_std.loc[:, '238U/206Pb'], chosen_drift_std.loc[:, '206Pb/204Pb'],
                                               chosen_drift_std.loc[:, '207Pb/206Pb'], chosen_drift_std.loc[:, '238U/235U'], chosen_drift_std.loc[:, '202Hg':'238U']], axis=1)
                drift_analytes = list(drift_analytes_df.columns[1:])
                drift_corrections_slope, drift_corrections_intercept, f_parray, p_slopearray = calc_fncs.get_driftstats(drift_analytes_df, drift_analytes)
                drift_slope = drift_corrections_slope[drift_analytes.index(self.drift_analyte_dropdown)]
                drift_intercept = drift_corrections_intercept[drift_analytes.index(self.drift_analyte_dropdown)]
                drift_fp = f_parray[drift_analytes.index(self.drift_analyte_dropdown)]
                drift_pslope = p_slopearray[drift_analytes.index(self.drift_analyte_dropdown)]
                return calc_fncs.plot_drift(drift_analytes_df, self.drift_analyte_dropdown, drift_slope, drift_intercept, drift_fp, drift_pslope)

            else:
                pass
        else:
            pass

    @pn.depends('file_path_ellipse', 'ellipse_mode_selector', watch=True)
    def _uploadelllipsefile(self):
        if self.ellipse_mode_selector == 'Ellipses':
            if self.file_path_ellipse != 'Insert Ellipse File Path':
                df = pd.read_excel(self.file_path_ellipse, sheet_name='Sheet1')
                col_bdl_condn = df[(df['206Pb/238U'] == 'bdl') | (df['207Pb/206Pb']
                                                                  == 'bdl') | (df['207Pb/235U'] == 'bdl')].index
                # drop rows that have 'bdl' in the 206/238 or 207/206 columns
                df.drop(col_bdl_condn, inplace=True)
                df = df.reset_index(drop=True)
                # having 'bdl's in columns makes them not numeric and thus they need to be changed to be used in calculations
                df['206Pb/238U'] = pd.to_numeric(df['206Pb/238U'])
                df['207Pb/206Pb'] = pd.to_numeric(df['207Pb/206Pb'])
                df['207Pb/235U'] = pd.to_numeric(df['207Pb/235U'])
                self.input_data_ellipse = df
                self.input_data_ellipse['207Pb/206Pb'] = self.input_data_ellipse['207Pb/206Pb']
                self.input_data_ellipse['238U/206Pb'] = 1 / self.input_data_ellipse['206Pb/238U']
                self.input_data_ellipse['206/238U_age_init'] = np.log(self.input_data_ellipse['206Pb/238U'] + 1) / lambda_238
            self.output_data_ellipse = pd.DataFrame([np.zeros(4)], columns=['Ellipse Center', 'Ell. Width', 'Ell. Height', 'Ell. Rotation'])
        else:
            pass

    @pn.depends('input_data', 'input_data_ellipse', 'common_206204_input', 'common_207204_input', 'ellipse_mode_selector', 'text_sample_selector', watch=True)
    def _updateCommonPb(self):
        if self.text_sample_selector != 'Input Sample ID':
            if self.ellipse_mode_selector == 'Point Estimates':
                if self.common_206204_input != 0:
                    self.input_data['SK 206Pb/204Pb'] = self.common_206204_input
                elif self.common_206204_input == 0:
                    self.input_data['SK 206Pb/204Pb'] = 11.152 + 9.74*(np.exp(
                        lambda_238*3.7e9)-np.exp(lambda_238*self.input_data['206/238U_age_init']))
                if self.common_207204_input != 0:
                    self.input_data['SK 207Pb/204Pb'] = self.common_207204_input
                elif self.common_207204_input == 0:
                    self.input_data['SK 207Pb/204Pb'] = 12.998 + 9.74/137.82*(np.exp(
                        lambda_235*3.7e9)-np.exp(lambda_235*self.input_data['206/238U_age_init']))
                self.input_data['SK 207Pb/206Pb'] = self.input_data['SK 207Pb/204Pb'] / \
                    self.input_data['SK 206Pb/204Pb']
            elif self.ellipse_mode_selector == 'Ellipses':
                if self.common_206204_input != 0:
                    self.input_data_ellipse['SK 206Pb/204Pb'] = self.common_206204_input
                elif self.common_206204_input == 0:
                    self.input_data_ellipse['SK 206Pb/204Pb'] = 11.152 + 9.74*(np.exp(
                        lambda_238*3.7e9)-np.exp(lambda_238*self.input_data_ellipse['206/238U_age_init']))
                if self.common_207204_input != 0:
                    self.input_data_ellipse['SK 207Pb/204Pb'] = self.common_207204_input
                elif self.common_207204_input == 0:
                    self.input_data_ellipse['SK 207Pb/204Pb'] = 12.998 + 9.74/137.82*(np.exp(
                        lambda_235*3.7e9)-np.exp(lambda_235*self.input_data_ellipse['206/238U_age_init']))
                self.input_data_ellipse['SK 207Pb/206Pb'] = self.input_data_ellipse['SK 207Pb/204Pb'] / \
                    self.input_data_ellipse['SK 206Pb/204Pb']

    @pn.depends('input_data', 'input_data_ellipse', 'text_sample_selector', 'y_axis_TW', 'label_toggle',
                'x_axis_TW_min', 'x_axis_TW_max', 'ellipse_mode_selector', 'power')
    def call_TW(self):
        if self.ellipse_mode_selector == 'Point Estimates':
            if self.text_sample_selector != 'Input Sample ID':
                data_toplot = self.input_data[self.input_data['SampleLabel'].str.contains(
                    self.text_sample_selector)]
                return calc_fncs.plot_TW(data_toplot,
                                         self.x_axis_TW_min, self.x_axis_TW_max,
                                         self.y_axis_TW[0], self.y_axis_TW[1],
                                         self.label_toggle, self.ellipse_mode_selector, self.power)
        elif self.ellipse_mode_selector == 'Ellipses':
            if self.text_sample_selector != 'Input Sample ID':
                data_toplot = self.input_data_ellipse[self.input_data_ellipse['SampleLabel'].str.contains(
                    self.text_sample_selector)]
                return calc_fncs.plot_TW(data_toplot,
                                         self.x_axis_TW_min, self.x_axis_TW_max,
                                         self.y_axis_TW[0], self.y_axis_TW[1],
                                         self.label_toggle, self.ellipse_mode_selector, self.power)

    @pn.depends('input_data', 'text_sample_selector', 'text_standard_selector', 'label_toggle', 'ThU_zrn_input', 'ThU_magma_input', 'Pb_Th_std_crct_selector', 'regression_selector',
                'ellipse_mode_selector', 'power', 'UTh_std_norm')
    def call_boxplot(self):
        if self.ellipse_mode_selector == 'Point Estimates':
            if self.text_sample_selector != 'Input Sample ID':
                data_toplot = self.input_data[self.input_data['SampleLabel'].str.contains(
                    self.text_sample_selector)]
                chosen_std = self.input_data[self.input_data['SampleLabel'].str.contains(
                    self.text_standard_selector)]
                if self.text_sample_selector == 'Input Sample ID':
                    return 'Placeholder'
                else:
                    ages = calc_fncs.correct_sample_ages(data_toplot, chosen_std, self.text_standard_selector, self.ThU_zrn_input, self.ThU_magma_input, self.Pb_Th_std_crct_selector, self.regression_selector,
                                                         self.ellipse_mode_selector, self.power, self.UTh_std_norm)
                    return calc_fncs.plot_boxplot(ages['206Pb/238U_correctedage']/(1e6), ages['SampleLabel'], self.label_toggle, self.ellipse_mode_selector)
        else:
            pass

    def add_output_data(self, event=None):
        if self.ellipse_mode_selector == 'Point Estimates':
            data_to_update = self.input_data[self.input_data['SampleLabel'].str.contains(
                self.text_sample_selector)]
            chosen_std = self.input_data[self.input_data['SampleLabel'].str.contains(
                self.text_standard_selector)]
            ages = calc_fncs.correct_sample_ages(data_to_update, chosen_std, self.text_standard_selector, self.ThU_zrn_input, self.ThU_magma_input, self.Pb_Th_std_crct_selector, self.regression_selector,
                                                 self.ellipse_mode_selector, self.power, self.UTh_std_norm)
            if self.output_data is None:
                self.output_data = ages
            else:
                self.output_data = self.output_data.append(
                    ages, ignore_index=True)
        elif self.ellipse_mode_selector == 'Ellipses':
            chosen_std = self.input_data_ellipse[self.input_data_ellipse['SampleLabel'].str.contains(
                self.text_standard_selector)]
            data_to_update = self.input_data_ellipse[self.input_data_ellipse['SampleLabel'].str.contains(
                self.text_sample_selector)]
            for i in data_to_update.SampleLabel.unique():
                data = calc_fncs.correct_sample_ages(data_to_update[data_to_update['SampleLabel'] == i], chosen_std, self.text_standard_selector, self.ThU_zrn_input, self.ThU_magma_input, self.Pb_Th_std_crct_selector, self.regression_selector,
                                                     self.ellipse_mode_selector, self.power, self.UTh_std_norm)
                if self.output_data_ellipse is None:
                    self.output_data_ellipse = data
                else:
                    self.output_data_ellipse = self.output_data_ellipse.append(
                        data, ignore_index=True)

    @pn.depends('output_data', watch=True)
    def _update_data_widget(self):
        if self.output_data is not None:
            self.output_data_widget = self.output_data
            self.output_data_widget.height = 40
            self.output_data_widget.heightpolicy = 'Fixed'
            return pn.widgets.Tabulator(self.output_data_widget, width=800)

    @pn.depends('output_data')
    def export_data(self, event=None):
        if self.ellipse_mode_selector == 'Point Estimates':
            self.output_data.to_excel('output_lasertramZ_ages.xlsx')
        elif self.ellipse_mode_selector == 'Ellipses':
            self.output_data_ellipse.to_excel(
                'output_lasertramZ_ellipses.xlsx')

    def export_plot(self, event=None):
        if self.ellipse_mode_selector == 'Point Estimates':
            data_toplot = self.input_data[self.input_data['SampleLabel'].str.contains(
                self.text_sample_selector)]
            plot = calc_fncs.plot_TW(data_toplot,
                                     self.x_axis_TW_min, self.x_axis_TW_max,
                                     self.y_axis_TW[0], self.y_axis_TW[1],
                                     self.label_toggle, self.ellipse_mode_selector, self.power)
            plot.savefig('LaserTRAMZ_TW.pdf', format='pdf', dpi=250)
        elif self.ellipse_mode_selector == 'Ellipses':
            data_toplot = self.input_data_ellipse[self.input_data_ellipse['SampleLabel'].str.contains(
                self.text_sample_selector)]
            plot = calc_fncs.plot_TW(data_toplot,
                                     self.x_axis_TW_min, self.x_axis_TW_max,
                                     self.y_axis_TW[0], self.y_axis_TW[1],
                                     self.label_toggle, self.ellipse_mode_selector, self.power)
            plot.savefig('LaserTRAMZ_TW.pdf', format='pdf', dpi=250)


reduce_ages = finalize_ages(name='Reduce Ages')

# %%
pn.extension('tabulator', 'mathjax')

width_ratios = [10, 5]
grid_layout = pn.GridSpec(sizing_mode='scale_both')
# grid_layout = pn.GridSpec()

grid_layout[0, 0] = pn.Column(reduce_ages.call_TW)
grid_layout[2, 0] = pn.Row(reduce_ages._update_data_widget)

grid_layout[0, 2] = pn.Row(pn.WidgetBox(pn.Param(reduce_ages.param,
                                                 widgets={'label_toggle': pn.widgets.CheckBoxGroup,
                                                          'export_data_button': pn.widgets.Button(name='DDDT!', button_type='success'),
                                                          'regression_selector': pn.widgets.RadioButtonGroup,
                                                          'ellipse_mode_selector': pn.widgets.RadioButtonGroup,
                                                          'Pb_Th_std_crct_selector': pn.widgets.RadioButtonGroup,
                                                          'mass_bias_nist_selector': pn.widgets.RadioButtonGroup,
                                                          'UTh_std_norm': pn.widgets.RadioBoxGroup}
                                                 )
                                        )
                           )
grid_layout[0, 1] = pn.Column(reduce_ages.call_boxplot)
grid_layout[1, 0] = pn.Column(reduce_ages.call_drift_plot)

grid_layout.show()
