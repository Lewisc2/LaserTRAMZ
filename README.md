# LaserTRAMZ
## PROGRAM UNDER CONSTRUCTION
Old files on another branch. New files are fully functional and run faster than old program. More features in terms of plot navigation and data reduction were also added. Keep an eye out here for more concise updates.

## Purpose
LaserTRAMZ was written to handle time resolved U-Pb zircon analyses from LA-ICP-Quadrupole-MS and to provide a free, non-automated open-source program, with special interest given to 'young' (≤10Ma Zircons). This is effectively the sister to [LaserTRAM-DB](https://github.com/jlubbersgeo/laserTRAM-DB#readme) ([Lubbers et al., 202x](https://doi.org/10.5281/zenodo.7826697)) which was built for handling trace elements using the same instrumentation. The program was born from a need to measure young zircons on the quadrupole and is built to handle isotopic ratios for the analyst that wants to count as long as possible during analysis (i.e., no trace elements). A variety of reduction techniques are available to reduce isotopic ratios, as described below.

## Data Input Format
### Part 1: Analyte Reduction

Similar to LaserTRAM-DB, data are in the following format:

| SampleLabel        | Time           | 202Hg    | ... | 238U|
| ------------------ | -------------- | -------- | --- | --- |
| FishCanyon         | 12.74          | 100.0004 | ... | 0   |
| FishCanyon         | 158.31         | 400.0064 | ... | 0   |
| ...                | ...            | ...      | ... | ... |

To make the data prepared for LaserTRAMZ, and to save a little bit of time, the [multifiler](https://github.com/jlubbersgeo/multifiler) tool can be used. Simply delete the 'Timestamp' header that is output in the LT_ready file. 

### Part 2: Standard reference materials

Primary and validation standards have a required input format in order to get the proper age reduction. Adding standards can be done with ease (see contact information below). Currently, standards are from the PlasmAge Consortium except for NIST which is used for a Pb-Pb mass bias correction and session-wide drift correction (see below):

| Standard          | Reference                 | Input Name Required |
| ----------------- | ------------------------- | ------------------- |
| Fish Canyon Tuff  | Schmitz and Bowring, 2001 | FishCanyon          |
| 94-35             | Klepeis et al., 1998      | 94-35               |
| Plešovice         | Sláma et al., 2008        | Plesovice           |
| Temora2           | Black et al., 2004        | Temora              |
| R33               | Black et al., 2004        | R33                 |
| 91500             | Wiedenbeck et al., 1995   | 91500               |
| FC-1              | Paces and Miller, 1993    | FC1                 |
| Oracle            | Bowring, unpublished      | Oracle              |
| Tan-BrA           | Pecha, unpublished        | Tan-Bra             |
| OG-1              | Stern et al., 2009        | OG1                 |
| NIST-614          | Woodhead and Hergt, 2007  | NIST614             |

The output excel file from part one, which is output into the folder where the script is, will have a column that needs to be deleted (column A) and a row that needs to be deleted (Row 2). Once this is done the file can be uploaded into Part Two and should look like:
| SampleLabel | t start | t end | ... | SE% 238/235 |
| ----------- | ------- | ------| --- | ----------- |
| FishCanyon  | 32      | 60    | ... | 1.04        |
| FishCanyon  | 32      | 57    | ... | 0.85        |
| ...         | ...     | ...   | ... | ...         |

## Part One: Analyte Reduction
Running the script for Part One will open a page in your default web browser with a grey column and some default selections. `Copy + Paste` the file path to your LT_ready file (see above) into the bottom blank string input. After some loading time, an empty data table will pop up and the dropdown titled `Sample Subset` will be populated with the analyses in the file. Select one of these (it does not auto populate) and four plots as well as some regression statistics will pop up that looks like the following screenshot. Labelled plots correspond to the letters as follows:

<ol type="A">
  <li>Ablation intensites</li>
  <li>Time resolved ratios</li>
  <li>Confidence Ellipse derived from the eigen vectors of the covariance matrix</li>
  <li>Regression statistics</li>
  <li>Residuals for the 206Pb/238U fit</li>
  <li>Output data table</li>
  <li>GUI tools</li>
</ol>

![Screenshot 2023-06-27 at 3 21 47 PM](https://github.com/Lewisc2/LaserTRAMZ/assets/65908927/3a56b040-a7f0-4a20-aad4-aba46f3c16a5)

### Detailed Descriptions


#### A. 
Analytes can be toggled and put on a log scale using the buttons in the GUI tools. The y axis scale can be changed using hte slider titled 'Ablation plot ylim slider: x'.

#### B. 
Time resolved ratios can be toggled using the buttons in the GUI tools. Y-axes limit can be changed using the Ratio plot ylim min / max inputs. Background interval is chosen with the 'Background slider:'. Background begins at the lower number and ends at the higher number, and both are viewed as dashed vertical black lines on the ablation intensity plot (A) and the ratio plot. Ablation interval is the same except the interval is viewed as two solid black vertical lines. Estimating 206Pb/238U with a regression (Kosler et al., 2002) is based on the assumption that fractionation at the ablation site is nill at the *start* of the ablation. The regression can thus be projected back to the start by inputting a value into 'Ablation Start True' in the GUI tools. If attempting to reduce an interval that starts and ends deep downhole, clearly it is better to not project to the start of the analysis. Regression of 207/235 is coming in the near future and will simply employ the same paramters from the 6/38 regression assuming that fractionation is element dependent and not abundance dependent (which imo is not true based on my observations... I digress).

LaserTRAMZ also allows one to choose which regression to use. Current options are 1st order, 2nd order, and exponential (weighted regressions are coming soon). Fitting an exponential regression virtually always requires an initial 'guess'. Default parameters can be forced into the intial guess by toggling 'Force Exp. Params'. If the optimizer fails to provide a decent fit based on heuristic observation and/or the residuals (E), toggling 'Joggle Exp. Params' allows the user to manually input values into the paramters A, B, and C for the exponential fit equation y=A*exp(-b*t)+C, where t is the input value for 'Ablation Start True'. You may wish to take the log of this equation to prove to yourself that A is the y-intercept, b is the slope, and c is an additive term that moves the regression vertically (this should help you choose parameters). Regression statistics

Other relevant isotopic ratios (6/4, 7/35, 7/6, 38/35) are chosen using one of the three large depressable buttons in the GUI tools (Total Counts, Poisson, Means and Regression). Total counts integrates the counts and then takes the ratio. Means and regression is the standard community practice, where regression is for 206/238 and the means of the time resovled ratios and their standard deviation is taken. Poisson reduces 206 as a poisson statistic and 207 as a zero-inflated poisson statistics and is for the use of young zircons (Lewis et al., in prep). This also provides 204 (Hg corrected) as a zero-inflated poisson though we have found this to be mostly unreliable. Dwell times for integration should be input into the Arrayofdwelltimes input near the bottom of the GUI tools. These are in the same order as the analytes in the datatable. Do not delete the brackets or commas. Errors are output in 1SE.

#### C. 
Confidence ellipses are derived from every pass of the detector (Paton et al., 2010). LaserTRAMZ derives them directly from the covariance matrix. One can view the confidence ellipse for Tera-Wasserburg (TW) or Wetherhill Concordia (Weth.) by using the selection ribbon above the plot. Every pass (data point) contributing to the confidence ellipse is treated based on the selections for the reduction of isotopic ratios (see explanation above). The power associated with the confidence of the ellipse can be changed by changing the value in the 'Power' input. Note that choosing to keep the confidence ellipse will slow down the application in the browser. You can toggle this off by toggling the 'Generate Ellipse' button.

#### D. 
Regression statistics include R-squared values for each regression and the standard error (%) of the value for 'Ablation Start true (Ramsey and Schafer, 2013; Kosler et al., 2002).

#### E. 
Residual plots for assistance in selecting and analyzing regressions are plotted together.

#### F. 
Output datatable is displayed in the bottom. *The values in the table are editable*

#### G. 
Gui tools. Clicking the Approval Interval button will record the analysis. Clicking Export All Plots will export all plots currently displayed. Clicking DDDT! will output the datafile.

## Part Two: Age Reduction
Running the script for Part Two will open up a page in your default web browser. Before inputting a file, select 1) which 206Pb/238U regression you will be using and 2) whether you are getting point estimates (ages) or the entire confidence ellipse.  Copy and Paste the filepath from part one into the appropriate input (File path for point estimates, File path ellipse for confidence ellipses). Again, a datatable will pop up with a column of GUI tools. 

ype in any string into the 'Text sample selector' input to populate Tera-Wasserburg Concordia. Typing a standard name into the 'Text standard selector' input will populate the boxplot. The GUI should now contain:

<ol type="A">
  <li>Concordia plot</li>
  <li>Boxplot of sample ages (if using point estimates)</li>
  <li>Data Table</li>
  <li>GUI Tools</li>
  <li>Drift Correction</li>

</ol>

![Screenshot 2023-07-06 at 4 34 29 PM](https://github.com/Lewisc2/LaserTRAMZ/assets/65908927/b0009092-fd8f-4ada-8c69-7b6f1a5ee3b2)

### Detailed Descriptions

#### A. 
Tera-Wasserburg Concordia is shown as the solid black line (Wetherhill option coming soon). The y and x axes limits may be adjusted using hte X axis min/max inputs and y-axis slider. Measured data are plotted as black diamonds and green dots along concordia are the concordant age without the common Pb component. Black lines passing through these points are the projection from common Pb (Stacey and Kramers, 1975) through the measured data point and onto concordia. LaserTRAMZ outputs the concordant age.

Common Pb composition may be manually input if measured externally (e.g., feldspars) into the Common 206204 / 207204 inputs. These must be input before inputting the datafile. Default of the program is to correct only for Common Pb and not Th-disequilibrium.

Th-disequilbrium correction (Schaer, 1984) may be chosen by clicking the depressable Common Pb+Th Disequil. button. [Th/U] (ug/g) in the coexisting melt and [Th/U] (ug/g) in the zircon may be input in the ThU magma input and ThU zrn input, respectively. Should the user be confident in the measured Th/U in the standard and the value is mostly invariant (e.g., Plesovice) Th/U in the measured zircon can be reduced by clicking the Calc U/Th from Std. button and inputting a value into ThU magma input. The calculation follows from Kent et al. (2008).

#### B. 
Boxplot of samples ages with some basic statistics.

#### C. 
Output data table.

#### D. 
GUI tools. Click approve data to reduce the current set of samples. Save Plot button will output the currently displayed plots. DDDT! will output the dataframe.

#### E. 
Drift Correction. Our practice is generally to use NIST-614 in order to monitor session wide drift, as this allows simple and highly reproducible drift corrections rather than monitoring drift in natural standards that frequently have documented heterogeneity. Drift corrections are only made for those analytes and ratios that have a statistically significant slope and f-test for the linear regression. Drift correction is carried out by shifting all values for a given analyte based on the slope of the regression. If the drift is significant, both measured and observed values will be visible in the drift correction plot. Type in NIST614 (see input format section above) into the appropriate string input in order to correct for drift. Note you may use any given standard, however.

We also use NIST-614 to correct for Pb-Pb isotope mass bias on 206/204 and 207/206. Select the depressable button in the GUI tools to correct for mass bias using NIST614. Corrected values will be in the output data.

## Final Output Information
* t start, t end, t project: The ablation start, end, and projected value.
* Analytes (e.g., 206Pb, 238U) and analyte errors (e.g., 206Pb_1SE, 238U_1SE): Mean background subtracted intensities and 1*Standard errors for each analyte.
* 206/238 1st, 2nd, exp. and SE/SE% 206/238 1st, 2nd, exp.: Value of the reduced 206/238 value determined at 'Ablation start true' in part one as well as the 1SE and 1SE%
* Other isotopic ratios and their SE: Values estimated according to the option in part one as well as 1 SE and 1 SE% for each
* R2 2st, 2nd, exp.: R-squared values for the 206/238 regression from part one
* 20iPb counts / SE/SE% 20iPb: Estimated counts and 1SE / 1SE% for lead isotope analytes when using the 'Poisson' option in part one. 
* 238U/206Pb / 206/238 Reg. Err: 1/value regression reduction chosen for part two and the error for the value.
* 206/238U_age_init: 'Guess' age to get the common Pb from common Pb. We avoid use of the 207/206 age (e.g., Pullen et al., 2018) as the program was written to deal with young zircons specifically.
* SK 20i/20i: Estiamted common Pb component from the Stacey-Kramers model. If inputting common Pb manually, this column should contain the value from that input.
* 207/206r: Radiogenic 207/206 composition after the common Pb correction
* f, counts_pb206r: fraction of common Pb estimated from the 207 method (e.g., Andersen et al., 2019) and the counts of 206Pb that are common Pb corrected
* 206Pb/238Upbc_numerical and 206Pb/238Upbc: Calculation of common Pb corrected 206Pb/238U and the concordant 206Pb/238U. These should be the same value and were put here mostly so I could double check my work when writing the code and I haven't gotten around to taking one out yet.
* 206Pb/238UPbThc: Pb and Th corrected 206/238
* 206Pb/238Uc_age: Common Pb corrected 206Pb/238U age
* 206Pb/238U_corrected age: If using Th correction, this will be the Pb and Th corrected age. Otherwise it will only be common Pb corrected
* ∆206/238 age (meas.): (mostly) just the analytical error on the age (only assumes 206/238 error and the decay constant error) *do not use this error when reporting ages* 
* ∆206/238 age (tot.): Fully propagated error. See Lewis et al. (in prep)
  
## Installation and Use
We recommend running the program through [Anaconda](https://www.anaconda.com/download). You may also need to download [Git](https://github.com/git-guides/install-git).  After downloading, one may clone the repository by opening their terminal or Anconda prompt and running the following lines (one by one). It is best to create a virtual environment, which is included in the code block below:

```
git clone https://github.com/Lewisc2/LaserTRAMZ.git
cd /path/to/LaserTRAMZ
conda create -n LaserTRAMZ python==3.9.17
conda activate LaserTRAMZ
pip install -r localrequirements.txt
```

where /path/to/LaserTRAMZ is the file path to the cloned repository. You will need to accept the install by typing y then pressing enter when prompted. Once this is complete, the program can be run by opening your IDE of choice and running the scripts or

```
cd /path/to/LaserTRAMZ
conda activate LaserTRAMZ
python LaserTRAMZ_Analyte_Redcution.py
```

for part one and 

```
cd /path/to/LaserTRAMZ
conda activate LaserTRAMZ
python LaserTRAMZ_Concordia.py
```
for part two.

To shut down the virtual environemnt, run the following:
```
conda deactivate LaserTRAMZ
```

## Demos
(https://www.youtube.com/watch?v=683U5F1hsdM&list=PLs1w2r4kCIlWRb5ARSq0PzIfl4XBUQuaG)

## Feedback and Questions
Feecback and suggestions may be made by opening an [issue](https://github.com/Lewisc2/LaserTRAMZ/issues) or emailing Chuck Lewis (lewisc2 _at_ oregonstate.edu). I'm happy to answer any questions that may come up.

## Citing and Related Documentation
If you use LaserTRAMZ in your work, please citte it! See citation.cff file for citing information. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8097362.svg)](https://doi.org/10.5281/zenodo.8097362)

- Publication coming `soon`!
  
## Publications utilizing LaserTRAMZ

- Lewis et al., (_in prep_)

## To-Do List
* ~~Wetherhill Concordia~~
* ~~Auto Session wide Drift Correction~~
* ~~Pb-Pb Mass Bias Correction (NIST614)~~
* Fit Pb loss dates on Wetherhill
* Eruption age estimates
* Demo Videos
* Get this GD paper published
* ~~207/235 Regression~~
* Weighted Regressions



## References
