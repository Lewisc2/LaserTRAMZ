# LaserTRAMZ
## Citing and Related Documentation
#### See citation.cff file for citing information. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8097362.svg)](https://doi.org/10.5281/zenodo.8097362)
#### Publication coming 'soon'!

## Purpose
#### This program was written to handle time resolved U-Pb zircon data using LA-ICP-Quadrupole-MS and to provide a free, non-automated open-source program, with special interest given to 'young' (≤10Ma Zircons). This is effectively the sister to [LaserTRAM-DB](https://github.com/jlubbersgeo/laserTRAM-DB#readme) ([Lubbers et al., 202x](https://doi.org/10.5281/zenodo.7826697)) which was built for handling trace elements using the same instrumentation. The program was born from a need to measure young zircons on the quadrupole and is built to handle isotopic ratios for the analyst that wants to count as long as possible during analysis (i.e., no trace elements). A variety of reduction techniques are available to reduce isotopic ratios, as described below.

## Data Input Format
**Part One**
____________
#### Similar to LaserTRAM-DB, the [multifiler](https://github.com/jlubbersgeo/multifiler) repository  make the data prepared for LaserTRAMZ after deleting the 'Timestamp' header that is output in the LT_ready file. The input format should look like:
| SampleLabel        | Time           | 202Hg    | ... | 238U|
| ------------------ | -------------- | -------- | --- | --- |
| FishCanyon         | 12.74          | 100.0004 | ... | 0   |
| FishCanyon         | 158.31         | 400.0064 | ... | 0   |
| ...                | ...            | ...      | ... | ... |

**Part Two**
____________
#### Primary and validation standards have a required input format in order to get the proper age reduction. Adding standards can be done with ease (see contact information below). Currently, standards are from the PlasmAge Consortium:
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

#### The output excel file from part one, which is output into the folder where the script is, will have a column that needs to be deleted (column A) and a row that needs to be deleted (Row 2). Once this is done the file can be uploaded into Part Two and should look like:
| SampleLabel | t start | t end | ... | SE% 238/235 |
| ----------- | ------- | ------| --- | ----------- |
| FishCanyon  | 32      | 60    | ... | 1.04        |
| FishCanyon  | 32      | 57    | ... | 0.85        |
| ...         | ...     | ...   | ... | ...         |

## Part One: Analyte Reduction
#### Running the script for Part One will open a page in your default web browser with a grey column and some default selections. Copy + Paste the file path to your LT_ready file (see above) into the bottom blank string input. After some loading time, an empty data table will pop up and the dropdown titled 'Sample Subset' will be populated with the analyses in the file. Select one of these (it does not auto populate) and four plots as well as some regression statistics will pop up that looks like the following screenshot. Labelled plots correspond to the letters as follows:
#### A. Ablation intensites, B. Time resolved ratios., C. Confidence Ellipse derived from the eigen vectors of the covariance matrix, D. Regression statistics, E. Residuals for the 206Pb/238U fit, F. Output data table, G. GUI tools
#### More information is given about each of these below.
![Screenshot 2023-06-27 at 3 21 47 PM](https://github.com/Lewisc2/LaserTRAMZ/assets/65908927/4cf8a84d-87c0-4adb-8129-fcad73683de3)

#### A. Analytes can be toggled and put on a log scale using the buttons in the GUI tools. The y axis scale can be changed using hte slider titled 'Ablation plot ylim slider: x'.

#### B. Time resolved ratios can be toggled using the buttons in the GUI tools. Y-axes limit can be changed using the Ratio plot ylim min / max inputs. Background interval is chosen with the 'Background slider:'. Background begins at the lower number and ends at the higher number, and both are viewed as dashed vertical black lines on the ablation intensity plot (A) and the ratio plot. Ablation interval is the same except the interval is viewed as two solid black vertical lines. Estimating 206Pb/238U with a regression (Kosler et al., 2002) is based on the assumption that fractionation at the ablation site is nill at the *start* of the ablation. The regression can thus be projected back to the start by inputting a value into 'Ablation Start True' in the GUI tools. If attempting to reduce an interval that starts and ends deep downhole, clearly it is better to not project to the start of the analysis. Regression of 207/235 is coming in the near future and will simply employ the same paramters from the 6/38 regression assuming that fractionation is element dependent and not abundance dependent (which imo is not true based on my observations... I digress).
#### LaserTRAMZ also allows one to choose which regression to use. Current options are 1st order, 2nd order, and exponential (weighted regressions are coming soon). Fitting an exponential regression virtually always requires an initial 'guess'. Default parameters can be forced into the intial guess by toggling 'Force Exp. Params'. If the optimizer fails to provide a decent fit based on heuristic observation and/or the residuals (E), toggling 'Joggle Exp. Params' allows the user to manually input values into the paramters A, B, and C for the exponential fit equation y=A*exp(-b*t)+C, where t is the input value for 'Ablation Start True'. You may wish to take the log of this equation to prove to yourself that A is the y-intercept, b is the slope, and c is an additive term that moves the regression vertically (this should help you choose parameters). Regression statistics
#### Other relevant isotopic ratios (6/4, 7/35, 7/6, 38/35) are chosen using one of the three large depressable buttons in the GUI tools (Total Counts, Poisson, Means and Regression). Total counts integrates the counts and then takes the ratio. Means and regression is the standard community practice, where regression is for 206/238 and the means of the time resovled ratios and their standard deviation is taken. Poisson reduces 206 as a poisson statistic and 207 as a zero-inflated poisson statistics and is for the use of young zircons (Lewis et al., in prep). This also provides 204 (Hg corrected) as a zero-inflated poisson though we have found this to be mostly unreliable. Dwell times for integration should be input into the Arrayofdwelltimes input near the bottom of the GUI tools. These are in the same order as the analytes in the datatable. Do not delete the brackets or commas. Errors are output in 1SE.

#### C. Confidence ellipses are derived from every pass of the detector (Paton et al., 2010). LaserTRAMZ derives them directly from the covariance matrix. One can view the confidence ellipse for Tera-Wasserburg (TW) or Wetherhill Concordia (Weth.) by using the selection ribbon above the plot. Every pass (data point) contributing to the confidence ellipse is treated based on the selections for the reduction of isotopic ratios (see explanation above). The power associated with the confidence of the ellipse can be changed by changing the value in the 'Power' input. Note that choosing to keep the confidence ellipse will slow down the application in the browser. You can toggle this off by toggling the 'Generate Ellipse' button.

#### D. Regression statistics include R-squared values for each regression and the standard error (%) of the value for 'Ablation Start true (Ramsey and Schafer, 2013; Kosler et al., 2002).

#### E. Residual plots for assistance in selecting and analyzing regressions are plotted together.

#### F. Output datatable is displayed in the bottom. *The values in the table are editable*

#### G. Gui tools. Clicking the Approval Interval button will record the analysis. Clicking Export All Plots will export all plots currently displayed. Clicking DDDT! will output the datafile.
## Part Two: Age Reduction
#### Running the script for Part Two will open up a page in your default web browser. Before inputting a file, select 1) which 206Pb/238U regression you will be using and 2) whether you are getting point estimates (ages) or the entire confidence ellipse.  Copy and Paste the filepath from part one into the appropriate input (File path for point estimates, File path ellipse for confidence ellipses). Again, a datatable will pop up with a column of GUI tools. 
#### Type in any string into the 'Text sample selector' input to populate Tera-Wasserburg Concordia. Typing a standard name into the 'Text standard selector' input will populate the boxplot. The GUI should now contain: A. Concordia plot, B. Boxplot of sample ages (if using point estimates), C. Data Table, D. GUI Tools. All are discussed in more detail below.
![Screenshot 2023-06-27 at 3 25 33 PM](https://github.com/Lewisc2/LaserTRAMZ/assets/65908927/df0d65f9-74f6-4c84-aa5b-fb006c74c1e4)


#### A. Tera-Wasserburg Concordia is shown as the solid black line (Wetherhill option coming soon). The y and x axes limits may be adjusted using hte X axis min/max inputs and y-axis slider. Measured data are plotted as black diamonds and green dots along concordia are the concordant age without the common Pb component. Black lines passing through these points are the projection from common Pb (Stacey and Kramers, 1975) through the measured data point and onto concordia. LaserTRAMZ outputs the concordant age.
#### Common Pb composition may be manually input if measured externally (e.g., feldspars) into the Common 206204 / 207204 inputs. These must be input before inputting the datafile. Default of the program is to correct only for Common Pb and not Th-disequilibrium.
#### Th-disequilbrium correction (Schaer, 1984) may be chosen by clicking the depressable Common Pb+Th Disequil. button. [Th/U] (ug/g) in the coexisting melt and [Th/U] (ug/g) in the zircon may be input in the ThU magma input and ThU zrn input, respectively. Should the user be confident in the measured Th/U in the standard and the value is mostly invariant (e.g., Plesovice) Th/U in the measured zircon can be reduced by clicking the Calc U/Th from Std. button and inputting a value into ThU magma input. The calculation follows from Kent et al. (2008).

#### B. Boxplot of samples ages with some basic statistics.

#### C. Output data table.

#### D. GUI tools. Click approve data to reduce the current set of samples. Save Plot button will output the currently displayed plots. DDDT! will output the dataframe.

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
#### We recommend running the program through [Anaconda](https://www.anaconda.com/download). You may also need to download [Git](https://github.com/git-guides/install-git).  After downloading, one may clone the repository by opening their terminal or Anconda prompt and running the following lines (one by one):
```
git clone https://github.com/Lewisc2/LaserTRAMZ.git
cd /path/to/LaserTRAMZ
pip install -r localrequirements.txt
```
#### where /path/to/LaserTRAMZ is the file path to the cloned repository. Once this is complete, the program can be run by opening Spyder from the Anaconda navigator and running the scripts, or,
```
cd /path/to/LaserTRAMZ
python LaserTRAMZ_Analyte_Redcution.py
```
#### for part one and 
```
cd /path/to/LaserTRAMZ
python LaserTRAMZ_Concordia.py
```
#### for part two

## Demos
* Coming soon!

## Feedback and Questions
Feecback and suggestions may be made by opening an [issue](https://github.com/Lewisc2/LaserTRAMZ/issues) or emailing Chuck Lewis (<lewisc2@oregonstate.edu>). I'm happy to answer any questions that may come up.

## To-Do List
* Wetherhill Concordia
* Auto Session wide Drift Correction
* Demo Videos
* Get this GD paper published
* 207/235 Regression
* Weighted Regressions

## References
