# LaserTRAMZ
____________
## Purpose
LaserTRAMZ was written to provide a free, non-automated open-source program to handle time resolved U-Pb zircon analyses gathered via LA-Q-ICP-MS. This is effectively the sister to [LaserTRAM-DB](https://github.com/jlubbersgeo/laserTRAM-DB#readme) ([Lubbers et al., 202x](https://doi.org/10.5281/zenodo.7826697)) which was built for handling trace elements using the same instrumentation.
____________
## Reference flowchart for repeat users that don't necessarily need to read the details below
[LaserTRAMZ_Reduction_Flowchart.pdf](https://github.com/user-attachments/files/19507746/LaserTRAMZ_Reduction_Flowchart.pdf)
____________
## Data Input Format
### Part 1: Analyte Reduction

Data need to be uploaded in the following format:

| SampleLabel        | Time           | 202Hg    | ... | 238U|
| ------------------ | -------------- | -------- | --- | --- |
| FishCanyon         | 12.74          | 100.0004 | ... | 0   |
| FishCanyon         | 158.31         | 400.0064 | ... | 0   |
| ...                | ...            | ...      | ... | ... |

To make the data prepared for LaserTRAMZ the [multifiler](https://github.com/jlubbersgeo/multifiler) tool can be used. Simply delete the 'Timestamp' header that is output in the LT_ready file before uploading. Upload by copy+pasting the file path into the appropriate box (see detail below).

### Part 2: Age Calculation

Part 1 will output an excel file into the LaserTRAMZ folder which is output into the folder. :
| SampleLabel | t start | t end | ... | TW rho      |
| ----------- | ------- | ------| --- | ----------- |
| FishCanyon  | 32      | 60    | ... | 0.05        |
| FishCanyon  | 32      | 57    | ... | 0.85        |
| ...         | ...     | ...   | ... | ...         |

### Standard reference materials

Primary and validation standards have a required naming scheme in order to be read by LaserTRAMZ properly. Adding in-house standards can be done with ease if requested. Standards currently supported are those available the PlasmAge Consortium as well as NIST-61x standard glasses:

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
| NIST-614          | Woodhead and Hergt, 2007 / Duffin et al., 2015  | NIST-614             |
| NIST-612          | Woodhead and Hergt, 2007 / Duffin et al., 2015  | NIST-614             |
| NIST-610          | Woodhead and Hergt, 2007 / Duffin et al., 2015  | NIST-614             |
____________
## Part One: Analyte Reduction
Running the script for Part One will open a page in your default web browser with a grey column and some default selections. `Copy + Paste` the file path to your LT_ready file into the bottom string input. After some loading time, an empty data table will pop up and the dropdown titled `Sample Subset` will be populated with the analyses in the file. Select one of these (it does not auto populate) and the program will load the data onto the screen:

<img width="1061" alt="Screenshot 2025-03-20 at 3 34 25 PM" src="https://github.com/user-attachments/assets/19848f12-c55c-4b3f-8a5d-05fb1987bd07" />

The top graph shows the intensities for each analyte. The bottom left graph shows the time resolved ratios selected in the left banner. The bottom right graph will always only show the 207Pb/206Pb ratio. Note that the plots are interactable (e.g., panning, zooming, downloading, etc). Buttons at the top of the plot toggle the interaction options.

### Descriptions of GUI visualization and tools
Background and ablation sliders are used to choose the respective intervals. These are plotted as dashed and solid lines, respectively. The Ablation Start True input allows the user to project the Pb/U ratios back to the start of ablation to deal wtih LIEF (Kosler et al., 2002). Experience leads me to recommended that you bind this to the ablation start by selecting the ablation start true button - the uncertainty will be smaller and the ratios will be more consistent.
Sliders for the Y-axis limits are also shown. The default limits of the graphs can be changed by highlighting the numbers above the sliders and typing in a number.
Buttons next to ratios and isotopes are for toggling these on the graphs. Ratios only change the bottom left graph. You may also put the intensities plot on a log scale.
The Total Counts, Poisson, and Means and Regression buttons allow the user to change how the isotope ratios are handled.
* Total Counts: Individual signals are integrated and a ratio of these 'total counts' is taken. See below for changing integration/dwell times.
* Poisson: 207/206 ratios are dealt with according to the manuscript currently under review here: https://eartharxiv.org/repository/dashboard/8699/. Pb/U ratios are dealt with by regression (c.a. Paton et al., 2010).
* Means and Regression: Pb/U ratios are dealt with by regression (c.a. Paton et al., 2010). Other ratios use the mean and SE of time-resolved ratios.

Power: Changes the alpha value for the confidence ellipsoid. Leave at 0.05 for the typical 95% confidence ellipsoid.
Arrayofdwelltimes: These are the dwell times for each analyte from low mass to high mass. By default everything is 10 ms. To change them simply type the dwell times into their respective position. Do NOT change any commas, brackets, etc.

### Evaluating an Interval and Recording Data
When you are ready to inspect the intervals you've chosen, click the Evaluate Interval button at the top left. A pop-up window will appear:

<img width="720" alt="Screenshot 2025-03-20 at 3 34 34 PM" src="https://github.com/user-attachments/assets/a4a01de3-c042-4339-9122-e1e394c6e639" />

The fitted regressions for the Pb/U ratios are shown in the top left graph. Residuals for the regressions are shown on the right. Confidence ellipsoids are shown in the bottom left. Regression statistics are shown in the lower right. If you are satisfied with the ablation interval, click Accept Interval in the top left. If not, close the pop-up window and adjust as desired, then re-evaluate.

### Exporting Data
Click the DDDT! button to export.
____________
## Part Two: Age Reduction
### Uploading and Setting reduction options
Running the script for Part Two will open up a page in your default web browser. Copy and Paste the filepath from part one into the file path input at the top of the sidebar. A window will pop up that looks like the following:

<img width="464" alt="Screenshot 2025-03-20 at 3 35 40 PM" src="https://github.com/user-attachments/assets/1258ef27-4e5c-4486-bd2e-d15b885377da" />

#### Input Options are in the order that follows:
* Primary Standard: Start typing the standard you wish to use, which should have a name the same as the table above. A dropdown will pop-up and you can either click the standard name or finish typing.
* Secondary Standard Excess Errors: Choose the standard you wish to use to evaluate excess variance according to the normalized ratios, c.a. Horstwood et al. (2016)
* List of buttons with sample names: Check those you want to use as validation standards. You need at least one for the program to run
* By Age, etc: Pb-Pb mass bias correction standard. If by Age, uses primary ZRM. Currently uses exponential correction of Woodhead (2002).
* Pb-Pb ratios list: Which measured Pb-Pb ratio to use for correcting Pb-Pb mass bias
* Primary 238U/235U, ...: Which standard to use for correcting the 238/235 mass bias. Primary 238U/235U uses primary ZRM and assumes 137.818 (Hiess et al., 2012)
* 238U/235U: Which ratio to use for correcting 38/35 bias. Currently only allows 38/35
* By Age, None: Drift correct according to reduced standard ages or not at all.
* Nearest Number: Nearest number of standards to use for correcting data. I.e., this is the sliding window of Gehrels et al. (2008)
* Decay Series Corrections: Whether to correct just for Common Pb or to correct for Common Pb and Th Diseuqilibrium
* Estimate Zircon [U/Th] By: How to estimate the [U/Th] concentration in zircon. For primary and selected secondary (same one as that used for excess variance), simply gets a factor according to the 238U and 232Th signals. Selecting Fitted regression standard will find all standards that you've selected (primary + those in the list of buttons) that have accepted concentrations and a simple linear regression will be used to make a calibration curve. If you have estimated the [U/Th] by some other way (e.g., split stream) you MUST manually enter these into the output from part one yourself. The column headers MUST be the following: [U] µg/g, [Th] µg/g, [U/Th].
* Calculate D [Th/U] From: First option uses inputs below to estimate D [Th/U] for the Th-disequilibrium correction (e.g., if you want to assume a constant value). Second option uses the estimated [U/Th] and the input melt [U/Th] in the input below.
* Th/U Zircon: Th/U ratio to use in zircon for Th disequilibrium correction
* Th/U Magma: Th/U ratio to use for host melt for Th disequilbrium correction
* Calculate Excess Variance From: First option increases 6/38 and 7/6 uncertainty equally until selected secondary ages have an MSWD of one. Second option increases 6/38 and 7/6 uncertainty until the normalized ratios of the preferred secondary standard (selected above) have an MSWD of one. Third option increases the uncertainty of the raw primary standard ratios until an MSWD of one is achieved.
* Exported Data Format: Simple is your best bet and has all the data you need to report. Annotated has additional output columns to help users assess sources of uncertainties. Full output is admittedly a bit of a mess and is mostly there for program development.
* Accept Reduction Parameters: Click this value to accept all your inputs

### Descriptions of GUI visualization and tools
Type a sample name into the Text sample selector box. The input only needs to be part of a sample name (e.g., typing 91500 will pull up all analyses that have "91500" in them). You can change the input [Th/U] concentrations of zircon and melt using the boxes on the left. After doing so, you should get a screen that looks like the following:

<img width="939" alt="Screenshot 2025-03-20 at 3 37 30 PM" src="https://github.com/user-attachments/assets/be6f18fc-b11b-40b1-95e1-868ff9ec5e50" />


* By default the common Pb correction uses the Stacey-Kramers model. If you have an estimate for the common 207Pb/206Pb ratio and its error you may input these in the appropriate boxes on the left.
Tera-Waserburg and Wetherhill Concordia can be visualized by selecting the tabs at the top. Limits to the Concordia plots can be changed using the slider or highlighting the numbers above the sliders and typing a number.
* The Concordia plot can be toggled from Tera-Waserburg to Wetherhill. By default the Stacey-Kramers model is used for the common Pb correction and the resultant ages are those projected onto Concordia from Common Pb through the data point. Common Pb ratios can be changed in the input on the left, as can the Th-disequilibrium parameters.
* The drift plot currently has some bugs and is mostly useful to load data first without drift correcting in order to assess how severe the drift is.
* The excess variance plot shows the ages (for secondary standard options) and ratios of standard data used to calculate excess variance. Black bars are 2s. Red is the additional excess variance.

#### The following standards have published concentrations for U and Th. Note Plesovice is not included due to the known heterogeneity (Slama et al., 2008). Choosing a standard(s) that does not have published concentrations and selecting an option that includes that (those) standard(s) to estimate concentrations will result in erroneous concentration estimations for unknowns (values set to 1e-7 in the script to avoid divide by zero errors).
* Temora
* Fish Canyon
* R33
* 91500
* FC1
* Oracle
* OG1
____________
## Final Output Information
* t start, t end, t project / bstart, bend: The ablation & background start, end, and projected value.
* Analytes (e.g., 206Pb, 238U) and analyte errors (e.g., 206Pb_1SE, 238U_1SE): Mean background subtracted intensities and 1*Standard errors for each analyte.
* Isotope ratios and their 1SE% are output. Pb/U ratios have their reduction method (i.e., regression type) included in the headers.
* Weth / TW C, Wid1, Wid2, rho: center, axes widths, and degree of rotation for confidence ellipsoids on Wetherhill and Tera-Wasserburg Concordia
* Isotope ratios with 'c' after them: Mass bias corrected isotope ratios. Also includes common Pb and Th-disequilibrium correction wehre present
* Pb/U ratio 'Age': Age in years
* Pb/U Age 1s (meas.): 1s analytical error on the age in years
* Pb/U Age 1s (tot): 1s full error on the age in year
____________
## Caveats
You need to have columns for 202Hg, 204-208Pb, 232Th, 235U, 238U for the program to run; largely due to the way I wrote the program initially and I haven't gotten around to changing this yet. If you need these columns to be input into the data there is a .ipynb file in the repository that can be used to input these columns provided you are familiar with python. I am also happy to do this for you if you reach out by email.
____________
## Installation and Use
We recommend running the program through [Anaconda](https://www.anaconda.com/download). You may also need to download [Git](https://github.com/git-guides/install-git).  After downloading, one may clone the repository by opening their terminal or Anconda prompt and running the following lines (one by one). It is best to create a virtual environment, which is included in the code block below:

```
git clone https://github.com/Lewisc2/LaserTRAMZ.git
cd /path/to/LaserTRAMZ
conda create -n LaserTRAMZ python==3.9.18
conda activate LaserTRAMZ
pip install -r localrequirements.txt
```

where /path/to/LaserTRAMZ is the file path to the cloned repository. You will need to accept the install by typing y then pressing enter when prompted. Once this is complete, the program can be run by opening your IDE of choice and running the scripts or

```
cd /path/to/LaserTRAMZ
conda activate LaserTRAMZ
python LaserTRAMZ_Quad_Analyte.py
```

for part one and 

```
cd /path/to/LaserTRAMZ
conda activate LaserTRAMZ
python LaserTRAMZ_Quad_Concordia.py
```
for part two.

To shut down the virtual environemnt, run the following:
```
conda deactivate LaserTRAMZ
```
____________
## Demos
(https://www.youtube.com/watch?v=683U5F1hsdM&list=PLs1w2r4kCIlWRb5ARSq0PzIfl4XBUQuaG)
____________
## Feedback and Questions
Feecback and suggestions may be made by opening an [issue](https://github.com/Lewisc2/LaserTRAMZ/issues) or emailing Chuck Lewis (lewisc2 _at_ oregonstate.edu). I'm happy to answer any questions that may come up.
____________
## Citing and Related Documentation
If you use LaserTRAMZ in your work, please citte it! See citation.cff file for citing information. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8097362.svg)](https://doi.org/10.5281/zenodo.8097362)

- Publication coming `soon`!
____________
## References
```
Black, L.P. et al., 2004, Improved 206Pb/238U microprobe geochronology by the monitoring of a trace-element-related matrix effect; SHRIMP, ID–TIMS, ELA–ICP–MS and oxygen isotope documentation for a series of zircon standards: Chemical Geology, v. 205, p. 115–140, doi:10.1016/j.chemgeo.2004.01.003.
Cheng, H., Edwards, R.L., Hoff, J., Gallup, C.D., Richards, D.A., and Asmerom, Y., 2000, The half-lives of uranium-234 and thorium-230: Chemical Geology, v. 169, p. 17–33, doi:10.1016/S0009-2541(99)00157-6.
Gehrels, G.E., Valencia, V.A., and Ruiz, J., 2008, Enhanced precision, accuracy, efficiency, and spatial resolution of U-Pb ages by laser ablation-multicollector-inductively coupled plasma-mass spectrometry: TECHNICAL BRIEF: Geochemistry, Geophysics, Geosystems, v. 9, p. n/a-n/a, doi:10.1029/2007GC001805.
Hiess, J., Condon, D.J., McLean, N., and Noble, S.R., 2012, 238U/235U Systematics in Terrestrial Uranium-Bearing Minerals: v. 335.
Horstwood, M.S.A. et al., 2016, Community‐Derived Standards for LA ‐ ICP ‐ MS U‐(Th‐)Pb Geochronology – Uncertainty Propagation, Age Interpretation and Data Reporting: Geostandards and Geoanalytical Research, v. 40, p. 311–332, doi:10.1111/j.1751-908X.2016.00379.x.
Jaffey, A.H., Flynn, K.F., Glendenin, L.E., Bentley, W.C., and Essling, A.M., 1971, Precision Measurement of Half-Lives and Specific Activities of U 235 and U 238: Physical Review C, v. 4, p. 1889–1906, doi:10.1103/PhysRevC.4.1889.
Johnston, S., Gehrels, G., Valencia, V., and Ruiz, J., 2009, Small-volume U–Pb zircon geochronology by laser ablation-multicollector-ICP-MS: Chemical Geology, v. 259, p. 218–229, doi:10.1016/j.chemgeo.2008.11.004.
Klepeis, K.A., Crawford, M.L., and Gehrels, G., 1998, Structural history of the crustal-scale Coast shear zone north of Portland Canal, southeast Alaska and British Columbia: Journal of Structural Geology, v. 20, p. 883–904, doi:10.1016/S0191-8141(98)00020-0.
Košler, J., Fonneland, H., Sylvester, P., Tubrett, M., and Pedersen, R.-B., 2002, U–Pb dating of detrital zircons for sediment provenance studies—a comparison of laser ablation ICPMS and SIMS techniques: Chemical Geology, v. 182, p. 605–618, doi:10.1016/S0009-2541(01)00341-2.
Košler, J., and Sylvester, P.J., 2003, Present Trends and the Future of Zircon in Geochronology: Laser Ablation ICPMS: Reviews in Mineralogy and Geochemistry, v. 53, p. 243–275, doi:https://doi.org/10.2113/0530243.
Mattinson, J.M., 1987, UPb ages of zircons: A basic examination of error propagation: Chemical Geology: Isotope Geoscience section, v. 66, p. 151–162, doi:10.1016/0168-9622(87)90037-6.
Paces, J.B., and Miller, J.D., 1993, Precise U‐Pb ages of Duluth Complex and related mafic intrusions, northeastern Minnesota: Geochronological insights to physical, petrogenetic, paleomagnetic, and tectonomagmatic processes associated with the 1.1 Ga Midcontinent Rift System: Journal of Geophysical Research: Solid Earth, v. 98, p. 13997–14013, doi:10.1029/93JB01159.
Pullen, A., Ibáñez-Mejia, M., Gehrels, G.E., Giesler, D., and Pecha, M., 2018, Optimization of a Laser Ablation-Single Collector-Inductively Coupled Plasma-Mass Spectrometer (Thermo Element 2) for Accurate, Precise, and Efficient Zircon U-Th-Pb Geochronology: Geochemistry, Geophysics, Geosystems, v. 19, p. 3689–3705, doi:10.1029/2018GC007889.
Schmitz, M.D., and Bowring, S.A., 2001, U-Pb zircon and titanite systematics of the Fish Canyon Tuff: an assessment of high-precision U-Pb geochronology and its application to young volcanic rocks: Geochimica et Cosmochimica Acta, v. 65, p. 2571–2587, doi:10.1016/S0016-7037(01)00616-0.
Sláma, J. et al., 2008, Plešovice zircon — A new natural reference material for U–Pb and Hf isotopic microanalysis: Chemical Geology, v. 249, p. 1–35, doi:10.1016/j.chemgeo.2007.11.005.
Steiger, R.H., and Jäger, E., 1977, Subcommission on geochronology: Convention on the use of decay constants in geo- and cosmochronology: Earth and Planetary Science Letters, v. 36, p. 359–362, doi:10.1016/0012-821X(77)90060-7.
Stern, R.A., Bodorkos, S., Kamo, S.L., Hickman, A.H., and Corfu, F., 2009, Measurement of SIMS Instrumental Mass Fractionation of Pb Isotopes During Zircon Dating: Geostandards and Geoanalytical Research, v. 33, p. 145–168, doi:10.1111/j.1751-908X.2009.00023.x.
Wiedenbeck, M., Allé, P., Corfu, F., Griffin, W.L., Meier, M., Oberli, F., Quadt, A.V., Roddick, J.C., and Spiegel, W., 1995, THREE NATURAL ZIRCON STANDARDS FOR U-TH-PB, LU-HF, TRACE ELEMENT AND REE ANALYSES: Geostandards and Geoanalytical Research, v. 19, p. 1–23, doi:10.1111/j.1751-908X.1995.tb00147.x.
Woodhead, J., 2002, A simple method for obtaining highly accurate Pb isotope data by MC-ICP-MS: Journal of Analytical Atomic Spectrometry, v. 17, p. 1381–1385, doi:10.1039/b205045e.
Woodhead, J.D., and Hergt, J.M., 2001, Strontium, Neodymium and Lead Isotope Analyses of NIST Glass Certified Reference Materials: SRM 610, 612, 614: Geostandards and Geoanalytical Research, v. 25, p. 261–266, doi:10.1111/j.1751-908X.2001.tb00601.x.
```
