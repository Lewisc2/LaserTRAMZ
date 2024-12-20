# LaserTRAMZ

## Purpose
LaserTRAMZ was written to provide a free, non-automated open-source program to handle time resolved U-Pb zircon analyses gathered via LA-Q-ICP-MS. This is effectively the sister to [LaserTRAM-DB](https://github.com/jlubbersgeo/laserTRAM-DB#readme) ([Lubbers et al., 202x](https://doi.org/10.5281/zenodo.7826697)) which was built for handling trace elements using the same instrumentation.

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

## Part One: Analyte Reduction
Running the script for Part One will open a page in your default web browser with a grey column and some default selections. `Copy + Paste` the file path to your LT_ready file into the bottom string input. After some loading time, an empty data table will pop up and the dropdown titled `Sample Subset` will be populated with the analyses in the file. Select one of these (it does not auto populate) and the program will load the data onto the screen:

![Screenshot 2024-12-19 at 10 37 26 AM](https://github.com/user-attachments/assets/137a0702-1420-48d0-979e-777f54afe614)

### Descriptions of GUI visualization and tools


#### A. Isotope Ratio Graph
Ratios can be toggled on and off by using the check boxes next to the isotope ratios in the left sidebar. Dashed lines show background start and end. Solid lines show ablation interval start and end. These can be changed using the sliders in the left sidebar. *IMPORTANT: *Dotted line shows where the ablation is projected to in order to get the Pb/U ratio prior to any laser induced elemental fractionation (LIEF) at the ablation site. You have the option to project before the ablation interval. This is an uncommon implementation but does hold with the assumptions of downhole correcting for LIEF (Kosler et al., 2002). You can change where the projected start is by using hte 'Ablation start true' box in the sidebar. Alternatively you may wish to check hte 'Lock Back Projection' box in the left sidebar, which binds 'Ablation start true' to the start of the selected ablation interval.
Y-axes limits can be changed using the Ratio plot ylim min / mix.
You man zoom in and out on this graph. Zoom tools, graph download, and pan tool are in the small ribbon above the plot.

#### B. Residuals Graph
Residuals of the chosen regressions for Pb/U ratios.
You may zoom and pan on this graph as well.

#### C. Regression Statistics
Regression statistics for Pb/U ratios showing R-squared values and 1SE % errors.

#### D. Ablation Data
Time resolved analyte intensities. Vertical lines as in the ratio plot. Visualized analytes can be chagned by selecting them from the list in the sidebar. The y-axis can be put on the log scale using the corresponding checkbox as well.
Zooming around and panning is also possible on this graph. 

#### E. Confidence Ellipsoids
Confidence ellipsoids derived from the covariance matrix. One can view the confidence ellipse for Tera-Wasserburg (TW) or Wetherhill Concordia (Weth.) by using the tabs above the plot. The power associated with the confidence of the ellipse can be changed by changing the value in the 'Power' input (default = 95% confidence). If using means and regression, the Pb/U ratios are downhole corrected for each detector pass before getting the axes of the ellipsoid. If using Total Counts, the Pb/U ratios are not downhole corrected but the center of the ellipsoid is in the second part of LaserTRAMZ. Typically, the downhold corrected ellipsoid will result in more elongate ellipsoids on Wetherhill Concordia.



## Part Two: Age Reduction
### Uploading and Setting reduction options
Running the script for Part Two will open up a page in your default web browser. Copy and Paste the filepath from part one into the file path input at the top of the sidebar. A window will pop up that looks like the following screenshot. Annotations in bold text explain what the options are.

![Screenshot 2024-12-19 at 10 57 14 AM](https://github.com/user-attachments/assets/95e5dca2-d51a-4def-9ccc-a3983fd3913e)

If the reference materials are named appropriately (see above) they will appear in a dropdown when you begin typing in the first two input boxes. Selecting other secondary standards will cause them to be output in their own excel file with various graphs highlighting variation in analyte intensity, age, etc. The primary standard is the primary normalization standard. The Secondary standard selected to quantify reference material (RM) isotope ratio uncertainty is evaluated using excess variance (Horstwood et al., 2016). The first selection at the bottom of the window will add error iteratively to both the 206Pb/238U and 207Pb/206Pb ratios (in increments of 0.01%) until the ages of hte secondary RM have an MSWD of 1. The second option will add error iteratively the the normalized 206Pb/238U and 207Pb/206Pb separately until the ratios have MSWDs of 1. The third option will do the same but only for the primary (unnormalized) reference material.

The sliding window fractionation factor is handled via Gehrels et al. (2008).

### Descriptions of GUI visualization and tools
After selecting the Accept Reduction Parameters button, input a sample string into the 'Text Sample Selector' box in the sidebar. The input string just needs to be some part of the sample name in the dataset. For example, if your samples are named CTL24-CHX-xyz, then inputting CTL24 will bring up all samples with that pre-fix. More detailed strings will be accordingly more restrictive. A screen similar to the following two screenshots will appear.

![Screenshot 2024-12-19 at 11 01 18 AM](https://github.com/user-attachments/assets/148d7749-bf80-49e0-be55-040ffc17825d)
![Screenshot 2024-12-19 at 11 01 35 AM](https://github.com/user-attachments/assets/352bd9b1-751d-4191-b501-ce8287f74ee5)

#### A. Concordia
Tera-Wasserburg and Wetherhill Concordia can be visualized using the taps at the top of the graph. X and Y axis limits can be changed using the slider in the sidebar. Lines projecting through data points onto Tera-Wasserburg Concordia are projected from Common Pb onto Concordia and visually show the removal of common Pb using the 207Pb/206Pb Tera-Wasserburg method (i.e., Vermeesch, 2018; and references therein). Yellow diamonds show data. Beige circles show concordant ratios used to calculate ages. Naturally these latter points do not appear on Wetherhill Concordia. Sample labels can be toggled on and off using the checkbox at the bottom of the sidebar.

#### B. Boxplot
Standard boxplot showing distribution of data. Fliers shown as hollow circle. Mean shown as yellow diamond. Median shown as yellow line.

#### C. Drift Assessment Graph
Graph showing selected ratio or age versus measurement number in the analytical session for all secondary standards. You may change what is shown using the Drift analyte dropdown near the top of the sidebar.

#### D. Excess Variance Graph
Three graphs selectable by tab showing the secondary or primary RM (depending on your selections) vs measurement number in the analytical session. Black bars are 2s errors. Red bars are 2s errors with added excess variance. Shaded blue-gray bar with solid line are weighted error and mean, respectively. If you chose to correct for drift using hte sliding window, the dashed bold vertical line will separate the windows and each will have its own weighted error and mean.

#### Other important notes
If you have determined the common Pb ratio via external measurement or projection of data you may input the ratio and its error in the appropriate boxes in the sidebar. Ratios for Th/U in zircon and magma used to calculate the Th disequilibrum correction are also available. Use the selector to choose whether data should be corrected only with Common Pb or with Common Pb and the Th/U disequilibrum correction.
Click the approve data button to calculate ages for the currently selected data. Click the DDDT! button to output data


## Final Output Information
* t start, t end, t project / bstart, bend: The ablation & background start, end, and projected value.
* Analytes (e.g., 206Pb, 238U) and analyte errors (e.g., 206Pb_1SE, 238U_1SE): Mean background subtracted intensities and 1*Standard errors for each analyte.
* Isotope ratios and their 1SE% are output. Pb/U ratios have their reduction method (i.e., regression type) included in the headers.
* Weth / TW C, Wid1, Wid2, rho: center, axes widths, and degree of rotation for confidence ellipsoids on Wetherhill and Tera-Wasserburg Concordia
* Isotope ratios with 'c' after them: Mass bias corrected isotope ratios. Also includes common Pb and Th-disequilibrium correction wehre present
* Pb/U ratio 'Age': Age in years
* Pb/U Age 1s (meas.): 1s analytical error on the age in years
* Pb/U Age 1s (tot): 1s full error on the age in year
  
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

## Demos
(https://www.youtube.com/watch?v=683U5F1hsdM&list=PLs1w2r4kCIlWRb5ARSq0PzIfl4XBUQuaG)

## Feedback and Questions
Feecback and suggestions may be made by opening an [issue](https://github.com/Lewisc2/LaserTRAMZ/issues) or emailing Chuck Lewis (lewisc2 _at_ oregonstate.edu). I'm happy to answer any questions that may come up.

## Citing and Related Documentation
If you use LaserTRAMZ in your work, please citte it! See citation.cff file for citing information. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8097362.svg)](https://doi.org/10.5281/zenodo.8097362)

- Publication coming `soon`!
  
## Publications utilizing LaserTRAMZ

- Lewis et al., (_in prep_)



## References
