# LaserTRAMZ
## Citing and Related Documentation
#### doi
#### Publication coming 'soon'!

## Purpose
#### This program was written to handle time resolved U-Pb zircon data using LA-ICP-Quadrupole-MS. This is effectively the sister to [LaserTRAM-DB](https://github.com/jlubbersgeo/laserTRAM-DB#readme) ([Lubbers et al., 202x](https://doi.org/10.5281/zenodo.7826697)) which was built for handling trace elements using the same instrumentation. The program was born from a need to measure young zircons on the quadrupole and is built to handle isotopic ratios for the analyst that wants to count as long as possible during analysis (i.e., no trace elements). A variety of reduction techniques are available to reduce isotopic ratios, as described below.

## Data Input Format
**Part One**
#### Similar to LaserTRAM-DB, the [multifiler](https://github.com/jlubbersgeo/multifiler) repository  make the data prepared for LaserTRAMZ after deleting the 'Timestamp' header that is output in the LT_ready file. The input format should look like:
| SampleLabel        | Time           | 202Hg    | ... | 238U|
| ------------------ | -------------- | -------- | --- | --- |
| FishCanyon         | 12.74          | 100.0004 | ... | 0   |
| FishCanyon         | 158.31         | 400.0064 | ... | 0   |
| ...                | ...            | ...      | ... | ... |

**Part Two**
Primary and validation standards have a required input format in order to get the proper age reduction. Adding standards can be done with ease (see contact information below). Currently, standards are from the PlasmAge Consortium:
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
The output excel file from part one will have a column that needs to be deleted (column A) and a row that needs to be deleted (Row 2). Once this is done the file can be uploaded into Part Two and should look like:

| SampleLabel | t start | t end | ... | SE% 238/235 |
| ----------- | ------- | ------| --- | ----------- |
| FishCanyon  | 32      | 60    | ... | 1.04        |
| FishCanyon  | 32      | 57    | ... | 0.85        |
| ...         | ...     | ...   | ... | ...         |


## Part One: Analyte Reduction
Running the script for Part One will open a page in your default web browser

## Part Two: Age Reduction


## Installation and Use


## Demos


## Feedback and Questions
Feecback and suggestions may be made by opening an [issue](https://github.com/Lewisc2/LaserTRAMZ/issues) or emailing Chuck Lewis (<lewisc2@oregonstate.edu>). I'm happy to answer any questions that may come up.

## To-Do List
* Wetherhill Concordia
* Auto Session wide Drift Correction
* Demo Videos
* Get this GD paper published
