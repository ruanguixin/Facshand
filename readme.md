# How to use facshand

Are you bored to do the FACS data statistic post flowjo handling? 

**USE FACSHAND!**

## Quick start

Facshand script is based on R for statistic and plot of FACS data exported from Flowjo. You can easily caculate and plot your FACS data by typing the following code in your terminal.

```shelll
Rscript path/facshand.R path_dir/
```

or run the script on R console as the following codes:

```R
path <- "path_dir"
source("path/facshand.R", chdir = TRUE)
```

You can directly drag your script file and directory to the console to acheive the paths.

Make sure that R has been installed in your OS with **ggplot2** and **Hmisc** packages installed.

```R
install.packages("ggplot2")
install.packages("Hmisc")
install.packages("acepack")
```

## Input files

The input files are as follows, keep these files in your path_dir directory.

1. rawdata.csv 

   transform xlsx file derived from flowjo to csv;

2. celltype.csv 

   with the header of "celltype", input the **exact names** of the subsets to be calculated; 

3.  group.csv

   with the header of "group", input the **exact names** of your groups;

4. cellnumber.csv

   with two colunms of "tissue" and "number", tissue refers to the exact name of your FACS sample names, and number refers to the cell number (x 10<sup>4</sup>) for caculation; notice that the sample name should be **prefixed with group name**, eg WT_whatever, ko-dmso-whatever.

\* make sure that your subset names are not duplicated, even partially duplicated names are not allowed, ignore cases.

## Output files

The output files include

1. data.csv

   contains the population and cell number caculated data;

2. statistic.txt

   refers to the statistics for the compared groups using students' t test between two groups. If multi groups are inputted, one-way anova is used for statistic;

3. plots

   named by your celltype names in tiff format, including population and cell number caculated dot plots with p-values. The plots may not meet the submission demands, but are enough for the presentation in lab meeting (maybe) :p .



**20200923_update**: modify to make x axis display as the same sequence as group.csv shown, modify to tolerate the most frequently used symbol in FACS gating "+".

**20200812_update**: support the rawdata files exported from flowjo version <= 10.1, support direct run on R console, add notices.