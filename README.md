# Tooth Loss and Cancer History in Chile: A National Cross-Sectional Study

This repository contains the analytical code (Stata) used for the study: **"Tooth Loss and Cancer History in Chile: A National Cross-Sectional Study"**. 

This research analyzes the relationship between physiological dysregulation (Allostatic Load), oral health (Tooth Loss), and cancer history using data from the Chilean National Health Survey (ENS 2016-2017).

## Data Availability

The datasets analyzed during the current study are **publicly available** but are not hosted in this repository due to copyright and ownership regulations of the Chilean government.

Researchers can download the data directly from the Department of Epidemiology of the Chilean Ministry of Health (MINSAL) at:
**[https://epi.minsal.cl/bases-de-datos/](https://epi.minsal.cl/bases-de-datos/)**

* **Dataset required:** `ENS_2016_2017.dta` (or `.sav` converted to `.dta`)
* **Protocol:** Survey documentation and variable dictionaries are also available at the link above.

## Repository Contents

* `Analysis_Code.do`: The main Stata script that performs:
    1.  Data cleaning and variable generation (SVI, Allostatic Load, Tooth Loss).
    2.  Survey-weighted statistical analysis (Logistic Regression).
    3.  Sensitivity analyses (Splines, KHB Mediation, E-values).
    4.  Generation of Figures 1 and 2.

## Requirements

* **Software:** Stata SE or MP (Version 15 or higher is recommended).
* **Stata Packages:** The code automatically checks for and installs the following required community-contributed packages:
    * `coefplot`
    * `khb`
    * `evalue`
    * `mkspline`
    * `estout`

## Instructions for Use

1.  **Download the Data:** Go to the [MINSAL website](https://epi.minsal.cl/bases-de-datos/) and download the ENS 2016-2017 dataset. Save it as `ENS_2016_2017.dta` on your local machine.
2.  **Download the Code:** Clone this repository or download the `Analysis_Code.do` file.
3.  **Set Working Directory:** Open `Analysis_Code.do` in Stata and edit the `cd` command at the top of the script to point to the folder where you saved the dataset.
    ```stata
    * Example:
    cd "C:/MyDocuments/Research/Chile_Cancer_Study"
    ```
4.  **Run:** Execute the do-file. It will generate a log file (`Analysis_Log.log`) with the results and save the figures as PNG files in the same folder.

## Citation

If you use this code or methodology, please cite the original article:
> [Citation Placeholder: Authors. (Year). Tooth Loss and Cancer History in Chile: A National Cross-Sectional Study. Journal Name.]

## Contact

For questions regarding the code, please contact:
**Ignacio N. Retamal** Universidad Mayor  
Email: ignacio.retamal@umayor.cl
