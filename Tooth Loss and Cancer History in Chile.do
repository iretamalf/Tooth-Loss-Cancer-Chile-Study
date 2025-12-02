* ==============================================================================
* PROJECT: Tooth Loss and Cancer History in Chile: A National Cross-Sectional Study
* DATA: Third Chilean National Health Survey (ENS 2016-2017)
* ==============================================================================
* AUTHOR:  Ignacio N. Retamal
* EMAIL:   ignacio.retamal@umayor.cl
* DATE:    November 2025
* ==============================================================================
* DESCRIPTION:
* This do-file performs the complete analysis pipeline:
* 1. Data Cleaning & Variable Generation (SVI, Expanded Allostatic Load, Oral Health).
* 2. Descriptive Statistics (Table 1).
* 3. Survey-Weighted Logistic Regression Models (Table 2).
* 4. Sensitivity Analyses (Splines, Mediation, Robustness checks).
* 5. Figure Generation (Forest Plot, Margins Plot).
* ==============================================================================

capture log close
clear all
set more off
macro drop _all

* ------------------------------------------------------------------------------
* 1. SETUP & REQUIREMENTS
* ------------------------------------------------------------------------------
* IMPORTANT:
* Before running this script, you must set the working directory to the folder 
* where you have saved the "ENS_2016_2017.dta" dataset.
* ------------------------------------------------------------------------------

* ---> INSERT YOUR PATH BELOW <---
* Example: cd "C:/Users/Name/Documents/Research/ENS_Project"
* cd "..." 

* Initialize Log File
log using "Analysis_Log.log", replace text
display ">>> ANALYSIS STARTED: " c(current_date) " " c(current_time)

* Install required community-contributed packages
foreach pkg in coefplot khb evalue mkspline estout {
    capture which `pkg'
    if _rc != 0 {
        display "Installing package: `pkg'..."
        ssc install `pkg', replace
    }
}

* ------------------------------------------------------------------------------
* 2. DATA LOADING & VARIABLE GENERATION
* ------------------------------------------------------------------------------
* Load the dataset (Ensure the filename matches your local copy)
use "ENS_2016_2017.dta", clear

display ">>> GENERATING VARIABLES..."

* --- A. SOCIAL VULNERABILITY INDEX (SVI) ---
* Components: Income, Education, Housing (Overcrowding)
gen income_vuln = (as27 < 158145) if !missing(as27)
gen edu_vuln = (NEDU1_MINSAL_1 < 3) if !missing(NEDU1_MINSAL_1)
gen housing_vuln = (as36 >= 2.5) if !missing(as36)

* Calculate Index (0-1)
egen svi_score = rowmean(income_vuln edu_vuln housing_vuln)
egen svi_miss = rowmiss(income_vuln edu_vuln housing_vuln) 
replace svi_score = . if svi_miss > 1
label variable svi_score "Social Vulnerability Index"

* --- B. EXPANDED ALLOSTATIC LOAD (AL) ---
* Score range: 0-7. Includes Cardiovascular, Metabolic, Immune, Renal, and Autonomic markers.
* Adjusted for self-reported medication use to avoid misclassification.

* 1. Blood Pressure (Biomarker + Medication)
gen pas_final = (m2p9_1 + m2p10_1) / 2 
gen pad_final = (m2p9_2 + m2p10_2) / 2
gen meds_bp = 0
capture replace meds_bp = 1 if m2p4 == 1 // "Receiving HTN treatment"
gen risk_bp = (pas_final >= 140 | pad_final >= 90 | meds_bp == 1) if !missing(pas_final)

* 2. Cholesterol (Biomarker + Medication)
gen risk_col = (Colesterol_Total >= 200) if !missing(Colesterol_Total)
capture replace risk_col = 1 if m3p18 == 1 // "Receiving Cholesterol treatment"
gen risk_hdl = 0
replace risk_hdl = 1 if (Sexo == 1 & Colesterol_HDL < 40) | (Sexo == 2 & Colesterol_HDL < 50)
replace risk_hdl = . if missing(Colesterol_HDL)

* 3. Glucose (HbA1c + Medication)
gen risk_gluc = .
gen meds_dm = 0
capture replace meds_dm = 1 if m3p5 == 1 // "Receiving Diabetes treatment"
capture confirm variable aux_Hemoglobina_A1C
if _rc == 0 replace risk_gluc = (aux_Hemoglobina_A1C >= 6.5 | meds_dm == 1) if !missing(aux_Hemoglobina_A1C)
else {
    gen hba1c_num = real(m3p10)
    replace risk_gluc = (hba1c_num >= 6.5 | meds_dm == 1) if !missing(hba1c_num)
}

* 4. Waist Circumference (Abdominal Obesity)
gen risk_waist = 0
replace risk_waist = 1 if (Sexo == 1 & m4p3 >= 102) | (Sexo == 2 & m4p3 >= 88)
replace risk_waist = . if missing(m4p3)

* 5. Inflammation (C-Reactive Protein)
gen risk_pcr = .
capture confirm variable aux_Proteina_C_Reactiva_cuant
if _rc == 0 replace risk_pcr = (aux_Proteina_C_Reactiva_cuant >= 3) if !missing(aux_Proteina_C_Reactiva_cuant)

* 6. Autonomic System (Resting Heart Rate)
gen risk_pulse = .
gen pulse_mean = .
capture replace pulse_mean = (m2p5 + m2p11)/2 
if _rc == 0 replace risk_pulse = (pulse_mean >= 90) if !missing(pulse_mean)

* 7. Renal Function (Serum Creatinine)
gen risk_renal = .
capture confirm variable aux_Creatinina
if _rc == 0 {
    replace risk_renal = 0 if !missing(aux_Creatinina)
    replace risk_renal = 1 if (Sexo == 1 & aux_Creatinina > 1.2) | (Sexo == 2 & aux_Creatinina > 1.0)
    replace risk_renal = . if missing(aux_Creatinina)
}

* Calculate Final AL Score
egen al_score = rowtotal(risk_bp risk_col risk_hdl risk_gluc risk_waist risk_pcr risk_pulse risk_renal)
egen al_miss = rowmiss(risk_bp risk_col risk_hdl risk_gluc risk_waist risk_pcr risk_pulse risk_renal)
* Exclude participants with >3 missing biomarkers
replace al_score = . if al_miss > 3
label variable al_score "Allostatic Load Score" 

* --- C. ORAL HEALTH (Tooth Loss) ---
gen teeth_total = .
capture confirm variable Totalf2ajustado
if _rc == 0 replace teeth_total = Totalf2ajustado 
else replace teeth_total = m5p3 + m5p6 
replace teeth_total = . if teeth_total < 0
gen missing_teeth = 32 - teeth_total
label variable missing_teeth "Missing Teeth"

* --- D. COVARIATES & OUTCOME ---
* Demographics
gen sex_en = Sexo
gen smoking = 0 
capture confirm variable ta2
if _rc == 0 replace smoking = 1 if ta2 == 1 | ta2 == 2 
else { 
    capture confirm variable ta1
    if _rc == 0 replace smoking = 1 if ta1 == 1 
    else capture replace smoking = 1 if ta14_3a == 1 
}

* Outcome: History of Cancer
gen cancer_outcome = .
capture confirm variable m9p5A
if _rc == 0 {
    replace cancer_outcome = 1 if m9p5A == 1
    replace cancer_outcome = 0 if m9p5A == 2
}
else {
    capture confirm variable m9p5a
    if _rc == 0 {
        replace cancer_outcome = 1 if m9p5a == 1
        replace cancer_outcome = 0 if m9p5a == 2
    }
}

* --- E. LABELING (English) ---
label variable Edad "Age (years)"
label variable sex_en "Sex"
label define l_sex 1 "Male" 2 "Female"
label values sex_en l_sex
label variable smoking "Smoking Status"
label define l_smoke 0 "Non-Smoker" 1 "Current Smoker"
label values smoking l_smoke
label variable svi_score "Social Vulnerability Index"
label variable al_score "Allostatic Load Score"
label variable missing_teeth "Missing Teeth (Count)"
label variable cancer_outcome "Cancer History"
label define l_cancer 0 "No Cancer" 1 "Cancer History"
label values cancer_outcome l_cancer

* Additional vars for Table 1
gen education_level = .
capture replace education_level = 0 if NEDU1_MINSAL_1 < 3
capture replace education_level = 1 if NEDU1_MINSAL_1 >= 3
label define l_edu 0 "< Complete Secondary" 1 ">= Complete Secondary"
label values education_level l_edu
label variable education_level "Education Level"

capture label define l_zona 1 "Urban" 2 "Rural"
capture label values Zona l_zona
label variable Zona "Zone"

* ------------------------------------------------------------------------------
* 3. STATISTICAL ANALYSIS
* ------------------------------------------------------------------------------
* Define Complex Survey Design
svyset Conglomerado [pweight=Fexp_F1p_Corr], strata(Estrato) singleunit(centered)

display ""
display "=================================================================="
display ">>> GENERATING TABLE 1: DESCRIPTIVE STATISTICS <<<"
display "=================================================================="

dtable Edad i.sex_en i.education_level i.Zona i.smoking ///
       svi_score al_score missing_teeth, ///
    by(cancer_outcome, tests) svy ///
    sample(, statistics(freq) place(seplabels)) ///
    sformat("(%s)" se) ///
    nformat(%9.1f mean se) ///
    title("Table 1. Weighted Characteristics of the Study Population") ///
    export("Table1_Descriptive.docx", replace)

display ""
display "=================================================================="
display ">>> GENERATING TABLE 2: REGRESSION MODELS <<<"
display "=================================================================="

* Model A: Allostatic Load
svy: logistic cancer_outcome c.svi_score c.al_score i.smoking c.Edad i.sex_en
estimates store Model_AL

* Model B: Oral Health
* Note: Attempt to use dental-specific weights if available
capture svyset Conglomerado [pweight=Fexp_EX2p_Corr], strata(Estrato) singleunit(centered)
if _rc != 0 svyset Conglomerado [pweight=Fexp_F1p_Corr], strata(Estrato) singleunit(centered) 

svy: logistic cancer_outcome c.svi_score c.missing_teeth i.smoking c.Edad i.sex_en
estimates store Model_Teeth

* Export Table 2 to Word (RTF format)
esttab Model_AL Model_Teeth using "Table2_Regression.rtf", ///
    replace ///
    eform ci(2) label ///
    star(* 0.05 ** 0.01 *** 0.001) ///
    scalars(N) ///
    title("Table 2. Multivariable Logistic Regression Models (Weighted)") ///
    nonumbers mtitles("Model A (Physio)" "Model B (Oral)") ///
    onecell rtf ///
    coeflabels(c.svi_score "Social Vulnerability Index" ///
               c.al_score "Allostatic Load Score" ///
               c.missing_teeth "Missing Teeth" ///
               c.Edad "Age (Years)" ///
               1.sex_en "Sex: Female" ///
               1.smoking "Smoking: Current" ///
               _cons "Intercept") ///
    drop(0.smoking) 

* ------------------------------------------------------------------------------
* 4. SENSITIVITY ANALYSES
* ------------------------------------------------------------------------------
display ""
display ">>> RUNNING SENSITIVITY CHECKS..."
* Restore main weights
svyset Conglomerado [pweight=Fexp_F1p_Corr], strata(Estrato) singleunit(centered)

* A. Non-Linearity Check (Restricted Cubic Splines for AL)
capture drop al_spline*
capture mkspline al_spline = al_score, cubic nknots(3)
capture svy: logistic cancer_outcome al_spline* i.smoking c.Edad i.sex_en c.svi_score
capture testparm al_spline* if _rc != 0 display ">>> NOTE: Splines calculation skipped due to collinearity/sample size."

* B. Robustness Check (Excluding Edentulous Participants)
gen not_edentulous = (missing_teeth < 32)
svy, subpop(not_edentulous): logistic cancer_outcome c.missing_teeth c.svi_score i.smoking c.Edad i.sex_en

* C. Interaction Check (Age x Tooth Loss)
svy: logistic cancer_outcome c.missing_teeth##c.Edad i.smoking i.sex_en c.svi_score

* D. Bias Analysis (E-Values)
* Note: Values are based on the observed OR from Model B
evalue or 1.22, lcl(1.10) ucl(1.35)

* ------------------------------------------------------------------------------
* 5. GRAPHICS GENERATION
* ------------------------------------------------------------------------------
display ">>> GENERATING FIGURES..."

* Figure 1: Forest Plot
set scheme s1mono
coefplot (Model_AL, label("Allostatic Load Model")) ///
         (Model_Teeth, label("Oral Health Model")), ///
    drop(_cons *.sex_en *.smoking c.Edad) xline(1, lpattern(dash)) eform ///
    title("Factors Associated with Cancer History") name(Fig1_Forest, replace)
graph export "Figure1_Forest.png", replace width(2400)

* Figure 2: Margins Plot (Estimated Probability)
quietly svy: logistic cancer_outcome c.missing_teeth i.smoking c.Edad i.sex_en c.svi_score
margins, at(missing_teeth=(0(5)32))
marginsplot, recast(line) ciopts(recast(rarea) color(gs14)) ///
    title("Estimated Probability of Cancer History by Tooth Loss") ///
    ytitle("Estimated Probability") xtitle("Missing Teeth") ///
    name(Fig2_Margins, replace)
graph export "Figure2_Margins.png", replace width(2400)

display ">>> ANALYSIS COMPLETED SUCCESSFULLY."
display "Outputs created:"
display "1. Table1_Descriptive.docx"
display "2. Table2_Regression.rtf"
display "3. Figure1_Forest.png"
display "4. Figure2_Margins.png"

log close
