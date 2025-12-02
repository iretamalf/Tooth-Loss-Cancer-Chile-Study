* ==============================================================================
* PROJECT: Tooth Loss and Cancer History in Chile: A National Cross-Sectional Study
* DATA: Third Chilean National Health Survey (ENS 2016-2017)
* ==============================================================================
* DESCRIPTION: 
* This do-file performs data cleaning, variable generation (Social Vulnerability,
* Expanded Allostatic Load, Oral Health), and statistical analysis (Survey-weighted
* Logistic Regression, Splines, KHB Mediation, E-values).
* ==============================================================================

clear all
set more off
macro drop _all

* ------------------------------------------------------------------------------
* 1. SETUP & REQUIREMENTS
* ------------------------------------------------------------------------------
* IMPORTANT: Set your working directory here before running.
* The directory must contain the dataset "ENS_2016_2017.dta".
* Example: cd "C:/Users/Name/Documents/Research/ENS_Project"
* ------------------------------------------------------------------------------

* ---> INSERT YOUR CD COMMAND HERE <---

* Install required packages if missing
foreach pkg in coefplot khb evalue mkspline estout {
    capture which `pkg'
    if _rc != 0 {
        display "Installing package: `pkg'..."
        ssc install `pkg', replace
    }
}

* Start Logging
log using "Analysis_Log.log", replace text
display ">>> ANALYSIS STARTED: " c(current_date) " " c(current_time)

* Load Data (Ensure the file is in your working directory)
use "ENS_2016_2017.dta", clear

* ------------------------------------------------------------------------------
* 2. VARIABLE GENERATION
* ------------------------------------------------------------------------------
display ">>> GENERATING VARIABLES..."

* --- A. Social Vulnerability Index (SVI) ---
* Components: Income, Education, Housing (Overcrowding)
gen income_vuln = (as27 < 158145) if !missing(as27)
gen edu_vuln = (NEDU1_MINSAL_1 < 3) if !missing(NEDU1_MINSAL_1)
gen housing_vuln = (as36 >= 2.5) if !missing(as36)

* Calculate Index (Mean of components)
egen svi_score = rowmean(income_vuln edu_vuln housing_vuln)
egen svi_miss = rowmiss(income_vuln edu_vuln housing_vuln) 
replace svi_score = . if svi_miss > 1
label variable svi_score "Social Vulnerability Index"

* --- B. Expanded Allostatic Load (AL) ---
* Score range: 0-7. Adjusts for medication use.

* 1. Blood Pressure (Biomarker + Meds)
gen pas_final = (m2p9_1 + m2p10_1) / 2 
gen pad_final = (m2p9_2 + m2p10_2) / 2
gen meds_bp = 0
capture replace meds_bp = 1 if m2p4 == 1 // Treatment for Hypertension
gen risk_bp = (pas_final >= 140 | pad_final >= 90 | meds_bp == 1) if !missing(pas_final)

* 2. Cholesterol (Biomarker + Meds)
gen risk_col = (Colesterol_Total >= 200) if !missing(Colesterol_Total)
capture replace risk_col = 1 if m3p18 == 1 // Treatment for High Cholesterol
gen risk_hdl = 0
replace risk_hdl = 1 if (Sexo == 1 & Colesterol_HDL < 40) | (Sexo == 2 & Colesterol_HDL < 50)
replace risk_hdl = . if missing(Colesterol_HDL)

* 3. Glucose (HbA1c + Meds)
gen risk_gluc = .
gen meds_dm = 0
capture replace meds_dm = 1 if m3p5 == 1 // Treatment for Diabetes
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
replace al_score = . if al_miss > 3
label variable al_score "Allostatic Load Score"

* --- C. Oral Health (Tooth Loss) ---
gen teeth_total = .
capture confirm variable Totalf2ajustado
if _rc == 0 replace teeth_total = Totalf2ajustado 
else replace teeth_total = m5p3 + m5p6 
replace teeth_total = . if teeth_total < 0
gen missing_teeth = 32 - teeth_total
label variable missing_teeth "Missing Teeth"

* --- D. Covariates & Outcome ---
gen sex_en = Sexo
label define l_sex 1 "Male" 2 "Female"
label values sex_en l_sex

gen smoking = 0 
capture confirm variable ta2
if _rc == 0 replace smoking = 1 if ta2 == 1 | ta2 == 2 
else { 
    capture confirm variable ta1
    if _rc == 0 replace smoking = 1 if ta1 == 1 
    else capture replace smoking = 1 if ta14_3a == 1 
}

gen cancer_outcome = .
capture confirm variable m9p5A
if _rc == 0 {
    replace cancer_outcome = 1 if m9p5A == 1
    replace cancer_outcome = 0 if m9p5A == 2
}
label define l_ca 0 "No Cancer" 1 "History of Cancer"
label values cancer_outcome l_ca

* ------------------------------------------------------------------------------
* 3. STATISTICAL ANALYSIS
* ------------------------------------------------------------------------------
* Set Survey Design (Weights)
svyset Conglomerado [pweight=Fexp_F1p_Corr], strata(Estrato) singleunit(centered)

* --- Table 1: Descriptive Statistics ---
display ">>> GENERATING TABLE 1..."
dtable Edad sex_en smoking svi_score al_score missing_teeth, ///
    by(cancer_outcome, tests) svy nformat(%9.2f mean se) 

* --- Table 2: Multivariable Models ---
display ">>> RUNNING REGRESSION MODELS..."

* Model A: Allostatic Load
svy: logistic cancer_outcome c.svi_score c.al_score i.smoking c.Edad i.sex_en
estimates store Model_AL

* Model B: Oral Health
* Note: Attempt to use dental-specific weights if available
capture svyset Conglomerado [pweight=Fexp_EX2p_Corr], strata(Estrato) singleunit(centered)
if _rc != 0 svyset Conglomerado [pweight=Fexp_F1p_Corr], strata(Estrato) singleunit(centered) 

svy: logistic cancer_outcome c.svi_score c.missing_teeth i.smoking c.Edad i.sex_en
estimates store Model_Teeth

* Display Table Results
esttab Model_AL Model_Teeth, eform ci(2) label star(* 0.05 ** 0.01 *** 0.001) scalars(N)

* ------------------------------------------------------------------------------
* 4. SENSITIVITY ANALYSES
* ------------------------------------------------------------------------------
display ">>> SENSITIVITY ANALYSES..."
svyset Conglomerado [pweight=Fexp_F1p_Corr], strata(Estrato) singleunit(centered)

* A. Non-Linearity Check (Splines)
capture drop al_spline*
capture mkspline al_spline = al_score, cubic nknots(3)
capture svy: logistic cancer_outcome al_spline* i.smoking c.Edad i.sex_en c.svi_score
capture testparm al_spline* * B. Robustness (Excluding Edentulous)
gen not_edentulous = (missing_teeth < 32)
svy, subpop(not_edentulous): logistic cancer_outcome c.missing_teeth c.svi_score i.smoking c.Edad i.sex_en

* C. Interaction (Age x Tooth Loss)
svy: logistic cancer_outcome c.missing_teeth##c.Edad i.smoking i.sex_en c.svi_score

* D. E-Values (Bias Analysis for OR ~1.22)
evalue or 1.22, lcl(1.10) ucl(1.35)

* ------------------------------------------------------------------------------
* 5. GRAPHICS
* ------------------------------------------------------------------------------
display ">>> GENERATING FIGURES..."

* Figure 1: Forest Plot
set scheme s1mono
coefplot (Model_AL, label("Allostatic Load Model")) ///
         (Model_Teeth, label("Oral Health Model")), ///
    drop(_cons *.sex_en *.smoking c.Edad) xline(1, lpattern(dash)) eform ///
    title("Factors Associated with Cancer History") name(Fig1_Forest, replace)
graph export "Figure1_Forest.png", replace width(2400)

* Figure 2: Margins Plot (Probability)
quietly svy: logistic cancer_outcome c.missing_teeth i.smoking c.Edad i.sex_en c.svi_score
margins, at(missing_teeth=(0(5)32))
marginsplot, recast(line) ciopts(recast(rarea) color(gs14)) ///
    title("Estimated Probability of Cancer History by Tooth Loss") ///
    ytitle("Estimated Probability") xtitle("Missing Teeth") ///
    name(Fig2_Margins, replace)
graph export "Figure2_Margins.png", replace width(2400)

display ">>> ANALYSIS COMPLETED SUCCESSFULLY."
log close
