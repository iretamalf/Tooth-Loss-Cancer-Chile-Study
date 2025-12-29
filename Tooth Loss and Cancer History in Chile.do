* ==============================================================================
* PROJECT: Tooth Loss and Cancer History in Chile (ENS 2016-2017)
* JOURNAL: Oral Diseases
* FILE: Final_Analysis_English_Corrected.do
* AUTHOR: Research Team
* DATE: December 2025
* STATA VERSION: 19
* ==============================================================================

cls
clear all
set more off
set linesize 255 

* ------------------------------------------------------------------------------
* 1. ENVIRONMENT SETUP
* ------------------------------------------------------------------------------

* Attempt to set directory automatically
capture cd "/Users/ignacioretamalfarina/Box Sync/IR/Research/UMAYOR/COP/Ideas proyectos/Social_adversity/ENS/Data bruta ENS/Exploracion"
if _rc != 0 {
    capture cd "/Users/iretamalf/Library/CloudStorage/Box-Box/IR/Research/UMAYOR/COP/Ideas proyectos/Social_adversity/ENS/Data bruta ENS/Exploracion"
}

* Start Logging
capture log close
log using "Final_Analysis_Log.log", replace text
display ">>> LOGGING STARTED: Final_Analysis_Log.log"

* Install coefplot if missing
capture which coefplot
if _rc != 0 {
    ssc install coefplot, replace
}

* Load Data
use "ENS_2016_2017.dta", clear

* ==============================================================================
* 2. VARIABLE GENERATION & CLEANING
* ==============================================================================
display ">>> STEP 2: GENERATING VARIABLES..."

* --- A. Social Vulnerability Index (SVI) ---
gen income_vuln = (as27 < 158145) if !missing(as27)
gen edu_vuln = (NEDU1_MINSAL_1 < 3) if !missing(NEDU1_MINSAL_1)
gen housing_vuln = (as36 >= 2.5) if !missing(as36)

egen svi_score = rowmean(income_vuln edu_vuln housing_vuln)
egen svi_n = rownonmiss(income_vuln edu_vuln housing_vuln)
replace svi_score = . if svi_n < 2
label variable svi_score "Social Vulnerability Index (0-1)"

* --- B. Allostatic Load (AL) ---
* Blood Pressure
gen pas_final = (m2p9_1 + m2p10_1) / 2 
gen pad_final = (m2p9_2 + m2p10_2) / 2
gen risk_bp = (pas_final >= 140 | pad_final >= 90) if !missing(pas_final)

* Cholesterol
gen risk_col = (Colesterol_Total >= 200) if !missing(Colesterol_Total)

* HDL (Sex-specific)
gen risk_hdl = .
replace risk_hdl = 1 if (Sexo == 1 & Colesterol_HDL < 40) | (Sexo == 2 & Colesterol_HDL < 50)
replace risk_hdl = 0 if !missing(Colesterol_HDL) & missing(risk_hdl)

* Glucose (HbA1c)
capture confirm variable aux_Hemoglobina_A1C
if _rc == 0 {
    gen risk_gluc = (aux_Hemoglobina_A1C >= 6.5) if !missing(aux_Hemoglobina_A1C)
}
else {
    gen hba1c_num = real(m3p10)
    gen risk_gluc = (hba1c_num >= 6.5) if !missing(hba1c_num)
}

* Waist Circumference
gen risk_waist = .
replace risk_waist = 1 if (Sexo == 1 & m4p3 >= 102) | (Sexo == 2 & m4p3 >= 88)
replace risk_waist = 0 if !missing(m4p3) & missing(risk_waist)

* CRP (Inflammation)
capture confirm variable aux_Proteina_C_Reactiva_cuant
if _rc == 0 {
    gen risk_pcr = (aux_Proteina_C_Reactiva_cuant >= 3) if !missing(aux_Proteina_C_Reactiva_cuant)
}
else {
    gen risk_pcr = . 
}

* Final AL Score
egen al_score = rowtotal(risk_bp risk_col risk_hdl risk_gluc risk_waist risk_pcr)
egen al_miss = rowmiss(risk_bp risk_col risk_hdl risk_gluc risk_waist risk_pcr)

replace al_score = . if al_miss > 4
label variable al_score "Allostatic Load Score"

* --- C. Oral Health (Tooth Loss) ---
gen teeth_total = .
capture confirm variable Totalf2ajustado
if _rc == 0 { 
    replace teeth_total = Totalf2ajustado 
}
else { 
    replace teeth_total = m5p3 + m5p6 if !missing(m5p3) & !missing(m5p6) 
}

replace teeth_total = . if teeth_total < 0
gen missing_teeth = 32 - teeth_total
label variable missing_teeth "Number of Missing Teeth"

* --- D. Covariates (ENGLISH LABELS) ---
* Sex
gen sex_en = Sexo
label define l_sex_en 1 "Male" 2 "Female"
label values sex_en l_sex_en
label variable sex_en "Sex"

* Age (Fixing 'Edad' to 'Age' for English output)
clonevar age = Edad
label variable age "Age"

gen age_group = .
replace age_group = 1 if age < 45
replace age_group = 2 if age >= 45 & age < 65
replace age_group = 3 if age >= 65
label define l_age_en 1 "15-44y" 2 "45-64y" 3 "65+y"
label values age_group l_age_en
label variable age_group "Age Group"

* Smoking
gen smoking = 0 
capture replace smoking = 1 if ta2 == 1 | ta2 == 2 
label define l_smoke 0 "Non-Smoker" 1 "Current Smoker"
label values smoking l_smoke
label variable smoking "Tobacco Use"

* ==============================================================================
* 3. OUTCOME: COMPREHENSIVE CANCER HISTORY
* ==============================================================================
display ">>> STEP 3: CONSTRUCTING OUTCOME..."

gen cancer_history_total = 0

* WOMEN'S MODULE (M9)
capture replace cancer_history_total = 1 if m9p2A == 1
capture replace cancer_history_total = 1 if m9p3A == 1
capture replace cancer_history_total = 1 if m9p4A == 1
capture replace cancer_history_total = 1 if m9p5A == 1
capture replace cancer_history_total = 1 if m9p6A == 1
capture replace cancer_history_total = 1 if m9p7A == 1
capture replace cancer_history_total = 1 if m9p8A == 1

* MEN'S MODULE (M11)
capture replace cancer_history_total = 1 if m11p2 == 1
capture replace cancer_history_total = 1 if m11p3 == 1
capture replace cancer_history_total = 1 if m11p2A == 1

label define l_ca_total 0 "No Cancer History" 1 "History of Cancer"
label values cancer_history_total l_ca_total
label variable cancer_history_total "Cancer Diagnosis"

* CHECK
display "******************************************************"
display "*** REAL CASE COUNT (UNWEIGHTED N) ***"
tab cancer_history_total, missing
display "******************************************************"

* ==============================================================================
* 4. GENERATING TABLES 1, 2, 3, 4
* ==============================================================================
display ">>> STEP 4: GENERATING TABLES..."

* --- Table 1: Unweighted Sample Characteristics ---
dtable age_group sex_en smoking svi_score al_score missing_teeth, ///
    by(cancer_history_total, tests) ///
    title("Table 1. Unweighted Characteristics of the Sample (ENS 2016-2017)") ///
    note("Values are presented as raw counts (n) and unweighted column percentages (%). This table reflects the actual sample size observed.") ///
    export("Table1_Unweighted.docx", replace)

* --- Table 2: Weighted Population Characteristics ---
svyset Conglomerado [pweight=Fexp_F1p_Corr], strata(Estrato) singleunit(centered)

dtable age_group sex_en smoking svi_score al_score missing_teeth, ///
    by(cancer_history_total, tests) svy ///
    title("Table 2. Weighted Characteristics of the Study Population") ///
    note("Values are presented as weighted means (SE) and weighted population percentages. Estimates account for complex survey design.") ///
    export("Table2_Weighted.docx", replace)

* --- Table 3: Unweighted Logistic Regression (Sensitivity) ---
logistic cancer_history_total c.svi_score c.al_score i.smoking c.age i.sex_en
estimates store Unweighted_Bio
logistic cancer_history_total c.svi_score c.missing_teeth i.smoking c.age i.sex_en
estimates store Unweighted_Oral

etable, estimates(Unweighted_Bio Unweighted_Oral) ///
    column(index) cstat(_r_b) cstat(_r_ci) mstat(N) showstars showstarsnote ///
    title("Table 3. Unweighted Logistic Regression Models (Sensitivity Analysis)") ///
    note("Odds Ratios (OR) and 95% Confidence Intervals (CI) derived from logistic regression models without survey weights.") ///
    export("Table3_Unweighted_Reg.docx", replace)

* --- Table 4: Survey-Weighted Logistic Regression (Main Results) ---
* NOTE: In svy models, "N" refers to the Sample Size, not the Weighted Population.
svy: logistic cancer_history_total c.svi_score c.al_score i.smoking c.age i.sex_en
estimates store Weighted_Bio
svy: logistic cancer_history_total c.svi_score c.missing_teeth i.smoking c.age i.sex_en
estimates store Weighted_Oral

etable, estimates(Weighted_Bio Weighted_Oral) ///
    column(index) cstat(_r_b) cstat(_r_ci) mstat(N) showstars showstarsnote ///
    title("Table 4. Survey-Weighted Logistic Regression Models (Main Results)") ///
    note("Odds Ratios (OR) and 95% Confidence Intervals (CI) derived from complex survey-weighted logistic regression models.") ///
    export("Table4_Weighted_Reg.docx", replace)

* ==============================================================================
* 5. FIGURES (ENGLISH LABELS)
* ==============================================================================
display ">>> STEP 5: GENERATING FIGURES..."

set scheme s1mono

* Figure 1: Forest Plot (Based on Main Weighted Models)
* Dropping 'age' and categorical controls to focus on main predictors
coefplot (Weighted_Bio, label("Model A: Allostatic Load")) ///
         (Weighted_Oral, label("Model B: Oral Health")), ///
    drop(_cons *.sex_en *.smoking c.age 0.smoking 1.smoking 1.sex_en 2.sex_en) ///
    xline(1, lpattern(dash) lcolor(gray)) eform levels(95) ///
    xtitle("Odds Ratio (95% CI)") ///
    headings(c.al_score = "{bf:Biological Markers}" ///
             c.missing_teeth = "{bf:Oral Health}" ///
             c.svi_score = "{bf:Social Determinants}") ///
    name(ForestPlot, replace)
graph export "Figure1_ForestPlot.tif", replace width(2400)

* Figure 2: Bar Plot
svy: regress missing_teeth i.cancer_history_total
margins cancer_history_total

marginsplot, ///
    recast(bar) ///
    plotopts(barwidth(0.6) fcolor(gs12) lcolor(black)) ///
    ciopts(lcolor(black)) ///
    title("Average Tooth Loss by Cancer History") ///
    ytitle("Mean Number of Missing Teeth") ///
    xtitle("") ///
    xlabel(0 "No Cancer" 1 "Cancer History", noticks) ///
    name(BarPlot, replace)
    
graph export "Figure2_BarChart_Teeth.tif", replace width(2400)

display " "
display ">>> COMPLETED SUCCESSFULLY."
display ">>> Generated: Table1_Unweighted.docx, Table2_Weighted.docx"
display ">>> Generated: Table3_Unweighted_Reg.docx, Table4_Weighted_Reg.docx"
display ">>> Generated: Figure1_ForestPlot.tif, Figure2_BarChart_Teeth.tif"

log close
