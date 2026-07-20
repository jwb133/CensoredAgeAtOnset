clear

version 18.5
set seed 398268137
set sortseed 611887103

matrix M = (46, 46)
matrix COV = (81, 48.6 \ 48.6, 81)
* distributions for actual and parental onset

forvalues i = 1/10000{
qui{
clear
drawnorm actual_onset parental_onset, n(100) cov(COV) means(M)
label var actual_onset "age at actual onset"
label var parental_onset "age at parental onset"

gen age=40 + rnormal(0,10)
gen id=_n

gen time_rel_actual = age - actual_onset  
label var time_rel_actual "time relative to actual onset"

gen y = 100 + time_rel_actual + rnormal(0,27)
label var y "outcome"
gen y2 = y + age - 40
* generating outcomes
* parameters chosen to give 90% power to detect a statistically significant (p<0.05) slope

gen time_rel_parental = age - parental_onset
label var time_rel_parental "time relative to parental onset"

gen time_rel_hybrid = time_rel_actual
replace time_rel_hybrid = time_rel_parental if time_rel_actual < 0
label var time_rel_hybrid "time relative to actual onset (if observed), otherwise parental"

gen hybrid_onset = actual_onset
replace hybrid_onset = parental_onset if time_rel_actual < 0
label var time_rel_hybrid "actual onset (if observed), otherwise parental"

gen actual_onset_observed=actual_onset
replace actual_onset_observed = . if actual_onset > age
* in fact actual onset time is unknown if after current age 
label var actual_onset_observed "actual onset (if observed)"
}
save D`i', replace
}
