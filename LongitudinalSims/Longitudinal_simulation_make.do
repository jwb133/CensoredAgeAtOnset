* cd "H:\DRG\AgeOnset\2023\Longitudinal 2024"
* log using "H:\DRG\AgeOnset\2023\Longitudinal 2024\Make_HPC_08_10_25.log", replace

* program drop mycubicspline
program define mycubicspline, rclass 
gen u1=t
gen u2=term1 
replace u2=u2-(term3*(knot4-knot1)-term4*(knot3-knot1))/(knot4-knot3)  
replace u2=u2/(knot4-knot1)^2

gen u3=term2 
replace u3=u3-(term3*(knot4-knot2)-term4*(knot3-knot2))/(knot4-knot3)  
replace u3=u3/(knot4-knot1)^2
end

version 18.5
set seed 973686515 
set sortseed 413688647

* forvalues i = 1/100 {	
forvalues i = 1/10000 {
qui{

clear 
set obs 50
gen mutation=_n
gen onset_re=rnormal(0,9*sqrt(0.6))
* onset random mutation  effect
gen intercept_re=rnormal(0,5)
* biomarker random mutation effect (on intecept)
summ *re

gen c=rpoisson(6)
expand c
* Mean of 6 individuals for each mutation in the older generation (generation=0)
* All of these 6 individuals are carriers
drop c
gen actual_onset = 46+onset_re+rnormal(0, 9*sqrt(0.4))
* Actual onset times in older generation
gen parental_onset = actual_onset
* For use later after generating the individuals in the younger generation

gen c=1+rpoisson(2)
summ c
expand c, gen(generation)
drop c
* Generating a mean of 2 children for each individual in the older generation

gen carrier=rbinomial(1,0.5)
* Individuals in the younger generation have a 50% chance of being a carrier 
replace carrier=1 if generation==0
* All individuals in the older generation are carriers
replace parental_onset=. if generation==0
* Individuals in the older generation have no parental data
table generation carrier

replace actual_onset = 46+onset_re + rnormal(0, 9*sqrt(0.4)) if generation==1 & carrier==1
* Individuals in the younger generation have onset times if carriers  
replace actual_onset = . if carrier==0
* Individuals in the younger generation have no onset time if non-carriers  

label var actual_onset "age at actual onset"
label var parental_onset "age at parental onset"

regress actual_onset parental_onset if generation==1 & carrier==1
* Checking slope = 0.6 (0.6/(0.6+0.4))

gen slope=rnormal(0,0.5) if generation==1 & carrier==1
* Biomarker random subject effect (on slope)
gen intercept=intercept_re+rnormal(0,5) if generation==1 & carrier==1
* Adding biomarker random subject effect (on intercept) to random mutation effect 

label var intercept "intercept"
label var slope "slope"
drop *_re
bysort generation carrier: summ

gen knot1=-20
gen knot2=-10
gen knot3=0
gen knot4=10
* Knot positions for restricted cubic splines

gen c=runiform()
* To be used to generate between 1 & 5 measures for younger generation
replace c=0 if generation==0
* Older generation individuals only have a single measurement
summ c

gen age1 = rnormal(38,5)
gen age2 = age1+1 if c>0.2
gen age3 = age1+2 if c>0.4
gen age4 = age1+3 if c>0.6
gen age5 = age1+4 if c>0.8
gen final_visit_age=max(age1,age2,age3,age4,age5)
gen id=_n
drop c

reshape long age, i(mutation id actual_onset parental_onset slope intercept final_visit_age) j(visit)
drop if age==.
replace age=. if generation==0
* Now have between 1 & 5 biomarker measurements for younger generation 

gen time_rel_actual = age - actual_onset if generation==1 & carrier==1  
label var time_rel_actual "time relative to actual onset"

gen max_time_rel_actual = final_visit_age - actual_onset if generation==1 & carrier==1

gen actual_onset_observed=actual_onset if generation==1 & carrier==1
replace actual_onset_observed = . if max_time_rel_actual < 0 & generation==1 & carrier==1
label var actual_onset_observed "actual onset (if observed)"

gen t=time_rel_actual
gen term1 = max((t-knot1)^3,0)
gen term2 = max((t-knot2)^3,0)
gen term3 = max((t-knot3)^3,0)
gen term4 = max((t-knot4)^3,0)

mycubicspline

gen y = intercept + slope*(time_rel_actual + 25) + 5*rnormal(0,1)
* Generating outcomes 
* Note that intercept relates to time 25 years before actual onset 
gen y1 = y + 100 + 1*u1 + 0*u2 + 0*u3 
gen y2 = y + 100 + 0*u1 + 0*u2 + 1.5*u3 
label var y "outcome"
* Generating two outcomes

rename u1 u1_actual 
rename u2 u2_actual 
rename u3 u3_actual 


gen time_rel_parental = age - parental_onset
label var time_rel_parental "time relative to parental onset"

replace t=time_rel_parental
replace term1 = max((t-knot1)^3,0)
replace term2 = max((t-knot2)^3,0)
replace term3 = max((t-knot3)^3,0)
replace term4 = max((t-knot4)^3,0)

mycubicspline

rename u1 u1_parental
rename u2 u2_parental
rename u3 u3_parental 


gen time_rel_hybrid = time_rel_actual
replace time_rel_hybrid = time_rel_parental if max_time_rel_actual< 0
label var time_rel_hybrid "time relative to actual onset (if observed), otherwise parental"

replace t=time_rel_hybrid
replace term1 = max((t-knot1)^3,0)
replace term2 = max((t-knot2)^3,0)
replace term3 = max((t-knot3)^3,0)
replace term4 = max((t-knot4)^3,0)

mycubicspline

rename u1 u1_hybrid
rename u2 u2_hybrid
rename u3 u3_hybrid 

}
save NewData`i', replace
display `i'
}

summ slope*