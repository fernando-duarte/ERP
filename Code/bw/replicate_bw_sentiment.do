* This program takes the raw data in columns D through P in the DATA worksheet and produces the sentiment measures in columns B and C
version 15

global dir `c(pwd)'
global bw_dir = "${dir}/Input/Baker Wurgler"

import delimited using "${bw_dir}/bw_data_input.csv", delim(",") varnames(1) clear

gen year = int(yearmo/100) 
gen month = (yearmo-year*100)
gen mo = (year-1960)*12+(yearmo-year*100)

cap destring ripo, force replace

foreach varname of varlist indpro employ {
	gen g`varname' = `varname'/`varname'[_n-12]-1
	gen g`varname'1 = g`varname'[_n-12]
}
foreach varname of varlist consdur consnon consserv {
	gen g`varname' = (`varname'/`varname'[_n-12])/(cpi/cpi[_n-12])-1
	gen g`varname'1 = g`varname'[_n-12]
}
gen recess1 = recess[_n-1]

capture tsset mo
sort mo

* Define monthly equivalents to annual variables

gen nipo_sum = sum(nipo)
gen nipo_am = nipo_sum - nipo_sum[_n-12] if nipo~=.&nipo[_n-11]~=.
replace nipo_am = nipo_sum if nipo~=.&nipo[_n-11]~=.&nipo[_n-12]==.

replace ripo=0 if nipo==0 | ripo==.
gen ripo_am = .
sum mo
local lastmo = r(max)


forvalues t = 1/`lastmo' {
	qui sum nipo_am if mo == `t'
	local nipo_am = r(mean)
	capture reg ripo if (mo>`t'-12)&(mo<=`t') [aweight = nipo/`nipo_am']
	capture predict ripowtdavg if mo==`t'
	qui sum ripo nipo if (mo>`t'-12)&(mo<=`t') 
	capture replace ripo_am = ripowtdavg if mo==`t'&r(N)==12
	capture drop ripowtdavg
}

gen ripom = ripo
gen nipom = nipo
replace ripo = ripo_am
replace nipo = nipo_am

drop *_sum nipo_am ripo_am

sort mo
local 1 = 196507
local 2 = 201812
foreach varname of varlist pdnd ripo ripom {
	gen raw_lag`varname' = `varname'[_n-12]
	egen sraw_lag`varname' = std(raw_lag`varname') if yearmo>=`1' & yearmo<=`2' 
	reg raw_lag`varname' gindpro1 gconsdur1 gconsnon1 gconsserv1 gemploy1 recess1 if yearmo>=`1' & yearmo<=`2'
	predict e_lag`varname', resid
	egen se_lag`varname' = std(e_lag`varname') if yearmo>=`1' & yearmo<=`2'
}

foreach varname of varlist nipo cef s pdnd ripo nipom ripom {
	gen raw_`varname' = `varname'
	egen sraw_`varname' = std(`varname') if yearmo>=`1' & yearmo<=`2'
	reg raw_`varname' gindpro gconsdur gconsnon gconsserv gemploy recess if yearmo>=`1' & yearmo<=`2'
	predict e_`varname', resid
	egen se_`varname' = std(e_`varname') if yearmo>=`1' & yearmo<=`2'
}

foreach type in raw_ e_ {
	pca `type'cef `type'nipo `type'lagripo `type'lagpdnd `type's if yearmo>=`1' & yearmo<=`2'
	predict `type'f2 if yearmo>=`1' & yearmo<=`2'
	egen s`type'f2 = std(`type'f2)
	reg s`type'f2 s`type'cef s`type'nipo s`type'lagripo s`type'lagpdnd s`type's
}

list yearmo *f2
twoway line sraw_f2 se_f2 yearmo
gen ispos_sraw_f2_temp=(sraw_f2>0&yearmo==200012)
egen ispos_sraw_f2=max(ispos_sraw_f2_temp)
replace sraw_f2=-sraw_f2 if ispos_sraw_f2==0
gen ispos_se_f2_temp=(se_f2>0&yearmo==200012)
egen ispos_se_f2=max(ispos_se_f2_temp)
replace se_f2=-se_f2 if ispos_se_f2==0
twoway line sraw_f2 se_f2 yearmo
sum *f2

rename sraw_f2 SENT
rename se_f2 SENT_ORTH

order yearmo SENT_ORTH SENT pdnd ripom nipom cefd se sd s indpro consdur consnon consserv recess employ cpi  
keep yearmo SENT_ORTH SENT pdnd ripom nipom cefd se sd s indpro consdur consnon consserv recess employ cpi  

label var SENT "SENT"
label var SENT_ORTH "SENT_ORTH"

outsheet using "${bw_dir}/bw_data_final.csv", c replace


