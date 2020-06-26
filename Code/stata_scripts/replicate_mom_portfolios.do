cd ..
cd ..

global dir `c(pwd)'
global indir = "${dir}/Input"

version 16

*****************************************************************************
* Clean CRSP data
*****************************************************************************
import delimited using "${indir}/crsp_m.txt", delim(" ") varnames(1) clear

tostring date, replace
gen date2 = date(date,"YMD")
replace date2 = mofd(date2)
format date2 %tmCCYYNN
drop date
rename date2 date
order date, after(permno)

replace ret = "0" if ret == ""
replace ret = "" if inlist(ret,"E","D","C","B","A")
destring ret, replace

*****************************************************************************
* Generate ME and lagged ME
*****************************************************************************
gen prc_adj = prc
replace prc_adj = altprc if mi(prc_adj)

* Calculate ME using prices adjusted by altprc
gen ME = abs(prc_adj) * shrout if prc_adj ~= 0 & shrout ~= 0
duplicates drop permno date prc, force
tsset permno date
gen l_ME = l.ME

*****************************************************************************
* Generate average returns for firms with multiple stocks
*****************************************************************************
bys permco date: egen sum_ME = sum(l_ME)
gen ret_vw = ret * l_ME / sum_ME if sum_ME ~= 0
bys permco date: egen avg_ret_vw = sum(ret_vw) if ~mi(ret_vw)	// value-weighted
replace ret = avg_ret_vw

*****************************************************************************
* Keep one observation for each permco-date pair
*****************************************************************************
duplicates drop permco date, force
tsset permco date
sort permco date

*****************************************************************************
* Calculate prior 2-12 returns
*****************************************************************************
gen ret_12_2 = (1 + L12.ret)
forval i = 2 / 11 {
	replace ret_12_2 = ret_12_2 * (1 + L`i'.ret) 
}
replace ret_12_2 = (ret_12_2 - 1) * 100

* replace with missing value if the criterion aren't met
replace ret_12_2 = . if mi(l13.prc) | mi(l2.ret) | mi(l_ME) | mi(prc)

*****************************************************************************
* Keep common stocks in NYSE, NASDAQ and AMEX
*****************************************************************************
* keep if stock is ordinary common and either "not further defined" or "need not be further defined" (as in Fama-French)
keep if shrcd==10 | shrcd==11

* keep if traded in NYSE, AMEX or Nasdaq only
keep if exchcd==1 | exchcd==2 | exchcd==3

*****************************************************************************
* Use French's breakpoints to classify portfolios
*****************************************************************************
preserve
import delimited "${indir}/Prior_2-12_Breakpoints.txt", ///
delimiter(space, collapse) varnames(nonames) rowrange(4) clear 

rename v2 count
forval i = 3 / 22 {
	local percentile = (`i' - 2) * 5
	rename v`i' prct`percentile'
}

tostring v1, replace
gen date = date(v1,"YM")
replace date = ym(year(date),month(date))
format date %tm
drop v1

save temp.tmp, replace

restore
merge m:1 date using temp.tmp, keep(1 3) nogen

sort permco date
gen p1 = 1 if ret_12_2 <= prct10 & ~mi(prct10)
forval i = 2 / 10 {
	local percentile = `i' * 10
	local percentile10 = (`i' - 1) * 10
	gen p`i' = 1 if ret_12_2 > prct`percentile10' & ret_12_2 <= prct`percentile'
}

replace ret = ret * 100

*****************************************************************************
* Generate equal- and value-weighted portfolio returns
*****************************************************************************
* Equal-weighted
cap rm "${indir}/10_Portfolios_Prior_12_2_Equal_weighted.dta"
forval i = 1 / 10 {
	preserve
	collapse (mean) ret if p`i' == 1, by(date) cw
	rename ret ret_`i'_rep
	cap merge 1:1 date using "${indir}/10_Portfolios_Prior_12_2_Equal_weighted", nogen
	order _all, sequential
	save "${indir}/10_Portfolios_Prior_12_2_Equal_weighted", replace
	restore
}

* Value-weighted
cap rm "${indir}/10_Portfolios_Prior_12_2_Value_weighted.dta"
forval i = 1 / 10 {
	preserve
	collapse (mean) ret if p`i' == 1 [w = l_ME], by(date) cw
	rename ret ret_`i'_rep
	cap merge 1:1 date using "${indir}/10_Portfolios_Prior_12_2_Value_weighted", nogen
	order _all, sequential
	save "${indir}/10_Portfolios_Prior_12_2_Value_weighted", replace
	restore
}

use "${indir}/10_Portfolios_Prior_12_2_Equal_weighted", clear
outsheet using "${indir}/10_Portfolios_Prior_12_2_Equal_weighted.txt", replace

use "${indir}/10_Portfolios_Prior_12_2_Value_weighted", clear
outsheet using "${indir}/10_Portfolios_Prior_12_2_Value_weighted.txt", replace
rm temp.tmp
