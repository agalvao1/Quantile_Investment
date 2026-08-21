/*
Date Created:	20250113
Date Modified: 	20250113
Quotes from:	CheckWRDS2010DifferenceUsingcgj19table1420220613.do
Data:			WRDS archive 5
Updates:		Restrict obs to USA firms
				All outputs are in latex
*/

clear all
log close _all

*ssc install asdoc
*ssc install estout

*local suffix: display %tdCCYYNNDD =daily("`c(current_date)'", "DMY")

*log using "C:\Users\cuijing3\OneDrive - Michigan State University\Research\Summer_RA_task\Log\log_`suffix'.txt", text append

cd "/Users/WRDS_data"


*****************************************************************************************
***********************************WRDS part I data cleaning*****************************
*****************************************************************************************
*****************************************************************************************

use "/Users/WRDS_data/AnnualCPI25.dta", replace
ren year fyear
g cpi1971 = cpi[26]
g adjcpi = cpi/cpi1971
save temp1, replace

use "/Users/WRDS_data/d0vcffx1yhjs1mhc.dta",replace

gen year = year(datadate)
gen month = month(datadate)
gen adjust = (month < 6)
gen fyear = year - adjust

*gen fyear=fyr+1970
merg m:1 fyear using temp1

keep if _mer==3
drop _merge cpi*

ren sic dnum
destring dnum, replace
destring gvkey, replace
sort fyear 


***Industry-Time Filters***
**Use only manufacturers
keep if dnum>1999 & dnum<4000

**Adjusting all data for CPI
quiet for @ in var act ib dp dvp dvc capx at ch che csho ceq sale dlc dltt ni prstkc ppegt ppent prcc_f txdb xint pstk invt: quiet replace @ = @/adjcpi


compress


********************
***obs of interest**
********************
* a brief reminder for additional steps:
* Need to remove indfmt = FS or tsset gvkey year suffers from repeated time values within panel, n = 19, gvkey == 5860 | gvkey ==163678 |gvkey == 166368
*keep if fic == "USA"

** new code start: remove duplicates and use the correct time variable
duplicates report gvkey fyear
bysort gvkey fyear: gen dup = _N
preserve
drop if dup == 1
list
restore

drop if dup == 2
* Or, for future codes, you might want to use something like the following which drop the obs in spring/summer:
* drop if dup == 2 & month < 7

* to avoid confusion, remove the temp variables and rename year
drop year month adjust dup
rename fyear year
** new code end --20260117

gen SIC3 = int(dnum/10)
tsset gvkey year

drop if ppent==.
drop if ppent<5 
*drop if at<1
drop if year< 1969

*The following two lines are not mentioned in the paper but is included in cgj8.do
g KGrowth = ppent/L.ppent
drop if KGrowth>2 & KGrowth!=.

g AssetGrowth = at/L.at
drop if AssetGrowth>2 & AssetGrowth!=.
g SalesGrowth = sale/L.sale
drop if SalesGrowth>2 & SalesGrowth!=. 
drop AssetGrowth KGrowth SalesGrowth
sort gvkey year

*****VARIABLE DEFINITIONS*****
tsset gvkey year

g KexpK1 = capx/L.ppent 
label variable KexpK1 "Investment1"

g KexpK2 = at/L.at 
label variable KexpK2 "Investment2"

g PPE_A = ppent/at

g CF = ib+dp
*g CF = ib+dp-dvp-dvc
g CF_K = CF/L.ppent 
label variable CF_K "Cash Flow"

g Q1 = (at+(prcc_f*csho)-ceq-txdb)/at
label variable Q1 "Tobins Q1(WP)"

g E = ceq + pstk
label variable E "Equity"

gen D = dlc + dltt
label variable D "Total Debt"


gen Ch = ch/at
label variable D "Cash"

gen Che = che/at
label variable D "Cash"

*gen Q2 = (E+D-invt)/L.ppent
*label variable Q2 "Tobins Q2"
*not used: the book value of assets minus the book value of equity plus the marketvalue of equity divided by the book value of assets Smith and Watts 1992 (Hubbard and Palia 1999)


gen Q2 = (prcc_f*csho + pstk + D)/at 
label variable Q2 "Tobins Q2(Daines)"


**Dropping missings
for @ in var Q1 Q2 KexpK1 KexpK2 CF_K: drop if @==.


**Some Sorting Variables
g payout = (dvp+dvc+prstkc)/ni
drop if payout==.
*g TotalLeverage = (dltt/*+data34*/)/at
g assets = at
*g bondrates = splticrm
*g cprates = spsticrm

**************************************
***Financial Constraints Categories***
**************************************

****Size Categories (Book Values)****
** Categorize firms year-by-year
egen Pct40 = pctile(at), p(40) by (year)
egen Pct60 = pctile(at), p(60) by (year)
g Small = (at<=Pct40)
g Large = (at>=Pct60)
drop Pct40 Pct60


****Dividend payout ratio****
**Categorize firms year-by-year
egen Pct40 = pctile(payout), p(40) by (year)
egen Pct60 = pctile(payout), p(60) by (year)
g LowPayout = (payout<=Pct40)
g HighPayout = (payout>=Pct60)
drop Pct40 Pct60


****Debt****
**Categorize firms year-by-year
egen Pct40 = pctile(D), p(40) by (year)
egen Pct60 = pctile(D), p(60) by (year)
g LowDebt = (D<=Pct40)
g HighDebt = (D>=Pct60)
drop Pct40 Pct60


****Cash****
**Categorize firms year-by-year
egen Pct40 = pctile(Che), p(40) by (year)
egen Pct60 = pctile(Che), p(60) by (year)
g LowCash = (Che<=Pct40)
g HighCash = (Che>=Pct60)
drop Pct40 Pct60



preserve
label variable ppent "Capital Stock"
label variable at "Asset"

summarize KexpK1 KexpK2 Q1 Q2 ppent at, detail

egen idaux=group(gvkey)
bysort idaux : keep if _N >= 5
drop idaux

egen id=group(gvkey)

save WRDS_clean_constrained_2025n2p40, replace

*use WRDS_clean, clear

