*** Project: education expansion and health in US
*** Author: S Bauldry (based on R Frase)
*** Date: Dec 29, 2020

*** Extracting analysis variables
cd ~/dropbox/research/data/hrs/hrs-old-data

*** Tracker file
use HHID PN BIRTHYR using tracker2016, replace
gen hhidpn = HHID + PN
tempfile d1
save `d1', replace

*** Rand longitudinal
use hhidpn raeduc raracem ragender rahispan racohbyr rabyear radyear ///
  r1agey_e-r13agey_e r1cenreg-r13cenreg rameduc rafeduc ravetrn rabplace ///
  r1shlt-r13shlt r2adla-r13adla r1psyche-r13psyche r1livsib-r13livsib ///
  hacohort using randhrs1992_2016v1, replace
tostring hhidpn, replace
tempfile d2
save `d2', replace

*** 1995 Core A data
use HHID PN D718 using A95A_R, replace
gen hhidpn = HHID + PN
tempfile d3
save `d3', replace

*** 1996 Core A data
use HHID PN E718 using H96A_R, replace
gen hhidpn = HHID + PN
tempfile d4
save `d4', replace

*** 1998 Core A data
use HHID PN F992 F993 F994 F995 F996 F997HM F1038 using H98A_R, replace
gen hhidpn = HHID + PN
tempfile d5
save `d5', replace

*** 2000 Core A data
use HHID PN G1079 G1080 G1081 G1082 G1083 G1084M G1085 G1125 using H00A_R, ///
  replace
gen hhidpn = HHID + PN
tempfile d6
save `d6', replace

*** 2002 Core B data
use HHID PN HB019 HB020 HB021 HB022 HB023 HB024M HB025 HB049 using H02B_R, ///
  replace
gen hhidpn = HHID + PN
tempfile d7
save `d7', replace

*** 2004 Core B data
use HHID PN JB019 JB020 JB021 JB022 JB023 JB024M JB025 JB049 using H04B_R, ///
  replace
gen hhidpn = HHID + PN
tempfile d8
save `d8', replace

*** 2006 Core B data
use HHID PN KB019 KB020 KB021 KB022 KB023 KB024M KB025 KB088 KB049 using ///
  H06B_R, replace
gen hhidpn = HHID + PN
tempfile d9
save `d9', replace

*** 2008 Core B data
use HHID PN LB019 LB099 LB104 LB020 LB021 LB117 LB022 LB023 LB025 LB088 ///
  LB049 LB122 LB123 using H08B_R, replace
gen hhidpn = HHID + PN
tempfile d10
save `d10', replace

*** 2010 Core B data
use HHID PN MB019 MB099 MB104 MB020 MB021 MB117 MB022 MB023 MB025 MB088 ///
  MB049 MB122 MB123 using H10B_R, replace
gen hhidpn = HHID + PN
tempfile d11
save `d11', replace

*** 2012 HRS Core B
use HHID PN NB019 NB099 NB104 NB020 NB021 NB117 NB022 NB023 NB025 NB088 ///
  NB049 NB122 NB123 using H12B_R, replace
gen hhidpn = HHID + PN
tempfile d12
save `d12', replace

*** 2014 HRS Core B
use HHID PN OB019 OB099 OB104 OB020 OB021 OB117 OB022 OB023 OB025 OB088 ///
  OB049 OB122 OB123 using H14B_R, replace
gen hhidpn = HHID + PN
tempfile d13
save `d13', replace

*** 2016 HRS Core B
use HHID PN PB019 PB099 PB104 PB020 PB021 PB117 PB022 PB023 PB025 PB088 ///
  PB049 PB122 PB123 using H16B_R, replace
gen hhidpn = HHID + PN
tempfile d14
save `d14', replace

*** Merge extracted data
use `d1', replace

forval i = 2/14 {
  merge 1:1 hhidpn using `d`i''
  drop if _merge == 2
  drop _merge
}



*** Preparing variables for analysis

* renaming HRS sample variable
rename hacohort smp

* pre- and post-expansion cohorts
rename BIRTHYR byr
gen     coh = 1 if byr >= 1916 & byr <= 1920
replace coh = 2 if byr >= 1957 & byr <= 1961
lab var coh "pre- and post-cohort"

* college degree
gen ba = (raeduc >= 5) if !mi(raeduc)
lab var ba "college degree or higher"

* health outcomes
rename radyear dyr

gen srh = r1shlt
gen adl = .
forval i = 2/13 {
  replace srh = 6 - r`i'shlt if !mi(r`i'shlt) & mi(srh)
  replace adl = r`i'adla if !mi(r`i'adla) & mi(adl)
}

gen dsrh = (srh >= 4) if !mi(srh)
gen dadl = (adl == 0) if !mi(adl)

lab var srh "self-reported health at first available wave"
lab var dsrh "self-reported health 4+ at first available wave"
lab var adl "ADL count at first available wave"
lab var dadl "no ADL at first available wave"

* demographic and other variables
rename racohbyr hrsc 

gen age = r1agey_e
gen crg = r1cenreg
gen sib = r1livsib
forval i = 2/13 {
  replace age = r`i'agey_e if mi(age) & !mi(r`i'agey_e)
  replace crg = r`i'cenreg if mi(crg) & !mi(r`i'cenreg)
  replace sib = r`i'livsib if mi(sib) & !mi(r`i'livsib)
}

rename raracem rce

gen fem = (ragender == 2) if !mi(ragender)
gen fbr = (rabplace == 11) if !mi(rabplace)
rename ravetrn vet

lab var age "age at first available wave"
lab var crg "census region at first available wave"
lab var sib "n siblings at first available wave"
lab var rce "race"
lab var fem "female"
lab var fbr "foreign born"
lab var vet "veteran"

* childhood variables
rename (rameduc rafeduc) (med fed)

gen csrh = .
foreach x of varlist F992 G1079 HB019 JB019 KB019 LB019 MB019 NB019 OB019 ///
  PB019 {
  recode `x' (8 9 = .)
  replace csrh = 6 - `x' if mi(csrh)
}

gen msch = .
foreach x of varlist LB099 MB099 NB099 OB099 PB099 {
  recode `x' (8 9 = .)
  replace msch = (`x' == 1) if mi(msch) & !mi(`x')
}

gen cses = .
foreach x of varlist F993 G1080 HB020 JB020 KB020 LB020 MB020 NB020 OB020 ///
  PB020 {
  recode `x' (6 8 9 = .) (1 = 3) (3 = 2) (5 = 1)
  replace cses = `x' if mi(cses)
}

gen cdrg = .
foreach x of varlist LB117 MB117 NB117 OB117 PB117 {
  recode `x' (8 9 = .)
  replace cdrg = (`x' == 1) if mi(cdrg) & !mi(`x')
}

gen cmvf = .
foreach x of varlist F994 G1081 HB021 JB021 KB021 LB021 MB021 NB021 OB021 ///
  PB021 {
  recode `x' (8 9 = .)
  replace cmvf = (`x' == 1) if mi(cmvf) & !mi(`x')
}

gen cfhp = .
foreach x of varlist F995 G1082 HB022 JB022 KB022 LB022 MB022 NB022 OB022 ///
  PB022 {
  recode `x' (8 9 = .)
  replace cfhp = (`x' == 1) if mi(cfhp) & !mi(`x')
}

gen cump = .
foreach x of varlist F996 G1083 HB023 JB023 KB023 LB023 MB023 NB023 OB023 ///
  PB023 {
  recode `x' (7 8 9 = .)
  replace cump = (`x' == 1 | `x' == 6) if mi(cump) & !mi(`x')
}

gen clvg = .
foreach x of varlist G1085 HB025 JB025 KB025 LB025 MB025 NB025 OB025 PB025 {
  recode `x' (8 9 = .)
  replace clvg = (`x' == 1) if mi(clvg) & !mi(`x')
}

gen mwrk = .
foreach x of varlist KB088 LB088 MB088 NB088 OB088 PB088 {
  recode `x' (7 8 9 = .) (5 = 0) (3 = 1) (1 = 2) 
  replace mwrk = `x' if mi(mwrk)
}

gen crur = .
foreach x of varlist D718 E718 F1038 G1125 HB049 JB049 KB049 LB049 MB049 ///
  NB049 OB049 PB049 {
  recode `x' (7 8 9 = .)
  replace crur = (`x' == 1) if mi(crur) & !mi(`x')
}

gen psmk = .
foreach x of varlist LB104 MB104 NB104 OB104 PB104 {
  recode `x' (8 9 = .)
  replace psmk = (`x' == 1 | `x' == 2) if mi(psmk) & !mi(`x')
}

gen csmk = .
foreach x of varlist LB122 MB122 NB122 OB122 PB122 {
  recode `x' (8 9 = .)
  replace csmk = (`x' == 1) if mi(csmk) & !mi(`x')
}

lab var csrh "childhood self-rated health"
lab var msch "childhood missed school due to health"
lab var cses "childhood SES"
lab var cdrg "childhood drug use"
lab var cmvf "childhood move for financial reasons"
lab var cfhp "childhood financial help"
lab var cump "childhood father unemployed"
lab var clvg "childhood lived with grandparents"
lab var mwrk "mother worked"
lab var crur "childhood rural"
lab var psmk "parent smoked"
lab var csmk "childhood smoked"


*** Keeping analysis variables
*** Note: set based on preliminary analyses building BA models
order smp coh ba srh dsrh adl dadl age rce fem fbr crur sib clvg vet med ///
  fed cses cump mwrk csrh msch csmk psmk
keep smp-psmk
  
*** Setting analysis sample and checking missing data
keep if !mi(coh)
dis _N
keep if !mi(ba)
keep if !mi(srh, adl)
keep if !mi(rce)

* checking distribution of missing data
misstable summarize age fem fbr crur sib clvg vet
misstable summarize med fed cses cump mwrk
misstable summarize csrh msch csmk psmk

*** Running multiple imputation
recode vet med fed rce (.r .d .m = .)
mi set wide
mi reg imp rce crur sib clvg vet med fed cses cump mwrk csrh msch csmk psmk
mi impute chain (ologit) csrh cses mwrk (logit) crur clvg vet cump msch psmk ///
  csmk (reg) med fed sib = age fem i.rce fbr ba srh dsrh adl dadl i.smp, ///
  by(coh) add(5) aug force dots rseed(135711) 

save ~/desktop/heeh-hrs-data-mi, replace
