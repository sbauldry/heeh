*** Project: conduct education expansion and health analysis
*** Author: S Bauldry (based on R Frase)
*** Date: Dec 30, 2020


*** Load prepared data for descriptive statistics
cd ~/desktop
use heeh-hrs-data-mi, replace

*** Table 1: descriptive statistics
preserve
qui mi xeq: gen rce1 = (rce == 1)
qui mi xeq: gen rce2 = (rce == 2)
qui mi xeq: gen rce3 = (rce == 3)

qui mi xeq: gen cses1 = (cses == 1)
qui mi xeq: gen cses2 = (cses == 2)
qui mi xeq: gen cses3 = (cses == 3)

qui mi xeq: gen mwrk0 = (mwrk == 0)
qui mi xeq: gen mwrk1 = (mwrk == 1)
qui mi xeq: gen mwrk2 = (mwrk == 2)

* program to estimate mean, difference in means, and p-value
capture program drop emd
program emd, rclass
  args x c
  
  qui mi est: mean `x' if coh == `c' & ba == 0
  mat b0 = e(b_mi)
  mat v0 = e(V_mi)

  qui mi est: mean `x' if coh == `c' & ba == 1
  mat b1 = e(b_mi)
  mat v1 = e(V_mi)
  
  local m0  = b0[1,1]
  local m1  = b1[1,1]
  local dm  = b1[1,1] - b0[1,1]
  local dse = sqrt(v1[1,1] + v0[1,1])
  local pv  = 2*(1 - normal(abs(`dm'/`dse')))
  
  return scalar m0 = `m0'
  return scalar m1 = `m1'
  return scalar pv = `pv'

end

tab coh ba

* pre-cohort
foreach x of varlist dsrh dadl age rce1-rce3 fem fbr crur sib clvg vet med ///
  fed cses1-cses3 cump mwrk0-mwrk2 csrh msch csmk psmk {
  emd `x' 1
  dis "`x'," as res %5.2f r(m0) "," as res %5.2f r(m1) "," as res %5.3f r(pv)
}

* post-cohort
foreach x of varlist dsrh dadl age rce1-rce3 fem fbr crur sib clvg vet med ///
  fed cses1-cses3 cump mwrk0-mwrk2 csrh msch csmk psmk {
  emd `x' 2
  dis "`x'," as res %5.2f r(m0) "," as res %5.2f r(m1) "," as res %5.3f r(pv)
}
restore


*** Table A1: propensity score model estimates
eststo clear
mi est, post saving(m1, replace): logit ba age i.rce i.fem i.fbr i.crur sib ///
  i.clvg i.vet med fed i.cses i.cump i.mwrk csrh i.msch i.csmk i.psmk       ///
  if coh == 1
eststo pm1

mi est, post saving(m2, replace): logit ba age i.rce i.fem i.fbr i.crur sib ///
  i.clvg i.vet med fed i.cses i.cump i.mwrk csrh i.msch i.csmk i.psmk       ///
  if coh == 2
eststo pm2

esttab pm1 pm2, b(%5.2f) se(%5.2f) nogap wide nostar one


*** Diagnostic: check VIFs for covariates
*** no estimates above 2.5 in either cohort
preserve
mi convert flong
qui reg ba age i.rce i.fem i.fbr i.crur sib i.clvg i.vet med fed i.cses ///
  i.cump i.mwrk csrh i.msch i.csmk i.psmk if coh == 1 & _mi_m == 1
vif

qui reg ba age i.rce i.fem i.fbr i.crur sib i.clvg i.vet med fed i.cses ///
  i.cump i.mwrk csrh i.msch i.csmk i.psmk if coh == 2 & _mi_m == 1
vif
restore



*** Figure X1: distribution of NE, CE, and AE in primary analysis
use heeh-hrs-anly-est-1, replace
forval i = 2/5 {
  append using heeh-hrs-anly-est-`i'
}

collapse (mean) ntpr-atps, by(sam)
keep if _n == 1

sum 

foreach x in nt cp at {
  gen p`x'pr = `x'pr/2467
  gen p`x'ps = `x'ps/4189
}

rename (pntpr pcppr patpr pntps pcpps patps) (p1 p2 p3 p4 p5 p6)
gen id = 1
order id p1 p2 p3 p4 p5 p6
keep id-p6

reshape long p, i(id) j(grp)
gen coh = 1 if _n < 4
replace coh = 2 if _n >= 4
recode grp (4 = 1) (5 = 2) (6 = 3)

lab def g 1 "NE" 2 "CE" 3 "AE"
lab val grp g

lab def c 1 "pre-expansion" 2 "post-expansion"
lab val coh c

format p %5.2f

tempfile g1 g2
graph bar p if coh == 1, over(grp) scheme(s1color) ylab(0(0.2)1, angle(h) grid       ///
  gstyle(dot)) ytit("proportion respondents") blabel(total, format(%5.2g)) ///
  tit("pre-expansion") saving(`g1', replace)
  
graph bar p if coh == 2, over(grp) scheme(s1color) ylab(0(0.2)1, angle(h) grid       ///
  gstyle(dot)) ytit("proportion respondents") blabel(total, format(%5.2g)) ///
  tit("post-expansion") saving(`g2', replace)
  
graph combine "`g1'" "`g2'", scheme(s1color) rows(1)
graph export ~/desktop/heeh-figX1.pdf, replace


*** Figure X2: estimates from primary analysis
use heeh-hrs-anly-est-1, replace
forval i = 2/5 {
  append using heeh-hrs-anly-est-`i'
}

collapse (mean) cp-vcpat (sd) sd_vcp = vcp sd_vat = vat sd_vpc = vpc ///
  sd_vcppc = vcppc sd_vatpc = vatpc sd_vcpat = vcpat, by(y)

foreach x in cp at pc cppc atpc cpat {
  gen se_`x' = sqrt( v`x' + (1 + 1/500)*(sd_v`x')^2 )
  gen ub_`x' = `x' + 1.96*se_`x'
  gen lb_`x' = `x' - 1.96*se_`x'
}

gen id = _n

rename (cp at pc cppc atpc cpat) (es1 es2 es3 es4 es5 es6)
rename (ub_cp ub_at ub_pc ub_cppc ub_atpc ub_cpat) (ub1 ub2 ub3 ub4 ub5 ub6)
rename (lb_cp lb_at lb_pc lb_cppc lb_atpc lb_cpat) (lb1 lb2 lb3 lb4 lb5 lb6)

keep id es* ub* lb*
reshape long es ub lb, i(id) j(vr)

replace vr = 7 - vr

tempfile g1 g2
graph twoway (scatter vr es, mc(black)) (rspike ub lb vr, hor lc(blue)) ///
  if id == 1, scheme(s1color) ylab(6 "CE" 5 "AE" 4 "post" 3 "CE x post" ///
  2 "AE x post" 1 "CE vs AE", angle(h) grid gstyle(dot)) ytit("")       ///
  xtit("estimate") xlab(-0.4(0.2)0.4, grid gstyle(dot)) xline(0,        ///
  lc(black) lp(dash)) legend(off) tit("No ADL") saving(`g1')
  
graph twoway (scatter vr es, mc(black)) (rspike ub lb vr, hor lc(blue)) ///
  if id == 2, scheme(s1color) ylab(6 "CE" 5 "AE" 4 "post" 3 "CE x post" ///
  2 "AE x post" 1 "CE vs AE", angle(h) grid gstyle(dot)) ytit("")       ///
  xtit("estimate") xlab(-0.4(0.2)0.4, grid gstyle(dot)) xline(0,        ///
  lc(black) lp(dash)) legend(off) tit("VG/E SRH") saving(`g2')
  
graph combine "`g1'" "`g2'", scheme(s1color)
graph export ~/desktop/heeh-figX2.pdf, replace 


*** Table A2: estimates from primary regression models
use heeh-hrs-anly-est-1, replace
forval i = 2/5 {
  append using heeh-hrs-anly-est-`i'
}

collapse (mean) cp-vcon (sd) sd_vcp = vcp sd_vat = vat sd_vpc = vpc ///
  sd_vcppc = vcppc sd_vatpc = vatpc sd_vcpat = vcpat sd_vage = vage ///
  sd_vblk = vblk sd_votr = votr sd_vfem = vfem sd_vfbr = vfbr       ///
  sd_vrur = vrur sd_vsib = vsib sd_vlvg = vlvg sd_vvet = vvet       ///
  sd_vmed = vmed sd_vfed = vfed sd_vses1 = vses1 sd_vses2 = vses2   ///
  sd_vump = vump sd_vmwk1 = vmwk1 sd_vmwk2 = vmwk2 sd_vcsrh = vcsrh ///
  sd_vmsc = vmsc sd_vcsk = vcsk sd_vpsk = vpsk sd_vcon = vcon, by(y)

foreach x in cp at pc cppc atpc cpat age blk otr fem fbr rur sib lvg vet ///
  med fed ses1 ses2 ump mwk1 mwk2 csrh msc csk psk con {
  gen se_`x' = sqrt( v`x' + (1 + 1/500)*(sd_v`x')^2 )
}

gen id = _n

rename (cp at pc cppc atpc cpat age blk otr fem fbr rur sib lvg vet med fed ///
  ses1 ses2 ump mwk1 mwk2 csrh msc csk psk con) (e1 e2 e3 e4 e5 e6 e7 e8 e9 ///
  e10 e11 e12 e13 e14 e15 e16 e17 e18 e19 e20 e21 e22 e23 e24 e25 e26 e27)
  
rename (se_cp se_at se_pc se_cppc se_atpc se_cpat se_age se_blk se_otr      ///
  se_fem se_fbr se_rur se_sib se_lvg se_vet se_med se_fed se_ses1 se_ses2   ///
  se_ump se_mwk1 se_mwk2 se_csrh se_msc se_csk se_psk se_con) (se1 se2 se3  ///
  se4 se5 se6 se7 se8 se9 se10 se11 se12 se13 se14 se15 se16 se17 se18 se19 ///
  se20 se21 se22 se23 se24 se25 se26 se27)
  
keep id e* se*

reshape long e se, i(id) j(vr)

list e se if id == 1, clean noobs table
list e se if id == 2, clean noobs table

