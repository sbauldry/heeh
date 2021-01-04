*** Project: conduct education expansion and health analysis
*** Author: S Bauldry (based on R Frase)
*** Date: Dec 30, 2020


*** Load prepared data
cd ~/desktop
use heeh-hrs-data-mi, replace


*** generate predicted probabilities
mi est, saving(m1, replace): logit ba age i.rce i.fem i.fbr i.crur sib  ///
  i.clvg i.vet med fed i.cses i.cump i.mwrk csrh i.msch i.csmk i.psmk   ///
  if coh == 1
qui mi predict ba1 using m1
qui mi xeq: gen pba1 = invlogit(ba1)

mi est, saving(m2, replace): logit ba age i.rce i.fem i.fbr i.crur sib  ///
  i.clvg i.vet med fed i.cses i.cump i.mwrk csrh i.msch i.csmk i.psmk   ///
  if coh == 2
qui mi predict ba2 using m2
qui mi xeq: gen pba2 = invlogit(ba2)

  
*** simulation program
capture program drop SimDiD
program SimDiD
  args i j
  preserve
  
  * generate predicted degrees
  qui mi xeq: gen preba = rbinomial(1, pba1)
  qui mi xeq: gen pstba = rbinomial(1, pba2)
  
  * assign cases to NT, AT, and C
  qui mi xeq: gen     cat = 1 if ba == 0 & coh == 1 & pstba == 0 
  qui mi xeq: replace cat = 1 if ba == 0 & coh == 2
  qui mi xeq: replace cat = 2 if ba == 0 & coh == 1 & pstba == 1
  qui mi xeq: replace cat = 2 if ba == 1 & coh == 2 & preba == 0 
  qui mi xeq: replace cat = 3 if ba == 1 & coh == 1 
  qui mi xeq: replace cat = 3 if ba == 1 & coh == 2 & preba == 1

  qui mi xeq: gen nt = (cat == 1)
  qui mi xeq: gen cp = (cat == 2)
  qui mi xeq: gen at = (cat == 3) 
  
  foreach x of varlist nt cp at {
    qui sum `x' if `x' == 1 & coh == 1
	local `x'pr = r(N)
	
	qui sum `x' if `x' == 1 & coh == 2
	local `x'ps = r(N)
  }
  
  * fit primary DiD model and post results
  qui mi xeq: gen pstc = (coh == 2)
  foreach x of varlist dsrh dadl {
	
  	* full sample
    qui mi est (d: _b[1.cp#1.pstc] - _b[1.at#1.pstc]): reg `x' i.cp i.at     ///
      i.pstc i.cp#i.pstc i.at#i.pstc age i.rce i.fem i.fbr i.crur sib i.clvg ///
	  i.vet med fed i.cses i.cump i.mwrk csrh i.msch i.csmk i.psmk
    mat b1 = e(b_mi)
    mat b2 = e(b_Q_mi)
    mat v1 = e(V_mi)
    mat v2 = e(V_Q_mi)
  
    post pf (`i') (`j') ("all") ("`x'") (`ntpr') (`cppr') (`atpr') (`ntps') ///
	  (`cpps') (`atps') (b1[1,2]) (v1[2,2]) (b1[1,4]) (v1[4,4]) (b1[1,6])   ///
	  (v1[6,6]) (b1[1,10]) (v1[10,10]) (b1[1,14]) (v1[14,14]) (b2[1,1])     ///
	  (v2[1,1]) (b1[1,15]) (v1[15,15]) (b1[1,17]) (v1[17,17]) (b1[1,18])    ///
	  (v1[18,18]) (b1[1,20]) (v1[20,20]) (b1[1,22]) (v1[22,22]) (b1[1,24])  ///
	  (v1[24,24]) (b1[1,25]) (v1[25,25]) (b1[1,27]) (v1[27,27]) (b1[1,29])  ///
	  (v1[29,29]) (b1[1,30]) (v1[30,30]) (b1[1,31]) (v1[31,31]) (b1[1,33])  ///
	  (v1[33,33]) (b1[1,34]) (v1[34,34]) (b1[1,36]) (v1[36,36]) (b1[1,38])  ///
	  (v1[38,38]) (b1[1,39]) (v1[39,39]) (b1[1,40]) (v1[40,40]) (b1[1,42])  ///
	  (v1[42,42]) (b1[1,44]) (v1[44,44]) (b1[1,46]) (v1[46,46]) (b1[1,47])  ///
	  (v1[47,47])
  }  
	 
  
  
end


*** run simulation
forval i = 1/5 {
  postutil clear
  postfile pf sim1 sim2 str5 sam str5 y ntpr cppr atpr ntps cpps atps cp vcp ///
    at vat pc vpc cppc vcppc atpc vatpc cpat vcpat age vage blk vblk otr     ///
	votr fem vfem fbr vfbr rur vrur sib vsib lvg vlvg vet vvet med vmed fed  ///
	vfed ses1 vses1 ses2 vses2 ump vump mwk1 vmwk1 mwk2 vmwk2 csrh vcsrh     ///
	msc vmsc csk vcsk psk vpsk con vcon using heeh-hrs-anly-est-`i', replace
  
  forval j = 1/100 {
  	dis `i' `j'
  	SimDiD `i' `j'
  }
  
  postclose pf
}
