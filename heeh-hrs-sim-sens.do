*** Project: conduct education expansion and health analysis
*** Author: S Bauldry (based on R Frase)
*** Date: Dec 30, 2020


*** Load prepared data
cd ~/desktop
use heeh-hrs-data-mi, replace


*** generate predicted probabilities
qui mi est, saving(m1, replace): logit ba age i.rce i.fem i.fbr i.crur sib  ///
  i.clvg i.vet med fed i.cses i.cump i.mwrk csrh i.msch i.csmk i.psmk       ///
  if coh == 1
qui mi predict ba1 using m1
qui mi xeq: gen pba1 = invlogit(ba1)

qui mi est, saving(m2, replace): logit ba age i.rce i.fem i.fbr i.crur sib  ///
  i.clvg i.vet med fed i.cses i.cump i.mwrk csrh i.msch i.csmk i.psmk       ///
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
  qui mi xeq: gen pstc = (coh == 2)
  
  * sensitivity: fit DiD models across subgroups and post results
  foreach x of varlist dsrh dadl {
  	
	* women
	qui mi est (d: _b[1.cp#1.pstc] - _b[1.at#1.pstc]): reg `x' i.cp i.at   ///
      i.pstc i.cp#i.pstc i.at#i.pstc age i.rce i.fbr i.crur sib i.clvg     ///
	  i.vet med fed i.cses i.cump i.mwrk csrh i.msch i.csmk i.psmk if fem
    mat b1 = e(b_mi)
    mat b2 = e(b_Q_mi)
    mat v1 = e(V_mi)
    mat v2 = e(V_Q_mi)
  
    post pf (`i') (`j') ("fem") ("`x'") (b1[1,2]) (v1[2,2]) (b1[1,4]) ///
	  (v1[4,4]) (b1[1,6]) (v1[6,6]) (b1[1,10]) (v1[10,10]) (b1[1,14]) ///
	  (v1[14,14]) (b2[1,1]) (v2[1,1])
	  
	* men
	qui mi est (d: _b[1.cp#1.pstc] - _b[1.at#1.pstc]): reg `x' i.cp i.at   ///
      i.pstc i.cp#i.pstc i.at#i.pstc age i.rce i.fbr i.crur sib i.clvg     ///
	  i.vet med fed i.cses i.cump i.mwrk csrh i.msch i.csmk i.psmk if !fem
    mat b1 = e(b_mi)
    mat b2 = e(b_Q_mi)
    mat v1 = e(V_mi)
    mat v2 = e(V_Q_mi)
  
    post pf (`i') (`j') ("mal") ("`x'") (b1[1,2]) (v1[2,2]) (b1[1,4]) ///
	  (v1[4,4]) (b1[1,6]) (v1[6,6]) (b1[1,10]) (v1[10,10]) (b1[1,14]) ///
	  (v1[14,14]) (b2[1,1]) (v2[1,1])
	  
    * white
	qui mi est (d: _b[1.cp#1.pstc] - _b[1.at#1.pstc]): reg `x' i.cp i.at   ///
      i.pstc i.cp#i.pstc i.at#i.pstc age i.fem i.fbr i.crur sib i.clvg     ///
	  i.vet med fed i.cses i.cump i.mwrk csrh i.msch i.csmk i.psmk if rce == 1
    mat b1 = e(b_mi)
    mat b2 = e(b_Q_mi)
    mat v1 = e(V_mi)
    mat v2 = e(V_Q_mi)
  
    post pf (`i') (`j') ("wht") ("`x'") (b1[1,2]) (v1[2,2]) (b1[1,4]) ///
	  (v1[4,4]) (b1[1,6]) (v1[6,6]) (b1[1,10]) (v1[10,10]) (b1[1,14]) ///
	  (v1[14,14]) (b2[1,1]) (v2[1,1])
	  
	* black
	qui mi est (d: _b[1.cp#1.pstc] - _b[1.at#1.pstc]): reg `x' i.cp i.at   ///
      i.pstc i.cp#i.pstc i.at#i.pstc age i.fem i.fbr i.crur sib i.clvg     ///
	  i.vet med fed i.cses i.cump i.mwrk csrh i.msch i.csmk i.psmk if rce == 2
    mat b1 = e(b_mi)
    mat b2 = e(b_Q_mi)
    mat v1 = e(V_mi)
    mat v2 = e(V_Q_mi)
  
    post pf (`i') (`j') ("blk") ("`x'") (b1[1,2]) (v1[2,2]) (b1[1,4]) ///
	  (v1[4,4]) (b1[1,6]) (v1[6,6]) (b1[1,10]) (v1[10,10]) (b1[1,14]) ///
	  (v1[14,14]) (b2[1,1]) (v2[1,1])
  }
  
  * sensitivity: fit DiD models on continuous measures
  foreach x of varlist srh adl {
	
  	* full sample / regress
    qui mi est (d: _b[1.cp#1.pstc] - _b[1.at#1.pstc]): reg `x' i.cp i.at     ///
      i.pstc i.cp#i.pstc i.at#i.pstc age i.rce i.fem i.fbr i.crur sib i.clvg ///
	  i.vet med fed i.cses i.cump i.mwrk csrh i.msch i.csmk i.psmk
    mat b1 = e(b_mi)
    mat b2 = e(b_Q_mi)
    mat v1 = e(V_mi)
    mat v2 = e(V_Q_mi)
  
    post pf (`i') (`j') ("all") ("`x'r") (b1[1,2]) (v1[2,2]) (b1[1,4]) ///
	  (v1[4,4]) (b1[1,6]) (v1[6,6]) (b1[1,10]) (v1[10,10]) (b1[1,14])  ///
	  (v1[14,14]) (b2[1,1]) (v2[1,1])
  } 
  
  * sensitivity: full sample / ologit  
  capture qui mi est (d: _b[1.cp#1.pstc] - _b[1.at#1.pstc]): ologit srh     ///
	  i.cp i.at i.pstc i.cp#i.pstc i.at#i.pstc age i.rce i.fem i.fbr i.crur ///
	  sib i.clvg i.vet med fed i.cses i.cump i.mwrk csrh i.msch i.csmk i.psmk
    mat b1 = e(b_mi)
    mat b2 = e(b_Q_mi)
    mat v1 = e(V_mi)
    mat v2 = e(V_Q_mi)
  
    post pf (`i') (`j') ("all") ("srho") (b1[1,2]) (v1[2,2]) (b1[1,4]) ///
	  (v1[4,4]) (b1[1,6]) (v1[6,6]) (b1[1,10]) (v1[10,10]) (b1[1,14])  ///
	  (v1[14,14]) (b2[1,1]) (v2[1,1])
	  
  * sensitivity: full sample / count  
  capture qui mi est (d: _b[1.cp#1.pstc] - _b[1.at#1.pstc]): nbreg adl      ///
	  i.cp i.at i.pstc i.cp#i.pstc i.at#i.pstc age i.rce i.fem i.fbr i.crur ///
	  sib i.clvg i.vet med fed i.cses i.cump i.mwrk csrh i.msch i.csmk i.psmk
    mat b1 = e(b_mi)
    mat b2 = e(b_Q_mi)
    mat v1 = e(V_mi)
    mat v2 = e(V_Q_mi)
  
    post pf (`i') (`j') ("all") ("adln") (b1[1,2]) (v1[2,2]) (b1[1,4]) ///
	  (v1[4,4]) (b1[1,6]) (v1[6,6]) (b1[1,10]) (v1[10,10]) (b1[1,14])  ///
	  (v1[14,14]) (b2[1,1]) (v2[1,1])
  
end


*** run simulation
forval i = 1/5 {
  postutil clear
  postfile pf sim1 sim2 str5 sam str5 y cp vcp at vat pc vpc cppc vcppc atpc ///
    vatpc cpat vcpat using heeh-hrs-anly-sens-est-`i', replace
  
  forval j = 1/100 {
  	dis `i' `j'
  	SimDiD `i' `j'
  }
  
  postclose pf
}
