*! v.1.0.0 12Sept2014 Orsini N. 

capture program drop laplace_surv
program define laplace_surv
version 9.2
syntax [,  AT1(string) /*
		*/ AT2(string) AT3(string) AT4(string) AT5(string) /*
		*/ AT6(string) AT7(string) AT8(string) AT9(string) /*
		*/ AT10(string) /*
		*/ Range(numlist ascending min=2 max=2) Format(string)  /*
		*/ GENerate(namelist max=11) scatter line lowess LOptions(string) SHow LPred * ]

if "`e(cmd)'"!="laplacereg" {
	error 301
}

tempvar perc surv 
quietly gen `perc' = 0 in 1
local j = 2
local nq : word count `e(qlist)'
local currentn = `c(N)'
if (`currentn' < (`nq'+1)) qui set obs `=`nq'+1' 

if "`format'" == "" 	local fmt = "%2.1fc"
	else local fmt = "`format'" 
	
_get_gropts , graphopts(`options')
local options `"`s(graphopts)'"'
	
if `"`range'"'~="" {
		tokenize `range'
		local ranopt  "`1', `2'"
	}
	else {    
			qui summ `e(depvar)' if e(sample), meanonly
			local ranopt  "0, `r(max)'"
		}
		
quietly foreach p of any `e(qlist)' {
	replace `perc' = `p' in `j++'
}

local eqnames `e(eqnames)'
local fq "`: word 1 of `eqnames''"
tempname b c
mat `b'= get(_b)
mat `c'= `b'[1, "`fq':"]
local rhsorig: colnames(`c')					
local rhsorig: subinstr local rhsorig "_cons" "", all
local countat 0
				
if `"`at1'`at2'`at3'`at4'`at5'`at6'`at7'`at8'`at9'`at10'"' ~= "" {
		local i 1
		while `i'<=10 {			
		 if "`at`i''"~="" {
				local countat = `countat' + 1

				tempvar pred`i'
				qui gen `pred`i'' = .
                local p = 2
			    if "`lpred'" != ""  display _n "at`i'(`at`i'')"

								 
			foreach q of local eqnames {	
				local equation`i' 
				local fixedvar 
				
				tokenize "`at`i''" , parse(" =")
					local j = 1
					local k = 1
				
					while "``j''" != "" {
					if `k' > 3 local k = 1
						if (`k'==1) {
									 unab namev: ``j''
									 local fixedvar `fixedvar' `namev'
									 local equation`i' `equation`i'' [`q']`namev'
									 }
						if (`k'==3) local equation`i' `equation`i''*``j''+
						local `k++'
						local j = `j'+1
					}
				local restvar: subinstr local rhsorig "`fixedvar'" "", all
				foreach x of local restvar {
						qui su `x'  if e(sample), meanonly
						local equation`i' `equation`i'' [`q']`x'*`=r(mean)' +
					}
			    local equation`i' `equation`i'' [`q']_cons
				if "`lpred'" != ""  display "`equation`i''"
				qui lincom  `equation`i''
				qui replace `pred`i'' = r(estimate) in `p++'
			}
		}
		local i = `i' +1 
		}
}		

if `"`at1'`at2'`at3'`at4'`at5'`at6'`at7'`at8'`at9'`at10'"' == "" {

		tempvar pred
		qui gen `pred' = .
        local p = 2
				
			foreach q of local eqnames {
            local equation 
				if `"`rhsorig'"' != "" {
				   foreach x of local rhsorig {
						qui su `x'  if e(sample), meanonly
						local equation `equation' [`q']`x'*`=r(mean)' +
					}
				}
				
			    local equation `equation' [`q']_cons
				if "`lpred'" != "" display "`equation`i''"

				qui lincom  `equation'
				qui replace `pred' = r(estimate) in `p++'
			}			
}

qui replace `perc' = `perc'*100
qui gen `surv' = (100-`perc')
char  `surv'[varname] "S(t)"

if `"`at1'`at2'`at3'`at4'`at5'`at6'`at7'`at8'`at9'`at10'"' != "" {
local patterns "solid solid dot dash_dot shortdash shortdash_dot longdash longdash_dot dash_3dot dash_dot_dot"          
	tokenize `patterns'
	forv i = 1/`countat' {			
		char  `pred`i''[varname] "t at`i'"
		local leglabels `leglabels' label(`i' at`i') 
		qui replace `pred`i'' = 0 in 1
		local predlist `predlist' `pred`i''		
		   if `"`line'"' != ""     local plot`i' "(line `surv' `pred`i'' if inrange(`pred`i'', `ranopt'), lp(``i'')  lw(med) connect(stairstep))"
		   if `"`scatter'"' != ""  local plot`i' "(scatter `surv' `pred`i'' if inrange(`pred`i'', `ranopt'), lp(``i'')  lw(med)   c(l) ms(o) mc(black) msize(*2))"
		   if `"`lowess'"' != ""   local plot`i' "(lowess `surv' `pred`i'' if inrange(`pred`i'', `ranopt'), lp(``i'')  lw(med)   `loptions')"
		  local graph `graph' `plot`i''
	}
	
if "`show'" != "" {
					set more off
					format `fmt' `predlist'
					list `surv' `predlist' if `surv' != . , noobs clean subvarname
}
			
}
else {
		char  `pred'[varname] "t"
		qui replace `pred' = 0 in 1
		   if `"`line'"' != ""     local plot "(line `surv' `pred' if inrange(`pred', `ranopt'), lp(l)  lw(med) connect(stairstep)  )"
		   if `"`scatter'"' != ""  local plot "(scatter `surv' `pred' if inrange(`pred', `ranopt'), lp(l)  lw(med)   c(l) ms(o) mc(black) msize(*2))"
		   if `"`lowess'"' != ""   local plot "(lowess `surv' `pred' if inrange(`pred', `ranopt'), lp(l)  lw(med)   `loptions')"
			local graph `plot'
			local predlist `pred'		
			if "`show'" != "" {
					set more off
					format `fmt' `pred'
					list `surv' `pred' if `surv' != . , noobs clean subvarname
			}
}

// Save new variables containing the displayed results

if "`generate'" != "" {
		local listvarnames "`surv' `predlist'" 
		local nnv : word count `generate' 
		tokenize `generate'
		forv i = 1/`nnv' {	
				qui gen ``i'' = `: word `i' of `listvarnames''
		}
}

if "`graph'" != "" {

if "`: variable label `e(depvar)''"	== "" local xtitle "`e(depvar)'"
else local xtitle "`: variable label `e(depvar)''"

twoway `graph',  ///  
	ytitle("Survival Proportion (%)") ///
	xtitle("`xtitle'") ///
	ylabel(#10, angle(horiz)  )  ///
	xlabel(#10)  ///
	legend(`leglabels'  col(1) ring(0) pos(1) region(style(none))) ///
	scheme(s1rcolor) /// 
	plotregion(style(none)) ///
	`options' 
}

* back to original dataset if number of quantiles is great
if "`generate'" == "" {
if (`c(N)' > `currentn') qui keep in 1/`currentn'   
}
end

exit

sysuse cancer, clear
set more off
qui xi: laplacereg studytime i.drug age, q(20(10)50)

laplace_surv, scatter
				
laplace_surv, scatter at1(_Idrug_2=1 _Idrug_3=0) ///
                at2( _Idrug_2=0 _Idrug_3=1) ///
			    at3( _Idrug_2=0 _Idrug_3=0) lp show 
exit
qui xi: laplacereg studytime age i.drug, q(10(10)50)
laplace_surv, scatter at1(age=50 _Idrug_2=1 _Idrug_3=0) ///
               at2(age=50 _Idrug_2=0 _Idrug_3=1) ///
			    at3(age=50 _Idrug_2=0 _Idrug_3=0)
