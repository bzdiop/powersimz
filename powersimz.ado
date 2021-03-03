* Created: 	   01/20/2016 | bzdiop
* Last Modified:   02/23/2016 | bzdiop
* comments added:  01/23/2016 | bzdiop




capture program drop powersimz
program define powersimz, rclass
	syntax varname,  COVariates(varlist) [ITERations(integer 150) seed(integer 1) level(varlist) TREATmentshare(real 0.5) takeup(real 1) alpha(real 0.05) power(real 0.80)] 
	
	* Displaying parameters of simulations
	display "Iterations: `iterations'"
	display "Treated: `treatmentshare'"
	display "Take-up: `takeup' of sample"
	display "Covariates of reg: `covariates'"
	display "Alpha: `alpha'"
	display "Power: `power'"
	
	if "`level'" != ""{
		display "Simulations at the `level' level"
	}
	
	if "`level'" == ""{
		display "Simulations at the individual level"
	}
	
	* Pop-up when there is no temp global created
	if missing("${temp}") {
		window stopbox stop "You must create a global named temp"  "to do that type in the command box:" "global temp [path]"
	}
	 
	preserve 
		
		* Stores the seeds to loop over for each simulation
		local minseed `seed'
		local maxseed = `minseed' + `iterations'
		
		* loops over (maxseed-minseed) times
		forvalues seeds = `minseed'/`maxseed'{
		
		*************************************************************
		*  	RUNS IF RANDOMIZATION IS DONE AT THE `level LEVEL 		*
		*************************************************************
		
			if "`level'" != ""{
			
				* randomizes at the `level' level
				quietly: save "$temp/templevel.dta", replace
					* keeps one observation per level
					keep `level'
					quietly: duplicates drop
						
					* counts the number of observation to then generate the number of observation that should be in T and in C
					qui: count 
					quietly: gen nb_slots    = r(N) * `treatmentshare'
				
					* sets seed for randomization (from minseed to maxseed)
					set seed `seeds'
					
					sort `level'
					
					** Randomization **
					quietly: gen random_draw = runiform()
					sort random_draw
					quietly: egen app_num  = rank(random_draw)
					quietly: gen treatment = (app_num <= nb_slots)
					
					* keeps the relevant variables
					keep `level' treatment
					quietly: save "$temp/randfile.dta", replace
				
				
				use "$temp/templevel.dta", clear
				********************
			
				* merges randomized list of levels to full sample 
				qui: merge m:1 `level' using "$temp/randfile.dta", nogen
				erase "$temp/templevel.dta"
				erase "$temp/randfile.dta"
				********************
				
				* Simulataion of the outcome regression
				qui: reg `varlist' treatment `covariates', cluster (`level')
				
				
				* Stores the outcome of the regression: Betas and SEs
				local beta_level`seeds' = _b[treatment]
				local se_level`seeds'   = _se[treatment]
				
				
				* drops randomization variables
				drop treatment 
			}
			

		*****************************************************************
		*  				RUNS IF FOR THE individual LEVEL 				*
		*****************************************************************
		
			* counts the number of observation to then generate the number of observation that should be in T and in C
			qui: count 
			gen nb_slots    = r(N) * `treatmentshare'
			
			* sets seed for randomization (from minseed to maxseed)
			set seed `seeds'
			
			sort `covariates'
			
			** Randomization (individual level) **
			gen random_draw = runiform()
			sort random_draw
			egen app_num  = rank(random_draw)
			gen treatment = (app_num <= nb_slots)
			********************
			
			* Simulataion of the outcome regression
			qui: reg `varlist' treatment `covariates'
			
			* Stores the outcome of the regression: Betas and SEs
			local beta_ind`seeds' = _b[treatment]
			local se_ind`seeds'   = _se[treatment]
			
			* drops randomization variables
			drop nb_slots random_draw app_num treatment 
		} 

		
		*****************************************************************
		*  	      Creates dataset of estimates and then the MDEs 		*
		*****************************************************************
		
		* From here on, this section creates a datasets of all estimates and uses them to estimate the MDEs
			clear 
			quietly: set obs `maxseed' // creates a dataset large enough to to store all estimates
			foreach i in beta_ind se_ind beta_level se_level {
				quietly: gen `i' = .
			} 



			foreach stat in beta se{ // over beta and se
				forvalues seeds = `minseed'/`maxseed' { // over all iterations
					quietly: replace `stat'_ind  = ``stat'_ind`seeds'' in `seeds' 
					
					if "`level'" != "" { //estimations of MDEs at the level of randomization different from ind.
						quietly: replace `stat'_level  = ``stat'_level`seeds'' in `seeds'
					}
					}
			}
		
		
		

		if "`level'" != ""{ // portion ran if there is level randomization
		*rename (v3 v4) (level_estimates level_standard_errors)

		
		
		
		* Generates Empirical and theoretical MDE based on simulations
		
				
		/* i) Theoretical MDEs: Based on the SE's returned by lm (or glm or ivreg), which are based on modeling assumptions (e.g. normally
              distributed errors in the case of OLS, where the SE is just sqrt(sigma*(X'X)^-1)). */
		
		* levels
		quietly: sum se_level
		local level_theoretical_mde = r(mean) * (invnormal(1 - (.5 * `alpha')) + invnormal(`power'))  / `takeup'
		
		

		/* ii) Empirical MDEs: Permutation based, so there are no modeling assumptions, it just takes the SD of the permutation sample (i.e. the
				sample of beta-hats from the simulations). The empirical version is based on the logic of permutation tests, which shuffle the
				treatment vector in order to simulate the null distribution. */
		
		* levels
		quietly: sum beta_level
		local level_empirical_mde   = r(sd)   * (invnormal(1 - (.5 * `alpha')) + invnormal(`power')) / `takeup'
		
	
		* stores levels
		return local level_theoretical_mde `level_theoretical_mde' 
		return local level_empirical_mde   `level_empirical_mde'
		
		
		}
		
		
		else {
		drop beta_level se_level
		}
		


		/* i) Theoretical MDEs: Based on the SE's returned by lm (or glm or ivreg), which are based on modeling assumptions (e.g. normally
              distributed errors in the case of OLS, where the SE is just sqrt(sigma*(X'X)^-1)). */
			  
		* individuals
		quietly: sum se_ind
		local  indiv_theoretical_mde = r(mean) * (invnormal(1 - (.5 * `alpha')) + invnormal(`power'))  / `takeup'
		

		/* ii) Empirical MDEs: Permutation based, so there are no modeling assumptions, it just takes the SD of the permutation sample (i.e. the
				sample of beta-hats from the simulations). The empirical version is based on the logic of permutation tests, which shuffle the
				treatment vector in order to simulate the null distribution. */
				
		* individuals
		quietly: sum beta_ind
		local indiv_empirical_mde0   =  r(sd)   * (invnormal(1 - (.5 * `alpha')) + invnormal(`power')) 
		local indiv_empirical_mde   = `indiv_empirical_mde0' / `takeup'
		

		
		* stores the outcomes in return list
		return local indiv_theoretical_mde `indiv_theoretical_mde' 
		return local indiv_empirical_mde   `indiv_empirical_mde'

		
		
		quietly: save "$temp/simulations_`varlist'_${date}.dta", replace
	
	restore
	

end	
