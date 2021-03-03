{smcl}
{* *! version 1.0  02Feb2017}{...}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col:{manlink R powersimz} {hline 2}}Simulation-based minimum detectable effects (MDE) {p_end}
{p2colreset}{...}

{p 4 0 2}
Author: {opt bintazahra.diop@economics.ox.ac.uk}, last update: May, 2019
{p_end}

{p 4 0 2}
{opt Required}: have a $temp global that links to a folder path
{p_end}

{marker syntax}{...}
{title:Syntax}

{p 8 18 2}
{cmdab:powersimz} [{varlist}]
 {cmd:,} 
{cmdab:cov:ariates(}{it:varlist}{cmd:)}
[{it:additional options}
{cmdab:iter:ations(}{it:integer}{cmd:)}
{cmdab:seed(}{it:integer}{cmd:)}
{cmdab:level(}{it:variable}{cmd:)}
{cmdab:treat:mentshare(}{it:integer}{cmd:)}
{cmdab:takeup(}{it:real}{cmd:)}
{cmdab:alpha(}{it:real}{cmd:)}
{cmdab:power(}{it:power}{cmd:)}
]


{marker description}{...}
{title:Description}

{pstd}
{opt powersimz} calculates and stores theoretical and empirical simulation-based Minimum Detectable Effects (MDEs). 

{p 4 6 2}
i) Theoretical MDEs: Based on the SE's returned by the {opt regress} command, which are based on modeling assumptions (e.g. normally distributed errors in the case of OLS, where the SE is just sqrt(sigma*(X'X)^-1)). In this case the MDE is essentially: (tα/2 + t1−κ) x Mean(SEs)
{p_end}

{p 4 6 2}
ii) Empirical MDEs: Permutation based, so there are no modeling assumptions, it just takes the SD of the permutation sample (i.e. the sample of beta-hats from the simulations). The empirical version is based on the logic of permutation tests, which shuffle the treatment vector in order to simulate the null distribution. In this case the MDE is essentially: (t(α/2) + t(1−κ)) x SD(β-hat)
{p_end}

{p 4 0 2}
Note that MDEs are defined by the following: (t(α/2) + t(1−κ))σ(β-hat)
{p_end} 
{p 4 4 2}
The main difference between the two is how we define σ(β-hat), in the “theoretical” version I am assuming that it is the mean of the standard errors of the β-hats from all the simulated randomizations (this makes the underlying assumptions that errors are normally distributed), and in the “empirical” version I are assuming that it is the standard deviation of all the β-hats from the same simulations (this makes no modeling assumptions).
{p_end}




Below are the steps in which these MDEs and computed. 



{p 4 6 2}
{it:if} {opt level} is specified the regression of {varlist} on covariates is ran {opt iterations} times. The steps are the following: 
{p_end}

{p 4 6 2}
{it:   - step1} Randomize at the {opt level} level
{p_end}

{p 4 6 2}
{it:   - step2} Assign each individual to their respective randomization {opt level} 
{p_end}

{p 4 6 2}
{it:   - step3} Runs the regression at the individual level clustering SEs at the {opt level} level
{p_end}

{p 4 6 2}
{it:if} {opt level} is not specified the regression of {varlist} on covariates is ran {opt iterations} times. The steps are the following: 
{p_end}

{p 4 6 2}
{it:   - step1} Randomize all individuals into treatment or control (The number of individuals treated depends on the specified share of treated individuals)
{p_end}

{p 4 6 2}
{it:   - step2} Run regression at the individual level (no clustering).
{p_end}

{p 4 6 2}
If no {opt alpha} is specified, it is assumed at 5%, if no {opt power} is specified it is assumed at 80%, if no {opt takeup} is specified it is assumed at 100%. If no {opt level} is specified, it is assumed to be individual level randomization. 
{p_end}

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt covariates} The list of variables that will be used in the simulated regressions (in addition to the treatment variable). 

{phang}
{opt seed}/{opt iterations}, each regression is ran on an artificial randomization. {cmdab:powersimz} randomizes treatment {opt iterations} number of times. Each iteration is done setting seeds between the specified {opt seed} and {opt seed} + {opt iterations}.  

{dlgtab:Optional}

{phang}
{opt level}({it:varname}) specifies a greater level of randomization than the most disaggregated level in the dataset.  When level is specified, the randomization is done at the {it:varname} level 

{phang}
{opt takeup}({it:numlist})  allows to include sub-perfect take-up. {opt takeup} has to be between 0 and 1. It basically decides the beta by the specified {opt takeup} -- in a 2sls fashion.

{phang}
{opt alpha}({it:numlist})  significance level; default is {opt alpha(0.05)}

{phang}
{opt power}({it:numlist})  power; default is {opt power(0.8)}

{phang}
{opt iterations}({it:numlist})  number of simulations; default is {opt iterations(150)}

{phang}
{opt treatmentshare}({it:numlist})  share of the sample that will be randomized into treatment; default is {opt treatment(0.5)}

{phang}
{opt seed}({it:numlist}) seed for the randomization of the first simulation; default is {opt seed(1)}


{marker examples}{...}
{title:Examples}

{cmd:. sysuse} auto.dta, clear 

{cmd:. powersimz} price, cov(mpg headroom weight) seed(2019) iterations(100)
{cmd:. return list}

{cmd:. powersimz} price, cov(mpg headroom weight) level(trunk) seed(2019) iterations(100)
{cmd:. return list}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:powersimz} stores the following in {cmd:r()}:

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:r(indiv_theoretical_mde)}} stores the SE of the Betas from all regressions -- with randomizations the most disaggregated level in the dataset {p_end}
{synopt:{cmd:r(indiv_empirical_mde)}} stores the mean of the SE from all regressions -- with randomizations the most disaggregated level in the dataset {p_end}
{synopt:{cmd:r(level_theoretical_mde)}} stores the SE of the Betas from all regressions -- with randomizations the level specified in {opt level()} {p_end}
{synopt:{cmd:r(level_empirical_mde)}} stores the mean of the SE from all regressions -- with randomizations the level specified in {opt level()} {p_end}
{p2colreset}{...}
