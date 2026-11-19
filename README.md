## Age Structured SEIHRD model with vaccination

This is the same SEIHRD model as in my other repository, optimisingNPIs, but translated from Python to R and wrapped up in a package. I hope that this code is (or, will be) far more readable than its Python predecessor, with proper commenting and documentation.

Only the epidemiological model is included, no cost functions. The overall function of this package is to obtain plots (or numerical values) for the different compartments within a population, for given lockdown duration(s), the possibility of a second variant, and school closure. The epidemiological parameters (such as transmission, probability of death, etc.) can all be changed within the 'parameters.json' file. Note that a selection of parameters are age-dependent; here, elderly people (65+) have higher chance of hospitalisation and death, and are also vaccinated first - following what was seen during the SARS-CoV-2 pandemic. The default here is that once 40% of the elderly population are vaccinated, adults (20-64) may be vaccinated, and once 40% of that population are vaccinated, those aged 0-19 can be. 


### Usage
To install this package, simply run:
`devtools::install_github("rebeccamann/seihrdModelAgeStruct")`

To run the model, the following parameters must be defined
```{r}
days            # number of days to run simulation for
eta             # vaccine efficacy
vaccineStart    # day to start vaccinating the elderly
outbreak2Start  # day where a second variant emerges
lockdown1Start  # start day of lockdown 1
lockdown1End    # end day of lockdown 1
lockdown2Start  # start day of lockdown 2
lockdown2End    # end day of lockdown 2
schoolClosure   # TRUE (closed) or FALSE (open) according to whether schools should be closed during the lockdowns
```
Please note, that to set a certain lockdown to 'off', it should be set to `days+1` (NOT 0!). Similarly, if `lockdown1Start=days+1`, then it must also be true that `lockdown2Start=days+1`. 
To obtain results:
```{r}
# make an instance of the simAllAges class:
simulation <- simAllAges$new(eta, vaccineStart, lockdown1Start, lockdown1End, lockdown2Start, lockdown2End, outbreak2Start, days, schoolClosure)

# run the ODE solver:
results <- simulation$ODEsolver()  # this is a list of length 7, each entry containing a the outputs for S,E,I,H,R,D,V respectively.

# plot these results:
plotter(results, colourPal, listOfCompartments)
```
where colourPal is a colour palette of the user's choice, and listOfCompartments is a list of any combination of "S","E","I","H","R","D","V"  according to what the user wants to plot (for example, `list("I","H")` will plot two graphs: for the number of infections, and hospitalisations).
