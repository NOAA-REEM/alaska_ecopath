library(Rpath)

# Model file location
  eco_folder <- "https://raw.githubusercontent.com/kaydin/alaska_ecopath/main/EBS_ACLIM_72/"

# File names
  Ebase <- paste(eco_folder, "EBS_ACLIM_72_base.csv", sep='')  
  Ediet <- paste(eco_folder, "EBS_ACLIM_72_diet.csv", sep='')   
  Eped  <- paste(eco_folder, "EBS_ACLIM_72_pedigree.csv", sep='')    
  Estz  <- paste(eco_folder, "EBS_ACLIM_72_stanzas.csv", sep='')   
  Estg  <- paste(eco_folder, "EBS_ACLIM_72_stanza_groups.csv",sep='')  

  biomass.datafile  <- paste(eco_folder, "EBS_ACLIM_72_biomass_fitdata.csv",sep='')
  catch.datafile    <- paste(eco_folder, "EBS_ACLIM_72_catch_fitdata.csv",sep='')

  fit.years <- 1991:2017

# Setup Base Ecopath and Base Rsim scenario
  unbal <- rpath.stanzas(read.rpath.params(Ebase, Ediet, Eped, Estg, Estz)) 
  bal   <- rpath(unbal)
  basescene91 <- rsim.scenario(bal, unbal, years = fit.years) # Ecosim params
  scene0 <- basescene91

# Read in fitting data
# Biomass data (e.g. surveys)
  scene0 <- read.fitting.biomass(scene0, biomass.datafile)

# Read time series of catch data and re-arrange catch forcing
  scene0 <- read.fitting.catch(scene0, catch.datafile) 
  # Apply the fit catch as a forcing catch
  scene0 <- fitcatch.to.forcecatch(scene0)
  # Turn off fishing effort, freeze discards/offal using forced biomass
  scene0 <- adjust.fishing(scene0, "ForcedEffort", rpath.gears(bal), fit.years, value=0.0)
  scene0$forcing$ForcedBio[,"Discards.offal"] <- bal$Biomass["Discards.offal"]
 
  # For species without catch, reapply Ecopath F (originally through gears) to ForcedFRate
  F_equil <- (rowSums(bal$Landings) + rowSums(bal$Discards))/(bal$Biomass) 
  Equil_species <- c("Transient.killer.whales", "Toothed.whales", "Gray.whales", 
              "Baleen.whales", "Pinnipeds", "Walrus.bd.seals", "N.fur.seal_Juv", 
              "Ice.seals", "Oth.birds", "Murres.puffins", "Kittiwakes", "Auklets", 
              "Fulmars", "Albatross", "W.pollock_Juv", "P.cod_Juv", "Arrowtooth_Juv", 
              "Gr.Turbot_Juv", "P.halibut_Juv", "Pel.zooplankton", "Euphausiids", 
              "Copepods", "Pelagic.microbes", "Benthic.microbes", "Lg.phytoplankton", 
              "Sm.phytoplankton", "Discards.offal", "Pelagic.Detritus", "Benthic.Detritus")
  for (sp in Equil_species){
     scene0 <- adjust.fishing(scene0, 'ForcedFRate', sp, fit.years, value=F_equil[sp])
  }


# Run model
  run0 <- rsim.run(scene0, method='AB', years=fit.years)
  par(mfrow=c(1,2))
  rsim.plot.biomass(scene0, run0, "W.pollock_Adu")
  rsim.plot.catch(scene0, run0, "W.pollock_Adu")

# Some Diagnostics 
  rsim.fit.table(scene0,run0)
  rsim.fit.obj.species(scene0,run0,"W.pollock_Adu")

# Testing switching between absolute and index
  test_sp <- "Pandalidae" 
  par(mfrow=c(1,2))
  scene0$fitting$Biomass$Type[scene0$fitting$Biomass$Group==test_sp] <- "absolute"
  run0 <- rsim.run(scene0, method='AB', years=fit.years)
  rsim.plot.biomass(scene0, run0, test_sp)

  #rsim.fit.obj.species(scene0,run0,test_sp)$Biomass
  scene0$fitting$Biomass$Type[scene0$fitting$Biomass$Group==test_sp] <- "index"
  run0 <- rsim.run(scene0, method='AB', years=fit.years)
  rsim.plot.biomass(scene0, run0, test_sp)
  #rsim.fit.obj.species(scene0,run0,test_sp)$Biomass

#Functon Update
  # Species to test 
    test_sp <- c("W.pollock_Adu","Pandalidae")
    data_type <- "absolute"  #"index"
  # Set data weightings for all data input low (zeros not allowed)
    scene0$fitting$Biomass$wt[] <- 1e-36
    scene0$fitting$Catch$wt[]   <- 1e-36
  # Set data type for test species
    scene0$fitting$Biomass$Type[scene0$fitting$Biomass$Group %in% test_sp] <- data_type
  # Set data weighting for one species to fit to 1
    scene0$fitting$Biomass$wt[scene0$fitting$Biomass$Group %in% test_sp]   <- 1

# Different fitting examples
  # mzero only
  fit_values   <- rep(0,length(test_sp))
  fit_species  <- test_sp
  fit_vartype  <- rep("mzero",length(test_sp))

  # predvul only
  fit_values   <- rep(0,length(test_sp))
  fit_species  <- test_sp
  fit_vartype  <- rep("predvul",length(test_sp))

  # preyvul only
  fit_values   <- rep(0,length(test_sp))
  fit_species  <- test_sp
  fit_vartype  <- rep("preyvul",length(test_sp))

  # all combined
  fit_values   <- c(rep(0,length(test_sp)),rep(0,length(test_sp)),rep(0,length(test_sp))) 
  fit_species  <- c(test_sp,test_sp,test_sp)
  fit_vartype  <- c(rep("mzero",length(test_sp)),
                    rep("predvul",length(test_sp)),
                    rep("preyvul",length(test_sp)))

# Core fitting procedure here
data.frame(fit_vartype,fit_species,fit_values)

par(mfrow=c(2,2))

fit.initial  <- rsim.fit.run(fit_values, fit_species, fit_vartype, scene0, verbose=T,
                             run_method='AB', years=fit.years)
rsim.plot.biomass(scene0, fit.initial, test_sp[1])
rsim.plot.biomass(scene0, fit.initial, test_sp[2])

# Run optimization
fit.optim    <- optim(fit_values, rsim.fit.run, #lower=-3, upper=3, 
                species=fit_species, vartype=fit_vartype, scene=scene0,   
                run_method='AB', years=fit.years) 
out_values <- fit.optim$par

data.frame(fit_vartype,fit_species,fit_values,out_values)
fit.final  <- rsim.fit.run(out_values, fit_species, fit_vartype, scene0, verbose=T,
                           run_method='AB', years=fit.years) 
rsim.plot.biomass(scene0, fit.final, test_sp[1])
rsim.plot.biomass(scene0, fit.final, test_sp[2])


scene1 <- rsim.fit.update(out_values, fit_species, fit_vartype, scene0)
  run1 <- rsim.run(scene1, method='AB', years=fit.years)
rsim.plot.biomass(scene1, run1, test_sp[1])
rsim.plot.biomass(scene1, run1, test_sp[2])

# Comparing scenario vuls between scene0 and scene1
get.rsim.predprey(scene0,"Pandalidae","all")
get.rsim.predprey(scene0,"all","Pandalidae")

get.rsim.predprey(scene1,"Pandalidae","all")
get.rsim.predprey(scene1,"all","Pandalidae")
# > data.frame(fit_vartype,fit_species,fit_values,out_values)
#   fit_vartype   fit_species fit_values  out_values
# 1       mzero W.pollock_Adu          0  0.07841237
# 2       mzero    Pandalidae          0  0.32816386
# 3     predvul W.pollock_Adu          0 -3.08705448
# 4     predvul    Pandalidae          0 -0.66494780
# 5     preyvul W.pollock_Adu          0 -0.78028658
# 6     preyvul    Pandalidae          0  0.80811375


# Combining the objective function with sense

sense.scene <- scene1
# From ecosense vignette
# Note that fit_years is less than burn_years here!  How do we do sense
# with data (e.g. catch data) when less than 50 years of data?
sense.scene$params$BURN_YEARS <- 50
NUM_RUNS <- 1000 # how many ecosystem parameter sets to generate
parlist<-as.list(rep(NA,NUM_RUNS)) # create lists to store the generated parameters
kept<-rep(NA,NUM_RUNS) # object to keep track of kept systems
set.seed(666) # Optional, set seed so output can be replicated
for (i in 1:NUM_RUNS){
  EBSsense <- sense.scene # scenario object
  # INSERT SENSE ROUTINE BELOW
  parlist[[i]]<- sense.scene$params 		# Base ecosim params
  parlist[[i]]<- rsim.sense(sense.scene,unbal,Vvary = c(-4.5,4.5), Dvary = c(0,0))	# Replace the base params with Ecosense params
  EBSsense$start_state$Biomass <- parlist[[i]]$B_BaseRef # Apply the Ecosense starting biomass
  parlist[[i]]$BURN_YEARS <- 50			# Set Burn Years to 50
  EBSsense$params <- parlist[[i]] # replace base params with the Ecosense generated params
  EBStest <- rsim.run(EBSsense, method="AB", years=fit.years) # Run rsim with the generated system
  failList <- which(is.na(EBStest$end_state$Biomass))
  {if (length(failList)>0)
  {cat(i,": fail in year ",EBStest$crash_year,": ",EBStest$params$spname[failList],rsim.fit.obj(EBSsense,EBStest,F),"\n"); kept[i]<-F; flush.console()}
    else 
    {cat(i,": success!",rsim.fit.obj(EBSsense,EBStest,F),"\n"); kept[i]<-T;  flush.console()}} # output for the console
  parlist[[i]]$BURN_YEARS <- 1
}

KEPT <- which(kept==T) # the number associated with the kept system
nkept <- length(KEPT); nkept # how many were kept
1-(nkept/NUM_RUNS) # rejection rate
ecos <- as.list(rep(NA,length(KEPT))) # lists for simulated ecosystems
k <- 0  # counter for simulated ecosystems
for (i in KEPT) {
  EBSsense <- scene1 # set up the scenario object
  EBSsense$start_state$Biomass <- parlist[[i]]$B_BaseRef # set the starting Biomass to the generated values

  EBSsense$params <- parlist[[i]] # set the params in the scenario object equal to the generated params.
  EBSsense$BURN_YEARS <- -1 # no burn-in period
  k <- k + 1 # set the number for the simulated ecosystem
  ecos[[k]] <- rsim.run(EBSsense,method='AB', years=fit.years) # run rsim.run on the generated system
  cat("Ecosystem no.", k, "out of", nkept, rsim.fit.obj(EBSsense,ecos[[k]],F),"\n"); flush.console() # progress output to console
}

par(mfrow=c(3,3))
rsim.plot.biomass(scene1, run1, "W.pollock_Adu")
for(i in 1:8){
  EBSsense <- scene1 # set up the scenario object
  EBSsense$start_state$Biomass <- parlist[[i]]$B_BaseRef # set the starting Biomass to the generated values
  EBSsense$params <- parlist[[i]] # set the params in the scenario object equal to the generated params.
  EBSsense$BURN_YEARS <- -1 # no burn-in period
  k <- k + 1 # set the number for the simulated ecosystem
  eco.kept <- rsim.run(EBSsense,method='AB', years=fit.years) # run rsim.run on the generated system
  rsim.plot.biomass(EBSsense, eco.kept, "W.pollock_Adu")
}

