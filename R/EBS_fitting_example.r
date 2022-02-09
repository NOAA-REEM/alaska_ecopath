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
  # Load a fresh ecosystem (including re-loading csvs)
    source("EBS_scene0_setup.r")
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



########
# Mzero fitting example    
  fit.function <- rsim.fit.run.mzero
  fit.vector   <- scene0$params$MzeroMort[test_sp]

#######
# Predvul fitting example
  fit.function <- rsim.fit.run.predvuls
  fit.vector   <- rep(0, length(test_sp))
  names(fit.vector) <- test_sp

#######
# Preyvul fitting example
  fit.function <- rsim.fit.run.preyvuls
  fit.vector   <- rep(0, length(test_sp))
  names(fit.vector) <- test_sp

# Core fitting routine
  out.vector   <- fit.vector
  fit.initial  <- fit.function(fit.vector, scene=scene0, species=names(fit.vector), verbose=T,
                               run_method='AB', years=fit.years)

  # Brent is for 1-variable problems only
  #fit.result   <- optim(fit.vector, fit.function, method="Brent", lower=-3, upper=3, 
  #                scene=scene0, species=names(fit.vector), run_method='AB', years=fit.years)

  fit.optim    <- optim(fit.vector, fit.function, #lower=-3, upper=3, 
                  scene=scene0, species=names(fit.vector), run_method='AB', years=fit.years)    

  out.vector[names(fit.vector)] <- fit.optim$par
  fit.final   <- fit.function(out.vector, scene=scene0, species=names(out.vector), verbose=T,
                 run_method='AB', years=fit.years)
    
  par(mfrow=c(2,2))
  rsim.plot.biomass(scene0, fit.initial, test_sp[1])
  rsim.plot.biomass(scene0, fit.initial, test_sp[2])
  rsim.plot.biomass(scene0, fit.final, test_sp[1])
  rsim.plot.biomass(scene0, fit.final, test_sp[2])
  
  rsim.fit.table(scene0,fit.initial)
  rsim.fit.table(scene0,fit.final)
  rsim.fit.table(scene0,fit.final) - rsim.fit.table(scene0,fit.initial)


  predvuls <- rep(0,length(test_sp))


  #par(mfrow=c(2,1))
    #rsim.plot.biomass(scene0, fit.initial, test_sp)
    #rsim.plot.biomass(scene0, fit.final, test_sp)


    test <- optimize(rsim.fit.run, c(-1,1),  
                  scene=scene0, species=names(mz_test),run_method='AB', years=fit.years)

    par(mfrow=c(2,1))
    rsim.plot.biomass(scene0, run1, "W.pollock_Adu")
    rsim.plot.biomass(scene0, run1, test_sp)
    rsim.fit.obj(scene0, run1, FALSE)


  scene0$fitting$Biomass$Type[scene0$fitting$Biomass$Group==test_sp] <- "absolute"

  run0 <- rsim.run(scene0, method='AB', years=fit.years)
  rsim.plot.biomass(scene0, run0, test_sp)
  rsim.fit.obj.species(scene0,run0,test_sp)$Biomass


  scene0$fitting$Biomass$wt[scene0$fitting$Biomass$Group==test_sp] <- 0
  scene0$fitting$Biomass["999","wt"] <- 1
  scene0$fitting$Biomass$Type[scene0$fitting$Biomass$Group==test_sp] <- "index"
  run0 <- rsim.run(scene0, method='AB', years=fit.years)
  rsim.plot.biomass(scene0, run0, test_sp)
  rsim.fit.obj.species(scene0,run0,test_sp)$Biomass

  scene0$fitting$Biomass["999","wt"] <- 1
  scene0$fitting$Biomass$Type[scene0$fitting$Biomass$Group==test_sp] <- "index"
  run0 <- rsim.run(scene0, method='AB', years=fit.years)
  rsim.plot.biomass(scene0, run0, test_sp)
  rsim.fit.obj.species(scene0,run0,test_sp)$Biomass


mz <- scene0$params$MzeroMort
rsim.fit.function(mz,scene0,fit.years)

mz[] <- 0
scene1 <- scene0
scene1$params$MzeroMort <- mz
run1 <- rsim.run(scene1, method='AB', years=fit.years)



  qdat <- scene$fitting$Biomass[scene$fitting$Biomass$Group==species,]
  
  #survey_q <- 1
  mn   <- qdat$obs_scaled #/survey_q
  up   <- mn + 1.96*qdat$sd * qdat$survey_q #/survey_q
  dn   <- mn - 1.96*qdat$sd * qdat$survey_q #/survey_q 
  tot  <- 0 #tot <- sum(qdat$fit)
  plot(as.numeric(rownames(run$annual_Biomass)),run$annual_Biomass[,species],type="l",
       ylim=c(0,max(up,run$annual_Biomass[,species])),xlab=tot,ylab="")
  mtext(side=2, line=2.2, paste(species,"biomass"), font=2, cex=1.0)
  points(as.numeric(qdat$Year),mn)





# Temporary = convert grid catch to series catch (prior to read-fitting.catch)
catchdat <- as.matrix(read.csv("data/EBS_ACLIM_72_BIO_catch_2017_NA.csv", header=T, row.names=1, dec="."))
Value <- as.vector(catchdat)
Year  <- as.vector(rownames(catchdat)[row(catchdat)])
Group <- as.vector(colnames(catchdat)[col(catchdat)])
catchframe <- data.frame(Series="EBS_ACLIM_72_BIO_catch_2017", Type="catch", 
              Year, Group, Value, Stdev=0.1*Value, Scale=2.14E-06)
write.csv(catchframe,"data/ebs_catchout_EBS_ACLIM_72_2017_NA.csv",row.names=F)


# Equilibrium by-species F-rate
F_equil <- (rowSums(bal$Landings) + rowSums(bal$Discards))/(bal$Biomass); 
F_equil  <- F_equil[c(all_living,all_detritus)]
# F_equil was zero for KAM and Eelpouts. I'm setting KAM F_equil equal to the 1992
# exploitation rate in the 2018 SA (Bryan et al. 2018). Also, F_equil is 0 for Eelpouts.
# I'm assigning them the F_equil of Motile.epifauna, which is their trophic guild.
F_equil["Kamchatka"] <- 0.02
F_equil["Eelpouts"] <- F_equil["Motile.epifauna"]
B_equil  <- bal$Biomass; #names(B_equil) <- bal$Group; 
B_equil  <- B_equil[c(all_living,all_detritus)] 

capped_sp <- c("Arrowtooth_Adu","Atka.mackerel","FH.sole","Gr.Turbot_Adu","Kamchatka",
               "North.rockfish","Octopi","Other.flatfish","Oth.rockfish","P.cod_Adu",
               "AK.Plaice","POP","W.pollock_Adu","N.Rock.sole","Rougheye.rock","Sablefish",
               "Lg.Sculpins","Sharks","Shortraker.rock","Skates","Squids","YF.sole")

names(capped_sp)<-c("Arrowtooth", "Atka", "Flathead", "Greenland", "Kamchatka",
                    "Northern", "Octopus", "OtherFlat", "OtherRock", "PCod", 
                    "Plaice", "POP", "Pollock", "Rock", "Rougheye", "Sablefish", 
                    "Sculpin", "Shark", "Shortraker", "Skate", "Squid", "Yellowfin")

crab_sp <- c("Bairdi", "King.Crab", "Opilio") 
# Species for which we have no historical catch data.
F_species <- c("Toothed.whales","Baleen.whales","Pinnipeds","Walrus.bd.seals",
               "Ice.seals","Oth.birds","Murres.puffins","Kittiwakes","Auklets",
               "Fulmars","Albatross","Salmon.returning")
# Species for which we have some historical catch data.
F_partial <- c("Sharks","Kamchatka","Other.flatfish","Skates","Eelpouts","Deep.demersals",
               "Shallow.demersals","Lg.Sculpins","Octopi","Mycto.bathy",
               "Capelin","Sandlance","Other.forage","Pandalidae","Ben.zooplankton",
               "Motile.epifauna","Structural.epifauna","Infauna","Jellyfish")

# Year Ranges
#fit.years <- 1971:2016
fit.years91 <- 1991:2017
fore_years <- 2018:2039
all_years91  <- c(fit.years91,fore_years)

# Baseline scenario SETUP - Starting with 1991
basescene91 <- rsim.scenario(bal, unbal, years = all_years91) # Ecosim params

scene0 <- basescene91
scene1 <- read.fitting.biomass(scene0, "data/ebs_fitout.csv")

run1 <- rsim.run(scene1, method='AB', years=all_years91)
rsim.plot.biomass(scene1, run1, "W.pollock_Adu")
  
fit.species <- unique(scene1$fitting$Biomass$Group)
rsim.plot.biomass(scene1, run1, "Pandalidae")


# -----***** UPDATE *****-----#
# Because we are starting in 1991 in this script there is no need to tweak starting
# biomasses or mortality rates.

# Set up hindcast fishing, turn off fishing effort
scene.full91 <- basescene91

# Zero our effort, freeze Discards.offal
scene.full91 <- adjust.fishing(scene.full91, "ForcedEffort", all_gears, all_years91, value=0.0)

scene.full91$forcing$ForcedBio[,"Discards.offal"] <- B_equil["Discards.offal"]

# Read fisheries Catch file, convert units

catchdat91 <- catchdat[21:47,]
catchdat91 <- catchdat91/495218
# assumes first year in input file is 1970, copy hindcast catch to CATCH
catchdat91 <- as.matrix(catchdat91[1:(length(fit.years91)),])
scene.full91$fishing$ForcedCatch[1:(length(fit.years91)),] <- catchdat91 
# Add historical F for species that we do not have historical catch for
for (i in 1:(length(F_species))) {
  scene.full91$fishing$ForcedFRate[1:(length(fit.years91)),F_species[i]] <- F_equil[F_species[i]]
}
#' There are several (20) species that we only have a limited amount of historical
#' catch data for. The other years of the hindcast have zero catch. However, there
#' was catch. To remedy this I will attempt to apply the equilibrium F (F_equil)
#' for each species to those years for which we have no catch. Otherwise there is
#' no catch for those species in those years and that affects biomass dynamics.
for (i in 1:(length(F_partial))) {
  #print(i)
  scene.full91$fishing$ForcedFRate[1:27,F_partial[i]] <- 
    ifelse(scene.full91$fishing$ForcedCatch[1:27,F_partial[i]]<=0.0,
      F_equil[F_partial[i]],0.0)
}

# Run hindcast years
run.hind91 <- rsim.run(scene.full91,method='AB',years=fit.years91)
rsim.plot(run.hind91)


# Catch data

dat <- read.csv("EBS_ACLIM_72_BIO_catch_2017_condensed.csv",row.names="Year")