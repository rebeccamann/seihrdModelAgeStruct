

ODEsClass <- setRefClass("ODEsClass", 
                         fields=list(params="ANY", initialVals="ANY"),
                         methods=list(
                           initialize = function(ageGrpParams, eta, lockdown1Start,
                                                 lockdown1End, lockdown2Start,
                                                 lockdown2End, outbreak2Start, days,
                                                 schoolClosure){
                             .self$params <- ageGrpParams
                             .self$params$eta1 <- eta
                             .self$params$eta2 <- 0.8*eta
                             .self$params$schoolClosure <- schoolClosure
                             .self$params$lockdown1Start <- lockdown1Start
                             .self$params$lockdown1End <- lockdown1End
                             .self$params$lockdown2Start <- lockdown2Start
                             .self$params$lockdown2End <- lockdown2End
                             .self$params$outbreak2Start <- outbreak2Start
                             .self$params$days <- days
                             
                             .self$initialVals <- unlist(ageGrpParams$compartments,
                                                         use.names = FALSE)
                           },
                           betaFunc = function(time, var, ageGrp){
                             beta <- .self$params[[var]][ageGrp]
                             if((.self$params$schoolClosure == FALSE) && 
                                (ageGrp ==1)){
                               newBeta <- beta
                             } else if((time > .self$params$lockdown1Start) && 
                                       (time <= .self$params$lockdown1End)){
                               newBeta <- (beta / .self$params$compliance) * (1 + (.self$params$compliance - 1) * exp(-(time - .self$params$lockdown1Start)))
                             } else if((time > .self$params$lockdown2Start) &&
                                       (time <= .self$params$lockdown2End)){
                               newBeta <- (beta / .self$params$compliance) * (1 + (.self$params$compliance - 1) * exp(-(time - .self$params$lockdown2Start)))
                             } else if((time > .self$params$lockdown1End) && 
                                       (time <= .self$params$lockdown2Start)){
                               newBeta <- beta + (beta * (1 - .self$params$compliance) / .self$params$compliance) * exp(-.self$params$easing * (time - .self$params$lockdown1End))
                             } else if(time > .self$params$lockdown2End){
                               newBeta <- beta + (beta * (1 - .self$params$compliance) / .self$params$compliance) * exp(-.self$params$easing * (time - .self$params$lockdown2End))
                             } else{
                               newBeta <- beta
                             }
                             return(newBeta)
                           },
                           defineVars = function(IVS, FOI1, FOI2, state){
                             
                             S_U <- c(.self$params$omega1*IVS[31] - FOI1*IVS[1]/.self$params$N,
                                      .self$params$omega1*IVS[31] - FOI1*IVS[1]/.self$params$N - .self$params$sigma*IVS[1],
                                      .self$params$omega1*IVS[31] + .self$params$omega2*IVS[33] - (FOI1+FOI2)*IVS[1]/.self$params$N,
                                      .self$params$omega1*IVS[31] + .self$params$omega2*IVS[33] - (FOI1+FOI2)*IVS[1]/.self$params$N -
                                        .self$params$sigma*IVS[1])[state]
                             
                             S1_V <- c(0,
                                       .self$params$sigma*(1-.self$params$eta1)*IVS[1] + .self$params$omega1*IVS[32] - 
                                         FOI1*IVS[2] / .self$params$N,
                                       0,
                                       .self$params$sigma*(1-.self$params$eta1)*IVS[1] + .self$params$omega1*IVS[32] +
                                         .self$params$omega2*IVS[34] - (FOI1+FOI2)*IVS[2]/.self$params$N)[state]
                             
                             S1_VI <- c(0,
                                        .self$params$sigma*.self$params$eta1*IVS[1],
                                        0,
                                        .self$params$sigma*.self$params$eta1*IVS[1] - IVS[3])[state]
                             
                             S2_V <- c(0,
                                       0,
                                       0,
                                       (1-.self$params$eta2)*IVS[3] + .self$params$omega2*IVS[35] - FOI2*IVS[4]/.self$params$N)[state]
                             
                             S2_VI <- c(0,
                                        0,
                                        0,
                                        .self$params$eta2*IVS[3])[state]
                             
                             E1_U <- c(FOI1*IVS[1]/.self$params$N - .self$params$epsilon1*IVS[6],
                                       FOI1*IVS[1]/.self$params$N - .self$params$epsilon1*IVS[6],
                                       FOI1*IVS[1]/.self$params$N - .self$params$epsilon1*IVS[6],
                                       FOI1*IVS[1]/.self$params$N - .self$params$epsilon1*IVS[6])[state]
                             
                             E1_V <- c(0,
                                       FOI1*IVS[2]/.self$params$N -.self$params$epsilon1*IVS[7],
                                       0,
                                       FOI1*IVS[2]/.self$params$N -.self$params$epsilon1*IVS[7])[state]
                             
                             E2_U <- c(0,
                                       0,
                                       FOI2*IVS[1]/.self$params$N - .self$params$epsilon2*IVS[8],
                                       FOI2*IVS[1]/.self$params$N - .self$params$epsilon2*IVS[8])[state]
                             
                             E2_V <- c(0,
                                       0,
                                       0,
                                       FOI2*IVS[2]/.self$params$N -.self$params$epsilon2*IVS[9])[state]
                             
                             E2_VI <- c(0,
                                        0,
                                        0,
                                        FOI2*IVS[4]/.self$params$N -.self$params$epsilon2*IVS[10])[state]
                             
                             I1_SU <- c(.self$params$epsilon1*.self$params$d1*IVS[6] - .self$params$gamma1prime*IVS[11],
                                        .self$params$epsilon1*.self$params$d1*IVS[6] - .self$params$gamma1prime*IVS[11],
                                        .self$params$epsilon1*.self$params$d1*IVS[6] - .self$params$gamma1prime*IVS[11],
                                        .self$params$epsilon1*.self$params$d1*IVS[6] - .self$params$gamma1prime*IVS[11])[state]
                             
                             I1_AU <- c(.self$params$epsilon1*(1-.self$params$d1)*IVS[6] - .self$params$gamma1*IVS[12],
                                        .self$params$epsilon1*(1-.self$params$d1)*IVS[6] - .self$params$gamma1*IVS[12],
                                        .self$params$epsilon1*(1-.self$params$d1)*IVS[6] - .self$params$gamma1*IVS[12],
                                        .self$params$epsilon1*(1-.self$params$d1)*IVS[6] - .self$params$gamma1*IVS[12])[state]
                             
                             I1_RU <- c(.self$params$gamma1prime*(1 - .self$params$h1)*IVS[11] - .self$params$delta1prime*IVS[13],
                                        .self$params$gamma1prime*(1 - .self$params$h1)*IVS[11] - .self$params$delta1prime*IVS[13],
                                        .self$params$gamma1prime*(1 - .self$params$h1)*IVS[11] - .self$params$delta1prime*IVS[13],
                                        .self$params$gamma1prime*(1 - .self$params$h1)*IVS[11] - .self$params$delta1prime*IVS[13])[state]
                             
                             I1_SV <- c(0,
                                        .self$params$epsilon1*.self$params$d1*IVS[7] - .self$params$gamma1prime*IVS[14],
                                        0,
                                        .self$params$epsilon1*.self$params$d1*IVS[7] - .self$params$gamma1prime*IVS[14])[state]
                             
                             I1_AV <- c(0,
                                        .self$params$epsilon1*(1-.self$params$d1)*IVS[7] - .self$params$gamma1*IVS[15],
                                        0,
                                        .self$params$epsilon1*(1-.self$params$d1)*IVS[7] - .self$params$gamma1*IVS[15])[state]
                             
                             I1_RV <- c(0,
                                        .self$params$gamma1prime*(1 - .self$params$h1)*IVS[14] - .self$params$delta1prime*IVS[16],
                                        0,
                                        .self$params$gamma1prime*(1 - .self$params$h1)*IVS[14] - .self$params$delta1prime*IVS[16])[state]
                             
                             I2_SU <- c(0,
                                        0,
                                        .self$params$epsilon2*.self$params$d2*IVS[8] - .self$params$gamma2prime*IVS[17],
                                        .self$params$epsilon2*.self$params$d2*IVS[8] - .self$params$gamma2prime*IVS[17])[state]
                             
                             I2_AU <- c(0,
                                        0,
                                        .self$params$epsilon2*(1-.self$params$d2)*IVS[8] - .self$params$gamma2*IVS[18],
                                        .self$params$epsilon2*(1-.self$params$d2)*IVS[8] - .self$params$gamma2*IVS[18])[state]
                             
                             I2_RU <- c(0,
                                        0,
                                        .self$params$gamma2prime*(1 - .self$params$h2)*IVS[17] - .self$params$delta2prime*IVS[19],
                                        .self$params$gamma2prime*(1 - .self$params$h2)*IVS[17] - .self$params$delta2prime*IVS[19])[state]
                             
                             I2_SV <- c(0,
                                        0,
                                        0,
                                        .self$params$epsilon2*.self$params$d2*IVS[9] - .self$params$gamma2prime*IVS[20])[state]
                             
                             I2_AV <- c(0,
                                        0,
                                        0,
                                        .self$params$epsilon2*(1-.self$params$d2)*IVS[9] - .self$params$gamma2*IVS[21])[state]
                             
                             I2_RV <- c(0,
                                        0,
                                        0,
                                        .self$params$gamma2prime*(1 - .self$params$h2)*IVS[20] - .self$params$delta2prime*IVS[22])[state]
                             
                             I2_SVI <- c(0,
                                         0,
                                         0,
                                         .self$params$epsilon2*.self$params$d2*IVS[10] - .self$params$gamma2prime*IVS[23])[state]
                             
                             I2_AVI <- c(0,
                                         0,
                                         0,
                                         .self$params$epsilon2*(1-.self$params$d2)*IVS[10] - .self$params$gamma2*IVS[24])[state]
                             
                             I2_RVI <- c(0,
                                         0,
                                         0,
                                         .self$params$gamma2prime*(1 - .self$params$h2)*IVS[23] - .self$params$delta2prime*IVS[25])[state]
                             
                             H1_U <- c(.self$params$gamma1prime*.self$params$h1*IVS[11] - .self$params$delta1*IVS[26],
                                       .self$params$gamma1prime*.self$params$h1*IVS[11] - .self$params$delta1*IVS[26],
                                       .self$params$gamma1prime*.self$params$h1*IVS[11] - .self$params$delta1*IVS[26],
                                       .self$params$gamma1prime*.self$params$h1*IVS[11] - .self$params$delta1*IVS[26])[state]
                             
                             H1_V <- c(0,
                                       .self$params$gamma1prime*.self$params$h1*IVS[14] - .self$params$delta1*IVS[27],
                                       0,
                                       .self$params$gamma1prime*.self$params$h1*IVS[14] - .self$params$delta1*IVS[27])[state]
                             
                             H2_U <- c(0,
                                       0,
                                       .self$params$gamma2prime*.self$params$h2*IVS[17] - .self$params$delta2*IVS[28],
                                       .self$params$gamma2prime*.self$params$h2*IVS[17] - .self$params$delta2*IVS[28])[state]
                             
                             H2_V <- c(0,
                                       0,
                                       0,
                                       .self$params$gamma2prime*.self$params$h2*IVS[20] - .self$params$delta2*IVS[29])[state]
                             
                             H2_VI <- c(0,
                                        0,
                                        0,
                                        .self$params$gamma2prime*.self$params$h2*IVS[23] - .self$params$delta2*IVS[30])[state]
                             
                             R1_U <- c(.self$params$delta1*(1-.self$params$k1)*IVS[26] + .self$params$delta1prime*IVS[13] + .self$params$gamma1*IVS[12] - .self$params$omega1*IVS[31],
                                       .self$params$delta1*(1-.self$params$k1)*IVS[26] + .self$params$delta1prime*IVS[13] + .self$params$gamma1*IVS[12] - .self$params$omega1*IVS[31],
                                       .self$params$delta1*(1-.self$params$k1)*IVS[26] + .self$params$delta1prime*IVS[13] + .self$params$gamma1*IVS[12] - .self$params$omega1*IVS[31],
                                       .self$params$delta1*(1-.self$params$k1)*IVS[26] + .self$params$delta1prime*IVS[13] + .self$params$gamma1*IVS[12] - .self$params$omega1*IVS[31])[state]
                             
                             R1_V <- c(0,
                                       .self$params$delta1*(1-.self$params$k1)*IVS[27] + .self$params$delta1prime*IVS[16] + .self$params$gamma1*IVS[15] - .self$params$omega1*IVS[32],
                                       0,
                                       .self$params$delta1*(1-.self$params$k1)*IVS[27] + .self$params$delta1prime*IVS[16] + .self$params$gamma1*IVS[15] - .self$params$omega1*IVS[32])[state]
                             
                             R2_U <- c(0,
                                       0,
                                       .self$params$delta2*(1-.self$params$k2)*IVS[28] + .self$params$delta2prime*IVS[19] + .self$params$gamma2*IVS[18] - .self$params$omega2*IVS[33],
                                       .self$params$delta2*(1-.self$params$k2)*IVS[28] + .self$params$delta2prime*IVS[19] + .self$params$gamma2*IVS[18] - .self$params$omega2*IVS[33])[state]
                             
                             R2_V <- c(0,
                                       0,
                                       0,
                                       .self$params$delta2*(1-.self$params$k2)*IVS[29] + .self$params$delta2prime*IVS[22] + .self$params$gamma2*IVS[21] - .self$params$omega2*IVS[34])[state]
                             
                             
                             R2_VI <- c(0,
                                        0,
                                        0,
                                        .self$params$delta2*(1-.self$params$k2)*IVS[30] + .self$params$delta2prime*IVS[25] + .self$params$gamma2*IVS[24] - .self$params$omega2*IVS[35])[state]
                             
                             D1_U <- c(.self$params$delta1*.self$params$k1*IVS[26],
                                       .self$params$delta1*.self$params$k1*IVS[26],
                                       .self$params$delta1*.self$params$k1*IVS[26],
                                       .self$params$delta1*.self$params$k1*IVS[26])[state]
                             
                             D1_V <- c(0,
                                       .self$params$delta1*.self$params$k1*IVS[27],
                                       0,
                                       .self$params$delta1*.self$params$k1*IVS[27])[state]
                             
                             D2_U <- c(0,
                                       0,
                                       .self$params$delta2*.self$params$k2*IVS[28],
                                       .self$params$delta2*.self$params$k2*IVS[28])[state]
                             
                             D2_V <- c(0,
                                       0,
                                       0,
                                       .self$params$delta2*.self$params$k2*IVS[29])[state]
                             
                             D2_VI <- c(0,
                                        0,
                                        0,
                                        .self$params$delta2*.self$params$k2*IVS[30])[state]
                             
                             Y = c(S_U, S1_V, S1_VI, S2_V, S2_VI, E1_U, E1_V, E2_U, E2_V, E2_VI,I1_SU, I1_AU, I1_RU, I1_SV, I1_AV,
                                   I1_RV,I2_SU, I2_AU, I2_RU, I2_SV, I2_AV, I2_RV, I2_SVI, I2_AVI, I2_RVI,H1_U, H1_V, H2_U, H2_V,
                                   H2_VI,R1_U, R1_V, R2_U, R2_V, R2_VI,D1_U, D1_V, D2_U, D2_V, D2_VI)
                             return(Y)
                           }
                         )
)

simAllAges <- setRefClass("simAllAges",
                          fields = list(params = "ANY", youngSim = "ANY", adultSim = "ANY", 
                                        elderlySim = "ANY", outputs = "ANY", events = "ANY", ageGroups = "ANY",
                                        days = "ANY", outbreak2Start = "ANY", vaccineStart = "ANY"),
                          methods = list(
                            initialize = function(eta, vaccineStart, lockdown1Start,
                                                  lockdown1End, lockdown2Start, lockdown2End,
                                                  outbreak2Start, days, schoolClosure){
                              allAgeGrpParams <- .loadParameters()
                              youngAgeGrp <- allAgeGrpParams[[1]]
                              adultAgeGrp <- allAgeGrpParams[[2]]
                              elderlyAgeGrp <- allAgeGrpParams[[3]]
                              
                              .self$youngSim <- ODEsClass$new(youngAgeGrp, eta, lockdown1Start,
                                                              lockdown1End, lockdown2Start, lockdown2End,
                                                              outbreak2Start, days, schoolClosure)
                              .self$adultSim <- ODEsClass$new(adultAgeGrp, eta, lockdown1Start,
                                                              lockdown1End, lockdown2Start, lockdown2End,
                                                              outbreak2Start, days, schoolClosure = TRUE)
                              .self$elderlySim <- ODEsClass$new(elderlyAgeGrp, eta, lockdown1Start,
                                                                lockdown1End, lockdown2Start, lockdown2End,
                                                                outbreak2Start, days, schoolClosure = TRUE)
                              .self$outputs <- list()
                              
                              .self$events <- list('elderlyVaccStart' = NULL, 'adultVaccStart' = NULL,
                                                   'youngVaccStart' = NULL,'outbreak2Start' = NULL,
                                                   'simulationEnd' = NULL)
                              
                              .self$ageGroups <- list(.self$youngSim, .self$adultSim, .self$elderlySim)
                              .self$days <- days
                              .self$outbreak2Start <- outbreak2Start
                              .self$vaccineStart <- vaccineStart
                              
                            },
                            seihrvdModel = function(time,IVS, state){
                              IVSlength <- length(IVS)
                              IVSyoung <- IVS[1:floor(IVSlength/3)]
                              IVSadult <- IVS[(floor(IVSlength/3) + 1):floor(2*IVSlength/3)]
                              IVSelderly <- IVS[(floor(2*IVSlength/3) + 1):IVSlength]
                              
                              ageIVS <- list(IVSyoung, IVSadult, IVSelderly)
                              
                              diffEqs <- c()
                              for(j in 1 : length(.self$ageGroups)){
                                age <- .self$ageGroups[[j]]
                                FOI1 <- 0
                                FOI2 <- 0
                                for(i in 1 : length(.self$ageGroups)){
                                  b1 <- age$betaFunc(time, 'beta1', i)
                                  b2 <- age$betaFunc(time, 'beta2', i)
                                  FOI1 <- FOI1 + b1*(age$params$tau1*(ageIVS[[i]][12] + ageIVS[[i]][15]) + ageIVS[[i]][11] + ageIVS[[i]][13] + ageIVS[[i]][14] + ageIVS[[i]][16] + age$params$rho1*(ageIVS[[i]][26] + ageIVS[[i]][27]))
                                  FOI2 <- FOI2 + b2*(age$params$tau2*(ageIVS[[i]][18] + ageIVS[[i]][21] + ageIVS[[i]][24]) + ageIVS[[i]][17] + ageIVS[[i]][19] + ageIVS[[i]][20] + ageIVS[[i]][22] + ageIVS[[i]][23] + ageIVS[[i]][25] + age$params$rho2*(ageIVS[[i]][28] + ageIVS[[i]][29] + ageIVS[[i]][30]))
                                }
                                diffEqs <- append(diffEqs, age$defineVars(ageIVS[[j]], FOI1, FOI2, state[[j]]))
                              }
                              return(diffEqs)
                            },
                            ODEsolver = function(){
                              ageInitialValues <- lapply(.self$ageGroups, function(age) age$initialVals)
                              currentIVS <- unlist(ageInitialValues, use.names = FALSE)
                              
                              timeAtLastEvent <- 0
                              
                              elderlyVaccineBit <- 0
                              adultVaccineBit <- 0
                              youngVaccineBit <- 0
                              outbreak2Bit <- 0
                              
                              state <- list(1, 1, 1)
                              
                              eventActive <- list(elderlyVaccine = TRUE, adultVaccine = TRUE, youngVaccine = TRUE, outbreak2 = TRUE)
                              
                              valuesArray <- list()
                              adultVaccIndices <- c(82, 83, 84, 85, 87, 89, 90, 94, 95, 96, 100, 101, 102, 103, 104, 105, 107, 109, 110, 112, 114, 115, 117, 119, 120)
                              youngVaccIndices <- c(42, 43, 44, 45, 47, 49, 50, 54, 55, 56, 60, 61, 62, 63, 64, 65, 67, 69, 70, 72, 74, 75, 77, 79, 80)
                              outbreak2Indices <- c(8, 17, 18, 48, 57, 58, 88, 97, 98)
                              
                              gradFunc <- function(t, y, params){
                                list(.self$seihrvdModel(t, y, state))
                              }
                              
                              rootFunc <- function(t, y, params){
                                elderlyVaccRoot <- if(eventActive$elderlyVaccine){
                                  t - .self$vaccineStart
                                } else 1
                                
                                adultVaccRoot <- if(eventActive$adultVaccine){
                                  totalElderlyVaccinated <- sum(y[adultVaccIndices])
                                  targetElderlyVaccinated <- 0.4 * .self$elderlySim$params$N
                                  totalElderlyVaccinated - targetElderlyVaccinated
                                } else 1
                                
                                youngVaccRoot <- if(eventActive$youngVaccine){
                                  totalAdultsVaccinated <- sum(y[youngVaccIndices])
                                  targetAdultsVaccinated <- 0.4 * .self$adultSim$params$N
                                  totalAdultsVaccinated - targetAdultsVaccinated
                                } else 1
                                
                                outbreak2Root <- if(eventActive$outbreak2){
                                  t - .self$outbreak2Start
                                } else 1
                                
                                return(c(elderlyVaccRoot, adultVaccRoot, youngVaccRoot, outbreak2Root))
                              }
                              repeat{
                                times <- seq(timeAtLastEvent, .self$days, by = 1)
                                if(length(times) < 2) break
                                out <- deSolve::ode(y = currentIVS, times = times, func = gradFunc, parms = NULL, rootfun = rootFunc, method = "lsodar")
                                outMat <- as.matrix(out)
                                
                                valuesArray[[length(valuesArray) + 1]] <- outMat[-1, -1, drop = FALSE]
                                
                                newTime <- outMat[nrow(outMat), 1]
                                currentIVS <- outMat[nrow(outMat), -1]
                                timeAtLastEvent <- newTime
                                
                                tRoot <- attr(out, "troot")
                                indRoot <- attr(out, "iroot")
                                
                                if(!is.null(indRoot)){
                                  if(indRoot[1] == 1 && eventActive$elderlyVaccine){
                                    .self$events[["elderlyVaccStart"]] <- tRoot
                                    elderlyVaccineBit <- 1
                                    eventActive$elderlyVaccine <- FALSE
                                  }
                                  if(indRoot[2] == 1 && eventActive$adultVaccine){
                                    .self$events[["adultVaccStart"]] <- tRoot
                                    adultVaccineBit <- 1
                                    eventActive$adultVaccine <- FALSE
                                  }
                                  if(indRoot[3] == 1 && eventActive$youngVaccine){
                                    .self$events[["youngVaccStart"]] <- tRoot
                                    youngVaccineBit <- 1
                                    eventActive$youngVaccine <- FALSE
                                  }
                                  if(indRoot[4] == 1 && eventActive$outbreak2){
                                    .self$events[["outbreak2Start"]] <- tRoot
                                    outbreak2Bit <- 1
                                    eventActive$outbreak2 <- FALSE
                                    currentIVS[outbreak2Indices] <- currentIVS[outbreak2Indices] + 100
                                  }
                                }
                                
                                
                                youngState <- bitwOr(bitwShiftL(outbreak2Bit, 1), youngVaccineBit) + 1
                                adultState <- bitwOr(bitwShiftL(outbreak2Bit, 1), adultVaccineBit) + 1
                                elderlyState <- bitwOr(bitwShiftL(outbreak2Bit, 1), elderlyVaccineBit) + 1
                                
                                state <- list(youngState, adultState, elderlyState)
                                if(timeAtLastEvent >= .self$days || is.null(indRoot)){
                                  break
                                }
                              }
                              allRows <- do.call(rbind, valuesArray)
                              nComp <- ncol(allRows) / 3
                              
                              valuesYoung   <- allRows[, 1:nComp, drop = FALSE]
                              valuesAdult   <- allRows[, (nComp + 1):(2 * nComp), drop = FALSE]
                              valuesElderly <- allRows[, (2 * nComp + 1):(3 * nComp), drop = FALSE]
                              
                              totalValues <- list(valuesYoung, valuesAdult, valuesElderly)
                              totalAdderConsec <- function(j, lower, upper){
                                mat <- totalValues[[j]]
                                rowSums(mat[, lower:upper, drop = FALSE])
                              }
                              totalAdderIndices <- function(j, indices){
                                mat <- totalValues[[j]]
                                rowSums(mat[, indices, drop = FALSE])
                              }
                              for(i in 1:length(.self$ageGroups)){
                                age <- .self$ageGroups[[i]]
                                
                                totalSusceptible1 <- totalAdderConsec(i, 1, 3)
                                totalSusceptible2 <- totalAdderConsec(i, 4, 5)
                                totalExposed1 <- totalAdderConsec(i, 6, 7)
                                totalExposed2 <- totalAdderConsec(i, 8, 10)
                                totalInfected1 <- totalAdderConsec(i, 11, 16)
                                totalInfected2 <- totalAdderConsec(i, 17, 25)
                                totalHospitalised1 <- totalAdderConsec(i, 26, 27)
                                totalHospitalised2 <- totalAdderConsec(i, 28, 30)
                                totalRecovered1 <- totalAdderConsec(i, 31, 32)
                                totalRecovered2 <- totalAdderConsec(i, 33, 35)
                                totalDead1 <- totalAdderConsec(i, 36, 37)
                                totalDead2 <- totalAdderConsec(i, 38, 40)
                                totalVaccinated <- totalAdderIndices(i, c(2, 3, 4, 5, 7, 9, 10, 14, 15, 16, 20, 21, 22, 23, 24, 25, 27, 29, 30, 32, 34, 35, 37, 39, 40))
                                
                                totals <- list(totalSusceptible1, totalSusceptible2, totalExposed1, totalExposed2, totalInfected1, totalInfected2, totalHospitalised1, totalHospitalised2, totalRecovered1, totalRecovered2, totalDead1, totalDead2, totalVaccinated)
                                
                                vaccinatedImmune <- totalAdderIndices(i, c(2,4,9,22,23,24,29,34,39))
                                
                                totalInfectiousS1 <- totalAdderIndices(i, c(10,13))
                                dailyHosp1 <- age$params$gamma1prime * age$params$h1 * totalInfectiousS1
                                
                                totalInfectiousS2 <- totalAdderIndices(i, c(16,19,22))
                                dailyHosp2 <- age$params$gamma2prime * age$params$h2 * totalInfectiousS2
                                
                                dailyHosp <- list(dailyHosp1, dailyHosp2)
                                
                                .self$outputs[[i]] <- list(totals, vaccinatedImmune, dailyHosp)
                              }
                              return(.self$outputs)
                            }
                          )
)
