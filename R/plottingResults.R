

.dataFramedResults <- function(results){
  youngResults <- results[[1]]
  adultResults <- results[[2]]
  elderlyResults <- results[[3]]
  
  susceptiblePlotting <- data.frame(
    time = rep(seq_len(length(youngResults[[1]][[1]])), times = 6),
    value = c(youngResults[[1]][[1]], youngResults[[1]][[2]],
              adultResults[[1]][[1]], adultResults[[1]][[2]],
              elderlyResults[[1]][[1]], elderlyResults[[1]][[2]]),
    variant_age = rep(c("0-19 yrs, variant 1", "0-19 yrs, variant 2",
                        "20-64 yrs, variant 1", "20-64 yrs, variant 2",
                        "65+ yrs, variant 1", "65+ yrs, variant 2"), each = length(youngResults[[1]][[1]]))
  )
  exposedPlotting <- data.frame(
    time = rep(seq_len(length(youngResults[[1]][[3]])), times = 6),
    value = c(youngResults[[1]][[3]], youngResults[[1]][[4]],
              adultResults[[1]][[3]], adultResults[[1]][[4]],
              elderlyResults[[1]][[3]], elderlyResults[[1]][[4]]),
    variant_age = rep(c("0-19 yrs, variant 1", "0-19 yrs, variant 2",
                        "20-64 yrs, variant 1", "20-64 yrs, variant 2",
                        "65+ yrs, variant 1", "65+ yrs, variant 2"), each = length(youngResults[[1]][[3]]))
  )
  
  infectedPlotting <- data.frame(
    time = rep(seq_len(length(youngResults[[1]][[5]])), times = 6),
    value = c(youngResults[[1]][[5]], youngResults[[1]][[6]],
              adultResults[[1]][[5]], adultResults[[1]][[6]],
              elderlyResults[[1]][[5]], elderlyResults[[1]][[6]]),
    variant_age = rep(c("0-19 yrs, variant 1", "0-19 yrs, variant 2",
                        "20-64 yrs, variant 1", "20-64 yrs, variant 2",
                        "65+ yrs, variant 1", "65+ yrs, variant 2"), each = length(youngResults[[1]][[5]]))
  )
  
  hospitalisedPlotting <- data.frame(
    time = rep(seq_len(length(youngResults[[1]][[7]])), times = 6),
    value = c(youngResults[[1]][[7]], youngResults[[1]][[8]],
              adultResults[[1]][[7]], adultResults[[1]][[8]],
              elderlyResults[[1]][[7]], elderlyResults[[1]][[8]]),
    variant_age = rep(c("0-19 yrs, variant 1", "0-19 yrs, variant 2",
                        "20-64 yrs, variant 1", "20-64 yrs, variant 2",
                        "65+ yrs, variant 1", "65+ yrs, variant 2"), each = length(youngResults[[1]][[7]]))
  )
  
  recoveredPlotting <- data.frame(
    time = rep(seq_len(length(youngResults[[1]][[9]])), times = 6),
    value = c(youngResults[[1]][[9]], youngResults[[1]][[10]],
              adultResults[[1]][[9]], adultResults[[1]][[10]],
              elderlyResults[[1]][[9]], elderlyResults[[1]][[10]]),
    variant_age = rep(c("0-19 yrs, variant 1", "0-19 yrs, variant 2",
                        "20-64 yrs, variant 1", "20-64 yrs, variant 2",
                        "65+ yrs, variant 1", "65+ yrs, variant 2"), each = length(youngResults[[1]][[9]]))
  )
  
  deadPlotting <- data.frame(
    time = rep(seq_len(length(youngResults[[1]][[11]])), times = 6),
    value = c(youngResults[[1]][[11]], youngResults[[1]][[12]],
              adultResults[[1]][[11]], adultResults[[1]][[12]],
              elderlyResults[[1]][[11]], elderlyResults[[1]][[12]]),
    variant_age = rep(c("0-19 yrs, variant 1", "0-19 yrs, variant 2",
                        "20-64 yrs, variant 1", "20-64 yrs, variant 2",
                        "65+ yrs, variant 1", "65+ yrs, variant 2"), each = length(youngResults[[1]][[11]]))
  )
  
  vaccinatedPlotting <- data.frame(
    time = rep(seq_len(length(youngResults[[1]][[13]])), times = 3),
    value = c(youngResults[[1]][[13]], adultResults[[1]][[13]], elderlyResults[[1]][[13]]),
    variant_age = rep(c("0-19 yrs1","20-64 yrs", "65+ yrs"), each = length(youngResults[[1]][[13]]))
  )
  return(list(susceptiblePlotting, exposedPlotting, infectedPlotting, hospitalisedPlotting, recoveredPlotting, deadPlotting,  vaccinatedPlotting))
}

plotter<- function(results, colourPal, listOfCompartments){
  dfResults <- .dataFramedResults(results)
  if("S" %in% listOfCompartments){
    susceptiblePlotting <- dfResults[[1]]
    susceptiblePlot <- ggplot2::ggplot(susceptiblePlotting, ggplot2::aes(x = time, y = value, colour = variant_age)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::scale_colour_manual(values = colourPal) +
      ggplot2::labs(x = "Time (days)", y = "Susceptibles", colour = "Variant and age group") +
      ggplot2::annotate("rect", xmin=lockdown1Start, xmax=lockdown1End, ymin=-Inf, ymax=Inf, fill = alpha("#88E788", 0.2)) +
      ggplot2::annotate("rect", xmin=lockdown2Start, xmax=lockdown2End, ymin=-Inf, ymax=Inf, fill = alpha("#88E788", 0.2)) +
      ggplot2::geom_vline(xintercept = vaccineStart, colour = "black", linetype = "dashed")
    print(susceptiblePlot)
  }
  if("E" %in% listOfCompartments){
    exposedPlotting <- dfResults[[2]]
    exposedPlot <- ggplot2::ggplot(exposedPlotting, ggplot2::aes(x = time, y = value, colour = variant_age)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::scale_colour_manual(values = colourPal) +
      ggplot2::labs(x = "Time (days)", y = "Exposed", colour = "Variant and age group") +
      ggplot2::annotate("rect", xmin=lockdown1Start, xmax=lockdown1End, ymin=-Inf, ymax=Inf, fill = alpha("#88E788", 0.2)) +
      ggplot2::annotate("rect", xmin=lockdown2Start, xmax=lockdown2End, ymin=-Inf, ymax=Inf, fill = alpha("#88E788", 0.2)) +
      ggplot2::geom_vline(xintercept = vaccineStart, colour = "black", linetype = "dashed")
    print(exposedPlot)
  }
  if("I" %in% listOfCompartments){
  infectedPlotting <- dfResults[[3]]
  infectedPlot <- ggplot2::ggplot(infectedPlotting, ggplot2::aes(x = time, y = value, colour = variant_age)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_colour_manual(values = colourPal) +
    ggplot2::labs(x = "Time (days)", y = "Infected", colour = "Variant and age group") +
    ggplot2::annotate("rect", xmin=lockdown1Start, xmax=lockdown1End, ymin=-Inf, ymax=Inf, fill = alpha("#88E788", 0.2)) +
    ggplot2::annotate("rect", xmin=lockdown2Start, xmax=lockdown2End, ymin=-Inf, ymax=Inf, fill = alpha("#88E788", 0.2)) +
    ggplot2::geom_vline(xintercept = vaccineStart, colour = "black", linetype = "dashed")
  print(infectedPlot)
  }
  if("H" %in% listOfCompartments){
    hospitalisedPlotting <- dfResults[[4]]
    hospitalisedPlot <- ggplot2::ggplot(hospitalisedPlotting, ggplot2::aes(x = time, y = value, colour = variant_age)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::scale_colour_manual(values = colourPal) +
      ggplot2::labs(x = "Time (days)", y = "Hospitalisations", colour = "Variant and age group") +
      ggplot2::annotate("rect", xmin=lockdown1Start, xmax=lockdown1End, ymin=-Inf, ymax=Inf, fill = alpha("#88E788", 0.2)) +
      ggplot2::annotate("rect", xmin=lockdown2Start, xmax=lockdown2End, ymin=-Inf, ymax=Inf, fill = alpha("#88E788", 0.2)) +
      ggplot2::geom_vline(xintercept = vaccineStart, colour = "black", linetype = "dashed")
    print(hospitalisedPlot)
  }
  if("R" %in% listOfCompartments){
    recoveredPlotting <- dfResults[[5]]
    recoveredPlot <- ggplot2::ggplot(recoveredPlotting, ggplot2::aes(x = time, y = value, colour = variant_age)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::scale_colour_manual(values = colourPal) +
      ggplot2::labs(x = "Time (days)", y = "Recovered", colour = "Variant and age group") +
      ggplot2::annotate("rect", xmin=lockdown1Start, xmax=lockdown1End, ymin=-Inf, ymax=Inf, fill = alpha("#88E788", 0.2)) +
      ggplot2::annotate("rect", xmin=lockdown2Start, xmax=lockdown2End, ymin=-Inf, ymax=Inf, fill = alpha("#88E788", 0.2)) +
      ggplot2::geom_vline(xintercept = vaccineStart, colour = "black", linetype = "dashed")
    print(recoveredPlot)
  }
  if("D" %in% listOfCompartments){
    deadPlotting <- dfResults[[6]]
    deadPlot <- ggplot2::ggplot(deadPlotting, ggplot2::aes(x = time, y = value, colour = variant_age)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::scale_colour_manual(values = colourPal) +
      ggplot2::labs(x = "Time (days)", y = "Dead", colour = "Variant and age group") +
      ggplot2::annotate("rect", xmin=lockdown1Start, xmax=lockdown1End, ymin=-Inf, ymax=Inf, fill = alpha("#88E788", 0.2)) +
      ggplot2::annotate("rect", xmin=lockdown2Start, xmax=lockdown2End, ymin=-Inf, ymax=Inf, fill = alpha("#88E788", 0.2)) +
      ggplot2::geom_vline(xintercept = vaccineStart, colour = "black", linetype = "dashed")
    print(deadPlot)
  }
  if("V" %in% listOfCompartments){
    vaccinatedPlotting <- dfResults[[7]]
    vaccinatedPlot <- ggplot2::ggplot(vaccinatedPlotting, ggplot2::aes(x = time, y = value, colour = variant_age)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::scale_colour_manual(values = colourPal) +
      ggplot2::labs(x = "Time (days)", y = "Vaccinated", colour = "Variant and age group") +
      ggplot2::annotate("rect", xmin=lockdown1Start, xmax=lockdown1End, ymin=-Inf, ymax=Inf, fill = alpha("#88E788", 0.2)) +
      ggplot2::annotate("rect", xmin=lockdown2Start, xmax=lockdown2End, ymin=-Inf, ymax=Inf, fill = alpha("#88E788", 0.2)) +
      ggplot2::geom_vline(xintercept = vaccineStart, colour = "black", linetype = "dashed")
    print(vaccinatedPlot)
  }
}