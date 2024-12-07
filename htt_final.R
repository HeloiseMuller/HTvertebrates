# Packages
######
library(vioplot)
library(brms)
library(rethinking)
library(phytools)
library(tidyverse)
library(treeio)
library(ggtree)
library(ggnewscale)
library(gridExtra)
######

# Functions
######
# produces a bar plot showing the predicted number of HTT for different groups of species
barplot_comp <- function(meta, orders, order_colours, data_vector, title, xlab, seed, species) {
  order_vector <- c()
  order_n <- c()
  # orders is a vector with the clades to be considered
  # loop over the clades
  # order_vector will have the group names
  # with each group name repeated as many times as there are species
  # order_n will have the corresponding number of transfers (already provided
  # with all the necessary transformations done)
  for (order in orders) {
    curr_species <- meta$species[meta$group == order]
    curr_species <- curr_species[curr_species %in% species]
    order_vector <- c(order_vector, rep(order, length(curr_species)))
    order_n <- c(order_n, data_vector[species %in% curr_species])
  }
  # calculate the median number of transfers per group so that the groups could be sorted
  # by this value
  order_medians <- aggregate(order_n ~ order_vector, FUN = median)
  ordered_orders <- order(order_medians[,2], decreasing = FALSE)
  orders_sorted <- order_medians[,1][ordered_orders]
  # make a vector of HTT numbers (n_clean), of group names (orders_clean)
  # and of colours (col_vector), all in the final order of display
  n_clean <- c()
  orders_clean <- c()
  col_vector <- c()
  for (order in orders_sorted) {
    n_clean <- c(n_clean, order_n[order_vector == order])
    orders_clean <- c(orders_clean, rep(order, sum(order_vector == order)))
    col_vector <- c(col_vector, rep(order_colours[[order]], sum(order_vector == order)))
  }
  # set up plotting area
  par("xpd" = TRUE)
  default <- par("mar")
  par("mar" = c(5.1, 8, 4.1, 2.1))
  # draw bar plot
  bp <- barplot(n_clean, horiz = TRUE, col = col_vector, border = NA, 
                xlab = xlab, main = title, cex.main = 1.5)
  # for each group, draw a vertical line and write the name of the group
  for (order in 1:length(orders_sorted)) {
    all_points <- bp[orders_clean == orders[order]]
    segments(-sd(n_clean)/25, min(all_points), -sd(n_clean)/25, max(all_points), col = order_colours[[orders[[order]]]])
    if (length(all_points) >= 3) {
      text(x = -sd(n_clean) * 0.8, y = mean(c(min(all_points), max(all_points))), 
           labels = orders[order],
           col = order_colours[[orders[[order]]]], font = 2, cex = 0.85)
    }
  }
  par("mar" = default)
}

# fit a negative binomial regression model of the number of HTT (for one particular species)
fit_model <- function(curr, agg, agg_function, time_mean, time_sd, habitat = FALSE, habitat_control = FALSE, alternative = FALSE, iter = 2000) {
  # by default, model the raw number of HTT ("n")
  if (alternative == FALSE) {
    alternative <- "n"
  }
  # these are the outputs that will be returned from this function
  lower <- FALSE
  upper <- FALSE
  sample_size <- FALSE
  ppd <- FALSE
  test_times <- FALSE
  # if the partner species are to be averaged by spClade
  if (agg == TRUE) {
    # calculate median HTT numbers by the scaled divergence time, and the spClade and habitat of the partner species
    if (habitat == TRUE) {
      curr_temp <- aggregate(curr[, alternative] ~ curr[,"zdivTime"] + curr[,"spClade.2"] + curr[,"same_habitat"], FUN = agg_function)
      columns <- c("zdivTime", "spClade.2", "habitat", "n", "divTime")
      # add the raw divergence time to the data
      curr_temp2 <- aggregate(curr[, alternative] ~ curr[,"divTime"] + curr[,"spClade.2"] + curr[,"same_habitat"], FUN = agg_function)
    }
    else {
      # don't consider habitat
      curr_temp <- aggregate(curr[, alternative] ~ curr[,"zdivTime"] + curr[,"spClade.2"], FUN = agg_function)
      columns <- c("zdivTime", "spClade.2", "n", "divTime")
      curr_temp2 <- aggregate(curr[, alternative] ~ curr[,"divTime"] + curr[,"spClade.2"], FUN = agg_function)
    }
    curr <- cbind(curr_temp, curr_temp2[,1])
    colnames(curr) <- columns
    # after taking the median, the number of HTT may have decimal points so round to the
    # nearest integer
    curr$n <- round(curr$n)
  }
  # also returned by the function
  lower <- NA 
  upper <- NA
  ppd <- NA 
  expected_500 <- NA 
  expected_1300 <- NA
  sample_size <- NA
  test_times <- NA
  intercept <- NA
  coefficient <- NA
  habitat_same <- NA
  habitat_diff <- NA
  confidence <- NA
  # don't fit model if no HTT were observed
  if (sum(curr$n) > 0) {
    # set priors
    # effect of phylogenetic proximity
    prior_prox <- set_prior("normal(0,1)", class = "b", coef = "zdivTime")
    # effect of shared habitat
    prior_habitat <- set_prior("normal(0,3)", class = "b", coef = "habitatyes")
    prior_i <- set_prior("normal(0,10)", class = "Intercept")
    prior_shape <- set_prior("gamma(0.001,0.001)", class = "shape",
                             lb = 0)
    # model that includes habitat
    if (habitat == TRUE) {
      column_name <- "b_habitatyes"
      model <- brm(n ~ zdivTime + habitat, data = curr, family = negbinomial(),
                   iter = iter, control = list(adapt_delta = 0.99),
                   backend = "cmdstanr", cores = 4,
                   prior = c(prior_prox, prior_habitat, prior_i, prior_shape), silent = 2, refresh = 0,
                   warmup = 1000)
    }
    # model without habitat
    else {
      column_name <- "b_zdivTime"
      model <- brm(n ~ zdivTime, data = curr, family = negbinomial(),
                   iter = iter, control = list(adapt_delta = 0.99),
                   backend = "cmdstanr", cores = 4,
                   prior = c(prior_prox, prior_i, prior_shape), silent = 2, refresh = 0,
                   warmup = 1000)
    }
    if (habitat_control == TRUE) {
      column_name <- "b_zdivTime"
    }
    # check Rhats
    rhats <- rhat(model)
    if (max(rhats) > 1.01) {
      print("Warning! High Rhat!")
      print(rhats)
    }
    # check ESS
    # calculate number of postwarmup draws
    draws <- (iter - 1000) * 4
    ess <- neff_ratio(model) * draws
    if (min(ess) < 2000) {
      print("Warning! Low ESS!")
      print(ess)
    }
    # get the MCMC draws
    samples <- as.matrix(model)
    # get regression coefficient of divergence time/habitat and the intercept
    coefficient <- median(samples[, column_name])
    intercept <- median(samples[, "b_Intercept"])
    # get lower and upper bound of HPDI (either for the effect of the divergence time or of the habitat)
    hpdi <- HPDI(samples[, column_name], 0.95)
    lower <- as.numeric(hpdi[1])
    upper <- as.numeric(hpdi[2])
    confidence <- mean(samples[, column_name] > 0)
    # to get minimum divTime
    times <- min(curr$divTime)
    # reduce to a multiple of 100
    minimum <- min(times) - min(times)%%100
    # these are the divergence times that we will make predictions for
    test_times <- (seq(minimum, 1400, by = 100) - time_mean)/time_sd
    # if we are studying habitats, we also want to make predictions for the two habitats 
    if (habitat == TRUE) {
      test_data <- data.frame(zdivTime = rep(test_times, 2), habitat = rep(c("yes", "no"), each = length(test_times)))
      test_times <- rep(test_times, 2) 
    }
    else {
      test_data <- data.frame(zdivTime = test_times)
      
    }
    # calculate predicted mean number of HTT for each divergence time (and potentially habitat)
    # based on each retained MCMC draw
    ppd <- posterior_epred(model, newdata = test_data)
    # predictions for 500 and 1300 MYA for other visualizations
    test_times_two <- (c(500, 1300) - time_mean)/time_sd
    if (habitat == TRUE) {
      ppd_three <- posterior_epred(model, newdata = data.frame(zdivTime = rep(test_times_two, 2), habitat = rep(c("yes", "no"), each = length(test_times_two))))
    }
    else {
      ppd_three <- posterior_epred(model, newdata = data.frame(zdivTime = test_times_two))
    }
    # take median across MCMC draws
    expected_500 <- median(ppd_three[,1])
    expected_1300 <- median(ppd_three[,2])
    sample_size <- dim(curr)[1]
    if (habitat == TRUE) {
      habitat_same <- median(ppd_three[,1])
      habitat_diff <- median(ppd_three[,3])
    }
  }
  return(list("lower"=lower, "upper"=upper, "ppd"=ppd, "expected_500"=expected_500, "expected_1300"=expected_1300, "sample_size"=sample_size, "test_times"=test_times, "intercept"=intercept, "coefficient"=coefficient, "habitat_same" = habitat_same, "habitat_diff"=habitat_diff, "confidence" = confidence))
}

# fit a negative binomial regression model of the number of HTT (for one particular species),
# for Class I/Class II elements separately and visualize
# see function fit_model() for more detailed comments
fit_model_classes <- function(HTT, species, subsample = FALSE) {
  # prepare colours for the plot
  desert <- rgb(1, 0.89, 0.77, alpha = 0.5)
  lightcoral <- rgb(250/255, 168/255, 68/255, alpha = 0.5)
  lightblue <- rgb(0.68, 0.85, 0.90, alpha = 0.7)
  pink <- rgb(255, 182, 193, alpha = 255/2, maxColorValue = 255)
  # calculate mean and SD of the divergence time so you could convert to and from Z-scores
  time_mean <- mean(HTT$divTime)
  time_sd <- sd(HTT$divTime)
  curr <- HTT[HTT$species.1 == species,]
  # assemble data for eithe class separately
  columns <- c("zdivTime", "spClade.2", "habitat", "n", "divTime")
  curr_temp <- aggregate(curr[, "type1"] ~ curr[,"zdivTime"] + curr[,"spClade.2"] + curr[,"same_habitat"], FUN = median)
  curr_temp2 <- aggregate(curr[, "type1"] ~ curr[,"divTime"] + curr[,"spClade.2"] + curr[,"same_habitat"], FUN = median)
  curr1 <- cbind(curr_temp, curr_temp2[,1])
  colnames(curr1) <- columns
  print(curr1)
  curr1$n <- round(curr1$n)
  curr_temp <- aggregate(curr[, "type2"] ~ curr[,"zdivTime"] + curr[,"spClade.2"] + curr[,"same_habitat"], FUN = median)
  curr_temp2 <- aggregate(curr[, "type2"] ~ curr[,"divTime"] + curr[,"spClade.2"] + curr[,"same_habitat"], FUN = median)
  curr2 <- cbind(curr_temp, curr_temp2[,1])
  colnames(curr2) <- columns
  curr2$n <- round(curr2$n)
  sample_string <- ""
  main <- gsub("_", " ", species)
  # if you want to subsample the Class 2 elements to be of the same number as Class 1
  if (subsample == TRUE) {
    probs <- curr2$n/sum(curr2$n)
    sampled <- sample(curr2$spClade.2, size = sum(curr1$n), prob = probs, replace = TRUE)
    new_n <- c()
    for (spc2 in curr2$spClade.2) {
      new_n <- c(new_n, sum(sampled == spc2))
    }
    curr2$n <- new_n
    sample_string <- "_subsampled"
    main <- paste(species, " (subsampled)", sep = "")
  }
  # fit one model per class
  # effect of phylogenetic proximity
  prior_prox <- set_prior("normal(0,1)", class = "b", coef = "zdivTime")
  # effect of shared habitat
  prior_habitat <- set_prior("normal(0,3)", class = "b", coef = "habitatyes")
  prior_i <- set_prior("normal(0,10)", class = "Intercept")
  prior_shape <- set_prior("gamma(0.001,0.001)", class = "shape",
                           lb = 0)
  column_name <- "b_zdivTime"
  model1 <- brm(n ~ zdivTime + habitat, data = curr1, family = negbinomial(),
                 iter = 10000, control = list(adapt_delta = 0.99),
                 backend = "cmdstanr", cores = 4,
                 prior = c(prior_prox, prior_habitat, prior_i, prior_shape), silent = 2, refresh = 0,
                warmup = 1000)
  model2 <- brm(n ~ zdivTime + habitat, data = curr2, family = negbinomial(),
                iter = 10000, control = list(adapt_delta = 0.99),
                backend = "cmdstanr", cores = 4,
                prior = c(prior_prox, prior_habitat, prior_shape), silent = 2, refresh = 0,
                warmup = 1000)
  # predict for the same divergence times, for the two classes separately
  times <- min(c(curr1$divTime, curr2$divTime))
  minimum <- min(times) - min(times)%%100
  test_times_raw <- seq(minimum, 1400, by = 100) 
  test_times <- (test_times_raw - time_mean)/time_sd
  test_data <- data.frame(zdivTime = test_times, habitat = "yes")
  ppd1 <- posterior_epred(model1, newdata = test_data)
  ppd2 <- posterior_epred(model2, newdata = test_data)
  # open plot
  pdf(paste("figures/", species, "_class_comp", sample_string, ".pdf", sep = ""), width = 5.7, height = 5)
  par(mfrow = c(1,1))
  maximum <- max(c(max(curr1$n), max(curr2$n)) + 1)
  # make empty plot just to create the plot area
  plot(ppd1[100,] ~ test_times_raw, pch = "", xlab = "Divergence time (MY)",
       ylab = "Predicted # of HTT", ylim = c(0, maximum),
       main = main, cex.lab = 1, cex.axis = 1, cex.main = 1.5)
  # add prediction lines from the first 1000 MCMC draws
  for (draw in 1:1000) {
    lines(ppd2[draw,] ~ test_times_raw, lwd = 2, col = pink)
    lines(ppd1[draw,] ~ test_times_raw, lwd = 2, col = lightblue)
  }
  # add points with the observed data
  # plot out Class 1 and Class 2 in alternation so one wouldn't
  # cover up the other
  breakpoints <- as.integer(seq(5, dim(curr1)[1], length.out = 35))
  start <- 1
  # to make sure even empty levels get plotted out
  for (bp in breakpoints) {
    divtimes <- curr1[start:bp,]$divTime + rnorm(bp - start + 1, 0, 4)
    points(curr1[start:bp,]$n ~ divtimes, pch = 20,
           col = rgb(0, 0, 255, 60, maxColorValue = 255))
    divtimes <- curr2[start:bp,]$divTime + rnorm(bp - start + 1, 0, 4)
    points(curr2[start:bp,]$n ~ divtimes, pch = 20,
      col = rgb(255, 0, 0, 80, maxColorValue = 255))
    start <- bp
  }
  # add legend
  legend("topright", pch = 19, col = c("blue", "red"),
         legend = c("Class 1", "Class 2"), cex = 1)
  dev.off()  
}

# fit the model to just a single species and visualize. To be used when there is no wrapper.
# see fit_model() and fit_model_classes() for more detailed comments
fit_model_single <- function(HTT, species) {
  desert <- rgb(1, 0.89, 0.77, alpha = 0.5)
  time_mean <- mean(HTT$divTime)
  time_sd <- sd(HTT$divTime)
  # filter to only the species of interest
  curr <- HTT[HTT$species.1 == species,]
  curr_temp <- aggregate(curr[, "n"] ~ curr[,"zdivTime"] + curr[,"spClade.2"] + curr[,"same_habitat"], FUN = median)
  columns <- c("zdivTime", "spClade.2", "habitat", "n", "divTime")
  # add the raw divergence time to the data
  curr_temp2 <- aggregate(curr[, "n"] ~ curr[,"divTime"] + curr[,"spClade.2"] + curr[,"same_habitat"], FUN = median)
  curr <- cbind(curr_temp, curr_temp2[,1])
  colnames(curr) <- columns
  curr$n <- round(curr$n)
  print("% above 0 at <1000 Myrs divergence:")
  print(mean(curr$n[curr$divTime < 1000] > 0) * 100)
  print("% above 0 at >1000 Myrs divergence:")
  print(mean(curr$n[curr$divTime > 1000] > 0) * 100)
  # in species name, replace underscore with space
  main <- gsub("_", " ", species) 
  # effect of phylogenetic proximity
  prior_prox <- set_prior("normal(0,1)", class = "b", coef = "zdivTime")
  # effect of shared habitat
  prior_habitat <- set_prior("normal(0,3)", class = "b", coef = "habitatyes")
  prior_i <- set_prior("normal(0,10)", class = "Intercept")
  prior_shape <- set_prior("gamma(0.001,0.001)", class = "shape",
                           lb = 0)
  column_name <- "b_zdivTime"
  model <- brm(n ~ zdivTime + habitat, data = curr, family = negbinomial(),
                iter = 10000, control = list(adapt_delta = 0.99),
                backend = "cmdstanr", cores = 4,
                prior = c(prior_prox, prior_habitat, prior_i, prior_shape), silent = 2, refresh = 0,
               warmup = 1000)
  print(summary(model))
  plot(model)
  # to get minimum divTime
  times <- min(curr$divTime)
  minimum <- min(times) - min(times)%%100
  test_times_raw <- seq(minimum, 1400, by = 100) 
  test_times <- (test_times_raw - time_mean)/time_sd
  test_data <- data.frame(zdivTime = test_times, habitat = "yes")
  ppd <- posterior_epred(model, newdata = test_data)
  pdf(paste("figures/", species, "_single.pdf", sep = ""), width = 5.7, height = 5)
  par(mfrow = c(1,1))
  maximum <- max(curr$n) + 1
  plot(ppd[100,] ~ test_times_raw, pch = "", xlab = "Additive divergence time (MY)",
       ylab = "Predicted # of HTT", ylim = c(0, maximum),
       main = main, cex.lab = 1, cex.axis = 1, cex.main = 1.5)
  for (draw in 1:1000) {
    lines(ppd[draw,] ~ test_times_raw, lwd = 2, col = desert)
  }
  stripchart(curr$n ~ curr$divTime, pch = 20, vertical = TRUE,
             col = rgb(112, 128, 144, 60, maxColorValue = 255), 
             method = "jitter", add = TRUE,
             at = sort(unique(curr$divTime)), jitter = 10)
  dev.off()  
  print(test_times_raw)
  print(HPDI(ppd[,1], 0.95))
  print(HPDI(ppd[,length(test_times_raw)], 0.95))
}

# wrapper around fit_model() to loop over all the species, apply fit_model()
# to each one and then make one plot per species group
fit_models_wrapper <- function(orders, meta, HTT, agg, agg_function, agg_string, with_type = FALSE, habitat = FALSE, habitat_control = FALSE, alternative = FALSE, iter = 2000) {
  # all the stuff that will be returned as output by the function (one value per species)
  negative <- c()
  positive <- c()
  species <- c()
  coefficients <- c()
  sample_sizes <- c()
  intercepts <- c()
  check1_big <- c()
  check2_big <- c()
  upper_big <- c()
  lower_big <- c()
  expected_500 <- c()
  expected_1300 <- c()
  habitat_same <- c()
  habitat_diff <- c()
  confidence <- c()
  # loop over the species groups
  for (order in orders) {
    # the stuff that will only be used to make the plot but that will not be returned as output
    print(order)
    lower <- c()
    upper <- c()
    check1 <- c()
    check2 <- c()
    # get the species in the species group
    to_check <- meta[meta$group == order, "species"]
    to_check <- to_check[to_check %in% HTT$species.1]
    # print out how many species there are in the group
    print("total:")
    print(length(to_check))
    counter <- 1
    # will contain the posterior estimates per MCMC draw
    curves <- list()
    sample_sizes_local <- c()
    species_local <- c()
    test_times <- list()
    time_mean <- mean(HTT$divTime)
    time_sd <- sd(HTT$divTime)
    # loop over the species in the species group
    for (sp in to_check) {
      print(counter)
      print(sp)
      counter <- counter + 1
      # filter data to only the current species
      # fit model
      curr <- HTT[HTT$species.1 == sp,]
      result <- fit_model(curr, agg, agg_function, time_mean, time_sd, 
                          habitat = habitat, habitat_control = habitat_control, 
                          alternative = alternative, iter = iter) 
      # store the various outputs from the model
      species <- c(species, sp)
      species_local <- c(species_local, sp)
      upper <- c(upper, result[["upper"]])
      lower <- c(lower, result[["lower"]])
      intercepts <- c(intercepts, result[["intercept"]])
      coefficients <- c(coefficients, result[["coefficient"]])
      sample_sizes <- c(sample_sizes, result[["sample_size"]])
      sample_sizes_local <- c(sample_sizes_local, result[["sample_size"]])
      habitat_same <- c(habitat_same, result["habitat_same"])
      habitat_diff <- c(habitat_diff, result["habitat_diff"])
      confidence <- c(confidence, result[["confidence"]])
      curves[[sp]] <- result[["ppd"]]
      if (length(result[["test_times"]]) == 1) {
        test_times[[sp]] <- NA
      }
      else {
        test_times[[sp]] <- (result[["test_times"]] * time_sd) + time_mean
      }
      expected_500 <- c(expected_500, result[["expected_500"]])
      expected_1300 <- c(expected_1300, result[["expected_1300"]])
      check1 <- c(check1, sum(curr$type1))
      check2 <- c(check2, sum(curr$type2))
    }
    # get species where the entire HPDI for the coefficient of interest is either above or below 0
    negative <- c(negative, species_local[upper < 0 & !is.na(upper)])
    positive <- c(positive, species_local[lower > 0 & !is.na(lower)])
    check1_big <- c(check1_big, check1)
    check2_big <- c(check2_big, check2)
    lower_big <- c(lower_big, lower)
    upper_big <- c(upper_big, upper)
    # make pretty graph
    plot_order(order, lower, upper, ppd, species_local, sample_sizes_local, 
               curves, test_times, agg, agg_function, agg_string, check1, 
               check2, with_type = with_type, habitat = habitat, 
               habitat_control = habitat_control, alternative = alternative)
  }
  return(list("positive" = positive, "negative" = negative, "intercepts" = intercepts, "coefficients" = coefficients, "species" = species, "sample_sizes" = sample_sizes, "check1" = check1_big, "check2" = check2_big, "upper" = upper_big, "lower" = lower_big, "expected_500" = expected_500, "expected_1300" = expected_1300, "habitat_same" = habitat_same, "habitat_diff" = habitat_diff, "confidence" = confidence))
}

# takes as input one vector with names of species that show a negative effect
# and another with the species that show a positive effect
# and returns a new vector that says for each species whether it is in one of the 
# two lists
get_affected_leaves <- function(pos, neg, species, binary = FALSE) {
  affected_leaves <- c()
  counter <- 0
  for (sp in species) {
    counter <- counter + 1
    if (binary == TRUE) {
      if (sp %in% c(pos, neg)) {
        affected_leaves <- c(affected_leaves, "effect present")
      }
      else {
        affected_leaves <- c(affected_leaves, "no evidence")
      }      
    }
    else {
      if (sp %in% pos) {
        affected_leaves <- c(affected_leaves, "positive")
      }
      else if (sp %in% neg) {
        affected_leaves <- c(affected_leaves, "negative")
      }
      else {
        affected_leaves <- c(affected_leaves, "no evidence")
      }      
    }
  }
  return(affected_leaves)
}

# reads and cleans up a file with HTT data
get_HTT <- function(file_name) {
  # read in file
  HTT <- read.csv(file_name, sep = "\t")
  # make a copy of the data with the columns ordered
  # differently to capture the possibility that the 
  # transfer could have been from species 2 to species 1
  # rather than 1 to 2
  HTT2 <- HTT[, c(1, 2, 4, 3, 5, 6, 7, 8, 10, 9, 11)]
  colnames(HTT2) <- colnames(HTT)
  # bind the two data frames together, so that each HTT event is
  # represented twice, once for each direction
  HTT <- rbind(HTT, HTT2)
  # scale the divergence time
  HTT$zdivTime <- scale(HTT$divTime)  
  return(HTT)
}

# fit the NB regression model to all the species, make per-group plots 
# and write the results to a series of files
model_and_write <- function(ID_string, orders, meta, HTT, agg, agg_function, agg_string, with_type = FALSE, habitat = FALSE, habitat_control = FALSE, alternative = FALSE, iter = 2000) {
  # fit models
  affected <- fit_models_wrapper(orders, meta, HTT, agg = agg, agg_function = agg_function, agg_string = agg_string, with_type = with_type, habitat = habitat, habitat_control = habitat_control, alternative = alternative, iter = iter)
  hab_string = ""
  if (habitat == TRUE) {
    hab_string <- "_habitat"
  }
  if (habitat_control == TRUE) {
    hab_string <- "_habitat_control"
  }
  # write results to files
  write.table(affected[["positive"]], file = paste(ID_string, "_positive_", agg_string, hab_string, ".txt", sep = ""), quote = FALSE,
              row.names = FALSE)
  write.table(affected[["negative"]], file = paste(ID_string, "_negative_", agg_string, hab_string, ".txt", sep = ""), quote = FALSE,
              row.names = FALSE)
  write.table(affected[["intercepts"]], file = paste(ID_string, "_intercepts_", agg_string, hab_string, ".txt", sep = ""), quote = FALSE,
              row.names = FALSE)
  write.table(affected[["coefficients"]], file = paste(ID_string, "_coefficients_", agg_string, hab_string, ".txt", sep = ""), quote = FALSE,
              row.names = FALSE)
  write.table(affected[["species"]], file = paste(ID_string, "_species_", agg_string, hab_string, ".txt", sep = ""), quote = FALSE,
              row.names = FALSE)
  write.table(affected[["sample_sizes"]], file = paste(ID_string, "_sample_sizes_", agg_string, hab_string, ".txt", sep = ""), quote = FALSE,
              row.names = FALSE)
  write.table(affected[["upper"]], file = paste(ID_string, "_upper_", agg_string, hab_string, ".txt", sep = ""), quote = FALSE,
              row.names = FALSE)
  write.table(affected[["lower"]], file = paste(ID_string, "_lower_", agg_string, hab_string, ".txt", sep = ""), quote = FALSE,
              row.names = FALSE)
  write.table(affected[["expected_500"]], file = paste(ID_string, "_expected_500_", agg_string, hab_string, ".txt", sep = ""), quote = FALSE,
              row.names = FALSE)
  write.table(affected[["expected_1300"]], file = paste(ID_string, "_expected_1300_", agg_string, hab_string, ".txt", sep = ""), quote = FALSE,
              row.names = FALSE)
  write.table(affected[["check1"]], file = paste(ID_string, "_check1_", agg_string, hab_string, ".txt", sep = ""), quote = FALSE,
                row.names = FALSE)
  write.table(affected[["check2"]], file = paste(ID_string, "_check2_", agg_string, hab_string, ".txt", sep = ""), quote = FALSE,
                row.names = FALSE) 
  write.table(unlist(affected[["habitat_same"]]), file = paste(ID_string, "_habitat_same_", agg_string, hab_string, ".txt", sep = ""), quote = FALSE,
              row.names = FALSE) 
  write.table(unlist(affected[["habitat_diff"]]), file = paste(ID_string, "_habitat_diff_", agg_string, hab_string, ".txt", sep = ""), quote = FALSE,
              row.names = FALSE)   
  write.table(unlist(affected[["confidence"]]), file = paste(ID_string, "_confidence_", agg_string, hab_string, ".txt", sep = ""), quote = FALSE,
              row.names = FALSE)   
}

# make a plot of # HTT as a function of divergence time for a specific species group
# see fit_model_classes for more detailed comments
plot_order <- function(order_name, lower, upper, ppd, species, sample_sizes, curves, test_times, agg, agg_function, agg_string, check1, check2, with_type = FALSE, habitat = FALSE, habitat_control = FALSE, alternative = FALSE) {
  desert <- rgb(1, 0.89, 0.77, alpha = 0.5)
  red_trans <- rgb(0.38, 0.16, 0.06, alpha = 0.5)
  blue_trans <- rgb(0.1, 0.1, 1, alpha = 0.5)
  # additions to the file name depending on the exact type of model
  type_string <- ""
  if (with_type == TRUE) {
    type_string <- "_with_type"
  }
  if (habitat == TRUE) {
    hab_string <- "_habitat"
    xlab <- "Effect of shared habitat"
  }
  else {
    hab_string <- ""
    xlab <- "Effect of divTime"
  }
  if (habitat_control == TRUE) {
    hab_string <- "_habitat_control"
    xlab <- "Effect of divTime"
  }
  if (alternative == FALSE) {
    alt_string <- ""
  }
  else {
    alt_string <- paste("_", alternative, "_")
  }
  # don't make a plot if none of the species had enough data to model
  species <- species[!is.na(upper)]
  if (length(species) > 0) {
    if (agg == TRUE) {
      pdf(paste("figures/", agg_string, "_", order_name, type_string, hab_string, alt_string, ".pdf", sep = ""), width = 9, height = 3 + (length(lower) * 0.1))
    }
    else {
      pdf(paste("figures/", order_name, type_string, hab_string, alt_string, ".pdf", sep = ""), width = 28, height = 3 + (length(lower) * 0.1))
    }
    par(mfrow = c(1,3))
    lower <- lower[!is.na(upper)]
    sample_sizes <- sample_sizes[!is.na(upper)]
    check1 <- check1[!is.na(upper)]
    check2 <- check2[!is.na(upper)]
    upper <- upper[!is.na(upper)]
    # order the species by the upper limit of the HPDI
    correct_order <- order(upper)
    maximum <- max(upper + 22, 22)
    shift <- 0
    if (with_type == TRUE) {
      shift <- 1
    }
    # make a plot with a segment for the 95% HPDI of the regression coefficient
    # for each of the species in the group
    plot(seq(1,length(lower)) ~ lower[correct_order], pch = "", xlim = c(min(lower) - 1, maximum),
         col = "red4", main = order_name, xlab = xlab,
         yaxt = "n", ylab = "", ylim = c(0, length(lower) + shift))
    segments(lower[correct_order], seq(1, length(lower)), upper[correct_order], seq(1, length(lower)),
             col = "black")
    # add a point to the either end of each segment because it's pretty
    points(seq(1,length(lower)) ~ upper[correct_order], pch = 19, col = "red4")
    points(seq(1,length(lower)) ~ lower[correct_order], pch = 19, col = "red4")
    # annotate with the names of the species
    text(maximum - 10, y = seq(1, length(lower)), labels = species[correct_order],
         cex = 1)
    # draw a straight vertical line at 0
    abline(v = 0, lty = 3, lwd = 2)
    # add a blue point if the species has Class I elements
    # and an orange point if it has Class II elements
    if (with_type == TRUE) {
      check1 <- check1 > 0
      check2 <- check2 > 0
      points(seq(1,length(lower))[check1[correct_order]] ~ rep(min(lower) - 0.8, sum(check1)), col = "steelblue2", pch = 19)
      points(seq(1,length(lower))[check2[correct_order]] ~ rep(min(lower) - 0.6, sum(check2)), col = "coral1", pch = 19)
      legend("topleft", legend = c("class I", "class II"), pch = 19, col = c("steelblue2", "coral1"),
             horiz = TRUE)
    }
    # pick the spcies with the most negative and the least negative effect
    lowest <- species[upper == min(upper)]
    highest <- species[upper == max(upper)]
    # make a pretty plot for either of these species (see fit_model_classes() for more detailed comments)
    if (alternative == FALSE) {
      alternative <- "n"
    }
    for (sp in c(lowest, highest)) {
      curr <- HTT[HTT$species.1 == sp,]
      if (agg == TRUE) {
        real <- aggregate(curr[, alternative] ~ curr[, "divTime"] + curr[, "spClade.2"],
                          FUN = agg_function)
      }
      else {
        real <- HTT[HTT$species.1 == sp,]
      }
      if (habitat == TRUE) {
        real <- aggregate(curr[, alternative] ~ curr[, "divTime"] + curr[, "spClade.2"] + curr[, "same_habitat"],
                          FUN = agg_function)
        colnames(real) <- c("divTime", "spClade.2", "same_habitat", "n")
      }
      maximum <- max(c(max(curves[[sp]][100,]) + 1, max(real[, dim(real)[2]])))
      plot(curves[[sp]][100,] ~ test_times[[sp]], pch = "", xlab = "Divergence time",
           ylab = "Predicted # of HTT", ylim = c(0, maximum),
           main = sp)
      for (draw in 1:1000) {
        if (habitat == TRUE) {
          original <- as.integer(length(test_times[[sp]])/2)
          lines(curves[[sp]][draw,1:original] ~ test_times[[sp]][1:original], lwd = 2, col = alpha("coral1", alpha = 0.5))
          lines(curves[[sp]][draw,(original + 1):length(test_times[[sp]])] ~ test_times[[sp]][(original + 1):length(test_times[[sp]])], lwd = 2, col = alpha("lightblue", alpha = 0.5))
        }
        else {
          lines(curves[[sp]][draw,] ~ test_times[[sp]], lwd = 2, col = desert)
        }
      }
      if (habitat == TRUE) {
        points(real[, dim(real)[2]][real$same_habitat == "yes"] ~ real[,1][real$same_habitat == "yes"], pch = 19, col = "gray1")
        points(real[, dim(real)[2]][real$same_habitat == "no"] ~ real[,1][real$same_habitat == "no"], pch = 1, col = "gray1")
      }
      else {
        points(real[, dim(real)[2]] ~ real[,1], pch = 19, col = "gray1")
      }
    }
    dev.off()
  }
}

# display results on a phylogeny (small version)
plot_trees <- function(ID_string, agg_string, tree, habitat = FALSE, habitat_control = FALSE, numbers = c()) {
  hab_string <- ""
  if (habitat == TRUE) {
    hab_string <- "_habitat"
  }
  if (habitat_control == TRUE) {
    hab_string <- "_habitat_control"
  }
  affected_pos <- read.csv(paste(ID_string, "_positive_", agg_string, hab_string, ".txt", sep = ""))$x
  affected_neg <- read.csv(paste(ID_string, "_negative_", agg_string, hab_string, ".txt", sep = ""))$x
  coefficients <- read.csv(paste(ID_string, "_coefficients_", agg_string, hab_string, ".txt", sep = ""))$x
  confidence <- read.csv(paste(ID_string, "_confidence_", agg_string, hab_string, ".txt", sep = ""))$x
  intercepts <- read.csv(paste(ID_string, "_intercepts_", agg_string, hab_string, ".txt", sep = ""))$x
  species <- read.csv(paste(ID_string, "_species_", agg_string, hab_string, ".txt", sep = ""))$x
  # prepare colours based on positive vs negative coefficients (95% HPDI)
  affected_leaves <- c()
  affected_cols <- c()
  divTime_coefs <- c()
  divTime_intercepts <- c()
  counter <- 0
  for (sp in species) {
    counter <- counter + 1
    if (sp %in% affected_pos) {
      affected_leaves <- c(affected_leaves, "positive")
      affected_cols <- c(affected_cols, "steelblue2")
    }
    else if (sp %in% affected_neg) {
      affected_leaves <- c(affected_leaves, "negative")
      affected_cols <- c(affected_cols, "red2")
    }
    else {
      affected_leaves <- c(affected_leaves, "no evidence")
      affected_cols <- c(affected_cols, "gray92")
    }
    divTime_coefs <- c(divTime_coefs, coefficients[counter])
    divTime_intercepts <- c(divTime_intercepts, intercepts[counter])
  }
  nodes <- data.frame("species" = species, 
                      "Effect_of_time" = affected_leaves,
                      "Coefficient" = divTime_coefs,
                      "intercepts" = divTime_intercepts,
                      "cols" = affected_cols,
                      "Confidence_in_pos_effect" = confidence,
                      "numbers" = numbers)  
  if (habitat == TRUE) {
    nodes <- data.frame("species" = species, 
                        "Effect_of_shared_habitat" = affected_leaves,
                        "Coefficient" = divTime_coefs,
                        "intercepts" = divTime_intercepts,
                        "cols" = affected_cols,
                        "Confidence_in_pos_effect" = confidence,
                        "numbers" = numbers)
  }
  if (habitat_control == TRUE) {
    nodes <- data.frame("species" = species, 
                        "Effect_of_time" = affected_leaves,
                        "Coefficient" = divTime_coefs,
                        "intercepts" = divTime_intercepts,
                        "cols" = affected_cols,
                        "Confidence_in_pos_effect" = confidence,
                        "numbers" = numbers)
  }
  x <- full_join(as_tibble(tree), nodes, join_by(label == species))
  tree2 <- as.treedata(x)
  pdf(paste("figures/", ID_string, hab_string, "_tree_", agg_string, ".pdf", sep = ""), width = 23, height = 20)
  if (habitat_control == TRUE) {
    print(ggtree(tree2, layout = "circular") + geom_tiplab(aes(geom = "point")) +
            geom_tippoint(pch=16, size=5, aes(col=Effect_of_time)) +
            scale_color_manual(values=c("blue", "gray92","red2")) +
            theme(legend.key.size = unit(1, 'cm'),
                  legend.title = element_text(size=40),
                  legend.text = element_text(size=40)))
    dev.off()    
  } else if (habitat == TRUE) {
    print(ggtree(tree2, layout = "circular") + geom_tiplab(aes(geom = "point")) +
            geom_tippoint(pch=16, size=5, aes(col=Effect_of_shared_habitat)) +
            scale_color_manual(values=c("blue", "gray92","red2")) +
            theme(legend.key.size = unit(1, 'cm'),
                  legend.title = element_text(size=40),
                  legend.text = element_text(size=40)))
    dev.off()    
  }
  else {
    print(ggtree(tree2, layout = "circular") + geom_tiplab(aes(geom = "point")) +
            geom_tippoint(pch=16, size=5, aes(col=Effect_of_time)) +
            scale_color_manual(values=c("blue", "gray92","red2")) +
            theme(legend.key.size = unit(1, 'cm'),
                  legend.title = element_text(size=40),
                  legend.text = element_text(size=40)))
    dev.off()   
  }
  pdf(paste("figures/", ID_string, hab_string, "_tree_coefficients_", agg_string, ".pdf", sep = ""), width = 23, height = 20)
  minimum <- -5.5
  if (min(coefficients, na.rm = TRUE) < minimum) {
    print("Warning! Minimum not low enough!")
  }
  maximum <- 5.5
  if (max(coefficients, na.rm = TRUE) > maximum) {
    print("Warning! Maximum not high enough!")
  }
  print(ggtree(tree2, layout = "circular") + geom_tiplab(aes(geom = "point", label = numbers), offset = 5) +
          geom_tippoint(pch=16, size=5, aes(col=Coefficient)) +
          theme(legend.key.size = unit(1, 'cm'),
                legend.title = element_text(size=40),
                legend.text = element_text(size=40)) + scale_color_gradient2(midpoint=0, low="blue", mid="gray92",
                                                                             high="red2", limits = c(minimum, maximum)))
  dev.off()
  pdf(paste("figures/", ID_string, hab_string, "_tree_confidence_", agg_string, ".pdf", sep = ""), width = 23, height = 20)
  minimum <- 0
  maximum <- 1
  print(ggtree(tree2, layout = "circular") + geom_tiplab(aes(geom = "point")) +
          geom_tippoint(pch=16, size=5, aes(col=Confidence_in_pos_effect)) +
          theme(legend.key.size = unit(1, 'cm'),
                legend.title = element_text(size=40),
                legend.text = element_text(size=40)) + scale_color_gradient2(midpoint=0.5, low="blue", mid="gray92",
                                                                             high="red2", limits = c(minimum, maximum)))
  dev.off()
  pdf(paste("figures/", ID_string, hab_string, "_tree_intercepts_", agg_string, ".pdf", sep = ""), width = 23, height = 20)
  print(ggtree(tree2, layout = "circular") + geom_tiplab(aes(geom = "point")) +
          geom_tippoint(pch=16, size=5, aes(col=exp(intercepts))) +
          theme(legend.key.size = unit(1, 'cm'),
                legend.title = element_text(size=40),
                legend.text = element_text(size=40)) + scale_color_gradient(low="gray92", high="red2"))
  dev.off()
  return(tree2)
}

# make a massive figure where the results are plotted out on a phylogeny
plot_trees_large <- function(ID_string, ID_string_class1, ID_string_class2, agg_string, tree, order_colours, meta, confs, filter = NULL) {
  # read in the necessary data from file
  affected_pos <- read.csv(paste(ID_string, "_positive_", agg_string, "_habitat_control.txt", sep = ""))$x
  affected_neg <- read.csv(paste(ID_string, "_negative_", agg_string, "_habitat_control.txt", sep = ""))$x
  print("Positive:")
  print(length(affected_pos))
  print(affected_pos)
  print("Negative:")
  print(length(affected_neg))
  affected_pos1 <- read.csv(paste(ID_string_class1, "_positive_", agg_string, "_habitat_control.txt", sep = ""))$x
  affected_neg1 <- read.csv(paste(ID_string_class1, "_negative_", agg_string, "_habitat_control.txt", sep = ""))$x
  affected_pos2 <- read.csv(paste(ID_string_class2, "_positive_", agg_string, "_habitat_control.txt", sep = ""))$x
  affected_neg2 <- read.csv(paste(ID_string_class2, "_negative_", agg_string, "_habitat_control.txt", sep = ""))$x
  lower <- read.csv(paste(ID_string, "_lower_", agg_string, ".txt", sep = ""))$x
  upper <- read.csv(paste(ID_string, "_upper_", agg_string, ".txt", sep = ""))$x
  class1_lower <- read.csv(paste(ID_string_class1, "_lower_", agg_string, "_habitat_control.txt", sep = ""))$x
  class1_upper <- read.csv(paste(ID_string_class1, "_upper_", agg_string, "_habitat_control.txt", sep = ""))$x
  class2_lower <- read.csv(paste(ID_string_class2, "_lower_", agg_string, "_habitat_control.txt", sep = ""))$x
  class2_upper <- read.csv(paste(ID_string_class2, "_upper_", agg_string, "_habitat_control.txt", sep = ""))$x
  species <- read.csv(paste(ID_string, "_species_", agg_string, ".txt", sep = ""))$x
  check1 <- read.csv(paste(ID_string, "_check1_", agg_string, ".txt", sep = ""))$x
  check2 <- read.csv(paste(ID_string, "_check2_", agg_string, ".txt", sep = ""))$x
  intercepts1 <- read.csv(paste(ID_string_class1, "_expected_500_", agg_string, "_habitat_control.txt", sep = ""))$x
  intercepts2 <- read.csv(paste(ID_string_class2, "_expected_500_", agg_string, "_habitat_control.txt", sep = ""))$x
  total <- check1 + check2
  affected_leaves <- get_affected_leaves(affected_pos, affected_neg, species)
  affected_leaves1 <- get_affected_leaves(affected_pos1, affected_neg1, species, binary = TRUE)
  affected_leaves2 <- get_affected_leaves(affected_pos2, affected_neg2, species, binary = TRUE)
  # get colours for the species based on their species group
  col_vector <- c()
  for (sp in species) {
    col_vector <- c(col_vector, order_colours[[meta[meta$species == sp, "group"]]])  
  }
  # associate numbers to species
  number_vector <- c()
  for (sp in species) {
    number_vector <- c(number_vector, meta[meta$species == sp, "numbers"])  
  }  
  # make a data frame with all the data that will be necessary for the plot
  curr_data <- data.frame("class1_lower" = class1_lower, "class1_upper" = class1_upper,
                          "class2_lower" = class2_lower, "class2_upper" = class2_upper,
                          "lower" = lower,
                          "upper" = upper,
                          "species" = species,
                          "Effect_of_time" = affected_leaves,
                          "Class1" = affected_leaves1,
                          "Class2" = affected_leaves2,
                          "total" = total,
                          "check1" = check1,
                          "check2" = check2,
                          "colours" = col_vector,
                          "intercepts1" = intercepts1,
                          "intercepts2" = intercepts2,
                          "Stronger_effect" = confs,
                          "numbers" = number_vector)
  # create file name
  file_name <- paste("figures/", ID_string, "_tree_large_", agg_string, "_habitat_control.pdf", sep = "")
  # set up various parameters of the plotting area
  height = 50
  n <- length(curr_data$species)
  HPDI_limits <- c(1 - 2.5, n + 2.6)
  bar_limits <- c(-2, n + 2.8)
  tree_factor <- 0.02
  bar_expand <- 1
  # if only certain species are to be plotted, filter the data to just those
  # and adjust the dimensions of the graphs accordingly
  if (!is.null(filter)) {
    curr_data <- curr_data[curr_data$species %in% filter,]
    file_name <- paste("figures/", ID_string, "_tree_small_", agg_string, "_habitat_control.pdf", sep = "")
    height = 50 * (length(filter)/length(species))
    n <- length(curr_data$species)
    HPDI_limits <- c(1 - 3.3, n + 4.1)
    tree_factor <- 0.04
    bar_expand = 1
    bar_limits <- c(-2, n + 3.7)
  }
  # join in the data from the files with the structure of the phylogeny
  x <- full_join(as_tibble(tree), curr_data, join_by(label == species))
  tree2 <- as.treedata(x)
  # create figure
  pdf(file_name, width = 27, height = height)
  # start by drawing the phylogeny at the very left
  # annotate negative divergence time effects with a black dot, positive effects with an orange dot
  # and the rest in grey
  tree_plot <- ggtree(tree2, layout = "rectangular") +
            geom_tippoint(pch=16, size=5, aes(col=Effect_of_time)) +
            scale_color_manual(values=c("black", "gray92","orange")) +
            theme(legend.key.size = unit(1, 'cm'),
                  legend.title = element_text(size=40),
                  legend.text = element_text(size=40),
        legend.position="bottom") +
    theme(legend.position = "none") +
    scale_y_continuous(expand=expand_scale(tree_factor))
  # loop over species groups
  for (order in names(order_colours)) {
    # get only the species in current group
    curr_species <- meta[meta$group == order, "species"]
    if (!is.null(filter)) {
      curr_species <- curr_species[curr_species %in% filter]
    }
    # draw a rectangle over the part of the tree corresponding to the
    # species group
    if (length(curr_species > 0)) {
      mrca <- MRCA(tree2, curr_species)
      tree_plot <- tree_plot + geom_hilight(node=mrca, fill=order_colours[[order]])
      tree_data <- fortify(tree2)
      curr <- tree_data[tree_data$node == mrca,]
      # add species group names
      size <- 10
      if (!is.null(filter) & order %in% c("Squamata", "Afrotheria", "Coelacanthi", "Testudines")) {
        size <- 6
      }
      tree_plot <- tree_plot + annotate("text", x = max(min(curr$x - 200, 500), 100), y = curr$y, label = order,
                                        size = unit(size, "pt"))
    }
  }
  # get the order of the tips of the tree so that the rest of the plot
  # could show the species in the same order
  tip_order <- get_taxa_name(tree_plot)
  rownames(curr_data) <- curr_data$species
  # reorder the data in the same order
  curr_data <- curr_data[tip_order,]
  curr_data$indices <- rev(seq(1, n))
  # add a number for each species
  tree_plot <- tree_plot + annotate("text", x = 700, y = curr_data$indices, label = curr_data$numbers,
                                    size = unit(2, "pt"))
  # plot out HPDIs
  # if filtered plot, do only one HDPI for everything
  # if all data, do separate ones for Class 1 and Class 2
  if (!is.null(filter)) {
    HPDI_plot <- ggplot(curr_data, aes(x=lower, y=indices, col = Effect_of_time)) + geom_point() +
      geom_point(aes(x=upper, y=indices, col = Effect_of_time)) + 
      labs(col = "Effect of divergence time")
    HPDI_plot <- HPDI_plot + geom_segment(aes(x = lower, y = indices, xend = upper, yend = indices, col = Effect_of_time),
                                          inherit.aes = FALSE, data = curr_data) +
      geom_vline(xintercept = 0) +
      scale_color_manual(values=c("black", "gray70","orange")) + guides(col = guide_legend(order = 2)) +
      geom_vline(xintercept = 0, linewidth = 1) + theme(axis.line.y = element_blank(),
                                                        axis.title.y = element_blank(),
                                                        axis.text.y = element_blank(),
                                                        axis.ticks.y = element_blank(), 
                                                        panel.grid.major = element_blank(), 
                                                        panel.grid.minor = element_blank(),
                                                        panel.background = element_blank(),
                                                        plot.title = element_text(size = 30),
                                                        axis.text.x = element_text(size = 20),
                                                        axis.title = element_text(size = 20)) +
      scale_y_continuous(limits = HPDI_limits, expand=c(0,tree_factor)) + 
        xlab("Coefficient of divergence time") +
      ggtitle("Effect of divergence time") + guides(col = guide_legend(order = 1)) +
      theme(legend.position = "none")
  }
  else {
    HPDI_plot <- ggplot(curr_data, aes(x=class2_lower, y=indices, col = Class2)) + geom_point() +
      geom_point(aes(x=class2_upper, y=indices, col = Class2))
    HPDI_plot <- HPDI_plot + geom_segment(aes(x = class2_lower, y = indices, xend = class2_upper, yend = indices, col = Class2),
                                          inherit.aes = FALSE, data = curr_data) +
      geom_vline(xintercept = 0) +
      scale_color_manual(values=c("red", "pink")) + guides(col = guide_legend(order = 2))
    HPDI_plot <- HPDI_plot + new_scale_color() + geom_point(aes(x=class1_lower, y=indices - 0.2, col = Class1)) +
      geom_point(aes(x=class1_upper, y=indices - 0.2, col = Class1))
    HPDI_plot <- HPDI_plot + geom_segment(aes(x = class1_lower, y = indices - 0.2, xend = class1_upper, yend = indices - 0.2, col = Class1),
                                          inherit.aes = FALSE, data = curr_data) +
      geom_vline(xintercept = 0, linewidth = 1) + theme(axis.line.y = element_blank(),
                                                        axis.title.y = element_blank(),
                                                        axis.text.y = element_blank(),
                                                        axis.ticks.y = element_blank(), 
                                                        panel.grid.major = element_blank(), 
                                                        panel.grid.minor = element_blank(),
                                                        panel.background = element_blank(),
                                                        plot.title = element_text(size = 30),
                                                        axis.text.x = element_text(size = 20),
                                                        axis.title = element_text(size = 20)) +
      scale_y_continuous(limits = HPDI_limits, expand=c(0,tree_factor)) +
      scale_color_manual(values=c("blue", "lightblue")) + xlab("Coefficient of divergence time") +
      ggtitle("Effect of divergence time") + guides(col = guide_legend(order = 1))
    # annotate based on whether the subsampling analysis gave a stronger effect for Class 1
    # or Class 2
    HPDI_plot <- HPDI_plot + new_scale_color() + geom_point(aes(x = 5, y = indices - 0.1, col = Stronger_effect), shape = 19) +
      scale_color_manual(values=c("white", "lightblue", "blue", "pink", "red")) + guides(col = guide_legend(order = 3)) +
      theme(legend.position = "none")
  }
  colours <- curr_data$colours
  curr_data$species <- factor(curr_data$species, levels = rev(tip_order))
  curr_data$indices <- rev(seq(0, n + 0.5, length.out = n))
  x_breaks <- seq(-2, +3, by = 1)
  x_labels <- 10 ^ abs(x_breaks)
  x_labels[x_labels == 1] <- 0
  # make a bar plot with the observed number of HTT events per species,
  # with Class 1 and Class 2 separately
  raw_bar_plot <- ggplot(data = curr_data, mapping = aes(y = log10(check2 + 1), x = indices, 
                                                         fill = species)) + 
    geom_bar(stat="identity", width = 1) +
    coord_flip(ylim = c(-max(log10(curr_data$check1)), max(log10(curr_data$check2 + 1)))) +
    theme(axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 30),
          axis.text.x = element_text(size = 20),
          axis.title = element_text(size = 20)) +
    scale_fill_manual(values = rev(colours)) +
    scale_x_continuous(limits = bar_limits, expand=c(0,bar_expand)) +
    geom_bar(mapping = aes(x = indices, y = -log10(check1 + 1)), stat="identity", width = 1) +
    geom_hline(yintercept = 0, linewidth = 1) +
    scale_y_continuous("", breaks = x_breaks, labels = x_labels) +
    ggtitle("Total # of HTT")
  raw_bar_plot <- raw_bar_plot + annotate("text", x = -2, y = c(-1, +1), label = c("Class 1", "Class 2"),
                                    size = unit(10, "pt"))
  # raw_bar_plot <- raw_bar_plot + geom_vline(xintercept = curr_data$indices, linewidth = 0.5)
  x_breaks <- seq(-2, +3, by = 1)
  x_labels <- 10 ^ abs(x_breaks)
  x_labels[x_labels == 1] <- 0
  # another bar plot but this time with the model estimates for a fixed divergence time
  est_bar_plot <- ggplot(data = curr_data, mapping = aes(y = log10(intercepts2 + 1), x = indices, 
                                                         fill = species)) + 
    geom_bar(stat="identity", width = 1) +
    coord_flip() +
    theme(axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 30),
          axis.text.x = element_text(size = 20),
          axis.title = element_text(size = 20)) +
    scale_fill_manual(values = rev(colours)) +
    scale_x_continuous(limits = bar_limits, expand=c(0,bar_expand)) +
    geom_bar(mapping = aes(x = indices, y = -log10(intercepts1+1)), stat="identity", width = 1) +
    geom_hline(yintercept = 0, linewidth = 1) +
    ggtitle("Est. for 500 MY divergence") +
    scale_y_continuous("", breaks = x_breaks, labels = x_labels) +
    xlab("")
  est_bar_plot <- est_bar_plot + annotate("text", x = -2, y = c(-1, +1), label = c("Class 1", "Class 2"),
                                          size = unit(10, "pt"))
  # put the graphs side-by-side
  grid.arrange(tree_plot, HPDI_plot, raw_bar_plot, est_bar_plot, ncol=4)
  dev.off()
  return(tree2)
}

# subsampling analysis to compare the strength of the effect between 
# Class 1 and Class 2, with Class 2 events subsampled to be of the same number as Class 1 events
# the analysis is performed for one specific species
subsampling_analysis <- function(HTT, species, iter = 10, visual = FALSE) {
  # make a table with just the data for the species of interest
  curr <- HTT[HTT$species.1 == species,]
  columns <- c("divTime", "spClade.2", "n")
  curr1 <- aggregate(curr[, "type1"] ~ curr[,"divTime"] + curr[,"spClade.2"], FUN = median)
  colnames(curr1) <- columns
  curr1$n <- round(curr1$n)
  curr2 <- aggregate(curr[, "type2"] ~ curr[,"divTime"] + curr[,"spClade.2"], FUN = median)
  colnames(curr2) <- columns
  curr2$n <- round(curr2$n)
  # don't continue if there are fewer than 3 events for either Class 1 or Class 2
  if (sum(curr2$n) <= 3 | sum(curr1$n) <= 3) {
    return(NA)
  }
  # the probability of sampling at different divergence times for Class 2
  probs <- curr2$n/sum(curr2$n)
  # correlation coefficient between divergence time and # HTT for Class 1
  corr1 <- cor.test(log(curr1$n + 1), curr1$divTime, method = "spearman")$estimate
  # the full correlation test output for both Class 1 and Class 2, 
  # just for printing out
  corr1_all <- cor.test(log(curr1$n + 1), curr1$divTime, method = "spearman")
  corr2_all <- cor.test(log(curr2$n + 1), curr2$divTime, method = "spearman")
  # if visual is set to TRUE, make a scatter plot of the relationship with divergence
  # time, for Class 1 and Class 2 separately
  if (visual == TRUE) {
    par(mfrow = c(1,2))
    plot(log(curr1$n + 1) ~ curr1$divTime, main = "Class I", col = "blue")
    plot(log(curr2$n + 1) ~ curr2$divTime, col = "red", main = "Class II")
  }
  # now we get to the subsampling
  sim_values <- c()
  # perform iter iterations
  for (sim in 1:iter) {
    temp2 <- curr2
    # sample from Class 2, using the probabilities of real Class 2 data for each divergence
    # time but the sample size of Class 1 data
    sampled <- sample(curr2$spClade.2, size = sum(curr1$n), prob = probs, replace = TRUE)
    # make a new data frame with the subsampled data
    new_n <- c()
    for (spc2 in curr2$spClade.2) {
      new_n <- c(new_n, sum(sampled == spc2))
    }
    temp2$n <- new_n
    # correlation coefficient with subsampled data
    corr2 <- cor.test(log(temp2$n + 1), temp2$divTime, method = "spearman")$estimate
    # store corrleation coefficient
    sim_values <- c(sim_values, corr2)
  }
  # calculate proportion of simulations where the correlation coefficient 
  # for Class 1 is more negative than the coefficient obtained with the
  # subsampling (by at least 0.05)
  conf_class1 <- mean((sim_values - corr1) >= 0.05)
  # same thing for Class 2
  conf_class2 <- mean((corr1 - sim_values) >= 0.05)
  # print stuff out
  print(species)
  print("class 1:")
  print(corr1_all)
  print("class 2:")
  print(corr2_all)
  print("sim:")
  print(sim_values[1:50])
  print(conf_class1)
  print(conf_class2)
  return(c(conf_class1, conf_class2))
}

######

# Define colours
######
desert <- rgb(1, 0.89, 0.77, alpha = 0.5)
desert2 <- rgb(1, 0.89, 0.77)
red_trans <- rgb(0.38, 0.16, 0.06, alpha = 0.5)
######

# Read in and preprocess data (this subsection only needs to be run once, afterwards one can read the clean
# data file below)
######
meta <- read.csv("metadata_210623_dbBusco_orderish.tbl", sep = "\t")
HTT <- get_HTT("nbHTTevents_perPairAll_noNeg_greedyCorrected_withClades.txt")
# class 1 and class 2 TEs
HTT_type1 <- read.csv("nbHTTevents_perPairAll_noNeg_greedyCorrected_LINEandLTR.txt", sep = "\t")
HTT_type2 <- read.csv("nbHTTevents_perPairAll_noNeg_greedyCorrected_DNAandRC.txt", sep = "\t")
# add TE class information to main data frame
type1_counts <- c()
type2_counts <- c()
# HTT will have two copies of every HTT event with species 1 and species 2 switched,
# whereas the class1/class2 files only have one
for (line in 1:dim(HTT)[1]) {
  match <- HTT_type1[HTT_type1$species.1 == HTT[line,"species.1"] & HTT_type1$species.2 == HTT[line,"species.2"],]
  if (dim(match)[1] == 0) {
    match <- HTT_type1[HTT_type1$species.2 == HTT[line,"species.1"] & HTT_type1$species.1 == HTT[line,"species.2"],]
  }
  type1_counts <- c(type1_counts, match$n)
  match <- HTT_type2[HTT_type2$species.1 == HTT[line,"species.1"] & HTT_type2$species.2 == HTT[line,"species.2"],]
  if (dim(match)[1] == 0) {
    match <- HTT_type2[HTT_type2$species.2 == HTT[line,"species.1"] & HTT_type2$species.1 == HTT[line,"species.2"],]
  }
  type2_counts <- c(type2_counts, match$n)
}
HTT$type1 <- type1_counts
HTT$type2 <- type2_counts
# store version with class1/class2 information added to file
write.table(HTT, file = "HTT.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
HTT <- read.csv("HTT.txt", sep = "\t")
# check that class 1 + class 2 == n
summary(rowSums(cbind(HTT$type1, HTT$type2)) == HTT$n)
######
# read in clean data file
HTT <- read.csv("HTT.txt", sep = "\t")
# metadata
meta <- read.csv("metadata_210623_dbBusco_orderish.tbl", sep = "\t")
# phylogeny
nwk <- ("species_247leaves_noRateau.nwk")
tree <- read.newick(nwk)
# colours for each species group
order_cols <- read.csv("colorsTaxa.txt", sep = "\t")
order_colours <- as.list(order_cols$col)
names(order_colours) <- order_cols[,1]
# modify a couple of the colours to be more visible
order_colours[["Ephydroidea"]] <- "wheat2"
order_colours[["Amphibia"]] <- "lightgoldenrod"
# species numbers
sp_numbers <- read.csv("tableSA_colors.csv", sep = ";")
numbers <- c()
for (sp in meta$species) {
  numbers <- c(numbers, sp_numbers[sp_numbers$speciesName == sp, "speciesNumber"])
}
meta$numbers <- numbers
# the names of all the species groups
orders <- unique(meta$group)

# Modelling the relationship between divTime and n
######
# fit the negative binomial models
# take the median of each spClade for the partner species
model_and_write("agg", orders, meta, HTT, agg = TRUE, agg_function = median, agg_string = "median", iter = 10000)
# Class 1 only
model_and_write("class1", orders, meta, HTT, agg = TRUE, agg_function = median, agg_string = "median", alternative = "type1", iter = 10000)
# Class 2 only
model_and_write("class2", orders, meta, HTT, agg = TRUE, agg_function = median, agg_string = "median", alternative = "type2", iter = 10000)
# test the effect of shared habitat
model_and_write("agg", orders, meta, HTT, agg = TRUE, agg_function = median, agg_string = "median", iter = 10000, habitat = TRUE)
# the effect of divTime controlling for shared habitat
# this is the main analysis reported in the manuscript
model_and_write("agg", orders, meta, HTT, agg = TRUE, agg_function = median, agg_string = "median", iter = 10000, habitat = TRUE, habitat_control = TRUE)
model_and_write("class1", orders, meta, HTT, agg = TRUE, agg_function = median, agg_string = "median", alternative = "type1", iter = 10000, habitat = TRUE, habitat_control = TRUE)
model_and_write("class2", orders, meta, HTT, agg = TRUE, agg_function = median, agg_string = "median", alternative = "type2", iter = 10000, habitat = TRUE, habitat_control = TRUE)

# plot Class 1 and Class 2 together
fit_model_classes(HTT, "Dendrolimus_kikuchii")
fit_model_classes(HTT, "Dendrolimus_kikuchii", subsample = TRUE)

# plots of single species
fit_model_single(HTT, "Dendrolimus_kikuchii")
fit_model_single(HTT, "Amia_calva")
fit_model_single(HTT, "Scaeva_pyrastri")
fit_model_single(HTT, "Macaca_fascicularis")
fit_model_single(HTT, "Erpetoichthys_calabaricus")
fit_model_single(HTT, "Metaphire_vulgaris")
fit_model_single(HTT, "Hymenochirus_boettgeri")
fit_model_single(HTT, "Batillaria_attramentaria")
fit_model_single(HTT, "Portunus_trituberculatus")
fit_model_single(HTT, "Argiope_bruennichi")
fit_model_single(HTT, "Nymphon_striatum")
fit_model_single(HTT, "Arion_vulgaris")
fit_model_single(HTT, "Trichechus_manatus_latirostris")

# Would a non-linear model be more appropriate?
######
# pick one species group as an example
order <- "Neuropterida"
# filter data to only this species
to_check <- meta[meta$group == order, "species"]
to_check <- to_check[to_check %in% HTT$species.1]
# mean and SD of divergence time to be able to convert to and from
# Z-scores
time_mean <- mean(HTT$divTime)
time_sd <- sd(HTT$divTime)
# store the different possibel models here
models <- list()
# loop over the species
for (sp in to_check) {
  # linear negative binomial regression with log link
  jpeg(paste("figures/", sp, "_nonlinearity.jpeg", sep = ""), width = 600, height = 300)
  curr <- HTT[HTT$species.1 == sp,]
  curr_temp <- aggregate(n ~ zdivTime + spClade.2, data = curr, FUN = median)
  curr_temp$n <- round(curr_temp$n)
  curr_temp2 <- aggregate(n ~ divTime + spClade.2, data = curr, FUN = median)
  curr <- cbind(curr_temp, curr_temp2$divTime)
  colnames(curr) <- c("zdivTime", "spClade.2", "n", "divTime")
  prior_b <- set_prior("normal(0,2)", class = "b")
  prior_i <- set_prior("normal(0,10)", class = "Intercept")
  prior_shape <- set_prior("gamma(1, 1)", class = "shape",
                           lb = 0)
  model1 <- brm(n ~ zdivTime, data = curr, family = negbinomial(),
               iter = 2000, control = list(adapt_delta = 0.99),
               backend = "cmdstanr", cores = 2,
               prior = c(prior_b, prior_i, prior_shape), silent = 2, refresh = 0)
  models[[sp]][["model1"]] <- model1
  # plot model estimates, as well as real data
  par(mfrow = c(1,2))
  ppd <- posterior_epred(model1)
  real <- aggregate(n ~ divTime + spClade.2, data = curr,
                    FUN = median)
  plot(colMeans(ppd) ~ curr$divTime, col = desert, main = paste(sp, "\nlinear", sep = ""),
       ylim = c(0,max(real[,3]) + 1), xlab = "divTime", ylab = "# HTT events",
       pch = "")
  for (example in 1:1000) {
    lines(ppd[example,] ~ curr$divTime, col = desert, type = "l")
  }
  points(real[,3] ~ real[,1], pch = 20)
  # negative binomial GAM
  model2 <- brm(n ~ s(zdivTime, k = 8), data = curr, family = negbinomial(),
                iter = 2000, control = list(adapt_delta = 0.99),
                backend = "cmdstanr", cores = 2,
                prior = c(prior_b, prior_i), silent = 2, refresh = 0)
  models[[sp]][["model2"]] <- model2
  ppd <- posterior_epred(model2)
  plot(colMeans(ppd) ~ curr$divTime, col = desert2, main = paste(sp, "\nnon-linear", sep = ""), type = "l", lwd = 2,
       ylim = c(0,max(real[,3]) + 1))
  for (example in 1:1000) {
    lines(ppd[example,] ~ curr$divTime, col = desert, type = "l")
  }
  points(real[,3] ~ real[,1], pch = 20)
  dev.off()
  print(sp)
  print(loo(model1, model2))
  # make plots with the posterior predictions from the different models
  jpeg(paste("figures/", sp, "_PPD_sim.jpeg", sep = ""), width = 900, height = 600)
  par(mfrow = c(2, 3))
  stripchart(real[,3] ~ real[,1], pch = 20, main = paste(sp, "\nobserved", sep = ""),
       ylim = c(0,max(real[,3]) + 5), xlab = "divTime", ylab = "# HTT events",
       method = "jitter", vertical = TRUE)
  ppd <- posterior_predict(models[[sp]][["model1"]])
  simulated <- sample(seq(1,dim(ppd)[1]), 5)
  for (sim in simulated) {
    stripchart(ppd[sim,] ~ real[,1], pch = 20, main = paste(sp, "\nsim ", sim, sep = ""),
         ylim = c(0,max(real[,3]) + 5), xlab = "divTime", ylab = "# HTT events",
         col = "coral1", method = "jitter", vertical = TRUE)
    
  }
  dev.off()
  jpeg(paste("figures/", sp, "_PPD_sim_nonlinear.jpeg", sep = ""), width = 900, height = 600)
  par(mfrow = c(2, 3))
  stripchart(real[,3] ~ real[,1], pch = 20, main = paste(sp, "\nobserved", sep = ""),
       ylim = c(0,max(real[,3]) + 5), xlab = "divTime", ylab = "# HTT events",
       method = "jitter", vertical = TRUE)
  ppd <- posterior_predict(models[[sp]][["model2"]])
  simulated <- sample(seq(1,dim(ppd)[1]), 5)
  for (sim in simulated) {
    stripchart(ppd[sim,] ~ real[,1], pch = 20, main = paste(sp, "\nsim ", sim, sep = ""),
         ylim = c(0,max(real[,3]) + 5), xlab = "divTime", ylab = "# HTT events",
         col = "steelblue2", method = "jitter", vertical = TRUE)
    
  }
  dev.off()
}
######

# Read in modelling results from various files
######
check1 <- read.csv("agg_check1_median.txt", sep = "\t")$x
check2 <- read.csv("agg_check2_median.txt", sep = "\t")$x
species <- read.csv("agg_species_median.txt", sep = "\t")$x
b01 <- read.csv("class1_intercepts_median.txt", sep = "\t")$x
b02 <- read.csv("class2_intercepts_median.txt", sep = "\t")$x
b0 <- read.csv("agg_intercepts_median.txt", sep = "\t")$x
expected_500 <- read.csv("agg_expected_500_median_habitat_control.txt", sep = "\t")$x
expected_1300 <- read.csv("agg_expected_1300_median_habitat_control.txt", sep = "\t")$x
check1_500 <- read.csv("class1_expected_500_median.txt", sep = "\t")$x
check1_1300 <- read.csv("class1_expected_1300_median.txt", sep = "\t")$x
check2_500 <- read.csv("class2_expected_500_median.txt", sep = "\t")$x
check2_1300 <- read.csv("class2_expected_1300_median.txt", sep = "\t")$x
coef1 <- read.csv("class1_coefficients_median.txt", sep = "\t")$x
coef2 <- read.csv("class2_coefficients_median.txt", sep = "\t")$x
coef_combined <- read.csv("agg_coefficients_median.txt", sep = "\t")$x
coef_habitat <- read.csv("agg_coefficients_median_habitat.txt", sep = "\t")$x
lower_shared_habitat <- read.csv("agg_lower_median_habitat.txt", sep = "\t")$x
upper_shared_habitat <- read.csv("agg_upper_median_habitat.txt", sep = "\t")$x
coef_combined_habitat_control <- read.csv("agg_coefficients_median_habitat_control.txt", sep = "\t")$x
habitat_same <- read.csv("agg_habitat_same_median_habitat.txt", sep = "\t")$x
habitat_diff <- read.csv("agg_habitat_diff_median_habitat.txt", sep = "\t")$x
habitat_confidence <- read.csv("agg_confidence_median_habitat.txt", sep = "\t")$x
lower_no_habitat_control <- read.csv("agg_lower_median.txt", sep = "\t")$x
upper_no_habitat_control <- read.csv("agg_upper_median.txt", sep = "\t")$x
lower_habitat_control <- read.csv("agg_lower_median_habitat_control.txt", sep = "\t")$x
upper_habitat_control <- read.csv("agg_upper_median_habitat_control.txt", sep = "\t")$x
######

# effect of shared habitat
######
# filter only to species that could be modelled
to_keep <- complete.cases(lower_shared_habitat)
species_shared_habitat <- species[to_keep]
lower_shared_habitat_complete <- lower_shared_habitat[to_keep]
upper_shared_habitat_complete <- upper_shared_habitat[to_keep]
coef_habitat_complete <- coef_habitat[to_keep]
# effect of shared habitat
mean(coef_habitat_complete > 0)
mean(lower_shared_habitat_complete > 0)
mean(upper_shared_habitat_complete < 0)
sum(lower_shared_habitat_complete > 0)
sum(upper_shared_habitat_complete < 0)
# for the plot, order the numbers the same as the species
ordered_numbers <- c()
for (sp in species) {
  ordered_numbers <- c(ordered_numbers, meta[meta$species == sp, "numbers"])
}
plot_trees("agg", "median", tree, habitat = TRUE, numbers = ordered_numbers)
plot_trees("agg", "median", tree, habitat_control = TRUE, numbers = ordered_numbers)

######
pdf(file = "figures/shared_habitat_effect.pdf", width = 10, height = 20)
yaxis <- seq(1, length(habitat_same))
maximum <- max(log10(c(habitat_diff, habitat_same) + 1))
par(xpd = TRUE)
plot(yaxis ~ log10(habitat_same + 1), pch = "",
     xlim = c(0, maximum), ylab = "",
     xlab = "Log10 of HTT at 500 MY",
     yaxt = "n", bty = "n", xaxt = "n")
labels <- 10 ^ seq(0, 2)
labels[1] <- 0
axis(1, at = seq(0, 2), labels = labels)
starting_point <- 1
for (order in orders) {
  curr_species <- meta$species[meta$group == order]
  curr_indices <- species_shared_habitat %in% curr_species
  species_number <- sum(curr_indices)
  curr_yaxis <- seq(starting_point, starting_point + species_number - 1)
  segments(log10(habitat_diff[curr_indices] + 1), curr_yaxis, log10(habitat_same[curr_indices] + 1), curr_yaxis)
  points(curr_yaxis ~ log10(habitat_same[curr_indices] + 1), pch = 19)
  points(curr_yaxis ~ log10(habitat_diff[curr_indices] + 1), pch = 1)  
  rect(0, starting_point - 0.5, maximum, starting_point + species_number - 1 + 0.5,
       col = alpha(order_cols[order_cols$group == order, "col"], alpha = 0.5),
       border = NA)
  text(-0.1, y = mean(c(starting_point, (starting_point + species_number - 1))),
       col = "black",
       labels = order)
  starting_point <- starting_point + species_number
}
legend("topright", legend = c("shared habitat", "different habitat"),
       pch = c(19, 1))
dev.off()
######
######

# calculate fold change
fc <- expected_500/expected_1300
species[fc == median(fc, na.rm = TRUE)]

# how many species have 95% HPDI of divTime entirely below 0?
sum(complete.cases(upper_habitat_control))
sum(upper_habitat_control[complete.cases(upper_habitat_control)] < 0)
# how many species have 95% HPDI of divTime entirely above 0?
sum(lower_habitat_control[complete.cases(lower_habitat_control)] > 0, na.rm = TRUE)


# how many times could the model be fit for Class 1/Class 2/everything?
sum(!is.na(coef1))
sum(!is.na(coef2))
sum(!is.na(coef_combined))

# check if flies and fish show more of a phylogenetic proximity effect 
######
fly_species <- meta[meta$group %in% c("Syrphoidea", "Nematocera", "Ephydroidea"), "species"]
fly_coefs <- coef_combined[species %in% fly_species]
nonfly_coefs <- coef_combined[!species %in% fly_species]
vioplot(fly_coefs, nonfly_coefs)
wilcox.test(fly_coefs, nonfly_coefs)
median(fly_coefs, na.rm = TRUE)
median(nonfly_coefs, na.rm = TRUE)

fish_species <- meta[meta$group == "Actinopterygii", "species"]
fish_species <- meta[meta$group == "Amphibia", "species"]
fish_coefs <- coef_combined[species %in% fish_species]
nonfish_coefs <- coef_combined[!species %in% fish_species]
vioplot(fish_coefs, nonfish_coefs)
wilcox.test(fish_coefs, nonfish_coefs)
median(fish_coefs, na.rm = TRUE)
median(nonfish_coefs, na.rm = TRUE)
######

# Compare the # of transfer events between different clades
######
# barplot of the observed number of transfers vs the number of transfers expected at 
# the mean divergence time, based on the modelling results
pdf("figures/n_vs_b0_barplot.pdf", width = 11, height = 7)
par(mfrow = c(1,2))
barplot_comp(meta, orders, order_colours, (check1 + check2), "All HTT events", 
             "# HTT events", seed, species)
barplot_comp(meta, orders, order_colours, exp(b0), "Expected for mean div. time", 
             "# HTT events", seed, species)
dev.off()
# barplot of the observed number of transfers, + expectation at 500 MYA and 1300 MYA
pdf("figures/n_vs_b0_barplot_three.pdf", width = 17, height = 7)
par(mfrow = c(1,3))
barplot_comp(meta, orders, order_colours, (check1 + check2), "All HTT events", 
                               "# HTT events", seed, species)
barplot_comp(meta, orders, order_colours, expected_500, "Expected for 500 MYA", 
             "# HTT events", seed, species)
barplot_comp(meta, orders, order_colours, expected_1300, "Expected for 1300 MYA", 
             "# HTT events", seed, species)
dev.off()

######

# Check for habitat effect
######
# a vector that will contain the habitat of each species
habitats <- c()
# a vector that will contain the phylogenetic group
# of each species
groups <- c()
# a vector that will contain the total number of HTT for each species
n <- c()
# construct the two vectors
for (sp in species) {
  habitats <- c(habitats, HTT[HTT$species.1 == sp, "habitat.1"][1])
  n <- c(n, sum(HTT[HTT$species.1 == sp, "n"]))
  curr_group <- meta[meta$species == sp, "group"]
  groups <- c(groups, curr_group)
}
vioplot(coef_combined ~ habitats)

# Habitat violin plot
######
# violin of the difference in HTT number
# between aquatic and terrestrial
# per species group
######
pdf("figures/effect_of_habitat_vioplot.pdf", width = 8, height = 5.5)
groups_keep <- c()
diffs_keep <- c()
aquatic <- c()
terrestrial <- c()
groups_long_aquatic <- c()
groups_long_terrestrial <- c()
# only consider species group where there are both aquatic and 
# terrestrial species
for (group in unique(meta$group)) {
  curr_species <- meta$species[meta$group == group]
  terra <- expected_500[species %in% curr_species & habitats == "terrestrial"]
  terra <- terra[is.finite(terra)]
  if (length(terra) > 0) {
    mare <- expected_500[species %in% curr_species & habitats == "aquatic"]
    mare <- mare[is.finite(mare)]
    if (length(mare) > 0) {
      aquatic <- c(aquatic, mare)
      terrestrial <- c(terrestrial, terra)
      groups_long_aquatic <- c(groups_long_aquatic, rep(group, length(mare)))
      groups_long_terrestrial <- c(groups_long_terrestrial, rep(group, length(terra)))
      diffs_keep <- c(diffs_keep, log(median(mare) + 1) - log(median(terra) + 1))
      groups_keep <- c(groups_keep, group)
    }
  }
}
diff_order <- groups_keep[order(diffs_keep, decreasing = TRUE)]
groups_long_aquatic <- factor(groups_long_aquatic, levels = diff_order)
groups_long_terrestrial <- factor(groups_long_terrestrial, levels = diff_order)
par(mfrow = c(1,1))
par(xpd = TRUE)
x <- seq(1, length(groups_keep) * 2, by = 2)
# cols <- rep("steelblue", length(diffs_keep))
# cols[diffs_keep < 0] <- desert2
stripchart(log2(aquatic+1) ~ groups_long_aquatic,
           col = "skyblue3", at = x - 0.9,
           vertical = TRUE, pch = "", 
           xlab = "", las = 2,
           ylab = "Est. # HTT at 500 MY divergence",
           yaxt = "n", xlim = c(0, max(x)), ylim = c(0,6),
           xaxt = "n")
rect_starts <- seq(-0.5, x[length(x)], by = 4)
rect_ends <- seq(1.5, x[length(x)] + 0.5, by = 4)
for (slot in 1:length(rect_starts)) {
  rect(rect_starts[slot], 0, rect_ends[slot], 6, col = "burlywood1",
       border = NA)  
}
axis(1, at = x-0.5, labels = rep("", length(x)))
text(x = x-0.5, y = -1, labels = levels(groups_long_aquatic),
     cex = 0.55, srt = 50)
axis(2, at = seq(0,6), labels = 2 ^ (seq(0,6)) - 1, las = 2)
vioplot(log2(aquatic+1) ~ groups_long_aquatic,
           col = "white", at = x - 0.9, add = TRUE)
stripchart(log2(aquatic+1) ~ groups_long_aquatic,
           col = "skyblue3", at = x - 0.9,
           vertical = TRUE, pch = 20, add = TRUE)
vioplot(log2(terrestrial+1) ~ groups_long_terrestrial,
        col = "white", at = x, add = TRUE)
stripchart(log2(terrestrial+1) ~ groups_long_terrestrial,
        col = "burlywood3", add = TRUE, at = x,
        vertical = TRUE, pch = 20)
# barplot(diffs_keep ~ factor(groups_keep, levels = groups_keep), ylab = "Log-fold excess in aquatic (500 MY divergence)", xlab = "",
#         las = 2, cex.names = 0.5, col = cols,
#         at = x + 0.5)
dev.off()
# test to compare between habitats
for_test <- aggregate(expected_500 ~ habitats + groups,
                      FUN = median)
for_test <- for_test[for_test$groups %in% groups_keep,]
t.test(log(for_test$expected_500)[for_test$habitats == "aquatic"],
       log(for_test$expected_500)[for_test$habitats == "terrestrial"],
       paired = TRUE)
######

# Class 1 vs Class 2 subsampling
######
confs1 <- c()
confs2 <- c()
counter <- 0
# loop over the species
for (sp in species) {
  counter <- counter + 1
  if (counter %% 10 == 0) {
    print(counter)
  }
  # apply the subsampling analysis to each species
  curr_confs <- subsampling_analysis(HTT, sp, iter = 10000)
  # confidence that Class 1 shows stronger effect
  confs1 <- c(confs1, curr_confs[1])
  # confidence that Class 2 shows stronger effect
  confs2 <- c(confs2, curr_confs[2])
}

# write the results into file as this is slow
conf_table <- data.frame("species" = species, "confs_class1" = confs1, "confs_class2" = confs2)
write.table(conf_table, file = "class1_class2_confs.csv", sep = ",", quote = FALSE, col.names = TRUE,
            row.names = FALSE)
# read from file
confs_all <- read.csv("class1_class2_confs.csv", sep = ",")
confs1 <- confs_all[,2]
confs2 <- confs_all[,3]
# turn into categorical variable
confs_categorical <- rep("no evidence", length(confs1))
confs_categorical[(confs1 > 0.8 & !is.na(coef1)) & !is.na(confs1)] <- "weak evidence for Class I"
confs_categorical[(confs1 > 0.95 & !is.na(coef1)) & !is.na(confs1)] <- "strong evidence for Class I"
confs_categorical[(confs2 > 0.8 & !is.na(coef2)) & !is.na(confs2)] <- "weak evidence for Class II"
confs_categorical[(confs2 > 0.95 & !is.na(coef2)) & !is.na(confs2)] <- "strong evidence for Class II"
confs_categorical <- factor(confs_categorical, levels = c("no evidence", "weak evidence for Class I",
                                                          "strong evidence for Class I", 
                                                          "weak evidence for Class II",
                                                          "strong evidence for Class II"))
# count the cases
sum(complete.cases(confs_all))
sum(confs_categorical == "strong evidence for Class I")
sum(confs_categorical == "weak evidence for Class I")
sum(confs_categorical == "strong evidence for Class II")
sum(confs_categorical == "weak evidence for Class II")
######

# Make massive figure
######
# all species
tree2 <- plot_trees_large("agg", "class1", "class2", "median", tree, order_colours, meta, confs_categorical)
# smaller subset
small_tree <- read.newick("species_153leaves_repSpecies.nwk")
species_sampled <- read.csv("spRep_clade.txt", sep = "\t")[,1]
tree3 <- plot_trees_large("agg", "class1", "class2", "median", small_tree, order_colours, meta, confs_categorical, filter = species_sampled)

######