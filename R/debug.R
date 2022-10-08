



# Loop over effect size
for (eff_size in seq(

  from = l$params$min_effect_size,
  to = l$params$max_effect_size,
  by = l$params$effect_size_step)

) {


  # Pre-allocate memory to results
  m_pvals <- matrix(

    nrow = l$params$max_vil_sampled,
    # row number is equal to number of villager per treatment group

    ncol = l$params$max_hh_sample_per_vil)
  # column number is equal to number of HHs per village


  # Matrix must not be NA or operations will only yield NA
  m_pvals[is.na(m_pvals)] <- 0


  # Simulate results many times
  for (iter in 1:l$params$nu_simulations) {

    # Set a fixed seed for reproducibility
    set.seed(1944 + iter) # Year of Bretton woods


    # Template of village-household universe
    vil_hh_u <- data.table()


    # Creating the universe - loop over villages
    for (village in 1:l$params$vil_in_universe) {


      # Select a village size
      vil_size <- sample(l$params$min_vil_size:l$params$max_vil_size, size = 1)


      # Select if village will be treated or not
      treatment <- sample(0:1,
                          size = 1,
                          prob = c(1 - l$params$prob_vil_is_treated,
                                   l$params$prob_vil_is_treated)
      )


      # Create a village data set
      tmp_vil <- data.table(

        # Village ID
        "village" = village,

        # Within-village Household ID
        "household" = 1:vil_size,

        # Create baseline score
        "baseline" = rnorm(vil_size,

                           # Mean depends on the village
                           mean = runif(1,
                                        min = l$params$bl_min_mean_val,
                                        max = l$params$bl_max_mean_val),

                           # SD depends on village
                           sd = runif(1,
                                      min = l$params$bl_min_sd_val,
                                      max = l$params$bl_max_sd_val)),


        # Treatment status - village-wide effect
        "treatead" = treatment %>% rep(vil_size))


      # Concatenate village data set with all previous villages
      vil_hh_u <- rbind(vil_hh_u, tmp_vil)
    }


    # # Number of villages in treated and control groups
    # vil_hh_u[, .("Nu. of villages" = uniqueN(.SD, by = "village")), by = .(treatead)]
    # # commented because adds little value inside the loop


    # Homogeneous effects of treatment
    vil_hh_u[, outcome := baseline + rnorm(.N)]
    vil_hh_u[treatead == 1, outcome := baseline + rnorm(.N, mean = eff_size)]
    # Treatment magnitude is the mean value


    # New seed for stage-specific reproducibility
    set.seed(2005 + iter) # Year DIME was created


    ## Parallelize

    # Runs for different village sample sizes the effect of varying
    # household sample size
    p_result <- clusterApply(clusters,
                             x = 1:l$params$max_vil_sampled,
                             fun = change_hh_size,
                             universe = vil_hh_u)


    # Unlist results into a matrix
    tmp_mat <- p_result %>% unlist %>%
      matrix(ncol = l$params$max_hh_sample_per_vil, byrow = T)
    # columns represent number HHs sampled per village,
    # rows represent the number of villages per treatment group


    # Warn if there are NAs, then replace NA's with 0
    if (any(is.na(tmp_mat))) {
      warning(
        paste0("NAs detected. effect size: ", eff_size,
               " Iteration: ", iter)
      )
    }
    tmp_mat[is.na(tmp_mat)] <- 0


    # Add results from previous iterations
    m_pvals <- m_pvals + tmp_mat


    # Verbalize iteration and time
    cat(
      paste0(
        "\n\nFinished cycle: ", iter, " - for effect size: ", eff_size,
        "\n", Sys.time())
    )

    flush.console()


  }


  # Calculate average p-value
  avg_pvals <- m_pvals / l$params$nu_simulations


  # Save full list of hyper parameters?
  if (l$opts$full_params_in_file) {

    file_details <-
      paste(
        paste0(names(l), "_", l),
        collapse = " ")

  } else {

    file_details <-
      paste(
        paste0(names(l$opts$selct_params_in_file),
               "_",

               # Evaluate expressions in current setting
               lapply(l$opts$selct_params_in_file, eval)),
        collapse = " ")

  }


  # Save result
  fwrite(avg_pvals, paste0("Outputs/HH - Village surface/",
                           "Homogeneous effects ",
                           file_details,
                           ".csv"), yaml = T)

}


# Close parallelization cluster
stopCluster(clusters)


# Clean up after loop
rm()

