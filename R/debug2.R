

# Hyper-parameter dashboard                                                 ----


# Hyper-parameters list template
l <- list()


# Maximum number of households sampled per village
l$params$max_hh_sample_per_vil <- 4


# Maximum number of villages sampled
l$params$max_vil_sampled <- 5


# Minimum effect size
l$params$min_effect_size <- 0.1


# Maximum effect size
l$params$max_effect_size <- 0.1


# Effect size step
l$params$effect_size_step <- 0


# Number of times to run simulation over the same sample parameters
l$params$nu_simulations <- 3


# Maximum number of villages in the population
l$params$vil_in_universe <- 20


# Minimum number of HH per village in population
l$params$min_vil_size <- 10


# Maximum number of HH per village in population
l$params$max_vil_size <- 20


# Independent probability of a village being treated
l$params$prob_vil_is_treated <- 0.5


# Minimum value of the mean of baseline score
l$params$bl_min_mean_val <- 2


# Maximum value of the mean of baseline score
l$params$bl_max_mean_val <- 100


# Minimum value of the sd of baseline score
l$params$bl_min_sd_val <- 0.5


# Maximum value of the sd of baseline score
l$params$bl_max_sd_val <- 1.5


## Not Simulation parameters but options
# Save full list of hyper-parameters in file name?
l$opts$full_params_in_file <- F


# Select parameters to include in file name in case
# l$params$full_params_in_file == F
l$opts$selct_params_in_file <- list("eff" = quote(eff_size))



# Warnings and error handling coming in the future








# Parallelization function
change_hh_size <-
  function(nu_of_villages, universe,
           max_nu_hh_per_vil = l$params$max_hh_sample_per_vil) {

    # This functions takes the number of villages to be sampled per treatment
    # group (nu_of_villages) and a data set of the universe of villages and
    # households (universe). Optionally one can also specify the maximum number
    # of hh sampled per village. In this script this value is determined by  the
    # max_hh_sample_per_vil parameter.

    # The function then estimates the effect of varying the number of hh sampled
    # per village and returns is a 1x50 (default) vector of the p-values of
    # each estimation.

    # Required packages for the function
    require(data.table)
    require(fixest)

    # Pre-allocate space for results into memory
    tmp <- vector(length = max_nu_hh_per_vil)


    # Loop over HHs sampled in each village
    for (nu_hh in 1:max_nu_hh_per_vil) {

      # Sample of treated villages
      v_sample <- universe[treatead == 1, unique(.SD), .SDcols = "village"
      ][sample(.N, size = nu_of_villages)]


      # Sample of untreated villages
      v_sample <- rbind(
        v_sample,
        universe[treatead != 1, unique(.SD), .SDcols = "village"
        ][sample(.N, size = nu_of_villages)])


      # From a simple OLS - get the p-value for treated dummy
      tmp[nu_hh] <-

        # OLS
        feols(outcome ~ baseline + treatead,
              data = universe[

                # Only villages in sample
                v_sample, on = "village"

                # Sample a of hh in each village
              ][, .SD[sample(.N, size = min(.N, nu_hh))],
                by = village]

              # Grab p-value for treated dummy
        )[["coeftable"]][["Pr(>|t|)"]][3]

    }


    # Result is the vector of p-values - for that Nu. of villages in each
    # treatment group
    return(tmp)
  }






# change_hh_size(nu_of_villages = 1, universe = vil_hh_u)








# furrr::plan()











