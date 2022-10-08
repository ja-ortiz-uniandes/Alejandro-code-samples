
# Code sample R - Simulation

# By: Alejandro Ortiz - ja.ortiz@uniandes.edu.co

# This is an exercise that models a question from DIME.


# IMPORTANT !!
# Please look at the HTML report before looking at this script





# Set up                                                                    ----


# Clean - R's environment
# .rs.restartR()
cat("\f")
graphics.off()
remove(list = ls())
gc(full = T)



# Publish working directory
getwd()


# Set options
# options(java.parameters = "-Xmx8000m")
options(max.print = 200)


# Update and load packages
# update.packages(ask = F)

# Plot results
library(plotly)

# Parallelization
library(furrr)

# Estimation
library(fixest)

# Core
library(tidyverse)
library(data.table)



# Hyper-parameter dashboard                                                 ----


# Hyper-parameters list template
l <- list()


# Maximum number of households sampled per village
l$params$max_hh_sample_per_vil <- 50


# Maximum number of villages sampled per treatment group
l$params$max_vils_sampled <- 50


# Minimum effect size
l$params$min_effect_size <- 0.1


# Maximum effect size
l$params$max_effect_size <- 0.3


# Effect size step
l$params$effect_size_step <- 0.05


# Number of times to run simulation over the same sample parameters
l$params$nu_simulations <- 10^4


# Minimum number of villages in the population
l$params$min_vils_in_universe <- 200


# Minimum number of HH per village in population
l$params$min_vil_size <- 50


# Maximum number of HH per village in population
l$params$max_vil_size <- 200


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
l$opts$select_params_in_file <- list("eff" = quote(eff_size))



# in future - Warnings and error handling




# Simulations                                                               ----


# Function to generate a random village
gen_village <- function() {


  # Select a village size
  vil_size <- sample(l$params$min_vil_size:l$params$max_vil_size, size = 1)


  # Select if village will be treated or not
  treatment <- sample(0:1,
    size = 1,
    prob = c(
      1 - l$params$prob_vil_is_treated,
      l$params$prob_vil_is_treated
    )
  )


  # Create a village data set
  tmp_vil <- data.table(

    # Village ID
    "village" = village,

    # Within-village Household ID
    "household" = seq_len(vil_size),

    # Create baseline score
    "baseline" = rnorm(vil_size,

      # Mean depends on the village
      mean = runif(1,
        min = l$params$bl_min_mean_val,
        max = l$params$bl_max_mean_val
      ),

      # SD depends on village
      sd = runif(1,
        min = l$params$bl_min_sd_val,
        max = l$params$bl_max_sd_val
      )
    ),


    # Treatment status - village-wide effect
    "treated" = treatment %>% rep(vil_size)
  )


  # Concatenate village data set with all previous villages
  return(tmp_vil)
}


## Parallelization


# Parallelization function
change_hh_size <-
  function(nu_of_villages, universe,
           max_nu_hh_per_vil = l$params$max_hh_sample_per_vil,
           return_messages = F, return_warnings = F) {

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
    pval_results <- vector("numeric", length = max_nu_hh_per_vil)


    # Loop over HHs sampled in each village
    for (nu_hh in seq_len(max_nu_hh_per_vil)) {

      # Sample of treated villages
      v_sample <- universe[treated == 1, unique(.SD), .SDcols = "village"][sample(.N, size = nu_of_villages)]


      # Sample of untreated villages
      v_sample <- rbind(
        v_sample,
        universe[treated != 1, unique(.SD), .SDcols = "village"][sample(.N, size = nu_of_villages)]
      )


      # Temporary result template
      tmp <- 0
      # Used so furrr - can transport the object


      # From a simple OLS - get the p-value for treated dummy


      tryCatch(
        expr = {
          tmp <<-
            feols(outcome ~ baseline + treated,
              data = universe[

                # Only villages in sample
                v_sample,
                on = "village"

                # Sample a of hh in each village
              ][, .SD[sample(.N, size = min(.N, nu_hh))],
                by = village
              ]

              # Grab p-value for treated dummy
            )[["coeftable"]][["Pr(>|t|)"]][3]
        },
        warning = function(war) {

          # Return warning if necessary
          if (return_warnings) {
            war
          }
        },
        message = function(mes) {

          # Return message if necessary
          if (return_messages) {
            mes
          }
        }
      )

      # Save to result vector
      pval_results[nu_hh] <- tmp
    }

    # Result is the vector of p-values - for that Nu. of villages in each
    # treatment group
    return(pval_results)
  }


# plan for future processes
plan(multisession)


# Loop over effect size
for (
  eff_size in seq(
    from = l$params$min_effect_size,
    to = l$params$max_effect_size,
    by = l$params$effect_size_step
  )
) {


  # Pre-allocate memory to results
  m_pvals <- matrix(
    nrow = l$params$max_vils_sampled,
    # row number is equal to number of villager per treatment group

    ncol = l$params$max_hh_sample_per_vil
  )
  # column number is equal to number of HHs per village


  # Matrix must not be NA or operations will only yield NA
  m_pvals[is.na(m_pvals)] <- 0


  # Simulate results many times
  for (iter in seq_len(l$params$nu_simulations)) {

    # Set a fixed seed for reproducibility
    set.seed(1944 + iter) # Year of Bretton woods


    # Template of village-household universe
    vil_hh_u <- data.table()


    # Creating the universe - loop over villages
    for (village in seq_len(l$params$min_vils_in_universe)) {
      vil_hh_u <- rbind(vil_hh_u, gen_village())
    }


    ## Contingency in case treatment or control groups are too small

    # Size of the smallest treatment group
    smallest_group <-
      vil_hh_u[, .(nu_vils = uniqueN(.SD, by = "village")), by = treated][, min(nu_vils)]


    # Guarantee that there are enough villages in both treatment groups
    while (smallest_group < l$params$max_vils_sampled) {
      vil_hh_u <- rbind(vil_hh_u, gen_village())

      # Size of the smallest treatment group
      smallest_group <-
        vil_hh_u[, .(nu_vils = uniqueN(.SD, by = "village")), by = treated][, min(nu_vils)]
    }


    # Homogeneous effects of treatment
    vil_hh_u[, outcome := baseline + rnorm(.N)]
    vil_hh_u[treated == 1, outcome := baseline + rnorm(.N, mean = eff_size)]
    # Treatment magnitude is the mean value


    # New seed for stage-specific reproducibility
    set.seed(2005 + iter) # Year DIME was created


    ## Parallelize

    # Runs for different village sample sizes the effect of varying
    # household sample size
    p_result <- future_map(seq_len(l$params$max_vils_sampled),
      .f = change_hh_size,
      universe = vil_hh_u,
      .progress = T,
      .options = furrr_options(seed = T)
    )


    # Unlist results into a matrix
    tmp_mat <- p_result %>%
      unlist() %>%
      matrix(ncol = l$params$max_hh_sample_per_vil, byrow = T)
    # columns represent number HHs sampled per village,
    # rows represent the number of villages per treatment group


    # Warn if there are NAs, then replace NA's with 0
    if (any(is.na(tmp_mat))) {
      warning(
        paste0(
          "NAs detected for effect size: ", eff_size,
          " Iteration: ", iter
        )
      )
    }
    tmp_mat[is.na(tmp_mat)] <- 1


    # Add results from previous iterations
    m_pvals <- m_pvals + tmp_mat


    # Verbalize iteration and time
    cat(
      paste0(
        "Finished cycle: ", iter, " - for effect size: ", eff_size,
        "\n", Sys.time(), "\n\n"
      )
    )

    flush.console()
  }


  # Calculate average p-value
  avg_pvals <- m_pvals / l$params$nu_simulations


  # Sample size of 2 is too small
  avg_pvals[1, 1] <- NA


  # Save full list of hyper parameters?
  if (l$opts$full_params_in_file) {
    file_details <-
      paste(
        paste0(names(l), "_", l),
        collapse = " "
      )
  } else {
    file_details <-
      paste(
        paste0(
          names(l$opts$select_params_in_file),
          "_",

          # Evaluate expressions in current setting
          lapply(l$opts$select_params_in_file, eval)
        ),
        collapse = " "
      )
  }


  # Transform to data.table
  avg_pvals <- avg_pvals %>% as.data.table()


  # Bulk rename
  names(avg_pvals) <- paste("HH_per_vil", seq_len(l$params$max_hh_sample_per_vil),
    sep = "."
  )


  # Add variable indicating the number of villages
  avg_pvals[, nu_villages := seq_len(.N)]


  # Save result
  fwrite(avg_pvals,
    paste0(
      "Outputs/HH - Village surface/",
      "Homogeneous effects ",
      file_details,
      ".csv"
    ),
    yaml = T
  )
}


# End parallelization
plan(sequential)


# Clean up after loop
remove(list = ls()[!ls() == "l"])




# Interactive plotly results                                                ----


# Plot list template
pval_plots <- list()


# Generate plots by looping over effect size
for (file in list.files("Outputs/HH - Village surface/",
  pattern = "Homogeneous effects",
  full.names = T
)) {


  # Import data
  avg_pvals <- fread(file, drop = "nu_villages", yaml = T)
  # nu_villages is not imported because it's not needed


  # Get hyper parameters form the data
  local.max_hh_sample_per_vil <- NCOL(avg_pvals)
  local.max_vils_sampled <- NROW(avg_pvals)


  # Add an empty row and column
  avg_pvals <- rbind(NA, avg_pvals, fill = T)
  # this is done because the plotting function starts at 0 by default not 1
  # The first value is estimated using 1 hh in 1 village.


  # List of x values
  xvals <-
    rep(
      0:local.max_vils_sampled,
      each = (local.max_hh_sample_per_vil + 1)
    ) %>%
    matrix(ncol = local.max_hh_sample_per_vil + 1, byrow = T)


  # List of y values
  yvals <-
    rep(
      0:local.max_hh_sample_per_vil,
      (local.max_vils_sampled + 1)
    ) %>%
    matrix(ncol = local.max_hh_sample_per_vil + 1, byrow = T)


  ## Matrix of significance

  # Template
  sig_star <- matrix(
    rep("", local.max_hh_sample_per_vil * local.max_vils_sampled),
    ncol = local.max_hh_sample_per_vil
  )


  # Include significance
  sig_star[as.matrix(avg_pvals) < 0.10] <- "*"
  sig_star[as.matrix(avg_pvals) < 0.05] <- "**"
  sig_star[as.matrix(avg_pvals) < 0.01] <- "***"


  # Extract effect size from title
  local.eff_size <-
    gsub(
      pattern = ".*eff_(\\d*\\.\\d*).*",
      replacement = "\\1",
      x = file,
      perl = T
    )


  # Actual plot
  pval_plots[[paste0("eff_", local.eff_size)]] <-
    plot_ly(z = as.matrix(avg_pvals)) %>%
    # Surface plot
    add_surface(

      # Custom information on hover
      hovertext = paste(

        paste0(
          "Nu. of Villages: ", xvals
        ),
        paste0(
          "<br>H.H. per Village: ", yvals
        ),
        paste0(
          "<br>p-value: ", paste0(
            round(avg_pvals %>% as.matrix(), 3),
            sig_star
          )
        ),
        paste0(
          "<br>Sample size: ", xvals * yvals * 2
        ),
        sep = " "
      ) %>%
        matrix(ncol = local.max_hh_sample_per_vil + 1) %>% t(),
      # plotly assigns values in a different order in hovertext than in a surface
      # trace. To account for this the matrix must be transposed.

      hovertemplate = "%{hovertext}<extra></extra>",
      showscale = F
    ) %>%
    # Layout options
    layout(
      hoverlabel = list(namelength = 10L),

      # Add margins to the title
      margin = list(
        l = 50,
        r = 50,
        b = 50,
        t = 50,
        pad = 20
      ),

      # Graph title
      title = list(
        text = paste0(
          "P-value for different sample size distributions",
          " - effect size: ", local.eff_size
        )
      ),
      scene = list(

        # x-axis options
        xaxis = list(
          title = list(
            text = "Households per village",
            font = list(
              size = 12
            )
          )
        ),

        # y-axis options
        yaxis = list(
          title = list(
            text = "Villages per treatment group",
            font = list(
              size = 12
            )
          )
        ),

        # z-axis options
        zaxis = list(title = "P-value"),

        # Camera options
        camera = list(
          center = list(
            x = 0,
            y = 0,
            z = -0.3
          ),
          eye = list(
            x = 1.4,
            y = 1.4,
            z = 1
          )
        )
      )
    )
}



# Clean up after loop
remove(list = ls()[!ls() == "pval_plots"])
