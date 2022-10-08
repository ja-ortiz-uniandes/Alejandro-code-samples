
# Code sample R

# By: Alejandro Ortiz - ja.ortiz@uniandes.edu.co


# IMPORATANT
# Please look at the HTML report before looking at this script

# This project uses renv to so it is portable and reproducible


# This is the an exercise that models a question from DIME.

# Clustered vs. Stratified vs. Systematic Sampling


# Set up                                                                    ----


# Clean - R's environment
# .rs.restartR()
cat("\f")
# dev.off()
remove(list = ls())
gc(full = T)


# Publish working directory
getwd()


# Set options
# options(java.parameters = "-Xmx8000m")
options(max.print = 200)


# Update and load packages
# update.packages(ask = F)
library(plotly)
library(furrr)
library(fixest)
library(tidyverse)
library(data.table)



# Hyper-parameter dashboard                                                 ----


# Hyper-parameters list template
l <- list()


# Maximum number of households sampled per village
l$params$max_hh_sample_per_vil <- 50


# Maximum number of villages sampled per treatment group
l$params$max_vil_sampled <- 50


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
# l$params$full_params_in_file == F
l$opts$selct_params_in_file <- list("eff" = quote(eff_size))



# Warnings and error handling coming in the future




# Simulations                                                               ----


# Function to generate a random village
gen_village <- function() {


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
  return(tmp_vil)

}


## Parallelization

# plan for future processes
plan(multisession)


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
    pval_results <- vector(length = max_nu_hh_per_vil)


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


      tryCatch(
        expr = {
          tmp <<-
            feols(outcome ~ baseline + treatead,
                  data = universe[

                    # Only villages in sample
                    v_sample, on = "village"

                    # Sample a of hh in each village
                  ][, .SD[sample(.N, size = min(.N, nu_hh))],
                    by = village]

                  # Grab p-value for treated dummy
            )[["coeftable"]][["Pr(>|t|)"]][3]
        },

        warning = function(war) {}
      )

      # Save to result vector
      pval_results[nu_hh] <- tmp

    }

    # Result is the vector of p-values - for that Nu. of villages in each
    # treatment group
    return(pval_results)
  }


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
    for (village in 1:l$params$min_vils_in_universe) {

      vil_hh_u <- rbind(vil_hh_u, gen_village())

    }


    ## Contingency in case treatment or control groups are too small

    # Size of the smallest treatment group
    smallest_group <-
      vil_hh_u[, .(nu_vils = uniqueN(.SD, by = "village")), by = treatead
      ][, min(nu_vils)]


    # Guarantee that there are enough villages in both treatment groups
    while (smallest_group < l$params$max_vil_sampled) {

      vil_hh_u <- rbind(vil_hh_u, gen_village())

      # Size of the smallest treatment group
      smallest_group <-
        vil_hh_u[, .(nu_vils = uniqueN(.SD, by = "village")), by = treatead
        ][, min(nu_vils)]

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
    p_result <- future_map(1:l$params$max_vil_sampled,
                           .f = change_hh_size,
                           universe = vil_hh_u,
                           .progress = T,
                           .options = furrr_options(seed = T))


    # Unlist results into a matrix
    tmp_mat <- p_result %>% unlist %>%
      matrix(ncol = l$params$max_hh_sample_per_vil, byrow = T)
    # columns represent number HHs sampled per village,
    # rows represent the number of villages per treatment group


    # Warn if there are NAs, then replace NA's with 0
    if (any(is.na(tmp_mat))) {
      warning(
        paste0("NAs detected for effect size: ", eff_size,
               " Iteration: ", iter)
      )
    }
    tmp_mat[is.na(tmp_mat)] <- 0


    # Add results from previous iterations
    m_pvals <- m_pvals + tmp_mat


    # Verbalize iteration and time
    cat(
      paste0(
        "\nFinished cycle: ", iter, " - for effect size: ", eff_size,
        "\n", Sys.time(), "\n")
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
  fwrite(avg_pvals %>% as.data.table,
         paste0("Outputs/HH - Village surface/",
                "Homogeneous effects ",
                file_details,
                ".csv"),
         yaml = T)

}


# End parallelization
plan(sequential)


# Clean up after loop
remove(list = ls())



# Interactive plotly results                                                ----


# loop over effect size
for (file in list.files("Outputs/HH - Village surface/",
                        pattern = "Homogeneous effects",
                        full.names = T)) {


  # Import data
  avg_pvals <- fread(paste0("Outputs/HH - Village surface/",
                            "Homogeneous effects ",
                            file_details,
                            ".csv"), yaml = T)


  # # Transform matrix into a long-format data.table
  # # Add variable of number of villages per group
  # avg_pvals[, nu_of_villages_per_group := .I]
  #
  #
  # # Melt to long
  # avg_pvals <- melt(avg_pvals,
  #      id.vars = "nu_of_villages_per_group",
  #      variable.name = "nu_hh_per_vil",
  #      value.name = "p_val",
  #      variable.factor = F)
  #
  #
  # # Transform nu. hh per village into numeric
  # avg_pvals[, nu_hh_per_vil := gsub(pattern = "\\D",
  #                           replacement = "",
  #                           x = nu_hh_per_vil,
  #                           perl = T) %>%
  #     as.numeric]


  # Add an empty row and column
  avg_pvals <- rbind(NA, avg_pvals, fill = T)
  # this is done because the plotting function starts at 0 by default not 1
  # The first value is estimated using 1 hh in 1 village.



  # List of values to be placed row-wise
  xvals <-
    rep(
      0:l$params$max_hh_sample_per_vil,
      each = (l$params$max_vil_sampled + 1)
    )


  # List of values to be placed column-wise
  yvals <-
    rep(
      0:l$params$max_vil_sampled,
      (l$params$max_hh_sample_per_vil + 1)
    )


  # Matrix of significance
  sig_star <- matrix(rep("", dim(avg_pvals)[1] * dim(avg_pvals)[2])) %>%
    matrix(ncol = dim(avg_pvals)[1])

  sig_star[as.matrix(avg_pvals) < 0.10] <- "*"
  sig_star[as.matrix(avg_pvals) < 0.05] <- "**"
  sig_star[as.matrix(avg_pvals) < 0.01] <- "***"


  # Actual plot
  plot_ly(z = ~as.matrix(avg_pvals)) %>%

    # Surface plot
    add_surface(

      hovertext = paste(

        # HH per village text
        "H.H. per Village:", xvals,

        # Nu. of Villages text
        "<br>Nu. of Villages:", yvals,

        # P-value and significance text
        "<br>p-value:", paste0(
          round(avg_pvals %>% as.matrix, 3), sig_star),

        # Hover on sample size
        "<br>Sample size: ", paste0(

          # Multiply both values to get sample per group
          xvals * yvals * 2
          # multiply by 2 because there are 2 groups - treatment & control
        )
      ) %>%

        matrix(ncol = (l$params$max_hh_sample_per_vil + 1)) %>%

        # Plotly maps values in reverse order so one has to transpose
        # The resulting matrix for proper mapping
        t,

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
        text = paste0("P-value for different sample size distributions",
                      " - effect size: ", eff_size)
      ),

      scene = list(

        # x-axis options
        xaxis = list(
          title = list(
            text = 'Households per village',
            font = list(
              size = 12
            )
          )
        ),

        # y-axis options
        yaxis = list(
          title = list(
            text = 'Villages per treatment group',
            font = list(
              size = 12
            )
          )
        ),

        # z-axis options
        zaxis = list(title = 'P-value'),

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




























