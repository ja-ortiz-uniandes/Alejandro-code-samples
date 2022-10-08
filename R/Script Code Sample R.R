
# Code sample R

# By: Alejandro Ortiz - ja.ortiz@uniandes.edu.co


# Make project Portable and reproducible
renv::activate()



# Set up                                                                    ----


# Clean - R's enviornment
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
library(fixest)
library(tidyverse)
library(data.table)



# Clustered vs. Stratified vs. Systematic Sampling                          ----


# This is the an excercice that models a question from DIME.


# Create population




tot.m <- matrix(nrow = 50, ncol = 20)


mat <- matrix(nrow = 50, ncol = 20)


ncores <- detectCores(logical = F)
cl <- makeCluster(ncores)


for (i in 1:1000) {


  set.seed(1944 + i) # Year of Bretton woods

  vil <- data.table()

  for (v in 1:150) {


    v.size <- sample(30:200, size = 1)

    treatment <- sample(c(T, F), 1)

    tmp_vil <- data.table("village" = v,

                          "household" = 1:v.size,

                          "baseline" = rnorm(v.size,
                                             mean = runif(1, 0, 100),
                                             sd = runif(1, 0.5, 3)),

                          "treatead" = treatment %>% rep(v.size))


    vil <- rbind(vil, tmp_vil)


  }

  print("fin crear muestra")
  print(Sys.time())

  # vil[, .("Nu. of villages" = uniqueN(.SD, by = "village")), by = .(treatead)]



  # vil[, outcome := baseline]
  # vil[treatead == 1, outcome := baseline + rnorm(.N, mean = 0.1)]

  vil[, baseline := baseline + rnorm(.N, mean = 0.1)]


  set.seed(2005 + i) # Year DIME was created



  hh.size.village <- function(nu.vil, db) {

    require(data.table)
    require(fixest)

    tmp <- vector(length = 20)


    for (nu.hh in 1:20) {

      # Sample of treated villages
      v.sample <- db[treatead == 1, unique(.SD), .SDcols = "village"
      ][sample(.N, size = nu.vil)]


      # Sample of untreated villages
      v.sample <- rbind(
        v.sample,
        db[treatead != 1, unique(.SD), .SDcols = "village"
        ][sample(.N, size = nu.vil)])





      tmp[nu.hh] <-
        feols(baseline ~ treatead,
              data = db[v.sample, on = "village"
              ][, .SD[sample(.N, size = min(.N, nu.hh))],
                by = village]
        )[["coeftable"]][["Pr(>|t|)"]][2]



    }

    return(tmp)
  }


  p.result <- clusterApply(cl, x = 1:50,
                           fun = hh.size.village, db = vil)

  mat <- p.result %>% unlist %>% matrix(nrow = 50)

  mat[is.na(mat)] <- 0

  tot.m <- tot.m + mat

  print(paste("fin ciclo", i))
  print(Sys.time())

}

k <- tot.m / 102

k <- cbind(NA, k)
k <- rbind(NA, k)

library(plotly)

# volcano is a numeric matrix that ships with R

fig <- plot_ly(z = ~k)

fig <- fig %>% add_surface() %>%
  layout(
    scene = list(
      xaxis = list(title = 'Nu. of households per village'),
      yaxis = list(title = 'Nu. of villages per group'),
      zaxis = list(title = 'P-value')
    )
  )


fig


tot.m[is.na(tot.m)] <- 0
