#!/usr/bin/env Rscript
#
# Wim Otte (w.m.otte@umcutrecht.nl)
# Oct 1, 2024
#
################################################################################
library( 'readr' ) # read and write (large) tsv-files
library( 'brms' )

# set readr progress bar to silent mode
options( readr.show_progress = FALSE )

################################################################################
# CUSTOM FUNCTIONS                                                             #
################################################################################

###
# Helper function plot.
##
number_ticks <- function( n ) { function( limits ) pretty( limits, n ) }

###
# Get data
##
get_data <- function()
{
    # read first tab sheet
    infile <- 'datasets/Tellez-Zenteno.csv'
    df <- read.csv( infile, stringsAsFactors = TRUE, row.names = 1 )
    df$time <- 1:nrow( df )

    rownames( df ) <- NULL

    return( df )
}

###
# Get BF
##
get_bayes_factor <- function( data )
{
    # truncate normal priors for the intercept and slope within reasonable seizure loads
    priors <- c(
        prior( normal( 0, 100 ), class = "Intercept", lb = -10, ub = 10 ),
        prior( normal( 0, 100 ), class = "b", lb = -5, ub = 5 )
    )

    # fit Bayesian Poisson model (treatment)
    bmodel_H1 <- brms::brm( seizures ~ treatment, data = data,
                            iter = 25000, warmup = 3000,
                            family = poisson(),
                            prior = priors,
                            save_pars = save_pars( all = TRUE ), refresh = 0 )

    # update Bayesian Poisson model [H0], with intercept only
    bmodel_H0 <- update( bmodel_H1, formula = seizures ~ 1, refresh = 0 )

    # Bayes Factor, through bridge sampling
    BF_brms_bridge <- brms::bayes_factor( bmodel_H1, bmodel_H0, log = TRUE )

    print( BF_brms_bridge )

    # output, last row of data
    output <- data[ nrow( data ), ]
    output$log_bf <- BF_brms_bridge$bf

    return( output )

}

get_bayes_factor2 <- function( data )
{
  # truncate normal priors for the intercept and slope within reasonable seizure loads
  priors <- c(
    prior( normal( 0, 100 ), class = "Intercept", lb = -10, ub = 10 ),
    prior( normal( 0, log(2)/2 ), class = "b" )  # log(2)/2 ascertains that 95$ of the posterior probability is within the 1/2 to 2x rate reduction range
  )

  # fit Bayesian Poisson model (treatment)
  bmodel_H1 <- brms::brm( seizures ~ treatment, data = data,
                          iter = 25000, warmup = 3000,
                          family = poisson(),
                          prior = priors,
                          save_pars = save_pars( all = TRUE ), refresh = 0 )

  # update Bayesian Poisson model [H0], with intercept only
  bmodel_H0 <- update( bmodel_H1, formula = seizures ~ 1, refresh = 0 )

  # Bayes Factor, through bridge sampling
  BF_brms_bridge <- brms::bayes_factor( bmodel_H1, bmodel_H0, log = TRUE )

  print( BF_brms_bridge )

  # output, last row of data
  output <- data[ nrow( data ), ]
  output$log_bf <- BF_brms_bridge$bf

  return( output )

}

################################################################################
# CUSTOM FUNCTIONS                                                             #
################################################################################

set.seed( 565 )

# out dir
outdir <- 'out.12.BF.hippocampal.stim'
dir.create( outdir, showWarnings = FALSE )


# get data
head( df <- get_data() )

# make subsets
data_set_1 <- df[ df$month <= 6, ]
data_set_2 <- df[ df$month <= 12, ]
data_set_3 <- df[ df$month <= 24, ]
data_set_4 <- df[ df$month <= 28, ]
data_set_5 <- df[ df$month <= 48, ]

# get BFs
bf1 <- get_bayes_factor( data_set_1 ) #  X log(BF)
bf2 <- get_bayes_factor( data_set_2 ) #  X log(BF)
bf3 <- get_bayes_factor( data_set_3 ) #  X log(BF)
bf4 <- get_bayes_factor( data_set_4 ) #  X log(BF)
bf5 <- get_bayes_factor( data_set_5 ) #  X log(BF)

# get BFs
bf1b <- get_bayes_factor2( data_set_1 ) #  X log(BF)
bf2b <- get_bayes_factor2( data_set_2 ) #  X log(BF)
bf3b <- get_bayes_factor2( data_set_3 ) #  X log(BF)
bf4b <- get_bayes_factor2( data_set_4 ) #  X log(BF)
bf5b <- get_bayes_factor2( data_set_5 ) #  X log(BF)

# combine
container1 <- rbind( bf1, bf2 )
container2 <- rbind( bf3, bf4 )
container3 <- rbind( container1, container2 )
container <- rbind( container3, bf5 )

# combine
container1b <- rbind( bf1b, bf2b )
container2b <- rbind( bf3b, bf4b )
container3b <- rbind( container1b, container2b )
containerb <- rbind( container3b, bf5b )
# the evidence here is actually lower than with the previous model as the effect is much more extreme than expected
# this might however change if the treatment was less strong, i.e., only rate of 1.5 or 2


# print
print( round( container$log_bf, 1 ) )
rownames( container ) <- NULL
# write to file
readr::write_tsv( container, file = paste0( outdir, '/bfs_cycles.tsv' ), quote = 'all' )


