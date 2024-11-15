#!/usr/bin/env Rscript
#
# Wim Otte (w.m.otte@umcutrecht.nl)
# 9 February 2023
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
    infile <- 'datasets/SMS.csv'
    df <- read.csv( infile, stringsAsFactors = TRUE, row.names = 1 )

    df[ df$group == 'injection', 'group' ] <- 'baseline'

    rownames( df ) <- NULL

    return( df )
}


###
# Get BF
##
get_bayes_factor <- function( data )
{
    # truncate gamma and normal priors for the intercept and slope within reasonable
    # severity ranges (0-1 for intercept and -1 to +1 for group)
    priors <- c(
        prior( gamma( 2, 3 ), class = "Intercept", lb = 0, ub = 1 ),
        prior( normal( 0, 5 ), class = "b", lb = -1, ub = 1 )
    )

    # fit Bayesian model (treatment)
    # refresh = 0 turns off iteration info print
    bmodel_H1 <- brms::brm( delta_severity ~ group, data = data,
                            iter = 25000, warmup = 3000,
                            prior = priors,
                            save_pars = save_pars( all = TRUE ), refresh = 0 )

    # update Bayesian model [H0], with intercept only
    bmodel_H0 <- update( bmodel_H1, formula = delta_severity ~ 1, refresh = 0 )

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
  # truncate gamma and normal priors for the intercept and slope within reasonable
  # severity ranges (0-1 for intercept and -1 to +1 for group)
  priors <- c(
    prior( gamma( 2, 3 ), class = "Intercept", lb = 0, ub = 1 ),
    prior( normal( 0, 5 ), class = "b", lb = -1, ub = 1 )
  )

  # fit Bayesian model (treatment)
  # refresh = 0 turns off iteration info print
  data$group2 <- ifelse(data$group == "placebo", -0.5, 0.5)
  bmodel_H1 <- brms::brm( delta_severity ~ group2, data = data,
                          iter = 25000, warmup = 3000,
                          prior = priors,
                          save_pars = save_pars( all = TRUE ), refresh = 0 )

  # update Bayesian model [H0], with intercept only
  bmodel_H0 <- update( bmodel_H1, formula = delta_severity ~ 1, refresh = 0 )

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

# set seed
set.seed( 3333 )

outdir <- 'out.10.BF.headache'
dir.create( outdir, showWarnings = FALSE )

# get data
head( df <- get_data() )

# remove NA's
df <- df[ !is.na( df$delta_severity ), ]

# make subsets
data_day_1 <- df[ df$day <= 1, ]
data_day_2 <- df[ df$day <= 2, ]
data_day_3 <- df[ df$day <= 3, ]
data_day_4 <- df[ df$day <= 4, ]

# get BFs
bf1 <- get_bayes_factor( data_day_1 )
bf2 <- get_bayes_factor( data_day_2 )
bf3 <- get_bayes_factor( data_day_3 )
bf4 <- get_bayes_factor( data_day_4 )

# centered parameterization BFs (verifying that the results are equivalent)
bf1b <- get_bayes_factor2( data_day_1 )
bf2b <- get_bayes_factor2( data_day_2 )
bf3b <- get_bayes_factor2( data_day_3 )
bf4b <- get_bayes_factor2( data_day_4 )
# this is OK

# BIC Bayes factors
(BICbf1 <- bayestestR::bic_to_bf(BIC(lm(delta_severity ~ group, data = data_day_1)), BIC(lm(delta_severity ~ 1, data = data_day_1)), log = FALSE))
(BICbf2 <- bayestestR::bic_to_bf(BIC(lm(delta_severity ~ group, data = data_day_2)), BIC(lm(delta_severity ~ 1, data = data_day_2)), log = FALSE))
(BICbf3 <- bayestestR::bic_to_bf(BIC(lm(delta_severity ~ group, data = data_day_3)), BIC(lm(delta_severity ~ 1, data = data_day_3)), log = FALSE))
(BICbf4 <- bayestestR::bic_to_bf(BIC(lm(delta_severity ~ group, data = data_day_4)), BIC(lm(delta_severity ~ 1, data = data_day_4)), log = FALSE))
# these show much more evidence in favor of the alternative hypothesis, this is worrying as an informed prior should generally provide more evidence
# for the alternative hypothesis (if the alternative is true) than the unit information prior

# another common default option would be JZS (which uses standardized effect sizes)
summary(BAS::bas.lm(delta_severity ~ group, data = data_day_1, prior = "JZS"))
summary(BAS::bas.lm(delta_severity ~ group, data = data_day_2, prior = "JZS"))
summary(BAS::bas.lm(delta_severity ~ group, data = data_day_3, prior = "JZS"))
summary(BAS::bas.lm(delta_severity ~ group, data = data_day_4, prior = "JZS"))
# the JZS prior however does not provide that much evidence in favor of the alternative hypothesis either (I would've expected stronger support for the model 1)

# similarly to the BayesFactor package
BayesFactor::lmBF(delta_severity ~ group, data = data_day_1)
BayesFactor::lmBF(delta_severity ~ group, data = data_day_2)
BayesFactor::lmBF(delta_severity ~ group, data = data_day_3)
BayesFactor::lmBF(delta_severity ~ group, data = data_day_4)

# these information together might indicate that the selected prior distribution is not as bad
# however, it is definitelly much less informative than it could be

container1 <- rbind( bf1, bf2 )
container2 <- rbind( bf3, bf4 )

container <- rbind( container1, container2 )

# print
# log(BF): -0.4, -0.2, 2.1, 4.3
print( round( container$log_bf, 1 ) )

# write to file
readr::write_tsv( container, file = paste0( outdir, '/bfs_cycles.tsv' ), quote = 'all' )

