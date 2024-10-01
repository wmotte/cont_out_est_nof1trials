#!/usr/bin/env Rscript
#
# Wim Otte (w.m.otte@umcutrecht.nl)
# 9 February 2024
#
# Nicotine n-of-1 trial updating
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
    infile <- 'datasets/nicotine.csv'   
    df <- read.csv( infile, stringsAsFactors = TRUE, row.names = 1 )
    df$time <- 1:nrow( df )

    # only n-of-1 part
    df <- df[ df$type == 'blinded', ]
    rownames( df ) <- NULL
    df$type <- NULL
    
    df$treatment <- relevel( as.factor( df$treatment ), 'placebo' )
    
    df$period <- 'first'
    df[ 43:nrow( df ), 'period' ] <- 'second'
    
    return( df )
}

###
# BF
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
    
    
    # output, last row of data
    output <- data[ nrow( data ), ]
    output$log_bf <- BF_brms_bridge$bf
    
    return( output )
}

################################################################################
# CUSTOM FUNCTIONS                                                             #
################################################################################

# output dir
outdir <- 'out.11.BF.nicotine'
dir.create( outdir, showWarnings = FALSE )

# get data
head( df <- get_data() )

#######################################################
# get Bayes Factors
#######################################################

set.seed( 111 )

# get subsets for first and second periods
data_first <- df[ df$period == 'first', ]
data_second <- df

# get bf's
bf1 <- get_bayes_factor( data_first )
bf2 <- get_bayes_factor( data_second )

container <- rbind( bf1, bf2 )

# print [4.3, 18.0]
print( round( container$log_bf, 1 ) )

# write to file
readr::write_tsv( container, file = paste0( outdir, '/bfs_cycles.tsv' ), quote = 'all' )
