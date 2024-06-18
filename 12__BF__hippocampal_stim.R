#!/usr/bin/env Rscript
#
# Wim Otte (w.m.otte@umcutrecht.nl)
# 9 February 2024
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
    # fit Bayesian Poisson model (treatment)
    bmodel_H1 <- brms::brm( seizures ~ treatment, data = data, 
                            iter = 25000, warmup = 3000,
                            family = poisson(), 
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


outdir <- 'out.12.BF.hippocampal.stim'
dir.create( outdir, showWarnings = FALSE )


# get data
head( df <- get_data() )

# make subsets
data_set_1 <- df[ df$month <= 12, ]
data_set_2 <- df[ df$month <= 24, ]
data_set_3 <- df[ df$month <= 48, ]

# get BFs
bf1 <- get_bayes_factor( data_set_1 )
bf2 <- get_bayes_factor( data_set_2 )
bf3 <- get_bayes_factor( data_set_3 )

# combine
container1 <- rbind( bf1, bf2 )
container <- rbind( container1, bf3 )

# print
print( round( container$log_bf, 1 ) )

# write to file
readr::write_tsv( container, file = paste0( outdir, '/bfs_cycles.tsv' ), quote = 'all' )


