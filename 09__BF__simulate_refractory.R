#!/usr/bin/env Rscript
#
# Wim Otte (w.m.otte@umcutrecht.nl)
# March 2024
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
    infile <- 'out.02.simulate.refractory/data__weekly-seizures-baseline-35__placebo-verum-risk-0.1--0.6.tsv.gz' 
    df <- as.data.frame( readr::read_tsv( infile, show_col_types = FALSE ) )
    
    # make null_group (with only 'baseline' and 'post-baseline'
    df$groupnull <- 'postbaseline'
    df[ df$group == 'baseline', 'groupnull' ] <- 'baseline'
    
    return( df )
}


###
# Get bayes factor
##
get_bayes_factor <- function( data )
{
    
    # fit Bayesian Poisson model [H1]
    # refresh = 0 turns off iteration info print
    bmodel_H1 <- brms::brm( y ~ group, data = data, 
                            iter = 50000, warmup = 3000,
                            family = poisson(), 
                            save_pars = save_pars( all = TRUE ), refresh = 0 )
    
    # update Bayesian Poisson model [H0]
    bmodel_H0 <- update( bmodel_H1, formula = y ~ groupnull, newdata = data, refresh = 0 ) 
    
    # get summary and parameters
    prior_summary( bmodel_H1 )
    prior_summary( bmodel_H0 )
    
    posterior_summary( bmodel_H1 )
    posterior_summary( bmodel_H0 )
    
    summary( bmodel_H1 )
    summary( bmodel_H0 )
    
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
set.seed( 999 )

outdir <- 'out.09.BF.simulate.refractory'
dir.create( outdir, showWarnings = FALSE )

# get data
head( df <- get_data() )

####################################

# loop over data and run model as soon as cycle is about to switch
current_cycle <- 'cycle-1'
i <- 110 # somewhere at end of cycle-1

# BF container
container <- NULL

for( i in 110:( nrow( df ) - 1 ) )
{
    dfcurrent <- df[ i, ]
    dfnext <- df[ i + 1, ]
    
    # check if next day is new cycle
    if( ( dfcurrent$cycle != dfnext$cycle ) | ( i + 1 == nrow( df ) ) )
    {
        print( paste0( "Finalizing day before new cycle: ", i ) )
        
        # evaluation data
        data <- df[1:i, ]
        
        # get Bayes Factor
        current_bf <- get_bayes_factor( data )
        
        # add to container
        container <- rbind( container, current_bf )

        print( container )
    }
}

print( round( container$log_bf, 1 ) )
readr::write_tsv( container, file = paste0( outdir, '/bfs_cycles.tsv' ), quote = 'all' )
 


## additional Bayes Factor for first day the posterior probability is clinically relevant.
data <- df[ 1:152, ]
bf <- get_bayes_factor( data )
readr::write_tsv( bf, file = paste0( outdir, '/bf_at_day_152.tsv' ), quote = 'all' )
