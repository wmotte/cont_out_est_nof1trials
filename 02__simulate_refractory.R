#!/usr/bin/env Rscript
#
# Wim Otte (w.m.otte@umcutrecht.nl)
# 16 Jan 2024
#
# Simulate refractory epilepsy subject
################################################################################
library( 'readr' )

# set readr progress bar to silent mode
options( readr.show_progress = FALSE )

################################################################################
# CUSTOM FUNCTIONS                                                             #
################################################################################

###
# get study design
##
get_study_design <- function( ncycles, ndays_baseline, ndays_verum, ndays_placebo )
{
    data <- NULL
    
    # first add baseline
    single <- data.frame( 
        cycle = rep( 'prior', ndays_baseline ),
        group = c( rep( 'baseline', ndays_baseline ) ) )

    # add to container
    data <- rbind( data, single )
    
    # now add cycles
    for( ncycle in 1:ncycles )
    {
        single <- data.frame( 
            cycle = rep( paste0( 'cycle-', ncycle ), length( single ) ),
            group = c( rep( 'placebo', ndays_placebo ), rep( 'verum', ndays_verum ) ) )
        
        # add to container
        data <- rbind( data, single )
    }
    
    # add day index for full design
    data$day <- 1:nrow( data )
    
    return( data )
}

###
# Add lambda (used for Poisson seizure simulation) to each day
##
add_lambda_to_design <- function( design, weekly_seizures_baseline, prob_drop_verum, prob_drop_placebo )
{
    # get daily seizures
    daily_seizures_baseline <- weekly_seizures_baseline / 7
    
    # get daily seizures after drop in risk
    daily_seizures_placebo <- daily_seizures_baseline - ( daily_seizures_baseline * prob_drop_placebo )
    daily_seizures_verum <- daily_seizures_baseline - ( daily_seizures_baseline * prob_drop_verum )

    # prepare
    design$lambda <- NA     
    
    design[ design$group %in% 'baseline', 'lambda' ] <- daily_seizures_baseline
    design[ design$group %in% 'placebo', 'lambda' ] <- daily_seizures_placebo
    design[ design$group %in% 'verum', 'lambda' ] <- daily_seizures_verum
    
    return( design )
}

###
# Simulate y
##
simulate_y <- function( design, nsubjects )
{
    
    data <- NULL
    
    for( nsubject in 1:nsubjects )
    {
        sdesign <- design
        sdesign$subject <- paste0( 'subject-', nsubject )
        sdesign$y <- NA
        
        
        for( i in 1:nrow( sdesign ) )
        {
            sdesign[ i, 'y' ] <- rpois( 1, lambda = sdesign[ i, 'lambda' ] ) 
            
        }
        
        data <- rbind( data, sdesign )
    }
    
    return( data )
}

################################################################################
# END CUSTOM FUNCTIONS                                                         #
################################################################################

# specify seed for random number generation 
set.seed( 777 )

# output directory
outdir <- 'out.02.simulate.refractory'
dir.create( outdir, showWarnings = FALSE )

# meta
nsubjects <- 1
ncycles <- 5
ndays_baseline <- 7 * 8
ndays_verum <- 7 * 4
ndays_placebo <- 7 * 4

# get study design (for single n-of-1 trial)
design <- get_study_design( ncycles, ndays_baseline, ndays_verum, ndays_placebo )

# meta info on effect sizes [5 per day, 35 per week]
cweekly_seizures_baseline <- c( 35 )

# probabilities of drop in seizure risk
cprob_drop_verum <- c( 0.6 )
cprob_drop_placebo <- c( 0.1 )

# loop over 3 option sets
for( weekly_seizures_baseline in cweekly_seizures_baseline )
{
    for( prob_drop_verum in cprob_drop_verum )
    {
        for( prob_drop_placebo in cprob_drop_placebo )
        {
            designl <- NULL
            data <- NULL
            
            print( paste0( '*** simulation ***: ', weekly_seizures_baseline, ' - ', prob_drop_verum, ' - ', prob_drop_placebo ) )
            designl <- add_lambda_to_design( design, weekly_seizures_baseline, prob_drop_verum, prob_drop_placebo )
            
            # simulate seizures based on design
            data <- simulate_y( designl, nsubjects )
            
            # save data to file
            readr::write_tsv( data, file = gzfile( paste0( outdir, '/data__weekly-seizures-baseline-', 
                                                   weekly_seizures_baseline, '__placebo-verum-risk-', 
                                                   prob_drop_placebo, '--', prob_drop_verum, '.tsv.gz' ) ), 
                                                    quote = 'all' )            
        }
    }
}

