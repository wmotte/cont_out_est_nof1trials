#!/usr/bin/env Rscript
#
# Wim Otte (w.m.otte@umcutrecht.nl)
# 16 Jan 2024
#
# CBD-trial simulation
################################################################################
library( 'readr' ) # read and write (large) tsv-files
library( 'rstanarm' ) # pre-compiled Bayesian models
library( 'bayestestR' ) # ROPE function

# set readr progress bar to silent mode
options( readr.show_progress = FALSE )

################################################################################
# CUSTOM FUNCTIONS                                                             #
################################################################################

###
# Process model
##
process_model <- function( datafile, outfile_base, subject = 1, max_cycle = 1 )
{
    df <- readr::read_tsv( datafile, show_col_types = FALSE )
    
    # select single subject
    df <- df[ df$subject == paste0( 'subject-', subject ), ]
    df$subject <- NULL
    
    # check input boundaries
    if( max_cycle < 1 | max_cycle > 5 )
        stop( "*** error ***: max_cycle should be 1,2 or 3!" )
    
    # select baseline up to max_cycle
    df <- df[ df$cycle %in% c( 'prior', paste0( 'cycle-', 1:max_cycle ) ), ]
    
    # convert to factor
    df$cycle <- as.factor( df$cycle )
    df$group <- as.factor( df$group )
    
    # set reference 'prior'
    df <- within( df, cycle <- relevel( cycle, ref = 'prior' ) )
    
    # fit model and plot output
    output <- fit_bayesian_poisson_model( df, max_cycle, outfile_base )
    return( output )
}

###
# Get certainty of drop in outcome (ROPE method)
# https://easystats.github.io/bayestestR/articles/region_of_practical_equivalence.html
##
get_drop_certainty <- function( model, outfile_base, prob_drop = 0.5 )
{
    # get baseline coefficient (lambda, i.e., average daily seizures)
    baseline_value <- model$coefficients[ '(Intercept)' ]
    
    # get lower ROPE value [100% reduction - 50%/30% reduction]
    upper_rope_value <- -1 * prob_drop * baseline_value
    lower_rope_value <- -1 * baseline_value
    
    # get perc in ROPE area
    perc_in_rope <- bayestestR::rope( model, ci = 1.0, range = c( lower_rope_value, upper_rope_value ) )
    
    return( perc_in_rope )
}

###
# Fit Bayesian Poisson model
##
fit_bayesian_poisson_model <- function( df, max_cycle = 1, outfile_base )
{
    # fit Bayesian Poisson model (group)
    model <- rstanarm::stan_glm( y ~ group, data = df,
                           chains = 1, cores = 1, iter = 5000, refresh = 0 )

    output <- NULL
    
    # calculate and save drop percentages with ROPE method
    for( prob_drop in c( 0.3 ) )
    {
        # get ROPE
        perc_in_rope <- get_drop_certainty( model, outfile_base, prob_drop )
        
        # get table with percentage certainty
        tt <- as.data.frame( perc_in_rope )
        
        perc_certainty_placebo <- round( 100 * tt[ tt$Parameter == 'groupplacebo', 'ROPE_Percentage' ], 1 )
        perc_certainty_verum <- round( 100 * tt[ tt$Parameter == 'groupverum', 'ROPE_Percentage' ], 1 )
        
        # output
        res <- data.frame( prob_drop, perc_certainty_placebo, perc_certainty_verum )
        
        output <- rbind( output, res )
    }
    
    return( output )
}


################################################################################
# END CUSTOM FUNCTIONS                                                         #
################################################################################

# specify seed for random number generation 
set.seed( 888 )

# data input directory
indir <- 'out.02.simulate.refractory'

# output directory
outdir <- 'out.03.model.refractory'
dir.create( outdir, showWarnings = FALSE )

# which subject to process? [for now single subject n-of-1]
csubject <- 1

# meta info on effect sizes
cweekly_seizures_baseline <- c( 35 )

# probabilities of drop in seizure risk
cprob_drop_verum <- c( 0.6 )
cprob_drop_placebo <- c( 0.1 )

# maximum cycle to process (if == 1, only first cycle is processed)
cmax_cycle <- c( 1, 2, 3, 4, 5 )

all <- NULL

for( subject in csubject )
{
    # loop over 4 option sets
    for( max_cycle in cmax_cycle ) 
    {
        print( paste0( 'subject: ', subject, ' -- max_cycle: ', max_cycle ) )
        
        for( weekly_seizures_baseline in cweekly_seizures_baseline )
        {
            for( prob_drop_verum in cprob_drop_verum )
            {
                for( prob_drop_placebo in cprob_drop_placebo )
                {
                    datafile <- paste0( indir, '/data__weekly-seizures-baseline-', 
                                     weekly_seizures_baseline, '__placebo-verum-risk-', 
                                     prob_drop_placebo, '--', prob_drop_verum, '.tsv.gz' )
                    
                    # only proceed if datafile is available
                    if( file.exists( datafile ) )
                    {
                        # process
                        single <- process_model( datafile, outfile_base, subject, max_cycle )
                        single$subject <- subject
                        single$max_cycle <- max_cycle
                        single$weekly_seizures_baseline <- weekly_seizures_baseline
                        single$prob_drop_verum <- prob_drop_verum
                        single$prob_drop_placebo <- prob_drop_placebo
                        
                        all <- rbind( all, single )
                    }
                }
            }
        }
    }
}

# write to disk
readr::write_tsv( all, file = gzfile( paste0( outdir, '/all_stats.tsv.gz' ) ), quote = 'all' )
