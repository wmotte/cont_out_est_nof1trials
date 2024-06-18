#!/usr/bin/env Rscript
#
# Wim Otte (w.m.otte@umcutrecht.nl)
# 16 January 2023
#
# CBD-trial simulation
################################################################################
library( 'readr' ) # read and write (large) tsv-files
library( 'ggplot2' ) # plot
library( 'rstanarm' ) # pre-compiled Bayesian models
library( 'bayestestR' ) # ROPE function
library( 'see' ) # plotting

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
# Process plot
##
process_plot <- function( datafile, outfile_plot, subject = 1, max_cycle = 1 )
{
    df <- readr::read_tsv( datafile, show_col_types = FALSE )
    
    # select single subject
    df <- df[ df$subject == paste0( 'subject-', subject ), ]
    
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
    
    # get breaks/labels for x-axis at weekly basis
    lab_y <- seq( from = 0, to = max( df$day ), by = 7 )
    
    # plot outcome over time
    p <- ggplot( data = df, aes( group = group, fill = group, x = day, y = y ) ) + 
        geom_point() + # for the fat-black strokes
        geom_point( shape = 21, colour = 'gray40', stroke = 0.2 ) +
        scale_x_continuous( breaks = lab_y, labels = lab_y ) +
        scale_y_continuous( breaks = number_ticks( 6 ) ) +
        facet_wrap( ~cycle, ncol = 3, scale = 'free_x' ) +
        xlab( 'Time' ) +
        ylab( 'Seizures (per day)' ) +
        scale_colour_manual( values = c( "#999999", "#E69F00", "#56B4E9" ) ) +
        scale_fill_manual( values = c( "#999999", "#E69F00", "#56B4E9" ) ) +
 
        theme_classic() +
        theme( legend.position = 'top', strip.text.x = element_text( face = 'bold' ), 
               strip.text.y = element_text( face = 'bold' ) )
    
    # save to disk
    ggsave( plot = p, file = outfile_plot, dpi = 600, width = 8, height = 4 )
}

###
# Process model
##
get_process_df <- function( datafile, outfile_base, subject = 1, max_cycle = 1 )
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
    
    return( df )
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

    # get posterior distribution samples, post warm-up
    draws <- as.data.frame( model )
    
    # only keep 2500 (should be sufficient for plotting)
    draws <- draws[ 1:2500, ]
    
    # convert to vectors
    vec_baseline <- draws$`(Intercept)`
    vec_placebo <- draws$groupplacebo
    vec_verum <- draws$groupverum
    
    # get percentage change relative to baseline
    perc_drop_placebo <- 100 * ( vec_placebo / vec_baseline )
    perc_drop_verum <- 100 * ( vec_verum / vec_baseline )
    
    # drop, corrected for placebo effect    
    perc_drop_verum__minus_placebo <- perc_drop_verum - perc_drop_placebo
    
    # collect output    
    out <- data.frame( max_cycle, vec_baseline, vec_placebo, vec_verum, 
                       perc_drop_placebo, perc_drop_verum, 
                       perc_drop_verum__minus_placebo )
    
    # save to disk
    readr::write_tsv( out, file = gzfile( paste0( outfile_base, '__posteriors.tsv.gz' ) ), quote = 'all' )
    
    # make plot of placebo-effect-corrected verum drop in outcome
    pdata <- data.frame( value = out$perc_drop_verum__minus_placebo, x = 1:nrow( out ), name = 'verum' )
    
    # get custom ticks breaks at -100, -90, ..., 0
    lab_x <- seq( from = 0, to = 50, by = 1 )
    lab_y <- seq( from = 200, to = -300, by = -10 )

    # plot
    p <- ggplot( pdata ) + 
        geom_histogram( aes( y = value, x = after_stat( 100 * ( count / sum( count ) ) ) ), colour = 'gray20', fill = '#56B4E9' ) +
        geom_hline( yintercept = 0, linetype = "dotted", colour = 'gray20' ) +
        geom_hline( yintercept = -30, linetype = "dashed", colour = 'gray40', size = 0.8 ) +
        geom_hline( yintercept = -50, linetype = "dashed", colour = 'gray40', size = 0.8 ) +
        scale_x_continuous( breaks = lab_x, labels = lab_x ) +
        scale_y_continuous( breaks = lab_y, labels = lab_y ) +
        xlab( 'certainty (%)' ) +
        ylab( 'drop in outcome (%)' ) +
        theme( axis.text.x = element_text( angle = -90, vjust = 0.5, hjust = 0 ) )
    
    # save to disk
    # output file data
    outfile_p <- paste0( outfile_base, '__drop_certainty-placebo-corrected.png' ) 
    ggsave( p, file = outfile_p, dpi = 600, width = 5, height = 5 )
    
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
        
        # output file data
        outfile_res <- paste0( outfile_base, '__drop_certainty-', prob_drop, '.tsv' ) 
        
        # write to tsv
        readr::write_tsv( res, file = outfile_res, quote = 'all' )
    }
}


################################################################################
# END CUSTOM FUNCTIONS                                                         #
################################################################################

# specify seed for random number generation 
set.seed( 888 )

# data input directory
indir <- 'out.02.simulate.refractory'

# output directory
outdir <- 'out.04.plot.refractory'
dir.create( outdir, showWarnings = FALSE )

# which subject to process? [for now single subject n-of-1]
csubject <- c( 1 )

# meta info on effect sizes
cweekly_seizures_baseline <- c( 35 )

# probabilities of drop in seizure risk
cprob_drop_verum <- c( 0.6 )
cprob_drop_placebo <- c( 0.1 )

# maximum cycle to process (if == 1, only first cycle is processed)
cmax_cycle <- c( 5 )

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
                        # process plot, output file
                        outfile_plot <- paste0( outdir, '/subject-', subject, 
                                    '__weekly-seizures-baseline-', 
                                    weekly_seizures_baseline, '__placebo-verum-risk-', 
                                    prob_drop_placebo, '--', prob_drop_verum, '__max-cycle-', max_cycle, '.png' )
                        
 
                        process_plot( datafile, outfile_plot, subject, max_cycle )
                        
                        # process output file, base
                        outfile_base <- gsub( ".png$", "", outfile_plot )
                        
                        # process
                        df <- get_process_df( datafile, outfile_base, subject, max_cycle )
                        
                        # fit model and plot output
                        fit_bayesian_poisson_model( df, max_cycle, outfile_base )
                    }
                }
            }
        }
    }
}

