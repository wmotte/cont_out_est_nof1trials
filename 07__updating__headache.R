#!/usr/bin/env Rscript
#
# Wim Otte (w.m.otte@umcutrecht.nl)
# 6 February 2023
#
################################################################################
library( 'readr' ) # read and write (large) tsv-files
library( 'ggplot2' ) # plot
library( 'rstanarm' ) # pre-compiled Bayesian models
library( 'bayestestR' ) # ROPE function
library( 'see' ) # plotting
library( 'parameters' )

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
    df$day <- paste0( "day ", df$day )

    df[ df$group == 'injection', 'group' ] <- 'baseline'

    rownames( df ) <- NULL
    
    return( df )
}

###
# get seizure reduction draws
##
get_pain_reduction <- function( bmodel )
{   
    # posterior post-warmup samples [5000]
    draws <- as.data.frame( bmodel )
    
    # get mean pain reduction
    draws$pain_reduction <- draws$groupSMS
    
    values <- round( quantile( draws$pain_reduction, c( 0.025, 0.5, 0.975 ) ), 3 )
    
    out <- data.frame( median = values[ 2 ], CI_low = values[ 1 ], CI_high = values[ 3 ] )
    
    return( out )
}


################################################################################
# CUSTOM FUNCTIONS                                                             #
################################################################################

outdir <- 'out.07.updating.headache'
dir.create( outdir, showWarnings = FALSE )

# get data
head( df <- get_data() )

####################################

all <- NULL
bmodel <- NULL

# first non-zero
start_row <- 14
end_row <- 20
end_row <- nrow( df )

for( i in start_row:end_row )
{
    # set seed
    set.seed( 123 + i )
    
    # get data upto time == i
    data <- df[ 1:i, ]

    day <- df[ i, 'day' ]
    time <- df[ i, 'time' ]
    measurement <- df[ i, 'measurement' ]
    group <- df[ i, 'group' ]
    
    delta_severity <- df[ i, 'delta_severity' ]
    
    if( !is.na( delta_severity ) )
    {
        print( i )
        
        if( i == start_row )
        {
            # fit Bayesian Poisson model (treatment)
            # refresh = 0 turns off iteration info print
            bmodel <- rstanarm::stan_glm( delta_severity ~ group, data = data, iter = 5000, refresh = 0 )
        } else
        {
            # update Bayesian Poisson model
            bmodel <- update( bmodel, formula = delta_severity ~ group, data = data, iter = 5000, refresh = 0 )  
        }
        
        # get PD
        #p <- parameters::parameters( bmodel )
        
        single <- get_pain_reduction( bmodel )
        single$day <- day
        single$time <- time
        single$measurement <- measurement
        single$group <- group
        
        # add to container
        all <- rbind( all, single )
    }
}

# plot update in seizure reduction
p <- 
    ggplot( data = all ) + 
    geom_hline( yintercept = 0, linetype = 2, colour = 'orange' ) + 
    geom_ribbon( aes( x = time, y = median, ymin = CI_low, ymax = CI_high ), fill = 'gray70', alpha = 0.5 ) + 
    geom_line( aes( x = time, y = median ), size = 1.1 ) +
    facet_wrap( ~day, nrow = 1, scale = 'free_x' ) +
    ylab( "SMS-induced pain reduction" ) +
    xlab( 'Time (hour)' ) +
    scale_x_continuous( breaks = c( 10, 12, 14, 16, 18, 20, 22, 24 ) ) +
    scale_y_continuous( breaks = number_ticks( 10 ) ) +
    theme_classic() +
    coord_cartesian( expand = FALSE, xlim = c( 9, 24 ), ylim = c( -0.75, 0.15 ) ) +
    theme( legend.position = 'top', strip.text.x = element_text( face = 'bold' ), 
           strip.text.y = element_text( face = 'bold' ) )

# save
ggsave( plot = p, file = paste0( outdir, '/plot_update_SMS_on_headache.png' ), dpi = 600, height = 3, width = 6 )

# save raw data
write.csv( all, file = paste0( outdir, '/plot_update_SMS_on_headache.csv' ), quote = TRUE )
