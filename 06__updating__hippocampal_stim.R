#!/usr/bin/env Rscript
#
# Wim Otte (w.m.otte@umcutrecht.nl)
# 5 February 2023
#
#
# The Probability of Direction (pd) is an index of effect existence, ranging from 50% to 100%
#
# https://easystats.github.io/blog/posts/bayestestr_pd/
# https://easystats.github.io/bayestestR/articles/probability_of_direction.html
#
# Makowski, D., Ben-Shachar, M. S., & Lüdecke, D. (2019). 
# bayestestR: Describing Effects and their Uncertainty, 
# Existence and Significance within the Bayesian Framework. 
# Journal of Open Source Software, 4(40), 1541. 
# https://doi.org/10.21105/joss.01541
#
# Hipp. stimulation n-of-1 trial updating
#
# rr <- readxl::read_xlsx( 'datasets/Tellez-Zenteno.xlsx' )
# write.csv( rr, file = 'datasets/Tellez-Zenteno.csv', quote = TRUE )
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
    infile <- 'datasets/Tellez-Zenteno.csv'   
    df <- read.csv( infile, stringsAsFactors = TRUE, row.names = 1 )
    df$time <- 1:nrow( df )

    rownames( df ) <- NULL
    
    return( df )
}

###
# Get certainty of drop in outcome (ROPE method)
# https://easystats.github.io/bayestestR/articles/region_of_practical_equivalence.html
##
get_drop_certainty <- function( i, model, prob_drop = 0.3 )
{
    # get baseline coefficient (lambda, i.e., average daily seizures)
    baseline_value <- model$coefficients[ '(Intercept)' ]
    
    # get lower ROPE value [100% reduction - 50%/30% reduction]
    upper_rope_value <- -1 * prob_drop * baseline_value
    lower_rope_value <- -1 * baseline_value
    
    # get perc in ROPE area
    perc_in_rope <- bayestestR::rope( model, ci = 1.0, range = c( lower_rope_value, upper_rope_value ) )
    
    # get table with percentage certainty
    tt <- as.data.frame( perc_in_rope )
    
    perc_certainty_verum <- round( 100 * tt[ tt$Parameter == 'treatmentnicotine', 'ROPE_Percentage' ], 1 )
    
    # output
    output <- data.frame( i, prob_drop, perc_certainty_verum )
    
    return( output )
}

###
# get seizure reduction draws
##
get_seizure_reduction <- function( bmodel )
{   
    # posterior post-warmup samples [5000]
    draws <- as.data.frame( bmodel )
    
    # get mean seizure reduction [i.e., -1.0 = one daily seizure less in treatment group]
    draws$seizure_reduction <- exp( draws$`(Intercept)` ) + exp( draws$treatmenton )
    
    values <- round( quantile( draws$seizure_reduction, c( 0.025, 0.5, 0.975 ) ), 1 )
    
    out <- data.frame( median = values[ 2 ], CI_low = values[ 1 ], CI_high = values[ 3 ] )
    
    return( out )
}


################################################################################
# CUSTOM FUNCTIONS                                                             #
################################################################################


outdir <- 'out.06.updating.hippocampal.stim'
dir.create( outdir, showWarnings = FALSE )


# get data
head( df <- get_data() )


# colors
#custom_colours <- c( "#E69F00", "#56B4E9" )

####################################

all <- NULL
bmodel <- NULL

# first non-zero
start_row <- 5
end_row <- nrow( df )

for( i in start_row:end_row )
{
    # set seed
    set.seed( 123 + i )
    
    # get data upto time == i
    data <- df[ 1:i, ]

    month <- df[ i, 'month' ]
    treatment <- df[ i, 'treatment' ]
    
    print( paste0( 'month: ', month ) )
    
    if( i == start_row )
    {
        # fit Bayesian Poisson model (treatment)
        # refresh = 0 turns off iteration info print
        bmodel <- rstanarm::stan_glm( seizures ~ treatment, data = data, iter = 5000, refresh = 0, family = poisson() )
    } else
    {
        # update Bayesian Poisson model
        bmodel <- update( bmodel, formula = seizures ~ treatment, data = data, iter = 5000, refresh = 0, family = poisson() )  
    }
    
    # get PD
    #p <- parameters::parameters( bmodel )
    
    #single <- data.frame( i = i, pd_treatment = p$pd[ 2 ], pd_to_pvalue = pd_to_p( p$pd[ 2 ] ) )
    single <- get_seizure_reduction( bmodel )
    single$month <- month
    single$treatment <- treatment
    
    # add to container
    all <- rbind( all, single )
}

# make marker to plot period of stimulation (on-off)
all$marker <- 0.01
all[ all$treatment == 'on', 'marker' ] <- 0.3


# plot update in seizure reduction
p <- ggplot( data = all ) + 
    geom_ribbon( aes( x = month, y = median, ymin = CI_low, ymax = CI_high ), fill = 'gray80' ) + 
    geom_line( aes( x = month, y = median ), size = 1.1 ) +
    geom_ribbon( aes( x = month, ymin = 0, ymax = marker ), fill = '#E69F00', alpha = 0.5 ) +
    geom_line( aes( x = month, y = marker ), colour = 'orange', linetype = 2 ) +
    ylab( "Reduced seizures" ) +
    xlab( 'Time (month)' ) +
    scale_x_continuous( breaks = number_ticks( 20 ) ) +
    scale_y_continuous( breaks = number_ticks( 10 ) ) +
    theme_classic() +
    coord_cartesian( expand = FALSE, xlim = c( 0, 48 ), ylim = c( 0, 12.5 ) ) +
    theme( legend.position = 'top', strip.text.x = element_text( face = 'bold' ), 
           strip.text.y = element_text( face = 'bold' ) )

# save
ggsave( plot = p, file = paste0( outdir, '/plot_update_hippocampal_stimulation.png' ), dpi = 600, height = 3, width = 6 )

# save raw data
write.csv( all, file = paste0( outdir, '/plot_update_hippocampal_stimulation.csv' ), quote = TRUE )


#################################

# colors
custom_colours <- c( "#56B4E9", "#E69F00" )

# plot
p2 <- ggplot( df, aes( x = month, y = seizures, fill = treatment, group = treatment ) ) + 
    geom_vline( xintercept = 0, linetype = 2, colour = "gray60" ) +
    geom_point() + # for the fat-black strokes
    geom_point( shape = 21, colour = 'gray40', stroke = 0.2 ) +
    ylab( "Seizures" ) +
    xlab( 'Time (month)' ) +
    scale_x_continuous( breaks = number_ticks( 30 ) ) +
    scale_y_continuous( breaks = number_ticks( 16 ) ) +
    scale_color_manual( values = custom_colours ) +
    scale_fill_manual( values = custom_colours ) +
    coord_cartesian( expand = FALSE, xlim = c( -5, 49 ), ylim = c( -0.3, 16.5 ) ) +
    theme_classic() +
    theme( legend.position = 'top', strip.text.x = element_text( face = 'bold' ), 
           strip.text.y = element_text( face = 'bold' ) )

# save
ggsave( plot = p2, file = paste0( outdir, '/plot_nr_seizures_hippocampal_stimulation.png' ), dpi = 600, height = 3, width = 7 )
