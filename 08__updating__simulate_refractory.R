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
    infile <- 'out.02.simulate.refractory/data__weekly-seizures-baseline-35__placebo-verum-risk-0.1--0.6.tsv.gz' 
    df <- readr::read_tsv( infile, show_col_types = FALSE )
    
    return( df )
}

###
# get seizure reduction draws
##
get_seizure_reduction <- function( bmodel )
{   
    # posterior post-warmup samples [5000]
    draws <- as.data.frame( bmodel )
    
    # get mean seizure reduction
    draws$seizure_reduction <- exp( draws$`(Intercept)` + draws$groupverum ) - exp( draws$`(Intercept)` + draws$groupplacebo )
    
    values <- round( quantile( draws$seizure_reduction, c( 0.025, 0.5, 0.975 ) ), 2 )
    
    out <- data.frame( median = values[ 2 ], CI_low = values[ 1 ], CI_high = values[ 3 ] )
    
    return( out )
}

################################################################################
# CUSTOM FUNCTIONS                                                             #
################################################################################

outdir <- 'out.08.updating.simulate.refractory'
dir.create( outdir, showWarnings = FALSE )

# get data
head( df <- get_data() )

####################################

all <- NULL
bmodel <- NULL

# first row with data from both groups
start_row <- 85
end_row <- nrow( df )

for( i in start_row:end_row )
{
    # set seed
    set.seed( 123 + i )
    
    # get data upto time == i
    data <- df[ 1:i, ]

    day <- df[ i, 'day' ]
    cycle <- df[ i, 'cycle' ]
    
    print( i )
        
    if( i == start_row )
    {
        # fit Bayesian Poisson model (treatment)
        # refresh = 0 turns off iteration info print
        bmodel <- rstanarm::stan_glm( y ~ group, data = data, iter = 5000, refresh = 0, family = poisson() )
    } else
    {
        # update Bayesian Poisson model
        bmodel <- update( bmodel, formula = y ~ group, data = data, iter = 5000, refresh = 0, family = poisson() )  
    }
    
    # get PD
    #p <- parameters::parameters( bmodel )
    
    single <- get_seizure_reduction( bmodel )
    single$day <- as.character( day )
    single$cycle <- as.character( cycle )

    
    # add to container
    all <- rbind( all, single )
    
}

all$day <- NULL
all$cycle <- NULL
all$time <- start_row - 1 + 1:nrow( all )

# reduction
baseline <- 5
yintercept <- -1 * baseline * 0.3 # -1.5

# plot update in seizure reduction
p <- ggplot( data = all ) + 
    # 30% reduction lines (horizontal)
    geom_hline( yintercept = yintercept, linetype = 2, colour = 'orange' ) + 
    # cycle boundaries (vertical)
    geom_vline( xintercept = ( 7 * 8 ) + 0.0, linetype = 2, colour = 'gray50' ) + 
    geom_vline( xintercept = ( 7 * 16 ) + 0.0, linetype = 2, colour = 'gray50' ) + 
    geom_vline( xintercept = ( 7 * 24 ) + 0.0, linetype = 2, colour = 'gray50' ) + 
    geom_vline( xintercept = ( 7 * 32 ) + 0.0, linetype = 2, colour = 'gray50' ) + 
    geom_vline( xintercept = ( 7 * 40 ) + 0.0, linetype = 2, colour = 'gray50' ) + 
    geom_ribbon( aes( x = time, y = median, ymin = CI_low, ymax = CI_high ), fill = 'gray70', alpha = 0.5 ) + 
    geom_line( aes( x = time, y = median ), size = 1.1 ) +
    ylab( "Seizure reduction (verum-induced)" ) +
    xlab( 'Time (day)' ) +
    scale_x_continuous( breaks = seq( from = 0, to = 336, by = 14 ) ) +
    scale_y_continuous( breaks = number_ticks( 10 ) ) +
    theme_classic() +
    coord_cartesian( expand = FALSE, xlim = c( 49, 338 ), ylim = c( -4.6, -0.4 ) ) +
    theme( legend.position = 'top', strip.text.x = element_text( face = 'bold' ), 
           strip.text.y = element_text( face = 'bold' ) )

p2 <- p + 
    annotate( "text", label = "cycle 1", x = 70, y = -0.5, size = 3, colour = "#56B4E9", fontface = "bold" ) +
    annotate( "text", label = "cycle 2", x = 126, y = -0.5, size = 3, colour = "#56B4E9", fontface = "bold" ) +
    annotate( "text", label = "cycle 3", x = 182, y = -0.5, size = 3, colour = "#56B4E9", fontface = "bold" ) +
    annotate( "text", label = "cycle 4", x = 238, y = -0.5, size = 3, colour = "#56B4E9", fontface = "bold" ) +
    annotate( "text", label = "cycle 5", x = 294, y = -0.5, size = 3, colour = "#56B4E9", fontface = "bold" )

# save
ggsave( plot = p2, file = paste0( outdir, '/plot_update_simulate_refractory.png' ), dpi = 600, height = 3, width = 6 )

# save raw data
write.csv( all, file = paste0( outdir, '/plot_update_simulate_refractory.csv' ), quote = TRUE )
