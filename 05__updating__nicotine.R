#!/usr/bin/env Rscript
#
# Wim Otte (w.m.otte@umcutrecht.nl)
# 16 January 2023
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
# Nicotine n-of-1 trial updating
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
    infile <- 'datasets/nicotine.csv'   
    df <- read.csv( infile, stringsAsFactors = TRUE, row.names = 1 )
    df$time <- 1:nrow( df )

    df <- df[ df$type == 'blinded', ]
    rownames( df ) <- NULL
    
    
    df$treatment <- relevel( as.factor( df$treatment ), 'placebo' )
    
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

################################################################################
# CUSTOM FUNCTIONS                                                             #
################################################################################


outdir <- 'out.05.updating.nicotine'
dir.create( outdir, showWarnings = FALSE )


# get data
head( df <- get_data() )


# colors
#custom_colours <- c( "#E69F00", "#56B4E9" )

####################################

all <- NULL
model <- NULL

# first non-zero
start_row <- 22
end_row <- 85

data <- df[1:344,]

for( i in start_row:end_row )
{
    print( i )
    
    # set seed
    set.seed( 123 + i )
    
    # get data upto time == i
    data <- df[ 1:i, ]

    print( paste0( 'nrows: ', nrow( data ) ) )
    
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
    p <- parameters( bmodel )
    
    single <- data.frame( i = i, pd_treatment = p$pd[ 2 ], pd_to_pvalue = pd_to_p( p$pd[ 2 ] ) )
    
    # add to container
    all <- rbind( all, single )
}

library( 'ggplot2' )

# get plot of probability direction over time
p <- ggplot( data = all, aes( x = i, y = pd_treatment * 100 ) ) + 
    geom_line( colour = "gray30", size = 1 ) +
    geom_vline( xintercept = 15, linetype = 2, colour = "gray60" ) +
    geom_vline( xintercept = 43, linetype = 2, colour = "gray60" ) +
    geom_vline( xintercept = 71, linetype = 2, colour = "gray60" ) +
    ylab( "Probability of Direction (%)" ) +
    xlab( 'Time' ) +
    scale_x_continuous( breaks = number_ticks( 16 ) ) +
    scale_y_continuous( breaks = number_ticks( 10 ) ) +
    theme_classic() +
    coord_cartesian( expand = FALSE, xlim = c( 1, 86 ), ylim = c( 99, 100.3 ) ) +
    theme( legend.position = 'top', strip.text.x = element_text( face = 'bold' ), 
           strip.text.y = element_text( face = 'bold' ) )

# remove y-axis information >100% (because not possible)
p_clean <- p + scale_y_continuous(
                breaks = seq( 99, 100, by = 0.1 ),  # Breaks only up to 100%
                labels = function( x ) ifelse( x > 100, "", x ) ) # Remove labels above 100%

# save
ggsave( plot = p_clean, file = paste0( outdir, '/plot_update_nicotine_PD.png' ), dpi = 600, height = 3, width = 6 )





# Compute indices
pd <- p_direction( model )
print( percentage_in_rope <- rope(model, ci = 1) )

# Visualise the pd
plot(pd)
pd

# Visualise the percentage in ROPE
plot(percentage_in_rope)
percentage_in_rope



pp_check( model, type="bars", nsamples = 100)

glm1 <- glm( seizures ~ treatment, data = data, family = poisson )



# TODO

library( "parameters" )
library( "bayestestR" )
library( "rstanarm" )
library( "see" )

# frequentist: get p-value
fmodel <- lm( disp ~ carb, data = mtcars )
parameters( fmodel )

# bayesian: probability of direction [pd] closest to p-value
bmodel <- rstanarm::stan_glm( disp ~ carb, data = mtcars, priors = NULL, prior_intercept = NULL, refresh = 0 )
p <- parameters( bmodel )

# plot
plot( bayestestR::p_direction( bmodel ) )

