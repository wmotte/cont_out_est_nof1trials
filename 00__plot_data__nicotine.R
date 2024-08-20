#!/usr/bin/env Rscript
#
# Visualize n-of-1 trial data - epilepsy nicotine patches (2003 study).
#
# Wim Otte (w.m.otte@umcutrecht.nl)
#
################################################################################
library( 'ggplot2' )        # for plotting

###
# Get data
##
get_data <- function()
{
    # read first tab sheet
    infile <- 'datasets/nicotine.csv'   
    df <- read.csv( infile, stringsAsFactors = TRUE, row.names = 1 )
    df$time <- 1:nrow( df )
    #df$day <- paste0( "day ", df$day )
    
    return( df )
}

### 
# Helper function plot.
##
number_ticks <- function( n ) { function( limits ) pretty( limits, n ) }

################################################################################
# END FUNCTIONS
################################################################################

# output dir
outdir <- 'out.00.plot.nicotine'
dir.create( outdir, showWarnings = FALSE )

# get data
head( df <- get_data() )

# colors
custom_colours <- c( "#E69F00", "#56B4E9" )

####################################

# get start and end of n-of-1 trial
startd <- as.numeric( rownames( df[ df$type == 'blinded', ] )[ 1 ] )
endd <- nrow( df )

# plot
p <- ggplot( df, aes( x = as.Date( day ), y = seizures, fill = treatment, group = treatment ) ) + 
        geom_point() + # for the fat-black strokes
        geom_point( shape = 21, colour = 'gray40', stroke = 0.2 ) +
        geom_vline( xintercept = as.numeric( as.Date( df$day[ c( startd, endd ) ] ) ), linetype = 2, colour = "gray60" ) +
        ylab( "Seizures (per day)" ) +
        xlab( 'Time' ) +
        scale_y_continuous( breaks = number_ticks( 10 ) ) +
        scale_color_manual( values = custom_colours ) +
        scale_fill_manual( values = custom_colours ) +
        theme_classic() +
        theme( legend.position = 'top', strip.text.x = element_text( face = 'bold' ), 
                                    strip.text.y = element_text( face = 'bold' ) )


# save
ggsave( plot = p, file = paste0( outdir, '/plot_nicotine.png' ), dpi = 600, height = 4, width = 8 )


####################################################################

# n-of-1 only (so remove 'open' pre-phase)
df_blinded <- df[ df$type == 'blinded', ]
df_blinded$time <- 1:nrow( df_blinded )

# plot
p_blind <- ggplot( data = df_blinded, aes( x = time, y = seizures, fill = treatment, group = treatment ) ) +
    geom_point() + # for the fat-black strokes
    geom_point( shape = 21, colour = 'gray40', stroke = 0.2 ) +
    
    #geom_point( size = 2 ) +
    ylab( "Seizures (per day)" ) +
    xlab( 'Time' ) +
    scale_x_continuous( breaks = number_ticks( 16 ) ) +
    scale_y_continuous( breaks = number_ticks( 7 ) ) +
    scale_color_manual( values = custom_colours ) +
    scale_fill_manual( values = custom_colours ) +
    theme_classic() +
    coord_cartesian( expand = FALSE, xlim = c( 1, 86 ), ylim = c( -0.14, 6 ) ) +
    theme( legend.position = 'top', strip.text.x = element_text( face = 'bold' ), 
           strip.text.y = element_text( face = 'bold' ) )

# save
ggsave( plot = p_blind, file = paste0( outdir, '/plot_nicotine_blinded.png' ), dpi = 600, height = 2, width = 6 )

p2 <- p_blind +
    geom_vline( xintercept = 0.5 + 7 * 2, linetype = 2, colour = 'gray50' ) + 
    geom_vline( xintercept = 0.5 + 7 * 4, linetype = 2, colour = 'gray50' ) +
    geom_vline( xintercept = 0.5 + 7 * 6, linetype = 2, colour = 'gray50' ) + 
    geom_vline( xintercept = 0.5 + 7 * 8, linetype = 2, colour = 'gray50' ) +
    geom_vline( xintercept = 0.5 + 7 * 10, linetype = 2, colour = 'gray50' ) +
    geom_vline( xintercept = 0.5 + 7 * 12, linetype = 2, colour = 'gray50' )

# save
ggsave( plot = p2, file = paste0( outdir, '/plot_nicotine_blinded__with_lines.png' ), dpi = 600, height = 2, width = 6 )

