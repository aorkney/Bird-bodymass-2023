# This script is a prelude to 'Integration_v_body_mass_11_08_2022.R', 
# Which will then be used to run analyses to generate Figure 1 from Orkney & Hedrick 2023.
# After running 'Integration_v_body_mass_11_08_2022.R', the plotting script
# 'Plot_integration_v_body_mass_14_11_2022.R' will generate the plot

# The purpose of this script is to load 
# centroid size information for 13 skeletal elements across a sample
# of 149 birds from across crown Aves. 
# The landmark constellations from which these data were extracted 
# are available at: https://doi.org/10.18563/journal.m3.125
# (Bjarnason & Benson, 2021)
# Thereafter, the birds' masses will be used to perform
# a phylogenetically gnostic Generalised Least Squares regression, 
# and the residuals will be extracted.
# This will remove allometric scaling in centroid size. 
# The residual centroid sizes may then be investigated by other scripts. 
# The phylogeny of Prum et al., 2015 was used:
# https://doi.org/10.1038/nature15697
# Residual centroid sizes of skeletal elements across birds are known 
# to exhibit a modular organisation of evolutionary covariances. 
# https://doi.org/10.1038/s41559-021-01509-w
# (Orkney, Bjarnason, et al., 2021)

# Load required analytical R packages
library( geomorph )
library( ape )
library( nlme )
# Done

# Set work directory as appropriate
setwd('D:/Documents/Alex_birds/Ideas_2022/SICBY')
# The user will need to customise

# Load centroid sizes for study birds. 
load('Csize.22.10.2022.RData')
# In our study, we defined this object as an array called 'GPA.Csize'
# Be aware the names do not yet match the names in the phylogeny we will use.

# Load phylogenetic tree of birds
load('tree.22.10.2022.RData')
# (Pruned from Prum et al., 2015)
# The object is a phylogeny called 'pruned.tree'.

# Load the tree names required to match birds to the closest genus on the tree.
load('tree.names.22.10.2022.RData')
# The object is a vector called 'tree_names'.

# Load bird masses
load('masses.22.10.2022.Rdata'))
# The object is a named numeric vector called 'masses'.
# The names will not yet match those on the phylogeny.

# Let's now ensure names match the phylogeny. 
names(masses) <- tree_names
for( i in 1:length(GPA.Csize) ){
	names(GPA.Csize[[i]]) <- tree_names 
}
# Done


# We may wish to investigate a specific taxonomic subset of 
# the available birds. 

# Set work directory as appropriate
setwd("D:/Documents/Alex_birds/")
# The user will need to customise

# Load a vector of birds belonging to clade Telluraves
Telluraves<-read.csv('Telluraves.csv')
# This is a pre-prepared list of birds belonging to Telluraves.

# Define a vector of birds belonging to clade Apodiformes
Apodiformes<- c('Archilochus_colubris',
'Topaza_pella',
'Chaetura_brachyura',
'Streptoprocne_zonaris',
'Hemiprocne_comata')
# This is a pre-prepared list of birds belonging to Apodiformes.

# Set work directory as appropriate
setwd('D:/Documents/Alex_birds/Ideas_2022/SICBY')
# The user will need to customise

# Define a sundry function to interract with data objects such as GPA.coords
get.item <- function( X , item ) { X[[ item ]] }
# This is a simple indexing function

# Define a function to compute the residual skeletal element centroid cizes, after the affects of allometry are removed
get.residual.Csize <- function( array, masses, phylogeny, taxa ){ # The function takes a shape array, mass vector, phylogeny and list of taxa
	species<-taxa
	allometry.Csize <- list() # Dumby variable to receive allometric models
	newphy <- drop.tip( phylogeny , phylogeny$tip.label[ !phylogeny$tip.label %in% taxa ] ) # Phylogeny pruned to desired taxa
	for(i in 1:length(array) ){ # For each skeletal element 
		df <- data.frame( mass= log10(masses[taxa]), Csize = log10( array[[ i ]][taxa] ), species = taxa  ) # Define a data frame
		allometry.Csize[[ i ]] <- gls( Csize ~ mass, correlation = corPagel( 0, phy=newphy, form= ~ species ), data=df ) # Compute an allometric model
	}
	names(allometry.Csize) <- names(array) # Ensure names of allometric models match the skeletal elements they describe 
	residual_Csize <- lapply( lapply( allometry.Csize , get.item , item = "residuals" ) , c ) # Extract model residuals
	return(residual_Csize) # Return residuals
} # Conclude function
# This function uses phylogenetic Generalised Least Squares to compute the allometrically adjusted centroid sizes of skeletal elements.
# The user must provide the array of original bird skeletal element centroid sizes, bodymasses, phylogeny and taxa of interest. 

# Compute residual centroid sizes by applying function:
residual.Csize <- get.residual.Csize( array = GPA.Csize, masses = masses, phylogeny = pruned.tree, taxa=names(masses) )
# Which taxa do we wish to investigate? 
# In this example, we will investigate all birds. 

# An array of residual centroid sizes is now available for investigation. 

# Script concludes
