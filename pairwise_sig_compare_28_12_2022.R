# This script produces Figure 2 in Orkney & Hedrick, 2023

# The purpose of this script is to load 
# centroid size information for 13 skeletal elements across a sample
# of 149 birds from across crown Aves...
# The landmark constellations from which these data were extracted 
# are available at: https://doi.org/10.18563/journal.m3.125
# (Bjarnason & Benson, 2021)
# Thereafter, the birds' masses will be used to perform
# a phylogenetically gnostic Generalised Least Squares regression, 
# and the residuals will be extracted...
# This will remove allometric scaling in centroid size...
# The phylogeny of Prum et al., 2015 was used:
# https://doi.org/10.1038/nature15697
# Residual centroid sizes of skeletal elements across birds are known 
# to exhibit a modular organisation of evolutionary covariances...
# https://doi.org/10.1038/s41559-021-01509-w
# (Orkney, Bjarnason, et al., 2021)

# Having prepared the dataset, 
# phylogenetic two block Partial Least Squares analyses will be performed,
# in order to assess the degree of integration between different pairwise 
# combinations of bones (carpometacarpus-humerus, sternum-scapula,
# sternum-carpometacarpus); These integration statistics will be computed
# in different cohorts of 30 birds of varying mean body mass...
# The significance of difference between the integration statistics will 
# then be computed between the cohorts, using the compare.pls function
# from the R Geomorph package... 
# This will clarify whether integration within the wing, within the trunk
# and between the wing and trunk, varies significantly as a function of body mass,
# across birds...
# The resultant plots will illustrate:

# 1) Whether integration is more or less intense in the heavier cohort of birds.
# 2) Only fields where p < 0.05 will be displayed.

# Load required analytical R packages
library( geomorph )
library( ape )
library( nlme )
# If you do not have a package installed, run the command 'install.packages('name')'
# and proceed as guided

# Set work directory as appropriate
setwd('D:/Documents/Alex_birds/Ideas_2022/SICBY')
# The user will need to customise this

# Load required data
load('Csize.22.10.2022.RData')
# Load centroid sizes for study birds
 
# In this study, this object was defined as an array called 'GPA.Csize'
# Be aware the names do not yet match the names in the phylogeny will be used

# Load the phylogeny
load('tree.22.10.2022.RData')
# (Pruned from Prum et al., 2015)
# The object is a phylogeny called 'pruned.tree'

# Load the tree names required to match birds to the closest genus on the tree
load('tree.names.22.10.2022.RData')
# The object is a vector called 'tree_names'

# Load bird masses
load('masses.22.10.2022.RData')
# The object is a named numeric vector called 'masses'
# The names will not yet match those on the phylogeny

# Ensure names match the phylogeny
names(masses) <- tree_names
for( i in 1:length(GPA.Csize) ){
	names(GPA.Csize[[i]]) <- tree_names 
}
# Done



# The user ,au wish to investigate a specific taxonomic subset of 
# the available birds

# Set work directory as appropriate
setwd("D:/Documents/Alex_birds/")
# The user will need to customise this path

# Read a list of birds belonging to clade 'Telluraves'
Telluraves<-read.csv('Telluraves.csv')
# This is a pre-prepared vector of birds belonging to Telluraves
# This will not actually be used in analysis, but may be helpful
# for other researchers who wish to perform custom analyses

# Define a vector of birds belonging to clade Apodiformes
Apodiformes<- c('Archilochus_colubris',
'Topaza_pella',
'Chaetura_brachyura',
'Streptoprocne_zonaris',
'Hemiprocne_comata')
# This is a pre-prepared list of birds belonging to clade Apodiformes

# Set work directory as appropriate
setwd('D:/Documents/Alex_birds/Ideas_2022/SICBY')
# The user will need to customise this path

# Define a sundry function for extracting data from pre-defined data objects
get.item <- function( X , item ) { X[[ item ]] }
# This is a simple indexing function

# Define a function to compute the residual centroid sizes of landmark constellations in a phylogenetic context, accounting for allometric scaling
get.residual.Csize <- function( array, masses, phylogeny, taxa, bones ){ # The function takes a landmark array, mass vector, phylogeny, taxa and bone names as arguments
	allometry.Csize <- list() # Define a dumby list to receive data 
	newphy <- drop.tip( phylogeny , phylogeny$tip.label[ !phylogeny$tip.label %in% taxa ] ) # Prune tree to taxa
	for(i in 1:length(array[bones]) ){ # For each bone
		df <- data.frame( mass= log10(masses[taxa]), Csize = log10( array[[ bones[i] ]][taxa] ), species = taxa ) # Organise data into data frame
		allometry.Csize[[ i ]] <- gls( Csize ~ mass, correlation = corPagel( 0, phy=newphy, form= ~ species  ), data=df ) # Compute allometric model in phylogenetic context; centroid size depending on mass
	}
	names(allometry.Csize) <- names(array[bones]) # Make sure computed models conform to the specified bones being requested 
	residual_Csize <- lapply( lapply( allometry.Csize , get.item , item = "residuals" ) , c ) # Get the residuals of the allometrid model fit
	return(residual_Csize) # Return the desired data
} # Function concludes
# This function uses phylogenetic Generalised Least Squares to compute the allometrically adjusted centroid sizes of skeletal elements
# The user must provide the array of original bird skeletal element centroid sizes, bodymasses, phylogeny and taxa of interest

# Define the number of birds to be included in any cohort (size of 'taxa')
width <- 30
# Done

# The following is a function that will compute significance of difference of pairwise integration
# between two bones among various combinations of lighter and heavier cohorts of birds, 
# within a specific taxonomic subset
# Comments within the function will generally follow lines, unless they describe multiple lines' actions, in which case they will be placed on the line immediately underneath the associated block
pair.sig <- function( mass, phylogeny, bone1, bone2){ This function takes a mass vector, phylogeny and two bones as input 
	lights <- c(1:width) # Begin by defining the lightest cohort of birds within the birds of interest 
	bin.mass<-matrix(NA,(length( mass )-2*(length(lights))),(length( mass )-2*(length(lights))))  
	p.dff<-matrix(NA,(length( mass )-2*(length(lights))),(length( mass )-2*(length(lights))))
	heavy.int<-matrix(NA,(length( mass )-2*(length(lights))),(length( mass )-2*(length(lights))))
	light.int<-matrix(NA,(length( mass )-2*(length(lights))),(length( mass )-2*(length(lights))))
	# Define matrices to receive data resulting from analyses 

	for(k in 1: dim(bin.mass)[1]){ # For all possible combinations of cohorts of birds. 
		lighter<-which( mass[order(mass)] < mass[order(mass)][min(lights)] ) # Find all birds less massive than the 'light' cohort
		heavies<-c(1:length(mass))[-c(lighter,lights)] # Define a cohort of birds heavier than the 'light' cohort
		for(i in 1:(length(heavies)-width) ){ # For all possible cohorts of birds within the heavier group
			start <- heavies[i] 
			end <- heavies[i+(width-1)]
			# Define the vector of heavier birds to be considered
			bin.mass[k,i+length(lighter)] <- mean(mass[order(mass)][start:end]) # Compute the mean mass of the heavy cohort
			taxa.heavies <- names(mass[order(mass)][start:end]) # Name the birds in the heavy cohort
			newphy.heavies <- drop.tip( phylogeny, phylogeny$tip.label[ !phylogeny$tip.label  %in% taxa.heavies] ) # Prune the phylogeny as required
			residual.Csize <-tryCatch( expr= { get.residual.Csize( array = GPA.Csize, masses = masses, phylogeny = pruned.tree, taxa=taxa.heavies, bones<-c(bone1,bone2) )},
			error=function(e) e
			 )
			if(inherits(residual.Csize, "error")) next
			# Compute the residual centroid sizes for the heavy cohort

			int.heavies <- phylo.integration(as.matrix(residual.Csize[[bone1]][taxa.heavies]),as.matrix(residual.Csize[[bone2]][taxa.heavies]),phy=newphy.heavies,print.progress=F)# Compute integration within the heavier cohort 
			taxa.lights <- names(mass[order(mass)][lights ]) # Name the birds in the light cohort
			newphy.lights  <- drop.tip( phylogeny, phylogeny$tip.label[ !phylogeny$tip.label  %in% taxa.lights ] ) # Prune the phylogeny as required 
			residual.Csize <-tryCatch( expr= { get.residual.Csize( array = GPA.Csize, masses = masses, phylogeny = pruned.tree, taxa=taxa.lights, bones<-c(bone1,bone2) )},
			error=function(e) e
			 )
			if(inherits(residual.Csize, "error")) next
			# Compute the residual centroid sizes for the light cohort

			int.lights  <- phylo.integration(as.matrix(residual.Csize[[bone1]][taxa.lights ]),as.matrix(residual.Csize[[bone2]][taxa.lights ]),phy=newphy.lights,print.progress=F )# Compute integration within the light cohort 
			p.dff[k,i+length(lighter)] <- compare.pls(int.lights,int.heavies)$pairwise.P[2]# Compute the significance of difference between the two integration analyses
			heavy.int[k,i+length(lighter)] <- int.heavies$Z
			light.int[k,i+length(lighter)] <- int.lights$Z
			# Assign values to the output matrices
		}
	print(round(100*(k/dim(bin.mass)[1]))) # Update on progress; this loop will take a long time to run
	lights<-lights+1 # Shift the light cohort to exclude the lightest bird within, and include the lightest bird of the remaining heavier birds
	}

	lights<-c(1:width)
	light.masses<-list()
	for(k in 1: dim(bin.mass)[1]){
		light.masses[[k]]<-mean(mass[order(mass)][lights])
		lights<-lights+1
	}
	# Compute the masses of the light cohorts

	colnames(p.dff)<-bin.mass[1,]
	rownames(p.dff)<-unlist(light.masses)
	colnames(heavy.int)<-bin.mass[1,]
	rownames(heavy.int)<-unlist(light.masses)
	colnames(light.int)<-bin.mass[1,]
	rownames(light.int)<-unlist(light.masses)
	# Tidy the data output

	output <- list(p.dff,bin.mass,heavy.int,light.int)
	return(output)
	# Produce output
} # Loop concludes

# The defined functions can now be applied to perform analyses
wing.all <- pair.sig(masses, pruned.tree, 'carpometacarpus', 'humerus') Figure 2 subplot a
trunk.all <- pair.sig(masses,  pruned.tree, 'scapula', 'sternum') Figure 2 subplot b 
cross.all <- pair.sig(masses,  pruned.tree, 'carpometacarpus', 'sternum') Figure 2 subplot c
# Analyses have been performed asking how the significance of integration between the carpometacarpus and humerus, 
# scapula and sternum, and carpometacarpus and sternum, evolve with body mass across birds 

# Check that the working directory is a satisfactory location to save the output; it took a long time to produce it
getwd()
# The user must verify they are happy
# save.image(file='pairwise_sig_compare_28_12_2022.RData')
# The analysis takes a long time to run, so the user may wish to save the output; uncomment the above line to do so
# load(file='pairwise_sig_compare_28_12_2022.RData')
# Uncomment the above line to load previously saved output

# The following packages are necessary to reshape and interpolate data
library(reshape2)
library(MBA)
# Packages loaded

# The following lines will define colour bar limits for plots:
zlims<-c( min( c(
(wing.all[[3]]-wing.all[[4]]),
(trunk.all[[3]]-trunk.all[[4]]),
(cross.all[[3]]-cross.all[[4]]) 
),na.rm=T ), 
max( c(
(wing.all[[3]]-wing.all[[4]]),
(trunk.all[[3]]-trunk.all[[4]]),
(cross.all[[3]]-cross.all[[4]]) 
),na.rm=T ) )
zlims[1]<-zlims[1]-0.1
zlims[2]<-zlims[2]+0.1
# Done

# The following lines will compute the difference in integration between light and heavy cohorts fo birds 
# and re-organise into a format that can be supplied to plotting functions 
difference <- wing.all[[3]]-wing.all[[4]]
difference <- melt(difference)
difference <- difference[complete.cases(difference),]
df<-melt(wing.all[[1]])
df <- df[complete.cases(df),]
# Done 

# Organise the data into a format that can be provided to interpolation functions
orig.df<-difference
colnames(orig.df)<-c('x','y','z')
orig.df$x<-log10(orig.df$x)
orig.df$y<-log10(orig.df$y)
# Done 

# Possible definitions for limits for subplot Figure 2 a
wing.lims.x<-c(min(c(orig.df$y)),max(c(orig.df$y)))
wing.lims.y<-c(min(c(orig.df$x)),max(c(orig.df$x)))
# Defined

# Interpolate the difference of integration statistics onto a dense grid
grid_mba<-mba.surf(orig.df, no.X = 100, no.Y = 100, extend = T)
dimnames(grid_mba$xyz.est$z) <- list(grid_mba$xyz.est$x, grid_mba$xyz.est$y)
grid_mba <- melt(grid_mba$xyz.est$z, varnames = c('x', 'y'), value.name = 'z')
grid_mba[grid_mba$x>grid_mba$y,3]<-NA
# Done 

# Define plotting limits
wing.lims.z<-c(min(c(grid_mba$z),na.rm=T),max(c(grid_mba$z),na.rm=T))
# 

# The following lines will repeat the above protocol for associated p-values:
orig.df<-df
orig.df<-log10(orig.df)
colnames(orig.df)<-c('x','y','z')
pgrid_mba<-mba.surf(orig.df, no.X = 100, no.Y = 100, extend = T)
dimnames(pgrid_mba$xyz.est$z) <- list(pgrid_mba$xyz.est$x, pgrid_mba$xyz.est$y)
pgrid_mba <- melt(pgrid_mba$xyz.est$z, varnames = c('x', 'y'), value.name = 'z')
pgrid_mba[pgrid_mba$x>pgrid_mba$y,3]<-NA
# Interpolate the significance statistics onto a dense grid

# Restrict the interpolated difference of integration statistics to only those with associated
# significant P-values
grid_mba <- grid_mba[pgrid_mba$z < log10(0.05),]
# Done

# Load package required for plotting
library(ggplot2)
# Done

# The code for the first subplot will be extensively commented, and comments will then largely be desisted for following subplots that repeat existing patterns of code
mp.wing <- ggplot(grid_mba, aes(y, x, fill = z, z=z)) + Define Figure 2 subplot a
lims(x=wing.lims.x,y=wing.lims.y)+ # Define axial limits
geom_raster( )+ # Create a raster layer
labs(y='mass light cohort', x='mass heavy cohort',fill=expression(paste(Delta,italic('Z') )))+ # Label the axes approrpiately
scale_fill_gradient2(high='red',low='blue',mid='white',na.value='white',midpoint=0,limits=zlims)+ # Fill the raster with a colour gradient representing difference in integration 
# between lighter and heavier cohorts of birds, within the pairwise combination of bones of intereset (in this case it was carpometacarpus-humerus) 
# Sundry thematic plotting instructions follow:
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=10,colour='black'),
      axis.text.y=element_text(size=10,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.pos='bottom',
	legend.key.width=unit(1.4,'cm'),
	legend.key.height=unit(0.8,'cm'),
	legend.title=element_text(size=15),
	legend.text=element_text(size=12),
	)
# Plotting instructions concluded

# Extract a legend
mp.legend<-cowplot::get_legend(mp.wing)
# (package 'cowplot' must be installed)

# Truncate the plotting limits
wing.lims.y[2]<-2
# Done

mp.wing<- ggplot(grid_mba, aes(y, x, fill = z, z=z)) + # Produce a version of Figure 2 subplot a without a legend
lims(x=wing.lims.x,y=wing.lims.y)+
geom_raster( )+
labs(y=expression(paste(log[10],'(',italic(M[Light]),') [g]')), x=expression(paste(log[10],'(',italic(M[Heavy]),') [g]')),fill=expression(paste(Delta,italic('Z') )))+
scale_fill_gradient2(high='red',low='blue',mid='white',na.value='white',midpoint=0,limits=zlims)+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=12,colour='black'),
      axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
	legend.position='none',
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
# Plotting concludes


# The above pattern of data organisation, interpolation, curation, and plotting is now computed for the trunk region (sternum-scapula combination);
# Comments will be minimal
difference <- trunk.all[[3]]-trunk.all[[4]]
difference <- melt(difference)
difference <- difference[complete.cases(difference),]
df<-melt(trunk.all[[1]])
df <- df[complete.cases(df),]
# Differences of integration between lighter and heavier cohorts of birds computed

orig.df<-difference
colnames(orig.df)<-c('x','y','z')
orig.df$x<-log10(orig.df$x)
orig.df$y<-log10(orig.df$y)
# Data organised so as to foster interpolation

trunk.lims.x<-c(min(c(orig.df$y)),max(c(orig.df$y)))
trunk.lims.y<-c(min(c(orig.df$x)),max(c(orig.df$x)))
trunk.lims.z<-c(min(c(grid_mba$z),na.rm=T),max(c(grid_mba$z),na.rm=T))
# Define plotting limits

grid_mba<-mba.surf(orig.df, no.X = 100, no.Y = 100, extend = T)
dimnames(grid_mba$xyz.est$z) <- list(grid_mba$xyz.est$x, grid_mba$xyz.est$y)
grid_mba <- melt(grid_mba$xyz.est$z, varnames = c('x', 'y'), value.name = 'z')
grid_mba[grid_mba$x>grid_mba$y,3]<-NA
# Interpolate the difference of integration statistics onto a dense grid 

orig.df<-df
orig.df<-log10(orig.df)
colnames(orig.df)<-c('x','y','z')
pgrid_mba<-mba.surf(orig.df, no.X = 100, no.Y = 100, extend = T)
dimnames(pgrid_mba$xyz.est$z) <- list(pgrid_mba$xyz.est$x, pgrid_mba$xyz.est$y)
pgrid_mba <- melt(pgrid_mba$xyz.est$z, varnames = c('x', 'y'), value.name = 'z')
pgrid_mba[pgrid_mba$x>pgrid_mba$y,3]<-NA
# Interpolate the associated significance statistics onto a dense grid 

grid_mba <- grid_mba[pgrid_mba$z < log10(0.05),]
# Restrict interpolation to only significant cells


mp.trunk<- ggplot(grid_mba, aes(y, x, fill = z, z=z)) + # Produce Figure 2 subplot b
lims(x=wing.lims.x,y=wing.lims.y)+
geom_raster( )+
labs(y=expression(paste(log[10],'(',italic(M[Light]),') [g]')), x=expression(paste(log[10],'(',italic(M[Heavy]),') [g]')),fill=expression(paste(Delta,italic('Z') )))+
scale_fill_gradient2(high='red',low='blue',mid='white',na.value='white',midpoint=0,limits=zlims)+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=12,colour='black'),
      axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
	legend.position='none',
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
# Plot the trunk analysis output; plotting concludes


# The above pattern of data organisation, interpolation, curation, and plotting is now computed between the trunk and wing region (sternum-carpometacarpus combination);
# Comments will be minimal
difference <- cross.all[[3]]-cross.all[[4]]
difference <- melt(difference)
difference <- difference[complete.cases(difference),]
df<-melt(cross.all[[1]])
df <- df[complete.cases(df),]
# Differences of integration between lighter and heavier cohorts of birds computed

orig.df<-difference
colnames(orig.df)<-c('x','y','z')
orig.df$x<-log10(orig.df$x)
orig.df$y<-log10(orig.df$y)
# Data organised so as to foster interpolation

cross.lims.x<-c(min(c(orig.df$y)),max(c(orig.df$y)))
cross.lims.y<-c(min(c(orig.df$x)),max(c(orig.df$x)))
cross.lims.z<-c(min(c(grid_mba$z),na.rm=T),max(c(grid_mba$z),na.rm=T))
# Define plotting limits

grid_mba<-mba.surf(orig.df, no.X = 100, no.Y = 100, extend = T)
dimnames(grid_mba$xyz.est$z) <- list(grid_mba$xyz.est$x, grid_mba$xyz.est$y)
grid_mba <- melt(grid_mba$xyz.est$z, varnames = c('x', 'y'), value.name = 'z')
grid_mba[grid_mba$x>grid_mba$y,3]<-NA
# Interpolate the difference of integration statistics onto a dense grid 

orig.df<-df
orig.df<-log10(orig.df)
colnames(orig.df)<-c('x','y','z')
pgrid_mba<-mba.surf(orig.df, no.X = 100, no.Y = 100, extend = T)
dimnames(pgrid_mba$xyz.est$z) <- list(pgrid_mba$xyz.est$x, pgrid_mba$xyz.est$y)
pgrid_mba <- melt(pgrid_mba$xyz.est$z, varnames = c('x', 'y'), value.name = 'z')
pgrid_mba[pgrid_mba$x>pgrid_mba$y,3]<-NA
# Interpolate the associated significance statistics onto a dense grid 

grid_mba <- grid_mba[pgrid_mba$z < log10(0.05),]
# Restrict interpolation to only those cells associated with significant differences


mp.cross<- ggplot(grid_mba, aes(y, x, fill = z, z=z)) + # Produce Figure 2 subplot c
lims(x=wing.lims.x,y=wing.lims.y)+
geom_raster( )+
labs(y=expression(paste(log[10],'(',italic(M[Light]),') [g]')), x=expression(paste(log[10],'(',italic(M[Heavy]),') [g]')),fill=expression(paste(Delta,italic('Z') )))+
scale_fill_gradient2(high='red',low='blue',mid='white',na.value='white',midpoint=0,limits=zlims)+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=12,colour='black'),
      axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
	legend.position='none',
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
# Plot the cross-wing-trunk analysis output; plotting concludes

# The following lines produce thematic axis decorations that will be added below the plot; they
# are not integral to analysis and will not be extensively commented:

df<-as.data.frame(cbind(c(-.2,-.2,.2,.2),c(-.2,.2,-.2,.2)))
colnames(df)<-c('x','y')

arrow_plot<-ggplot(data=df, aes(x, y)) +
geom_point(col='white')+
lims(x=c(-2,2),y=c(0.5,0.8))+
  geom_segment(
    x = 1, y = 0.75,
    xend = 1.8, yend = 0.75,
    lineend = "round",    linejoin = "round",size = 2, 
    arrow = arrow(length = unit(0.3, "cm")),
    colour = 'black'   )+
  geom_segment(
    x = -1, y = 0.75,
    xend = -1.8, yend = 0.75,
    lineend = "round",    linejoin = "round",size = 2, 
    arrow = arrow(length = unit(0.3, "cm")),
    colour = 'black'   )+
geom_text(x=1.5,y=0.6,aes(label='Higher integration \n in heavy birds'),size=3)+
geom_text(x=-1.5,y=0.6,aes(label='Higher integration \n in light birds'),size=3)+
  theme(axis.line.x.bottom=element_blank(),
	axis.line.y.left=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
	axis.title.x.bottom=element_blank(),
	axis.title.y.left=element_blank(),
      axis.ticks=element_blank(),
	legend.position='none',
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
# A subsidiary plot is made to label the colour bar

# Load a package to arrange multiple subplots into a single combined plot
library(ggpubr)
# Done


dev.new(width=9,height=18,unit='cm')
ggarrange('',mp.wing,mp.trunk,mp.cross,'',mp.legend,arrow_plot,ncol=1,nrow=7,
labels=c('','a','b','c','','',''),label.x=0,label.y=1.160,font.label=list(size=25),heights=c(0.1,1,1,1,0.1,0.275,0.3))
# Assemble the subplots

# Check working directory
getwd() 
# The user should establish that they are satisfied this is a suitable location to save the output
# ggsave(filename='pairwise_plot_28_12_2022.pdf')
# ggsave(filename='pairwise_plot_28_12_2022.png',dpi=300)
# ggsave(filename='pairwise_plot_28_12_2022.jpg',dpi=300)
# Uncomment the above lines to save .PDF and high quality .PNG and .JPG output


# Should the user desire an alternative representation:
dev.new(width=9,height=18,unit='cm')
ggarrange('',mp.wing,mp.trunk,mp.cross,'','','',ncol=1,nrow=7,
labels=c('','a','b','c','','',''),label.x=0,label.y=1.160,font.label=list(size=25),heights=c(0.1,1,1,1,0.1,0.275,0.3))
# Assemble plots
ggsave(filename='pairwise_plot_17_01_2022_subplot_labels.pdf')
ggsave(filename='pairwise_plot_17_01_2022_subplot_labels.jpeg',dpi=300)
# Save output as high quality .PDF and .JPEG output

# Should the user desire an alternative representation:
dev.new(width=9,height=18,unit='cm')
ggarrange('',mp.wing,mp.trunk,mp.cross,'','','',ncol=1,nrow=7,
labels=c('','','','','','',''),label.x=0,label.y=1.160,font.label=list(size=25),heights=c(0.1,1,1,1,0.1,0.275,0.3))
# Assemble plots
ggsave(filename='pairwise_plot_17_01_2022_no_subplot_labels.pdf')
ggsave(filename='pairwise_plot_17_01_2022_no_subplot_labels.jpeg',dpi=300)
# Save output as high quality .PDF and .JPEG output

# Should the user wish to save only the labelled legend for subsequent use in illustration applications:
dev.new(width=4,height=1,unit='cm')
ggarrange(mp.legend)
ggsave(filename='pairwise_plot_17_01_2022_just_legend.pdf')
ggsave(filename='pairwise_plot_17_01_2022_just_legend.jpeg',dpi=300)
# Save output as high quality .PDF and .JPEG output

# Script concludes

