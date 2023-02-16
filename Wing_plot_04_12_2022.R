# This script will generate Figure 3 from Orkney & Hedrick 2023

# The purpose of this script is to produce a series of plots
# illustrating how integration within the avian wing is structured
# by body mass, and how this structure over body mass is influenced
# by ecological relationships to flight style variety

# Call some relevant packages
library(phytools)
library(caTools)
library(geomorph)
library(nlme)
library(cluster)
# Done

# Set the work directory
setwd('D:/Documents/Alex_birds/Ideas_2022/SICBY')
# The user will need to customise this

# Load some requisite datasets
load('tree.22.10.2022.RData')
load('tree.names.22.10.2022.RData')
load('masses.22.10.2022.RData')
load('Csize.22.10.2022.RData')
# Done 

# Ensure names match the phylogeny. 
names(masses) <- tree_names
for( i in 1:length(GPA.Csize) ){
	names(GPA.Csize[[i]]) <- tree_names 
}
# Done

# Define a sundry function to interrogate data objects
get.item <- function( X , item ) { X[[ item ]] }
# This is a simple indexing function

# Define a function to compute landmark constellation centroid sizes, corrected for allometric scaling 
get.residual.Csize <- function( array, masses, phylogeny, taxa, bones ){ # Input arguments of shape-data array, mass vector, phylogeny, taxa and bones of interest
	allometry.Csize <- list() # Dumby variable defined to receive allometric models
	newphy <- drop.tip( phylogeny , phylogeny$tip.label[ !phylogeny$tip.label %in% taxa ] ) # Phylogeny pruned to taxa of interest
	for(i in 1:length(array[bones]) ){ # For each bone of interest 
		df <- data.frame( mass= log10(masses[taxa]), Csize = log10( array[[ bones[i] ]][taxa] ), species = taxa ) # Organise the data into a dataframe
		allometry.Csize[[ i ]] <- gls( Csize ~ mass, correlation = corPagel( 0, phy=newphy, form= ~ species  ), data=df ) # Compute allometric model
	}
	names(allometry.Csize) <- names(array[bones]) # Check correspondance of allometric models with bone names that describe them 
	residual_Csize <- lapply( lapply( allometry.Csize , get.item , item = "residuals" ) , c ) # Extract allometric model residuals 
	return(residual_Csize) # Return residuals 
}
# This function uses phylogenetic Generalised Least Squares to compute
# the allometrically adjusted centroid sizes of skeletal elements.
# The user must provide the array of original bird skeletal element centroid sizes, bodymasses,
# phylogeny, bones and taxa of interest. 

# The analyses in this script will require categorial flight style information that describes the bird species' ecologies;
# These are represented as a matrix of binary scores that will need to be transformed into a continuous representation
flight<- read.csv('flight_masses_22_10_2022_plus_A.csv') # Load the flight style data
flight.rownames<-flight$X
flight<-flight[,-c(1,2,14)]
rownames(flight)<-flight.rownames
flightdist <- daisy(flight, metric='gower', stand=F, type=list(asymm=c(1:11))) # Compute distance matrix
flightdist[which(is.na(flightdist)==T)]<-1 # Replace Na values with 1 
flightpcoa <- cmdscale(flightdist,eig=T) # Compute a Principal Coordinate Analysis to make the data continuously distributed 
k <- max(which(flightpcoa$eig/sum(flightpcoa$eig)>=0.05)) # Discover the lowest axis at which 95% of original variance is explained 
full.flightpcoa <- cmdscale(flightdist, k=k, eig=T) # Perform a truncation 
# Done 

# The following lines describe a bi-variate function that takes two vectors of skeletal element centroid sizes and computes
# the integration of their allometrically-adjusted residuals. The distance of individual taxa from the major axis is then extracted, 
# so that the relationship of this derived variable to body mass can be investigated: 
m.axis.resid <- function( array, masses, bone1, bone2, phylogeny, taxa ){ # This function takes shape data, a mass vector, bone names, phylogeny and taxa of interest inputs
	bones<-c(bone1,bone2) # Define a pairwise combination of bones 
	r.Csize <- get.residual.Csize( array = GPA.Csize, masses = masses, phylogeny = phylogeny, taxa= taxa,bones= bones ) # Compute residual centroid size
	newphy <- drop.tip( phylogeny , phylogeny$tip.label[ !phylogeny$tip.label %in% taxa ] ) # Prune phylogeny to taxa of interest 
	two.block <- phylo.integration(r.Csize[[bone1]],r.Csize[[bone2]],phy=newphy) # Compute 2 block partial least squares (integration) 
	residuals<- prcomp(cbind(two.block$XScores[,1],two.block$YScores[,1]))$x[,2] # Find residuals from the major axis of covaraince 
	list<-list()
	list[[1]]<-two.block
	list[[2]]<-residuals
	return(list) # Return residuals 
}
# This is a function which finds the major axis of covariance between two unique arrays of bone centroid sizes
# using phylogenetic 2-Block Partial Least Squares analysis...
# Taxonomic subsets and phylogenetic hypotheses can be input...
# Centroid sizes are adjusted for inter-specific allometric scaling in a phylogenetic context, and 
# body masses are therefore required...
# Euclidean distances of each species from the axis of major covariance are found and returned...  

# Input for the above function will now be defined
bone1<-'carpometacarpus'
bone2<-'humerus'
# Specify the bones of interest
# Let's explore integration within the wing

# Underttake analysis
output<-m.axis.resid( GPA.Csize, masses, bone1, bone2, pruned.tree, names(masses) )
# Find the distances from the major axis of covariance

# Construct a linear model between the residual distances and body masses.
model <- lm(output[[2]] ~log10(masses) )
# Done 

# Perform a Breusch-Pagan test, depending on body mass, to test for Heteroskedacity
car::ncvTest(model, ~log10(masses))$p
# A p-value is returned
# A value <0.05 may be taken as evidence of significance; that body mass variety causes a change in integration within the wing in birds 

# Load packages to be required to plot the data
library(ggplot2)
library(ggpubr)
# Done 
# Define a colour palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Colour Blind Palette

# Sundry plotting variables to be defined
 axis.title <- 15
 axis.text <- 10
 point.size <- 2
 line.width <- 1/2
 genus.text<-2.5
# Done 

# Set work directory as appropriate
setwd("D:/Documents/Alex_birds/")
# The user will need to customise 

# Define a vector of birds belonging to different taxonomic groups
Telluraves<-read.csv('Telluraves.csv')
# This is a pre-prepared vector of birds belonging to Telluraves
Apodiformes<- c('Archilochus_colubris',
'Topaza_pella',
'Chaetura_brachyura',
'Streptoprocne_zonaris',
'Hemiprocne_comata')
# This is a pre-prepared vector of birds belonging to Apodiformes

# Set work directory as appropriate
setwd('D:/Documents/Alex_birds/Ideas_2022/SICBY')
# The user will need to customise

# Define a vector to store the identities of taxonomic groups
groups <- rep('0',length=length(masses))
groups[match(Telluraves$x,names(masses))]<-'1'
groups[match(Apodiformes,names(masses))]<-'2'
cols <- c(cbbPalette[6],cbbPalette[7],cbbPalette[4]) # associated plotting point colours 
shapes <- c(15,18,17) # shapes 
labels <- c('Non-Telluraves','Telluraves','Apodiformes') # legend labels
# Done 

# Organise data for the first analysis into a dataframe. 
df <- as.data.frame(cbind(log10(masses),(output[[2]]), groups))
df$V1<-as.numeric(as.character(df$V1))
df$V2<-as.numeric(as.character(df$V2))
# Done 

# Compute 1 and 2 sigma confidence intervals for the data distribution
ribbons <- as.data.frame(cbind( rep(0,length(masses)),
log10(masses)[order(masses)],
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.95, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.68, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.32, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses)], k=30, probs= 0.05, endrule='quantile')
))
# Done 

# The names of outlying species across a range of body masses are desired
# for illustration purposes
smallest<-which(log10(masses)<1)
small<-which(log10(masses)>1 & log10(masses) <2)
big<-which(log10(masses)>2 & log10(masses) <3)
biggest<-which(log10(masses)>3)
smallest.out<-which(abs(output[[2]][smallest]) == max(abs(output[[2]][smallest])))
smallest[smallest.out]
small.out<-which(abs(output[[2]][small]) == max(abs(output[[2]][small])))
small[small.out]
big.out<-which(abs(output[[2]][big]) == max(abs(output[[2]][big])))
big[big.out]
biggest.out<-which(abs(output[[2]][biggest]) == max(abs(output[[2]][biggest])))
biggest[biggest.out]
names.out<-c(smallest[smallest.out],
small[small.out],big[big.out],biggest[biggest.out])
name.data<-as.data.frame(cbind(log10(masses)[names.out],output[[2]][names.out]))
# Done

# Species names must be italic
labs <- sapply(
  sub("\\_.*", "", names(names.out)), 
  function(x) parse(text = paste0("italic('", x[1], "')"))
)
# Done 

# The code for the first subplot will be commented extensively, and repetivive subplots thereafter will be commented sparingly:
# The following is Figure 3 subplot b 
mp <- ggplot()+ # Make plot
lims(x= c( min(log10(masses)),max(log10(masses)))+c(-0.1,0.1))+ # Axial limits
geom_ribbon(aes(x=V2,ymin=V6,ymax=V3),fill='gray90',data=ribbons )+ # Confidence envelope 1
geom_ribbon(aes(x=V2,ymin=V5,ymax=V4),fill='gray80',data=ribbons )+ # Confidence envelope 2 
geom_point( aes( x=V1, y=V2, col=groups,shape=groups ), size=point.size, stroke=1, data=df )+ # Species scatterplot
scale_color_manual( values= cols, labels = labels)+ # Colour species by taxonomic group
scale_shape_manual( values= shapes, labels = labels )+ # Point shape according to taxonomic group 
labs(y=expression(paste( '',italic('D'[m]),'' )),x=expression(paste(log[10],'(',italic(mass),')',' [g]')))+ # Axial labels 
# Sundry thematic instructions that regulate plot appearance:
theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=axis.text,colour='black'),
      axis.text.y=element_text(size=axis.text,colour='black'),
	axis.title.x.bottom=element_text(size=axis.title,colour="black"),
	axis.title.y.left=element_text(size=axis.title,colour="black"),
      axis.ticks=element_line(size=2),
      legend.position="right",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.key.width = unit(2, "cm"))
# Conclude plot

# Extract the legend
mp.legend<-cowplot::get_legend(mp)
# The package 'cowplot' must be installed

# The first subplot is going to be re-made with an overlay of silhouettes of outlying taxa
# This will require some more packages to be loaded:
library(png)
library(grid)
# Done 
# Phylopic files, from https://www.phylopic.org will be used as silhouettes:
Archilocus<-readPNG('Archilochus_colubris.png')
Archilocus <- rasterGrob(Archilocus, interpolate=TRUE)
Chaetura<-readPNG('Chaetura_brachyura.png')
Chaetura<- rasterGrob(Chaetura, interpolate=TRUE)
Puffinus<-readPNG('Puffinus_griseus.png')
Puffinus<-rasterGrob(Puffinus, interpolate=TRUE)
Phoebastria<-readPNG('Phoebastria_nigripes.png')
Phoebastria<-rasterGrob(Phoebastria, interpolate=TRUE)
# Loaded 

# The lines responsible for generating the overlay will be indicated with comments 
mp <- ggplot()+
lims(x= c( min(log10(masses)),max(log10(masses)))+c(-0.1,0.1))+
geom_ribbon(aes(x=V2,ymin=V6,ymax=V3),fill='gray90',data=ribbons )+
geom_ribbon(aes(x=V2,ymin=V5,ymax=V4),fill='gray80',data=ribbons )+
geom_point( aes( x=V1, y=V2, col=groups,shape=groups ), size=point.size, stroke=1, data=df )+
geom_text(data=name.data, aes(x=V1,y=V2) , label= labs , nudge_x=c(0.2,0.4,0.5,0.0),nudge_y=c(-0.003,0.005,0.002,-0.003), size=genus.text )+ # Add outlying taxa names
scale_color_manual( values= cols, labels = labels)+
scale_shape_manual( values= shapes, labels = labels )+
labs(y=expression(paste( '',italic('D'[m]),'' )),x=expression(paste(log[10],'(',italic(mass),')',' [g]')))+
# The following lines add silhouettes of outlying taxa as annotations: 
annotation_custom(Archilocus,xmin=name.data$V1[1]-0.2,xmax=name.data$V1[1]+0.2,
ymin=name.data$V2[1]-0.1,ymax=name.data$V2[1]+0.1)+
annotation_custom(Chaetura,xmin=name.data$V1[2]-0.2,xmax=name.data$V1[2]+0.2,
ymin=name.data$V2[2]-0.1,ymax=name.data$V2[2]+0.1)+
annotation_custom(Puffinus,xmin=name.data$V1[3]-0.2,xmax=name.data$V1[3]+0.2,
ymin=name.data$V2[3]-0.09,ymax=name.data$V2[3]+0.1)+
annotation_custom(Phoebastria,xmin=name.data$V1[4]-0.2,xmax=name.data$V1[4]+0.2,
ymin=name.data$V2[4]-0.1,ymax=name.data$V2[4]+0.1)+
# Done 
theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=axis.text,colour='black'),
      axis.text.y=element_text(size=axis.text,colour='black'),
	axis.title.x.bottom=element_text(size=axis.title,colour="black"),
	axis.title.y.left=element_text(size=axis.title,colour="black"),
      axis.ticks=element_line(size=2),
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.key.width = unit(2, "cm"))




# The following lines will prepare data to be plotted for Figure 3 subplot a 
df<- as.data.frame(cbind(log10(masses),output[[1]]$XScores[,1],output[[1]]$YScores[,1]))
# Compile data for first analysis into a dataframe.
pca<-prcomp(df[,2:3])
df<-cbind(df,pca$x)
# Compute principal component decomposition.
pca$rotation
svd<-(svd(df[,2:3])$v %*% diag(svd(df[,2:3])$d)) 
factor<-max(svd)/max(pca$x)
svd1<-(svd/factor)*0.8
# Arbitrary scaling 

# Ascribe new coordinates to the names of outlying taxa
name.data<-as.data.frame(cbind(df$PC1[names.out],df$PC2[names.out]))
rownames(name.data)<-names(names.out)
# Done 



mp.one<- # Produce plot Figure 3 subplot a; the first instance will be commented extensively and further subplots of the same type more sparingly
ggplot()+ # Make plot 
geom_point(data=df, aes(x=PC1,y=PC2), shape=21, size=3, fill='black')+ # Add a scatter of points representing individual species
geom_abline(intercept=0,slope= 0 )+ # Add a horizontal line representing the major axis 
geom_segment( aes(x=PC1,xend=PC1,y=0,yend=PC2, col= V1), data=df,lwd=3/2 )+ # Link the scatter of points to the major axis
scale_colour_gradient2(high='red',low='blue',mid='white',na.value='black',midpoint=2.25,limits=c(min(log10(masses)),max(log10(masses))) )+ # Colour the above lines according to body mass 
labs(x=expression(paste( '',italic('m'),'' )), y=expression(paste( '',italic('D'[m]),'' )),colour=expression(paste(log[10],'(',italic(mass),')',' [g]')))+ # Label axes and colour field
geom_segment(aes(x=c(0,0),xend=svd1[,1], y=c(0,0),yend=svd1[,2]))+ # Create vectors indicating direction of original skeletal variables 
geom_text(aes(x=svd1[,1]*1.1,y=svd1[,2]*1.1), label= c('carpometacarpus size','humerus size'))+ # Label those original vectors 
# Sundry thematic instructions that regulate plot appearance:
theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
	axis.title.x.bottom=element_text(size=20,colour="black"),
	axis.title.y.left=element_text(size=20,colour="black"),
      axis.ticks=element_line(size=2),
      legend.position="right",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.key.width = unit(2/3, "cm"))
# Plot concludes 


# Extract legend
mass.legend<-cowplot::get_legend(mp.one)
# package 'cowplot' must be installed



mp.one<- # Remake plot without legend
ggplot()+
geom_point(data=df, aes(x=PC1,y=PC2), shape=21, size=point.size-1, fill='black')+
geom_abline(intercept=0,slope= 0 )+
geom_segment( aes(x=PC1,xend=PC1,y=0,yend=PC2, col= V1), data=df,lwd=line.width )+
scale_colour_gradient2(high='red',low='blue',mid='white',na.value='black',midpoint=2.25,limits=c(min(log10(masses)),max(log10(masses))) )+
labs(x=expression(paste( '',italic('m'),'' )), y=expression(paste( '',italic('D'[m]),'' )),colour=expression(paste(log[10],'(',italic(mass),')',' [g]')))+
geom_segment(aes(x=c(0,0),xend=svd1[,1], y=c(0,0),yend=svd1[,2]))+
geom_text(aes(x=svd1[,1]*1.1,y=svd1[,2]*1.1), label= c('carpometacarpus size','humerus size'), size=genus.text, nudge_x=-max(svd1)*0.3)+
geom_text(data=name.data, aes(x=V1,y=V2) , label=  labs  , nudge_x=-max(svd1)*0.4, nudge_y=c(0), size=genus.text )+
theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
	axis.title.x.bottom=element_text(size=axis.title,colour="black"),
	axis.title.y.left=element_text(size=axis.title,colour="black"),
      axis.ticks=element_line(size=2),
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.key.width = unit(2/3, "cm"))



# The above process will now be repeated to investigate the way that body mass structures
# the integration of carpometacarpus size and a multivariate decomposition representing avian flight style variety 

# The above task will require the 'm.axis.resid' function to be altered subtly;
m.axis.resid <- function( array, masses, bone1, ecology, phylogeny, taxa ){
	r.Csize <- get.residual.Csize( array = GPA.Csize, masses = masses, phylogeny = phylogeny, taxa= taxa,bones= bone1 )
	newphy <- drop.tip( phylogeny , phylogeny$tip.label[ !phylogeny$tip.label %in% taxa ] )
	two.block <- phylo.integration(r.Csize[[bone1]],ecology,phy=newphy) # The function has been adapted to incorporate ecology
	residuals<- prcomp(cbind(two.block$XScores[,1],two.block$YScores[,1]))$x[,2] 
	list<-list()
	list[[1]]<-two.block
	list[[2]]<-residuals
	return(list)
}
# This time the function is adapted to compute relationship between bone and ecology

# The following lines perform the essential analysis:
output2<-m.axis.resid(GPA.Csize, masses, 
	'carpometacarpus', full.flightpcoa$points, pruned.tree, names(masses))
# Find the distances from the major axis of covariance
model <- lm(output2[[2]] ~log10(masses) )
# Construct a linear model between the residual distances and body masses 
car::ncvTest(model, ~log10(masses))$p
# Perform a Breusch-Pagan test, depending on body mass, to test for Heteroskedacity
# A p-value is returned.
# A value <0.05 may be taken as evidence of significance

# Data are prepared for plotting as before
df <- as.data.frame(cbind(log10(masses),(output2[[2]]), groups))
df$V1<-as.numeric(as.character(df$V1))
df$V2<-as.numeric(as.character(df$V2))
# Organise data for the first analysis into a dataframe

# Confidence envelopes are computed as before
ribbons <- as.data.frame(cbind( rep(0,length(masses)),
log10(masses)[order(masses)],
runquantile( x=(output2[[2]])[order(masses)], k=30, probs= 0.95, endrule='quantile'),
runquantile( x=(output2[[2]])[order(masses)], k=30, probs= 0.68, endrule='quantile'),
runquantile( x=(output2[[2]])[order(masses)], k=30, probs= 0.32, endrule='quantile'),
runquantile( x=(output2[[2]])[order(masses)], k=30, probs= 0.05, endrule='quantile')
))
# Compute 1 and 2 sigma confidence intervals for the data distribution

# Outlying taxa are found as before
smallest<-which(log10(masses)<1)
small<-which(log10(masses)>1 & log10(masses) <2)
big<-which(log10(masses)>2 & log10(masses) <3)
biggest<-which(log10(masses)>3)
smallest.out<-which(abs(output2[[2]][smallest]) == max(abs(output2[[2]][smallest])))
smallest[smallest.out]
small.out<-which(abs(output2[[2]][small]) == max(abs(output2[[2]][small])))
small[small.out]
big.out<-which(abs(output2[[2]][big]) == max(abs(output2[[2]][big])))
big[big.out]
biggest.out<-which(abs(output2[[2]][biggest]) == max(abs(output2[[2]][biggest])))
biggest[biggest.out]
names.out<-c(smallest[smallest.out],
small[small.out],big[big.out],biggest[biggest.out])
name.data<-as.data.frame(cbind(log10(masses)[names.out],output2[[2]][names.out]))
# Done

# Species names are italicised as before
labs <- sapply(
  sub("\\_.*", "", names(names.out)), 
  function(x) parse(text = paste0("italic('", x[1], "')"))
)
# Done

# The image files for outlying taxa are loaded, as before
Malurus<-readPNG('Malurus_splendens.png')
Malurus <- rasterGrob(Malurus, interpolate=TRUE)
Neodrepanis<-readPNG('Neodrepanis_coruscans.png')
Neodrepanis<- rasterGrob(Neodrepanis, interpolate=TRUE)
Coracias<-readPNG('Coracias_cyanogaster.png')
Coracias<-rasterGrob(Coracias, interpolate=TRUE)
Spheniscus<-readPNG('Spheniscus_humboldti.png')
Spheniscus<-rasterGrob(Spheniscus, interpolate=TRUE)
# Done

mp2 <- ggplot()+ # Figure 3 subplot c
lims(x= c( min(log10(masses)),max(log10(masses)))+c(-0.1,0.1))+
geom_ribbon(aes(x=V2,ymin=V6,ymax=V3),fill='gray90',data=ribbons )+
geom_ribbon(aes(x=V2,ymin=V5,ymax=V4),fill='gray80',data=ribbons )+
geom_point( aes( x=V1, y=V2, col=groups,shape=groups ), size=point.size, stroke=1, data=df )+
geom_text(data=name.data, aes(x=V1,y=V2) , label= labs , nudge_x=c(0.0,0.3,-0.5,-0.2),nudge_y=c(0.007,0.0075,0,-0.007), size=genus.text )+
scale_color_manual( values= cols, labels = labels)+
scale_shape_manual( values= shapes, labels = labels )+
labs(y=expression(paste( '',italic('D'[m]),'' )),x=expression(paste(log[10],'(',italic(mass),')',' [g]')))+
annotation_custom(Malurus,xmin=name.data$V1[1]-0.1,xmax=name.data$V1[1]+0.3,
ymin=name.data$V2[1]-0.1,ymax=name.data$V2[1]+0.1)+
annotation_custom(Neodrepanis,xmin=name.data$V1[2]+0.0,xmax=name.data$V1[2]+0.4,
ymin=name.data$V2[2]-0.1,ymax=name.data$V2[2]+0.1)+
annotation_custom(Coracias,xmin=name.data$V1[3]-0.05,xmax=name.data$V1[3]+0.35,
ymin=name.data$V2[3]-0.1,ymax=name.data$V2[3]+0.1)+
annotation_custom(Spheniscus,xmin=name.data$V1[4]-0.1,xmax=name.data$V1[4]+0.1,
ymin=name.data$V2[4]-0.1,ymax=name.data$V2[4]+0.1)+
theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=axis.text,colour='black'),
      axis.text.y=element_text(size=axis.text,colour='black'),
	axis.title.x.bottom=element_text(size=axis.title,colour="black"),
	axis.title.y.left=element_text(size=axis.title,colour="black"),
      axis.ticks=element_line(size=2),
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.key.width = unit(2, "cm"))
# Plot concludes

# Data is now prepared for Figure 3 subplot d, as before 
df<- as.data.frame(cbind(log10(masses),output2[[1]]$XScores[,1],output2[[1]]$YScores[,1]))
# Compile data for first analysis into a dataframe
pca<-prcomp(df[,2:3])
df<-cbind(df,pca$x)
# Compute principal component decomposition
pca$rotation
svd<-(svd(df[,2:3])$v %*% diag(svd(df[,2:3])$d)) 
factor<-max(svd)/max(pca$x)
svd2<-(svd/factor)*0.4
# Arbitrary scaling 

# Coordinates are ascribed for the outlying taxa names
name.data<-as.data.frame(cbind(df$PC1[names.out],df$PC2[names.out]))
rownames(name.data)<-names(names.out)
# Done 


mp.two <- ggplot()+ # Figure 3 subplot d 
geom_point(data=df, aes(x=PC1,y=PC2), shape=21, size=point.size-1, fill='black')+
geom_abline(intercept=0,slope= 0 )+
geom_segment( aes(x=PC1,xend=PC1,y=0,yend=PC2, col= V1), data=df,lwd=line.width )+
scale_colour_gradient2(high='red',low='blue',mid='white',na.value='black',midpoint=2.25,limits=c(min(log10(masses)),max(log10(masses))) )+
labs(x=expression(paste( '',italic('m'),'' )), y=expression(paste( '',italic('D'[m]),'' )),colour=expression(paste(log[10],'(',italic(mass),')',' [g]')))+
geom_segment(aes(x=c(0,0),xend=svd2[,1], y=c(0,0),yend=svd2[,2]))+
geom_text(aes(x=svd2[,1]*1.1,y=svd2[,2]*1.1), label= c('carpometacarpus size','flight style'), size=genus.text, nudge_x=c(max(svd)*0.3,max(svd2)*0.8), nudge_y=c(0,-max(svd2)*0.1))+
geom_text(data=name.data, aes(x=V1,y=V2) , label=  labs  , nudge_x=-max(svd2)*0.7, nudge_y=c(0), size=genus.text )+
theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
	axis.title.x.bottom=element_text(size=axis.title,colour="black"),
	axis.title.y.left=element_text(size=axis.title,colour="black"),
      axis.ticks=element_line(size=2),
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.key.width = unit(2/3, "cm"))
# Plot concludes 


# An investigation shall now be performed for carpometacarpus 3-D shape verses flight style ecology; does
# it reveal the same story as that evident from carpometacarpus size?
# This will require new data to be loaded and managed, and some functions to be re-defined...
# A simpler implementation is presented here, compared to the code to execute Figure 4; This implementation 
# is less flexible for user customistion (such as removing individual taxa)

# Load allometry-adjusted 3-D landmark constellation shape data and make sure names conform to phylogeny
load('allom_coords.29.10.2022.RData')
for( i in 1:length(allometry.residuals) ){
	dimnames(allometry.residuals[[i]])[[3]] <- tree_names
}
# Be aware that removing any outlying taxa will require the allometric models to be fitted again from scratch

# A new sundry function is defined, that will be used to generalise the 'm.axis.resid' function to higher dimensions for 3-D shape data
add <- function(x) Reduce("+", x)
# Done 

# The m.axis.resid function is redefined; comments will identify the specific changes
m.axis.resid <- function( array, masses, bone, ecology, phylogeny, taxa ){
	newphy <- drop.tip( phylogeny , phylogeny$tip.label[ !phylogeny$tip.label %in% taxa ] )
	two.block <- phylo.integration(array[[bone]],ecology,phy=newphy)
	residuals<-list()
	eig<-two.block$svd$d^2/sum(two.block$svd$d^2)
	for(i in 1:dim(two.block$XScores)[2]){ # For each pair of axes 
		residuals[[i]]<- prcomp(cbind(two.block$XScores[,i],two.block$YScores[,i]))$x[,2]*eig[i] # find the major axis of covariance and scale residuals by explained variance
	}
	residuals.multivar<-add(residuals) # Combine the scaled residuals of all axes into a single representation 
	list<-list()
	list[[1]]<-two.block
	list[[2]]<-residuals.multivar
	list[[3]]<-residuals
	return(list)
}
# This is a multivariate version of 'm.axis.resid', which was previously bivariate in nature 

# The following lines will compute an analysis asking how integration between 3-D carpometacarpus shape and flight style variety is influenced by
# body mass across the study birds
output3<- m.axis.resid(allometry.residuals, masses, 
	'carpometacarpus', full.flightpcoa$points, pruned.tree, names(masses))
# Residuals from major axis of covariance found
model <- lm(output3[[2]] ~log10(masses) )
# Construct a linear model between the residual distances and body masses
car::ncvTest(model, ~log10(masses))$p
# Perform a Breusch-Pagan test, depending on body mass, to test for Heteroskedacity
# A p-value is returned.
# A value <0.05 may be taken as evidence of significance

# The organisation, curation and plotting of data now follows as before:
df <- as.data.frame(cbind(log10(masses),(output3[[2]]), groups))
df$V1<-as.numeric(as.character(df$V1))
df$V2<-as.numeric(as.character(df$V2))
# Organise data for the first analysis into a dataframe

# Compute confidence envelopes
ribbons <- as.data.frame(cbind( rep(0,length(masses)),
log10(masses)[order(masses)],
runquantile( x=(output3[[2]])[order(masses)], k=30, probs= 0.95, endrule='quantile'),
runquantile( x=(output3[[2]])[order(masses)], k=30, probs= 0.68, endrule='quantile'),
runquantile( x=(output3[[2]])[order(masses)], k=30, probs= 0.32, endrule='quantile'),
runquantile( x=(output3[[2]])[order(masses)], k=30, probs= 0.05, endrule='quantile')
))
# Compute 1 and 2 sigma confidence intervals for the data distribution

# Identify outlying taxa
smallest<-which(log10(masses)<1)
small<-which(log10(masses)>1 & log10(masses) <2)
big<-which(log10(masses)>2 & log10(masses) <3)
biggest<-which(log10(masses)>3)
smallest.out<-which(abs(output3[[2]][smallest]) == max(abs(output3[[2]][smallest])))
smallest[smallest.out]
small.out<-which(abs(output3[[2]][small]) == max(abs(output3[[2]][small])))
small[small.out]
big.out<-which(abs(output3[[2]][big]) == max(abs(output3[[2]][big])))
big[big.out]
biggest.out<-which(abs(output3[[2]][biggest]) == max(abs(output3[[2]][biggest])))
biggest[biggest.out]
names.out<-c(smallest[smallest.out],
small[small.out],big[big.out],biggest[biggest.out])
name.data<-as.data.frame(cbind(log10(masses)[names.out],output3[[2]][names.out]))
# Done

# Italicise taxa names
labs <- sapply(
  sub("\\_.*", "", names(names.out)), 
  function(x) parse(text = paste0("italic('", x[1], "')"))
)
# Done

# Call outlying taxa silhouette images
Malurus<-readPNG('Malurus_splendens.png')
Malurus <- rasterGrob(Malurus, interpolate=TRUE)
Upupa<-readPNG('Upupa_epops.png')
Upupa<- rasterGrob(Upupa, interpolate=TRUE)
Porphyrio<-readPNG('Porphyrio_porphyrio.png')
Porphyrio<-rasterGrob(Porphyrio, interpolate=TRUE)
Spheniscus<-readPNG('Spheniscus_humboldti.png')
Spheniscus<-rasterGrob(Spheniscus, interpolate=TRUE)
# Done

mp3 <- ggplot()+ # Produce Figure 3 subplot f
lims(x= c( min(log10(masses)),max(log10(masses)))+c(-0.1,0.1))+
geom_ribbon(aes(x=V2,ymin=V6,ymax=V3),fill='gray90',data=ribbons )+
geom_ribbon(aes(x=V2,ymin=V5,ymax=V4),fill='gray80',data=ribbons )+
geom_point( aes( x=V1, y=V2, col=groups,shape=groups ), size=point.size, stroke=1, data=df )+
scale_color_manual( values= cols, labels = labels)+
scale_shape_manual( values= shapes, labels = labels )+
geom_text(data=name.data, aes(x=V1,y=V2) , label= labs , nudge_x=c(-0.3,-0.5,-0.5,-0.6), nudge_y=c(-0.001,0,0,0), size=genus.text )+
labs(y=expression(paste( '',italic('D'[m]),'' )),x=expression(paste(log[10],'(',italic(mass),')',' [g]')))+
annotation_custom(Malurus,xmin=name.data$V1[1]-0.2,xmax=name.data$V1[1]+0.2,
ymin=name.data$V2[1]-0.1,ymax=name.data$V2[1]+0.1)+
annotation_custom(Upupa,xmin=name.data$V1[2]-0.2,xmax=name.data$V1[2]+0.2,
ymin=name.data$V2[2]-0.098,ymax=name.data$V2[2]+0.1)+
annotation_custom(Porphyrio,xmin=name.data$V1[3]-0.2,xmax=name.data$V1[3]+0.2,
ymin=name.data$V2[3]-0.1,ymax=name.data$V2[3]+0.1)+
annotation_custom(Spheniscus,xmin=name.data$V1[4]-0.1,xmax=name.data$V1[4]+0.1,
ymin=name.data$V2[4]-0.099,ymax=name.data$V2[4]+0.1)+
theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=axis.text,colour='black'),
      axis.text.y=element_text(size=axis.text,colour='black'),
	axis.title.x.bottom=element_text(size=axis.title,colour="black"),
	axis.title.y.left=element_text(size=axis.title,colour="black"),
      axis.ticks=element_line(size=2),
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.key.width = unit(2, "cm"))
# Done 

# Organise data in preparation for production of Figure 3 subplot e 
df<- as.data.frame(cbind(log10(masses),output3[[1]]$XScores[,1],output3[[1]]$YScores[,1]))
# Compile data for first analysis into a dataframe
pca<-prcomp(df[,2:3])
df<-cbind(df,pca$x)
# Compute principal component decomposition
pca$rotation
svd<-(svd(df[,2:3])$v %*% diag(svd(df[,2:3])$d)) 
factor<-max(svd)/max(pca$x)
svd3<-(svd/factor)*0.1
# Arbitrary scaling
# Done

# Ascribe coordinates in new space for the outlying taxa names
name.data<-as.data.frame(cbind(df$PC1[names.out],df$PC2[names.out]))
rownames(name.data)<-names(names.out)
# Done 

mp.three <- ggplot()+ # Produce Figure 3 subplot e 
geom_point(data=df, aes(x=PC1,y=PC2), shape=21, size=point.size-1, fill='black')+
geom_abline(intercept=0,slope= 0 )+
geom_segment( aes(x=PC1,xend=PC1,y=0,yend=PC2, col= V1), data=df,lwd=line.width )+
scale_colour_gradient2(high='red',low='blue',mid='white',na.value='black',midpoint=2.25,limits=c(min(log10(masses)),max(log10(masses))) )+
labs(x=expression(paste( '',italic('m'),'' )), y=expression(paste( '',italic('D'[m]),'' )),colour=expression(paste(log[10],'(',italic(mass),')',' [g]')))+
geom_segment(aes(x=c(0,0),xend=svd3[,1], y=c(0,0),yend=svd3[,2]))+
geom_text(aes(x=svd3[,1]*1.1,y=svd3[,2]*1.1), label= c('carpometacarpus shape','flight style'), size=genus.text, nudge_x=-max(svd3)*0.3, nudge_y=c(0,-max(svd3)*0.1))+
geom_text(data=name.data, aes(x=V1,y=V2) , label=  labs  , nudge_x=c(max(svd3)*0.5,-max(svd3)*1.75,0,max(svd3)*2.9), nudge_y=c(-max(svd3)*0.07,0,max(svd3)*0.1,0), size=genus.text )+
theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
	axis.title.x.bottom=element_text(size=axis.title,colour="black"),
	axis.title.y.left=element_text(size=axis.title,colour="black"),
      axis.ticks=element_line(size=2),
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.key.width = unit(2/3, "cm"))
# Done

# Package 'ggpubr' will be required if it has not been loaded already

# Organise like subplots into columns
col.1<-ggarrange(mp.one,mp.two,mp.three,nrow=3,
labels=c('a','c','e'),label.x=-0.010,label.y=1.035,font.label=list(size=20))
col.2<-ggarrange(mp,mp2,mp3,nrow=3,align='hv',
labels=c('b','d','f'),label.x=-0.010,label.y=1.035,font.label=list(size=20))
col.3<-ggarrange(mp.legend,mass.legend,nrow=2,align='hv')
# Done

# Produce a plot space
dev.new(width=16.9,height=16.9,unit='cm')
# Done
# Arrange columns in plot space
ggarrange(col.1,col.2,col.3,ncol=3,widths=c(0.5,0.6,0.35))
# Done

getwd() # Check working directory
# ggsave(filename='wing_plot_17_01_2022.pdf')
# save the output

# Script concludes

