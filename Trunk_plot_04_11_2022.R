# This script produces Figure 4 from Orkney & Hedrick 2023

# The purpose of this script is to produce a series of plots
# illustrating how integration within the avian trunk is structured
# by body mass, and how this structure over body mass is influenced
# by ecological relationships to flight style variety...


# Call required packages
library(phytools)
library(caTools)
library(geomorph)
library(nlme)
library(cluster)
# Done

# Set the work directory
setwd('D:/Documents/Alex_birds/Ideas_2022/SICBY')
# The user will need to customise as appropriate

# Load some requisite datasets
load('tree.22.10.2022.RData')
load('tree.names.22.10.2022.RData')
load('masses.22.10.2022.RData')
load('Csize.22.10.2022.RData')
# Done

# Ensure taxa names match the phylogeny 
names(masses) <- tree_names
for( i in 1:length(GPA.Csize) ){
	names(GPA.Csize[[i]]) <- tree_names 
}
# Done

# Define a sundry function to interrogate the data objects:
get.item <- function( X , item ) { X[[ item ]] }
# This is a simple indexing function

# The following is a function that will compute residual centroid sizes from phylogenetically-informed 
# allometric models of skeletal element centroid size depending on body mass:
get.residual.Csize <- function( array, masses, phylogeny, taxa, bones ){ # The function takes a shape data array, mass vector, phylogeny, taxa and vector of skeletal element names as arguments
	allometry.Csize <- list() # Define a dumby variable to receive the allometric models
	newphy <- drop.tip( phylogeny , phylogeny$tip.label[ !phylogeny$tip.label %in% taxa ] ) # Prune phylogeny to supplied taxa 
	for(i in 1:length(array[bones]) ){ # For each of the supplied skeletal elements 
		df <- data.frame( mass= log10(masses[taxa]), Csize = log10( array[[ bones[i] ]][taxa] ), species = taxa ) # Define a data frame
		allometry.Csize[[ i ]] <- gls( Csize ~ mass, correlation = corPagel( 0, phy=newphy, form= ~ species  ), data=df ) # Compute allometric model
	}
	names(allometry.Csize) <- names(array[bones]) # Ensure names in dumby variable match supplied skeletal elements
	residual_Csize <- lapply( lapply( allometry.Csize , get.item , item = "residuals" ) , c ) # Extract residual centroid sizes from allometric models
	return(residual_Csize) # Return the residuals
} # Function concludes
# This function uses phylogenetic Generalised Least Squares to compute
# the allometrically adjusted centroid sizes of skeletal elements....
# The user must provide the array of original bird skeletal element centroid sizes, bodymasses,
# phylogeny, bones and taxa of interest


# Define a function to compute the residuals from the major axis of covariance that describes integration between 2 bones' centroid sizes, corrected for allometry
m.axis.resid <- function( array, masses, bone1, bone2, phylogeny, taxa ){ # The user must supply a shape data array, a vector of masses, skeletal element names, a phylogeny and list of taxa
	bones<-c(bone1,bone2) # Define the pairwise combination of bones 
	r.Csize <- get.residual.Csize( array = GPA.Csize, masses = masses, phylogeny = phylogeny, taxa= taxa,bones= bones ) # Accommodate allometry
	newphy <- drop.tip( phylogeny , phylogeny$tip.label[ !phylogeny$tip.label %in% taxa ] ) # Prune the phylogeny to the requested taxa 
	two.block <- phylo.integration(r.Csize[[bone1]],r.Csize[[bone2]],phy=newphy) # Compute two block partial least squares for allometry-corrected centroid sizes for the pairwise combination of bones
	residuals<- prcomp(cbind(two.block$XScores[,1],two.block$YScores[,1]))$x[,2] # Extract residuals from the major axis of covariance 
	list<-list()
	list[[1]]<-two.block
	list[[2]]<-residuals
	return(list) # Define a return the output as a list
}
# This is a function which finds the major axis of covariance between two unique arrays of bone centroid sizes 
# using phylogenetic 2-Block Partial Least Squares analysis... 
# Taxonomic subsets and phylogenetic hypotheses can be input...
# Centroid sizes are adjusted for inter-specific allometric scaling in a phylogenetic context, and 
# body masses are therefore required...
# Euclidean distances of each species from the axis of major covariance are found and returned.  


# The user may be interested in disqualifying outlying taxa such as Apodiformes 
#nApods <- names(masses)[-match(Apodiformes,names(masses))] # non Apodiformes

# Carpometacarpus-Humerus outliers include
Outliers<-c('none')
# no Outliers have been defined this time
nOutliers<-names(masses)

# Specify the bones of interest
bone1<-'sternum'
bone2<-'scapula'
# Integration within the wing

# Compute the residuals from the major axis of covariance
output<-m.axis.resid( GPA.Csize, masses[nOutliers], bone1, bone2, pruned.tree, nOutliers )
# Done
# Construct a linear model between the residual distances and body masses
model <- lm(output[[2]] ~log10(masses[nOutliers]) )
# How does residual distance from the major axis depend on body mass?
car::ncvTest(model, ~log10(masses[nOutliers]))$p
# Perform a Breusch-Pagan test, depending on body mass, to test for Heteroskedacity
# A p-value is returned
# A value <0.05 may be taken as evidence of significance; that variation in body mass causes a gradient in integration within the trunk across birds

# Load packages for plotting functions
library(ggplot2)
library(ggpubr)
# These packages will be required to plot the data

# Define a colour palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Colour Blind Palette

# Define some sundry plot parameters
 axis.title <- 15
 axis.text <- 10
 point.size <- 2
 line.width <- 1/2
 genus.text<-2.5
# Done

# Set work directory as appropriate
setwd("D:/Documents/Alex_birds/")
# The user must specify this as they see fit

# Load a vector of clade Telluraves 
Telluraves<-read.csv('Telluraves.csv')
# This is a pre-prepared list of birds belonging to Telluraves

# Define a vector of Apodiform birds
Apodiformes<- c('Archilochus_colubris',
'Topaza_pella',
'Chaetura_brachyura',
'Streptoprocne_zonaris',
'Hemiprocne_comata')
# This is a pre-prepared list of birds belonging to Apodiformes

# Set work directory as appropriate
setwd('D:/Documents/Alex_birds/Ideas_2022/SICBY')
# The user must specify this as they see fit

# Define a vector to store the identities of taxonomic groups
groups <- rep('0',length=length(masses[nOutliers]))
groups[match(Telluraves$x,names(masses[nOutliers]))]<-'1'
groups[match(Apodiformes,names(masses[nOutliers]))]<-'2'
# Done

# Define corresponding plot parameters for point colour, shape and legend labels
cols <- c(cbbPalette[6],cbbPalette[7],cbbPalette[4])
shapes <- c(15,18,17)
labels <- c('Non-Telluraves','Telluraves','Apodiformes')
# Done 

# Organise data for the first analysis into a dataframe
df <- as.data.frame(cbind(log10(masses[nOutliers]),(output[[2]]), groups))
df$V1<-as.numeric(as.character(df$V1))
df$V2<-as.numeric(as.character(df$V2))
# Done

# The following lines will compute a confidence envelope:
ribbons <- as.data.frame(cbind( rep(0,length(masses[nOutliers])),
log10(masses[nOutliers])[order(masses[nOutliers])],
runquantile( x=(output[[2]])[order(masses[nOutliers])], k=30, probs= 0.95, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses[nOutliers])], k=30, probs= 0.68, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses[nOutliers])], k=30, probs= 0.32, endrule='quantile'),
runquantile( x=(output[[2]])[order(masses[nOutliers])], k=30, probs= 0.05, endrule='quantile')
))
# Compute 1 and 2 sigma confidence intervals for the data distribution

# Obtain names of outlying species across a range of body masses
# for illustration purposes:
smallest<-which(log10(masses[nOutliers])<1)
small<-which(log10(masses[nOutliers])>1 & log10(masses[nOutliers]) <2)
big<-which(log10(masses[nOutliers])>2 & log10(masses[nOutliers]) <3)
biggest<-which(log10(masses[nOutliers])>3)
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
name.data<-as.data.frame(cbind(log10(masses[nOutliers])[names.out],output[[2]][names.out]))
# Done 

# Sundry label definition: make labels italic
labs <- sapply(
  sub("\\_.*", "", names(names.out)), 
  function(x) parse(text = paste0("italic('", x[1], "')"))
)
# Done


# The first instance of subplot in its respective class will be heavily commented, 
# and thereafter comments will desist
mp <- ggplot()+ # Produce Figure 4 subplot b 
lims(x= c( min(log10(masses)),max(log10(masses)))+c(-0.1,0.1))+ # Define axial limits 
geom_ribbon(aes(x=V2,ymin=V6,ymax=V3),fill='gray90',data=ribbons )+ # Plot first sigma confidence interval
geom_ribbon(aes(x=V2,ymin=V5,ymax=V4),fill='gray80',data=ribbons )+ # Plot second sigma confidence interval
geom_point( aes( x=V1, y=V2, col=groups,shape=groups ), size=point.size, stroke=1, data=df )+ # Add points for species 
scale_color_manual( values= cols, labels = labels)+ # Set point colour to taxanomic groups
scale_shape_manual( values= shapes, labels = labels )+ # Set point shape to taxonomic group 
labs(y=expression(paste( '',italic('D'[m]),'' )),x=expression(paste(log[10],'(',italic(mass),')',' [g]')))+ # Label axes 
# Sundry thematic instructions:
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
# Plot concludes

# Extract the legend
mp.legend<-cowplot::get_legend(mp)
# Package 'cowplot' must be installed

# Load packages to process images and arrange plots
library(png)
library(grid)
# Done 

# Read phylopic images for outlying taxa: (These can be sourced and downloaded independently from the phylopic website: https://www.phylopic.org
Archilocus<-readPNG('Archilochus_colubris.png')
Archilocus <- rasterGrob(Archilocus, interpolate=TRUE)
Hymenops<-readPNG('Hymenops_perspicillatus.png')
Hymenops<- rasterGrob(Hymenops, interpolate=TRUE)
Rollandia<-readPNG('Rollandia_rolland.png')
Rollandia<-rasterGrob(Rollandia, interpolate=TRUE)
Gavia<-readPNG('Gavia_immer.png')
Gavia<-rasterGrob(Gavia, interpolate=TRUE)
# Done 

# The following lines will reproduce Figure 4 subplot b, with the exception that images of outlying taxa will be superimposed 
# Only those lines will be commented 
mp <- ggplot()+
lims(x= c( min(log10(masses)),max(log10(masses)))+c(-0.1,0.1))+
geom_ribbon(aes(x=V2,ymin=V6,ymax=V3),fill='gray90',data=ribbons )+
geom_ribbon(aes(x=V2,ymin=V5,ymax=V4),fill='gray80',data=ribbons )+
geom_point( aes( x=V1, y=V2, col=groups,shape=groups ), size=point.size, stroke=1, data=df )+
geom_text(data=name.data, aes(x=V1,y=V2) , label= labs , nudge_x=c(0.5,-0.2,0.5,0.5), nudge_y=c(0.003,-0.003,-0.002,0), size=genus.text )+ # Add text for outlying taxa
scale_color_manual( values= cols, labels = labels)+
scale_shape_manual( values= shapes, labels = labels )+
labs(y=expression(paste( '',italic('D'[m]),'' )),x=expression(paste(log[10],'(',italic(mass),')',' [g]')))+
# The following lines add the silhuoettes of outlying taxa 
annotation_custom(Archilocus,xmin=name.data$V1[1]-0.1,xmax=name.data$V1[1]+0.3,
ymin=name.data$V2[1]-0.1,ymax=name.data$V2[1]+0.1)+
annotation_custom(Hymenops,xmin=name.data$V1[2]-0.1,xmax=name.data$V1[2]+0.3,
ymin=name.data$V2[2]-0.1,ymax=name.data$V2[2]+0.1)+
annotation_custom(Rollandia,xmin=name.data$V1[3]-0.2,xmax=name.data$V1[3]+0.2,
ymin=name.data$V2[3]-0.1,ymax=name.data$V2[3]+0.1)+
annotation_custom(Gavia,xmin=name.data$V1[4]-0.2,xmax=name.data$V1[4]+0.2,
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
# plot remade with outlying taxa and without the legend 


# The following lines will prepare data for Figure 4 subplot a 
df<- as.data.frame(cbind(log10(masses[nOutliers]),output[[1]]$XScores[,1],output[[1]]$YScores[,1]))
# Compile data for first analysis into a dataframe
pca<-prcomp(df[,2:3])
df<-cbind(df,pca$x)
# Compute principal component decomposition
pca$rotation
svd<-(svd(df[,2:3])$v %*% diag(svd(df[,2:3])$d)) 
factor<-max(svd)/max(pca$x)
svd1<-(svd/factor)*0.8
# Sundry decisions such as scaling 

# Define coordinates of outlying taxa
name.data<-as.data.frame(cbind(df$PC1[names.out],df$PC2[names.out]))
rownames(name.data)<-names(names.out)
# Dome 


mp.one<- # The following produces a plot for Figure 4 subplot a; salient differences will be highlighted 
ggplot()+
geom_point(data=df, aes(x=PC1,y=PC2), shape=21, size=3, fill='black')+
geom_abline(intercept=0,slope= 0 )+
geom_segment( aes(x=PC1,xend=PC1,y=0,yend=PC2, col= V1), data=df,lwd=3/2 )+ # Vectors connecting taxa to PC1
scale_colour_gradient2(high='red',low='blue',mid='white',na.value='black',midpoint=2.25,limits=c(min(log10(masses)),max(log10(masses))) )+ # Line colours indicating taxa masses
labs(x=expression(paste( '',italic('m'),'' )), y=expression(paste( '',italic('D'[m]),'' )),colour=expression(paste(log[10],'(',italic(mass),')',' [g]')))+
geom_segment(aes(x=c(0,0),xend=svd1[,1], y=c(0,0),yend=svd1[,2]))+ # Vectors describing original input variable directions in PC-space 
geom_text(aes(x=svd1[,1]*1.1,y=svd1[,2]*1.1), label= c('sternum size','scapula size'))+ # Label these vectors 
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


mass.legend<-cowplot::get_legend(mp.one)



mp.one<- # Figure 4 subplot a will now be re-compiled without its legend
ggplot()+
geom_point(data=df, aes(x=PC1,y=PC2), shape=21, size=point.size-1, fill='black')+
geom_abline(intercept=0,slope= 0 )+
geom_segment( aes(x=PC1,xend=PC1,y=0,yend=PC2, col= V1), data=df,lwd=line.width )+
scale_colour_gradient2(high='red',low='blue',mid='white',na.value='black',midpoint=2.25,limits=c(min(log10(masses)),max(log10(masses))) )+
labs(x=expression(paste( '',italic('m'),'' )), y=expression(paste( '',italic('D'[m]),'' )),colour=expression(paste(log[10],'(',italic(mass),')',' [g]')))+
geom_segment(aes(x=c(0,0),xend=svd1[,1], y=c(0,0),yend=svd1[,2]))+
geom_text(aes(x=svd1[,1]*1.1,y=svd1[,2]*1.1), label= c('sternum size','scapula size'), size=genus.text, nudge_x= -max(svd1)*0.4)+
geom_text(data=name.data, aes(x=V1,y=V2) , label=  labs  , nudge_x=c(max(svd1)*0.2,max(svd1)*0.4,max(svd1)*0.5,max(svd1)*0.30), nudge_y=c(max(svd1)*0.055,-max(svd1)*0.001,0,0), size=genus.text )+
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


# This process will now be repeated for the remaining subplots, as such comments will be lightly distributed

nOutliers<-names(masses)

bone1<-'carpometacarpus'
bone2<-'sternum'
# A new pairwise combination of bones is being considered

output2<-m.axis.resid( GPA.Csize, masses[nOutliers], bone1, bone2, pruned.tree, nOutliers )

# Find the distances from the major axis of covariance. 
model <- lm(output2[[2]] ~log10(masses[nOutliers]) )
# Construct a linear model between the residual distances and body masses. 
car::ncvTest(model, ~log10(masses[nOutliers]))$p
# Perform a Breusch-Pagan test, depending on body mass, to test for Heteroskedacity. 
# A p-value is returned.
# A value <0.05 may be taken as evidence of significance. 


setwd("D:/Documents/Alex_birds/")
# Set work directory as appropriate

Telluraves<-read.csv('Telluraves.csv')
# This is a pre-prepared list of birds belonging to Telluraves

Apodiformes<- c('Archilochus_colubris',
'Topaza_pella',
'Chaetura_brachyura',
'Streptoprocne_zonaris',
'Hemiprocne_comata')
# This is a pre-prepared list of birds belonging to Apodiformes.

setwd('D:/Documents/Alex_birds/Ideas_2022/SICBY')
# Set work directory as appropriate

groups <- rep('0',length=length(masses[nOutliers]))
# Define a vector to store the identities of taxonomic groups
groups[match(Telluraves$x,names(masses[nOutliers]))]<-'1'
groups[match(Apodiformes,names(masses[nOutliers]))]<-'2'
cols <- c(cbbPalette[6],cbbPalette[7],cbbPalette[4])
shapes <- c(15,18,17)
labels <- c('Non-Telluraves','Telluraves','Apodiformes')

df <- as.data.frame(cbind(log10(masses)[nOutliers],(output2[[2]]), groups))
df$V1<-as.numeric(as.character(df$V1))
df$V2<-as.numeric(as.character(df$V2))
# Organise data for the first analysis into a dataframe. 

ribbons <- as.data.frame(cbind( rep(0,length(masses[nOutliers])),
log10(masses[nOutliers])[order(masses[nOutliers])],
runquantile( x=(output2[[2]])[order(masses[nOutliers])], k=30, probs= 0.95, endrule='quantile'),
runquantile( x=(output2[[2]])[order(masses[nOutliers])], k=30, probs= 0.68, endrule='quantile'),
runquantile( x=(output2[[2]])[order(masses[nOutliers])], k=30, probs= 0.32, endrule='quantile'),
runquantile( x=(output2[[2]])[order(masses[nOutliers])], k=30, probs= 0.05, endrule='quantile')
))
# Compute 1 and 2 sigma confidence intervals for the data distribution. 


# We also want the names of outlying species across a range of body masses
# for illustration purposes

smallest<-which(log10(masses[nOutliers])<1)
small<-which(log10(masses[nOutliers])>1 & log10(masses[nOutliers]) <2)
big<-which(log10(masses[nOutliers])>2 & log10(masses[nOutliers]) <3)
biggest<-which(log10(masses[nOutliers])>3)

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

name.data<-as.data.frame(cbind(log10(masses[nOutliers])[names.out],output2[[2]][names.out]))

labs <- sapply(
  sub("\\_.*", "", names(names.out)), 
  function(x) parse(text = paste0("italic('", x[1], "')"))
)


Acanthisitta<-readPNG('Acanthisitta_chloris.png')
Acanthisitta <- rasterGrob(Acanthisitta, interpolate=TRUE)
Neodrepanis<-readPNG('Neodrepanis_coruscans.png')
Neodrepanis<- rasterGrob(Neodrepanis, interpolate=TRUE)
Crypturellus<-readPNG('Crypturellus_tataupa.png')
Crypturellus<-rasterGrob(Crypturellus, interpolate=TRUE)
Cathartes<-readPNG('Cathartes_burrovianus.png')
Cathartes<-rasterGrob(Cathartes, interpolate=TRUE)


mp2 <- ggplot()+ # Figure 4 subplot d
lims(x= c( min(log10(masses)),max(log10(masses)))+c(-0.1,0.1))+
geom_ribbon(aes(x=V2,ymin=V6,ymax=V3),fill='gray90',data=ribbons )+
geom_ribbon(aes(x=V2,ymin=V5,ymax=V4),fill='gray80',data=ribbons )+
geom_point( aes( x=V1, y=V2, col=groups,shape=groups ), size=point.size, stroke=1, data=df )+
geom_text(data=name.data, aes(x=V1,y=V2) , label= labs , nudge_x=c(0.0,-0.4,-0.6,0.5), nudge_y=c(0.003,-0.0015,0,0), size=genus.text )+
scale_color_manual( values= cols, labels = labels)+
scale_shape_manual( values= shapes, labels = labels )+
labs(y=expression(paste( '',italic('D'[m]),'' )),x=expression(paste(log[10],'(',italic(mass),')',' [g]')))+
annotation_custom(Acanthisitta,xmin=name.data$V1[1]-0.1,xmax=name.data$V1[1]+0.1,
ymin=name.data$V2[1]-0.1,ymax=name.data$V2[1]+0.1)+
annotation_custom(Neodrepanis,xmin=name.data$V1[2]+0.0,xmax=name.data$V1[2]+0.4,
ymin=name.data$V2[2]-0.1,ymax=name.data$V2[2]+0.1)+
annotation_custom(Crypturellus,xmin=name.data$V1[3]-0.05,xmax=name.data$V1[3]+0.35,
ymin=name.data$V2[3]-0.099,ymax=name.data$V2[3]+0.101)+
annotation_custom(Cathartes,xmin=name.data$V1[4]-0.2,xmax=name.data$V1[4]+0.2,
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


df<- as.data.frame(cbind(log10(masses[nOutliers]),output2[[1]]$XScores[,1],output2[[1]]$YScores[,1]))
# Compile data for first analysis into a dataframe.
pca<-prcomp(df[,2:3])
df<-cbind(df,pca$x)
# Compute principal component decomposition.
pca$rotation
svd<-(svd(df[,2:3])$v %*% diag(svd(df[,2:3])$d)) 
factor<-max(svd)/max(pca$x)
svd2<-(svd/factor)*0.7

name.data<-as.data.frame(cbind(df$PC1[names.out],df$PC2[names.out]))
rownames(name.data)<-names(names.out)


mp.two <- ggplot()+ # Figure 4 subplot c
geom_point(data=df, aes(x=PC1,y=PC2), shape=21, size=point.size-1, fill='black')+
geom_abline(intercept=0,slope= 0 )+
geom_segment( aes(x=PC1,xend=PC1,y=0,yend=PC2, col= V1), data=df,lwd=line.width )+
scale_colour_gradient2(high='red',low='blue',mid='white',na.value='black',midpoint=2.25,limits=c(min(log10(masses)),max(log10(masses))) )+
labs(x=expression(paste( '',italic('m'),'' )), y=expression(paste( '',italic('D'[m]),'' )),colour=expression(paste(log[10],'(',italic(mass),')',' [g]')))+
geom_segment(aes(x=c(0,0),xend=svd2[,1], y=c(0,0),yend=svd2[,2]))+
geom_text(aes(x=svd2[,1]*1.1,y=svd2[,2]*1.1), label= c('carpometacarpus size','sternum size'), size=genus.text, nudge_x=c(-max(svd2)*0.35,0), nudge_y=c(max(svd2)*0.03,-max(svd2)*0.03))+
geom_text(data=name.data, aes(x=V1,y=V2) , label=  labs  , nudge_x=c(0,0,-max(svd2)*0.5,-max(svd2)*0.5), nudge_y=c(max(svd2)*0.05,max(svd2)*0.1,0,-max(svd2)*0.05), size=genus.text )+
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



# A plot will now be produced for 3-D coracoid shape verses flight style ecology, and the way in which this relationship is structured by body mass 

nOutliers <- names(masses)

setwd("D:/Documents/Alex_birds/")
# Set work directory as appropriate

Telluraves<-read.csv('Telluraves.csv')
# This is a pre-prepared list of birds belonging to Telluraves.

Apodiformes<- c('Archilochus_colubris',
'Topaza_pella',
'Chaetura_brachyura',
'Streptoprocne_zonaris',
'Hemiprocne_comata')
# This is a pre-prepared list of birds belonging to Apodiformes.

setwd('D:/Documents/Alex_birds/Ideas_2022/SICBY')
# Set work directory as appropriate

groups <- rep('0',length=length(masses[nOutliers]))
# Define a vector to store the identities of taxonomic groups
groups[match(Telluraves$x,names(masses[nOutliers]))]<-'1'
groups[match(Apodiformes,names(masses[nOutliers]))]<-'2'
cols <- c(cbbPalette[6],cbbPalette[7],cbbPalette[4])
shapes <- c(15,18,17)
labels <- c('Non-Telluraves','Telluraves','Apodiformes')

# A denser commenting style resumes here: Figure 4 subplots e and f are computed in a slightly different way:

# A dataset of flight style binary scores will be loaded, transformed by Principal Coordinate analysis, and truncated at 95% explained variance
flight<- read.csv('flight_masses_22_10_2022_plus_A.csv') # The document 'flight_masses_22_10_2022_plus_A.csv' is required 
flight.rownames<-flight$X
flight<-flight[match(nOutliers,flight.rownames),]
flight<-flight[,-c(1,2,14)]
rownames(flight)<-flight.rownames[match(nOutliers,flight.rownames)]
flightdist <- daisy(flight, metric='gower', stand=F, type=list(asymm=c(1:11))) # Produce a distance matrix 
flightdist[which(is.na(flightdist)==T)]<-1 # Replace Na values with '1'
flightpcoa <- cmdscale(flightdist,eig=T) # Compute Principal Coordinate analysis
k <- max(which(flightpcoa$eig/sum(flightpcoa$eig)>=0.05)) # Find lowest axis at with 95% of variance is explained 
full.flightpcoa <- cmdscale(flightdist, k=k, eig=T) # Truncate
# Done 

# Original 3-D landmark constellation data may be required for the following analyses;
# Sometimes the metadata includes spelling mistakes, so a bespoke document that corrects them is required
name.matches<-read.csv('name_matches_12_10_2022.csv') # The document 'name_matches_12_10_2022.csv' is needed
# Loaded

# Find matches between original taxa and phylogeny tip labels
match(nOutliers,name.matches$tree_names)
# Done

# In order to process the 3-D shape data we will need to load Alex Bjarnason's landmarks from scratch, so that
# the GPA and allometric models are not distorted by the outliers we have removed (if a decision is made to do this) 
load("D:/Documents/Alex_birds/NEW_analysis/2021_Bjarnason_125_SIdata/Supp 3 Datasets for Bjarnason & Benson 2021 revised submission 26 January 2021/Bird_landmarks_clean_26January2021_Bjarnason&Benson taxon set.RData")
# Done

# Load 3-D landmark constellation metadata:
landmark.info<-read.csv("D:\\Documents\\Alex_birds\\Bird_landmark_info_23Nov2019.csv" , row.names = 1)
# Done 

# Load package 'ape' to process phylogenies
library(ape)
# Done

# We are now loading the Avian Phylogeny of Prum et al., (2015)
setwd("D:/Documents/Alex_birds/")
# Path will need to be customed by the user 
PrumTree <- read.nexus(file = "Avian-TimeTree.tre")
# Done 

setwd('D:/Documents/Alex_birds/NEW_analysis/')
# Path will need to be customed by the user
# Load required meta data
metadata <- read.csv( "Eco_meta_data.csv" , row.names = 1 )
# Done (this contains information about bird masses and ecological properties)

# Set number of landmarks in semilandmark curves to minimum 
SL.counts <- "min"
# Done 

# Define list of skeletal elements
elements <- c( "skull" , "mandible" , "scapula", "coracoid", "sternum" , "humerus", "ulna", "radius", "carpometacarpus", "synsacrum" , "femur", "tibiotarsus", "tarsometatarsus" )
# Done

# The following lines will correct errors in the original metadata list, remove duplicate entries and standardise label spelling
name_list<-list()
for(i in 1:length(elements)){
	name_list[[i]]<-dimnames(compiled.bird.landmarks.list[[ elements[ i ] ]][[ SL.counts ]][[ "output.landmarks.array" ]] )[[ 3 ]]
}
AA<-table(unlist(name_list))
birds <- names( AA )[ AA == length( elements ) ]
birds <- birds[ birds != "Falco_sparverius_UMMZ_154452" ]

birds<-birds[match(nOutliers,name.matches$tree_names)]
# Done

# The following lines will perform Generalised Procrustes Alignments on the landmark constellations for each unique skeletal element
GPA.results <- list()
	for( i in 1:length( elements ) )	{
		element <- elements [ i ]
			GPA.results[[ i ]] <- gpagen( compiled.bird.landmarks.list[[ element ]][[ SL.counts ]][[ "output.landmarks.array" ]][ , , birds ] ,  curves = compiled.bird.landmarks.list[[ element ]][[ SL.counts ]][[ "sliders" ]] , ProcD = T )
		}
names( GPA.results ) <- elements
# Done


# We have performed the rotation to align the skeletons; We have rejected birds that do not appear in the Meta data 
# All of the rotated constellations host nOutliers birds (Birds that are not considered as outliers)

# Define a sundry function to interrogate existing data objects
get.item <- function( X , item ) { X[[ item ]] }
# Done

# Extract aligned landmark coordinates and centroid sizes	
GPA.coords <- lapply( GPA.results , get.item , item = "coords" )
GPA.Csize <- lapply( GPA.results , get.item , item = "Csize" )
#Done


# The following lines will compute phylogenetically-informed allometric models and correct shape data for allometry
# The package 'abind' will be needed as part of this process
library(abind)

allometry_list <- list()
allometry.residuals <- list()

for(i in 1:length(GPA.coords)){
	shape.temp<-GPA.coords[[i]]
	dimnames(shape.temp)[[3]]<-nOutliers #tree_names
	nOutliers.tree<-drop.tip( pruned.tree , pruned.tree$tip.label[ !pruned.tree$tip.label %in% nOutliers ] )
	gdf.temp<-geomorph.data.frame(shape=shape.temp,size=log10(masses[nOutliers]),phy=nOutliers.tree)
	allometry_list[[i]]<-procD.pgls(shape~size,phy=nOutliers.tree,data=gdf.temp)
	names(allometry_list)[[i]]<-elements[i]
	for(j in 1:nrow(allometry_list[[i]]$pgls.residuals)){
		if(j==1){
			allometry.residuals[[i]]<-matrix(allometry_list[[elements[i]]]$pgls.residuals[j,],ncol=3,byrow=T)+matrix(t(GPA.results[[elements[i]]]$consensus),ncol=3,byrow=T)
			rownames(allometry.residuals[[i]])<-dimnames(shape.temp)[[1]]
			colnames(allometry.residuals[[i]])<-c('x','y','z')
		} else {
		allometry.residuals[[i]]<-abind(allometry.residuals[[i]],matrix(allometry_list[[elements[i]]]$pgls.residuals[j,],ncol=3,byrow=T)+matrix(t(GPA.results[[elements[i]]]$consensus),ncol=3,byrow=T),along=3)
 		}
	}
}
names(allometry.residuals)<-names(allometry_list)<-names(GPA.coords)
# The 3-D landmark constellations have been corrected for allometric shape changes

The following lines will ensure that taxa names acros these original and derives shape data objects are standardised and conform to the phylogeny
for(i in 1:length(GPA.coords)){
	dimnames(allometry.residuals[[i]])<-dimnames(GPA.coords[[i]])
}

for( i in 1:length(allometry.residuals) ){
	dimnames(allometry.residuals[[i]])[[3]] <- nOutliers
}
# Done 

# A sundry function will be defined to help generalise an existing bivariate function to an n dimensional function
add <- function(x) Reduce("+", x)
# Done

# The 'get taxa residuals from the major axis of covariance' function is not suitable in its current state to interrogate data that has a dimension above 2,
# Therefore it requires generalisation; the lines that achieve this will be highlighted 

m.axis.resid <- function( array, masses, bone, ecology, phylogeny, taxa ){
	newphy <- drop.tip( phylogeny , phylogeny$tip.label[ !phylogeny$tip.label %in% taxa ] )
	two.block <- phylo.integration(array[[bone]],ecology,phy=newphy)
	residuals<-list()
	eig<-two.block$svd$d^2/sum(two.block$svd$d^2)
	# Begin generalisation 
	for(i in 1:dim(two.block$XScores)[2]){ # For each pair of axes in an n dimensional space 
		residuals[[i]]<- prcomp(cbind(two.block$XScores[,i],two.block$YScores[,i]))$x[,2]*eig[i] # Compute principal axes and scale by explained variance
	}
	residuals.multivar<-add(residuals) # Combine the multiple axes pairs into a single metric
	# Generalisation complete 
	list<-list()
	list[[1]]<-two.block
	list[[2]]<-residuals.multivar
	list[[3]]<-residuals
	return(list)
}
# This is a multivariate version of a previous function 'm.axis.resid'

# The above function will now be applied to discover how the strength of integration between 3-D coracoid shape and flight style ecology changes with body mass across birds
output3<- m.axis.resid(allometry.residuals, masses[nOutliers], 
	'coracoid', full.flightpcoa$points, pruned.tree, nOutliers)
# Analysis completed

model <- lm(output3[[2]] ~log10(masses[nOutliers]) )
# Construct a linear model between the residual distances and body masses 
car::ncvTest(model, ~log10(masses[nOutliers]))$p
# Perform a Breusch-Pagan test, depending on body mass, to test for Heteroskedacity
# A p-value is returne
# A value <0.05 may be taken as evidence of significance

# The organisation of resultant analytical output for plotting now proceeds as before:

df <- as.data.frame(cbind(log10(masses),(output3[[2]]), groups))
df$V1<-as.numeric(as.character(df$V1))
df$V2<-as.numeric(as.character(df$V2))
# Organise data for the first analysis into a dataframe. 

ribbons <- as.data.frame(cbind( rep(0,length(masses[nOutliers])),
log10(masses)[order(masses)],
runquantile( x=(output3[[2]])[order(masses)], k=30, probs= 0.95, endrule='quantile'),
runquantile( x=(output3[[2]])[order(masses)], k=30, probs= 0.68, endrule='quantile'),
runquantile( x=(output3[[2]])[order(masses)], k=30, probs= 0.32, endrule='quantile'),
runquantile( x=(output3[[2]])[order(masses)], k=30, probs= 0.05, endrule='quantile')
))
# Compute 1 and 2 sigma confidence intervals for the data distribution. 


smallest<-which(log10(masses[nOutliers])<1)
small<-which(log10(masses[nOutliers])>1 & log10(masses[nOutliers]) <2)
big<-which(log10(masses[nOutliers])>2 & log10(masses[nOutliers]) <3)
biggest<-which(log10(masses[nOutliers])>3)

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

name.data<-as.data.frame(cbind(log10(masses[nOutliers])[names.out],output3[[2]][names.out]))

labs <- sapply(
  sub("\\_.*", "", names(names.out)), 
  function(x) parse(text = paste0("italic('", x[1], "')"))
)

setwd('D:/Documents/Alex_birds/Ideas_2022/SICBY')

Archilochus<-readPNG('Archilochus_colubris.png')
Archilochus <- rasterGrob(Archilochus, interpolate=TRUE)
Rhynchocyclus<-readPNG('Rhynchocyclus_olivaceus.png')
Rhynchocyclus<- rasterGrob(Rhynchocyclus, interpolate=TRUE)
Puffinus<-readPNG('Puffinus_griseus.png')
Puffinus<-rasterGrob(Puffinus, interpolate=TRUE)
Phoebastria<-readPNG('Phoebastria_nigripes.png')
Phoebastria<-rasterGrob(Phoebastria, interpolate=TRUE)


mp3 <- ggplot()+ # Figure 4 subplot f 
lims(x= c( min(log10(masses)),max(log10(masses)))+c(-0.1,0.1))+
geom_ribbon(aes(x=V2,ymin=V6,ymax=V3),fill='gray90',data=ribbons )+
geom_ribbon(aes(x=V2,ymin=V5,ymax=V4),fill='gray80',data=ribbons )+
geom_point( aes( x=V1, y=V2, col=groups,shape=groups ), size=point.size, stroke=1, data=df )+
scale_color_manual( values= cols, labels = labels)+
scale_shape_manual( values= shapes, labels = labels )+
geom_text(data=name.data, aes(x=V1,y=V2) , label= labs , nudge_x=c(0.3,-0.45,-0.5,-0.0), nudge_y=c(-0.004,-0.003,0,-0.0025), size=genus.text )+
labs(y=expression(paste( '',italic('D'[m]),'' )),x=expression(paste(log[10],'(',italic(mass),')',' [g]')))+
annotation_custom(Archilochus,xmin=name.data$V1[1]-0.2,xmax=name.data$V1[1]+0.2,
ymin=name.data$V2[1]-0.102,ymax=name.data$V2[1]+0.1)+
annotation_custom(Rhynchocyclus,xmin=name.data$V1[2]-0.15,xmax=name.data$V1[2]+0.05,
ymin=name.data$V2[2]-0.0109,ymax=name.data$V2[2]+0.01)+
annotation_custom(Puffinus,xmin=name.data$V1[3]-0.2,xmax=name.data$V1[3]+0.2,
ymin=name.data$V2[3]-0.098,ymax=name.data$V2[3]+0.1)+
annotation_custom(Phoebastria,xmin=name.data$V1[4]-0.2,xmax=name.data$V1[4]+0.2,
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

df<- as.data.frame(cbind(log10(masses[nOutliers]),output3[[1]]$XScores[,1],output3[[1]]$YScores[,1]))
# Compile data for first analysis into a dataframe.
pca<-prcomp(df[,2:3])
df<-cbind(df,pca$x)
# Compute principal component decomposition.
pca$rotation
svd<-(svd(df[,2:3])$v %*% diag(svd(df[,2:3])$d)) 
factor<-max(svd)/max(pca$x)
svd3<-(svd/factor)*0.7

name.data<-as.data.frame(cbind(df$PC1[names.out],df$PC2[names.out]))
rownames(name.data)<-names(names.out)

mp.three <- ggplot()+ # Figure 4 subplot e
geom_point(data=df, aes(x=PC1,y=PC2), shape=21, size=point.size-1, fill='black')+
geom_abline(intercept=0,slope= 0 )+
geom_segment( aes(x=PC1,xend=PC1,y=0,yend=PC2, col= V1), data=df,lwd=line.width )+
scale_colour_gradient2(high='red',low='blue',mid='white',na.value='black',midpoint=2.25,limits=c(min(log10(masses)),max(log10(masses))) )+
labs(x=expression(paste( '',italic('m'),'' )), y=expression(paste( '',italic('D'[m]),'' )),colour=expression(paste(log[10],'(',italic(mass),')',' [g]')))+
geom_segment(aes(x=c(0,0),xend=svd3[,1], y=c(0,0),yend=svd3[,2]))+
geom_text(aes(x=svd3[,1]*1.1,y=svd3[,2]*1.1), label= c('coracoid shape','flight style'), size=genus.text, nudge_x=c(max(svd3)*0.2,max(svd3)*0.1), nudge_y=c(max(svd3)*0.01,-max(svd3)*0.015))+
geom_text(data=name.data, aes(x=V1,y=V2) , label=  labs  , nudge_x=c(max(svd3)*0,-max(svd3)*0.0,max(svd3)*0.4,max(svd3)*0.4), nudge_y=c(-max(svd3)*0.02,max(svd3)*0.04,max(svd3)*0.0,0), size=genus.text )+
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


# The subplots will now be arranged into like columns, and then placed in a larger combined plot to be saved as a .PDF file
# The package 'ggpubr' is required if it has not already been loaded:

col.1<-ggarrange(mp.one,mp.two,mp.three,nrow=3,
labels=c('a','c','e'),label.x=-0.010,label.y=1.035,font.label=list(size=20))
col.2<-ggarrange(mp,mp2,mp3,nrow=3,align='hv',
labels=c('b','d','f'),label.x=-0.010,label.y=1.035,font.label=list(size=20))
col.3<-ggarrange(mp.legend,mass.legend,nrow=2,align='hv')
# Like columns 

dev.new(width=16.9,height=16.9,unit='cm')
# New plot
ggarrange(col.1,col.2,col.3,ncol=3,widths=c(0.5,0.6,0.35))
# Arrange columns 

getwd() # Check working directory
# ggsave(filename='trunk_plot_17_01_2022.pdf')
# save the output

# Script concludes 

