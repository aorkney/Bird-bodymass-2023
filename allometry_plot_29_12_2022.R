# This script will produce Figure 5 from Orkney & Hedrick 2023


# A plot of the residual log carpometacarpus centroid size
# as a function of mass, 
# with flap-gliding and migrating birds identified in the point colour and shape, will be produced...
# Outliers will be excluded so that the allometric model is tightly fit...

# Call relevant packages
library(phytools)
library(geomorph)
library(nlme)
# Done

# Set work directory
setwd('D:/Documents/Alex_birds/Ideas_2022/SICBY')
# The user must modify as appropriate

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

# Sundry function defined to interrogate data objects
get.item <- function( X , item ) { X[[ item ]] }
# This is a simple indexing function


# Define a function to compute residual centroid size from a allometric model fir
get.residual.Csize <- function( array, masses, phylogeny, taxa, bones ){ # Shape data array, mass vector, phylogeny, taxa and bones of interest required
	allometry.Csize <- list() # Dumby variable to receive allometric models
	newphy <- drop.tip( phylogeny , phylogeny$tip.label[ !phylogeny$tip.label %in% taxa ] ) # Prune phylogeny to taxa of interest
	for(i in 1:length(array[bones]) ){ # For each bone 
		df <- data.frame( mass= log10(masses[taxa]), Csize = log10( array[[ bones[i] ]][taxa] ), species = taxa ) # Define a dataframe
		allometry.Csize[[ i ]] <- gls( Csize ~ mass, correlation = corPagel( 0, phy=newphy, form= ~ species  ), data=df ) # Compute an allometric model
	}
	names(allometry.Csize) <- names(array[bones]) # Ensure allometric model names reflect bones of interest
	residual_Csize <- lapply( lapply( allometry.Csize , get.item , item = "residuals" ) , c ) # Extract model residuals
	return(residual_Csize) # Return residuals
} # Function concluides
# This function uses phylogenetic Generalised Least Squares to compute
# the allometrically adjusted centroid sizes of skeletal elements...
# The user must provide the array of original bird skeletal element centroid sizes, bodymasses,
# phylogeny, bones and taxa of interest...

# Define a vector of bird taxa that are not unusual outliers, such as penguins
Outliers<-c('Malurus_splendens','Neodrepanis_coruscans','Coracias_cyanogaster','Spheniscus_humboldti')
nOutliers <- names(masses)[-match(Outliers,names(masses))]
# Done

# Compute residual centroid size of carpometacarpus
rCsize<-get.residual.Csize(GPA.Csize, masses, pruned.tree, nOutliers, 'carpometacarpus')
# Done 

# Load flight style data matrix (binary scores)
flight<- read.csv('flight_masses_22_10_2022_plus_A.csv')
flight.rownames<-flight$X
flight<-flight[,-c(1,2,14)]
rownames(flight)<-flight.rownames
# Done

# Extract migrants and flap-gliders
migrant<-flight[nOutliers,]$Migratory_flight
flap_glider<-flight[nOutliers,]$Flap_gliding
# Done

# Create categorisation of different groups, based on these flight styles
flap_glider[which(flap_glider==0 & migrant==0)]<-'neither'
flap_glider[which(flap_glider==1 & migrant==1)]<-'both'
migrant[which(flap_glider=='both' & migrant==1)]<-''
flap_glider[which(flap_glider==1)]<-'flap_glider'
migrant[which(migrant==1)]<-'migrant'
migrant[which(migrant==0)]<-''
flap_glider[which(flap_glider==0)]<-''
# Done

# Create categorical variable for flight style
flight_style<- paste(migrant,flap_glider,sep='')
# Done

# Organise the data into a dataframe
df<- as.data.frame(cbind(
log10(masses[nOutliers]),
rCsize$carpometacarpus,
flight_style
))
# Done

# Ensure data is correct class
df$V1<-as.numeric(as.character(df$V1))
df$V2<-as.numeric(as.character(df$V2))
# Done


# Plotting package to be loaded
library(ggplot2)
# Done

# Sundry plotting variables
 axis.title <- 15
 axis.text <- 10
# Defined


# Create a vector in order to determine the x-axial order of violin plots of carpometacarpus residual centroid size as a function of flight style group
order_by_mean<-list()
for(i in 1:length(unique(df$flight_style))){
	order_by_mean[[i]]<- mean((as.numeric(df$V2[grep(pattern=unique(df$flight_style)[i],df$flight_style)])),na.rm=T)
}
order_by_mean<-(unique(df$flight_style)[order(unlist(order_by_mean))])
# Done

# Factorise the flight style variable in the data frame
df$flight_style<-factor(df$flight_style)
# Done

plot1<- # Produce Figure 5 subplot a 
ggplot()+ # Make plot 
geom_point(aes(x=V1, y=V2, colour=factor(flight_style, levels= order_by_mean), shape=factor(flight_style, levels= order_by_mean) ), data=df, size=4, stroke=2)+ # Scatter plot of mass and residual centroid size
labs(y=expression(Carpometacarpus["rCsize"]),x=expression(paste(log[10],'(',italic(mass),')',' [g]')),colour='flight style',shape='flight style')+ # Labels and legend 
# Thematic instructions for plot appearance follow:
theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=axis.text,colour='black'),
      axis.text.y=element_text(size=axis.text,colour='black'),
	axis.title.x.bottom=element_text(size=axis.title,colour="black"),
	axis.title.y.left=element_text(size=axis.title,colour="black"),
      axis.ticks=element_line(size=2),
      legend.position="bottom",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.key.width = unit(1, "cm"))
# Plot concludes

plot2<- # Produce Figure 5 subplot b 
ggplot()+ # Make plot 
geom_violin(aes( x=factor(flight_style, levels= order_by_mean), y=V2, group=factor(flight_style, levels= order_by_mean), col=factor(flight_style, levels= order_by_mean)), data=df, trim=F, lwd=2)+ # Produce violin plots
geom_boxplot(width=1/3,aes( x=factor(flight_style, levels= order_by_mean), y=V2, group=factor(flight_style, levels= order_by_mean)), data=df)+ # Overlay boxplots
labs(y=expression(Carpometacarpus["rCsize"]))+ # Label y axis 
# Add summary statistics onto plot: a point to represent the mean for each violin plot
stat_summary(aes( x=factor(flight_style, levels= order_by_mean), y=V2, group=factor(flight_style, levels= order_by_mean)), data=df, fun = "mean",
               geom = "point",
               color = "red")+
# Done 
# Thematic instructions for plot appearance follow:
theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=axis.text,colour='black',angle=90),
      axis.text.y=element_text(size=axis.text,colour='black'),
	axis.title.x.bottom=element_blank(),
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

# Load package to combine subplots
library(ggpubr)
# Done

# Define new plot space
dev.new(width=6,height=9,unit='cm')
# Done 

# Arrange subplots and save
ggarrange(plot1,plot2,nrow=2, labels=c('a','b'), font.label=list(size=30))
ggsave(filename='allometry_flight_style_29_12_2022.pdf')
# Done

# Script concludes 


