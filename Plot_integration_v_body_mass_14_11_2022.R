# This script will produce plots demonstrating how integration with 
# avian body modules changes with increasing body mass
# These plots correspond to 'Figure 1' in Orkney & Hedrick 2023

# This script requires that previous analysis output have been produced and saved
# This can be achieved by...
# Running 'Residual_Csize_10_08_2022.R' 
# to first obtain the array of residual centroid sizes,
# for the desired selection of taxa
# Then running 'Integration_v_body_mass_11_08_2022.R'
# to compute the pairwise integration across the skeleton 
# over several intervals of body mass
# These individual scripts should be consulted for furhter details
# or customisation. 

# Set the working directory
setwd('D:/Documents/Alex_birds/Ideas_2022/SICBY')
# The path will need to be customised by the user

# Load packages that are needed to organise and plot data;
# If the user has not already installed a required package then run 
# 'install.packages('package name') and proceed as guided
library(ggplot2) # essential plotting functions
library(reshape2) # data organisation
library(ggpubr) # functions to combine plots
library(ggnewscale) # functions to overlay multiple fields of continuous data in plots
# Packages loaded

# Read the .csv file describing skeletal integration and its evolution with body mass
size.int<-read.csv('size_scale_integration_11_08_2022.csv')
# File read

# Remove the first row and set negative effect sizes to zero
size.int<-size.int[,-1]
size.int[size.int < 0 ]<-0
# Done

# Read the .csv file describing the P-values associated with skeletal integration and its evolution with body mass
size.int.p<-read.csv('size_scale_integration_p_11_08_2022.csv')
# File read

# Remove the first row and set negative effect sizes to zero
size.int.p<-size.int.p[,-1]
size.int.p[size.int.p < 0 ]<-0
# Done 

# The following lines define a vector of bones belonging to the 'head' module, 
# and retrieve the required indices to access the integration information
head.bones<-
combn(c('skull','mandible'),2) 
within.head<-matrix(NA,1,dim(head.bones)[2])
for(i in 1:dim(head.bones)[2] ){
	within.head[,i]<-intersect(grep(head.bones[1,i],colnames(size.int)),
	grep(head.bones[2,i],colnames(size.int)) )
}
within.head<-as.vector(within.head)
# Done

# The following lines define a vector of bones belonging to the 'wing' module, 
# and retrieve the required indices to access the integration information
wing.bones<-
combn(c('humerus','radius','ulna','carpometacarpus'),2) 
within.wing<-matrix(NA,1,dim(wing.bones)[2])
for(i in 1:dim(wing.bones)[2] ){
	within.wing[,i]<-intersect(grep(wing.bones[1,i],colnames(size.int)),
	grep(wing.bones[2,i],colnames(size.int)) )
}
within.wing<-as.vector(within.wing)
# Done

# The following lines define a vector of bones belonging to the 'trunk' module (body), 
# and retrieve the required indices to access the integration information
trunk.bones<-
combn(c('scapula','coracoid','sternum','synsacrum'),2) 
within.trunk<-matrix(NA,1,dim(trunk.bones)[2])
for(i in 1:dim(trunk.bones)[2] ){
	within.trunk[,i]<-intersect(grep(trunk.bones[1,i],colnames(size.int)),
	grep(trunk.bones[2,i],colnames(size.int)) )
}
within.trunk<-as.vector(within.trunk)
# Done

# The following lines define a vector of bones belonging to the 'leg' module, 
# and retrieve the required indices to access the integration information
leg.bones<-
combn(c('femur','tibiotarsus','tarsometatarsus'),2) 
within.leg<-matrix(NA,1,dim(leg.bones)[2])
for(i in 1:dim(leg.bones)[2] ){
	within.leg[,i]<-intersect(grep(leg.bones[1,i],colnames(size.int)),
	grep(leg.bones[2,i],colnames(size.int)) )
}
within.leg<-as.vector(within.leg)
# Done

# The following lines define a matrix to receive integration values binned by body
# module and body mass
inc<-unique(size.int$bin.mass)
mean.size.int<-matrix(NA,length(inc),21)
mean.size.int[,1]<-log10(inc)
colnames(mean.size.int)<-c('bin.mass',
'head','headmin','headmax','headminmin','headmaxmax',
'wing','wingmin','wingmax','wingminmin','wingmaxmax',
'trunk','trunkmin','trunkmax','trunkminmin','trunkmaxmax',
'leg','legmin','legmax','legminmin','legmaxmax')
# Done

# The commenting style will be simplified within the following loop, 
# which will produce boot-strapped probability intervals for change in integration within
# modular body regions, over body mass

for(i in 1: length(inc) ){ # For all unique mean body mass values
	ind<-which(size.int$bin.mass==inc[i]) # Find the associated indices
	mean.size.int[i,2]<-mean( size.int[ind,within.head] ) # Compute mean within-head integration

	boot<-list() # Define a dumby variable
	for(j in 1:1000){ # For 1000 replicates
		 boot[[j]]<-sample(size.int[ind,within.head],replace=T) # Sample integration with replacement
	}
	boot<-unlist(boot)
	mean.size.int[i,3]<-quantile(boot,probs=0.32) # Use the quantiles to extract probability intervals at 1 and 2 standard deviations
	mean.size.int[i,4]<-quantile(boot,probs=0.68)
	mean.size.int[i,5]<-quantile(boot,probs=0.05)
	mean.size.int[i,6]<-quantile(boot,probs=0.95)
	# Comments will not be repeated unnecessarily

	# Repeat process for within-wing integration:
	mean.size.int[i,7]<-mean( rowMeans( size.int[ind,within.wing] ))

	boot<-list()
	for(j in 1:1000){
		 boot[[j]]<-sample(rowMeans(size.int[ind,within.wing]),replace=T)
	}
	boot<-unlist(boot)
	mean.size.int[i,8]<-quantile(boot,probs=0.32)
	mean.size.int[i,9]<-quantile(boot,probs=0.68)
	mean.size.int[i,10]<-quantile(boot,probs=0.05)
	mean.size.int[i,11]<-quantile(boot,probs=0.95)

	# Repeat process for within-trunk integration:
	mean.size.int[i,12]<-mean( rowMeans( size.int[ind,within.trunk] ))

	boot<-list()
	for(j in 1:1000){
		 boot[[j]]<-sample(rowMeans(size.int[ind,within.trunk]),replace=T)
	}
	boot<-unlist(boot)
	mean.size.int[i,13]<-quantile(boot,probs=0.32)
	mean.size.int[i,14]<-quantile(boot,probs=0.68)
	mean.size.int[i,15]<-quantile(boot,probs=0.05)
	mean.size.int[i,16]<-quantile(boot,probs=0.95)

	# Repeat process for within-leg integration:
	mean.size.int[i,17]<-mean( rowMeans( size.int[ind,within.leg] ))

	boot<-list()
	for(j in 1:1000){
		 boot[[j]]<-sample( rowMeans(size.int[ind,within.leg]),replace=T)
	}
	boot<-unlist(boot)
	mean.size.int[i,18]<-quantile(boot,probs=0.32)
	mean.size.int[i,19]<-quantile(boot,probs=0.68)
	mean.size.int[i,20]<-quantile(boot,probs=0.05)
	mean.size.int[i,21]<-quantile(boot,probs=0.95)

} # Loop concludes
# The output matrix has now been populated with mean within-module integration statistics, 
# as well as 2 and 1 sigma confidence envelope statistics

# The following lines will re-organise the data for plotting
df<-as.data.frame(mean.size.int)
molten.df<-melt(df,id='bin.mass')
# Done

# Define a dataframe for the confidence envelopes
ribbons<-c(grep('min',molten.df$variable),grep('max',molten.df$variable))
# Done

# Incorporate confidence envelopes into the re-organised data
melt.df<-melt(df[-ribbons,],id='bin.mass')
# Done 

# Define a colour blind friendly palette for plotting
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Done

mp<- # This plot will be produced so that its legend can be extracted should the user desire it for their own uses
ggplot()+ # Make plot 
geom_line(aes(x=bin.mass,y=value,group=variable,linetype=variable,col=variable),size=1,data=molten.df[-ribbons,])+ # Plot lines of average within-module integration over body mass
scale_colour_manual(values=c('head'=cbbPalette[7],'wing'=cbbPalette[4],
'trunk'=cbbPalette[6],'leg'=cbbPalette[8]))+ # Define line colours of body modules 
labs(colour='module',linetype='module')+ # Label the line colours by corresponding body modules
# Sundry thematic settings follow that regulat plot appearance:
theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=15,colour='black'),
      axis.text.y=element_text(size=15,colour='black'),
	axis.title.x.bottom=element_text(size=20,colour="black"),
	axis.title.y.left=element_text(size=20,colour="black"),
      axis.ticks=element_line(size=2),
      legend.position="right",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.key.width = unit(2, "cm"))
# A line plot has been produced (the purpose of this plot is just to extract its legend)

# Extract the plot's legend
mp.legend<-cowplot::get_legend(mp)
# Make sure the cowplot package is installed. 


# The first subplot will be commented for illustrative purposes; after this further subplot code will be identified at its start, but not commented explicitly
wing.mp<- # Figure 1 subplot a; evolution of within-wing integration as a function of body mass
ggplot()+ # Make plot 
geom_ribbon(aes(x=bin.mass,ymin=wingminmin,ymax=wingmaxmax),alpha=1,data=df,bg="#80FEF3")+ # Plot the 2 sigma confidence envelope of within-wing integration with body mass
geom_ribbon(aes(x=bin.mass,ymin=wingmin,ymax=wingmax),alpha=1,data=df,bg="#40CEB3")+ # Plot the 1 sigma confidence envelope of within-wing integration with body mass
geom_line(aes(x=bin.mass,y=wing),size=1,col=cbbPalette[4], data=df,linetype='22')+ # Plot a line tracking the mean within-wing integration with body mass
geom_text(x=2.0, y=5.25, aes(label=c('Wing')), size=6 )+ # Add a text label to represent the body module being considered; the wing
lims(y=c(2.5,5.5))+ # Define axis limits
labs(colour='module',linetype='module')+ # This line is not required; it was included as a result of re-using code, but it does not have any undesireable effect
labs(x=expression(paste(log[10],'(',italic(mass),')',' [g]')), y=expression(italic(Z)))+ # Label the x and y axes of the plot
# Sundry instructions follow for the plot's thematic appearance:
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=12,colour='black'),
      axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
# Plot produced

# Figure 1 subplot a is not finished yet; the plot in the manuscript contains illustrative microplots overlaid upon it;

# The following 'yvals' define the arbitrary y-axial position of the microplots
yvals<-c(2.5,2.9,3.3,3.5,3.6,3.7)
#yvals<-c(0,3.2,0,3.5,3.6,3.7) # for Telluraves
#yvals<-c(3.2,0,0,3.2,0,3.2) # for non-Telluraves
#yvals<-c(3.1,0,3.3,3.5,3.6,3.7) # for non-Apodiformes
# Different y-axial positions may be desired if the user is investigating a bespokely chosen subset of birds

# 'kseq' is the vector of microplot identities that will be plotted, corresponding to x-axial positions
# The user may make bespoke choices about which microplots to represent, but two possibilities are presented here:
kseq<-inc[seq(1,10,by=2)]
kseq<-inc[c(1,2,3,5,7,9)]
# Done


mp.comp<-wing.mp # Make a new plot, overlaying the microplots upon the existing wing plot (Figure 1 subplot a)
for(k in 1:length(kseq)){ # For each microplot
	recompile <- matrix(NA,13,13) # The microplot will be 13 by 13 cells wide
	rownames(recompile)<-c('skull','mandible','carpometacarpus','radius','ulna','humerus',
	'sternum','coracoid','scapula','synsacrum','femur','tibiotarsus','tarsometatarsus')
	colnames(recompile)<-rownames(recompile) # The cells in the microplot correspond to average integration between pairwise combinations of bones across the bird body plan

	for(i in 1:length(size.int[1,-c(1:2)])  ){ # This loop will extract the relevant integration data from a previously produced and saved .csv file
		names<-colnames(size.int[,-c(1:2)])
		x<-  sub("\\..*","", names)[i]
		y<-  sub(".*\\.","", names)[i]
		ind<-which(size.int$bin.mass==kseq[k])
		recompile[which(rownames(recompile)==x),which(colnames(recompile)==y)]<-mean( size.int[ind,-c(1:2)][,i] )
		recompile[which(rownames(recompile)==y),which(colnames(recompile)==x)]<-mean( size.int[ind,-c(1:2)][,i] )
	} # Loop concludes

	result.rev<-apply(t(recompile),2,rev) # The resultant data for the microplot will be transposed
	longData<-melt(result.rev) # The data will be re-organised into a format suitable to supply to plotting functions
	longData$value[which(longData$value<0)]<-0 # Cells which have a value below 0 will be set to 0 
	longData$Var1<-(as.numeric(longData$Var1)/(100/(5.5-2.5)))+yvals[k] # The y-axial position of the microplot will be specified
	longData$Var2<-(as.numeric(longData$Var2)/(100/2))+log10(kseq[k]) # The x-axial position of the microplot will be specified

	mp.comp<- mp.comp+ # The overlaying operation will be performed
	coord_fixed(ratio=2/3)+ # A set coordinate ratio is defined to conserve geometry
  	geom_raster(aes(fill=value,x = Var2, y = Var1),data=longData)+ # The overlaid microplot is defined as a raster object
	scale_fill_gradientn(colours = c('yellow','red','black') ,na.value="white",lim=c(0,6))+ # Stronger integration within the microplot will be represented by a darker colour
	new_scale("fill") + # The new plotting scale of the overlaid microplot is defined as a colour fill
	theme(legend.position='none') # No legend is required for the new scale; it is for illustrative purposes alone and does not require explicit quantification
}# Loop condludes
# Microplots of matrices of pairwise integration at different
# representative ranges of body mass have been overlaid on Figure 1 subplot a 


head.mp<- # Figure 1 subplot d will be produced; it represents the evolution of within-head integration as a function of body mass
ggplot()+
coord_fixed(ratio=(2/3)/(1.75/3))+
geom_ribbon(aes(x=bin.mass,ymin=headminmin,ymax=headmaxmax),alpha=1,data=df,bg="#F5CE80")+
geom_ribbon(aes(x=bin.mass,ymin=headmin,ymax=headmax),alpha=1,data=df,bg="#F58E40")+
geom_line(aes(x=bin.mass,y=head),size=1,col=cbbPalette[7], data=df,linetype='solid')+
lims(y=c(4.25,6.0))+ 
geom_text(x=2.0, y=5.7, aes(label=c('Head')), size=6 )+
labs(colour='module',linetype='module')+
labs(x=expression(paste(log[10],'(',italic(mass),')',' [g]')), y=expression(italic(Z)))+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=12,colour='black'),
      axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

trunk.mp<- # Figure 1 subplot b will be produced; it represents the evolution of within-trunk integration as a function of body mass
ggplot()+
coord_fixed(ratio=(2/3)/(3.25/3))+
geom_ribbon(aes(x=bin.mass,ymin=trunkminmin,ymax=trunkmaxmax),alpha=1,data=df,bg="#80E2F2")+
geom_ribbon(aes(x=bin.mass,ymin=trunkmin,ymax=trunkmax),alpha=1,data=df,bg="#40A2F2")+
geom_line(aes(x=bin.mass,y=trunk),size=1,col=cbbPalette[6], data=df,linetype='31')+
lims(y=c(0.75,4))+
geom_text(x=2.0, y=3.45, aes(label=c('Trunk')), size=6 )+
labs(colour='module',linetype='module')+
labs(x=expression(paste(log[10],'(',italic(mass),')',' [g]')), y=expression(italic(Z)))+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=12,colour='black'),
      axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

leg.mp<- # Figure 1 subplot e will be produced; it represents the evolution of within-leg integration as a function of body mass
ggplot()+
coord_fixed(ratio=(2/3)/(2.25/3))+
geom_ribbon(aes(x=bin.mass,ymin=legminmin,ymax=legmaxmax),alpha=1,data=df,bg="#FCE9F7")+
geom_ribbon(aes(x=bin.mass,ymin=legmin,ymax=legmax),alpha=1,data=df,bg="#FCA9E7")+
geom_line(aes(x=bin.mass,y=leg),size=1,col=cbbPalette[8], data=df,linetype='33')+
lims(y=c(2.5,4.75))+
geom_text(x=2.0, y=4.4, aes(label=c('Leg')), size=6 )+
labs(colour='module',linetype='module')+
labs(x=expression(paste(log[10],'(',italic(mass),')',' [g]')), y=expression(italic(Z)))+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=12,colour='black'),
      axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())


# The following lines of code will compute Ordinary Least Squares fits for change in pairwise-Z scores
# of integration across all combinations of skeletal elements as a function of mean body mass
# within each bin
# This will produce Figure 1 subplot c, which is a plot produced for illustrative purposes; be aware that the bins of body mass
# may overlap, causing statistical non-independence among the data

# These lines define a matrix to receive the slopes of OLS fits:
recompile <- matrix(NA,13,13)
rownames(recompile)<-c('skull','mandible','carpometacarpus','radius','ulna','humerus',
'sternum','coracoid','scapula','synsacrum','femur','tibiotarsus','tarsometatarsus')
colnames(recompile)<-rownames(recompile)
# Done

# The following loop will populate the matrix with OLS slope fits:
inc<-unique(size.int[,1])
for(i in 1:length(size.int[1,-c(1:2)])  ){
	names<-colnames(size.int[,-c(1:2)])
	x<-  sub("\\..*","", names)[i]
	y<-  sub(".*\\.","", names)[i]
	ind<-which(size.int$bin.mass==kseq[k])

	temp<-list()
	for(j in 1:length(inc)){
		ind<-which(size.int[,1]==inc[j])
		temp[[j]]<-mean(size.int[,-c(1:2)][,i][ind])
		
	}
	temp<-unlist(temp)

	coef<-lm(temp~log10(inc))

	if(summary(coef)$coef[2,4] <1){ # The user may exclude cells above a specific significance threshold if they wish, by replacing '<1'
		recompile[which(rownames(recompile)==x),which(colnames(recompile)==y)]<-coef$coef[2]
		recompile[which(rownames(recompile)==y),which(colnames(recompile)==x)]<-coef$coef[2]
	} 
}
# Done

# Re-organise the output for plotting
result.rev<-apply(t(recompile),2,rev)
longData<-melt(result.rev)
# Done

# Define colours for different body modules of the avian skeleton 
bone_colours<-rev(c(rep(cbbPalette[7],2),rep(cbbPalette[4],4),rep(cbbPalette[6],4),rep(cbbPalette[8],3)))
# Done

mp.change<- # Produce Figure 1 subplot c
ggplot(longData, aes(x = Var2, y = Var1)) + # Make plot, specify data sources and axial variables
	coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+ # Define geometry 
  geom_tile(aes(fill=value)) + # Treat the plot as a tile plot 
   scale_fill_gradient2(low='blue',mid='white',high='red' ,na.value="white",midpoint=0) + # Fill the tiles with colour corresponding to the slope of OLS fit
	labs(fill = "slope")+ # Label the legend
  labs(x="", y="") + # No axial labels are needed
# Sundry thematic appearance instructions follow
  theme_bw() + theme(axis.text.x=element_text(size=8, angle=90, vjust=0.3,colour=rev(bone_colours),face = "bold"),
                     axis.text.y=element_text(size=8,colour=(bone_colours),face = "bold"),
				legend.key.width=unit(1/2,'cm'),
				legend.key.height=unit(1,'cm'), 
				legend.title=element_text(size=15),
				legend.text=element_text(size=12),
				plot.margin=unit(c(-0.5,-1,0.5,-1),'cm'),
                     plot.title=element_text(size=11,hjust=.5), legend.position="right",panel.background = element_blank())+
	labs(fill = expression(paste(Delta,italic('Z') )))+
geom_hline(yintercept = 3.5,size=1,colour="black")+
geom_hline(yintercept = 11.5,size=1,colour="black")+
geom_hline(yintercept = 7.5,size=1,colour="black")+
geom_vline(xintercept = 10.5,size=1,colour="black")+
geom_line(data=as.data.frame(cbind(c(2.5,2.5),c(0.5,13.5))),aes(x=V1,y=V2),size=1,colour='black')+
geom_line(data=as.data.frame(cbind(c(6.5,6.5),c(0.5,13.5))),aes(x=V1,y=V2),size=1,colour='black')+
geom_line(data=as.data.frame(cbind(c(10.5,10.5),c(0.5,13.5))),aes(x=V1,y=V2),size=1,colour='black')+
scale_x_discrete(position = "top") 
# Plot complete
# This has produced a plot that shows general trends in pairwise-Z of integration across the skeleton as a function
# of increasing body mass

# The following lines will assemble all the subplots into one combined plot 
mp.modules<-ggarrange(head.mp,trunk.mp,leg.mp,ncol=1,nrow=3,align='hv',
labels=c('b','d','e'),label.x=-0.025,label.y=1,font.label=list(size=25))
# Subplots b, d and e have been combined and aligned
mp.left<-ggarrange(mp.comp,mp.change,ncol=1,
labels=c('a','c'),label.x=0.075,label.y=1,font.label=list(size=25))
# Subplots a and c have been combined and aligned
dev.new(height=16,width=18,unit='cm')
# A new overall plot space has been created
ggarrange(mp.left,mp.modules,ncol=3,widths=c(1,0.8,0.05))
# All plots have been added to this overall space

# An alternative arrangement: 
mp.modules<-ggarrange(head.mp,trunk.mp,leg.mp,ncol=1,nrow=3,align='hv',
labels=c('b','d','e'),label.x=-0.025,label.y=1,font.label=list(size=25))
mp.left<-ggarrange(mp.comp,mp.change,ncol=1,
labels=c('a','c'),label.x=0.15,label.y=1,font.label=list(size=25))
dev.new(height=22,width=18,unit='cm')
ggarrange(mp.left,'',mp.modules,ncol=4,widths=c(1,0.1,0.8,0.05))


# Alternatively, to arrange these plots to leave a gap for an illustration,
# and save the plot as a large PNG file which can be used as an illustration
# background... 

# Define a left column of subplots
mp.left<-ggarrange(mp.comp,trunk.mp,ncol=1,nrow=2,
labels=c('a','b'),label.x=0,label.y=1,font.label=list(size=25))
blank<-ggplot()+ geom_blank()+
  theme(
	axis.text.x=element_text(size=12,colour='black'),
  	axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
# Done

# Define a central column of subplots
mp.centre<-ggarrange(blank,mp.change,ncol=1,nrow=2,
labels=c('','c'),label.x=0,label.y=1,font.label=list(size=25))
# Done

# Define a right column of subplots
mp.right<-ggarrange(head.mp,leg.mp,ncol=1,nrow=2,
labels=c('d','e'),label.x=0,label.y=1,font.label=list(size=25))
# Done

# Produce and populate a new plot space
dev.new(height=10,width=16,unit='cm')
ggarrange(mp.left,'',mp.centre,'',mp.right,ncol=6,widths=c(1,0.075,1.2,0.075,1,0.05))
# Done

# Save the plot as a pdf and as a high quality .PNG file
ggsave(filename='mass_v_integration_26_11_2022b.pdf')
ggsave('mass_v_integration_26_11_2022b.png',dpi=300)
# Done

# Script concludes