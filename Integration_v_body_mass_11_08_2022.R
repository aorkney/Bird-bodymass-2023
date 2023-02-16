# This script is required to perform the analyses that underlie 'Figure 1' in Orkney & Hedrick 2023

# This script will produce '.csv' files for 10 different bins of body mass, 
# showing how integration among pairwise combinations of bones changes with mass
# For the studied birds (149 species), these bins will be overlapping, 
# which will generate statistical non-independence between some bins
# This can be regarded as being similar to a rolling-average
# Because the birds are not uniformly distributed in frequency over log10(bodymass), 
# some bins will have wider or narrower ranges of mass represented than others

# The user will need to run 'Residual_Csize_10_08_2022.R' 
# to first obtain the array of residual centroid sizes used by this script

# Load package 'abind'; the 'array bind' function will be needed;
# If the user has not installed this package run 'install.packages('abind')' and proceed as guided
# no further comments advising package installation will be provided
library(abind)
# Package loaded


# Retrieve the names of the skeletal elements to be considered
bones<-names(residual.Csize)
# Names defined


# The following lines specify a function 
# that computes pairwise Z-scores (effect sizes) of integration between 
# different combinations of bones, for a stated group of taxa and a provided phylogeny relating the taxa

pair.int.phy <- function( input, phylogeny, taxa ){ # Define function names and required arguments; landmark coordinates (input), phylogeny and taxa

	newphy<-drop.tip( phylogeny, phylogeny$tip.label[ !phylogeny$tip.label  %in% taxa] ) # Prune the source phylogeny as appropriate
	
	bones <- names(input) # Define the bone names in the input data
	output <- matrix(NA,length(bones),length(bones)) # Define a matrix to receive pairwise combinations of Z-scores between different bones
	output.p <- matrix(NA,length(bones),length(bones)) # Define  matrix to receive corresponding P-values
	rownames(output)<-bones # Name the rows and columns appropriately after the bones they represent
	colnames(output)<-rownames(output)
	colnames(output.p)<-colnames(output)
	rownames(output.p)<-rownames(output)
	bone_combinations<- combn(bones,2,simplify=T) # Define a vector of all possible combinations of bones that need to be considered
	
	int_list<-list() # Define a dumby list to receive integration estimates 
	for(i in 1:dim(bone_combinations)[2] ){	# For each pairwise combination 
		x<-paste(bone_combinations[,i][1]) 
		y<-paste(bone_combinations[,i][2])
		int_list[[i]] <- phylo.integration(as.matrix(input[[which(bones==x)]][taxa]),as.matrix(input[[which(bones==y)]][taxa]),phy=newphy,iter=999,print.progress=F) # Compute the Z-score and store it
	}

	comparison<-compare.pls(int_list) # Compute the difference of integration statistics between the pairwise Z-scores
	
	z.score<-unlist(comparison$sample.z) # Turn the resultant effect sizes of difference of integration into a vector

	p.val<-unlist( lapply(int_list,function(x) x$P.value ) ) # Extract the corresponding p-values 

	for(i in 1:dim(bone_combinations)[2] ){ # For each indiviual pairwise combination of bones
		x<-paste(bone_combinations[,i][1]) 
		y<-paste(bone_combinations[,i][2])
		output[which(rownames(output)==x),which(colnames(output)==y)]<-z.score[i] # Assign the relevant Z-score to the output matrix
		output[which(rownames(output)==y),which(colnames(output)==x)]<-z.score[i] # "
	}

	for(i in 1:dim(bone_combinations)[2] ){ # For each indiviual pairwise combination of bones
		x<-paste(bone_combinations[,i][1])
		y<-paste(bone_combinations[,i][2])
		output.p[which(rownames(output.p)==x),which(colnames(output.p)==y)]<-p.val[i] # Assign the corresponding P-value to the output matrix
		output.p[which(rownames(output.p)==y),which(colnames(output.p)==x)]<-p.val[i] # "
	}
		
	return(abind( output, output.p,along=3 )) # Return the Z-scores and corresponding p-values as a bound array
} # Function concludes



# Define the upper-bound of the body mass bins to be studied
breaks<-round(seq(40,length(masses),length.out=10))
# Bounds defined

# Define the width of the body mass bins to be studied, by number of taxa
bin.width<-40
# Width defined

# Define the number of taxa to be sampled in any individual replicate within a bin
n<-30 
subsample.size<-30
# Each replicate will be computed by subsampling 30 taxa within each bin

# Extract the number of pairwise combinations of skeletal elements
k <- dim(combn(names(residual.Csize),2,simplify=2))[2]
# Done

# Define matrices to receive the Z-scores, subsessetted by binned body mass, and associated p-values
size.scale.integration<-matrix( NA, length(breaks)*n, k+4 ) 
size.scale.integration.p<-matrix( NA, length(breaks)*n, k+4 ) 
# Matrices defined

# Define a vector that corresponds to the breaks between the body mass bins
break.vector<-rep(breaks,each=n)
# Done

# Comments will follow each line directly within the next loop, unless the function of plural lines is described, in which case the comment will be underneath the appropriate lines

for( i in 1 : dim(size.scale.integration)[1] ){ # For each replicate, within each body mass bin...

	end<-break.vector[i]
	start<-end-bin.width
	# Define the bin of body mass.

	bin.mass <- mean(masses[order(masses)][start:end])
	bin.birds <- names(masses[order(masses)][start:end])
	# Compute the mean mass and name the birds within the bin 

	selection<-sample(bin.birds,subsample.size)# Extract a subsample of taxa.

	subsample.mass <- mean(masses[order(masses)][start:end][selection])
	min.mass <- min(masses[order(masses)][start:end][selection])
	max.mass <- max(masses[order(masses)][start:end][selection])
	# Compute the mean, min and max body masses within the subsample of taxa. 

	size.scale.integration[i,1]<-bin.mass
	size.scale.integration[i,2]<-subsample.mass
	size.scale.integration.p[i,1]<-bin.mass
	size.scale.integration.p[i,2]<-subsample.mass
	size.scale.integration[i,3]<-min.mass
	size.scale.integration[i,4]<-max.mass
	size.scale.integration.p[i,3]<-min.mass
	size.scale.integration.p[i,4]<-max.mass
	# Store these information to the output matrix 

	result<-pair.int.phy( residual.Csize, pruned.tree, selection)# Compute pairwise integration between all combinations of skeletal elements

	size.scale.integration[i,-c(1:4)] <- result[,,1][upper.tri(result[,,1])]
	size.scale.integration.p[i,-c(1:4)] <- result[,,2][upper.tri(result[,,2])]
	# Store the resultant test statistics in the data output matrix

	print(round(i/dim(size.scale.integration)[1]*100,digit=2))# Print the progress to the screen (this code will take some time to run)
}
# This loop has produced a matrix of Z-scores and p-values that change with increasing body mass of the taxanomic subsets

# Define the index for the upper triangle of pairwise integration scores
ind<- which( upper.tri(result[,,1],diag=F) , arr.ind = TRUE )
# The lower triangle is identical, so is not needed 

# Define column names for the data output matrix
colpaste<-paste(bones[ind[,1]],bones[ind[,2]],sep='.')
# Done

# Set the column names for the objects recording Z-scores and p-values
colnames(size.scale.integration) <- c('bin.mass','subsample.mass','min.mass','max.mass',colpaste)
colnames(size.scale.integration.p) <- c('bin.mass','subsample.mass','min.mass','max.mass',colpaste)
# Assigned



# Check the working directory is as expected
getwd()
# User action may be warranted; this is where the files are to be saved 


# Save the output for subsequent analyses
write.csv( x=size.scale.integration, file='size_scale_integration_11_08_2022.csv')
write.csv( x=size.scale.integration.p, file='size_scale_integration_p_11_08_2022.csv')
# argument 'x' describes the object, 'file' is the resultant file name

# Script concludes; the saved data may be interrogated or plotted by other scripts
