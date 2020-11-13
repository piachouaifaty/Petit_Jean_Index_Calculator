#the following commands will make sure BiocManager and igraph are present on the user's system
#we are using ChemmineR in BiocManager to read in molfiles and sdf files into the R environment
#igraph is used for graph manipulation

petitjeanindex = function(molsdfpath)
  
{
  #Loading Libraries
  
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
  if (!requireNamespace("BiocManager", quietly=TRUE))
    BiocManager::install("ChemmineR")
  library("ChemmineR")
  if (!require(igraph)) install.packages("igraph")
  library(igraph)
  
  #the read.SDFset function reads an sdf file and stores each molecule from the sdf in a list which we have called sdfset here
  #the number of items in the sdfset indicates the number of molecules documented in the sdf file
  
  sdfset = read.SDFset(molsdfpath) 
  l=length(sdfset) #this is the number of molecules stores in the sdf file
  
  #we need to generate a petit jean index for each molecule
  #loop through the number of molecules
  #if the file is a single molfile with one molecule, it will iterate only once
  
  for (i in 1:l) {
    
    #the indeces of each molecule object inside the list are:
    #[i] corresponds to the molecule in the list of molecules
    #[i][1] corresponds to the first two blocks of the molfile
    #[i][2] corresponds to the third block, the xyz coordinates of the atoms
    #[i][3] is the block of interest, it represents the connections between atoms, and we will use it to build our graph
    bblock = sdfset[[i]][[3]] #[3] will output the bond block of [1] the first molecule in the list as a matrix
    
    #we extract the edge pairs from the bond block matrix 
    #we use transpose (t) to get the output by row rather than by column
    edgevect = c(t(bblock[,c(1,2)]))
    #the edgevect vector contains the edges of our graph so we use it as input and create a graph g
    g=graph(edgevect, n=max(edgevect), directed=FALSE)
    
    #the eccentricity function will generate a vector of the eccentricity of each atom in our graph for all vectors 
    #we specify vids=V(g) to indicate that we want the eccentricities of all atoms
    ecc = eccentricity(g, vids = V(g))
    
    #the radius is the smallest eccentricity
    rad=min(ecc)
    
    #the diameter is the largest eccentricity
    dia=max(ecc)
    
    print(sdfset[[i]][[1]]) # print the first two blocks of the molfile
    print(paste("radius =", rad)) #print the radius
    print(paste("diameter =", dia)) # print the diameter
    pj=((dia-rad)/rad) #calculate the petit jean index according to the formula
    print(paste("Petit Jean Index =", pj)) #print the petit jean index
  }}