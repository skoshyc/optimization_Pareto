#Code to populate the ordinal list of biofilm thickness
#will use a latin hypercube sample 
library(rPref) 
library(dplyr) 
library(igraph) 
library(ggplot2)
#processing of the ouput files using iDynoR package
library(vegan)
library(XML)
library(gtools)
library(xml2)
library(iDynoR)
library(stringr)
library(dplyr)
library(lhs)#for the latin hypercube sampling
library(parallel) #for parallel computing

source(file='../paretofrontierhelper.R') # loads the relevant functions
path1 = "Insert working directory"
setwd

#Latin hypercube sampling
#Draws a Latin Hypercube Sample from a set of unit uniform distributions for use in creating a Latin
#Hypercube Design. This sample is taken in a random manner without regard to optimization.
set.seed(1234) # we set the seed so as to get the same results
X <- randomLHS(150, 2) #first input is the number of points and second input is the number of variables

#Suppose we want the first set of points from the interval (a,b) then we would take
#newX[,1] = a + (b-a) * X[,1]

#for the multi-heterotroph model
 newg1=matrix(0,150,3,dimnames = list(c(),c("Oxygen_conc","Nitrate_conc","Mean_biofilm_thickness")));
 newg1[,1] = 0.1 + 1 * X[,1] #Range for oxygen substration (0,1)
 newg1[,2] = 0.1 + 1 * X[,2] #Range for nitrate concentration (0,1)

#gingivalis-gordonii model
#newg1=matrix(0,150,3,dimnames = list(c(),c("Protein_conc","Toxin_conc","Mean_biofilm_thickness")));
#newg1[,1] = 0.1 + 20 * X[,1] #Range for protein substration (0.1,20)
#newg1[,2] = 0.1 + 5 * X[,2] #Range for toxin concentration (0.1,5)


#insert relevant iDynoMiCS protocol folder for xmlfolder
xmlfolder="../iDynoMiCS/protocol/multi_heterotroph_ordinal/"

for(i in 1:nrow(newg1)){
  setwd(xmlfolder)
  
  #for the multi_heterotroph model
  new_inputfile=paste0("multi_heterotroph_species_","step_",i,".xml")
  xmlparser_multiheterotroph("multi_heterotroph_species.xml",newg1[i,"Oxygen_conc"],newg1[i,"Nitrate_conc"],new_inputfile)
  
  #for the gingivalis-gordonii model
  #new_inputfile=paste0("gingivalis_gordonii_step_",i,".xml")
  #xmlparser_gingivalis("gingivalis_gordonii_Case3.xml",newg1[i,"Protein_conc"],newg1[i,"Toxin_conc"],new_inputfile)
  }  
  
cl <- makeCluster(3)


new_input=mixedsort(dir(xmlfolder,pattern = "multi_heterotroph_species_step_"))

len=1

idx=seq(1,150)

# run 3 instances of the code at the same time
while(len<=(length(idx))){
  if(len==(length(idx)-1)){
    parSapply(cl , new_input[len:(len+1)], runpytho1, xmlfolder, 1)
    len=len+2;
  }
  else if(len==length(idx)){
    parSapply(cl , new_input[len], runpytho1, xmlfolder, 1)
    len=len+1;
  }
  else{
    parSapply(cl , new_input[len:(len+2)], runpytho1, xmlfolder, 1)
    len=len+3;
  }
  
}

# close cluster object
stopCluster(cl)  

#insert relevant iDynoMiCS result folder  
resultfolder="../iDynoMiCS/results/multi_heterotroph_ordinal/"

  diversity_index=matrix(0,150,1,dimnames = list(c(),c("Diversity_index")))
  speciesabud=matrix(0,150,3) #number of columns the number of species
  
  for(i in 1:150){
    #multi-heterotroph model
    new_inputfile=paste0("multi_heterotroph_species_step_",i,"[(]")
     # parenthesis is a special character. 
    #For it to be recognized in a regular expression, you have to place it within [].
    #this is done so list1 gives step_1 only and not step_1, step_11, step_12 and so on. 
    
    #gingivalis-gordonii model
    #new_inputfile=paste0("gingivalis_gordonii_","step_",i,"[(]")
    
    list1=dir(resultfolder,pattern = new_inputfile);#To get the multiple simulations
    
    #for(k in 1:length(list1)){
    k=1
    xmlpath1=paste0(resultfolder,list1[k])
    setwd(xmlpath1)
    unzip('agent_Sum.zip',exdir = paste0(getwd(),'/agent_Sum'))
    unzip('env_State.zip',exdir = paste0(getwd(),'/env_State'))
    list2=mixedsort(dir(paste0(xmlpath1,'/agent_Sum')))
    n=as.numeric(str_extract(list2[length(list2)], "[0-9]+"))
    simpsonData=simpsonIndex_noEPS(xmlpath1,n,n,xmlpath1)
    
    diversity_index[i,1]=simpsonData[length(simpsonData)]
    speciesabud[i,] =as.matrix(speciesnumber_noEPS(xmlpath1,n,n)[2,]) # to get the species count at last iteration
    colnames(speciesabud)=colnames(speciesnumber_noEPS(xmlpath1,n,n))
    
    l=getEnvData("last") #l will be a named list of the environment data from the last simulation
    # heights is an array with the mean, stddev, and max interface biofilm thickness data
    newg1[i,"Mean_biofilm_thickness"]=l$heights[1]
    #}
  }

  #newg1=cbind.data.frame(newg1,diversity_index,speciesabud)
  

setwd(path1)
write.csv(newg1,file="Ordinal list of thickness_original model.csv")
