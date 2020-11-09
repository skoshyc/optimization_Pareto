#Code to find the list of biofilm thickness of the reduced model
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
library(irr)#for Cohen's kappa
library(parallel) #for parallel computing
source(file='../paretofrontierhelper.R') # loads the relevant functions
path1="Insert working directory" #The working directory should contain the ordinal list of controls...
#for the original model
setwd(path1)
newg1=read.csv(file="Ordinal list of thickness_original model.csv",header = TRUE,row.names = 1)




#Relevant iDynoMiCS protocol folder
xmlfolder="../iDynoMiCS/protocol/multi_heterotroph_ordinal_reduced/Resolution1/"


for(i in 1:150){
  setwd(xmlfolder)
  
  #multi-heterotroph model
  new_inputfile=paste0("multi_heterotroph_species_","step_",i,".xml")
  xmlparser_multiheterotroph("multi_heterotroph_species.xml",newg1[i,"Oxygen_conc"],newg1[i,"Nitrate_conc"],new_inputfile)
  
  #gingivalis-gordonii model
  #new_inputfile=paste0("gingivalis_gordonii_step_",i,".xml")
  #xmlparser_gingivalis("gingivalis_gordonii_Case3.xml",newg1[i,"Protein_conc"],newg1[i,"Toxin_conc"],new_inputfile)
  }
  cl <- makeCluster(3)
  
  new_input=mixedsort(dir(xmlfolder,pattern = "multi_heterotroph_species_step_"))
  #new_input=mixedsort(dir(xmlfolder,pattern = "gingivalis_gordonii_step_"))
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
  
  #insert relevant iDynoMiCS results folder
  resultfolder="../iDynoMiCS/results/multi_heterotroph_ordinal_reduced/Resolution1/"

  
  for(i in 1:150){
  #multi-heterotroph model
  new_inputfile=paste0("multi_heterotroph_species_","step_",i,"[(]") 
    # parenthesis is a special character. 
  #For it to be recognized in a regular expression, you have to place it within [].
  #this is done so list1 gives step_1 only and not step_1, step_11, step_12 and so on. 
    
    #gingivalis-gordonii model  
  #new_inputfile=paste0("gingivalis_gordonii_","step_",i,"[(]")
  list1=dir(resultfolder,pattern = new_inputfile);#To get the multiple simulations
  
  
    k=1
    xmlpath1=paste0(resultfolder,list1[k])
    setwd(xmlpath1)
    unzip('agent_Sum.zip',exdir = paste0(getwd(),'/agent_Sum'))
    unzip('env_State.zip',exdir = paste0(getwd(),'/env_State'))
    list2=mixedsort(dir(paste0(xmlpath1,'/agent_Sum')))
    n=as.numeric(str_extract(list2[length(list2)], "[0-9]+"))
    simpsonData=simpsonIndex_noEPS(xmlpath1,n,n,xmlpath1)
    newg1[i,"Diversity_index"]=simpsonData[length(simpsonData)]
    newg1[i,5:ncol(newg1)] =speciesnumber_noEPS(xmlpath1,n,n)[2,] # to get the species count at last iteration
    l=getEnvData("last")
    newg1[i,"Mean_biofilm_thickness"]=l$heights[1]
  
}
setwd(path1)


write.csv(newg1,file="Ordinal list of thickness_(Insert filename for reduced model).csv")
###comparing the two lists using the cohen's kappa
list1=read.csv(file = "Ordinal list of thickness_original model.csv",header = TRUE,row.names = 1)

list2=read.csv(file="Ordinal list of thickness_(Insert filename for reduced model).csv",header = TRUE,row.names = 1)
comparelist=cbind.data.frame(rank(-list1$Mean_biofilm_thickness),rank(-list2$Mean_biofilm_thickness))
colnames(comparelist)=c("Ranking of original model","Ranking of reduced model")
kappa2(comparelist,"squared")


