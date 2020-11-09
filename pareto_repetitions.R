####Pareto optimization algorithm
#### pareto frontier of the idynomics model
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
library(nsga2R) #to apply the NSGA2 (non-dominated sorting genetic algorithm)
library(parallel) #for parallel computing

source(file='../paretofrontierhelper.R') # loads the relevant functions
path1="Insert working directory" #The working directory should contain the ordinal list of controls
setwd(path1)
set.seed(1234)
gen_max=50#maximum number of generations
popSize=20 #population size
tour_size=2 #tournament size
varNo=2 # number of decision variables 

objDim=2 #the number of quantities to be optimized
#for the multi-heterotroph model the bounds of oxygen and nitrate concentrations was (0.1,1)

#for the gingivalis-gordonii model bounds of protein and toxin are (0.1,20) and (0.1,5) resp.
#we have 2 variables, hence bounds are of size 2.
lowerbounds=c(0.1,0.1) #lower bound for solute concentrations 
upperbounds=c(1,1) # upper bound

multiples1=1 #the number of simulations to be run for each input xml file
mutprob=0.2 #mutation probability 

###the data should have the variables as the first set of columns, then the quantities to be optimized
#third column has to be diversity index/species count, and 4th mean biofilm thickness



 d1=read.csv("Insert ordinal list filename.csv",header=TRUE,as.is=TRUE
             ,row.names = 1)

#the procedure is based on https://rdrr.io/cran/nsga2R/src/R/nsga2R.R
parent=d1[sample(nrow(d1),popSize),] #initial population




#we are sorting the 3rd and 4th column resp
#the default is minimization, if we want maximization we add a minus sign. min(-f(x))=max(f(x))
ranking <- fastNonDominatedSorting(parent[,(varNo+1):(varNo+objDim)]);
# Rank index for each chromosome
rnkIndex <- integer(popSize);
i <- 1;
while (i <= length(ranking)) {
  rnkIndex[ranking[[i]]] <- i;
  i <- i + 1;
} 
parent <- cbind(parent,rnkIndex);


#crowding distance calculation

objRange <- apply(parent[,(varNo+1):(varNo+objDim)], 2, max) - apply(parent[,(varNo+1):(varNo+objDim)], 2, min);
cd <- crowdingDist4frnt(parent,ranking,objRange);
parent <- cbind(parent,apply(cd,1,sum));
colnames(parent)[ncol(parent)]="crowdingDist"
setwd(path1)
write.csv(parent,file=paste("Parent_generation_",1,".csv",sep=""))

for(gen in 1:gen_max){
  
  child=boundedPolyMutation(parent[,1:varNo],lowerBounds=lowerbounds,upperBounds = upperbounds,mutprob,mum=10) #mutate the oxygen and nitrate concentrations
  
  child=cbind(child,parent[,(varNo+1):(varNo+objDim)]) #to add the quantities to be optimized
  colnames(child)[(varNo+1):(varNo+objDim)]=colnames(parent)[(varNo+1):(varNo+objDim)]
  #-1 is added so we get the rows that are different from parent
  #get the mutations so as to run the code only for the different values.
  child[which(child[,1]!=parent[,1]),(varNo+1):(varNo+objDim)]=-1
  child[which(child[,2]!=parent[,2]),(varNo+1):(varNo+objDim)]=-1 
  idx=which(child[,(varNo+1)]==-1)
  
  #write xml inputs only for the mutations, run the simulations and obtain new values of objective function
  #insert relevant iDynoMiCS protocol folder
  xmlfolder="../iDynoMiCS/protocol/multi_heterotroph_Pareto/"
  
  
  #make sure to strip comments from the xml file which you wish to change
  
  idx_repeats=rep(idx,each=10)
  
  for(i in 1:length(idx_repeats)){
    setwd(xmlfolder)
    
    #for the multi-heterotroph model
    new_inputfile=paste0("multi_heterotroph_species_gen_",gen,"_step_",i,".xml")
    xmlparser("multi_heterotroph_species.xml",
              child$Oxygen_conc[idx_repeats[i]],child$Nitrate_conc[idx_repeats[i]],new_inputfile)
    
    #for the gingivalis-gordonii model
    # new_inputfile=paste0("gingivalis_gordonii_gen_",gen,"_step_",i,".xml")
    # xmlparser_gingivalis("gingivalis_gordonii_Case3.xml",
    #     child$Protein_conc[idx_repeats[i]],child$Toxin_conc[idx_repeats[i]],new_inputfile)
  }
  # create cluster object
  cl <- makeCluster(3)
  
  new_input=mixedsort(dir(xmlfolder,pattern = paste("gen_",as.character(gen),sep="")))
  #So as to get the correct generation
  len=1
  # run 3 instances of the code at the same time
  while(len<=(length(idx_repeats))){
    if(len==(length(idx_repeats)-1)){
      parSapply(cl , new_input[len:(len+1)], runpytho1, xmlfolder, 1)
      len=len+2;
    }
    else if(len==length(idx_repeats)){
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
  child_repeats=matrix(0,length(idx_repeats),(ncol(child)+1),dimnames = list(c(),c(colnames(child),"repeats")))
  child_repeats[,ncol(child_repeats)]=idx_repeats
  for(i in 1:length(idx_repeats)){
    #insert relevant iDynoMiCS results folder
    resultfolder="../iDynoMiCS/results/multi_heterotroph_Pareto/"
    
    #for the multi-heterotroph model
    new_inputfile=paste0("multi_heterotroph_species_gen_",gen,"_step_",i,"[(]")
    
    #for the gingivalis-gordonii model
    #new_inputfile=paste0("gingivalis_gordonii_gen_",gen,"_step_",i,"[(]")
    
    list1=dir(resultfolder,pattern = new_inputfile);
    
    
    k=1  
    xmlpath1=paste0(resultfolder,list1[k])
      setwd(xmlpath1)
      unzip('agent_Sum.zip',exdir = paste0(getwd(),'/agent_Sum'))
      unzip('env_State.zip',exdir = paste0(getwd(),'/env_State'))
      list2=mixedsort(dir(paste0(xmlpath1,'/agent_Sum')))
      n=as.numeric(str_extract(list2[length(list2)], "[0-9]+"))
      simpsonData=simpsonIndex_noEPS(xmlpath1,n,n,xmlpath1) #for diversity index_multi-heterotroph model
      #speciesabud =matrix(speciesnumber_noEPS(xmlpath1,n,n)[2,],nrow=1,ncol=2) 
      # to get the species count at last iteration for gingivalis-gordonii model
      colnames(speciesabud)=colnames(speciesnumber_noEPS(xmlpath1,n,n))
      child_repeats[i,(varNo+1)]=speciesabud[which(colnames(speciesabud)==colnames(child_repeats)[(varNo+1)])]
      l=getEnvData(n)
      child_repeats[i,(varNo+objDim)]=l$heights[1]
     
    
  }
  
  child_thickness=data_summary(as.data.frame(child_repeats),varname = "Mean_biofilm_thickness",groupnames = c("repeats"))
  child_species=data_summary(as.data.frame(child_repeats),varname = "Diversity_index",groupnames = c("repeats"))
  #child_species=data_summary(as.data.frame(child_repeats),varname = "Mygordonii",groupnames = c("repeats"))
  for (i in 1:length(idx)) {
    rep_row=idx[i]
    child[rep_row,(varNo+1)]=floor(child_species[i,2])#the 3rd column will give the mean diversity index/
    #species number for gingivalis-gordonii model
    child[rep_row,(varNo+objDim)]=child_thickness[i,2]
  }
  
  
  #seed the next generation
  parentNext=rbind(parent[,-(which(colnames(parent)=="rnkIndex"):which(colnames(parent)=="crowdingDist"))],
                   child) #removing the rnkIndex and cd 
  parentNext=unique(parentNext[,1:(varNo+objDim)]) #to avoid repeated rows
  #ranking again
  ranking <- fastNonDominatedSorting(parentNext[,(varNo+1):(varNo+objDim)]);
  rnkIndex=integer(nrow(parentNext))
  i <- 1;
  while (i <= length(ranking)) {
    rnkIndex[ranking[[i]]] <- i;
    i <- i + 1;
  } 
  parentNext <- cbind(parentNext,rnkIndex, row.names = NULL);#so a warning message doesn't come about the row names
  #crowded comparison again
  objRange <- apply(parentNext[,(varNo+1):(varNo+objDim)], 2, max) - apply(parentNext[,(varNo+1):(varNo+objDim)], 2, min);
  cd <- crowdingDist4frnt(parentNext,ranking,objRange);
  parentNext <- cbind(parentNext,apply(cd,1,sum));
  colnames(parentNext)[ncol(parentNext)]="crowdingDist"
  parentNext.sort <- parentNext[order(parentNext[,which(colnames(parentNext)=="rnkIndex")],
                                      -parentNext[,which(colnames(parentNext)=="crowdingDist")]),];
  # choose the first 'popSize' rows for next generation
  parent <- parentNext.sort[1:popSize,]
  setwd(path1)
  write.csv(parent,file=paste("Parent_generation_",gen+1,".csv",sep=""))
}

#Plot the points in the parent file according to pareto dominance ranking 
#and the controls according to pareto dominance ranking 

setwd(path1)
file1=read.csv(paste("Parent_generation_",gen+1,".csv",sep=""),header=TRUE,row.names=1) #reads the parent of the last generation
pdf(file=paste("Pareto_Frontier after ",gen," generations.pdf",sep=""))
#In the multi-heterotroph model, x=Diversity_index and
#gingivalis gordonii it is Mygordonii
ggplot(file1,aes(x="Insert appropriate quantity",
                 y=Mean_biofilm_thickness,color=factor(rnkIndex)))+geom_point()+
  labs(color="Ranking of points",caption="Ranking of 1 indicates points on the Pareto frontier")


##In the multi-heterotroph model, Appropriate control1=Oxygen_conc
#Appropriate control2=Nitrate_conc
#gingivalis gordonii, Appropriate control1=Protein_conc
#Appropriate control2=Toxin_conc
p1=ggplot(file1,aes(x="Insert appropriate quantity",
                    y=Mean_biofilm_thickness,color=factor(rnkIndex),
                    shape=factor("Appropriate control1")))+
  geom_point()+
  scale_shape_manual(values=seq(0,15,1))
p2=ggplot(file1,aes(x="Insert appropriate quantity",
                    y=Mean_biofilm_thickness,color=factor(rnkIndex),
                    shape=factor("Appropriate control2")))+
  geom_point()+
  scale_shape_manual(values=seq(0,19,1))
library(gridExtra) #needed to get multiple ggplots in one frame
grid.arrange(p1,p2, ncol=1)
dev.off()
