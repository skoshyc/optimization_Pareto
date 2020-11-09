#functions needed for the pareto frontier code

#theme for ggplot
theme_pm <- function () {
  theme_bw(base_size=12) +
    theme(
      panel.grid=element_line(linetype="dashed", color="light grey", size=0.2),
      axis.ticks.length=unit(-0.25, "cm"),
      axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5),
                                             "cm")),
      axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5),
                                             "cm"))
    )
}
#From http://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization#barplot-with-error-bars
#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summarized
# groupnames : vector of column names to be used as
# grouping variables 
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}



####to parse the xml file and change the relevant input of the solute
### xmlfile is the xml file we wish to read/modify
#The solutes we are considering for the Example model from idynomics is oxygen and nitrate
#Also the concentrations have to be within quotes
#resultfile is the xml file we wish to write into
xmlparser_multiheterotroph=function(xmlfile,O2concentration,Nconcentration,resultfile){
  #Before using the as_list function make sure to remove all comments in the xml document as
  #as_list reads the comments in as a list element and that causes problems when we use the function as_xml_document.
  xDoc=as_list(read_xml(xmlfile))#converting into list makes it easier to manipulate
  for(i in 1:length(xDoc$idynomics$world$bulk)){
    if(names(xDoc$idynomics$world$bulk)[[i]]=="solute"){
      if(attr(xDoc$idynomics$world$bulk[[i]],"name")=="o2d"){ #Changing the oxygen concentration
        #the first [[1]] denotes the first solute, second [[1]] denotes the first param which is Sbulk, third [[1]] denotes the value
        xDoc$idynomics$world$bulk[[i]][[1]][[1]]=O2concentration
        #the second param of the solute is Sin which we also have to change
        xDoc$idynomics$world$bulk[[i]][[2]][[1]]=O2concentration
      }
      if(attr(xDoc$idynomics$world$bulk[[i]],"name")=="MyNO3"){ #changing the nitrate concentration
        #the first [[1]] denotes the first solute, second [[1]] denotes the first param which is Sbulk, third [[1]] denotes the value
        xDoc$idynomics$world$bulk[[i]][[1]][[1]]=Nconcentration
        #the second param of the solute is Sin which we also have to change
        xDoc$idynomics$world$bulk[[i]][[2]][[1]]=Nconcentration
      }
    }
  }
  x1=as_xml_document(xDoc)
  write_xml(x1,file=resultfile)

}

#for the gingivalis-gordonii model
xmlparser_gingivalis=function(xmlfile,proteinconcentration,toxinconcentration,resultfile){
  #Before using the as_list function make sure to remove all comments in the xml document as
  #as_list reads the comments in as a list element and that causes problems when we use the function as_xml_document.
  xDoc=as_list(read_xml(xmlfile))#converting into list makes it easier to manipulate
  for(i in 1:length(xDoc$idynomics$world$bulk)){
    if(names(xDoc$idynomics$world$bulk)[[i]]=="solute"){
      if(attr(xDoc$idynomics$world$bulk[[i]],"name")=="MyProtein"){ #Changing the protein concentration
        #the first [[1]] denotes the first solute, second [[1]] denotes the first param which is Sbulk, third [[1]] denotes the value
        xDoc$idynomics$world$bulk[[i]][[1]][[1]]=proteinconcentration
        #the second param of the solute is Sin which we also have to change
        xDoc$idynomics$world$bulk[[i]][[2]][[1]]=proteinconcentration
      }
      if(attr(xDoc$idynomics$world$bulk[[i]],"name")=="MyToxin"){ #changing the toxin concentration
        #the first [[1]] denotes the first solute, second [[1]] denotes the first param which is Sbulk, third [[1]] denotes the value
        xDoc$idynomics$world$bulk[[i]][[1]][[1]]=toxinconcentration
        #the second param of the solute is Sin which we also have to change
        xDoc$idynomics$world$bulk[[i]][[2]][[1]]=toxinconcentration
      }
    }
  }
  x1=as_xml_document(xDoc)
  write_xml(x1,file=resultfile)
  
}




####Run the RunIdyno.py from the command prompt
#inputfile is the path of the xmlfile to be run. Default number of multiples=5
runpytho= function(xmlfolder,inputfile,multiples){
  #the scripts_to_start_idynomics folder has to be the working directory
  python_folder="../iDynoMiCS/scripts_to_start_idynomics" 
  setwd(python_folder)
  l1=paste("python RunIdyno.py ",xmlfolder,inputfile," --multiples=",multiples,sep="")
  system(l1)
}





####Run the RunIdyno.py from the command prompt
#inputfile is the path of the xmlfile to be run. Default number of multiples=5
#Changing the order of the runpytho inputs so as to accommodate parallel runs.
runpytho1= function(inputfile,xmlfolder,multiples){
  #the scripts_to_start_idynomics folder has to be the working directory
  python_folder="../iDynoMiCS/scripts_to_start_idynomics" 
  setwd(python_folder)
  l1=paste("python RunIdyno.py ",xmlfolder,inputfile," --multiples=",multiples,sep="")
  system(l1)
}




simpsonIndex_noEPS = function(resultFileFolder,numTimepoints,outputPeriod, folderForGraphOut){
  # Set up the graph file location 
  GRAPHFILE = paste(folderForGraphOut,"/Species_Diversity_Simpson.pdf",sep="")
  pdf(GRAPHFILE,width=10,height=7)
  
  speciesAbundance<-getSpeciesSpecificAbundance(resultFileFolder,numTimepoints,outputPeriod)
  #to remove the EPS species from the dataframe
  if(any(grepl('EPS',colnames(speciesAbundance)))=='TRUE'){
    speciesAbundance=speciesAbundance[,-which(grepl('EPS',colnames(speciesAbundance)))]
  }
  div_time <- NULL
  specAbunT<-t(speciesAbundance)
  for (i in 1:dim(specAbunT)[2]) 
  {
    # get the diversity of species
    div_time[i] <- vegan::diversity(specAbunT[,i], index = "simpson", MARGIN = 2)
    maxc <- max(div_time, na.rm = TRUE)
  }
  
  plot(div_time,  xlab = "Time course", ylab = "Species diversity (Simpson)", ylim= c(0,maxc), main = "Species diversity", las = 1, type = "b")
  
  dev.off()
  
  return(div_time)
}


speciesnumber_noEPS = function(resultFileFolder,numTimepoints,outputPeriod){
 
  
  speciesAbundance<-getSpeciesSpecificAbundance(resultFileFolder,numTimepoints,outputPeriod)
  
  #to remove the EPS species from the dataframe
  if(any(grepl('EPS',colnames(speciesAbundance)))=='TRUE'){
    speciesAbud=as.data.frame(speciesAbundance[,-which(grepl('EPS',colnames(speciesAbundance)))])
    colnames(speciesAbud)=colnames(speciesAbundance)[-which(grepl('EPS',colnames(speciesAbundance)))]
  }
  else{
    speciesAbud=speciesAbundance
    colnames(speciesAbud)=colnames(speciesAbundance)
  }
  
  return(speciesAbud)
}



