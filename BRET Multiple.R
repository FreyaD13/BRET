BRETMultiple<-function(Experiments,
                       Figure,
                       export.means=TRUE, #Exporting the means is useful for combining independant experiments
                       #into one graph, treating each experiment as a data point. Set import.means to TRUE if
                       #this is what you want to do
                       export.plot=FALSE,
                       subset.ligands,
                       set.control.well,
                       set.control.row,
                       Normalise=FALSE){
  output<-list()

  if (!missing(Figure)){
    source("Figure.r",local=TRUE)
    assign("Experiments",get(Figure))
    if (exists((paste0(Figure,".subsetligands")))){
      assign("subset.ligands",get(paste0(Figure,".subsetligands")))
    }
    if (exists((paste0(Figure,".controlwell")))){
      assign("set.control.well",get(paste0(Figure,".controlwell")))
    }
    if (exists((paste0(Figure,".controlrow")))){
      assign("set.control.row",get(paste0(Figure,".controlrow")))
    }
  }
  
  Experiments.comp<-Experiments
  dflt.norm.factor<-"Missing"
  
    for (i in 1:length(Experiments)){
      if (export.means==TRUE){
        BRETexp<-BRET(Experiments[i],
                      save.plot = FALSE,
                      save.means=TRUE,
                      subset.ligands = subset.ligands,
                      set.control.well=set.control.well,
                      set.control.row=set.control.row,
                      nested.BRET=TRUE,
                      nested.run.no=i,
                      Normalise=Normalise,
                      dflt.norm.factor=dflt.norm.factor
                      )
        }
      if (export.means==FALSE){
        BRETexp<-BRET(Experiments[i],
                      save.plot = FALSE,
                      save.processed =TRUE,
                      subset.ligands = subset.ligands,
                      set.control.well=set.control.well,
                      set.control.row=set.control.row,
                      nested.BRET=TRUE,
                      nested.run.no=i,
                      Normalise=Normalise,
                      dflt.norm.factor=dflt.norm.factor
                      )
      }
      
      if (export.plot==TRUE){
        BRETexp<-BRET(Experiments[i],
                      save.plot = TRUE,
                      save.processed =FALSE,
                      subset.ligands = subset.ligands,
                      set.control.well=set.control.well,
                      set.control.row=set.control.row,
                      nested.BRET=TRUE,
                      nested.run.no=i,
                      Normalise=Normalise,
                      dflt.norm.factor=dflt.norm.factor
                      )
      }
      
      
      if (i==1 & !missing(Normalise)){
        dflt.norm.factor<-BRETexp$dflt.norm.factor
        BRETexp<-BRETexp[-1]
      }
      
      if (length(BRETexp$Warning)==0){
        output[[i]]<-BRETexp
        }
      
      if (length(BRETexp$Warning)==1){
        Experiments.comp<-Experiments.comp[Experiments.comp !=Experiments[i]]
        }
        
  }
  names(output)<-Experiments.comp
  output<-output[lapply(output,length)>0]
  return(output)
}
