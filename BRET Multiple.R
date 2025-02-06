BRETMultiple<-function(Experiments,
                       export.means=FALSE, #Exporting the means is useful for combining independant experiments
                       #into one graph, treating each experiment as a data point. Set import.means to TRUE if
                       #this is what you want to do
                       export.plot=FALSE){
  output<-list()

    for (i in 1:length(Experiments)){
      if (export.means==TRUE){
        BRETexp<-BRET(Experiments[i],
                      save.plot = FALSE,
                      save.means=TRUE)
        }
      if (export.means==FALSE){
        BRETexp<-BRET(Experiments[i],
                      save.plot = FALSE,
                      save.processed =TRUE)
      }
      
      if (export.plot==TRUE){
        BRETexp<-BRET(Experiments[i],
                      save.plot = TRUE,
                      save.processed =FALSE)
      }
      
    output[[i]]<-BRETexp
    
  }
  names(output)<-Experiments
  output<-output[lapply(output,length)>0]
  return(output)
}
