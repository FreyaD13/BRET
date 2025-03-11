BRET<-function(Experiment,
               Figure,
               import.data,
               Directory="Raw data/",
               lum.threshold=1000,
               Samples, #Cells in each row. Can be a list of 8 samples in the 
               #correct (A-H) order, or a single sample if they're 
               #consistent across the plate. See S2 for what to do 
               #with 2-4 repeating samples.
               Ligands,   #Treatment of the cells in each row. Use like Samples.
               Acceptor="Acceptor",  #Wavelength or identifier of acceptor channel.
               #Defaults to 'acceptor'
               Donor="Donor",     #Wavelength or identifier of donor channel
               #Defaults to 'donor'
               GProt=FALSE,   #Is this a g-protein BRET? does the data need inverted
               Concentrations=c("30", "10", "3", "1","0.3","0.1" ,"0.03" ,"0.01" ,"0.003","0.001","0.0003","0"), 
               #MUST be in uM and end in 0. Default is 
               #30,10,3,1 etc. but can be replaced by a list 
               #of 12 concentrations.
               seperate=FALSE,#Default false, analyse each plate seperately
               save.plot=TRUE,#Default true, save any plots made
               save.raw=FALSE, #Default false, save raw data for each plate
               save.means=FALSE,
               ignore.lig=FALSE, #do you want to ignore sorting by ligands?
               #for example if ligands only vary due to
               #different receptors
               save.processed=FALSE, #default is false, save tidy raw data
               colour, #if both sample+ligand are variables, which should
               #be indicated by colour (the other will be shape)
               find.ec50=FALSE,
               save.ec50.lines=FALSE,
               save.curve.fit=FALSE,
               subset.ligands,
               subset.samples,
               highlight.ec50,
               #CHANGED DEFAULT DESCRIPTION IS OLD;import.means and compare.exp are false by default meaning that
               #each data point will be treated as one data point in the final result
               #if compare.exp is set to true, all data points will still be
               #shown but results will be seperated by experiment
               #if import.means is set to true, the means from each experiment become
               #the new data points. this is necessary for avoiding pseudoreplication;
               #each experiment= 1 data point/1 replicate
               import.means=TRUE,
               compare.exp=FALSE,
               data.points=FALSE,
               set.line.resolution=0.001,
               constrain.lims=TRUE,
               ec_f, #to calculate eg. ec75
               error.bars=TRUE,
               subset.output=TRUE,
               set.control.well,
               set.control.row,
               constrain.hill=FALSE,
               Normalise=FALSE,
               nested.BRET=FALSE, #this just tells BRET whether it is nested within a larger BRET so can be used
               #to do things like turn off normalisation,
               nested.run.no,
               dflt.norm.factor
){#Set universal defaults
  
  if (!missing(Figure)){
    import.data<-BRETMultiple(Figure=Figure,
                              Normalise=Normalise)
    import.means=TRUE
  }
  
  if ((find.ec50==FALSE)&(save.plot==TRUE)){
    find.ec50=TRUE
  }
  
  output<-list()
  
  #first work out the data source; check there isn't more than one
  if ((!missing(Experiment))&(!missing(import.data))){
    warning("Multiple data sources provided, only raw data will be used")
    rm(import.data)
  }
  
  #check there is at least one
  if ((missing(Experiment))&(missing(import.data))){
    warning("NO DATA PROVIDED")
  }
  
  #Normal BRET Processing from source code
  if (!missing(Experiment)){
    Experiment<-as.character(Experiment)
    
    #create a list of all files containing the experiment stem in their name
    list<-grep(paste0(Experiment),list.files(Directory),value = TRUE)
    
    #import the layout data locally
    source("BRETLayouts.R",local=TRUE)
    
    #Check data is available for this experiment
    #If a file with the experiment name does not exist
    if (!exists(paste0(Experiment))&
        #plus no data has been otherwise supplied
        ((missing(GProt)) | (missing(Samples)) | (missing(Ligands)))) {
      #first look for a reason if one exists and paste a warning
      if (exists(paste0(Experiment,"_Missing"))){
        warning("Plate Layout Data Missing for ",Experiment,": ",get(paste0(Experiment,"_Missing")))
        #or print a generic explanation
      } else {warning("Plate Layout Data Missing for ",Experiment,": Reason not found.")}
    }
    
    #print any other warnings about the data
    if (exists(paste0(Experiment,"_Warning"))){
      warning(get(paste0(Experiment,"_Warning")))
    }
    
    #otherwise can proceed
    if (exists(paste0(Experiment))|
        ((exists("GProt")) & (exists("Samples")) & (exists("Ligands")))) {
    
      #check if there is a concentration modification. if so update the concentrations
      if ("Concentrations" %in% names(plate.layouts[[Experiment]])){
        Concentrations<-plate.layouts[[Experiment]]$Concentrations
      }
      
      #Import GProt, Samples and Ligands from layout data unless they have been
      #specified
      if (missing(GProt)){
        GProt<-plate.layouts[[Experiment]]$GProt
      }
      
      if (missing(Samples)){
        Samples<-plate.layouts[[Experiment]]$Samples
      }
      
      if (missing(Ligands)){
        Ligands<-plate.layouts[[Experiment]]$Ligands
      }
      
      #apply any additional settings from the plate data
      if (!is.null(plate.layouts[[Experiment]]$set.control.well)){
        set.control.well<-plate.layouts[[Experiment]]$set.control.well
      }
      
      if (!is.null(plate.layouts[[Experiment]]$set.control.row)){
        set.control.row<-plate.layouts[[Experiment]]$set.control.row
      }
      
      if (!is.null(plate.layouts[[Experiment]]$lum.threshold)){
        lum.threshold<-plate.layouts[[Experiment]]$lum.threshold
      }
      
      #set a Samples and Ligands source that wont be changed by looping
      SamplesSource<-Samples
      SamplesDims<-dim(data.frame(Samples))[2]
      LigandsSource<-Ligands
      LigandsDims<-dim(data.frame(Ligands))[2]
      
      #for every file in that list...
      for (file in 1:length(list)) {
        
        #.... first stick Raw data/ in front of it
        Data_Path<-paste0(Directory,list[file])
        
        #Then check set the samples and ligands for each plate if they're different
        
        #first check whether more than one list of samples has been provided. If only
        #one is provided it should be applied to each plate
        if (SamplesDims>1){
          #if there are more samples provided than the list of experiments it will
          #be an errror
          if (SamplesDims>length(list)){
            stop("More samples provided than experiments found")
          }
          #The opposite will also produce an error
          if (SamplesDims<length(list)){
            stop("More experiments found than samples provided")
          }
          #If the number of samples is the same as the length of the list it'll be
          #fine
          if (SamplesDims==length(list)){
            Samples<-t(SamplesSource[file])
          }
        }# end of set samples
        
        #Exactly the same thing is done for ligands
        if (LigandsDims>1){
          #error if there are more ligands than experiments 
          if (LigandsDims>length(list)){
            stop("More ligands provided than experiments found")
          } #end of error
          #error if more experiments than ligands
          if (LigandsDims<length(list)){
            stop("More experiments found than ligands provided")
          } #end of error
          #If no of ligands = number of experiments, the ligand list is set as
          #the number of that experiment (file)
          if (LigandsDims==length(list)){
            Ligands<-t(LigandsSource[file])
          }
        } #end of set ligands
        
        
        #then run the whole BRET function
        { #Save entire .csv as Bys
          Bys<-read_csv(Data_Path)
          
          #Extracts the specific rows with useful data in them
          Bys<-bind_rows(Bys[36:43,],Bys[47:54,],Bys[59:66,],Bys[70:77,])
          
          
          names(Bys)[2]<-"Char"
          
          #The second column is char by default, this changes it to a double
          Bys$Char<-as.double(Bys$Char)
          

            
          #first checks the correct number of concentrations have been supplied and
          #if not ends the function and prints error message
          if (length(Concentrations)==2){
            conc_1<-as.numeric(Concentrations[1])
            conc_2<-as.numeric(Concentrations[2])
            Concentrations[3]<-conc_1*10^-1
            Concentrations[4]<-conc_2*10^-1
            Concentrations[5]<-conc_1*10^-2
            Concentrations[6]<-conc_2*10^-2
            Concentrations[7]<-conc_1*10^-3
            Concentrations[8]<-conc_2*10^-3
            Concentrations[9]<-conc_1*10^-4
            Concentrations[10]<-conc_2*10^-4
            Concentrations[11]<-conc_1*10^-5
            Concentrations[12]<-0
          }
          
          if (!length(Concentrations)==12){
            stop("Wrong number of concentrations supplied")
          } else {names(Bys)[1]<-"Row"
           names(Bys)[2:13]<-Concentrations}
          
          #Establish what the vehicle concentration has been written as
          
          
          #Specifying intervals, first creates a new column and fills it with "Two"
          Bys$Interval <- "Two"
          #Then changes the value of the second half to "three"
          Bys[17:32,14]<-"Three"
          #If donor and acceptor channel names are unavailable they're just referred to
          #as donor and acceptor
          if (missing(Donor)){Donor<-"donor"}
          if (missing(Acceptor)){Acceptor<-"acceptor"}
          
          #Specifying channels, creates a new column and fills it with whichever
          #acceptor was specified
          Bys$Channel<-Acceptor
          #Then replaces every other 8 rows with Donor
          Bys[c(9:16,25:32),15]<-Donor
          #Creates new columns to fill with the Sample and Ligand for each row
          Bys$Sample<-"unspecified"
          Bys$Ligand<-"unspecified"
          
          #########################
          #Assign the sample/ligand labels to the correct rows. This is likely too 
          #convoluted and could be simplified
          
          #If statement turns the Samples into a list if they're not already
          if (!is.list(Samples)){
            #if not a blank list of x's is created
            Samples<-rep(Samples,8)
          }
          #Also if Samples is a list but there is only one element the same
          #must be done
          if (length(Samples)<8){
            Samples<-rep(Samples,ceiling(8/length(Samples)))
            Samples<-Samples[1:8]
          }
          
          #Exactly the same thing is done for ligands
          if (length(Ligands)<8){
            Ligands<-rep(Ligands,ceiling(8/length(Ligands)))
            Ligands<-Ligands[1:8]
          }
          
          #go through the data and assign the correct sample and ligand name to each row
          #by looking at which letter is in the 'Row' column.
          
          for (row in 1:dim(Bys)[1]) {
            if (Bys$Row[row]=="A") {Bys$Sample[row]=Samples[1]
            Bys$Ligand[row]=Ligands[1]}
            if (Bys$Row[row]=="B") {Bys$Sample[row]=Samples[2]
            Bys$Ligand[row]=Ligands[2]}
            if (Bys$Row[row]=="C") {Bys$Sample[row]=Samples[3]
            Bys$Ligand[row]=Ligands[3]}
            if (Bys$Row[row]=="D") {Bys$Sample[row]=Samples[4]
            Bys$Ligand[row]=Ligands[4]}
            if (Bys$Row[row]=="E") {Bys$Sample[row]=Samples[5]
            Bys$Ligand[row]=Ligands[5]}
            if (Bys$Row[row]=="F") {Bys$Sample[row]=Samples[6]
            Bys$Ligand[row]=Ligands[6]}
            if (Bys$Row[row]=="G") {Bys$Sample[row]=Samples[7]
            Bys$Ligand[row]=Ligands[7]}
            if (Bys$Row[row]=="H") {Bys$Sample[row]=Samples[8]
            Bys$Ligand[row]=Ligands[8]}}
          
          
          #Now everything has been done that requires specifying the actual row,
          #so rows can be removed if they contain values under threshold
          rows.removed<-FALSE
          
          #the automatic threshold is 1000
          if (missing(lum.threshold)){
            lum.threshold<-1000
          }
          
          
          #for the length of the Bys data frame
          for (row in 1:dim(Bys)[1]){
            #if any rows are under the set threshold, or contain NAs...
            if ((any((Bys[row,2:13]<lum.threshold)))|(any(is.na(Bys[row,2:13])))){
              #..first if this is the first row identified to be removed
              if (rows.removed==FALSE){
                #the first four elements of a rows.to.remove list are created and set as NAs
                rows.to.remove<-rep(NA,4)
                #a counter is set to 0
                counter<-0
                #rows.removed is changed to TRUE
                rows.removed<-TRUE
              }#end of if this is the first iteration
              
              #then for every row in the Bys data frame
              for (subrow in 1:dim(Bys)[1]){
                #if the row letter matches the one that has been identified as needing
                #to be removed
                if (Bys[row,1]==Bys[subrow,1]){
                  #the counter is updated
                  counter<-counter+1
                  #the row is added to the list
                  rows.to.remove[counter]<-subrow
                }#end of correct letter match
              }#end of checking letter matches across Bys
            }#end of correctly found NA in rows
          }#end of loop
          
          if (rows.removed==TRUE){
            #remove the rows.to.remove list from Bys
            Bys<-Bys[-rows.to.remove,]
            warning(paste0("Some rows contain values under threshold, set at ",
                           lum.threshold,
                           ". Adjust lum.threshold to include."))
          }
          
          
          
          ########
          #Tidying the data
          
          #Take all the emission values per concentration and put them in new columns
          #called emission and concentration
          Bys<-
            Bys |> 
            pivot_longer(
              #select all columns apart from the following, ie select only ones containing
              #concentration data
              cols = -c("Row", "Interval", "Channel","Sample","Ligand"),
              names_to = "Concentration", 
              values_to = "Emission",
              values_drop_na = TRUE,
              #make the concentration column numeric rather than a character so it will
              #display properly on the x-axis
              names_transform=list(Concentration=as.numeric))
          
          
          #!!!It may be useful to add a step here to identify anomalous results (ie
          #raw emission under a certain value)
          
          #Get the experiment ID from the name of the file. It is useful to tag the
          #data with this in case it needs further processing with other experiments
          Experiment_ID<-sub(".csv","", (sub("Raw data/","", Data_Path)))
          #add a new column containing the experiment ID
          Bys$Exp_ID<-Experiment_ID
          
          if (file==1){
            BysComb<-Bys
          } #If its the first plate it is called "BysComb"
          
          if (file>1){
            BysComb<-bind_rows(BysComb,Bys)
          } #if it isnt its added to byscomb
        }#end of this BRET function
      } #end of for loop for all files
      Bys<-BysComb
      veh<-(as.double(min(Bys$Concentration)))
      #Save to the global environment
      if (save.raw==TRUE){
        Bys$Exp_MasterID<-paste0(Experiment)
        output[["Raw"]]<-Bys
  #      assign((paste0(Experiment,'_Raw')),Bys,envir = .GlobalEnv)
        }
      
      #First split the emission values by interval and take an average, ending up
      #with a table containing Row, Channel, Sample, Concentration, Average Well
      #Emission and Exp_ID. Exp_ID is maintained here although it isn't used again
      #because i feel like it might be useful (if I want to save the intermediate
      #processed plots for example). Log_Conc is introduced here. It is in M rather
      #than uM so uM values are divided by 10^6.
      
      Int_Bys<-Bys|>
        pivot_wider(
          names_from=Interval,
          values_from=Emission) |>
        group_by(Row, Channel, Sample, Concentration,Ligand,Exp_ID) |> 
        summarise(AvWell_Emission=((Two+Three)/2),
                  Log_Conc=log10((Concentration/10^6)))
      
      #Then sort the average values by channel. Take a raw ratio of Acceptor
      #wavelength divided by Donor wavelength. Again, Exp_ID is maintained here in
      #case this data needs to be extracted. Get() is a very useful function here
      #as it allows you to refer to the specific value of an object, in this case
      #the channel names of the donor and acceptor.
      Chan_Bys<-Int_Bys|> 
        pivot_wider(
          names_from=Channel,
          values_from=AvWell_Emission) |>
        group_by(Row, Sample, Concentration,Log_Conc,Ligand,Exp_ID) |> 
        summarise(Raw_Ratio=(get(Acceptor)/get(Donor)))
      
      #A new column called ratio is made and filled with zeros.
      Chan_Bys$Ratio<-0
      #sort by Experiment_ID. This solves previous issue of the row vehicle being
      #set the same over the two plates
      Chan_Bys<-(Chan_Bys[order(Chan_Bys$Exp_ID),])
      
      #Each well is normalised to the control The code runs over every
      #row in the Chan_Bys dataframe.

      #if the data needs to be normalised to a specific row or well first the 
      #data will be filtered to find those wells if a specific well, the data 
      #is filtered based first off the row (ie the Ligand or Sample specified) 
      #and then the control well(s) is/are taken and averaged
      
      if (!missing(set.control.well)){
        control.ratio<-
          mean((filter(Chan_Bys,Ligand==set.control.well|Sample==set.control.well, 
                       Concentration==0))$Raw_Ratio)
      }
      
      #if the data should be normalised to an entire row, that row (based off sample/ligand name) is found and averaged 
      if (!missing(set.control.row)){
        control.ratio<-
          mean((filter(Chan_Bys,Ligand==set.control.row|Sample==set.control.well))$Raw_Ratio)
      }
      
      for (row in 1:dim(Chan_Bys)[1]){
        if (exists("control.ratio")){
          if (GProt==FALSE){
            Chan_Bys$Ratio[row]=((Chan_Bys$Raw_Ratio[row])-(control.ratio))
          } else if (GProt==TRUE){
            Chan_Bys$Ratio[row]=((control.ratio)-(Chan_Bys$Raw_Ratio[row]))
          } #end of g prot specific well
        } #end of normalise to specific row/well
        
        if (!exists("control.ratio")){
          if (GProt==FALSE){
            if ((Chan_Bys$Concentration[row]==veh)){
              Row_Vehicle<-Chan_Bys$Raw_Ratio[row]
            }
            #...each row raw_ratio is subtracted from the "Row_Vehicle"
            Chan_Bys$Ratio[row]=((Chan_Bys$Raw_Ratio[row])-(Row_Vehicle))
          } else if (GProt==TRUE){
            if ((Chan_Bys$Concentration[row]==veh)){
              Row_Vehicle<-Chan_Bys$Raw_Ratio[row]
            }
            #...each row raw_ratio is subtracted from the "Row_Vehicle"
            Chan_Bys$Ratio[row]=((Row_Vehicle)-(Chan_Bys$Raw_Ratio[row]))
          }
        }
      }
      
      
      #NORMALISE
      
  #    if ((!missing(Normalise))&nested.BRET==FALSE){
      if (!missing(Normalise)&(!Normalise==FALSE)){  
        if ((missing(import.data)&missing(Figure))|nested.BRET==TRUE){
          Normalisation<-filter(Chan_Bys,Sample==Normalise|Ligand==Normalise)
          if (dim(unique(Normalisation[c('Sample','Ligand')]))[1]>1){
            output[["Warning"]]<-"Normalisation Failed"
            warning("Multiple sample types per ligand or ligands per sample for normalisation.No solution developed")
          } else if (dim(Normalisation)[1]==0){
            output[["Warning"]]<-"Normalisation Failed"
            warning(paste0(Normalise," (required for normalisation) was not run for experiment ", Experiment," & experiment will be excluded from analysis."))
          } else if (dim(Normalisation)[1]>0) {
            Normalisation_Curve<-drm(
              #looking at ratio and concentration
              formula=Ratio~Concentration,
              #looking in the subs dataframe (just created as subset of Chan_Bys)
              data=Normalisation,
              fct = LL.4(names=c('hill','min_value','max_value','ec_50'))
              )
            norm.max<-Normalisation_Curve$coefficients[3]
            norm.factor<-100/norm.max
            Chan_Bys$Ratio<-Chan_Bys$Ratio*norm.factor}
          }
        }
  
      #the average and standard error are taken to plot the error bars. New
      #Experiment ID is introduced. Important in case these averages are used to
      #combine again, needs to be consistent with mean output of BRET() function.
      Av_Bys<-Chan_Bys|>
        group_by(Sample,Concentration,Log_Conc,Ligand)|>
        summarise(Exp_ID=paste0(Experiment),
                  m_ratio=mean(Ratio),
                  sem_ratio=sd(Ratio)/sqrt(n()))
      
    }
  } #end of normal BRET processing from raw data
  
  #IMPORT DATA from previous BRET
  if (!missing(import.data)){
    Experiment<-"BRET"
    
    for (p in 1:length(import.data)){
      if (p==1){
        Chan_Bys<-import.data[[1]]
        Chan_Bys<-Chan_Bys[[1]]
      }
      if (p>1){
        temp<-import.data[[p]]
        Chan_Bys<-bind_rows(Chan_Bys,temp[[1]])
      }
    }
    
    if (import.means==TRUE){
      Chan_Bys<-rename(.data=Chan_Bys,
                       Ratio=m_ratio,
                       Error=sem_ratio)
      
      if (!missing(Normalise)){
        Normalisation<-filter(Chan_Bys,Sample==Normalise|Ligand==Normalise)
        if (dim(unique(Normalisation[c('Sample','Ligand')]))[1]>1){
          warning("Multiple sample types per ligand or ligands per sample for normalisation.No solution developed")
        }
        
        Normalisation_Curve<-drm(
          #looking at ratio and concentration
          formula=Ratio~Concentration,
          #looking in the subs dataframe (just created as subset of Chan_Bys)
          data=Normalisation,
          fct = LL.4(names=c('hill','min_value','max_value','ec_50'))
        )
        
        norm.max<-Normalisation_Curve$coefficients[3]
        norm.factor<-100/norm.max
        Chan_Bys$Ratio<-Chan_Bys$Ratio*norm.factor
        
      }
      
      Av_Bys<-Chan_Bys|>
        group_by(Sample,Concentration,Log_Conc,Ligand)|>
        summarise(m_ratio=mean(Ratio),
                  sem_ratio=sd(Ratio)/sqrt(n()))
    }
    
    if (import.means==FALSE){
      if (compare.exp==FALSE){ ##DEFAULT!
        
        if (!missing(Normalise)){
          Normalisation<-filter(Chan_Bys,Sample==Normalise|Ligand==Normalise)
          if (dim(unique(Normalisation[c('Sample','Ligand')]))[1]>1){
            warning("Multiple sample types per ligand or ligands per sample for normalisation.No solution developed")
          }
          
          Normalisation_Curve<-drm(
            #looking at ratio and concentration
            formula=Ratio~Concentration,
            #looking in the subs dataframe (just created as subset of Chan_Bys)
            data=Normalisation,
            fct = LL.4(names=c('hill','min_value','max_value','ec_50'))
          )
          
          norm.max<-Normalisation_Curve$coefficients[3]
          norm.factor<-100/norm.max
          Chan_Bys$Ratio<-Chan_Bys$Ratio*norm.factor
          
        }
        
        Av_Bys<-Chan_Bys|>
          group_by(Sample,Concentration,Log_Conc,Ligand)|>
          summarise(m_ratio=mean(Ratio),
                    sem_ratio=sd(Ratio)/sqrt(n()))
      }
      if (compare.exp==TRUE){
        
        if (!missing(Normalise)){
          warning("Cannot yet normalise data when comparing experiments")
        }
        
        Av_Bys<-Chan_Bys|>
          group_by(Sample,Concentration,Log_Conc,Ligand,Exp_MasterID)|>
          summarise(m_ratio=mean(Ratio),
                    sem_ratio=sd(Ratio)/sqrt(n()))
        Av_Bys<-Av_Bys|>
          rename(Exp_ID=Exp_MasterID)
      }
    }
  }
  

  if (subset.output==FALSE){
    if (save.means==TRUE){
      Av_Bys$Exp_MasterID<-paste0(Experiment)
      output[["Means"]]<-Av_Bys
      
    }
    if (save.processed==TRUE){
      Chan_Bys$Exp_MasterID<-paste0(Experiment)
      output[["Processed"]]<-Chan_Bys
    }
  }
  
  if (!missing(subset.ligands)){
    Av_Bys<-subset(Av_Bys,Ligand %in% subset.ligands)
    Chan_Bys<-subset(Chan_Bys, Ligand %in% subset.ligands)
  }
  
  if (!missing(subset.samples)){
    Av_Bys<-subset(Av_Bys,Sample %in% subset.samples)
    Chan_Bys<-subset(Chan_Bys, Sample %in% subset.samples)
  }
  
  if (subset.output==TRUE){
    if (save.means==TRUE){
      Av_Bys$Exp_MasterID<-paste0(Experiment)
      output[["Means"]]<-Av_Bys
      
    }
    if (save.processed==TRUE){
      Chan_Bys$Exp_MasterID<-paste0(Experiment)
      output[["Processed"]]<-Chan_Bys
    }
  }
  
  ########THIS POINT TO NORMALISE
  
  
  if (find.ec50==TRUE){
    #EC50 CALCULATION
    #part1:estimation of EC50, this creates df of estimated ec50s 
    #set defaults
    if (missing(save.ec50.lines)){
      save.ec50.lines=FALSE
    }
    if (missing(save.curve.fit)){
      save.curve.fit=FALSE
    }
    
    #first create a dataframe of all unique sample and ligand combinations
    if (compare.exp==FALSE){
      VariableUnq<-as.data.frame(
        unique(Chan_Bys[c('Sample','Ligand')]))
    } else if (compare.exp==TRUE) {
      VariableUnq<-as.data.frame(
        unique(Chan_Bys[c('Sample','Ligand','Exp_ID')]))
    }
    
    
    #this for loop does the entire thing, it loops over the list of unique variables
    #and analyses the data for them one by one
    for (var in 1:dim(VariableUnq)[1]){
     
      if (compare.exp==FALSE){
        subs<-filter(Chan_Bys,Sample==(VariableUnq[var,1]),Ligand==(VariableUnq[var,2]))}
      
      if (compare.exp==TRUE){
        subs<-filter(Chan_Bys,Sample==(VariableUnq[var,1]),Ligand==(VariableUnq[var,2]),Exp_ID==(VariableUnq[var,3]))
      }
      
      #create the ec50 estimate
      {try(curve_fit<-drm(
        #looking at ratio and concentration
        formula=Ratio~Concentration,
        #looking in the subs dataframe (just created as subset of Chan_Bys)
        data=subs,
        fct = LL.4(names=c('hill','min_value','max_value','ec_50'))
      ))
        }
      
      
      if (exists("curve_fit")){
        
        if (compare.exp==FALSE){
          VariableUnq[var,3:6]<-curve_fit$coefficients
        }
        
        if (compare.exp==TRUE){
          VariableUnq[var,4:7]<-curve_fit$coefficients
        }
        
        
        ############################################     
        #IF EC50 LABELLING LINES ARE NEEDED
        #this will create a tribble of data for the addition of lines indicating the
        #ec50 on the plot. Default is FALSE but will also run if highlight.ec50
        #values are provided
        
        if ((save.ec50.lines==TRUE)|(!missing(highlight.ec50))){
          #find the midpoint (ie y axis value) of the predicted values
          
          if (constrain.lims==TRUE){
            curve_fit$coefficients["min_value:(Intercept)"]<-0
          }
          
          
          mid<-(curve_fit$coefficients["max_value:(Intercept)"]
                +curve_fit$coefficients["min_value:(Intercept)"])/2
          
          #find the pEC50: the negative log of the ec50
          ec50<-(curve_fit$coefficients["ec_50:(Intercept)"])
          pec50<-log10((ec50/10^6))
          
          #put these values into a tribble
          ec50_plot_lines <- tribble(
            ~x,   ~xend, ~y,   ~yend,
            -Inf, pec50,  mid, mid,
            pec50,pec50,mid,-Inf
          )
          
          #label these rows with what the current variables are
          
          ec50_plot_lines$Sample<-var2Samp
          ec50_plot_lines$Ligand<-var2Lig
          if (compare.exp==TRUE){
            ec50_plot_lines$Exp_MasterID<-var2Exp
            
          }
          
          ec50_plot_lines$pEC50<-pec50
          ec50_plot_lines$EC50<-ec50
          
          #if this is the first variable a master list is created
          if (var==1){ec50_lines_master<-ec50_plot_lines}
          #if this is the second variable the values are added to the master list
          if (var>1){ec50_lines_master<-bind_rows(ec50_lines_master,ec50_plot_lines)}
        }
        rm(curve_fit)
      }#end of loop for curve_fit works & exists
    } #end of for loop over variables
    
    if (compare.exp==FALSE){
      names(VariableUnq)[3:6]<-c("hill","min_value","max_value","ec_50")
    }
      
    if (compare.exp==TRUE){
      names(VariableUnq)[4:7]<-c("hill","min_value","max_value","ec_50")
    }
    
    #after the loop is complete first subset if relevant and then
    #save relevant data to the global environment
    
    if (!missing(subset.ligands)){
      VariableUnq<-subset(VariableUnq, Ligand %in% subset.ligands)
      if ((save.ec50.lines==TRUE)){
        ec50_lines_master<-subset(ec50_lines_master, Ligand %in% subset.ligands)
      }
    }
    
    if (!missing(subset.samples)){
      VariableUnq<-subset(VariableUnq, Sample %in% subset.samples)
      if ((save.ec50.lines==TRUE)){
        ec50_lines_master<-subset(ec50_lines_master, Sample %in% subset.ligands)
      }
    }
    
    VariableUnq$pec_50<-log10(VariableUnq$ec_50/(10^6))
    
    #Add Ec75 or other details
    if (!missing(ec_f)){
      VariableUnq[paste0("EC",ec_f)]<-VariableUnq$ec_50*((ec_f/(100-ec_f))^(1/-VariableUnq$hill))
      VariableUnq[paste0("pEC",ec_f)]<-log10(VariableUnq$ec_50*((ec_f/(100-ec_f))^(1/-VariableUnq$hill))/(10^6))
    }
    
    output[["CurveParams"]]<-VariableUnq
    
    if ((save.ec50.lines==TRUE)){
      output[["EC50Lines"]]<-ec50_lines_master
    }
    
    #CREATING GRAPHS
    if (save.plot==TRUE){
      
      
      #set parameters if there is no curve available (ie. if the line 
      #should be flat). important to do this before setting constraints
      for (f in 1:dim(VariableUnq)[1]){
        if (sum(is.na(VariableUnq[f,]))>0|(round(VariableUnq$hill[f],3)==0)){
          VariableUnq$hill[f]<-1
          VariableUnq$min_value[f]<-0
          VariableUnq$max_value[f]<-0
          VariableUnq$ec_50[f]<-1
        }
        #    VariableUnq$hill[f]<-(-1)
      }
      
      if (constrain.lims==TRUE){
        mean.max<-mean(filter(VariableUnq,hill>0)$max_value)
        for (f in 1:dim(VariableUnq)[1]){
          #if hill is over 0 for some reason this is a negative curve so the max must
          #be constrained
          if (VariableUnq$hill[f]>0){
            VariableUnq$max_value[f]<-mean.max
            #but if its under then its a positive curve so the max should be constrained
          } else if (VariableUnq$hill[f]<0){
            VariableUnq$min_value[f]<-0
          }
          
          if (!missing(Normalise)){
            if (VariableUnq$max_value[f]>100){
              VariableUnq$max_value[f]<-100
            }
        }
        }
      }
      
      if (constrain.hill==TRUE){
        for (f in 1:dim(VariableUnq)[1]){
          #if hill is over 0 for some reason this is a negative curve so the max must
          #be constrained
          if (VariableUnq$hill[f]>0){
            VariableUnq$hill[f]<-1
            #but if its under then its a positive curve so the max should be constrained
          } else if (VariableUnq$hill[f]<0){
            VariableUnq$hill[f]<-(-1)
          }
        }
      }
      
      
      # 
      # if (exists("control.ratio")){
      #   constrain.lims=FALSE
      #   VariableUnq$max_value<-0
      # }
      # 
      # if (constrain.lims==TRUE){
      #   VariableUnq$min_value<-0
      # }
      # 
      
     
      
      #standard graphs, no experiment comparisons
      if (compare.exp==FALSE){
        #there is only one ligand or you want to ignore extra ligands
        if ((length(unique(VariableUnq$Ligand))==1)|(ignore.lig==TRUE)) {
          #the agonist is set to the name of the ligand
          if (ignore.lig==TRUE){
            Agonist<-"Agonist"}
          else {Agonist <- VariableUnq$Ligand[1]}
          
          #If there is also only one sample
          if ((length(unique(VariableUnq$Sample))==1)){
            bys_plot <- ggplot(data=Av_Bys,
                               mapping=aes(x=Log_Conc,
                                           y=m_ratio)
            )+
              lapply(1:dim(VariableUnq)[1], function(i){
                eq=function(x){
                  VariableUnq$min_value[i]+(
                    (VariableUnq$max_value[i]-VariableUnq$min_value[i])/
                      (1+10^((x-(log10(((VariableUnq$ec_50[i])/10^6))))*VariableUnq$hill[i]))
                  )
                }
                line<-geom_function(fun=eq,
                                    size=1.5)})
            
            if (data.points==TRUE){
              bys_plot<-bys_plot+
              #Raw data points
              geom_point(data=Chan_Bys,
                         mapping= aes(x=Log_Conc,
                                      y=Ratio),
                         shape=21,fill="white")
            }
            
            
            
          } #otherwise if there are multiple samples but only one ligand
          else {
            bys_plot <- ggplot(data=Av_Bys,
                               mapping=aes(x=Log_Conc,
                                           y=m_ratio,
                                           colour=Sample)
            ) +
              lapply(1:dim(VariableUnq)[1], function(i){
                eq=function(x){
                  VariableUnq$min_value[i]+(
                    (VariableUnq$max_value[i]-VariableUnq$min_value[i])/
                      (1+10^((x-(log10(((VariableUnq$ec_50[i])/10^6))))*VariableUnq$hill[i]))
                  )
                }
                line<-geom_function(fun=eq,
                                    aes(colour=VariableUnq$Sample[i]),
                                    size=1.5)})
            
            
            if (data.points==TRUE){
              bys_plot<-bys_plot+
              #Raw data points
              geom_point(data=Chan_Bys,
                         mapping= aes(x=Log_Conc,
                                      y=Ratio,
                                      colour=Sample),
                         shape=21,fill="white")
            }
            
          }  #end of multiple samples, only one ligand
          
          #Set standard colours BY SAMPLE
          Samples<-VariableUnq$Sample
          source("StandardColours.r",local=TRUE)
          for (i in 1:length(Samples)){
            if (i==1){
              extra.colours<-c("#A6CEE3","#FB9A99","#66A61E","#A6761D","#666666","#E41A1C","#377EB8","#7570B3","#A65628","#F781BF","#E6AB02")
              coloursumm<-c()
            }
            sampcol<-filter(Sample.Colours,Sample==Samples[i])
            if (dim(sampcol)[1]==0){
              coloursumm[i]<-extra.colours[1]
              extra.colours<-extra.colours[-1]
            }
            if (dim(sampcol)[1]>0){
              coloursumm[i]<-sampcol[1,2]
            }
          }
          
          bys_plot<-bys_plot+
            scale_colour_manual(breaks=Ligands,
                                values=coloursumm)
          
        }#end of only one ligand
        else #there is more than one ligand but.....
        {
          Agonist <- "Agonist"
          ##..... there is only one sample
          if (length(unique(VariableUnq$Sample))==1){
            #The colour is by ligand
            bys_plot <- ggplot(data=Av_Bys,
                               mapping=aes(x=Log_Conc,
                                           y=m_ratio,
                                           colour=Ligand)
            )+
              lapply(1:dim(VariableUnq)[1], function(i){
                eq=function(x){
                  VariableUnq$min_value[i]+(
                    (VariableUnq$max_value[i]-VariableUnq$min_value[i])/
                      (1+10^((x-(log10(((VariableUnq$ec_50[i])/10^6))))*VariableUnq$hill[i]))
                  )
                }
                line<-geom_function(fun=eq,
                                    aes(colour=VariableUnq$Ligand[i]),
                                    size=1.5)
              }
              )
            
            if (data.points==TRUE){
              bys_plot<-bys_plot+
              geom_point(data=Chan_Bys,
                         mapping= aes(x=Log_Conc,
                                      y=Ratio,
                                      colour=Ligand),
                         shape=21,fill="white"
              )
            }
            
            #Set standard colours
            Ligands<-VariableUnq$Ligand
            source("StandardColours.r",local=TRUE)
            for (i in 1:length(Ligands)){
              if (i==1){
                extra.colours<-c("#A6CEE3","#FB9A99","#66A61E","#A6761D","#666666","#E41A1C","#377EB8","#7570B3","#A65628","#F781BF","#E6AB02")
                coloursumm<-c()
              }
              ligcol<-filter(Ligand.Colours,Ligand==Ligands[i])
              if (dim(ligcol)[1]==0){
                coloursumm[i]<-extra.colours[1]
                extra.colours<-extra.colours[-1]
              }
              if (dim(ligcol)[1]>0){
                coloursumm[i]<-ligcol[1,2]
              }
            }
            
            bys_plot<-bys_plot+
              scale_colour_manual(breaks=Ligands,
                                  values=coloursumm)
            
            }#end of there is only one sample
          else #else there are multiple samples and multiple ligands
          {
            #The default is that colour will indicate ligand and shape sample
            if (missing(colour)){colour<-"Ligand"}
            #if you want colour to indicate the ligand
            if (colour=="Ligand"){
              bys_plot <- ggplot(data=Av_Bys,
                                 mapping=aes(x=Log_Conc,
                                             y=m_ratio,
                                             colour=Ligand,
                                             shape=Sample
                                 )
              )+
                #if colour indicates ligand, sample will be indicated by shape
                #and linetype
                lapply(1:dim(VariableUnq)[1], function(i){
                  eq=function(x){
                    VariableUnq$min_value[i]+(
                      (VariableUnq$max_value[i]-VariableUnq$min_value[i])/
                        (1+10^((x-(log10(((VariableUnq$ec_50[i])/10^6))))*VariableUnq$hill[i]))
                    )
                  }
                  line<-geom_function(fun=eq,
                                      aes(colour=VariableUnq$Ligand[i],
                                          linetype=VariableUnq$Sample[i]),
                                      size=1.5)})
              
              
              if (data.points==TRUE){
                bys_plot<-bys_plot+
                geom_point(data=Chan_Bys,
                           mapping= aes(x=Log_Conc,
                                        y=Ratio,
                                        colour=Ligand,
                                        shape=Sample)
                )
              } 
            }#end of colour = ligand, shape = sample
            else
            {
              bys_plot <- ggplot(data=Av_Bys,
                                 mapping=aes(x=Log_Conc,
                                             y=m_ratio,
                                             colour=Sample,
                                             shape=Ligand
                                 )
              )+
                #if colour indicates sample, ligand will be indicated by shape
                #and linetype
                lapply(1:dim(VariableUnq)[1], function(i){
                  eq=function(x){
                    VariableUnq$min_value[i]+(
                      (VariableUnq$max_value[i]-VariableUnq$min_value[i])/
                        (1+10^((x-(log10(((VariableUnq$ec_50[i])/10^6))))*VariableUnq$hill[i]))
                    )
                    }
                  line<-geom_function(fun=eq,
                                      aes(colour=VariableUnq$Sample[i],
                                          linetype=VariableUnq$Linetype[i]),
                                      size=1.5)})
              
              
              if (data.points==TRUE){
                bys_plot<-bys_plot+
                geom_point(data=Chan_Bys,
                           mapping= aes(x=Log_Conc,
                                        y=Ratio,
                                        colour=Sample,
                                        shape=Ligand)
                ) 
              }
              
            }#end of colour = sample, shape =ligand
          }#end of there are multiple samples and multiple ligands
        }#end of there is more than one ligand
      }
      
      #if you're comparing between experiments this code will always set
      #exp to colour and the other parameter to shape. It will not allow comparisons
      #between all three parameters...
      if (compare.exp==TRUE){
        #First check you're not comparing three parameters. If ignore.ligand is 
        #true this is not an issue
        if (ignore.lig==FALSE){
          if ((length(unique(VariableUnq$Ligand))>1)&(length(unique(VariableUnq$Sample))>1)){
            warning("This function can not create graphs with more than two parameters.
                Subset samples or ligands to compare experiments.")
          }
        }#end of warning if you're trying to make a graph comparing experiments,
        #ligands and samples
        
        #there is only one ligand or you want to ignore extra ligands
        if ((length(unique(VariableUnq$Ligand))==1)|(ignore.lig==TRUE)) {
          #if there are multiple ligands but you're ignorning that
          #the agonist is set to agonist
          if (ignore.lig==TRUE){
            Agonist<-"Agonist"}
          
          #otherwise the agnoist is set to the actual name of the agonist
          else {Agonist <- VariableUnq$Ligand[1]}
          
          #If there is also only one sample
          if ((length(unique(VariableUnq$Sample))==1)){
            bys_plot <- ggplot(data=Av_Bys,
                               mapping=aes(x=Log_Conc,
                                           y=m_ratio,
                                           #the colour is set to ID
                                           colour=Exp_ID)
            ) +
              lapply(1:dim(VariableUnq)[1], function(i){
                eq=function(x){
                  VariableUnq$min_value[i]+(
                    (VariableUnq$max_value[i]-VariableUnq$min_value[i])/
                      (1+10^((x-(log10(((VariableUnq$ec_50[i])/10^6))))*VariableUnq$hill[i]))
                  )
                }
                line<-geom_function(fun=eq,
                                    #colour is by experiment ID
                                    aes(colour=VariableUnq$Exp_ID[i]),
                                    size=1.5)})
            
            
            
            if (data.points==TRUE){
              bys_plot<-bys_plot+#Raw data points
                geom_point(data=Chan_Bys,
                           mapping= aes(x=Log_Conc,
                                        y=Ratio,
                                        #raw data points added by Exp_MasterID
                                        colour=Exp_MasterID),
                           shape=21,fill="white")
            }
              
            # end of [compare.exp] one ligand one sample
          } else if ((length(unique(Av_Bys$Sample))>1)) {   #otherwise if there are MULTIPLE SAMPLES 
            bys_plot <- ggplot(data=Av_Bys,
                               mapping=aes(x=Log_Conc,
                                           y=m_ratio,
                                           colour=Exp_ID,
                                           shape=Sample)
            ) +
              lapply(1:dim(VariableUnq)[1], function(i){
                eq=function(x){
                  VariableUnq$min_value[i]+(
                    (VariableUnq$max_value[i]-VariableUnq$min_value[i])/
                      (1+10^((x-(log10(((VariableUnq$ec_50[i])/10^6))))*VariableUnq$hill[i]))
                  )
                }
                line<-geom_function(fun=eq,
                                    #colour is by experiment ID, linetype is by sample
                                    aes(colour=VariableUnq$Exp_ID[i],
                                        linetype=VariableUnq$Sample[i]),
                                    size=1.5)})
            
            if (data.points==TRUE){
              bys_plot<-bys_plot+
              #Raw data points
              geom_point(data=Chan_Bys,
                         mapping= aes(x=Log_Conc,
                                      y=Ratio,
                                      colour=Exp_MasterID,
                                      shape=Sample))
            }
            
            
          } #end of [compare.exp] one ligand multiple samples
          #end of [compare.exp] only one ligand
        } else if ((length(unique(Av_Bys$Ligand))>1)) #there are MUTLIPLE LIGANDS (therefore only one sample)
        {Agonist <- "Agonist"
        #The colour is by ligand
        bys_plot <- ggplot(data=Av_Bys,
                           mapping=aes(x=Log_Conc,
                                       y=m_ratio,
                                       colour=Exp_ID,
                                       shape=Ligand)
        )+
          lapply(1:dim(VariableUnq)[1], function(i){
            eq=function(x){
              VariableUnq$min_value[i]+(
                (VariableUnq$max_value[i]-VariableUnq$min_value[i])/
                  (1+10^((x-(log10(((VariableUnq$ec_50[i])/10^6))))*VariableUnq$hill[i]))
              )
            }
            line<-geom_function(fun=eq,
                                #colour is by experiment ID, linetype is by ligand
                                aes(colour=VariableUnq$Exp_ID[i],
                                    linetype=VariableUnq$Ligand[i]),
                                size=1.5)})
        
        if (data.points==TRUE){
          bys_plot<-bys_plot+
          geom_point(data=Chan_Bys,
                     mapping= aes(x=Log_Conc,
                                  y=Ratio,
                                  shape=Exp_MasterID,
                                  colour=Ligand)
          )
        }
        
        }#end of [compare.exp] there are multiple ligands therefore one sample
      }# end of compare.exp
      
      
      #format the graph
      
      #find min& max concentrations for scaling
      min.conc<-min(Av_Bys$Concentration[Av_Bys$Concentration>0])
      max.conc<-max(Av_Bys$Concentration)
      
      bys_plot<-bys_plot+
        #Labels, the x-axis can be labelled with the specific agonist 
        xlab(paste0('Log[',Agonist,'] M')) +
        ylab('BRET Ratio (Relative to Control)') +
        #title is the experiment ID
        #ggtitle((paste0('Experiment: ',Experiment)))+
        #scale x axis
        scale_x_continuous(n.breaks=0.5*length(unique(Av_Bys$Concentration)))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text.x=element_text(size=12, face="bold"),
              axis.text.y=element_text(size=12, face="bold"))+
        xlim((round((log10(((as.numeric(min.conc))/10^6))),1)-0.25),
             (round((log10(((as.numeric(max.conc))/10^6))),1)+0.25))+ 
        theme(legend.position = c(0.2,0.75))+
        ylim(-0.03,0.20)
      
      if (data.points==FALSE){
        bys_plot<-bys_plot+
          geom_point(size=3,shape=21,fill="white")
      }
      
      if (exists("control.ratio")){
        bys_plot<-bys_plot+
          theme(legend.position = c(0.2,0.5))
      }
      
      ##set y lims based on whether g prot or bystander
      if (!missing(GProt)){
        if (GProt==TRUE){
          bys_plot<-bys_plot+
            ylim(-0.03,0.20)
          } else if (GProt==FALSE){
            bys_plot<-bys_plot+
              ylim(-0.015,0.10)
          }
      }
      
      
      
      if (missing(GProt)){
        if (max(Av_Bys$m_ratio>0.1)){
          bys_plot<-bys_plot+
            ylim(-0.03,0.20)
        } else if (max(Av_Bys$m_ratio<0.1)) {
          bys_plot<-bys_plot+
            ylim(-0.015,0.10)
        }
      }
      
      if (!missing(Normalise)){
        bys_plot<-bys_plot+
          ylim(-25,150)
      }
      
      if (error.bars==TRUE){
        bys_plot<-bys_plot+
          #error bars are standard error
          geom_errorbar(aes(ymin = m_ratio - sem_ratio,
                            ymax = m_ratio + sem_ratio),
                        width = 0.3)
      }
      
      if (!missing(highlight.ec50)){
        #create a subset list of the EC50 lines. first split the sample and ligand
        for (lin in 1:length(highlight.ec50)){
          #split the string and extract it as a normal list
          to.highlight<-strsplit(highlight.ec50[lin],"[+]")
          to.highlight<-to.highlight[[1]]
          
          #reset each highlight to the full list of variables
          samp.highlight<-as.list(ec50_lines_master$Sample)
          lig.highlight<-as.list(ec50_lines_master$Ligand)
          if (compare.exp==TRUE){
          exp.highlight<-as.list(ec50_lines_master$Exp_MasterID)}
          
          #work out what each variable is and overwrite the full list with
          #the specific if it applies
          for (lin2 in 1:length(to.highlight)){
            if (to.highlight[lin2] %in% ec50_lines_master$Sample){
              samp.highlight<-to.highlight[lin2]
            }
            if (to.highlight[lin2] %in% ec50_lines_master$Ligand){
              lig.highlight<-to.highlight[lin2]
            }
            if (compare.exp==TRUE){
            if (to.highlight[lin2] %in% ec50_lines_master$Exp_MasterID){
              exp.highlight<-to.highlight[lin2]
            }
            }
          }
          
          #apply the filter
          if (compare.exp==FALSE){
            active.lines.temp<-ec50_lines_master|>
              filter(
                Sample==samp.highlight,
                Ligand==lig.highlight
              )
          } else if (comoare.exp==TRUE){
            active.lines.temp<-ec50_lines_master|>
              filter(
                Sample==samp.highlight,
                Ligand==lig.highlight,
                Exp_MasterID==exp.highlight
              ) 
            }
          
          #if this is the first loop it is called the master list
          if (lin==1){
            active.lines<-active.lines.temp}
          #if this is the second loop this list is added to the master list
          if(lin==2){
            active.lines<-bind_rows(active.lines,active.lines.temp)}
          if (lin==3){
          }
        }
        
        #the segments are added to the plot
        bys_plot<-bys_plot+
          geom_segment(
            data=active.lines,
            aes(x=x,y=y,xend=xend,yend=yend),
            colour='grey',linetype='dashed',size=1
          )
        #a shorter list is created to eliminate duplicates (ie where the x points
        #start and the ys end)
        labs<-unique(active.lines[c('y','Sample','Ligand','EC50')])
        
        #for each row of this list....
        for (lab in 1:dim(labs)[1]){
          bys_plot<-bys_plot+
            #... the point on the y axis the line is drawn is found, and the label
            #is put there
            annotate('text',x=-10,y=as.double(labs$y[lab]),
                     #the EC50 value is used, rounded to 2dp.
                     label=paste0("EC50: ",round(labs$EC50[lab],2),"μM"),
                     hjust=0,
                     vjust=-1)
        }
      }#end of add ec50 labels
      #save the graph
 #     assign((paste0(Experiment,'_Plot')),bys_plot,envir = .GlobalEnv)
      output[["Plot"]]<-bys_plot
    }#end of graphs
  }#end of ec50 calculations
  return(output)
}#end of function



