###Function to import plate setups
# 
# output<-list()
# output[["One"]]<-1
# output[[2]]<-c("one","two")
# output[[2]]
# output

#08
GProtein08<-list(
  GProtein=TRUE,
  Samples=c("D2S","D2L","D4.4","D4.7","D2S","D2L","D4.4","D4.7"),
  Ligands=data.frame(rep("Dopamine",8),rep("Aripiprizole",8),
                     rep("Dopamine",8),rep("Aripiprizole",8)),
  Concentrations=c(300,100,30,10,3,1,0.3,0.1,0.03,0.01,0.003,0)
)
#21
GProtein21<-list(
  GProtein=TRUE,
  Samples=c("D2S","D2L","D3","D4.4","D4.7","pcDNA","3ug","6ug"),
  Ligands="Dopamine"
)
#22
GProtein22<-list(
  GProtein=TRUE,
  Samples=c("D2S","D2L","D3","D4.4","D4.7","pcDNA","3ug","6ug"),
  Ligands="Dopamine"
)
#23 there is an issue with this one
Lig23<-c("Dopamine",
         "Quinpirole",
         "MLS-1547",
         "UNC-9994",
         "UNC-9975",
         "MLS-1547",
         "UNC-9994",
         "UNC-9975")
Ligs2627<-data.frame(Lig23,Lig23,Lig23,Lig23)
Lig23<-data.frame(Lig23,Lig23,Lig23,Lig23)
Lig23[6,4]<-"Dopamine"
Lig23[7,3:4]<-"Quinpirole"
Lig23[8,3:4]<-c("UNC-9994","UNC-9975")

GProtein25<-list(
  GProtein=TRUE,
  Samples=c("D2L","D2L","D2L","D2L","D2L","pcDNA","pcDNA","pcDNA"),
  Ligands=Lig23,
  Concentrations=c(30,10,1,3,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0)
)

#24 dodgy??
GProtein24_Missing<-"Experiment didn't work"

#25
GProtein25<-list(
  GProtein=TRUE,
  Samples=c("D2L","D2L","D2L","D2L","D2L","pcDNA","pcDNA","pcDNA"),
  Ligands=c("Dopamine",
            "Quinpirole",
            "MLS-1547",
            "UNC-9994",
            "UNC-9975",
            "MLS-1547",
            "UNC-9994",
            "UNC-9975")
)
#26
#data inherited from 23
Ligs2627[6,3:4]<-"Dopamine"
Ligs2627[7,3:4]<-"Quinpirole"
Ligs2627[8,3:4]<-c("UNC-9994","UNC-9975")

GProtein26<-list(
  GProtein=TRUE,
  Samples=c("D2L","D2L","D2L","D2L","D2L","pcDNA","pcDNA","pcDNA"),
  Ligands=Ligs2627
)

#27-29
GProtein29<-GProtein28<-GProtein27<-GProtein26

#30
GProtein30_Warning<-"GProtein30 compared dopamine prepared at different time points before the experiment"
GProtein30<-list(
  GProtein=TRUE,
  Samples=c("D2L","pcDNA","D2L","pcDNA","D2L","pcDNA","D2L","pcDNA"),
  Ligands=c("2-Hours",
            "2-Hours",
            "1-Hours",
            "1-Hours",
            "0.5-Hours",
            "0.5-Hours",
            "0-Hours",
            "0-Hours")
)


#31
GProtein31<-list(
  GProtein=TRUE,
  Samples=c("D2L","D2L","D2L","D2L","D2L","FAULT","FAULT","FAULT"),
  Ligands=c("Dopamine",
            "Quinpirole",
            "MLS-1547",
            "UNC-9994",
            "UNC-9975",
            "FAULT",
            "FAULT",
            "FAULT")
)
#32
GProtein32_Warning<-"GProtein32 compared different G-Protein amplification constructs"

GProtein32<-list(
  GProtein=TRUE,
  Samples=c("D2L"),
  Ligands=c("Quinpirole",
            "Quinpirole",
            "Dopamine",
            "Dopamine",
            "Quinpirole",
            "Quinpirole",
            "Dopamine",
            "Dopamine")
)
#33
GProtein33<-list(
  GProtein=TRUE,
  Samples=c("pcDNA","D2L","pcDNA","D2L","pcDNA","D2L","pcDNA","D2L"),
  Ligands=c("Quinpirole",
            "Quinpirole",
            "Dopamine",
            "Dopamine",
            "Quinpirole",
            "Quinpirole",
            "Dopamine",
            "Dopamine")
)
#34
GProtein34<-list(
  GProtein=TRUE,
  Samples=data.frame(rep("D2L",8),rep("pcDNA",8)),
  Ligands=c("Quinpirole",
            "Dopamine",
            "Aripiprizole",
            "Quinpirole",
            "Dopamine",
            "Aripiprizole",
            "Dopamine",
            "Aripiprizole")
)

Bystander34<-list(
  GProtein=FALSE,
  Samples=data.frame(rep("D2L",8),rep("pcDNA",8)),
  Ligands=c("Quinpirole",
            "Dopamine",
            "Aripiprizole",
            "Quinpirole",
            "Dopamine",
            "Aripiprizole",
            "Dopamine",
            "Aripiprizole")
)

#35
GProtein35<-list(
  GProtein=TRUE,
  Samples="D2L",
  Ligands=c("Quinpirole","Dopamine","Aripiprizole","Quinpirole","Dopamine","Aripiprizole","Dopamine","Aripiprizole")
)

Bystander35<-list(
  GProtein=FALSE,
  Samples="D2L",
  Ligands=c("Quinpirole","Dopamine","Aripiprizole","Quinpirole","Dopamine","Aripiprizole","Dopamine","Aripiprizole")
)

#36
GProtein36<-list(
  GProtein=TRUE,
  Samples="D2L",
  Ligands=c("Quinpirole","Dopamine","Aripiprizole","Quinpirole","Dopamine","Aripiprizole","Quinpirole","Dopamine")
)

Bystander36<-list(
  GProtein=FALSE,
  Samples="D2L",
  Ligands=c("Quinpirole",
            "Dopamine",
            "Aripiprizole",
            "Quinpirole",
            "Dopamine",
            "Aripiprizole",
            "Quinpirole",
            "Dopamine"))


#38
Bystander38<-list(
  GProtein=FALSE,
  Samples="D2L",
  Ligands=c("Dopamine",
            "Quinpirole",
            "Aripiprizole",
            "Dopamine",
            "Quinpirole",
            "Aripiprizole",
            "Dopamine",
            "Aripiprizole")
)

#39
GProtein39<-list(
  GProtein=TRUE,
  Samples="D2L",
  Ligands=c("Dopamine",
            "Quinpirole",
            "Aripiprizole",
            "Dopamine",
            "Quinpirole",
            "Aripiprizole",
            "Dopamine",
            "Aripiprizole")
)

Bystander39<-list(
  GProtein=FALSE,
  Samples="D2L",
  Ligands=c("Dopamine",
            "Quinpirole",
            "Aripiprizole",
            "Dopamine",
            "Quinpirole",
            "Aripiprizole",
            "Dopamine",
            "Aripiprizole")
)

#40
GProtein40<-list(
  GProtein=TRUE,
  Samples="D2L",
  Ligands=c("Dopamine",
            "MLS-1547",
            "UNC-9994",
            "UNC-9975")
)

Bystander40<-list(
  GProtein=FALSE,
  Samples="D2L",
  Ligands=c("Dopamine",
            "MLS-1547",
            "UNC-9994",
            "UNC-9975")
)

#41
GProtein41<-list(
  GProtein=TRUE,
  Samples="D2L",
  Ligands=c("Dopamine",
            "MLS-1547",
            "UNC-9994",
            "UNC-9975")
)

Bystander41<-list(
  GProtein=FALSE,
  Samples="D2L",
  Ligands=c("Dopamine",
            "MLS-1547",
            "UNC-9994",
            "UNC-9975")
)

#42
GProtein42<-list(
  GProtein=TRUE,
  Samples="D2L",
  Ligands=c("MLS-1547",
            "UNC-9975",
            "Dopamine",
            "UNC-9994")
)

Bystander42<-list(
  GProtein=FALSE,
  Samples="D2L",
  Ligands=c("MLS-1547",
            "UNC-9975",
            "Dopamine",
            "UNC-9994")
)

#43
GProtein43<-GProtein42

Bystander43<-Bystander42

#44
GProtein44<-list(
  GProtein=TRUE,
  Samples="D2L",
  Ligands=c("UNC-9994",
            "Dopamine",
            "UNC-9975",
            "MLS-1547")
)

Bystander44<-list(
  GProtein=FALSE,
  Samples="D2L",
  Ligands=c("UNC-9994",
            "Dopamine",
            "UNC-9975",
            "MLS-1547")
)

plate.layouts<-list(GProtein08,
                    GProtein21,
                    GProtein22,
                    GProtein25,
                    GProtein26,
                    GProtein27,
                    GProtein28,
                    GProtein29,
                    GProtein30,
                    GProtein31,
                    GProtein32,
                    GProtein33,
                    GProtein34,
                    Bystander34,
                    GProtein35,
                    Bystander35,
                    GProtein36,
                    Bystander36,
                    Bystander38,
                    GProtein39,
                    Bystander39,
                    GProtein40,
                    Bystander40,
                    GProtein41,
                    Bystander41,
                    GProtein42,
                    Bystander42,
                    GProtein43,
                    Bystander43,
                    GProtein44,
                    Bystander44)

names(plate.layouts)<-c("GProtein08",
                        "GProtein21",
                        "GProtein22",
                        "GProtein25",
                        "GProtein26",
                        "GProtein27",
                        "GProtein28",
                        "GProtein29",
                        "GProtein30",
                        "GProtein31",
                        "GProtein32",
                        "GProtein33",
                        "GProtein34",
                        "Bystander34",
                        "GProtein35",
                        "Bystander35",
                        "GProtein36",
                        "Bystander36",
                        "Bystander38",
                        "GProtein39",
                        "Bystander39",
                        "GProtein40",
                        "Bystander40",
                        "GProtein41",
                        "Bystander41",
                        "GProtein42",
                        "Bystander42",
                        "GProtein43",
                        "Bystander43",
                        "GProtein44",
                        "Bystander44")


# namelist<-Filter(function(x) !any(grepl("Warning", x)), ls(pattern="GProtein"))
# namelist<-Filter(function(x) !any(grepl("Warning", x)), ls(pattern="Bystander"))
# 
# names(plate.layouts)<-Filter(function(x) !any(grepl("Missing", x)), namelist)
