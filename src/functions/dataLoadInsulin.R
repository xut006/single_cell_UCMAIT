#### load insulin data ####
dataLoadInsulin <- function(){
  
  ## set data directory
  dataDir <- paste(baseDir, "data/", sep="")
  
  ### read in data
  ## age: 7 weeks (6 weeks)
  ## tissue: islets
  ## probes: 12-20 & 13-21
  ## mice: 1 & 2 (1)
  ctTable1 <- read.csv(paste(dataDir, "Islets_Insulin1_sort090415.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable1$mouse <- "1"
  ## mice: 3 & 4 (2)
  ctTable2 <- read.csv(paste(dataDir, "Islets_Insulin2_sort091015.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable2$mouse <- "2"
  
  ## age: 6 & 12 weeks
  ## tissue: Spleen & pLN
  ## probes: 12-20 & 13-21
  ## mice: A & B at 6 weeks = (3 & 4)
  ## mice: A & B at 12 weeks = (5 & 6)
  ctTable3 <- read.csv(paste(dataDir, "12-20A_6weeks.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable3$mouse <- "3"
  ctTable4 <- read.csv(paste(dataDir, "12-20A_12weeks.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable4$mouse <- "5"
  ctTable5 <- read.csv(paste(dataDir, "12-20B_6weeks.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable5$mouse <- "4"
  ctTable6 <- read.csv(paste(dataDir, "12-20B_12weeks.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable6$mouse <- "6"
  ctTable7 <- read.csv(paste(dataDir, "13-21A_6weeks.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable7$mouse <- "3"
  ctTable8 <- read.csv(paste(dataDir, "13-21A_12weeks.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable8$mouse <- "5"
  ctTable9 <- read.csv(paste(dataDir, "13-21B_6weeks.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable9$mouse <- "4"
  ctTable10 <- read.csv(paste(dataDir, "13-21B_12weeks.csv.trimmed", sep=""), 
                        header=TRUE, stringsAsFactors=FALSE)
  ctTable10$mouse <- "6"
  
  ## age: 7 weeks (6 weeks)
  ## tissue: islets & pLN
  ## probes: 12-20 & 13-21
  ## mice: 2 (7)
  ctTable11 <- read.csv(paste(dataDir, "insulin_7weeks_mouse2.csv.trimmed", sep=""), 
                        header=TRUE, stringsAsFactors=FALSE)
  ctTable11$mouse <- "7"
  ## mice: 1, 3 & 4 (8, 9, 10)
  ctTable12 <- read.csv(paste(dataDir, "314_insulin_7weeks.csv.trimmed", sep=""), 
                        header=TRUE, stringsAsFactors=FALSE)
  ctTable12 <- ctTable12[-(grep("empty", ctTable12$Name)), ]
  ctTable12$mouse <- substr(ctTable12$Name, 6, 6)
  ctTable12[which(ctTable12$mouse == "1"), "mouse"] <- "8"
  ctTable12[which(ctTable12$mouse == "3"), "mouse"] <- "9"
  ctTable12[which(ctTable12$mouse == "4"), "mouse"] <- "10"
  
  ## age: 13 weeks (12 weeks)
  ## tissue: islets & pLN
  ## probes: 12-20 & 13-21
  ## mice: 2, 3 & 4 (11, 12 & 13)
  ctTable13 <- read.csv(paste(dataDir, "Insulin13weeksMouse2.csv.trimmed", sep=""), 
                        header=TRUE, stringsAsFactors=FALSE)
  ctTable13$mouse <- "11"
  ctTable14 <- read.csv(paste(dataDir, "Insulin13weeksMouse3.csv.trimmed", sep=""), 
                        header=TRUE, stringsAsFactors=FALSE)
  ctTable14$mouse <- "12"
  ctTable15 <- read.csv(paste(dataDir, "Insulin13weeksMouse4.csv.trimmed", sep=""), 
                        header=TRUE, stringsAsFactors=FALSE)
  ctTable15$mouse <- "13"
  
  ## age: 7 weeks (6 weeks)
  ## tissue: islets & pLN
  ## probes: 12-20 & 13-21
  ## mice: 7, 10 & 11 (14, 15 & 16)
  ctTable16 <- read.csv(paste(dataDir, "insulin_7weeks_Mouse7_LN_and_islets.csv.trimmed", sep=""), 
                        header=TRUE, stringsAsFactors=FALSE)
  ctTable16$mouse <- "14"
  ctTable17 <- read.csv(paste(dataDir, "insulin_7weeks_Mouse10_LN_and_islets.csv.trimmed", sep=""), 
                        header=TRUE, stringsAsFactors=FALSE)
  ctTable17$mouse <- "15"
  ctTable18 <- read.csv(paste(dataDir, "insulin_7weeks_Mouse11_LN_and_islets.csv.trimmed", sep=""), 
                        header=TRUE, stringsAsFactors=FALSE)
  ctTable18$mouse <- "16"
  
  ctTable19 <- read.csv(paste(dataDir, "8Wp1220M3LNSIsletsPMC.csv.trimmed", sep=""), 
                        header=TRUE, stringsAsFactors=FALSE)
  ctTable19$mouse <- "17"
  
  ## everything
  ctTable <- rbind(ctTable1, ctTable2, ctTable3, ctTable4, ctTable5, ctTable6, ctTable7, ctTable8, ctTable9, ctTable10, ctTable11, ctTable12, ctTable13, ctTable14, ctTable15, ctTable16, ctTable17, ctTable18, ctTable19)
  
  # ctTable <- rbind(ctTable1, ctTable2, ctTable3, ctTable4, ctTable5, ctTable6, ctTable7, ctTable8, ctTable9, ctTable10, ctTable11, ctTable12, ctTable13, ctTable14, ctTable15, ctTable16)
  ## pretty clean
  # ctTable <- rbind(ctTable1, ctTable2, ctTable11, ctTable12, ctTable13, ctTable14, ctTable15, ctTable16)
  
  ## change column names
  names(ctTable)[c(2, 5, 7)] <- c("cellType", "gene", "ct")
  
  ## create cellSource column
  ctTable$cellSource <- NA
  ctTable[grepl("S", ctTable$cellType, fixed=TRUE), "cellSource"] <- "S"
  ctTable[grepl("LN", ctTable$cellType, fixed=TRUE), "cellSource"] <- "LN"
  ctTable[grepl("i", ctTable$cellType, fixed=TRUE), "cellSource"] <- "I"
  ctTable[grepl("I", ctTable$cellType, fixed=TRUE), "cellSource"] <- "I"
  ctTable[grepl("b", ctTable$cellType, fixed=TRUE), "cellSource"] <- "PB"
  ctTable[grepl("empty", ctTable$cellType, fixed=TRUE), "cellSource"] <- "controlIslet"
  ctTable[grepl("zero", ctTable$cellType, fixed=TRUE), "cellSource"] <- "controlIslet"
  
  ctTable$probe <- substr(ctTable$cellType, 1, 2)
  ctTable[grepl("p12M", ctTable$cellType, fixed=TRUE), "probe"] <- "12"
  
#   ctTable$mouse <- NA
#   ctTable[which(ctTable$cellSource=="I"), "mouse"] <- 
#     substr(ctTable[which(ctTable$cellSource=="I"), "cellType"], 6, 6)
#   ctTable[which(ctTable$cellSource=="S" | ctTable$cellSource=="LN"), "mouse"] <-
#     substr(ctTable[which(ctTable$cellSource=="S" | ctTable$cellSource=="LN"), "cellType"], 3, 3)
#   ctTable[which(ctTable$mouse=="-"), "mouse"] <-
#     substr(ctTable[which(ctTable$mouse=="-"), "cellType"], 6, 6)
#   ctTable[grepl("M2", ctTable$cellType, fixed=TRUE), "mouse"] <- "2"
#   ctTable[grepl("M3", ctTable$cellType, fixed=TRUE), "mouse"] <- "3"
#   ctTable[grepl("M4", ctTable$cellType, fixed=TRUE), "mouse"] <- "4"
#   ctTable[grepl("M7", ctTable$cellType, fixed=TRUE), "mouse"] <- "7"
#   ctTable[grepl("M10", ctTable$cellType, fixed=TRUE), "mouse"] <- "10"
#   ctTable[grepl("M11", ctTable$cellType, fixed=TRUE), "mouse"] <- "11"
  
  ctTable$age <- NA
  ctTable[which(ctTable$cellSource=="I"), "age"] <- "6"
  ctTable[grepl("A6", ctTable$cellType, fixed=TRUE), "age"] <- "6"
  ctTable[grepl("B6", ctTable$cellType, fixed=TRUE), "age"] <- "6"
  ctTable[grepl("A12", ctTable$cellType, fixed=TRUE), "age"] <- "12"
  ctTable[grepl("B12", ctTable$cellType, fixed=TRUE), "age"] <- "12"
  ctTable[grepl("W7", ctTable$cellType, fixed=TRUE), "age"] <- "6"
  ctTable[grepl("W13", ctTable$cellType, fixed=TRUE), "age"] <- "12"
  ctTable[grepl("W8", ctTable$cellType, fixed=TRUE), "age"] <- "6"
  ctTable[which(is.na(ctTable$age)), "age"] <- "6"
  
  ## blood data
  ctTableBlood1 <- read.csv(paste(dataDir, "blood7W12-20M1andM2.csv.trimmed", sep=""), 
                        header=TRUE, stringsAsFactors=FALSE)
  names(ctTableBlood1)[c(2, 5, 7)] <- c("cellType", "gene", "ct")
  ctTableBlood1$mouse <- NA
  ctTableBlood1[grepl("M1", ctTableBlood1$cellType, fixed=TRUE), "mouse"] <- "12"
  ctTableBlood1[grepl("M2", ctTableBlood1$cellType, fixed=TRUE), "mouse"] <- "13"
  
  ctTableBlood2 <- read.csv(paste(dataDir, "blood8W12-20M1andM2.csv.trimmed", sep=""), 
                        header=TRUE, stringsAsFactors=FALSE)
  ctTableBlood3 <- read.csv(paste(dataDir, "blood8W12-20M3M4.csv.trimmed", sep=""), 
                        header=TRUE, stringsAsFactors=FALSE)
  
  ctTableBlood2and3 <- rbind(ctTableBlood2, ctTableBlood3)
  names(ctTableBlood2and3)[c(2, 5, 7)] <- c("cellType", "gene", "ct")
  ctTableBlood2and3$mouse <- NA
  ctTableBlood2and3[grepl("M1", ctTableBlood2and3$cellType, fixed=TRUE), "mouse"] <- "14"
  ctTableBlood2and3[grepl("M2", ctTableBlood2and3$cellType, fixed=TRUE), "mouse"] <- "15"
  ctTableBlood2and3[grepl("M3", ctTableBlood2and3$cellType, fixed=TRUE), "mouse"] <- "16"
  ctTableBlood2and3[grepl("M4", ctTableBlood2and3$cellType, fixed=TRUE), "mouse"] <- "17"
  
  ctTableBlood <- rbind(ctTableBlood1, ctTableBlood2and3)
  
  ctTableBlood$probe <- substr(ctTableBlood$cellType, 1, 2)
  ctTableBlood$age <- "6"
  ctTableBlood$cellSource <- "PB"
  ctTableBlood[grepl("empty", ctTableBlood$cellType, fixed=TRUE), "cellSource"] <- "controlBlood"
  ctTableBlood[grepl("zero", ctTableBlood$cellType, fixed=TRUE), "cellSource"] <- "controlBlood"
  
  ctTableBlood <- ctTableBlood[, c(names(ctTableBlood)[1:14], "cellSource", "probe", "mouse", "age")]
  
  ctTableFull <- rbind(ctTable, ctTableBlood)
  
  return(ctTableFull)
  
}
