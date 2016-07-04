
#' A validation function of ProBiS .lig files
#'
#' Use this function to check whether your files are of correct format.
#' @param Your .lig file - a ProBiS algorithm output.
#' @keywords chackFile
#' @export
#' @examples
#' ValidateProBiS(imported_file_name e.g. mol1.lig)

ValidateProBiS <- function (file1){
  interaction_file <- unique(file1)
  flen<-nrow(interaction_file)
  by(interaction_file, 1:flen, function(row) 
      if (length(row[3])>5){
      print(paste("Found suspicious entry! ", row))
        }
#       if (row[5]%%1!=0){
#         print(paste("invalid locus entry", row))
#       }
    )
  print("PDB codes validated!")
}



#' Function to display predicted area of interaction
#'
#' Use this function to obtain all loci within predicted PI
#' @param Your ProBiS interaction file
#' @keywords ProBiSregion
#' @export
#' @examples
#' predictedArea(imported_file_name e.g. mol1.lig)

predictedArea<- function (interaction_file){
  predictionRange <- unique(interaction_file$V5)
  print(predictionRange)
  
}

#' Function to display most common interactions present
#'
#' Use this function to find loci within most interactions
#' @param Your ProBiS interaction file
#' @keywords ProBiSmaxLoci
#' @export
#' @examples
#' mostCommonInteraction(imported_file_name e.g. mol1.lig, 15)

mostCommonInteraction<- function (interaction_file, topn){
  print("trying to generate table of interaction partners...")
  mci<-tail(sort((table(interaction_file$V8))),topn)
  print(mci)
}

#' Function to display single loci within most interactions
#'
#' Use this function to obtain most present loci in interaction set
#' @param Your ProBiS interaction file
#' @keywords ProBiSmaxLoci
#' @export
#' @examples
#' maxInteraction(imported_file_name e.g. mol1.lig)

maxInteraction <- function (interaction_file){
  maxval<-max(sort(table(interaction_file$V5)))
  print ("Maximum number of detected interactions on one loci>")
  print (maxval)
}

#' Function to plot interaction intensities for specific loci
#'
#' Use this function to plot frequeency of loci within interaction
#' @param Your ProBiS interaction file and number of desired loci
#' @keywords ProBiSplot
#' @export
#' @examples
#' plotbc(imported_file_name e.g. mol1.lig, 10)

plotbc <- function(ifile, top){
maxval<-max(sort(table(ifile$V5)))
colfunc <- colorRampPalette(c("black", "green"))
barplot(tail(sort(table(ifile$V5)),top), ylim=c(0, maxval+5), 
        col= colfunc(top), 
        main="Number of predicted interactions at specific site", 
        las=3,
        xlab=c("Locus on whole protein sequence"),
        ylab=c("Number of predicted interactions"))
text(top/2,maxval, paste("Loci with most interactions> ",maxval))
}

#' Function to display unique interactions present according to PDBs
#'
#' Use this function to obtain uniqe PDB to PDB interactions
#' @param Your ProBiS interaction file
#' @keywords ProBiSregion
#' @export
#' @examples
#' interactionList(imported_file_name e.g. mol1.lig)

interactionList <- function (interaction_file){
  pdbs<-paste(interaction_file$V3, interaction_file$V8, sep="-")
  print("unique interactions detected>")
  print(unique(pdbs))
}

#' Entire summary of your .lig file.
#'
#' Use this function obtain quick insight into your .lig file
#' @param Your ProBiS interaction file
#' @keywords ProBiSsummary
#' @export
#' @examples
#' probisSummary(imported_file_name e.g. mol1.lig)

probisSummary<- function(interaction_file){
  interaction_file<- unique(interaction_file)
  predictionRange <- unique(interaction_file$V5)
  print("Range of loci in predicted interaction sites> ")
  print(predictionRange)
  print("Most common interactions> ")
  mci<-sort(table(interaction_file$V5))
  print(mci)
  maxval<-max(sort(table(interaction_file$V5)))
  print ("Maximum number of detected interactions on one loci>")
  print (maxval)
  pdbs<-paste(interaction_file$V3, interaction_file$V8, sep="-")
  print("List of unique interactions detected>")
  print(unique(pdbs))
  par(mfrow=c(2,2))
  sizes<- c(length(unique(interaction_file$V5))/10,
            length(unique(interaction_file$V5))/5,
            length(unique(interaction_file$V5))*0.5,
            length(unique(interaction_file$V5)))
  for (h in sizes){
    plotbc(interaction_file, h)
  }

}

#' Find PDB specific interactions
#'
#' Use this function to obtain all PDB specific interactions
#' @param Your ProBiS interaction file and PDB name
#' @keywords ProBiSregion
#' @export
#' @examples
#' findbyPDB(imported_file_name e.g. mol1.lig, '3epoA')

findbyPDB <- function(dbname, lname){
  ldataset <- data.frame(dbname[which(dbname[,3]==lname),])
  View(unique(ldataset))
  return(unique(dbname[which(dbname[,3]==lname),]))
  
}

#' Obtain raw data to paste into PyMol for plotting and selection
#'
#' Paste output to selection function in PyMol
#' @param Your ProBiS interaction file
#' @keywords ProBiSPymoL
#' @export
#' @examples
#' pymolRegion(imported_file_name e.g. mol1.lig)

pymolRegion <- function (ligandFile){
  f<-noquote(paste0(as.vector(predictedArea(ligandFile)), sep="+"))
  f<-gsub(",","",toString(f))
  f<-gsub(" ","",f)
  return(f)
}
