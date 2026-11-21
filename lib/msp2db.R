## function for converting MSP files to DB files for RTLS analyses
# author: Chris McGann

#Builds DB
library(BiocParallel)
library(RSQLite)
library(dplyr)
library(compiler)
library(vroom)
library(data.table)


# main function of shiny app
DBbuilder<- function(Library, FragmentationMode, MassAnalyzer, CollisionEnergy, CompoundClass,
                     Filter, DBoutput, topX, cutoff, TMTPro, massOffset, IonTypes,deltaFragment, oldMod, newMod) {
#Get file type from input file
  dbkey_root <- Sys.getenv("DBKEY_HOME", unset = path.expand("~/Repos/MSPtoDB"))
fileType<- (Library$name[1])
if(grepl("msp", tolower(fileType))) {
  fileType <- "Prosit"
} else if(grepl("sptxt", tolower(fileType))) {
  fileType <- "SpectraST"
} else if(grepl("blib", tolower(fileType))) {
  fileType<-"Skyline"
}
  # Complete-record streaming avoids PSOCK worker copies and ambiguous line offsets.
  # Keep a maximum of one incomplete final spectrum between reads.
  if (fileType != "Skyline") {
    if (fileType == "SpectraST") {
      source(file.path(dbkey_root, "lib", "SpXLibraryParser.R"), local = TRUE)
      parser <- SpXLibraryParser
    } else {
      source(file.path(dbkey_root, "lib", "LibraryParserv2.R"), local = TRUE)
      parser <- LibraryParser
    }
    previous_threads <- data.table::getDTthreads()
    data.table::setDTthreads(1L)
    on.exit(data.table::setDTthreads(previous_threads), add = TRUE)
    results <- list()
    spectra_done <- 0L
    parse_complete <- function(lines) {
      starts <- grep("^Name: ", lines)
      if (!length(starts)) stop("Input block has no spectral entries")
      lines <- lines[seq.int(starts[1L], length(lines))]
      answer <- parser(
        Library = lines, FragmentationMode = FragmentationMode,
        MassAnalyzer = MassAnalyzer, CollisionEnergy = CollisionEnergy,
        CompoundClassArg = CompoundClass, Filter = Filter, TMTPro = TMTPro,
        Source = fileType, topX = topX, cutoff = cutoff, massOffset = massOffset,
        IonTypes = IonTypes, deltaFragment = deltaFragment,
        oldMod = oldMod, newMod = newMod)
      if (nrow(answer) != length(starts))
        stop("Parser returned a different spectrum count than the input")
      spectra_done <<- spectra_done + nrow(answer)
      message("Parsed ", spectra_done, " spectra")
      answer
    }
    for (library_path in Library$datapath) {
      input_connection <- file(library_path, open = "rt")
      pending <- character()
      tryCatch({
        repeat {
          block <- readLines(input_connection, n = 100000L, warn = FALSE)
          if (!length(block)) break
          pending <- c(pending, block)
          starts <- grep("^Name: ", pending)
          if (length(starts) > 1L) {
            final_start <- tail(starts, 1L)
            results[[length(results) + 1L]] <-
              parse_complete(pending[seq_len(final_start - 1L)])
            pending <- pending[seq.int(final_start, length(pending))]
          }
        }
        if (length(grep("^Name: ", pending)))
          results[[length(results) + 1L]] <- parse_complete(pending)
      }, finally = close(input_connection))
    }
    if (!length(results)) stop("No spectra found in the input file(s)")
    resultsTable <- data.table::rbindlist(results, use.names = TRUE, fill = TRUE)
    rm(results)
  }

  if(fileType!="Skyline"){
 # #Builds tables
    MasterCompoundTable <- data.table(
      CompoundId = seq(1, length(resultsTable$Names), by=1),
      Formula = unlist(resultsTable$Formula),
      Name = unlist(resultsTable$Names),
      Synonyms = "",
      Tag = unlist(resultsTable$Tags),
      Sequence = gsub('.{2}$', '',resultsTable$Names),
      CASId = "",
      ChemSpiderId= "",
      HMDBId= "",
      KEGGId= "",
      PubChemId= "",
      Structure= "",
      mzCloudId= as.integer(NA),
      CompoundClass= resultsTable$CompoundClass,
      SmilesDescription= "",
      InChiKey= "")


    MasterSpectrumTable <- data.table(
      SpectrumId = seq(1, length(resultsTable$Names), by=1),
      CompoundId = seq(1, length(resultsTable$Names), by=1),
      mzCloudURL = "",
      ScanFilter = "",
      RetentionTime = resultsTable$iRT,
      ScanNumber = 0,#
      PrecursorMass = resultsTable$PrecursorMasses,
      NeutralMass= 0,
      CollisionEnergy= resultsTable$CollisionEnergy,
      Polarity= "+",
      FragmentationMode= resultsTable$FragmentationMode,
      IonizationMode= "ESI",
      MassAnalyzer= MassAnalyzer,
      InstrumentName= "",
      InstrumentOperator= "",
      RawFileURL= "",
      blobMass= I(resultsTable$blobMass),
      blobIntensity= I(resultsTable$blobInt),
      blobAccuracy="",
      blobResolution="",
      blobNoises="",
      blobFlags="",
      blobTopPeaks="",
      Version="",
      CreationDate=NA,
      Curator="",
      CurationType="",
      PrecursorIonType="",
      Acession="")

    HeaderTable <- data.frame(
      version = 5,
      CreationDate = NA,
      LastModifiedDate = NA,
      Description = NA,
      Company = NA,
      ReadOnly = NA,
      UserAccess = NA,
      PartialEdits= NA)

    MaintenanceTable <- data.frame(
      CreationDate = gsub("-"," ",Sys.Date()),
      NoofCompoundsModified=NA,
      Description=NA
    )
    
    ##Writing new DB
  conn4 <- dbConnect(SQLite(),DBoutput)
  dbWriteTable(conn4,"CompoundTable", MasterCompoundTable, overwrite = TRUE, append = F, field.types = c(CompoundId="INTEGER PRIMARY KEY", Synonyms="BLOB_TEXT", Structure="BLOB_TEXT"))
  dbWriteTable(conn4, "SpectrumTable", MasterSpectrumTable, overwrite = TRUE, append = F, field.types = c( SpectrumId = "INTEGER PRIMARY KEY", CompoundId ="INTEGER REFERENCES [CompoundTable]", blobMass="BLOB", blobIntensity="BLOB", RetentionTime = "DOUBLE", ScanNumber ="INTEGER", PrecursorMass="DOUBLE", NeutralMass="DOUBLE", blobAccuracy = "BLOB",blobFlags = "BLOB",blobResolution = "BLOB",blobNoises = "BLOB",blobTopPeaks = "BLOB", Version = "INTEGER"))
  dbWriteTable(conn4, "HeaderTable", HeaderTable, overwrite = TRUE, append = F, field.types = c(version = "INTEGER NOT NULL DEFAULT 0", CreationDate="TEXT",LastModifiedDate="TEXT", Description="TEXT", Company="TEXT", ReadOnly = "BOOL", UserAccess="TEXT", PartialEdits = "BOOL"))
  dbWriteTable(conn4, "MaintenanceTable", MaintenanceTable, overwrite = TRUE, append = F, field.types = c(CreationDate="TEXT",NoofCompoundsModified = "INTEGER", Description="TEXT"))
  dbDisconnect(conn4)
  return(invisible(DBoutput))
  }

   if( fileType == "Skyline") {
    source(file.path(dbkey_root, "lib", "SkylineConvert.R"))
     LibraryPath<-Library$datapath
SkylineConvert(x=LibraryPath,CollisionEnergy=CollisionEnergy,CompoundClassArg=CompoundClass,FragmentationMode=FragmentationMode,MassAnalyzer=MassAnalyzer,topX=topX,cutoff=cutoff, Filter=Filter,massOffset=massOffset, DBoutput=DBoutput)

  }
}


