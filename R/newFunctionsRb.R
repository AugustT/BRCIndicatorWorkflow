##### create rdata files in TrendSummaries format from Charlie's posterior samples csv files #####

## This function is a bit bespoke because I just needed the csv files in a format I could use. Will probably
## just use Gary's sampBind which works on the sparta outputs (.rdata files) for the generic workflow

sampBind2 <- function(group, inPath, outPath, region, minYear, maxYear) {
  
  files <- list.files(inPath,
                         full.names=T)[which(
                           grepl(pattern=paste0("_", group, "_"), 
                                 x=list.files(inPath, full.names = T)))] ## list all files for group
  
  "%!in%" = Negate("%in%")  
  
  filterCountry <- function(file) { ## function to extract posteriors for the chosen country and format for TrendSummaries
    
    dat <- read.csv(file)
    
    if (region == "UK" & "UK" %!in% dat$Region) {
      
      dat <- dat[which(dat$Region == "GB"),]
      
    } else {
      
      dat <- dat[which(dat$Region == region),]
      
    }
    
    dat <- dat[, -c(1,3)]
    
    dat <- dat[, c(3:ncol(dat), 2,1)]
    
  }
  
  stacked_samps <- purrr::map_df(.x = files, .f= filterCountry) ## apply function over all files and rbind the outputs
  
  yrCols <- paste0("year_", 1970:(1970+ (ncol(stacked_samps)-3)))
  
  colnames(stacked_samps) <- c(yrCols, "iteration","species")
  
  allCols <- paste0("year_", minYear:maxYear)

  for (i in allCols) { # add columns of NAs for years with no data (older models only go up to 2015)
    
    if (i %!in% colnames(stacked_samps)) {
      
      if (i == "year_2016") {
      stacked_samps <- tibble::add_column(stacked_samps,  "year_2016" = rep(NA, nrow(stacked_samps)), .before = "iteration")
      } else if (i == "year_2017") {
        stacked_samps <- tibble::add_column(stacked_samps,  "year_2017" = rep(NA, nrow(stacked_samps)), .before = "iteration")  
      } else if (i == "year_2018") {
        stacked_samps <- tibble::add_column(stacked_samps,  "year_2018" = rep(NA, nrow(stacked_samps)), .before = "iteration")
      }
      
    }
    
  }
  
  save(stacked_samps, file=paste0(outPath, group, "stacked.rdata"))
  
  return(stacked_samps)
  
}





















## Function to extract metadata from sparta input files 

extractMeta <- function(inPath, group) {
  
  file <- list.files(inPath,
                      full.names=T)[which(
                        grepl(pattern=group, 
                              x=list.files(inPath)))]

  if (length(file) > 0) {
    
    load(file)
    
    colnames(taxa_data) <- toupper(colnames(taxa_data))
    
  getMeta <- function(spp) {
    
    nRec <- length(taxa_data$CONCEPT[taxa_data$CONCEPT == spp]) # total number of records
    
    first <- min(taxa_data$YEAR[taxa_data$CONCEPT == spp]) # year of first record
    
    last <- max(taxa_data$YEAR[taxa_data$CONCEPT == spp]) # year of last record

    yrs <- sort(unique(taxa_data$YEAR[taxa_data$CONCEPT == spp]), decreasing = F) # used to calc largest gap in data

    gaps <- NULL 
    
    if (length(yrs) > 1) { # some species only have data for one year
      
      for (i in (1:length(yrs) - 1)) {
      
        gaps <- c(gaps, yrs[i+1] - yrs[i])
      
      }
    }
    
    if (!is.null(gaps)) { # if more than one year of data
      
      gap <- max(gaps) # find the biggest gap (years) in the data
      
    } else {
      
      gap <- 1 # does not matter what the gap is
    }
    
    return(data.frame(spp, nRec, first, last, gap))
    
  }
    
  spp <- unique(taxa_data$CONCEPT)
  
  names(spp) <- spp
  
  taxa_meta <- purrr::map_df(.x = spp, .f = getMeta) # apply getMeta over all species in the group and rbind outputs
  
  save(taxa_meta, file=paste0("F:/BRCIndicatorWorkflow/metaData/", group, ".rdata"))
  
  }
  
}














## This function formats the metadata for use in TrendSummaries. 
## It is horribly bespoke because of the porblems with concept/ latin names, and because I didn't have
## the .rdata sparta inputs to run extractMeta (above) for some taxa. 

sampMeta <- function(group, region, inPath) {

  spp.info <- read.csv(paste0(inPath, "speciesInfo.csv"), stringsAsFactors = FALSE)
  
  metaData <- paste0(inPath, group, ".rdata")
  
  if (file.exists(metaData)) { ## didn't have files with which the run extractMeta for some taxa so needed to use first / last records from a different file (Charlie's)
   
    load(paste0(inPath, group, ".rdata"))
    metaInput <- "new"
    
  } else {
    
    taxa_meta <- read.csv(paste0(inPath, "charlieSpeciesInfo.csv"))
    metaInput <- "charlie"
    
  }

  load(paste0("F:/BRCIndicatorWorkflow/workflow1/", group, "stacked.rdata"))

  spp.names <- unique(stacked_samps$species)

  getMeta <- function(spp) {  # extracts the metadata for the taxa no matter what the source of name format 

    print(spp)
    
    if (tolower(spp) %in% tolower(spp.info$concept)) {
        
      if (tolower(spp) %in% tolower(spp.info$Species)) {
        inc <- spp.info$Reason_not_included[which(tolower(spp.info$Species) == tolower(spp))]
      } else {
        inc <- spp.info$Reason_not_included[which(tolower(spp.info$concept) == tolower(spp))]
      }

      nrecs <- taxa_meta$nRec[which(tolower(taxa_meta$spp) == tolower(spp))]

      first <- taxa_meta$first[which(tolower(taxa_meta$spp) == tolower(spp))]
      
      last <- taxa_meta$last[which(tolower(taxa_meta$spp) == tolower(spp))]
      
      gap <- taxa_meta$gap[which(tolower(taxa_meta$spp) == tolower(spp))]

    } else {

      if (metaInput == "charlie" | tolower(spp) %in% tolower(taxa_meta$spp)) {

        x <- tolower(spp)
        
      } else {
        
        x <- tolower(tolower(spp.info$concept)[tolower(spp.info$Species) == tolower(spp)])

      }

      if (tolower(spp) %in% tolower(spp.info$Species)) {
      inc <- spp.info$Reason_not_included[which(tolower(spp.info$Species) == tolower(spp))]
      } else {
        inc <- NA
      }

      nrecs <- taxa_meta$nRec[which(tolower(taxa_meta$spp) == x)]

      first <- taxa_meta$first[which(tolower(taxa_meta$spp) == x)]
      
      last <- taxa_meta$last[which(tolower(taxa_meta$spp) == x)]
      
      gap <- taxa_meta$gap[which(tolower(taxa_meta$spp) == x)]

    }
    
     if (length(first) < 1 | length(last) <  1 | length(inc) < 1) { # if I did not have the needed metadata, I remove the species by setting "gap" to 11 so it will be filtered out 
       gap <- 11
       first <- 2010
       last <- 2010
       nrecs <- 0
       inc <- NA

  }

    ## for species with fewer than 50 records, they do not have a start or end date. We give them one so that 
    ## stackFilter will work, but it does not matter which dates we provide as the species is removed by
    ## the minObs filter anyway 

    if (is.na(first) | is.na(last)) { 
      first <- 1970
      last <- 2015
    }

    if (!is.na(inc) & inc != "Didn't meet criteria") { ## i.e. species for which scheme adivsed removal are given "0 records" so they are filtered out
      
      nrecs <- 0
      
    }
    
    newRow <- cbind(spp, nrecs, first, last, 1970, 2018, 0, 0, gap)
    
    return(newRow)
    
  }

  names(spp.names) <- spp.names 
  
  out <- data.frame(t(purrr::map_df(.x = spp.names, .f = getMeta)), stringsAsFactors = F) # apply getMeta over all species and rbind outputs
  
  metaNames <- c("Species_r_","n_obs_r_","min_year_data_r_",
                 "max_year_data_r_","min_year_model_r_",
                 "max_year_model_r_",
                 "gap_start_r_","gap_end_r_",
                 "gap_middle_r_")
  
  colnames(out) <- paste0(metaNames, region) # trendSummaries needs these column names
  
  ## convert all columns other than "species" to numeric for use in StackFilter

  out[, 2:ncol(out)] <- sapply(out[, 2:ncol(out)], as.numeric)
  
  return(out)

}














## This function produces a list of species for inclusion or exclusion from the indicator (e.g. pollinators).

sampSubset <- function(subset) {
  
  inPath ="F:/BRCIndicatorWorkflow/metaData/"
  
  spp.list <- read.csv(paste0(inPath, "speciesInfo.csv"), stringsAsFactors = FALSE)

  prio.list <- read.csv(paste0(inPath, "PrioritySpeciesNames_v2.csv"), stringsAsFactors = FALSE)
  
  poll.list <- read.csv(paste0(inPath, "pollinators.csv"), stringsAsFactors = FALSE)

  
  poll2Drop.list <- read.csv(paste0(inPath, "poll2Drop.csv"))
  print(poll2Drop.list)
  if (subset == "priority") {
    
    x <- spp.list$Species[which(tolower(spp.list$Species) %in% tolower(prio.list$NBN_Name) |
                                    tolower(spp.list$Species) %in% tolower(prio.list$DesigName) |
                                    tolower(spp.list$Species) %in% tolower(prio.list$PREFERRED_NBN_NAME) |
                                    tolower(spp.list$concept) %in% tolower(prio.list$MatchName))] ## i.e. any type of latin name or concept name so we don't miss any 
    
    y <- spp.list$concept[which(tolower(spp.list$Species) %in% tolower(prio.list$NBN_Name) |
                                  tolower(spp.list$Species) %in% tolower(prio.list$DesigName) |
                                  tolower(spp.list$Species) %in% tolower(prio.list$PREFERRED_NBN_NAME) |
                                  tolower(spp.list$concept) %in% tolower(prio.list$MatchName))]
    spp <- c(x,y)
      
  } else if (subset == "pollinators") {
    
    x <- spp.list$Species[which(tolower(spp.list$Species) %in% tolower(poll.list$species))]
    
    y <- spp.list$concept[which(tolower(spp.list$Species) %in% tolower(poll.list$species))]
    
    z <- read.csv("F:/BRCIndicatorWorkflow/metaData/hovConc.csv")
    z <- unique(z)
    z <- z[,1]
    
    spp <- c(x,y, as.character(z))
    
  } else if (subset == "poll2Drop") {
    
    spp <- poll2Drop.list$species
    
  }
  
  return(spp)
  
}












## little function to make sure all outputs of sampSubset start have caps for first character only. Will incorporate into sampMeta when I get the chance

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}














## Function to convert data into an array so we can run lambda_indicator 

sampArray <- function(startYear, endYear, niter, inPath) {
  
  "%notin%" <- Negate("%in%")
  
  files <- list.files(inPath, full.names=TRUE, pattern = ".rdata")
  
  combined.df <- NULL
  
  for (i in 1:length(files)) {
    
    load(files[i])
    
    print(files[i])
    
    if ("year_2019" %in% names(stacked_samps)) {
      
      stacked_samps  <- stacked_samps[ , -which(names(stacked_samps) %in% c("year_2019"))]
      
    }
    
    if ("year_2016" %notin% names(stacked_samps)) {
      
      stacked_samps <- tibble::add_column(stacked_samps, year_2016 = rep(NA, nrow(stacked_samps)), .after= "year_2015")
      
    }
    
    if ("year_2017" %notin% names(stacked_samps)) {
      
      stacked_samps <- tibble::add_column(stacked_samps, year_2017 = rep(NA, nrow(stacked_samps)), .after= "year_2016")
      
    }
    
    if ("year_2018" %notin% names(stacked_samps)) {
      
      stacked_samps <- tibble::add_column(stacked_samps, year_2018 = rep(NA, nrow(stacked_samps)), .after= "year_2017")
      
    }
    
    combined.df <- rbind(combined.df, stacked_samps)
    
  }
  
  combined.df <- combined.df[,-ncol(combined.df)] 
  
  combined.df$iteration <- as.numeric(combined.df$iteration)
  
  arr <- simplify2array(by(combined.df, combined.df$iteration, as.matrix))
  
  print(str(arr))
  
  start <- (startYear - 1970) + 1
  end <- (endYear - 1970) + 1  
  
  arr <- arr[,start:end,]

  #tArr <- aperm(a=arr, perm=c(2,1,3), resize=TRUE)
  
  #mode(tArr) <- "numeric"
  
  dimnames(arr)[[1]] <- 1:length(dimnames(arr)[[1]])
  
  dimnames(arr)[[2]] <- 1:length(dimnames(arr)[[2]])
  
  return(arr)
  
}











## function to plot indicators as JNCC require 

plotIndicator <- function(indicator, minYear, maxYear, inPath, label) {
  
  load(paste0(inPath, indicator, "_ind.rdata"))
  
  ind <- ind
  
  ind$summary$indicator

  n <- max(ind$summary$Species_Number)
  
  years <- minYear:maxYear
  
  p1 <- ggplot(data = NULL, aes(x= years, y= ind$summary$indicator)) +
    geom_ribbon(data= NULL, aes(ymax=ind$summary$upper, ymin = ind$summary$lower), fill = "grey80") +
    geom_line() +
    geom_point() +
    theme_linedraw() +
    ylab("Occupancy index") +
    xlab("") +
    ylim(c(0, 180)) +
    ggtitle(label) +
    annotate("text", x=1985, y=30, label= paste(n, "species"))

  p1  
  
  load(paste0(inPath, indicator, "ST.rdata"))
  
  st <- ass$species_assessment$category
  
  st <- data.frame(st, rep(as.factor(1), length(st)))
  
  colnames(st) <- c("val","type")
  
  load(paste0(inPath, indicator, "LT.rdata")) 
  
  lt <- ass$species_assessment$category
  
  lt <- data.frame(lt, rep(as.factor(2), length(lt)))
  
  colnames(lt) <- c("val","type")
  
  dat <- rbind(st,lt)
  
  p2 <- ggplot(dat, aes(x = factor(type), fill = forcats::fct_rev(val))) +
    geom_bar(position="fill") +
    theme_linedraw() +
    ylab("Proportion of species") +
    xlab("") +
    scale_x_discrete(labels=c("Short term","Long term")) +
    guides(fill=guide_legend(title="")) 
  
  p2
  
  #jpeg(paste0(path,"prio_spp_ind.jpg"), width=8, height=3.5,units="in", res=600)
  gridExtra::grid.arrange(p1, p2, ncol=2)
  #dev.off()
}













## extract a taxonomic breakdown of the species contributing to the indicator. Options are nSpec which gives
## the number of species per group, and sppNames which gives the latin names of included species

sampBreakdown <- function(group, indName, inPath, output) {
  
  file <- paste0(inPath, "prio", group, "stacked_Filter.rdata")

  if (file.exists(file)) {

    load(file)
    
  if (output == "nSpecies") {
  
  nSpec <- nrow(stacked_samps) / 1000 # because there are 1000 samples per species
  
  out <- data.frame(group, nSpec)
  
  } else if (output == "sppNames") {
    
    out <- data.frame(unique(stacked_samps$species))
    
  }
  
  return(out)
    
  }

}












## Function to extract summary stats. Produces a list with elements: 1) indicator mean +/- CI; 2) table of trends in long term,
## short-term, and for final year; and 3) raw data (% change) for each species in long, short and final year 

sampSumStats <- function(inPath, indName, startYear) {
  
  load(paste0(inPath, indName, "_ind.rdata"))

  ind$summary$year <- ind$summary$year + (startYear - 1)

  out1 <- ind$summary

  changeLT <- data.frame(table(ind$species_change$category))
  
  rawLT <- ind$species_change[,1]
  
  L <- length(rawLT)

  load(paste0(inPath, indName, "ST.rdata"))
  
  changeST <- data.frame(table(ass$species_assessment$category))[,2]
  
  rawST <- ass$species_assessment[,1]
  
  rawST <- c(rawST, rep(NA, (L - length(rawST))))
  
  load(paste0(inPath, indName, "Final.rdata"))
  
  changeFinal <- data.frame(table(ass$species_assessment$category))[,2]
  
  rawFinal <- ass$species_assessment[,1]
  
  rawFinal <- c(rawFinal, rep(NA, (L - length(rawFinal))))
  
  change <- data.frame(changeLT, changeST, changeFinal)
  
  colnames(change) <- c("Trend", "Long term", "Short term", "Final year")
  
  raw <- data.frame(rawLT, rawST, rawFinal)
  
  colnames(raw) <- c("Long term", "Short term", "Final year")

  out <- list(out1, change, raw)

}




