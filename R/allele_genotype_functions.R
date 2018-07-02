#' @title A string coverage list
#'
#' @description A list of tibbles, one for every marker, used to contain the sequencing information of STR MPS data.
#' The tibbles should include columns with the following names: "Marker", "BasePairs", "Allele", "Type", "MotifLength", "ForwardFlank", "Region", "ReverseFlank", "Coverage", "AggregateQuality", and "Quality".
setClass("stringCoverageList")

#' String coverage coontrol object
#'
#' @details Control function for the 'stringCoverage' function. Sets default values for the parameters.
#'
#' @param motifLength The motif lengths of each marker.
#' @param Type The chromosome type of each marker (autosomal, X, or Y).
#' @param simpleReturn TRUE/FALSE: Should the returned object be simplified?
#' @param includeLUS TRUE/FALSE: Should the LUS of each region be calculated?
#' @param numberOfThreads The number of cores used for parallelisation.
#' @param includeAverageBaseQuality Should the average base quality of the region be included?
#' @param meanFunction The function used to average the base qualities.
#' @param trace TRUE/FALSE: Show trace?
#' @param uniquelyAssigned TRUE/FALSE: Should regions not uniquely assigned be removed?
#'
#' @return List of parameters used for the 'stringCoverage' function.
stringCoverage.control <- function(motifLength = 4, Type = "AUTOSOMAL", simpleReturn = TRUE, includeLUS = FALSE, numberOfThreads = 4L, meanFunction = mean,
                                   includeAverageBaseQuality = FALSE, trace = FALSE, uniquelyAssigned = TRUE) {
    list(motifLength = motifLength, Type = Type, simpleReturn = simpleReturn, includeLUS = includeLUS, numberOfThreads = numberOfThreads, meanFunction = meanFunction,
         includeAverageBaseQuality = includeAverageBaseQuality, trace = trace, uniquelyAssigned = uniquelyAssigned)
}

.extractedReadsList.stringCoverage <- function(extractedReadsListObject, control = stringCoverage.control()) {
    if (control$uniquelyAssigned) {
        extractedReads <- extractedReadsListObject$identifiedMarkersSequencesUniquelyAssigned
    }
    else {
        warning("Only uniquely assigned sequences extract list should be used")
        extractedReads <- extractedReadsListObject$identifiedMarkers
    }

    if (length(control$motifLength) != length(extractedReads)) {
        if (length(control$motifLength) == 1L) {
            motifLengths <- rep(control$motifLength, length(extractedReads))
        }
        else {
            stop("'motifLength' must have length 1 or the same as 'extractedReads'")
        }

    }
    else {
        motifLengths = control$motifLength
    }

    if (length(control$Type) != length(extractedReads)) {
        if (length(control$Type) == 1L) {
            Types <- rep(control$Type, length(extractedReads))
        }
        else {
            stop("'Type' must have length 1 or the same as 'extractedReads'")
        }

    }
    else {
        Types = control$Type
    }

    alleles <- mclapply(seq_along(extractedReads), function(i) {
        if ((is.null(extractedReads[[i]]$trimmed)) | (dim(extractedReads[[i]]$trimmed)[1] == 0)) {
            return(NULL)
        }

        matchedFlanks <- extractedReads[[i]]

        if (control$trace)
            cat(i, "/", length(extractedReads), ":: Marker:", as.character(matchedFlanks$name), "\n")

        marker = matchedFlanks$name

        stringCoverageQuality <- cbind(matchedFlanks$trimmed, Quality = matchedFlanks$trimmedQuality$Region) %>%
            group_by(ForwardFlank, Region, ReverseFlank) %>%
            summarise(Coverage = n(), AggregateQuality = .aggregateQuality(Quality), Quality = list(as.character(Quality))) %>%
            ungroup() %>%
            mutate(Marker = marker, MotifLength = motifLengths[i], Type = Types[i], BasePairs = nchar(Region),
                   Allele = BasePairs / MotifLength) %>%
            select(Marker, BasePairs, Allele, Type, MotifLength, ForwardFlank, Region, ReverseFlank, Coverage, AggregateQuality, Quality)

        if (control$simpleReturn) {
            stringCoverageQuality <- stringCoverageQuality %>% group_by(Marker, BasePairs, Allele, Type, MotifLength, Region) %>%
                summarise(Coverage = sum(Coverage), AggregateQuality = .aggregateQuality(AggregateQuality), Quality = list(unlist(Quality))) %>% ungroup()
        }

        if (control$includeLUS) {
            validLUS = stringCoverageQuality$Allele >= 1
            stringCoverageQuality$LUS <- sapply(seq_along(stringCoverageQuality$Region), function(ss) ifelse(validLUS[ss], LUS(stringCoverageQuality$Region[ss], motifLength = stringCoverageQuality$MotifLength[ss], returnType = "string"), NA))
        }

        return(stringCoverageQuality)
    }, mc.cores = control$numberOfThreads)

    names(alleles) <- names(extractedReads)
    class(alleles) <- "stringCoverageList"
    return(alleles)
}

#' Get string coverage STR identified objects.
#'
#' \code{stringCoverage} takes an extractedReadsList-object and finds the coverage of every unique string for every marker in the provided list.
#'
#' @param extractedReadsListObject an extractedReadsList-object, created using the \link{identifySTRRegions}-function.
#' @param control an \link{stringCoverage.control}-object.
#' @return Returns a list, with an element for every marker in extractedReadsList-object, each element contains the string coverage of all unique strings of a given marker.
#' @example inst/examples/stringCoverageAggregated.R
setGeneric("stringCoverage", signature = "extractedReadsListObject",
           function(extractedReadsListObject, control = stringCoverage.control())
               standardGeneric("stringCoverage")
)

#' Get string coverage STR identified objects.
#'
#' \code{stringCoverage} takes an extractedReadsList-object and finds the coverage of every unique string for every marker in the provided list.
#'
#' @param extractedReadsListObject an extractedReadsList-object, created using the \link{identifySTRRegions}-function.
#' @param control an \link{stringCoverage.control}-object.
#' @return Returns a list, with an element for every marker in extractedReadsList-object, each element contains the string coverage of all unique strings of a given marker.
#' @example inst/examples/stringCoverageAggregated.R
setMethod("stringCoverage", "extractedReadsList",
           function(extractedReadsListObject, control = stringCoverage.control())
               .extractedReadsList.stringCoverage(extractedReadsListObject, control)
)

#' Get string coverage STR identified objects.
#'
#' \code{stringCoverage} takes an extractedReadsList-object and finds the coverage of every unique string for every marker in the provided list.
#'
#' @param extractedReadsListObject an extractedReadsList-object, created using the \link{identifySTRRegions}-function.
#' @param control an \link{stringCoverage.control}-object.
#' @return Returns a list, with an element for every marker in extractedReadsList-object, each element contains the string coverage of all unique strings of a given marker.
#' @example inst/examples/stringCoverageAggregated.R
setMethod("stringCoverage", "extractedReadsListReverseComplement",
          function(extractedReadsListObject, control = stringCoverage.control())
              .extractedReadsList.stringCoverage(extractedReadsListObject, control)
)

#' Get string coverage STR identified objects.
#'
#' \code{stringCoverage} takes an extractedReadsList-object and finds the coverage of every unique string for every marker in the provided list.
#'
#' @param extractedReadsListObject an extractedReadsList-object, created using the \link{identifySTRRegions}-function.
#' @param control an \link{stringCoverage.control}-object.
#' @return Returns a list, with an element for every marker in extractedReadsList-object, each element contains the string coverage of all unique strings of a given marker.
#' @example inst/examples/stringCoverageAggregated.R
setMethod("stringCoverage", "extractedReadsListCombined",
          function(extractedReadsListObject, control = stringCoverage.control())
              .extractedReadsList.stringCoverage(extractedReadsListObject, control)
)

#' Get string coverage STR identified objects.
#'
#' \code{stringCoverage} takes an extractedReadsList-object and finds the coverage of every unique string for every marker in the provided list.
#'
#' @param extractedReadsListObject an extractedReadsList-object, created using the \link{identifySTRRegions}-function.
#' @param control an \link{stringCoverage.control}-object.
#' @return Returns a list, with an element for every marker in extractedReadsList-object, each element contains the string coverage of all unique strings of a given marker.
#' @example inst/examples/stringCoverageAggregated.R
setMethod("stringCoverage", "extractedReadsListNonCombined",
          function(extractedReadsListObject, control = stringCoverage.control())
              stop("'stringCoverage' not implemented for 'extractedReadsListNReveseComplementList'. Use lapply on the two elements on the list.")
)

#' Genotype list
#'
#' A reduced stringCoverageList restricted to the identified genotypes.
setClass("genotypeIdentifiedList")

#' Noise list
#'
#' Creates a flag to the sequences in a stringCoverageList which cloud be classified as noise.
setClass("noiseIdentifiedList")


.stringCoverageList.NoiseGenotype <- function(stringCoverageListObject, colBelief = "Coverage",
                                              thresholdSignal = 0, thresholdHeterozygosity = 0, thresholdAbsoluteLowerLimit = 1,
                                              trueGenotype = NULL, identified = "genotype") {
    if (length(thresholdSignal) == 1L) {
        if(thresholdSignal < 1 & thresholdSignal > 0) {
            thresholdSignal <- unlist(lapply(stringCoverageListObject, function(s) thresholdSignal*max(s[, colBelief])))
            thresholdSignal <- sapply(seq_along(thresholdSignal), function(s) max(thresholdSignal[s], thresholdAbsoluteLowerLimit))
        }
        else {
            thresholdSignal <- rep(max(thresholdSignal, thresholdAbsoluteLowerLimit), length(stringCoverageListObject))
        }
    }

    if (length(thresholdHeterozygosity) == 1L) {
        thresholdHeterozygosity <- rep(thresholdHeterozygosity, length(stringCoverageListObject))
    }

    if (length(thresholdSignal) != length(stringCoverageListObject)) {
        stop("'stringCoverageListObject' and 'thresholdSignal' must have the same length.")
    }

    if (length(thresholdHeterozygosity) != length(stringCoverageListObject)) {
        stop("'stringCoverageListObject' and 'thresholdHeterozygosity' must have the same length.")
    }

    res <- vector("list", length(stringCoverageListObject))
    for (i in seq_along(stringCoverageListObject)) {
        stringCoverage_i <- stringCoverageListObject[[i]]
        if (is.null(trueGenotype)) {
            belief <- unname(stringCoverage_i[, colBelief] %>% as_vector())
            beliefMax <- max(belief)
            beliefKeepers <- which(belief > thresholdSignal[i] & belief > thresholdHeterozygosity[i]*beliefMax)
        }
        else {
            beliefKeepers <- which(stringCoverage_i$Region %in% trueGenotype[[i]])
        }
        res[[i]] <- stringCoverage_i[beliefKeepers, ] %>% mutate(Indices = beliefKeepers)
    }

    names(res) <- names(stringCoverageListObject)
    class(res) <- if(tolower(identified) == "genotype") "genotypeIdentifiedList" else if(tolower(identified) == "noise") "noiseIdentifiedList"
    return(res)
}

#' Assigns genotype.
#'
#' \code{getGenotype} takes an stringCoverageList-object, assumes the sample is a reference file and assings a genotype, based on a heterozygote threshold, for every marker in the provided list.
#'
#' @param stringCoverageListObject an stringCoverageList-object, created using the \link{stringCoverage}-function.
#' @param colBelief the name of the coloumn used for identification.
#' @param thresholdSignal threshold applied to the signal (generally the coverage) of every string.
#' @param thresholdHeterozygosity threshold used to determine whether a marker is hetero- or homozygous.
#' @param thresholdAbsoluteLowerLimit a lower limit on the coverage for it to be called as an allele.
#'
#' @return Returns a list, with an element for every marker in stringCoverageList-object, each element contains the genotype for a given marker.
#' @example inst/examples/getGenotype.R
setGeneric("getGenotype", signature = "stringCoverageListObject",
           function(stringCoverageListObject, colBelief = "Coverage", thresholdSignal = 0, thresholdHeterozygosity = 0.35, thresholdAbsoluteLowerLimit = 1)
               standardGeneric("getGenotype")
)

#' Assigns genotype.
#'
#' \code{getGenotype} takes an stringCoverageList-object, assumes the sample is a reference file and assings a genotype, based on a heterozygote threshold, for every marker in the provided list.
#'
#' @param stringCoverageListObject an stringCoverageList-object, created using the \link{stringCoverage}-function.
#' @param colBelief the name of the coloumn used for identification.
#' @param thresholdSignal threshold applied to the signal (generally the coverage) of every string.
#' @param thresholdHeterozygosity threshold used to determine whether a marker is hetero- or homozygous.
#' @param thresholdAbsoluteLowerLimit a lower limit on the coverage for it to be called as an allele.
#'
#' @return Returns a list, with an element for every marker in stringCoverageList-object, each element contains the genotype for a given marker.
#' @example inst/examples/getGenotype.R
setMethod("getGenotype", "stringCoverageList",
          function(stringCoverageListObject, colBelief = "Coverage", thresholdSignal = 0, thresholdHeterozygosity = 0.35, thresholdAbsoluteLowerLimit = 1)
              .stringCoverageList.NoiseGenotype(stringCoverageListObject, colBelief, thresholdSignal, thresholdHeterozygosity,
                                                thresholdAbsoluteLowerLimit, NULL, "genotype")
)

#' Idenfities the noise.
#'
#' \code{identifyNoise} takes an stringCoverageList-object and identifies the noise based on a signal threshold for every marker in the provided list.
#'
#' @param stringCoverageListObject an stringCoverageList-object, created using the \link{stringCoverage}-function.
#' @param colBelief the name of the coloumn used for identification.
#' @param thresholdSignal threshold applied to the signal (generally the coverage) of every string.
#'
#' @return Returns a list, with an element for every marker in stringCoverageList-object, each element contains the genotype for a given marker.
#' @example inst/examples/getNoise.R
setGeneric("identifyNoise", signature = "stringCoverageListObject",
           function(stringCoverageListObject, colBelief = "Coverage", thresholdSignal = 0.01)
               standardGeneric("identifyNoise")
)

#' Idenfities the noise.
#'
#' \code{identifyNoise} takes an stringCoverageList-object and identifies the noise based on a signal threshold for every marker in the provided list.
#'
#' @param stringCoverageListObject an stringCoverageList-object, created using the \link{stringCoverage}-function.
#' @param colBelief the name of the coloumn used for identification.
#' @param thresholdSignal threshold applied to the signal (generally the coverage) of every string.
#'
#' @return Returns a list, with an element for every marker in stringCoverageList-object, each element contains the genotype for a given marker.
#' @example inst/examples/getNoise.R
setMethod("identifyNoise", "stringCoverageList",
          function(stringCoverageListObject, colBelief = "Coverage", thresholdSignal = 0.01)
              .stringCoverageList.NoiseGenotype(stringCoverageListObject, colBelief, thresholdSignal, 0, 0, NULL, "noise")
)


.noiseGenotypeIdentified.stringCoverageList.merge <- function(stringCoverageListObject, noiseGenotypeIdentifiedListObject, identified = "genotype") {
    stringCoverageListObjectMerged <- vector("list", length(stringCoverageListObject))
    indValue <- if(tolower(identified) == "genotype") TRUE else if(tolower(identified) == "noise") FALSE
    indCol <- if(tolower(identified) == "genotype") "AlleleCalled" else if(tolower(identified) == "noise") "Noise"

    for(i in seq_along(stringCoverageListObject)) {
        stringCoverageListObjectMerged[[i]] <- stringCoverageListObject[[i]] %>% mutate(tempName = !indValue, FLAGMoreThanTwoAlleles = FALSE)

        if (!is.null(noiseGenotypeIdentifiedListObject[[i]]) && nrow(noiseGenotypeIdentifiedListObject[[i]]) > 0L) {
            stringCoverageListObjectMerged[[i]]$tempName[noiseGenotypeIdentifiedListObject[[i]]$Indices] <- indValue
        }

        if (nrow(noiseGenotypeIdentifiedListObject[[i]]) > 2L) {
            stringCoverageListObjectMerged[[i]]$FLAGMoreThanTwoAlleles <- TRUE
        }

        names(stringCoverageListObjectMerged[[i]]) <- gsub("tempName", indCol, names(stringCoverageListObjectMerged[[i]]))
    }

    names(stringCoverageListObjectMerged) <- names(stringCoverageListObject)
    class(stringCoverageListObjectMerged) <- if(tolower(identified) == "genotype") "stringCoverageGenotypeList" else if(tolower(identified) == "noise") "stringCoverageNoiseList"
    return(stringCoverageListObjectMerged)
}

#' Merge genotypeIdentifiedList and stringCoverageList.
#'
#' \code{mergeGenotypeStringCoverage} merges genotypeIdentifiedList-objects and stringCoverageList-objects.
#'
#' @param stringCoverageListObject a stringCoverageList-object, created using the \link{stringCoverage}-function.
#' @param noiseGenotypeIdentifiedListObject a noiseGenotypeIdentifiedList-object, created using the \link{getGenotype}-function.
#'
#' @return Returns a list, with an element for every marker in extractedReadsList-object, each element contains the string coverage of all unique strings of a given marker.
#' @example inst/examples/mergeLists.R
setGeneric("mergeGenotypeStringCoverage", signature = "noiseGenotypeIdentifiedListObject",
           function(stringCoverageListObject, noiseGenotypeIdentifiedListObject)
               standardGeneric("mergeGenotypeStringCoverage")
)

#' Merge genotypeIdentifiedList and stringCoverageList.
#'
#' \code{mergeGenotypeStringCoverage} merges genotypeIdentifiedList-objects and stringCoverageList-objects.
#'
#' @param stringCoverageListObject a stringCoverageList-object, created using the \link{stringCoverage}-function.
#' @param noiseGenotypeIdentifiedListObject a noiseGenotypeIdentifiedList-object, created using the \link{getGenotype}-function.
#'
#' @return Returns a list, with an element for every marker in extractedReadsList-object, each element contains the string coverage of all unique strings of a given marker.
#' @example inst/examples/mergeLists.R
setMethod("mergeGenotypeStringCoverage", "genotypeIdentifiedList",
          function(stringCoverageListObject, noiseGenotypeIdentifiedListObject)
              .noiseGenotypeIdentified.stringCoverageList.merge(stringCoverageListObject, noiseGenotypeIdentifiedListObject, identified = "genotype")
)

#' Merge noiseIdentifiedList and stringCoverageList.
#'
#' \code{mergeNoiseStringCoverage} merges noiseIdentifiedList-objects and stringCoverageList-objects.
#'
#' @param stringCoverageListObject a stringCoverageList-object, created using the \link{stringCoverage}-function.
#' @param noiseGenotypeIdentifiedListObject a noiseGenotypeIdentifiedList-object, created using the \link{identifyNoise}-function.
#'
#' @return Returns a list, with an element for every marker in extractedReadsList-object, each element contains the string coverage of all unique strings of a given marker.
#' @example inst/examples/mergeLists.R
setGeneric("mergeNoiseStringCoverage", signature = "noiseGenotypeIdentifiedListObject",
           function(stringCoverageListObject, noiseGenotypeIdentifiedListObject)
               standardGeneric("mergeNoiseStringCoverage")
)

#' Merge noiseIdentifiedList and stringCoverageList.
#'
#' \code{mergeNoiseStringCoverage} merges noiseIdentifiedList-objects and stringCoverageList-objects.
#'
#' @param stringCoverageListObject a stringCoverageList-object, created using the \link{stringCoverage}-function.
#' @param noiseGenotypeIdentifiedListObject a noiseGenotypeIdentifiedList-object, created using the \link{identifyNoise}-function.
#'
#' @return Returns a list, with an element for every marker in extractedReadsList-object, each element contains the string coverage of all unique strings of a given marker.
#' @example inst/examples/mergeLists.R
setMethod("mergeNoiseStringCoverage", "noiseIdentifiedList",
          function(stringCoverageListObject, noiseGenotypeIdentifiedListObject)
              .noiseGenotypeIdentified.stringCoverageList.merge(stringCoverageListObject, noiseGenotypeIdentifiedListObject, identified = "noise")
)

#' Combined stringCoverage- and genotypeIdentifiedList
#'
#' Merges a stringCoverageList with a genotypeIdentifiedList.
setClass("stringCoverageGenotypeList")

#' Combined stringCoverage- and noiseIdentifiedList
#'
#' Merges a stringCoverageList with a noiseIdentifiedList
setClass("stringCoverageNoiseList")
