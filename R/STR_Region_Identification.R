#' Extract STR region information
#'
#' Identifies the marker of the read using flanking regions and trims the read to include what is between the flanking regions.
setClass("extractedReadsList")

#' Extract STR region information of the reverse complement DNA strand.
#'
#' Identifies the marker of the read using reverse complement flanking regions and trims the read to include what is between the flanking regions.
setClass("extractedReadsListReverseComplement")

#' Combined extract STR region information.
#'
#' Identifies the marker of the read for both the provided and reverse complement flanking regions. The resulting lists are then combined into a single list.
setClass("extractedReadsListCombined")

#' Combined extract STR region information.
#'
#' Identifies the marker of the read for both the provided and reverse complement flanking regions.
setClass("extractedReadsListNonCombined", representation(identifiedReads = "extractedReadsList", identifiedReverseComplementReads = "extractedReadsListReverseComplement"))

.getCols <- function(colNames) {
    structure(lapply(c("marker", "forwardflank", "reverseflank"), function(n) grep(n, tolower(colNames))), .Names = c("markerCol", "forwardCol", "reverseCol"))
}

## .mclapplyvmatchPattern and .vmatchMultiplePatternsAlternate are used for diagnostic purposes only
.mclapplyvmatchPattern <- function(flanks, seqs, max.mismatch, numberOfThreads, limitSeqRange = NULL) {
    if (is.null(limitSeqRange)) {
        resList <- mclapply(seq_along(flanks), function(f) vmatchPattern(flanks[[f]], seqs, with.indels = FALSE, max.mismatch = max.mismatch),
                            mc.cores = numberOfThreads)
    }
    else {
        resList <- mclapply(seq_along(flanks), function(f) {
            if (!is.null(limitSeqRange[[f]])) {
                res <- vmatchPattern(flanks[[f]], seqs[limitSeqRange[[f]]], with.indels = FALSE, max.mismatch = max.mismatch)
                return(res)
            }
            else {
                return(NULL)
            }
        }, mc.cores = numberOfThreads)
    }
    return(resList)
}

.identifyFlankingRegions <- function(seqs, flankingRegions, colList, numberOfMutation = 1, numberOfThreads = 4, removeEmptyMarkers = TRUE,
                                     remainingSequences = NULL) {
    if (!is.null(remainingSequences))
        if (!is.list(remainingSequences))
            remainingSequences = rep(list(remainingSequences), times = length(flankingRegions$Marker))

    forwardFlank <- DNAStringSet(unname(unlist(flankingRegions[, colList$forwardCol])))
    reverseFlank <- DNAStringSet(unname(unlist(flankingRegions[, colList$reverseCol])))

    # Forward match
    identifiedForward <- .mclapplyvmatchPattern(forwardFlank, seqs, numberOfMutation, numberOfThreads, limitSeqRange = remainingSequences)
    forwardMatch <- lapply(1:nrow(flankingRegions), function(i) which(elementNROWS(identifiedForward[[i]]) >= 1L))

    # Reverse match
    identifiedReverse <- .mclapplyvmatchPattern(reverseFlank, seqs, numberOfMutation, numberOfThreads, limitSeqRange = forwardMatch)
    reverseMatch <- lapply(1:nrow(flankingRegions), function(i) which(elementNROWS(identifiedReverse[[i]]) >= 1L))

    matchedSeq <- lapply(seq_along(reverseMatch), function(i) forwardMatch[[i]][reverseMatch[[i]]])

    if(removeEmptyMarkers) {
        nonEmptyEntries <- which(sapply(matchedSeq, length) != 0)
        rList <- list(markers = unname(unlist(flankingRegions[nonEmptyEntries, colList$markerCol])),
                      identifiedForward = identifiedForward[nonEmptyEntries], identifiedReverse = identifiedReverse[nonEmptyEntries],
                      matchedSeq = matchedSeq[nonEmptyEntries], reverseMatchedSeq = reverseMatch[nonEmptyEntries],
                      remainingSequences = seq_along(seqs)[-unique(unlist(matchedSeq, recursive = TRUE))])
    }
    else {
        rList <- list(markers = unname(unlist((flankingRegions[, colList$markerCol]))), identifiedForward = identifiedForward,
                      identifiedReverse = identifiedReverse, matchedSeq = matchedSeq, reverseMatchedSeq = reverseMatch,
                      remainingSequences = seq_along(seqs)[-unique(unlist(matchedSeq, recursive = TRUE))])
    }

    return(rList)
}


#####
.extractAndTrimMarkerIdentifiedReads.mclapply <- function(seqs, qual, flankingRegions, colList, numberOfMutation, removeEmptyMarkers,
                                                          flankSizes, numberOfThreads = 4, reversed = FALSE, reverseComplementRun = FALSE,
                                                          remainingSequences = NULL) {
    identifiedFlanksObj <- .identifyFlankingRegions(seqs = seqs, flankingRegions = flankingRegions, colList = colList,
                                                        numberOfMutation = numberOfMutation, numberOfThreads = numberOfThreads,
                                                        removeEmptyMarkers = removeEmptyMarkers, remainingSequences = remainingSequences)

    whichMarkers = which(unlist(flankingRegions[, colList$markerCol]) %in% identifiedFlanksObj$markers)
    identifiedMarkers <- structure(mclapply(seq_along(identifiedFlanksObj$matchedSeq), function(i) {
        marker <- identifiedFlanksObj$markers[i]
        identifiedForward <- identifiedFlanksObj$identifiedForward[[i]]
        identifiedReverse <- identifiedFlanksObj$identifiedReverse[[i]]
        matchedSeq <- identifiedFlanksObj$matchedSeq[[i]]
        reverseMatchedSeq <- identifiedFlanksObj$reverseMatchedSeq[[i]]

        # F: Take the first in the IRanges
        endForward <- unlist(lapply(endIndex(identifiedForward)[matchedSeq], function(v) v[1L]))
        startForward <- endForward - (flankSizes[whichMarkers[i], 1] - 1)

        # R: Take the last in the IRanges
        startReverse <- unlist(lapply(startIndex(identifiedReverse)[reverseMatchedSeq], function(v) v[length(v)]))
        endReverse <- startReverse + (flankSizes[whichMarkers[i], 2] - 1)

        # Removes reads where forward is NOT observed before reverse
        keepSeq <- which(startReverse > (endForward + 1) & endReverse <= nchar(seqs[matchedSeq]) & startForward > 0)
        matchedSeq <- matchedSeq[keepSeq]

        startForward <- startForward[keepSeq]
        endForward <- endForward[keepSeq]

        startReverse <- startReverse[keepSeq]
        endReverse <- endReverse[keepSeq]

        if (length(keepSeq) == 0L) {
            rList <- NULL
        }
        else {
            trimmedIncludingFlanks <- subseq(seqs[matchedSeq], start = startForward, end = endReverse)
            trimmedQualityIncludingFlanks <- subseq(quality(qual[matchedSeq]), start = startForward, end = endReverse)

            if (!(reversed && reverseComplementRun)) {
                forwardFlankingRegion <- subseq(seqs[matchedSeq], start = startForward, end = endForward + 1)
                reverseFlankingRegion <- subseq(seqs[matchedSeq], start = startReverse - 1, end = endReverse)

                trimmed <- subseq(seqs[matchedSeq], start = endForward + 1, end = startReverse - 1)

                forwardFlankingRegionQuality <- subseq(quality(qual[matchedSeq]), start = startForward, end = endForward + 1)
                reverseFlankingRegionQuality <- subseq(quality(qual[matchedSeq]), start = startReverse - 1, end = endReverse)

                trimmedQuality <- subseq(quality(qual[matchedSeq]), start = endForward + 1, end = startReverse - 1)
            }
            else {
                forwardFlankingRegion <- reverseComplement(subseq(seqs[matchedSeq], start = startReverse - 1, end = endReverse))
                reverseFlankingRegion <- reverseComplement(subseq(seqs[matchedSeq], start = startForward, end = endForward + 1))

                trimmed <- reverseComplement(subseq(seqs[matchedSeq], start = endForward + 1, end = startReverse - 1))

                forwardFlankingRegionQuality <- reverse(subseq(quality(qual[matchedSeq]), start = startReverse - 1, end = endReverse))
                reverseFlankingRegionQuality <- reverse(subseq(quality(qual[matchedSeq]), start = startForward, end = endForward + 1))

                trimmedQuality <- reverse(subseq(quality(qual[matchedSeq]), start = endForward + 1, end = startReverse - 1))
            }

            trimmedTibble <- tibble(ForwardFlank = as.character(forwardFlankingRegion), Region = as.character(trimmed),
                                    ReverseFlank = as.character(reverseFlankingRegion))
            trimmedQualityTibble <- tibble(ForwardFlank = as.character(forwardFlankingRegionQuality), Region = as.character(trimmedQuality),
                                           ReverseFlank = as.character(reverseFlankingRegionQuality))
            rList <- list(name = marker, matchedSeq = matchedSeq,
                          startForward = startForward, endForward = endForward,
                          startReverse = startReverse, endReverse = endReverse,
                          trimmedIncludingFlanks = trimmedIncludingFlanks, trimmedQualityIncludingFlanks = trimmedQualityIncludingFlanks,
                          trimmed = trimmedTibble, trimmedQuality = trimmedQualityTibble)
        }

        return(rList)
    }, mc.cores = numberOfThreads), .Names = as.character(identifiedFlanksObj$markers)) #

    # Remove sequences occuring at several loci
    identifiedSeqTable <- sort(table(unlist(lapply(identifiedMarkers, function(r) r$matchedSeq))), decreasing = TRUE)
    identifiedSeqNotUnique <- as.integer(names(identifiedSeqTable[identifiedSeqTable > 1L]))
    identifiedMarkersUniqueSequences <- lapply(identifiedMarkers, function(r) {
        keepSeq <- which(!(r$matchedSeq %in% identifiedSeqNotUnique))
        return(list(name = r$name, matchedSeq = r$matchedSeq[keepSeq], startForward = r$startForward[keepSeq], endForward = r$endForward[keepSeq],
                    startReverse = r$startReverse[keepSeq], endReverse = r$endReverse[keepSeq],
                    trimmedIncludingFlanks = r$trimmedIncludingFlanks[keepSeq], trimmedQualityIncludingFlanks = r$trimmedQualityIncludingFlanks[keepSeq],
                    trimmed = r$trimmed[keepSeq, ], trimmedQuality = r$trimmedQuality[keepSeq, ]))
    })

    extractedReadsList <- list(n_reads = sum(sapply(identifiedMarkersUniqueSequences, function(vv) length(vv[["matchedSeq"]]))),
                               reverseComplement = reverseComplementRun, identifiedMarkers = identifiedMarkers,
                               identifiedMarkersSequencesUniquelyAssigned = identifiedMarkersUniqueSequences,
                               remainingSequences = identifiedFlanksObj$remainingSequences)
    class(extractedReadsList) <- "extractedReadsList"
    return(extractedReadsList)
}

.extractAndTrimMarkerIdentifiedReadsReverseComplement <- function(seqs, qual, flankingRegions, colList = NULL,
                                                                  numberOfMutation, removeEmptyMarkers = TRUE,
                                                                  numberOfThreads = numberOfThreads, reversed = TRUE,
                                                                  matchPatternMethod = "mclapply",
                                                                  remainingSequences = NULL) {

    flankingRegionsReverseComplement <- structure(data.frame(flankingRegions[, colList$markerCol],
                                                             unname(as.character(reverseComplement(DNAStringSet(as_vector(flankingRegions[, colList$reverseCol]))))),
                                                             unname(as.character(reverseComplement(DNAStringSet(as_vector(flankingRegions[, colList$forwardCol]))))),
                                                             stringsAsFactors = FALSE),
                                                  .Names = c("Marker", "ForwardRC", "ReverseRC"))

    extractAndTrimMarkerIdentifiedReads = switch(matchPatternMethod,
                                                 "mclapply" = .extractAndTrimMarkerIdentifiedReads.mclapply)

    extractedTrimmedRC <- extractAndTrimMarkerIdentifiedReads(seqs, qual, flankingRegionsReverseComplement,
                                                              colList = list(markerCol = 1, forwardCol = 2, reverseCol = 3),
                                                              numberOfMutation = numberOfMutation, removeEmptyMarkers = removeEmptyMarkers,
                                                              flankSizes = apply(flankingRegionsReverseComplement[, -1], 2, nchar),
                                                              numberOfThreads = numberOfThreads,
                                                              reversed = reversed, reverseComplementRun = TRUE,
                                                              remainingSequences = remainingSequences)

    class(extractedTrimmedRC) <- "extractedReadsListReverseComplement"
    return(extractedTrimmedRC)
}

.combineMarkerIdentifiedReadsLists <- function(extractedReadsList1, extractedReadsList2) {
    uniqueMarkers <- unique(c(names(extractedReadsList1$identifiedMarkersSequencesUniquelyAssigned), names(extractedReadsList2$identifiedMarkersSequencesUniquelyAssigned)))
    combinedList <- list(n_reads = extractedReadsList1$n_reads + extractedReadsList2$n_reads)

    combinedList$identifiedMarkers <- structure(lapply(uniqueMarkers, function(j) {
        if (is.null(extractedReadsList2$identifiedMarkers[[j]])) return(extractedReadsList1$identifiedMarkers[[j]])
        else if (is.null(extractedReadsList1$identifiedMarkers[[j]])) return(extractedReadsList2$identifiedMarkers[[j]])
        else return(.appendExtractLists(extractedReadsList1$identifiedMarkers[[j]], extractedReadsList2$identifiedMarkers[[j]]))
    }), .Names = c(uniqueMarkers))

    combinedList$identifiedMarkersSequencesUniquelyAssigned <- structure(lapply(uniqueMarkers, function(j) {
        if (is.null(extractedReadsList2$identifiedMarkersSequencesUniquelyAssigned[[j]])) return(extractedReadsList1$identifiedMarkersSequencesUniquelyAssigned[[j]])
        else if (is.null(extractedReadsList1$identifiedMarkersSequencesUniquelyAssigned[[j]])) return(extractedReadsList2$identifiedMarkersSequencesUniquelyAssigned[[j]])
        else return(.appendExtractLists(extractedReadsList1$identifiedMarkersSequencesUniquelyAssigned[[j]], extractedReadsList2$identifiedMarkersSequencesUniquelyAssigned[[j]]))
    }), .Names = c(uniqueMarkers))

    class(combinedList) <- "extractedReadsListCombined"
    return(combinedList)
}

#' Control function for identifySTRRegions
#'
#' A list containing default parameters passed to the \link{identifySTRRegions} function.
#'
#' @param colList The position of the forward, reverse, and motifLength columns in the flanking region tibble/data.frame. If 'NULL' a function searches for the words 'forward', 'reverse', and 'motif' ot identify the columns.
#' @param numberOfThreads The number of threads used by mclapply (stuck at '2' on windows).
#' @param reversed TRUE/FALSE: In a revrse complementary run, should the strings/quality be reversed (recommended)?
#' @param includeReverseComplement TRUE/FALSE: Should the function also search for the reverse complement DNA strand (recommended)?
#' @param combineLists TRUE/FALSE: If 'includeReverseComplement' is TRUE, should the sets be combined?
#' @param removeEmptyMarkers TRUE/FALSE: Should markers returning no identified regions be removed?
#' @param matchPatternMethod Which method should be used to identify the flanking regions (only 'mclapply' implemented at the moment)?
#'
#' @return A control list setting default behaviour.
identifySTRRegions.control <- function(colList = NULL, numberOfThreads = 4L, reversed = TRUE,
                                       includeReverseComplement = TRUE, combineLists = TRUE, removeEmptyMarkers = TRUE,
                                       matchPatternMethod = "mclapply") {
    controlList <- list(colList = NULL, numberOfThreads = numberOfThreads, removeEmptyMarkers = removeEmptyMarkers,
                        reversed = reversed, includeReverseComplement = includeReverseComplement,
                        combineLists = combineLists, matchPatternMethod = matchPatternMethod)
    return(controlList)
}


.ShortReadQ.identifySTRRegions <- function(reads, flankingRegions, numberOfMutation, control = identifySTRRegions.control()) {
    seqs <- ShortRead::sread(reads)
    qual <- Biostrings::quality(reads)

    colID <- if(is.null(control$colList)) .getCols(names(flankingRegions)) else control$colList
    flankSizes <- apply(flankingRegions[, c(colID$forwardCol, colID$reverseCol)], 2, nchar)

    extractAndTrimMarkerIdentifiedReads = switch(control$matchPatternMethod,
                                                 "mclapply" = .extractAndTrimMarkerIdentifiedReads.mclapply)

    extractedSTRs <- extractAndTrimMarkerIdentifiedReads(seqs = seqs, qual = qual, flankingRegions = flankingRegions,
                                                         colList = colID, numberOfMutation = numberOfMutation,
                                                         removeEmptyMarkers = control$removeEmptyMarkers, flankSizes = flankSizes,
                                                         numberOfThreads = control$numberOfThreads, reversed = FALSE,
                                                         reverseComplementRun = FALSE, remainingSequences = NULL)

    if (control$includeReverseComplement & (length(extractedSTRs$remainingSequences) > 0)) {
        extractedSTRs_RC <- .extractAndTrimMarkerIdentifiedReadsReverseComplement(seqs, qual, flankingRegions, colList = colID,
                                                                                  numberOfMutation = numberOfMutation,
                                                                                  removeEmptyMarkers = control$removeEmptyMarkers,
                                                                                  numberOfThreads = control$numberOfThreads,
                                                                                  reversed = control$reversed,
                                                                                  matchPatternMethod = control$matchPatternMethod,
                                                                                  remainingSequences = NULL) #extractedSTRs$remainingSequences)

        if (control$combineLists) {
            combinedRegions <- .combineMarkerIdentifiedReadsLists(extractedSTRs[which(names(extractedSTRs) != "remainingSequences")],
                                                                  extractedSTRs_RC[which(names(extractedSTRs) != "remainingSequences")])
            return(combinedRegions)
        }
        else {
            extractedSTRs_Dual <- list(identifiedReads = extractedSTRs, identifiedReverseComplementReads = extractedSTRs_RC)
            class(extractedSTRs_Dual) <- "extractedReadsListNonCombined"
            return(extractedSTRs_Dual)
        }
    }
    else {
        return(extractedSTRs)
    }
}

.character.identifySTRRegions <- function(reads, flankingRegions, numberOfMutation, control = identifySTRRegions.control()) {
    if(!file.exists(reads))
        stop("File does not exist.")

    .ShortReadQ.identifySTRRegions(reads = readFastq(reads), flankingRegions = flankingRegions, numberOfMutation = numberOfMutation, control = control)
}

#' Identify the STR regions of a fastq-file or ShortReadQ-object.
#'
#' \code{identifySTRRegions} takes a fastq-file location or a ShortReadQ-object and identifies the STR regions
#' based on a directly adjacent flanking regions.
#' The function allows for mutation in the flanking regions through the numberOfMutation argument.
#'
#' @param reads either a fastq-file location or a ShortReadQ-object
#' @param flankingRegions containing marker ID/name, the directly adjacent forward and reverse flanking regions, used for identification.
#' @param numberOfMutation the maximum number of mutations (base-calling errors) allowed during flanking region identification.
#' @param control an \link{identifySTRRegions.control}-object.
#'
#' @return The returned object is a list of lists. If the reverse complement strings are not included or if the \code{control$combineLists == TRUE},
#' a list, contains lists of untrimmed and trimmed strings for each row in \code{flankingRegions}. If \code{control$combineLists == FALSE}, the function returns a list of two such lists,
#' one for forward strings and one for the reverse complement strings.
#' @example inst/examples/identify.R
setGeneric("identifySTRRegions", signature = "reads",
           function(reads, flankingRegions, numberOfMutation, control)
               standardGeneric("identifySTRRegions")
)

#' Identify the STR regions of a fastq-file or ShortReadQ-object.
#'
#' \code{identifySTRRegions} takes a fastq-file location or a ShortReadQ-object and identifies the STR regions
#' based on a directly adjacent flanking regions.
#' The function allows for mutation in the flanking regions through the numberOfMutation argument.
#'
#' @param reads a ShortReadQ-object
#' @param flankingRegions containing marker ID/name, the directly adjacent forward and reverse flanking regions, used for identification.
#' @param numberOfMutation the maximum number of mutations (base-calling errors) allowed during flanking region identification.
#' @param control an \link{identifySTRRegions.control}-object.
#'
#' @return The returned object is a list of lists. If the reverse complement strings are not included or if the \code{control$combineLists == TRUE},
#' a list, contains lists of untrimmed and trimmed strings for each row in \code{flankingRegions}. If \code{control$combineLists == FALSE}, the function returns a list of two such lists,
#' one for forward strings and one for the reverse complement strings.
#' @example inst/examples/identify.R
setMethod("identifySTRRegions", "ShortReadQ",
          function(reads, flankingRegions, numberOfMutation = 1, control = identifySTRRegions.control())
              .ShortReadQ.identifySTRRegions(reads = reads, flankingRegions = flankingRegions, numberOfMutation = numberOfMutation, control = control)
)

#' Identify the STR regions of a fastq-file or ShortReadQ-object.
#'
#' \code{identifySTRRegions} takes a fastq-file location or a ShortReadQ-object and identifies the STR regions
#' based on a directly adjacent flanking regions.
#' The function allows for mutation in the flanking regions through the numberOfMutation argument.
#'
#' @param reads path to fastq-file.
#' @param flankingRegions containing marker ID/name, the directly adjacent forward and reverse flanking regions, used for identification.
#' @param numberOfMutation the maximum number of mutations (base-calling errors) allowed during flanking region identification.
#' @param control an \link{identifySTRRegions.control}-object.
#'
#' @return The returned object is a list of lists. If the reverse complement strings are not included or if the \code{control$combineLists == TRUE},
#' a list, contains lists of untrimmed and trimmed strings for each row in \code{flankingRegions}. If \code{control$combineLists == FALSE}, the function returns a list of two such lists,
#' one for forward strings and one for the reverse complement strings.
#' @example inst/examples/identify.R
setMethod("identifySTRRegions", "character",
          function(reads, flankingRegions, numberOfMutation = 1, control = identifySTRRegions.control())
              .character.identifySTRRegions(reads = reads, flankingRegions = flankingRegions, numberOfMutation = numberOfMutation, control = control)
)
