#' A neighbour list
#'
#' A list of the identified neightbours of the called alleles in a stringCoverageGenotypeList
setClass("neighbourList")

.findNeighbourStrings <- function(strings, alleles_i, motifLength, searchDirection, gapOpeningPenalty, gapExtensionPenalty) {
    motifDifference <- motifLength*abs(searchDirection)

    trueStutters <- vector("list", length(alleles_i))
    for(j in seq_along(alleles_i)) {
        alleles_j <- alleles_i[j]

        neighbourRepeatLength <- strings$Allele[alleles_j] + searchDirection
        neighbours_j <- which(abs(strings$Allele - neighbourRepeatLength) < 1e-10)
        neighbours_j <- neighbours_j[which((strings$ForwardFlank[alleles_j] == strings$ForwardFlank[neighbours_j]) & (strings$ReverseFlank[alleles_j] == strings$ReverseFlank[neighbours_j]))]

        if (length(neighbours_j) == 0) {
            next
        }

        subMatrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -nchar(strings$Region[alleles_j]), baseOnly = FALSE)
        stutterAligned <- pairwiseAlignment(DNAStringSet(as.character(strings$Region[neighbours_j])), strings$Region[alleles_j], substitutionMatrix = subMatrix,
                                            gapOpening = -gapOpeningPenalty, gapExtension = -gapExtensionPenalty)

        trueStutters[[j]] <- which(stutterAligned@score == (nchar(strings$Region[alleles_j]) - motifDifference - (gapOpeningPenalty + motifDifference*gapExtensionPenalty)))
    }

    df <- vector("list", length(alleles_i))
    for(j in seq_along(alleles_i)) {
        df_j <- vector("list", length = length(trueStutters[[j]]))
        if (length(trueStutters[[j]]) > 0) {
            alleles_j <- alleles_i[j]
            entireParentRepeatStructure <- BLMM(as.character(strings$Region[alleles_j]), motifLength, returnType = "fullList")
            lusOfMotifs <- entireParentRepeatStructure %>% group_by(Motif) %>% filter(Repeats == max(Repeats)) %>% ungroup()
            lus <- which.max(lusOfMotifs$Repeats)

            alleleRepeatLength <- strings$Allele[alleles_j]
            neighbourRepeatLength <- strings$Allele[alleles_j] + searchDirection
            neighbours_j <- which(abs(strings$Allele - neighbourRepeatLength) < 1e-10)
            neighbours_j <- neighbours_j[which((strings$ForwardFlank[alleles_j] == strings$ForwardFlank[neighbours_j]) & (strings$ReverseFlank[alleles_j] == strings$ReverseFlank[neighbours_j]))]

            if (length(neighbours_j) == 0) {
                next
            }

            subMatrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -nchar(strings$Region[alleles_j]), baseOnly = FALSE)
            stutterAligned <- pairwiseAlignment(DNAStringSet(as.character(strings$Region[neighbours_j])), strings$Region[alleles_j], substitutionMatrix = subMatrix,
                                                gapOpening = -gapOpeningPenalty, gapExtension = -gapExtensionPenalty)

            calledNeighbours <- which(strings$AlleleCalled[neighbours_j[trueStutters[[j]]]])
            for (k in seq_along(trueStutters[[j]])) {
                if ((j > 1) & (k %in% calledNeighbours)) {
                    next
                }
                if (searchDirection == -1) {
                    missingRepeatUnitStartPosition <- which(unlist(strsplit(as.character(aligned(stutterAligned)[trueStutters[[j]][k]]), "")) == "-")[1]
                    entireParentRepeatStructure_k <- entireParentRepeatStructure[which((missingRepeatUnitStartPosition >= entireParentRepeatStructure$Start) & (missingRepeatUnitStartPosition < entireParentRepeatStructure$End)),]
                    endingMotif <- entireParentRepeatStructure_k$Motif[which(entireParentRepeatStructure_k$End == (missingRepeatUnitStartPosition + motifDifference))]
                    missingRepeatUnit <- entireParentRepeatStructure_k$Motif
                    occurenceInParent <- entireParentRepeatStructure_k$Repeats
                }
                else {
                    missingRepeatUnit = NA
                    occurenceInParent = NA
                }

                AlleleDifference = -1
                if (j == 1 && length(alleles_i) > 1L) {
                    AlleleDifference <- strings$Allele[alleles_i[j + 1]] - strings$Allele[alleles_i[j]]
                }
                else if (length(alleles_i) > 1L) {
                    AlleleDifference <- strings$Allele[alleles_i[j]] - strings$Allele[alleles_i[j - 1]]
                }

                FLAGStutterIdentifiedMoreThanOnce <- FALSE
                if (AlleleDifference == 0) {
                    stutterID <- if (j == 1) (trueStutters[[j]][k] %in% trueStutters[[j + 1]]) else (trueStutters[[j]][k] %in% trueStutters[[j - 1]])
                    if (stutterID) FLAGStutterIdentifiedMoreThanOnce <- TRUE
                }

                FLAGMoreThanOneBlock <- ifelse(length(missingRepeatUnit) > 1, TRUE, FALSE)
                FLAGBlocksWithDifferentLengths <- ifelse(length(unique(occurenceInParent)) > 1, TRUE, FALSE)

                FLAGMoreThanTwoAlleles <- ifelse(length(alleles_i) > 2, TRUE, FALSE)
                FLAGAlleleDifferenceOne <- ifelse(length(calledNeighbours) > 0, TRUE, FALSE)

                neighbourRatio <- strings$Coverage[neighbours_j[trueStutters[[j]][k]]] / strings$Coverage[alleles_j]
                neighbourProportion <- strings$Coverage[neighbours_j[trueStutters[[j]][k]]] / (strings$Coverage[neighbours_j[trueStutters[[j]][k]]] + strings$Coverage[alleles_j])

                motifCycles <- sapply(entireParentRepeatStructure_k$Motif, function(m) .cyclicRotation(endingMotif, m))
                setOccurenceInParent <- max(occurenceInParent[motifCycles])

                df_j[[k]] <- data.frame(Genotype = paste(strings$Allele[alleles_i], collapse = ",", sep = ""),
                                    ParentAllele = alleleRepeatLength,
                                    ParentString = strings$Region[alleles_j],
                                    ParentLUS = paste("[", lusOfMotifs$Motif[lus], "]", lusOfMotifs$Repeats[lus], sep=""),
                                    ParentLUSLength = lusOfMotifs$Repeats[lus],
                                    ParentCoverage = strings$Coverage[alleles_j],
                                    NeighbourAllele = neighbourRepeatLength,
                                    NeighbourString = strings$Region[neighbours_j[trueStutters[[j]][k]]],
                                    NeighbourCoverage = strings$Coverage[neighbours_j[trueStutters[[j]][k]]],
                                    Block = paste("[", missingRepeatUnit, "]", occurenceInParent, sep = "", collapse = "/"),
                                    MissingMotif = paste(missingRepeatUnit, sep = "", collapse = "/"),
                                    BlockLengthMissingMotif = setOccurenceInParent,
                                    NeighbourRatio = neighbourRatio,
                                    NeighbourProportion = neighbourProportion,
                                    FLAGStutterIdentifiedMoreThanOnce = FLAGStutterIdentifiedMoreThanOnce,
                                    FLAGMoreThanTwoAlleles = FLAGMoreThanTwoAlleles,
                                    FLAGAlleleDifferenceOne = FLAGAlleleDifferenceOne,
                                    FLAGMoreThanOneBlock = FLAGMoreThanOneBlock,
                                    FLAGBlocksWithDifferentLengths = FLAGBlocksWithDifferentLengths)
            }

            df[[j]] <- do.call(rbind, df_j)
        }
    }

    df_res <- do.call(rbind, df)
    if (is.null(df_res)) {
        df_res <- data.frame(Genotype = NA, ParentAllele = NA, ParentString = NA,
                         ParentLUS = NA, ParentLUSLength = NA, ParentCoverage = NA, NeighbourAllele = NA,
                         NeighbourString = NA, NeighbourCoverage = NA, Block = NA, MissingMotif = NA, BlockLengthMissingMotif = NA,
                         NeighbourRatio = NA, NeighbourProportion = NA,
                         FLAGStutterIdentifiedMoreThanOnce = FALSE, FLAGMoreThanTwoAlleles = FALSE, FLAGAlleleDifferenceOne = FALSE,
                         FLAGMoreThanOneBlock = FALSE,
                         FLAGBlocksWithDifferentLengths = FALSE)
    }

    df_res <- df_res %>% as_tibble()
    return(df_res)
}

# trace = T; searchDirection = -1; gapOpeningPenalty = 6; gapExtensionPenalty = 1; i = 1; j = 1; k = 1
.findNeighbours <- function(stringCoverageGenotypeListObject, searchDirection, gapOpeningPenalty = 6, gapExtensionPenalty = 1, trace = FALSE) {
    if (length(searchDirection) != length(stringCoverageGenotypeListObject)) {
        if (length(searchDirection) == 1) {
            searchDirection <- rep(searchDirection, length(stringCoverageGenotypeListObject))
        }
    }

    res <- vector("list", length(stringCoverageGenotypeListObject))
    for (i in seq_along(stringCoverageGenotypeListObject)) {
        if (trace) {
            cat("Marker:", names(stringCoverageGenotypeListObject[i]), "::", i, "/", length(stringCoverageGenotypeListObject), "\n")
        }

        strings <- stringCoverageGenotypeListObject[[i]]
        if (is.null(suppressWarnings(strings$ForwardFlank))) {
            strings <- strings %>% mutate(ForwardFlank = "A")
        }
        if (is.null(suppressWarnings(strings$ReverseFlank))) {
            strings <- strings %>% mutate(ReverseFlank = "B")
        }

        searchDirection_i = searchDirection[i]
        alleles_i <- which(strings$AlleleCalled)

        motifLength <- round(unique(strings$MotifLength))
        if (length(alleles_i) > 0) {
            df <- .findNeighbourStrings(strings, alleles_i, motifLength, searchDirection_i, gapOpeningPenalty, gapExtensionPenalty)
            res[[i]] <- bind_cols(tibble(Marker = rep(names(stringCoverageGenotypeListObject[i]), dim(df)[1])), df)
        }
    }

    class(res) <- "neighbourList"
    return(res)
}


#' @title Find stutters
#'
#' @description Given identified alleles it search for '-1' stutters of the alleles.
#'
#' @param stringCoverageGenotypeListObject A \link{stringCoverageGenotypeList-class} object.
#' @param trace Should a trace be shown?
#'
#' @return A 'neighbourList' with the stutter strings for the identified allele regions.
#' @example inst/examples/stutter.R
setGeneric("findStutter", signature = "stringCoverageGenotypeListObject",
           function(stringCoverageGenotypeListObject, trace = FALSE)
               standardGeneric("findStutter")
)

#' @title Find stutters
#'
#' @description Given identified alleles it search for '-1' stutters of the alleles.
#'
#' @param stringCoverageGenotypeListObject A \link{stringCoverageGenotypeList-class} object.
#' @param trace Should a trace be shown?
#'
#' @return A 'neighbourList' with the stutter strings for the identified allele regions.
#' @example inst/examples/stutter.R
setMethod("findStutter", "stringCoverageGenotypeList",
          function(stringCoverageGenotypeListObject, trace = FALSE)
              .findNeighbours(stringCoverageGenotypeListObject, searchDirection = -1, gapOpeningPenalty = 6, gapExtensionPenalty = 1, trace)
)


#' @title Find neighbours
#'
#' @description Generic function for finding neighbouring strings, given identified alleles.
#'
#' @param stringCoverageGenotypeListObject A \link{stringCoverageGenotypeList-class} object.
#' @param searchDirection The direction to search for neighbouring strings. Default is -1, indicating a search for '-1' stutters.
#' @param trace Should a trace be shown?
#'
#' @return A 'neighbourList' with the neighbouring strings, in the specified direction, for the identified allele regions.
setGeneric("findNeighbours", signature = "stringCoverageGenotypeListObject",
           function(stringCoverageGenotypeListObject, searchDirection, trace = FALSE)
               standardGeneric("findNeighbours")
)

#' @title Find neighbours
#'
#' @description Generic function for finding neighbouring strings, given identified alleles.
#'
#' @param stringCoverageGenotypeListObject A \link{stringCoverageGenotypeList-class} object.
#' @param searchDirection The direction to search for neighbouring strings. Default is -1, indicating a search for '-1' stutters.
#' @param trace Should a trace be shown?
#'
#' @return A 'neighbourList' with the neighbouring strings, in the specified direction, for the identified allele regions.
setMethod("findNeighbours", "stringCoverageGenotypeList",
          function(stringCoverageGenotypeListObject, searchDirection = -1, trace = FALSE)
              .findNeighbours(stringCoverageGenotypeListObject, searchDirection, gapOpeningPenalty = 6, gapExtensionPenalty = 1, trace)
)
