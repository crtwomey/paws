#
# score.R
# Copyright (c) 2019, 2020 Colin Twomey
# All rights reserved.
#
# This file is part of PAWS.
#
# PAWS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PAWS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PAWS.  If not, see <https://www.gnu.org/licenses/>.
#

#' Scoring pain from paw features
#' 
#' Returns a pain score based on Jones et al. 2020
#' 
#' @param paw.features the paw features returned by
#'        \code{\link{extract_features}}
#' @param strains a vector of containing strain information for each mouse
#' @param strain.standards z-scores (centering and scaling) information
#'        by strain
#' @param feature.set use either pre-peak or post-peak features
#' 
#' @return one or more pain scores as a vector
#' 
#' @export
pain_score <- function(paw.features,
	strains          = NULL,
	strain.standards = jones2020.standards,
	feature.set      = c("post.peak","pre.peak"))
{
	# sanity check the value passed to feature.set
	feature.set <- match.arg(feature.set)
	if (   length(feature.set) != 1
	    || !(feature.set %in% c("post.peak", "pre.peak")))
       	{
		stop(paste("feature.set must be either post.peak",
		           "(default) or pre.peak"))
	}

	# choose pre- pr post-peak standards if available
	standards <- strain.standards
	if (any(c("post.peak", "pre.peak") %in% names(strain.standards))) {
		if (!(feature.set %in% names(strain.standards))) {
			stop(paste(feature.set,
				"not provided in strain.standards", sep=" "
			))
		}
		standards <- strain.standards[[feature.set]]
	}

	# extract the correct data.frame if given a paw.features object
	# or a list of paw.features objects.  otherwise check that we
	# have at least been given a data.frame.
	if (class(paw.features) == 'paw.features') {
		paw.features <- paw.features[[feature.set]]
	} else if (is.list(paw.features)) {
		oftype.paw.features <- sapply(paw.features, function(r) {
			class(r) == 'paw.features'
	    	})
		if (all(oftype.paw.features)) {
			# extract the correct feature set
			features <- lapply(paw.features, function(f) {
				f[[feature.set]]
			})
			paw.features <- do.call("rbind", features)
		} else stop("not all entries in list are paw.features")
	} else if (class(paw.features) != 'data.frame') {
		stop("invalid paw.features data")
	}

	# check that strain information is provided
	if (is.null(strains)) {
		if (is.null(paw.features$strain)) {
			stop("missing strain information")
		} else strains <- paw.features$strain
	} else if (!is.null(paw.features$strain)) {
		# strain information provided twice
		if (any(strains != paw.features$strain)) {
			stop("inconsistent strain information provided")
		}
	}

	# check that we have strain information for every paw.features row
	if (length(strains) != nrow(paw.features)) {
		stop("need strain information for each paw.features row")
	}

	# check that we have centering and scaling information for
	# every strain requested
	if (!all(strains %in% names(standards))) {
		stop(paste("missing centering and scaling information for",
		           "some strains"))
	}
	
	# check that we have the correct set of features in the provided
	# paw.features data.frame and strain.standards
	use.features <- names(jones2020.standards[[feature.set]][[1]][[1]])
	if (!all(use.features %in% names(paw.features))) {
		stop("missing features in provided paw.features data")
	}
	if (!all(use.features %in% names(standards[[1]][[1]]))) {
		stop("missing features in provided strain.standards")
	}

	# compute paw features standardized by strain
	standardized.features <- lapply(unique(strains), function(strain) {
		scale(paw.features[strains == strain, use.features],
			center = standards[[strain]]$center,
			scale  = standards[[strain]]$scale
		)
	})
	standardized.features <- do.call("rbind", standardized.features)

	# untransformed pain scores using ordinal logistic regression
	# coefficients from Jones et al. (2020) stored in R/sysdata.rda
	X <- as.matrix(standardized.features)
	Y <- X %*% coef(olog[[feature.set]])

	# identify nonpain -> pain and low-pain -> high pain transition points
	nonpain_pain_threshold       <- olog[[feature.set]]$zeta[2]
	low.pain_high.pain_threshold <- olog[[feature.set]]$zeta[3]

	# offset and scale pain scores by pain threshold transition points
	delta      <- low.pain_high.pain_threshold - nonpain_pain_threshold
	pain.score <- (Y - nonpain_pain_threshold) / delta

	# return results in original order
	original.order  <- 1:nrow(paw.features)
	order.by.strain <- unlist(lapply(unique(strains), function(strain) {
		original.order[strains == strain]
	}))
	pain.score[order.by.strain] <- pain.score

	return(as.vector(pain.score))
}

#' Identify pain class based on pain score
#' 
#' @param score pain scores generated by \code{\link{pain_score}}
#' @param type binary or trinary (ternary) pain classifications
#' 
#' @return one or more pain classes
#' 
#' @export
pain_class <- function(score, type=c("binary","trinary")) {
	type <- match.arg(type)
	if (type == "binary") {
		ifelse(score < 0, "non.pain", "pain")
	} else if (type == "trinary") {
		ifelse( score < 0, "non.pain", ifelse(
			score < 1, "low.pain", "high.pain"
		))
	} else stop("type not recognized")
}

