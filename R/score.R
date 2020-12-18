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

# internal function for validating arguments
validate_model_arguments <- function(env=parent.frame()) {
	# validate the following variables
	stopifnot(all(
		exists('paw.features',     envir=env, inherits=FALSE),
		exists('strains',          envir=env, inherits=FALSE),
		exists('strain.standards', envir=env, inherits=FALSE),
		exists('feature.set',      envir=env, inherits=FALSE)
	)); with(env,
{
	# validate the value passed to feature.set
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
	if (class(paw.features) == "paw.features") {
		paw.features <- paw.features[[feature.set]]
	} else if (class(paw.features) == "list") {
		oftype.paw.features <- sapply(paw.features, function(r) {
			class(r) == "paw.features"
	    	})
		if (all(oftype.paw.features)) {
			# extract the correct feature set
			features <- lapply(paw.features, function(f) {
				f[[feature.set]]
			})
			paw.features <- do.call("rbind", features)
		} else stop("not all entries in list are paw.features")
	} else if (class(paw.features) != "data.frame") {
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
	standardized.features <- as.data.frame(standardized.features)
	# note: ordered by strain

	# mapping to recover original order for paw.features
	original.order  <- 1:nrow(paw.features)
	order.by.strain <- unlist(lapply(unique(strains), function(strain) {
		original.order[strains == strain]
	}))

	# return to original order
	standardized.features <- standardized.features[original.order,]
})}

#' Scoring pain from paw features
#' 
#' Returns a pain score based on Jones et al. 2020 or the given pain.model.
#' Pain scores are standardized so that increasingly positive (negative) values
#' correspond to increasingly likely to be painful (unlikely to be painful)
#' experiences. Pain scores based on pain models derived from trinary pain
#' stimuli will additionally be scaled such that the transition between low
#' and high pain is occurs at 1. Scores are not comparable across pain models.
#' 
#' @param paw.features the paw features returned by
#'        \code{\link{extract_features}}
#' @param strains a vector containing strain information for each mouse
#' @param strain.standards z-scores (centering and scaling) information
#'        by strain
#' @param feature.set use either pre-peak or post-peak features
#' @param pain.model a pain model returned by \code{\link{create_pain_model}}
#'        or NULL, in which case the model from Jones et al. 2020 is used.
#' 
#' @return one or more pain scores as a vector
#' 
#' @export
pain_score <- function(paw.features,
	strains          = NULL,
	strain.standards = jones2020.standards,
	feature.set      = c("post.peak","pre.peak"),
	pain.model       = NULL)
{
	# validate arguments
	validate_model_arguments()

	# validate pain.model if one is provided
	pain.type <- NA
	if (!is.null(pain.model)) {
		pain.types <- c("binary", "trinary")
		pain.type  <- attr(pain.model, "pain.type")
		if (!(pain.type %in% pain.types)) {
			stop("unknown pain model type")
		}

		# ensure pain model is for the requested feature.set
		if (attr(pain.model, "feature.set") != feature.set) {
			stop(paste("pain.model calibrated for ",
				attr(pain.model, "feature.set"), "features,",
				"but pain scores requested for",
				feature.set, "features"
			))
		}
	} else {
		# Jones et al. 2020 data supports trinary pain classifications
		pain.type <- "trinary"
	}

	# generate pain scores
	X <- as.matrix(standardized.features)
	pain.score <- NA
	if (is.null(pain.model)) {
		# untransformed pain scores using ordinal logistic regression
		# coefficients from Jones et al. (2020) stored in R/sysdata.rda
		Y <- X %*% coef(olog[[feature.set]])

		# identify pain transition points
		nonpain_pain_threshold       <- olog[[feature.set]]$zeta[2]
		low.pain_high.pain_threshold <- olog[[feature.set]]$zeta[3]

		# offset and scale pain scores by transition points
		delta <- low.pain_high.pain_threshold - nonpain_pain_threshold
		pain.score <- (Y - nonpain_pain_threshold) / delta
	} else {
		# use provided pain.model for scoring
		cc <- coef(pain.model)
		pain.score <- switch(pain.type,
			binary  = cbind(1,X) %*% cc,
			trinary = (cbind(1,X) %*% cc[-1]) / (cc[2]-cc[1])
		)
	}

	# return annotated pain scores
	pain.score <- as.vector(pain.score)
	attr(pain.score, "pain.type")   <- pain.type
	attr(pain.score, "feature.set") <- feature.set
	return(pain.score)
}

#' Create a new pain model
#'
#' Create a new pain model based on a set of paw features, corresponding pain
#' stimuli, and corresponding strain identity.
#' 
#' @param paw.features the paw features returned by
#'        \code{\link{extract_features}}
#' @param strains a vector containing strain information for each mouse
#' @param pain.stimulus a vector containing stimulus information for each mouse
#' @param strain.standards z-scores (centering and scaling) information
#'        by strain
#' @param feature.set use either pre-peak or post-peak features
#' 
#' @return pain model
#' 
#' @export
create_pain_model <- function(paw.features,
	strains          = NULL,
	pain.stimulus    = NULL,
	strain.standards = jones2020.standards,
	feature.set      = c("post.peak","pre.peak"))
{
	# validate arguments
	validate_model_arguments()

	# validate pain stimulus information
	if (   is.null(pain.stimulus)
	    || length(pain.stimulus) != nrow(paw.features))
	{
		stop("need pain stimulus for each paw.features entry")
	}

	# pain types (note: order of trinary ensures comparable axis to binary)
	binary    <- c("non.pain", "pain")
	trinary   <- rev(c("non.pain", "low.pain", "high.pain"))
	pain.type <- c("binary","trinary")[c(
		all(pain.stimulus %in% binary),
		all(pain.stimulus %in% trinary)
	)]
	if (length(pain.type) != 1) {
		stop("unknown pain type(s) in pain.stimulus")
	}

	# check that all pain levels are present for the given pain type
	if (pain.type == "binary" && length(unique(pain.stimulus)) != 2) {
		stop("need both non.pain and pain pain.stimulus data")
	} else if (   pain.type == "trinary"
	           && length(unique(pain.stimulus)) != 3)
	{
		stop(paste("need non.pain, low.pain, and high.pain", 
		           "pain.stimulus data"))
	}

	# ensure pain.stimulus is an ordinal factor
	pain.levels   <- switch(pain.type, binary=binary, trinary=trinary)
	pain.stimulus <- factor(pain.stimulus,
		levels = pain.levels,
		order  = TRUE
	)

	# combine standardized features with pain.stimulus information
	standardized.features$stimulus <- pain.stimulus

	# create model based on pain type
	model <- NULL
	if (pain.type == "binary") {
		# run logistic regression
		model <- glm(stimulus ~ .,
			family = binomial(logit),
			data   = standardized.features,
			method = "brglmFit"
		)
	} else {
		# run ordinal logistic regression
		model <- bracl(stimulus ~ .,
			parallel = TRUE,
			data     = standardized.features
		)
	}

	# return annotated model
	attr(model, 'feature.set') <- feature.set
	attr(model, 'pain.type')   <- pain.type
	return(model)
}

#' Identify pain class based on pain score
#' 
#' Convenience function to convert pain scores to binary or trinary pain
#' classifications, e.g. pain / non-pain or pain / low-pain / high-pain. Note
#' that trinary classifications are only valid for pain models fit with trinary
#' classes (e.g. Jones et al. 2020). Otherwise the scale of the score is
#' arbitrary, and the boundary between low and high pain is not scaled to be 1.
#'
#' @param score pain scores generated by \code{\link{pain_score}}
#' @param type binary or trinary (ternary) pain classifications
#' 
#' @return one or more pain classes
#' 
#' @export
pain_class <- function(score, type=c("binary","trinary")) {
	type <- match.arg(type)
	if (attr(score, "pain.type") == "binary" && type == "trinary") {
		stop(paste("score computed from pain model based on binary",
			   "pain classes only."))
	}
	if (type == "binary") {
		ifelse(score < 0, "non.pain", "pain")
	} else if (type == "trinary") {
		ifelse( score < 0, "non.pain", ifelse(
			score < 1, "low.pain", "high.pain"
		))
	} else stop("type not recognized")
}

