#
# standards.R
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

#' Default strain-based standards for paw features
#'
#' @return standards used for Jones et al. (2020)
#' 
#' @export
default_standards <- function() jones2020.standards

#' Create new strain standards
#' 
#' Creates a new reference centering and scaling for the paw features of the
#' given strain or strains.
#' 
#' @param paw.features list of extracted paw features
#' @param strain character string or a vector of strain names, each of which
#'        is a character string, matching in length and order the list of
#'        paw features.
#' 
#' @return list of pre-peak and post-peak strain standards, indexed by strain
#' 
#' @export
create_strain_standard <- function(paw.features, strain)
{
	# returns centering and scaling information for given features
	standardize.features <- function(features) {
		zscores  <- scale(features)
		return(list(
			center = attr(zscores, 'scaled:center'),
			scale  = attr(zscores, 'scaled:scale')
		))
	}

	# total number of paw.feature entries
	npaws <- length(paw.features)

	# check that strain information is available
	if (length(strain) == 1) strain <- rep(strain, length(paw.features))
	else if (length(strain) < npaws) {
		stop("missing strain information for one or more paw.features")
	}

	# combine features
	pre.peak.features  <- lapply(paw.features, function(p) p$pre.peak)
	pre.peak.features  <- do.call('rbind', pre.peak.features)

	post.peak.features <- lapply(paw.features, function(p) p$post.peak)
	post.peak.features <- do.call('rbind', post.peak.features)

	# centering and scaling for pre and post features
	pre.peak  <- tapply(1:npaws, strain, function(js) {
		standardize.features(pre.peak.features[js,])
	})
	post.peak <- tapply(1:npaws, strain, function(js) {
		standardize.features(post.peak.features[js,])
	})

	return(list(
		pre.peak  = pre.peak,
		post.peak = post.peak
	))
}

