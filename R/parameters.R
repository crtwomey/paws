#
# parameters.R
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

#' Default parameters for extracting paw features
#' 
#' @return parameters used for Jones et al. (2020)
#' 
#' @export
default_parameters <- function() jones2020.parameters

#' Convenience function for changing parameters
#' 
#' @param ... name and value of parameters to set (e.g. fps=1500)
#' @param based.on an existing complete set of parameters to modify
#'        (see \code{\link{default_parameters}}).
#' 
#' @return the modified list of parameters
#' 
#' @export
set_parameters <- function(..., based.on = default_parameters()) {
	params <- list(...)
	based.on[names(params)] <- params
	return(based.on)
}

# Round to the nearest odd integer (useful for calculating filter window size)
round_odd <- function(x) {
	v <- floor(x)
	ifelse(v %% 2 == 0, v+1, v)
}

# Convert from time to frames
convert_to_frames <- function(parameters) {
	# identify parameters that need to be converted
	# from time (seconds) to frames
	params  <- names(parameters)
	convert <- c(
		which(params == "projection.window"),
		grep("filter.size", params)
	)

	# convert to frames
	parameters[convert] <- with(parameters, sapply(
		parameters[convert],
		function(p) round_odd(p * fps / downsample))
	)
	return(parameters)
}

