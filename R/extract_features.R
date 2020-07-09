#
# extract_features.R
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

#' Jones et al. (2020) paw trajectory data
#' 
#' Paw trajectory time series, strain, and stimulus information
#'  for the cohort 1 data used in Jones et al.
#' 
#' @source Jones et al. (2020) A machine-vision approach for automated pain
#'         measurement at millisecond timescales. bioRxiv.
#'         \url{https://doi.org/10.1101/2020.02.18.955070}
#' @format A list of paw trajectories, each containing:
#' \describe{
#'   \item{\code{id}}{A unique id for each mouse}
#'   \item{\code{strain}}{The corresponding mouse strain}
#'   \item{\code{stimulus}}{The stimulus used}
#'   \item{\code{time.series}}{The paw trajectory when stimulus was applied}
#' }
"jones2020.tracks"

#' Extract features for paw time series
#' 
#' @param x time series of horizontal paw movement. Alternatively,
#'        a two column matrix of x and y time series, respectively.
#' @param y time series of vertical paw movement, or NULL if x is
#'        a two column matrix.
#' @param parameters contains information about frames per second,
#'        filtering, windowing, and thresholds, for paw features
#'        (see \code{\link{default_parameters}}, or use
#'        \code{\link{set_parameters}} to modify the defaults).
#' 
#' @return pre-peak and post-peak paw features
#' 
#' @export
extract_features <- function(x, y=NULL,
	parameters = default_parameters())
{
	# check that time series data is provided as either
	# 1. x : vector, y : vector, length(x) == length(y)
	# 2. x : matrix, with ncol == 2
	if (is.vector(x)) {
	       if (!is.vector(y) || length(x) != length(y)) {
			stop("x and y must be equal length vectors")
	       }
	} else if (is.matrix(x)) {
		if (ncol(x) == 2) {
			# standardize input
			y <- x[,2]
			x <- x[,1]
		} else stop("x is a matrix but does not have exactly 2 columns")
	}  else stop("invalid time series x y data")

	# put parameters in units convenient for processing
	parameters <- convert_to_frames(parameters)
	with(parameters, {

	# clip sequence based on activity
	x.window <- get_window(
		x, window.filter.size, window.filter.order, window.threshold
	)
	y.window <- get_window(
		y, window.filter.size, window.filter.order, window.threshold
	)
	window <- list(
		start = min(x.window$start, y.window$start),
		end   = max(x.window$end,   y.window$end)
	)

	# check that we found a valid window on paw activity
	if (is.na(window$start) || is.na(window$end)) {
		stop("failed to find paw activity window")
	}

	x.clipped <- with(window, x[start:end])
	y.clipped <- with(window, y[start:end])

	# standardize height baseline
	taus      <- 1:length(y.clipped)
	ntaus     <- length(taus)
	baseline  <- y.clipped[1] * (ntaus-(taus-1))/(ntaus-1) +
	             y.clipped[ntaus] * (taus-1)/(ntaus-1)
	y.clipped <- (y.clipped - baseline)

	# compute univariate projection that preserves as
	# much variability in paw movement as possible.
	up <- compute_univariate_projection(
		x.clipped, y.clipped, projection.window
	)
	u <- up$projection
	
	# determine the global peak in paw height
	global.peak <- find_global_peak(
		y.clipped, global.peak.filter.size, global.peak.filter.order
	)

	# determine the first local peak to occur in paw height
	local.peaks <- find_local_peaks(
		y.clipped, local.peak.filter.size, local.peak.filter.order,
		local.peak.threshold
	)

	# get first peak (or global peak if no first peak detected)
	first.peak <- if (length(local.peaks) == 0) {
		global.peak
	} else list(
		time   = local.peaks[1],
		height = y.clipped[local.peaks[1]]
	)

	# compute maximum local peak height *after* the first peak,
	# if a second peak exists
	max.height.post.peak <- max(y.clipped[-(1:local.peaks[1])])

	# compute features of paw movement
	time.post.peak   <- get_time_post_peak(y.clipped, first.peak)
	shakes.post.peak <- get_shakes_post_peak(
		u, global.peak$height, first.peak$time,
		shake.filter.size, shake.filter.order, shake.threshold
	)

	# estimated velocities
	gvx <- estimate_velocity(
		x.clipped, velocity.filter.size, velocity.filter.order
	)
	gvy <- estimate_velocity(
		y.clipped, velocity.filter.size, velocity.filter.order
	)

	# max velocities (pre/post first peak)
	max.x.velocity <- get_max_velocity(gvx, first.peak)
	max.y.velocity <- get_max_velocity(gvy, first.peak)

	# compute guarding duration
	guarding.duration <- time.post.peak - shakes.post.peak$total.duration

	# compute distance paw travels
	distance.traveled.prepeak <- get_distance_traveled(
		x.clipped[1:first.peak$time],
		y.clipped[1:first.peak$time]
	)
	distance.traveled.postpeak <- get_distance_traveled(
		x.clipped[first.peak$time:length(x.clipped)],
		y.clipped[first.peak$time:length(y.clipped)]
	)

	# features computed on pre-peak period only
	pre.peak.measurements <- data.frame(
		max.height        = first.peak$height,
		max.x.velocity    = max.x.velocity$pre * fps,
		max.y.velocity    = max.y.velocity$pre * fps,
		distance.traveled = distance.traveled.prepeak
	)

	# features computed on post-peak period only
	post.peak.measurements <- data.frame(
		max.height        = max.height.post.peak,
		max.x.velocity    = max.x.velocity$post * fps,
		max.y.velocity    = max.y.velocity$post * fps,
		distance.traveled = distance.traveled.postpeak,
		number.of.shakes  = shakes.post.peak$total.length,
		shaking.duration  = shakes.post.peak$total.duration / fps,
		guarding.duration = guarding.duration / fps
	)

	features <- list(
		pre.peak    = pre.peak.measurements,
		post.peak   = post.peak.measurements,
		time.series = list(
			x     = x.clipped,
			y     = y.clipped,
			u     = u,
			tstar = first.peak$time
		),
		parameters = parameters
	)
	attr(features, 'class') <- 'paw.features'
	return(features)

	}) # with parameters
}

