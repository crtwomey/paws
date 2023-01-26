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
#' @param diagnostics set to TRUE will record intermediate values
#'        used when computing paw features. This information can
#'        be helpful for debugging parameter choices. The default,
#'        FALSE, is to not record these values.
#' 
#' @return pre-peak and post-peak paw features (plus diagnostics,
#'         if enabled)
#' 
#' @export
extract_features <- function(x, y=NULL,
	parameters  = default_parameters(),
	diagnostics = FALSE)
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
	if (diagnostics) {
		features <- c(features, list(
				diagnostics        = list(
					x                = x,
					y                = y,
					x.window         = x.window,
					y.window         = y.window,
					window           = window,
					x.clipped        = x.clipped,
					y.clipped        = y.clipped,
					u.projection     = u,
					global.peak      = global.peak,
					local.peaks      = local.peaks,
					first.peak       = first.peak,
					x.velocity       = gvx,
					y.velocity       = gvy,
					shakes.post.peak = shakes.post.peak
				)
		))
	}
	attr(features, 'class') <- 'paw.features'
	return(features)

	}) # with parameters
}

plot_diagnostics <- function(features, clipped = FALSE,
	panel = c("all",
		"x displacement", "y displacement", "vx velocity", "vy velocity",
		"univariate projection", "scaled univariate projection",
		"kinematics", "displacement", "velocity", "projections"
	)
) {
	if (!exists("diagnostics", features)) {
		stop(paste0(
			"No diagnostics available. Diagnostic values need to be generated ",
			"by passing diagnostics = TRUE to the call to extract_features."
		))
	}
	panel <- match.arg(panel)

	with(features$diagnostics, {
		s <- window$start
		w <- if (clipped) {
			c(window$start, window$end)
		} else {
			c(1, length(x))
		}
		decorate <- function() {
			abline(h = 0, col = "gray60")
			if (!clipped) {
				abline(v = c(window$start, window$end), lty = 2, col = 4)
			}
			abline(v = s + first.peak$time, lty = 1, col = 2)
			axis(3, at = s + first.peak$time, labels = expression(t^"*"),
				col.axis = "red", tick = FALSE, line = -0.8
			)
		}

		mfrows <- switch(panel,
			all          = c(3, 2),
			kinematics   = c(2, 2),
			displacement = c(1, 2),
			velocity     = c(1, 2),
			projections  = c(1, 2),
			c(1, 1)
		)
		op <- par(mar = c(5, 5, 4, 1), mfrow = mfrows)

		if (panel %in% c("all", "kinematics", "displacement", "x displacement")) {
			plot(x, type = "l", las = 1, xaxs = "i", xlim = w,
				xlab = "time (frames)",
				main = "horizontal displacement"
			)
			decorate()
		}

		if (panel %in% c("all", "kinematics", "velocity", "vx velocity")) {
			plot(s + seq_along(x.velocity), x.velocity,
				type = "l", las = 1, xaxs = "i", xlim = w,
				xlab = "time (frames)",
				ylab = expression(v[x]),
				main = "estimated horizontal velocity"
			)
			decorate()
		}

		if (panel %in% c("all", "kinematics", "displacement", "y displacement")) {
			plot(y, type = "l", las = 1, xaxs = "i", xlim = w,
				xlab = "time (frames)",
				main = "vertical displacement"
			)
			points(s + global.peak$time, y[s + global.peak$time],
				col = 4, pch = 16, cex = 1.25
			)
			points(s + local.peaks, y[s + local.peaks], cex = 1.25)
			decorate()
		}

		if (panel %in% c("all", "kinematics", "velocity", "vy velocity")) {
			plot(s + seq_along(y.velocity), y.velocity,
				type = "l", las = 1, xaxs = "i", xlim = w,
				xlab = "time (frames)",
				ylab = expression(v[y]),
				main = "estimated vertical velocity"
			)
			decorate()
		}

		if (panel %in% c("all", "projections", "univariate projection")) {
			plot(s + seq_along(u.projection), u.projection,
				type = "l", las = 1, xaxs = "i", xlim = w,
				xlab = "time (frames)",
				ylab = "u",
				main = "univariate projection"
			)
			decorate()
		}

		if (panel %in% c("all", "projections", "scaled univariate projection")) {
			shakes <- shakes.post.peak
			ymin   <- min(shakes$scaled.projection)
			ymax   <- max(shakes$scaled.projection)
			bump   <- (ymax - ymin) * 0.05
			plot(NA, type = "l", las = 1, xaxs = "i", xlim = w, ylim = c(ymin, ymax),
				xlab = "time (frames)",
				ylab = "s",
				main = "scaled univariate projection"
			)
			decorate()
			for (i in seq_len(nrow(shakes$sequence.times))) {
				t0 <- s + shakes$sequence.times[i, 1]
				t1 <- s + shakes$sequence.times[i, 2]
				ss <- 1 + shakes$shaking.sequence[i]
				cl <- c("black", "red")
				rect(
					xleft = t0, ybottom = ymin - bump, xright = t1, ytop = ymax + bump,
					col    = adjustcolor(cl, alpha.f = 0.1)[ss],
					border = "black",
					lty    = "dotted"
				)
				axis(3, at = t0 + (t1 - t0) / 2, tick = FALSE, line = -0.8,
					labels   = shakes$sequence.lengths[i],
					col.axis = cl[ss]
				)
			}
			lines(s + seq_along(shakes$scaled.projection), shakes$scaled.projection)
			points(s + shakes$zero.crossing,
				shakes$scaled.projection[shakes$zero.crossing],
				col = 1 + shakes$past.threshold
			)
		}
		par(op)
	})
}

#' @method plot paw.features
#' @export
plot.paw.features <- function(features, ...) {
	if (exists("diagnostics", features)) {
		plot_diagnostics(features, ...)
	} else {
		with(features$time.series, {
			decorate <- function() {
				abline(h = 0, col = "gray60")
				abline(v = tstar, lty = 1, col = 2, lwd = 2)
			}
			op <- par(mar = c(5, 5, 1, 1), mfrow = c(3, 1))
			plot(x, type = "l", las = 1, xaxs = "i",
				xlab = "time (frames)"
			)
			decorate()
			plot(y, type = "l", las = 1, xaxs = "i",
				xlab = "time (frames)"
			)
			decorate()
			plot(u, type = "l", las = 1, xaxs = "i",
				xlab = "time (frames)",
			ylab = "univariate projection"
			)
			decorate()
			par(op)
		})
	}
}


