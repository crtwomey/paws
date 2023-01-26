#
# features.R
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

# Find period of activity in given time series.
get_window <- function(u,
	filter.size,
	filter.order,
	threshold)
{
	# compute smoothed first derivative
	v <- signal::sgolayfilt(u, p=filter.order, n=filter.size, m=1)

	# approximate the beginning and end of the sequence
	rough.start <- head(which(abs(u - head(u,1)) > threshold), 1)
	rough.end   <- tail(which(abs(u - tail(u,1)) > threshold), 1)

	# refine estimate of start and end based on speed
	delta <- diff(sign(v))
	zero.crossing <- which(delta != 0)

	candidates <- zero.crossing[zero.crossing < rough.start]
	start <- ifelse(length(candidates) > 0, max(candidates), rough.start)

	candidates <- zero.crossing[zero.crossing > rough.end]
	end <- ifelse(length(candidates) > 0, min(candidates), rough.end)

	return(list(start=start, end=end))
}

# Project a 2-dimensional time series onto its principal
# axis of variation.
compute_univariate_projection <- function(ux, uy, window) {
	# unit vector in the same direction as vector u
	unit <- function(u) {
		u / sqrt(sum(u^2))
	}
	dot_prod <- function(u,v) {
		sum(u * v)
	}

	U  <- cbind(ux,uy)
	tw <- min(window, floor(nrow(U)/2-1))
	ts <- 1:tw
	v  <- matrix(0, nrow(U)-2*tw, ncol(U))
	m  <- matrix(0, nrow(U)-2*tw, ncol(U))
	for (t in 1:nrow(v)) {
		pc    <- prcomp(U[ts+t-1,], center=TRUE)
		v[t,] <- pc$rotation[,1]
		m[t,] <- pc$center
	}

	w <- matrix(0, nrow(v), ncol(v))
	w[1,] <- v[1,]
	for (t in 2:nrow(w)) {
		s     <- sign(dot_prod(unit(w[t-1,]), unit(v[t,])))
		w[t,] <- unit(v[t,]) * ifelse(s==0,1,s)
	}

	W <- rbind(
		matrix(w[1,], tw, 2, byrow=TRUE),
		w,
		matrix(w[nrow(w),], tw, 2, byrow=TRUE)
	)
	M <- rbind(
		matrix(m[1,], tw, 2, byrow=TRUE),
		m,
		matrix(m[nrow(m),], tw, 2, byrow=TRUE)
	)

	return(list(
		projection = rowSums(W * (U - M)),
		direction  = W
	))
}

# Identify a global peak in the time series (point in time
# at which the overal maximum is achieved).
find_global_peak <- function(u,
	filter.size,
	filter.order)
{
	# compute smoothed first derivative
	v <- signal::sgolayfilt(u, p=filter.order, n=filter.size, m=1)

	# refine estimate of start and end based on speed
	delta <- diff(sign(v))
	zero.crossing <- which(delta != 0)

	global.peak <- head(which(
		u[zero.crossing] == max(u[zero.crossing])
	), 1)
	global.peak.time   <- zero.crossing[global.peak]
	global.peak.height <- u[global.peak.time]

	return(list(
		time   = global.peak.time,
		height = global.peak.height
	))
}

# Find all local peaks in univariate time series of paw motion.
find_local_peaks <- function(u,
	filter.size,
	filter.order,
	threshold)
{
	# compute smoothed first derivative
	v <- signal::sgolayfilt(u, p=filter.order, n=filter.size, m=1)

	# refine estimate of start and end based on speed
	delta <- diff(sign(v))
	zero.crossing <- which(delta != 0)
	if (zero.crossing[1] != 1) zero.crossing <- c(1, zero.crossing)

	d <- c(0,sapply(2:length(zero.crossing), function(i) {
		u[zero.crossing[i]] - u[zero.crossing[i-1]]
	}))
	
	# global peak height
	global.peak <- head(which(
		u[zero.crossing] == max(u[zero.crossing]
	)), 1)
	global.peak.time   <- zero.crossing[global.peak]
	global.peak.height <- u[global.peak.time]

	candidate.peaks <- which(d / global.peak.height > threshold)

	return(zero.crossing[candidate.peaks])
}

# Compute time from occurence of a peak to the end of the time series.
get_time_post_peak <- function(u, peak) {
	time.post.peak <- length(u) - peak$time
	return(time.post.peak)
}

# Estimate the number of shakes based on univariate time series of paw motion.
get_shakes_post_peak <- function(u,
	global.peak.height,
	first.peak.time,
	filter.size, 
	filter.order,
	threshold)
{
	# compute smoothed first derivative
	v <- signal::sgolayfilt(u, p=filter.order, n=filter.size, m=1)

	# identify frames where local minima and maxima occur (zero crossings in
	# the first derivative of the univariate projection time series)
	delta <- diff(sign(v))
	zero.crossing <- which(delta != 0)

	# rescale time series by the maximum height attained. The shake threshold
	# parameter is expressed in terms of this rescaled displacement
	s <- u / global.peak.height

	# include time of the first peak in the set of zero crossings. Used for
	# measuring the displacement between successive critical points
	zero.crossing <- sort(c(zero.crossing, first.peak.time))

	# find zero crossings that exceed the shake threshold, i.e. local minima and
	# maxima where displacement from the previous minimum or maximum (or the
	# first peak) exceeds the shake threshold. The first entry in zero.crossing
	# has no predecessor, so it can never cross the displacement threshold
	past.threshold  <- abs(diff(s[zero.crossing])) > threshold
	past.threshold  <- c(FALSE, past.threshold)
	
	# only retain information about zero crossings after (and excluding) the
	# first peak time
	crossings.past.peak <- zero.crossing > first.peak.time
	threshold.past.peak <- past.threshold[crossings.past.peak]

	# if there are no zero crossings after the first peak, then there are
	# no shakes to report
	if (length(threshold.past.peak) == 0) {
		return(list(
			total.length      = 0,
			total.duration    = 0,
			sequence.times    = c(),
			sequence.lengths  = c(),
			shaking.sequence  = c(),
			zero.crossing     = c(),
			past.threshold    = c(),
			scaled.projection = s
		))
	}

	# identify shaking sequences by consecutive zero crossings that exceed the
	# displacement threshold parameter. A shaking sequence is terminated by a
	# crossing that fails to exceed this threshold. break.points bookend these
	# sequences, starting at 0 (before the first zero crossing after the first
	# peak) and ending with the final break point after the last zero crossing
	break.points <- which(abs(diff(threshold.past.peak))==1)
	break.points <- c(0, break.points)
	if (tail(break.points,1) != length(threshold.past.peak)) {
		break.points <- c(break.points, length(threshold.past.peak))
	}

	# a shaking sequence requires at least 2 zero crossings. A non-shaking
	# sequence has fewer than two crossings or no displacements above threshold.
	# For each sequence (bookended by consecutive break.points), shaking.sequence
	# codes the sequence as either a shaking sequence (TRUE) or not (FALSE)
	shaking.sequence <- sapply(2:length(break.points), function(i) {
		sequence <- (break.points[(i-1)]+1):break.points[i]
		any(threshold.past.peak[sequence]) && (length(sequence) > 1)
	})
	
	# sequence.lengths counts the number of local minima and maxima contained
	# within each sequence (computed for both types: shaking and non-shaking)
	sequence.lengths <- sapply(2:length(break.points), function(i) {
		sequence <- (break.points[(i-1)]+1):break.points[i]
		length(sequence)
	})

	# compute the start and end times of each sequence (both types)
	t0 <- head(which(zero.crossing >= first.peak.time),1)
	sequence.times <- t(sapply(2:length(break.points), function(i) {
		t1 <- t0 + break.points[(i-1)]
		t2 <- t0 + break.points[i]
		return(c(zero.crossing[t1], zero.crossing[t2]))
	}))
	
	# compute the durations of each sequence (both types)
	t0 <- head(which(zero.crossing >= first.peak.time),1)
	sequence.durations <- sapply(2:length(break.points), function(i) {
		t1 <- t0 + (break.points[(i-1)]+1)
		t2 <- t0 + break.points[i]
		zero.crossing[t2] - zero.crossing[t1]
	})

	# compute the total number of shakes and the total duration of shaking
	# periods across all shaking sequences
	tot.length   <- sum(sequence.lengths[shaking.sequence])
	tot.duration <- sum(sequence.durations[shaking.sequence])

	return(list(
		total.length      = tot.length,
		total.duration    = tot.duration,
		sequence.times    = sequence.times,
		sequence.lengths  = sequence.lengths,
		shaking.sequence  = shaking.sequence,
		zero.crossing     = zero.crossing[crossings.past.peak],
		past.threshold    = threshold.past.peak,
		scaled.projection = s
	))
}

# Estimate velocity based on Savitsky-Golay filter.
estimate_velocity <- function(u,
	filter.size,
	filter.order)
{
	# Compute smoothed velocity
	signal::sgolayfilt(u, p=filter.order, n=filter.size, m=1)
}

# Compute max velocity pre and post first peak.
get_max_velocity <- function(gv, first.peak) {
	list(
		pre  = max(abs(gv[1:first.peak$time])),
		post = max(abs(gv[first.peak$time:length(gv)]))
	)
}

# Compute total distance traveled by paw.
get_distance_traveled <- function(x, y) {
	sum(sqrt(diff(x)^2 + diff(y)^2))
}

