
## paws: Pain Assessment at Withdrawal Speeds

Automated pain scoring from paw withdrawl tracking data based on
[Jones et al. (2020)](https://doi.org/10.1101/2020.02.18.955070)
A machine-vision approach for automated pain measurement at millisecond
timescales. This R package takes paw trajectory data in response to a
stimulus and provides an automated scoring of pain.

### Installing

Using the devtools package (the first line can be skipped if devtools is
already installed), in R run:

```R
install.packages("devtools")
devtools::install_github("crtwomey/paws")
```

### Compiling

To compile from the commandline, cd to the parent directory of the paws
repository and run[^1]

```console
$ R CMD build paws
$ R CMD install paws_1.0.0.tar.gz
```

[^1]: The `$` symbol just refers to your commandline's prompt. It may
  look different on different systems. When typing these commands into
  your commandline, don't include the `$` symbol.

### Usage

The package documentation explains the three main functions:
`extract_features`, `pain_score`, and `pain_class`. A simple example
is shown here using paw trajectory data from Jones et al.,
which is included with the package as `jones2020.tracks`.

```R
# compute features for Jones et al. (2020) paw trajectories
paw.features <- lapply(jones2020.tracks, function(track) {
	extract_features(track$time.series)
})

# get strain information for each track
strains <- sapply(jones2020.tracks, function(track) track$strain)

# compute pain scores
scores <- pain_score(paw.features, strains=strains)
```

Plotting the kernel density estimates of score distributions by stimulus
(excluding strain SJL; see the paper) reproduces Fig. 4I of Jones et al.

![Figure 4I](/man/figures/figure4I.svg)


### Diagnostics

It's a good idea to check the output of the feature extraction process. Your
experimental conditions may differ from those of Jones et al., and require
different parameters. At a minimum, check that the time series of paw
displacements returned by the `extract_features` method captures the window of
paw activity (from the time when the paw first leaves the ground to when it last
returns) and identifies the time of the first peak in paw height. This can be
done as follows:

```R
# using the first track from the Jones et al. (2020) data as an example
paw <- jones2020.tracks[[1]]$time.series

# extract features for this time series
features <- extract_features(paw)

# plotting the features will show the x (horizontal) and y (vertical)
# displacement of the paw over time, alongside the computed univariate
# projection used for estimating shaking and guarding sequences
plot(features)
```

![Paw trajectory timeseries](/man/figures/timeseries.svg)

The time of the first peak in vertical paw displacement is indicated by the red
vertical line shown on each panel. The position of this peak can be checked in
the second panel, which shows vertical displacement of the paw. The first panel
shows horizontal displacement relative to the starting position of the paw.  The
time series of vertical paw displacements is shown using a linearly adjusted
baseline that ensures the paw's vertical position starts and ends at zero. The
time (in frames) is relative to the estimated window of paw activity.

Note that `extract_features` assumes that displacement from the ground increases
as `y` increases. If this is not true for your data because your tracking
software uses a reference frame with the opposite orientation, you will need to
flip your `y` axis data before passing it to `extract_features`. A simple way to
do this is to pass `max(y) - y` instead.

More detailed information is available if requested when calling the
`extract_features` function by including the `diagnostics = TRUE` parameter, as
in the following example.

```R
# extract features with diagnostic information
features <- extract_features(paw, diagnostics = TRUE)

# plot diagnostics for kinematic data (x and y displacements and velocity)
plot(features, panel = "kinematics")
```

![Kinematics diagnostics panel](/man/figures/kinematics.svg)

Calling the `plot` function with diagnostic enriched extracted features will
provide more detailed diagnostic plots. Including the `panel = "kinematics"`
parameter in the above example will display panels for horizontal and vertical
(x and y) displacements and estimated velocities. Here, frames are shown
numbered identically to the original time series passed to `extract_features`,
and the estimated window of paw activity is shown as dashed blue vertical bars.
The time series can be clipped to just the window of activity by passing
`clipped = TRUE` to the `plot` function (shown in the next example). The
vertical displacement panel additionally shows estimated local peaks in paw
height (the first of which determines $t^\star$), and the maximum height
attained is indicated by a filled blue point.

If no panel or subset of panels is specified, then all diagnostic plots are
shown, as in the next example.

```R
# plot all diagnostics using the extracted paw features with diagnostics
# information from the previous example, clipped to the period of paw activity
plot(features, clipped = TRUE)
```

![Diagnostics panel](/man/figures/diagnostics.svg)

The `clipped = TRUE` property subsets the frames viewed in each panel according
to the identified window of activity (the period between the blue dashed lines
in the previous example). Other than this windowing, the first four panels are
the same as in the previous example. The fifth panel shows the the univariate
projection time series, as in the first example. The sixth and final panel shows
this same univariate time series scaled by the maximum paw height attained, and
annotated with periods after the time of the first peak, $t^\star$, identified
as either "guarding" (gray) or "shaking" (red).

In the sixth panel, counts are annotated above the panel for the estimated
number of potential consecutive shakes in each period. Displacements above the
shaking threshold (set by the `shake.filter.threshold` parameter; see below) are
annotated by red points (these are points where displacement from previous local
minima or maxima exceeds the threshold).  Displacements below threshold are
annotated by black points.  Sequential runs of all below or all above threshold
displacements are shown as regions bordered by black dashed vertical lines. A
region is identified as shaking (colored red) if and only if at least two
consecutive displacements are above threshold. This example shows one such
shaking period with four "shakes", three periods with no shakes above threshold
(having 1, 2, and 2 consecutive displacements, respectively), and one period
with only a single displacement above threshold.

Different diagnostic panels can be shown by passing one of the following
arguments to the `panel` parameter of the `plot` function (by default all panels
will be shown):

| panel     | output                                              |
| :-------: | :-------------------------------------------------- |
| *"x"*     | horizontal displacement                             |
| *"y"*     | vertical displacement                               |
| *"vx"*    | estimated horizontal velocity                       |
| *"vy"*    | estimated vertical velocity                         |
| *"u"*     | univariate projection                               |
| *"s"*     | scaled univariate projection (shaking and guarding) |
| *"d"*     | all displacement panels (i.e. both "x" and "y")     |
| *"ve"*    | all velocity panels ("vx" and "vy")                 |
| *"k"*     | all kinematics panels ("x", "y", "vx", and "vy")    |
| *"p"*     | all projection panels ("u", "s")                    |
| *"a"*     | all panels (default)                                |

The first six options show a single panel in isolation. The remaining five
options show some combination of panels, which can be useful when checking
different parameter settings.

### Parameters

There are a number of parameters associated with the feature extraction process.
The defaults can be inspected and changed using the `default_parameters` and
`set_parameters` functions, respectively, as shown below.

```R
# returns the parameters used in Jones et al. 2020 (the default parameters used
# by the extract_parameters function)
default_parameters()

# an example of how to modify the default parameters (modifications are passed
# as arguments to set_parameters)
params <- set_parameters(fps = 1000, shake.filter.threshold = 0.4)

# the modified parameters can be passed to the extract_features function in
# place of the defaults
features <- extract_features(paw, parameters = params)
```

The default parameters assume that the `x` and `y` time series data passed to
`extract_features` is in millimeters. If this is not the case, the
`window.threshold` parameter should be scaled accordingly (this parameter is
used for determing when the paw first lifts off the ground, and when it last
returns to the ground). The list of parameters, their default values, units, and
what they are used for, is given below, along with the most informative
diagnostic panel for checking each value.

| parameter                | value | units           | panel   | use                                                              |
| -----------------------: | ----: | :-------------: | ------- | :--------------------------------------------------------------- |
| fps                      | 2000  | frames / second | &mdash; | conversion from seconds to frames                                |
| window.filter.size       | 0.045 | seconds         | *"k"*   | activity window identification                                   |
| window.filter.order      | 3     | integer[^2]     | *"k"*   | activity window identification                                   |
| window.threshold         | 0.5   | mm[^3]          | *"k"*   | activity window identification displacement threshold            |
| projection.window        | 0.04  | seconds         | *"u"*   | univariate projection sliding window size                        |
| velocity.filter.size     | 0.005 | seconds         | *"ve"*  | velocity estimation                                              |
| velocity.filter.order    | 3     | integer         | *"ve"*  | velocity estimation                                              |
| global.peak.filter.size  | 0.015 | seconds         | *"y"*   | maximum paw height                                               |
| global.peak.filter.order | 3     | integer         | *"y"*   | maximum paw height                                               |
| local.peak.filter.size   | 0.015 | seconds         | *"y"*   | first peak time, $t^\star$                                       |
| local.peak.filter.order  | 3     | integer         | *"y"*   | first peak time, $t^\star$                                       |
| local.peak.threshold     | 0.2   | %               | *"y"*   | first peak time, $t^\star$, (percentage of maximum paw height)   |
| shake.filter.size        | 0.015 | seconds         | *"s"*   | shaking and guarding sequences                                   |
| shake.filter.order       | 3     | integer         | *"s"*   | shaking and guarding sequences                                   |
| shake.filter.threshold   | 0.35  | %               | *"s"*   | shaking and guarding sequences (percentage of scaled projection) |

Many of these parameters can safely be held fixed across datasets (e.g. the
filter order parameters). Others may need to be changed based on the
characteristics of the noise in the tracking data, but are likely to have the
same value as each other (e.g.  `global.peak.filter.size`,
`local.peak.filter.size`, and `shake.filter.size` are likely to all have the
same value as each other, even if different from the default listed above). When
modifying parameters, it will be best to start at the top of this list and work
down to the bottom (the choice of parameters early in the list may affect later
choices; less so vice versa).

[^2]: The order number controls the order of the polynomial used in the
  Savitzky-Golay filter.
[^3]: The units for the Jones et al. default parameters are in millimeters
  because the tracking data in `jones2020.tracks` is in millimeters. If your time
  series data is in units other than millimeters, you can either convert to
  millimeters before passing it to `extract_features`, or modify the
  `window.threshold` parameter to match the units of your time series data.

### Credit

Please cite [Jones et al. (2020)](https://doi.org/10.1101/2020.02.18.955070)
and include a link to this repository if you use this code in an academic
publication.

Testing, bug reports, and code contributions very welcome.

Copyright (c) 2019 &ndash; 2023 Colin Twomey.
Shared under a GNU GPLv3 license (see COPYING).

