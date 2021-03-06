# Copyright (c) 2022 Roelof Bart Toonen
# License: MIT license (spdx.org MIT)
# NldsPredExamples.R

require (NldsPred)
require (dplyr)

# compLV
# Simulate competitive Lotka Volterra system
# r1, a11 and a12 are the parameters for species 1.
# r2, r21 and a22 are the parameters for species 2.
# sd1m and sd2m are the standard deviations of model noise for species 1 and 2.
# sd1o and sd2o are the standard deviations of observational noise for species 1 and 2.

compLV <- function (x1, r1, a11, a12, y1, r2, a21, a22, sd1m, sd1o, sd2m, sd2o, n)
{
    set.seed (10)
    # Create time series for 2 variable model + stochastic influence.
    makeNoise <- function (s, n)
    {
        if (s == 0.0)
            return (rep (0.0, n))
        else
            return (rnorm (n, mean = 0.0, sd = s))
    }

    xmnoise = makeNoise (sd1m, n)
    ymnoise = makeNoise (sd2m, n)
    xonoise = makeNoise (sd1o, n)
    yonoise = makeNoise (sd2o, n)

    x  = rep(NA, n)
    y  = rep(NA, n)
    ti = 1:n

    x[1] = x1 + xmnoise[1]
    y[1] = y1 + ymnoise[1]

    for(t in 2:n)
    {
        x[t] = x[t-1]*(1 + r1 - a11*x[t-1] - a12*y[t-1]) + xmnoise[t]
        y[t] = y[t-1]*(1 + r2 - a21*y[t-1] - a22*x[t-1]) + ymnoise[t]
    }

    x = x + xonoise
    y = y + yonoise

    ts = data.frame(ti, x, y)

    return (ts)
}

# Create a dataset to work with
data1 = compLV (x1 = 0.2, r1 = 2.7, a11 = 3.7, a12 = 0.05,
                y1 = 0.4, r2 = 2.7, a21 = 3.7, a22 = 0.38,
                sd1m = 0.0, sd1o = 0.00, sd2m = 0.0, sd2o = 0.00, 300)

###############################################################################
# Example 1:
# Predict y from x, using a weighted mean of the y values of the neighbors in 
# the embedding for variable x.

# Create a list with parameters for the exponential prediction function.
fo <- NldsPredFnExpOpts (expK = 1.0, nnnAdd = 1, excl = "timeCoord", fnDenom = "min")

# Create a list with validation set options for bootstrapping.
so <- NldsPredSetBootstrap (nBootstrap = 5000, libSizeIsEmbSize = TRUE)

# Create normalized matrix for the time-series columns.
ds <- scale (data.matrix (data.frame(x = data1$x, y = data1$y), rownames.force = NA),
             center = TRUE,
             scale = TRUE)

fname <- NldsPredRun (data = ds, logFile = tempfile(), fnOpts = fo, setOpts = so,
                      time = data1$ti, embRangeDef = "y,0,3,1,trail:x,0")

# Extract coordinates of points used in each embedding.
embpoints <- NldsPredGetTable (logFileName = fname,
                 fieldNames = c("embpar1_emb_label", "embpar1_emb_num",
                                "embpar1_n_row", "embpar1_e",
                                "embpar1_n_pre_val", "embpar1_n_bundle_val",
                                "vid_vec_num", "val_copr", "val_idx", "val_t",
                                "val_val"),
                 newRecName = "VAL")


# Extract general statistics that were generated by NldsPredRun.
stats <- NldsPredGetTable(logFileName = fname,
             fieldNames = c("embpar1_emb_num", "embpar1_e", "spbt_boot_i",
                            "stat_n_pre_obs", "stat_avg_obs", "stat_var_pre",
                            "stat_var_obs", "stat_cov_pre_obs",
                            "stat_mae_pre_obs", "stat_rmse_pre_obs",
                            "stat_md_pre", "stat_md_obs", "stat_mdad_pre",
                            "stat_mdad_obs", "stat_mdae_pre_obs"),
             newRecName = "STAT")
View(stats)
