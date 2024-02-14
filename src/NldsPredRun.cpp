/*
 * Copyright (c) 2022 Roelof Bart Toonen
 * License: MIT license (spdx.org MIT)
 *
 * NldsPredRun.cpp
 *
 */

#include <Rcpp.h>
#include <exception>
#include <cstring>
#include <string.h>

extern "C" {
#include "traverse.h"
#include "tsdat.h"
#include "log.h"
#include "fn.h"
#include "fn_exp.h"
#include "fn_tls.h"
#include "point.h"
#include "sets.h"
#include "tstoembdef.h"
#include "bundle.h"
}

using namespace Rcpp;

char **
GetColNames (const NumericMatrix& m)
{
    SEXP dimnms = m.attr("dimnames");
    if ( !(!Rf_isNull(dimnms) && Rf_length(dimnms) > 1 && !Rf_isNull(VECTOR_ELT(dimnms, 1)) ) )
        throw (std::invalid_argument ("Column names of matrix not specified"));

    CharacterVector labs = colnames (m);

    char **clabs = new char * [labs.size()];

    for (int i = 0; i < labs.size(); i++)
    {
        clabs[i] = new char[(Rcpp::as<std::string>(labs[i])).size() + 1];
        strcpy (clabs[i], (Rcpp::as<std::string>(labs[i])).c_str());
    }

    return clabs;
}

double *
CopyNumericMatrix (const NumericMatrix& m)
{
    double * dm = new double[m.size()];

    for (int i = 0; i < m.nrow(); i++)
    {
        for (int j = 0; j < m.ncol(); j++)
        {
            if (R_IsNA (m[i, j]))
                dm[i * m.ncol() + j] = NAN;
            else
                dm[i * m.ncol() + j] = (double) m(i, j);
        }
    }

    return dm;
}

//' @name NldsPredFnTlsOpts
//'
//' @title Set function options for the total least square vector autoregression (TLS-VAR) function
//'
//' @description
//' Validate the parameter values and create a list containing parameters for the TLS-VAR function.
//' This list will be used as an input argument for NldsPredRun.
//'
//' @param thetaMin For the TLS-VAR function, theta controls the amount of ’localness’ 
//'                 of the VAR estimates With theta = 0.0, all library vectors
//'                 have an equal influence in the VAR estimate. For linear systems,
//'                 this should give the best results. For nonlinear systems, it may
//'                 be that the system shows local linear behavior In that case, in
//'                 the VAR parameter estimation, the nearest vectors should be
//'                 assigned a higher weight than more distant vectors. Higher
//'                 values of theta correspond to the estimates becoming more 
//'                 local. Usually, theta can range from 0.0 to values around 1.0 or
//'                 higher. In a validation procedure, different values of theta may
//'                 be used to find the best results. It is possible to let the program
//'                 go through a range of values for theta, ranging from thetaMin
//'                 to thetaMax, with steps of deltaTheta.
//'
//' @param thetaMax See thetaMin.
//'
//' @param deltaTheta See thetaMin.
//'
//' @param nnn      The TLS-VAR function normally uses all library vectors to 
//'                 compute the parameter values of the linear model. In the case of
//'                 very large time series this may take a lot of time. Using this 
//'                 argument it is possible to use only a specified amount of nearest
//'                 neighbor vectors.
//'
//' @param excl     Exclusion method Options: none, timeCoord, timeWin, self.
//'
//' @param varWin   When the exclusion method is set to timeWin, use varWin to
//'                 specify the size of the variable time window around the target.
//'                 The unit of measurement is the number of time series points
//'                 before or after the target.
//'
//' @param center   Center embedding around the predictee vector. The
//'                 coordinates of the predictee will all become 0. Options: TRUE, FALSE.
//'
//' @param          restrictPred If this argument is greater than zero then predictions that are
//'                 higher than this value will be discarded when calculating the
//'                 statistics. This only pertains to the statistics that are calculated
//'                 by NldsPred; not to additional statistics that are calculated afterwards 
//'                 by the user.
//'
//' @param warnIsError The TLS method may return warnings. If this argument is set
//'                 to TRUE, the computed values will be regarded as erroneous
//'                 and will be discarded. Options: TRUE, FALSE.
//'
//'
//' @param refMeth  Similarly to the exponential nearest neighbors function, TLS-VAR 
//'                 uses a reference distance in the denominator in an exponential
//'                 weight function. Values can be: mean - to use the average distance between
//'                 the predictee and all library vectors that
//'                 are used in the TLS-VAR procedure; and xnn - to use the average
//'                 distance between the predictee and xnn closest nearest neighbors.
//'
//' @param refXnn   If refMeth has been set to xnn then refXnn should be set to the
//'                 number of neighbors that have to be used to calculate the value
//'                 of the denominator.
//'
//' @return         This function returns a list of parameter settings.
//'
//' @examples
//' # Create a list of parameters for the TLS-VAR
//' # prediction function.
//' fo <- NldsPredFnTlsOpts(thetaMin = 0.0, thetaMax = 2.0,
//' deltaTheta = 0.1, excl = "timeCoord",
//' center = TRUE, warnIsError = TRUE,
//' refMeth = "mean")
//' @export
// [[Rcpp::export]]
List NldsPredFnTlsOpts (const double thetaMin = 0.0,
                        const double thetaMax = 2.0,
                        const double deltaTheta = 0.1,
                        const int nnn = 0,
                        const std::string excl = "timeCoord",
                        const int varWin = 0,
                        const bool center = true,
                        const double restrictPred = 0.0,
                        const bool warnIsError = false,
                        const std::string refMeth = "mean",
                        const int refXnn = 0,
                        const bool objectOnlyOnce = true
                       )
{
BEGIN_RCPP
    Texcl iexcl;
    Ttls_ref_meth iref_meth;

    if (thetaMin < 0.0)
        throw (std::invalid_argument ("Minimum value of theta cannot be less than 0. Current value " + std::to_string(thetaMin)));

    if (thetaMin > thetaMax)
        throw (std::invalid_argument ("Minimum value of theta cannot be greater than maximum value. Min = " +
                    std::to_string(thetaMin) + " max = " + std::to_string(thetaMax)));

    if (deltaTheta <= 0)
        throw (std::invalid_argument ("Delta theta cannot be less than or equal to 0. Current value " + 
                    std::to_string(deltaTheta)));

    if (nnn < 0)
        throw (std::invalid_argument ("Number of nearest neighbours cannot be negative. Current value " + std::to_string(nnn)));

    if (excl == "none")
        iexcl = T_EXCL_NONE;
    else if (excl == "timeCoord")
        iexcl = T_EXCL_TIME_COORD;
    else if (excl == "timeWin")
        iexcl = T_EXCL_TIME_WIN;
    else if (excl == "self")
        iexcl = T_EXCL_SELF;
    else
        throw (std::invalid_argument ("Exclusion type incorrect " + excl));

    if (iexcl == T_EXCL_TIME_WIN && varWin < 1)
        throw (std::invalid_argument ("Exclusion window should be larger than 0. " + std::to_string(varWin)));

    if (refMeth == "mean")
        iref_meth = KTLS_REFMETH_MEAN;
    else if (refMeth == "xnn")
        iref_meth = KTLS_REFMETH_XNN_GT_ZERO;
    else
        throw (std::invalid_argument ("Reference method for TLS estimation unknown. " + refMeth));

    if (iref_meth == KTLS_REFMETH_XNN_GT_ZERO && refXnn < 1)
        throw (std::invalid_argument ("TLS reference value should be based on at least one nearest neighbour. " + 
                    std::to_string(refXnn)));

    return List::create (Named("fnType") = "TLS",
                         Named("thetaMin") = thetaMin,
                         Named("thetaMax") = thetaMax,
                         Named("deltaTheta") = deltaTheta,
                         Named("nnn") = nnn,
                         Named("tlsCenter") = center,
                         Named("restrictPred") = restrictPred,
                         Named("warnIsError") = warnIsError,
                         Named("excl") = (int) iexcl,
                         Named("varWin") = varWin,
                         Named("refMeth") = (int) iref_meth,
                         Named("refXnn") = refXnn,
                         Named("objectOnlyOnce") = objectOnlyOnce);
END_RCPP
}

//' @name NldsPredFnExpOpts
//'
//' @title Set function options for exponential weighted nearest neighbor function
//'
//' @description
//' Validate the parameter values and create a list containing parameters for the exponential function that uses the nearest
//' neighbors to predict values. This list will be used as an input argument for NldsPredRun.
//'
//' @param expK     The constant value that is used in the exponential function.
//'
//' @param nnnAdd   To enclose a point in an n-dimensional space, at least n+1 nearest neighbors around it are needed.
//'                 This parameter specifies how many additional nearest neighbors to add to the dimension number n.
//'                 So, to use n+1 nearest neighbors, use a value of 1.
//'
//' @param excl     Exclusion method. Options: none, timeCoord, timeWin, self.
//'
//' @param varWin   When the exclusion method is set to timeWin, use varWin to specify the size of the variable time window 
//'                 around the target.
//'                 The unit of measurement is the number of time series points before or after the target.
//'
//' @param fnDenom  The method that is used to calculate the denominator in the fraction in the exponent of the function.
//'                 Options: min – the denominator is the distance between the closest nearest neighbor and the target;
//'                 avgNn – the denominator is the average of the distances of the nearest neighbors to the target;
//'                 avgLib – the denominator is the average distance between all points of the library set.
//'
//' @return         This function returns a list of parameter settings.
//'
//' @examples
//' # Create a list of parameters for the exponential
//' # prediction function.
//' fo <- NldsPredFnExpOpts (expK = 1.0, nnnAdd = 1, excl = "timeCoord",
//'                          fnDenom = "min")
//'
//' @export
// [[Rcpp::export]]
List NldsPredFnExpOpts (const double expK = 1.0,
                        const int nnnAdd = 1,
                        const std::string excl = "timeCoord",
                        const int varWin = 0,
                        const std::string fnDenom = "avgLib",
                        const bool objectOnlyOnce = true
                       )
{
BEGIN_RCPP
    Texcl iexcl;
    Tfn_denom ifn_denom;

    if (nnnAdd < 0)
        throw (std::invalid_argument ("Number of nearest neighbours to add cannot be negative. Current value " + 
                    std::to_string(nnnAdd)));

    if (excl == "none")
        iexcl = T_EXCL_NONE;
    else if (excl == "timeCoord")
        iexcl = T_EXCL_TIME_COORD;
    else if (excl == "timeWin")
        iexcl = T_EXCL_TIME_WIN;
    else if (excl == "self")
        iexcl = T_EXCL_SELF;
    else
        throw (std::invalid_argument ("Exclusion type incorrect " + excl));

    if (iexcl == T_EXCL_TIME_WIN && varWin < 1)
        throw (std::invalid_argument ("Exclusion window should be larger than 0. " + std::to_string(varWin)));

    if (fnDenom == "min")
        ifn_denom = FN_WEIGHT_DENOM_MINIMUM;
    else if (fnDenom == "avgNn")
        ifn_denom = FN_WEIGHT_DENOM_AVG_NN;
    else if (fnDenom == "avgLib")
        ifn_denom = FN_WEIGHT_DENOM_AVG_LIB;
    else
        throw (std::invalid_argument ("Invalid value for denominator of exponent. " + fnDenom));

    return List::create (Named("fnType") = "EXP",
                         Named("expK") = expK,
                         Named("excl") = (int) iexcl,
                         Named("varWin") = varWin,
                         Named("fnDenom") = (int) ifn_denom,
                         Named("nnnAdd") = nnnAdd,
                         Named("objectOnlyOnce") = objectOnlyOnce);
END_RCPP
}

//' @name NldsPredSetConvergent
//'
//' @title Create a list of options for the convergent library validation method
//'
//' @description
//' This function checks the arguments and creates a list of options for the convergent library validation method that is used in
//' Sugihara’s CCM method. This list is used as an input argument for the NldsPredRun function.
//'
//' @param libSizeMin   Minimum library size.
//'
//' @param libSizeMax   Maximum library size.
//'
//' @param libInc       Increment to use when increasing the library size.
//'
//' @param libIncFactor When the library size increases, it is often possible to use larger values for the increment that is used.
//'                     If libIncFactor is greater than 0.0, libInc will be multiplied by this factor after each increment.
//'
//' @param shiftMethod  The method that is used to create multiple library sets with the same size.
//'                     Use option "shift" to move the selection window through the complete set of embedding points.
//'                     If the right side of the window exceeds the embedding size, it will continue at the beginning.
//'                     Use option "random" to first put the embedding points in random order.
//'                     Moving is done similar to the "shift" option. Use option "bootstrap" to create library sets by sampling
//'                     from the embedding.
//'
//' @param libShift     If shiftMethod has been set to "shift" or "random", then this argument specifies the number of
//'                     points the selection window will move at each shift.
//'
//' @param nBootstrap   The number of bootstrap sets to create when shiftMethod has been set to "bootstrap".
//'
//' @return             This function returns a list of options.
//'
//' @examples
//' so <- NldsPredSetConvergent(libSizeMin = 10, libSizeMax = 85, libInc = 1,
//'                             libIncFactor = 1.1, shiftMethod = "random",
//'                             libShift = 1)
//'
//' so <- NldsPredSetConvergent(libSizeMin = 10,libSizeMax = 200, libInc = 1,
//'                             libIncFactor = 1.3, shiftMethod = "bootstrap",
//'                             nBootstrap = 3000)
//'
//' @export
// [[Rcpp::export]]
List NldsPredSetConvergent (const int libSizeMin,
                            const int libSizeMax,
                            const int libInc = 1,
                            const double libIncIncFactor = 0.0,
                            const std::string shiftMethod = "bootstrap",
                            const int libShift = 1,
                            const int nBootstrap = 100
                            )
{
BEGIN_RCPP
    Tlib_shift_meth ilib_shift_meth;

    if (libSizeMin < 1)
        throw (std::invalid_argument ("Minimum lib size should be greater than zero " + std::to_string(libSizeMin)));

    if (libSizeMax < 1)
        throw (std::invalid_argument ("Maximum lib size should be greater than zero " + std::to_string(libSizeMax)));

    if (libSizeMax < libSizeMin)
        throw (std::invalid_argument ("Maximum lib size should be greater than minimum lib size " + 
                    std::to_string(libSizeMax) + std::to_string(libSizeMin)));

    if (libInc < 1)
        throw (std::invalid_argument ("Increment should be greater than zero " + std::to_string(libInc)));

    if (libIncIncFactor < 0.0)
        throw (std::invalid_argument ("Increment increase factor cannot be negative. " + 
                    std::to_string(libIncIncFactor)));

    if (shiftMethod == "shift")
        ilib_shift_meth = LIB_SHIFT_SHIFT;
    else if (shiftMethod == "random")
        ilib_shift_meth = LIB_SHIFT_RANDOM;
    else if (shiftMethod == "bootstrap")
        ilib_shift_meth = LIB_SHIFT_BOOTSTRAP;
    else if (shiftMethod == "permut")
        ilib_shift_meth = LIB_SHIFT_BOOT_PERMUT;
    else
        throw (std::invalid_argument ("Shift method incorrect " + shiftMethod));

    return List::create (Named("setMethod") = "CONVERGENT",
                         Named("libSizeMin") = libSizeMin,
                         Named("libSizeMax") = libSizeMax,
                         Named("libInc") = libInc,
                         Named("libIncIncFactor") = libIncIncFactor,
                         Named("shiftMethod") = (int) ilib_shift_meth,
                         Named("libShift") = libShift,
                         Named("nBootstrap") = nBootstrap);
END_RCPP
}

//' @name NldsPredSetKFold
//'
//' @title Create a list of options for the k-fold validation method
//'
//' @description
//' This function checks the arguments and creates a list of options for the k-fold validation method.
//' This list is used as an input argument for the NldsPredRun function.
//'
//' @param kFold    The factor k for the k-fold validation method.
//'
//' @param nRep     Number of repetitions.
//'                 Before each repetition, the order of the points in the full set is reshuffled.
//'                 Then a new cycle of k-fold validations is carried out.
//'
//' @return         This function returns a list of options.
//'
//' @export
// [[Rcpp::export]]
List NldsPredSetKFold (const int kFold,
                       const int nRep
                      )
{
BEGIN_RCPP
    return List::create (Named("setMethod") = "KFOLD",
                         Named("kFold") = kFold,
                         Named("nRep") = nRep);
END_RCPP
}

//' @name NldsPredSetLooc
//'
//' @title Create a list of options for the LOOC validation method
//'
//' @description
//' This function creates a list of options for the LOOC validation method. This list is used as an input argument for the
//' NldsPredRun function.
//'
//' @return         This function returns a list of options.
//'
//' @export
// [[Rcpp::export]]
List NldsPredSetLooc ()
{
BEGIN_RCPP
    return List::create (Named("setMethod") = "LOOC");
END_RCPP
}

//' @name NldsPredSetBootstrap
//'
//' @title Create a list of options for the bootstrap validation method
//'
//' @description
//' This function checks the arguments and creates a list of options for the bootstrap validation method.
//' This list is used as an input argument for the NldsPredRun function.
//'
//' @param nBootstrap  Number of bootstrap sets to create.
//'
//' @param libSize     Number of points to sample into each bootstrap set.
//'
//' @param preSetIsLib If this argument is set to TRUE, the prediction set will contain the same points 
//'                    as the library set. Options: TRUE, FALSE.
//'
//' @param libSizeIsEmbSize If this argument is set to TRUE, the number of points in the library set will always 
//'                         be set equal to the number of points in the embedding. If bundle embeddings are used,
//'                         the sampled bootstrap set will contain as many vectors per bundle as the library set.
//'                         Options: TRUE, FALSE.
//'
//' @param perAdditGroup    When additional variables are used, and this argument is set to TRUE, the number of
//'                         vectors with a specific value for the additional variable will be the same as in the
//'                         library set.  Options: TRUE, FALSE.
//'
//' @return         This function returns a list of options.
//'
//' @examples
//' so <- NldsPredSetBootstrap(nBootstrap = 2000, preSetIsLib = FALSE,
//'                            libSize = 300)
//'
//' so <- NldsPredSetBootstrap(nBootstrap = 1000, preSetIsLib = FALSE,
//'                            libSizeIsEmbSize = TRUE)
//'
//' @export
// [[Rcpp::export]]
List NldsPredSetBootstrap (const int nBootstrap,
                           const int libSize = 0,
                           const bool preSetIsLib = false,
                           const bool libSizeIsEmbSize = true,
                           const bool perAdditGroup = false
                          )
{
BEGIN_RCPP
    return List::create (Named("setMethod") = "BOOTSTRAP",
                         Named("libSize") = libSize,
                         Named("nBootstrap") = nBootstrap,
                         Named("preSetIsLib") = preSetIsLib,
                         Named("libSizeIsEmbSize") = libSizeIsEmbSize,
                         Named("perAdditGroup") = perAdditGroup);
END_RCPP
}

//' @name NldsPredRun
//'
//' @title Starts the computations
//'
//' @description
//' This function starts the computations.
//'
//' @param data           Numeric matrix with time-series values.
//'                       Each column corresponds to one variable.
//'                       Time, id and bundle values are not included in the matrix.
//'                       These are supplied through the time, id or bundle argument respectively.
//'                       Column names should be set and correspond to the variable names that are 
//'                       used in the embLagDef or rangeDef arguments.
//'
//' @param logFile        Name of the log file where output will be written to.
//'
//' @param fnOpts         List of function options. Generated by NldsPredFnExpOpts or NldsPredFnTlsOpts.
//'
//' @param setOpts        List of validation set options. Generated by NldsPredSetBootstrap,
//'                       NldsPredSetConvergent, NldsPredSetKFold or NldsPredSetLooc.
//'
//' @param time           When user-supplied time values are required in the log file, use this argument
//'                       to supply a vector with integer values that can be transformed in real time
//'                       values afterwards by the user. The length of this vector must be equal to the
//'                       number of rows in the data matrix. Can be NULL.
//'
//' @param id             When dewdrop embeddings are used, use this argument to supply the ID strings of
//'                       each row in the data matrix. The length of this vector must be equal to the
//'                       number of rows in the data matrix. Can be NULL. If it is not NULL, then the
//'                       NldsPredRun will automatically use dewdrop embeddings.
//'
//' @param bundle         When bundles are used, use this argument to supply the bundle values for each
//'                       row in the data matrix. The number of rows in this numeric matrix must be equal
//'                       to the number of rows in the data matrix. Can be NULL.
//'
//' @param embLagDef      Embedding lag definition string. Can be NULL if embRangeDef is supplied.
//'
//' @param embRangeDef    Embedding range definition string. Can be NULL if embLagDef is supplied.
//'
//' @return               This function returns the name of the log file.
//'
//' @examples
//' # Create a list of parameters for the exponential prediction
//' # function.
//' fo <- NldsPredFnExpOpts (expK = 1.0, nnnAdd = 1, excl = "timeCoord", fnDenom = "min")
//'
//' # Create a list of validation set options for bootstrapping.
//' so <- NldsPredSetBootstrap (nBootstrap = 5000, libSizeIsEmbSize = TRUE)
//'
//' # Matrix ds has two columns, x and y, with normalized
//' # time-series data. ti is a vector with integer times,
//' # corresponding to the rows in ds. The rows in ds, and the
//' # values in ti, must be ordered by increasing time value.
//' # Here, x is predicted from y, using several embedding
//' # dimensions (defined by embRangeDef).
//' fname <- NldsPredRun (data = ds, logFile = tempfile(), fnOpts = fo,
//'                       setOpts = so, time = ti,
//'                       embRangeDef = "y,0,3,1,trail:x,0")
//'
//' # Extract coordinates of points used in each embedding.
//' embpoints <- NldsPredGetTable (logFileName = fname,
//'                                fieldNames = c("embpar1_emb_label",
//'                                               "embpar1_emb_num",
//'                                               "embpar1_n_row",
//'                                               "embpar1_e",
//'                                               "embpar1_n_pre_val",
//'                                               "embpar1_n_bundle_val",
//'                                               "vid_vec_num",
//'                                               "val_copr",
//'                                               "val_idx",
//'                                               "val_t",
//'                                               "val_val"),
//'                                newRecName = "VAL")
//'
//' # Extract general statistics that were generated by
//' # NldsPredRun.
//' stats <- NldsPredGetTable(logFileName = fname,
//'                           fieldNames = c("embpar1_emb_num",
//'                                          "embpar1_e",
//'                                          "spbt_boot_i",
//'                                          "stat_n_pre_obs",
//'                                          "stat_avg_obs",
//'                                          "stat_var_pre",
//'                                          "stat_var_obs",
//'                                          "stat_cov_pre_obs",
//'                                          "stat_mae_pre_obs",
//'                                          "stat_rmse_pre_obs",
//'                                          "stat_md_pre",
//'                                          "stat_md_obs",
//'                                          "stat_mdad_pre",
//'                                          "stat_mdad_obs",
//'                                          "stat_mdae_pre_obs"),
//'                           newRecName = "STAT")
//'
//' @export
// [[Rcpp::export]]
String NldsPredRun (const NumericMatrix & data,
                    const std::string logFile,
                    const List & fnOpts,
                    const List & setOpts,
                    const Rcpp::Nullable<Rcpp::IntegerVector> & time = R_NilValue,
                    const Rcpp::Nullable<Rcpp::StringVector> & id = R_NilValue,
                    const Rcpp::Nullable<Rcpp::NumericMatrix> & bundle = R_NilValue,
                    const std::string embLagDef = "",
                    const std::string embRangeDef = "",
                    const int log_level = 0x0FFF     // Only change when you know what you are doing.
                   )
{
BEGIN_RCPP
    Tnew_fn_params  new_fn_params;  // Pointer to function that initializes fn parameters.
    Tnext_fn_params next_fn_params; // Pointer to function to switch to the next set of fn parameters. 
    Tfn             fn;             // Pointer to the function that performs the computations.

    Texcl excl;
    Tfn_denom denom;

    Tnew_sets new_sets;
    Tnext_set next_set;
    Tfree_set free_set;

    Temb_lag_def * emb_lag_def;
    int          n_lag_def;

    Tfdat * fdat = new Tfdat();

    int ret_val;
    StringVector fnType;


    fdat->lab = GetColNames (data);
    fdat->in_lib = NULL; // We don't use it.
    fdat->in_pre = NULL; // We don't use it.



    // Check the parameters and convert to C type struct.
    fdat->t = NULL;
    if (time.isNotNull())
    {
        IntegerVector _time(time);
        if (_time.size() != data.nrow())
            throw (std::invalid_argument ("Size of time vector does not match data matrix"));

        fdat->t = new long[_time.size()];

        for (int i = 0; i < _time.size(); i++)
            fdat->t[i] = _time[i];
    }

    fdat->id = NULL;

    if (id.isNotNull())
    {
        StringVector _id(id);
        if (_id.size() != data.nrow())
            throw (std::invalid_argument ("Size of id vector does not match data matrix"));

        fdat->id = new char * [_id.size()];

        for (int i = 0; i < _id.size(); i++)
        {
            fdat->id[i] = new char[(Rcpp::as<std::string>(_id[i])).size() + 1];
            strcpy (fdat->id[i], (Rcpp::as<std::string>(_id[i])).c_str());
        }
    }

    fdat->n_bundle_col = 0;
    fdat->bundle = NULL;
    fdat->bundle_lab = NULL;

    if (bundle.isNotNull())
    {
        NumericMatrix _bundle(bundle);
        if (_bundle.nrow() != data.nrow())
            throw (std::invalid_argument ("Number of rows in bundle matrix does not match data matrix"));

        fdat->n_bundle_col = _bundle.ncol();
        fdat->bundle_lab   = GetColNames (_bundle);
        fdat->bundle       = CopyNumericMatrix (_bundle);
    }

    fdat->n_col = data.ncol();
    fdat->n_dat = data.nrow();
    fdat->dat   = CopyNumericMatrix (data);

    // Open the binary logfile where all output will be written to. 0x0FFF = log everything.
    if (open_bin_log_file (logFile.c_str(), log_level) < 0)
        throw (std::invalid_argument ("Could not open the binary logfile: " + logFile));


    // Initialize the method to predict values.
    fnType = fnOpts["fnType"];

    if (fnType[0] == "EXP")
    {
        IntegerVector nnnAdd, varWin, excl, fnDenom;
        NumericVector expK;
        LogicalVector objectOnlyOnce;

        nnnAdd         = fnOpts["nnnAdd"];
        excl           = fnOpts["excl"];
        varWin         = fnOpts["varWin"];
        fnDenom        = fnOpts["fnDenom"];
        expK           = fnOpts["expK"];
        objectOnlyOnce = fnOpts["objectOnlyOnce"];

        if ((ret_val = init_fn_exponential (&new_fn_params, &next_fn_params, &fn,
                        (int) nnnAdd[0], (Texcl) ((int) excl[0]), (int) varWin[0], (Tfn_denom) ((int) fnDenom[0]),
                        (double) expK[0], (bool) objectOnlyOnce[0] )) < 0)
            throw (std::invalid_argument ("Error when initializing the exponential fn " + std::to_string(ret_val)));
    } else if (fnType[0] == "TLS")
    {
        IntegerVector nnn, varWin, refMeth, refXnn, excl;
        NumericVector thetaMin, thetaMax, deltaTheta, restrictPred;
        LogicalVector tlsCenter, warnIsError, objectOnlyOnce;

        thetaMin       = fnOpts["thetaMin"];
        thetaMax       = fnOpts["thetaMax"];
        deltaTheta     = fnOpts["deltaTheta"];
        nnn            = fnOpts["nnn"];
        excl           = fnOpts["excl"];
        varWin         = fnOpts["varWin"];
        tlsCenter      = fnOpts["tlsCenter"];
        restrictPred   = fnOpts["restrictPred"];
        warnIsError    = fnOpts["warnIsError"];
        refMeth        = fnOpts["refMeth"];
        refXnn         = fnOpts["refXnn"];
        objectOnlyOnce = fnOpts["objectOnlyOnce"];

        if ((ret_val = init_fn_tls (&new_fn_params, &next_fn_params, &fn,
                        (double) thetaMin[0],
                        (double) thetaMax[0],
                        (double) deltaTheta[0],
                        (int) nnn[0],
                        (Texcl) ((int) excl[0]),
                        (int) varWin[0],
                        (bool) tlsCenter[0],
                        (double) restrictPred[0],
                        (bool) warnIsError[0],
                        (Ttls_ref_meth) ((int) refMeth[0]),
                        (int) refXnn[0],
                        (bool) objectOnlyOnce[0] )) < 0)
            throw (std::invalid_argument ("Error when initializing the total least squares  fn " + std::to_string(ret_val)));
    } else
        throw (std::invalid_argument ("Value of function type not correct: " + fnType[0]));

    // Initialize the method to generate library and predictee sets.
    StringVector setMethod = setOpts["setMethod"];

    if (setMethod[0] == "CONVERGENT")
    {
        IntegerVector libSizeMin, libSizeMax, libInc, shiftMethod, libShift, nBootstrap;
        NumericVector libIncIncFactor;

        libSizeMin      = setOpts["libSizeMin"];
        libSizeMax      = setOpts["libSizeMax"];
        libInc          = setOpts["libInc"];
        shiftMethod     = setOpts["shiftMethod"];
        libShift        = setOpts["libShift"];
        nBootstrap      = setOpts["nBootstrap"];
        libIncIncFactor = setOpts["libIncIncFactor"];

        if ((ret_val = init_set_convergent_lib (
                        (int) libSizeMin[0],
                        (int) libSizeMax[0],
                        (int) libInc[0],
                        (double) libIncIncFactor[0],
                        (Tlib_shift_meth) ((int) shiftMethod[0]),
                        (int) libShift[0],
                        (int) nBootstrap[0],
                        &new_sets, &next_set, &free_set)) < 0)
            throw (std::invalid_argument ("Error when initializing convergent sets."));

    } else if (setMethod[0] == "KFOLD")
    {
        IntegerVector kFold, nRep;

        kFold = setOpts["kFold"];
        nRep  = setOpts["nRep"];

        if ((ret_val = init_set_k_fold ((int) kFold[0], (int) nRep[0], &new_sets, &next_set, &free_set)) < 0)
            throw (std::invalid_argument ("Error when initializing K fold validation sets."));

    } else if (setMethod[0] == "LOOC")
    {
        if ((ret_val = init_set_looc (&new_sets, &next_set, &free_set)) < 0)
            throw (std::invalid_argument ("Error when initializing looc validation sets."));

    } else if (setMethod[0] == "BOOTSTRAP")
    {
        IntegerVector libSize, nBootstrap;
        LogicalVector preSetIsLib, libSzIsEmbedSz, perAdditGroup;

        libSize        = setOpts["libSize"];
        nBootstrap     = setOpts["nBootstrap"];
        preSetIsLib    = setOpts["preSetIsLib"];
        libSzIsEmbedSz = setOpts["libSizeIsEmbSize"];
        perAdditGroup  = setOpts["perAdditGroup"];

        if ((ret_val = init_set_bootstrap ((int) libSize[0], (int) nBootstrap[0],
                                           (bool) preSetIsLib[0], (bool) libSzIsEmbedSz[0],
                                           (bool) perAdditGroup[0],
                                           &new_sets, &next_set, &free_set)) < 0)
            throw (std::invalid_argument ("Error when initializing bootstrap validation sets."));

    } else
        throw (std::invalid_argument ("Value of set method not correct: " + setMethod[0]));

    if (!embLagDef.empty())
    {
        char * emb_str = new char [embLagDef.size() + 1];
        strncpy (emb_str, embLagDef.c_str(), embLagDef.size() + 1);
        if ((emb_lag_def = str_to_lag_def (emb_str, &n_lag_def)) == NULL)
            throw (std::invalid_argument ("Error when initializing embedding lag definitions: " + embLagDef));

    } else if (!embRangeDef.empty())
    {
        char * emb_str = new char [embRangeDef.size() + 1];
        strncpy (emb_str, embRangeDef.c_str(), embRangeDef.size() + 1);
        if ((emb_lag_def = str_range_to_lag_def (emb_str, &n_lag_def)) == NULL)
            throw (std::invalid_argument ("Error when initializing embedding range definitions: " + embRangeDef));

    } else
        throw (std::invalid_argument ("No embedding parameters defined."));

    ret_val = traverse_all (fdat, n_lag_def, emb_lag_def,
                         new_sets, next_set, free_set,
                         new_fn_params, next_fn_params, fn,
                         true, true);

    return logFile;
END_RCPP
}
