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
                        const int refXnn = 0
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
                         Named("refXnn") = refXnn);
END_RCPP
}

// [[Rcpp::export]]
List NldsPredFnExpOpts (const double expK = 1.0,
                        const int nnnAdd = 1,
                        const std::string excl = "timeCoord",
                        const int varWin = 0,
                        const std::string fnDenom = "avgLib"
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
                         Named("nnnAdd") = nnnAdd);
END_RCPP
}

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

// [[Rcpp::export]]
List NldsPredSetLooc ()
{
BEGIN_RCPP
    return List::create (Named("setMethod") = "LOOC");
END_RCPP
}

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

// [[Rcpp::export]]
String NldsPredRun (const NumericMatrix & data,
                    const std::string logFile,
                    const List & fnOpts,
                    const List & setOpts,
                    const Rcpp::Nullable<Rcpp::IntegerVector> & time = R_NilValue,
                    const Rcpp::Nullable<Rcpp::StringVector> & id = R_NilValue,
                    const Rcpp::Nullable<Rcpp::NumericMatrix> & bundle = R_NilValue,
                    const std::string embLagDef = "",
                    const std::string embRangeDef = ""
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

    // Open the binary logfile where all output will be written to. 0xFFFF = log everything.
    if (open_bin_log_file (logFile.c_str(), 0xFFFF) < 0)
        throw (std::invalid_argument ("Could not open the binary logfile: " + logFile));


    // Initialize the method to predict values.
    fnType = fnOpts["fnType"];

    if (fnType[0] == "EXP")
    {
        IntegerVector nnnAdd, varWin, excl, fnDenom;
        NumericVector expK;

        nnnAdd      = fnOpts["nnnAdd"];
        excl        = fnOpts["excl"];
        varWin      = fnOpts["varWin"];
        fnDenom     = fnOpts["fnDenom"];
        expK        = fnOpts["expK"];

        if ((ret_val = init_fn_exponential (&new_fn_params, &next_fn_params, &fn,
                        (int) nnnAdd[0], (Texcl) ((int) excl[0]), (int) varWin[0], (Tfn_denom) ((int) fnDenom[0]),
                        (double) expK[0] )) < 0)
            throw (std::invalid_argument ("Error when initializing the exponential fn " + std::to_string(ret_val)));
    } else if (fnType[0] == "TLS")
    {
        IntegerVector nnn, varWin, refMeth, refXnn, excl;
        NumericVector thetaMin, thetaMax, deltaTheta, restrictPred;
        LogicalVector tlsCenter, warnIsError;

        thetaMin     = fnOpts["thetaMin"];
        thetaMax     = fnOpts["thetaMax"];
        deltaTheta   = fnOpts["deltaTheta"];
        nnn          = fnOpts["nnn"];
        excl         = fnOpts["excl"];
        varWin       = fnOpts["varWin"];
        tlsCenter    = fnOpts["tlsCenter"];
        restrictPred = fnOpts["restrictPred"];
        warnIsError  = fnOpts["warnIsError"];
        refMeth      = fnOpts["refMeth"];
        refXnn       = fnOpts["refXnn"];

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
                        (int) refXnn[0] )) < 0)
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
