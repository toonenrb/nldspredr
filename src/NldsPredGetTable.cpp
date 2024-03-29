/*
 * Copyright (c) 2022 Roelof Bart Toonen
 * License: MIT license (spdx.org MIT)
 *
 * NldsPredGetTable.cpp
 *
 */

#include <Rcpp.h>
#include <exception>
#include <cstring>
#include <cstdio>

extern "C" {
#include "log_meta.h"
#include "logtotbl.h"
}

using namespace Rcpp;

StringVector
charArrayToStringVector (char *ca, long ni)
{
    StringVector sv(ni);
    char         s[2] = " ";

    for (int i = 0; i < ni; i++)
    {
        s[0]  = ca[i];
        sv[i] = s;
    }
    return sv;
}

StringVector
stringArrayToStringVector (char ** sa, long ni)
{
    StringVector sv(ni);

    for (int i = 0; i < ni; i++)
        sv[i] = sa[i];

    return sv;
}

IntegerVector
intArrayToIntegerVector (int * ia, long ni)
{
    IntegerVector iv(ni);

    for (int i = 0; i < ni; i++)
        iv[i] = ia[i];

    return iv;
}

IntegerVector
longArrayToIntegerVector (long * la, long ni)
{
    IntegerVector iv(ni);

    for (int i = 0; i < ni; i++)
        iv[i] = la[i];

    return iv;
}

IntegerVector
shortArrayToIntegerVector (short * sa, long ni)
{
    IntegerVector iv(ni);

    for (int i = 0; i < ni; i++)
        iv[i] = sa[i];

    return iv;
}

NumericVector
doubleArrayToNumericVector (double * da, long ni)
{
    NumericVector nv(ni);

    for (int i = 0; i < ni; i++)
        nv[i] = da[i];

    return nv;
}

//' @name NldsPredGetTable
//' 
//' @title Read the log file and return a data frame with the required fields.
//'
//' @description
//' This function reads the log file that was generated by NldsPredRun. It returns a data frame with columns that were specified
//' by the fieldNames argument.
//'
//' @param logFileName   Name of the log file that was created by NldsPredRun.
//'
//' @param fieldNames    String vector with names of fields to include in the data frame.
//'
//' @param newRecName    Name of the lowest level record to include in one data frame row.
//'                      After this record has been encountered in the log file and included
//'                      in the current data-frame row, a new data-frame row will be created.
//'
//' @return              This function returns a data frame object.
//'
//' @examples
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
//'                               newRecName = "VAL")
//'
//' @export
// [[Rcpp::export]]
DataFrame NldsPredGetTable (const std::string logFileName,
                            const StringVector & fieldNames,
                            const std::string newRecName)
{
BEGIN_RCPP
    FILE     * logfile = NULL;
    char     ** field_names;
    Ttbl_rec * tbl_rec;
    

    logfile = fopen (logFileName.c_str(), "r");

    if (!logfile)
        throw (std::invalid_argument ("Could not open the specified binary nlds log file: " + logFileName));
    // Convert StringVector to char * array.
    if (fieldNames.size() < 1)
        throw (std::invalid_argument ("No field names specified."));

    field_names = new char * [fieldNames.size() + 1];  /* One extra for "" string to close list of field names.*/

    for (int i = 0; i < fieldNames.size(); i++)
    {
        field_names[i] = new char [(Rcpp::as<std::string>(fieldNames[i])).size() + 1];
        strcpy ( field_names[i], (Rcpp::as<std::string>(fieldNames[i])).c_str() );
    }

    field_names[fieldNames.size()] = new char[1];
    *field_names[fieldNames.size()] = '\0';

    if (newRecName.size() < 2)
        throw (std::invalid_argument ("No rec name for new line specified: " + newRecName));

    tbl_rec = log_to_tbl (logfile, field_names, (char *) newRecName.c_str());

    if (!tbl_rec)
        throw (std::runtime_error ("No tables were created."));

    // Convert the returned C structure to an R list that appears as a data frame.
    // TODO: directly create R list in table convert program.

    int  ncol = tbl_rec->ncol;
    long  nrow = tbl_rec->nrow;
    CharacterVector colNames(ncol);
    List listOut(ncol);

    for (int i = 0; i < ncol; i++)
        colNames[i] = std::string(tbl_rec->col[i].name);

    listOut.attr("names") = colNames;
    listOut.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, nrow);
    listOut.attr("class") = "data.frame";

    for (int i = 0; i < ncol; i++)
    {
        switch (tbl_rec->col[i].field_type)
        {
          case FT_CHAR:
            listOut[i] = charArrayToStringVector ((char *) tbl_rec->col[i].val, nrow);
            break;
          case FT_STRING:
            listOut[i] = stringArrayToStringVector ((char **) tbl_rec->col[i].val, nrow);
            break;
          case FT_INT:
            listOut[i] = intArrayToIntegerVector ((int *) tbl_rec->col[i].val, nrow);
            break;
          case FT_SHORT:
            listOut[i] = shortArrayToIntegerVector ((short *) tbl_rec->col[i].val, nrow);
            break;
          case FT_LONG:
            listOut[i] = longArrayToIntegerVector ((long *) tbl_rec->col[i].val, nrow);
            break;
          case FT_DOUBLE:
            listOut[i] = doubleArrayToNumericVector ((double *) tbl_rec->col[i].val, nrow);
            break;
        }
    }
    free_tbl();

    return listOut;
END_RCPP
}
