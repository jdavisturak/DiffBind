// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// mergePeaks
Rcpp::DataFrame mergePeaks(Rcpp::DataFrame data, int maxGap);
RcppExport SEXP DiffBind_mergePeaks(SEXP dataSEXP, SEXP maxGapSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type maxGap(maxGapSEXP);
    __result = Rcpp::wrap(mergePeaks(data, maxGap));
    return __result;
END_RCPP
}
// mergeScores
Rcpp::List mergeScores(Rcpp::DataFrame sMerged, Rcpp::NumericVector sScore, Rcpp::DataFrame sPeaks);
RcppExport SEXP DiffBind_mergeScores(SEXP sMergedSEXP, SEXP sScoreSEXP, SEXP sPeaksSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type sMerged(sMergedSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sScore(sScoreSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type sPeaks(sPeaksSEXP);
    __result = Rcpp::wrap(mergeScores(sMerged, sScore, sPeaks));
    return __result;
END_RCPP
}
// peakOrder
Rcpp::RObject peakOrder(SEXP schrom, SEXP sleft, SEXP sright);
RcppExport SEXP DiffBind_peakOrder(SEXP schromSEXP, SEXP sleftSEXP, SEXP srightSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type schrom(schromSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sleft(sleftSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sright(srightSEXP);
    __result = Rcpp::wrap(peakOrder(schrom, sleft, sright));
    return __result;
END_RCPP
}
