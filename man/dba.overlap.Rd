\name{dba.overlap}
\alias{dba.overlap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute binding site overlaps (occupancy analysis)
}
\description{
Computes binding overlaps and co-occupancy statistics 
}
\usage{
dba.overlap(DBA, mask, mode=DBA_OLAP_PEAKS, 
            contrast, method=DBA$config$AnalysisMethod, th=DBA$config$th, 
            bUsePval=DBA$config$bUsePval, 
            report, byAttribute, bCorOnly=TRUE, CorMethod="pearson", 
            DataType=DBA$config$DataType)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\emph{Required arguments:}


\item{DBA}{
DBA object
}

%\emph{Optional/defaulted arguments:}
 
\item{mask}{
mask or vector of peakset numbers indicating a subset of peaksets to use (see \code{\link{dba.mask}}). When generating overlapping/unique peaksets, either two, three, or four peaksets may be specified. If the mode type is DBA_OLAP_ALL, and a contrast is specified, a value of TRUE (mask=TRUE) indicates that all samples should be included (otherwise only those present in one of the contrast groups will be included).
}

\item{mode}{
indicates which results should be returned (see MODES below). One of:
\itemize{
\item \link{DBA_OLAP_PEAKS}
\item \link{DBA_OLAP_ALL}
\item \link{DBA_OLAP_RATE}
}
}

  \item{contrast}{
contrast number to use. Only specified if contrast data is to be used when mode=DBA_OLAP_ALL. 
See \code{\link{dba.show}(DBA, bContrast=T)} to get contrast numbers. 
}
  \item{method}{
if contrast is specified and mode=DBA_OLAP_ALL, use data from method used for analysis:
\itemize{ 
\item \link{DBA_DESEQ2}
\item \link{DBA_DESEQ2_BLOCK}
\item \link{DBA_EDGER}
\item \link{DBA_EDGER_BLOCK}
}
}
\item{th}{
if contrast is specified and mode=DBA_OLAP_ALL, significance threshold; all sites with FDR (or p-values, see \code{bUsePval}) 
less than or equal to this value will be included. 
A value of 1 will include all binding sites, but only the samples included in the contrast.
}
\item{bUsePval}{
if contrast is specified and mode=DBA_OLAP_ALL, logical indicating whether to use FDR (FALSE) or p-value (TRUE) for thresholding.
}
\item{report}{
if contrast is specified and mode=DBA_OLAP_ALL, a report (obtained from \code{\link{dba.report}}) 
specifying the data to be used. 
If counts are included in the report (and a contrast is specified), 
the count data from the report will be used to compute correlations, 
rather than the scores in the global binding affinity matrix. 
If report is present, the \code{method}, \code{th}, and \code{bUsePval} parameters are ignored.
}
\item{byAttribute}{
when computing co-occupancy statistics (DBA_OLAP_ALL), limit comparisons to peaksets with the same value for a specific attribute, one of:
\itemize{
\item \link{DBA_ID}
\item \link{DBA_TISSUE}
\item \link{DBA_FACTOR}
\item \link{DBA_CONDITION}
\item \link{DBA_TREATMENT}
\item \link{DBA_REPLICATE}
\item \link{DBA_CONSENSUS}
\item \link{DBA_CALLER}
}
}

\item{bCorOnly}{
when computing co-occupancy statistics (DBA_OLAP_ALL), logical indicating that only correlations, and not overlaps, should be computed. 
This is much faster if only correlations are desired (e.g. to plot the correlations using  \code{\link{dba.plotHeatmap}}).
}

\item{CorMethod}{
when computing co-occupancy statistics (DBA_OLAP_ALL), method to use when computing correlations.
}

\item{DataType}{
if mode==DBA_OLAP_PEAKS, the class of object that peaksets should be returned as:
\itemize{
\item \link{DBA_DATA_GRANGES}
\item \link{DBA_DATA_RANGEDDATA}
\item \link{DBA_DATA_FRAME}
}
Can be set as default behavior by setting \code{DBA$config$DataType}.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~

MODE: Generate overlapping/unique peaksets:

dba.overlap(DBA, mask, mode=DBA_OLAP_PEAKS, minVal)

MODE: Compute correlation and co-occupancy statistics (e.g. for dba.plotHeatmap):

dba.overlap(DBA, mask, mode=DBA_OLAP_ALL, byAttribute, minVal,
            attributes, bCorOnly, CorMethod)

MODE: Compute correlation and co-occupancy statistics using significantly differentially bound sites (e.g. for dba.plotHeatmap):

dba.overlap(DBA, mask, mode=DBA_OLAP_ALL, byAttribute, minVal,
            contrast, method, th=, bUsePval,
            attributes, bCorOnly, CorMethod)
       
Note that the scores from the global binding affinity matrix will be used for correlations unless a report containing count data is specified.

MODE: Compute overlap rates at different stringency thresholds:

dba.overlap(DBA, mask, mode=DBA_OLAP_RATE, minVal)
}
\value{
Value depends on the mode specified in the mode parameter.

If mode = DBA_OLAP_PEAKS, Value is an overlap record: a list of three peaksets for an A-B overlap,  seven peaksets for a A-B-C overlap, and fifteen peaksets for a A-B-C-D overlap:
\item{inAll}{peaks in all peaksets}
\item{onlyA}{peaks unique to peakset A}
\item{onlyB}{peaks unique to peakset B}
\item{onlyC}{peaks unique to peakset C}
\item{onlyD}{peaks unique to peakset D}
\item{notA}{peaks in all peaksets except peakset A}
\item{notB}{peaks in all peaksets except peakset B}
\item{notC}{peaks in all peaksets except peakset C}
\item{notD}{peaks in all peaksets except peakset D}
\item{AandB}{peaks in  peaksets A and B but not in peaksets C or D}
\item{AandC}{peaks in  peaksets A and C but not in peaksets B or D}
\item{AandD}{peaks in  peaksets A and D but not in peaksets B or C}
\item{BandC}{peaks in  peaksets B and C but not in peaksets A or D}
\item{BandD}{peaks in  peaksets B and D but not in peaksets A or C}
\item{CandD}{peaks in  peaksets C and D but not in peaksets A or B}


If mode = DBA_OLAP_ALL, Value is a correlation record: a matrix with a row for each pair of peaksets and the following columns:
\item{A}{peakset number of first peakset in overlap}
\item{B}{peakset number of second peakset in overlap}
\item{onlyA}{number of sites unique to peakset A}
\item{onlyB}{number of sites unique to peakset B}
\item{inAll}{number of peaks in both peakset A and B (merged)}
\item{R2}{correlation value A vs B}
\item{Overlap}{percentage overlap (number of overlapping sites divided by number of peaks unique to smaller peakset}


If mode = DBA_OLAP_RATE, Value is a vector whose length is the number of peaksets, containing the number of overlapping peaks at the corresponding minOverlaps threshold (i.e., Value[1] is the total number of unique sites, Value[2] is the number of unique sites appearing in at least two peaksets, Value[3] the number of sites overlapping in at least three peaksets, etc.).
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Rory Stark
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{dba.plotVenn}},  \code{\link{dba.plotHeatmap}}
}

\examples{
data(tamoxifen_peaks)
# default mode: DBA_OLAP_PEAKS -- get overlapping/non overlapping peaksets
mcf7 <- dba.overlap(tamoxifen,tamoxifen$masks$MCF7&tamoxifen$masks$Responsive)
names(mcf7)
mcf7$inAll

# mode:  DBA_OLAP_ALL -- get correlation record
mcf7 <- dba(tamoxifen,tamoxifen$masks$MCF7)
mcf7.corRec <- dba.overlap(mcf7,mode=DBA_OLAP_ALL,bCorOnly=FALSE)
mcf7.corRec

# mode: DBA_OLAP_RATE -- get overlap rate vector
data(tamoxifen_peaks)
rate <- dba.overlap(tamoxifen, mode=DBA_OLAP_RATE)
rate
plot(rate,type='b',xlab="# peaksets",ylab="# common peaks",
     main="Tamoxifen dataset overlap rate")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
