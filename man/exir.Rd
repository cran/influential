% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{exir}
\alias{exir}
\alias{ExIR}
\title{Experimental data-based Integrated Ranking}
\usage{
exir(
  Desired_list = NULL,
  Diff_data,
  Diff_value,
  Regr_value = NULL,
  Sig_value,
  Exptl_data,
  Condition_colname,
  Normalize = FALSE,
  cor_thresh_method = "mr",
  r = 0.5,
  mr = 20,
  max.connections = 50000,
  alpha = 0.05,
  num_trees = 10000,
  mtry = NULL,
  num_permutations = 100,
  inf_const = 10^10,
  ncores = "default",
  seed = 1234,
  verbose = TRUE
)
}
\arguments{
\item{Desired_list}{(Optional) A character vector of your desired features. This vector could be, for
instance, a list of features obtained from cluster analysis, time-course analysis,
or a list of dysregulated features with a specific sign.}

\item{Diff_data}{A dataframe of all significant differential/regression data and their
statistical significance values (p-value/adjusted p-value). Note that the differential data
should be in the log fold-change (log2FC) format.
You may have selected a proportion of the differential data as the significant ones according
to your desired thresholds. A function, named \code{\link[influential]{diff_data.assembly}}, has also been
provided for the convenient assembling of the Diff_data dataframe.}

\item{Diff_value}{An integer vector containing the column number(s) of the differential
data in the Diff_data dataframe. The differential data could result from any type of
differential data analysis. One example could be the fold changes (FCs) obtained from differential
expression analyses. The user may provide as many differential data as he/she wish.}

\item{Regr_value}{(Optional) An integer vector containing the column number(s) of the regression
data in the Diff_data dataframe. The regression data could result from any type of regression
data analysis or other analyses such as time-course data analyses that are based on regression models.}

\item{Sig_value}{An integer vector containing the column number(s) of the significance values (p-value/adjusted p-value) of
both differential and regression data (if provided). Providing significance values for the regression data is optional.}

\item{Exptl_data}{A dataframe containing all of the experimental data including a column for specifying the conditions.
The features/variables of the dataframe should be as the columns and the samples should come in the rows.
The condition column should be of the character class. For example, if the study includes several replicates of
cancer and normal samples, the condition column should include "cancer" and "normal" as the conditions of different samples.
Also, the prior normalization of the experimental data is highly recommended. Otherwise,
the user may set the Normalize argument to TRUE for a simple log2 transformation of the data.
The experimental data could come from a variety sources such as transcriptomics and proteomics assays.}

\item{Condition_colname}{A string or character vector specifying the name of the column "condition" of the Exptl_data dataframe.}

\item{Normalize}{Logical; whether the experimental data should be normalized or not (default is FALSE). If TRUE, the
experimental data will be log2 transformed.}

\item{cor_thresh_method}{A character string indicating the method for filtering the correlation results, either
"mr" (default; Mutual Rank) or "cor.coefficient".}

\item{r}{The threshold of Spearman correlation coefficient for the selection of correlated features (default is 0.5).}

\item{mr}{An integer determining the threshold of mutual rank for the selection of correlated features (default is 20). Note that
higher mr values considerably increase the computation time.}

\item{max.connections}{The maximum number of connections to be included in the association network.
Higher max.connections might increase the computation time, cost, and accuracy of the results (default is 50,000).}

\item{alpha}{The threshold of the statistical significance (p-value) used throughout the entire model (default is 0.05)}

\item{num_trees}{Number of trees to be used for the random forests classification (supervised machine learning). Default is set to 10000.}

\item{mtry}{Number of features to possibly split at in each node. Default is the (rounded down) square root of the
number of variables. Alternatively, a single argument function returning an integer, given the number of independent variables.}

\item{num_permutations}{Number of permutations to be used for computation of the statistical significance (p-values) of
the importance scores resulted from random forests classification (default is 100).}

\item{inf_const}{The constant value to be multiplied by the maximum absolute value of differential (logFC)
values for the substitution with infinite differential values. This results in noticeably high biomarker values for features
with infinite differential values compared with other features. Having said that, the user can still use the
biomarker rank to compare all of the features. This parameter is ignored if no infinite value
is present within Diff_data. However, this is used in the case of sc-seq experiments where some genes are uniquely
expressed in a specific cell-type and consequently get infinite differential values. Note that the sign of differential
value is preserved (default is 10^10).}

\item{ncores}{Integer; the number of cores to be used for parallel processing. If ncores == "default" (default), the number of 
cores to be used will be the max(number of available cores) - 1. We recommend leaving ncores argument as is (ncores = "default").}

\item{seed}{The seed to be used for all of the random processes throughout the model (default is 1234).}

\item{verbose}{Logical; whether the accomplishment of different stages of the model should be printed (default is TRUE).}
}
\value{
A list of one graph and one to four tables including:

- Driver table: Top candidate drivers

- DE-mediator table: Top candidate differentially expressed/abundant mediators

- nonDE-mediator table: Top candidate non-differentially expressed/abundant mediators

- Biomarker table: Top candidate biomarkers

The number of returned tables depends on the input data and specified arguments.
}
\description{
This function runs the Experimental data-based Integrated Ranking (ExIR)
model for the classification and ranking of top candidate features. The input
data could come from any type of experiment such as transcriptomics and proteomics.
A shiny app has also been developed for Running the ExIR model, visualization of its results as well as computational
simulation of knockout and/or up-regulation of its top candidate outputs, which is accessible using
the `influential::runShinyApp("ExIR")` command.
You can also access the shiny app online at https://influential.erc.monash.edu/.
}
\examples{
\dontrun{
MyDesired_list <- Desiredlist
MyDiff_data <- Diffdata
Diff_value <- c(1,3,5)
Regr_value <- 7
Sig_value <- c(2,4,6,8)
MyExptl_data <- Exptldata
Condition_colname <- "condition"
My.exir <- exir(Desired_list = MyDesired_list,
               Diff_data = MyDiff_data, Diff_value = Diff_value,
               Regr_value = Regr_value, Sig_value = Sig_value,
               Exptl_data = MyExptl_data, Condition_colname = Condition_colname)
}
}
\seealso{
\code{\link[influential]{exir.vis}},
\code{\link[influential]{diff_data.assembly}},
\code{\link[coop]{pcor}},
\code{\link[stats]{prcomp}},
\code{\link[ranger]{ranger}},
\code{\link[ranger]{importance_pvalues}}

Other integrative ranking functions: 
\code{\link{comp_manipulate}()},
\code{\link{hubness.score}()},
\code{\link{ivi.from.indices}()},
\code{\link{ivi}()},
\code{\link{spreading.score}()}
}
\concept{integrative ranking functions}
\keyword{exir}
