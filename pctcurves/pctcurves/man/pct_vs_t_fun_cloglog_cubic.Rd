\name{pct_vs_t_fun_cloglog_cubic}
\alias{pct_vs_t_fun_cloglog_cubic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Percentiles generated from P_r(t), with glm and clolog link gml funtion
%%  ~~function to do ... ~~
}
\description{
A function to create percentile curves using a glm function and cloglog link function, with time since infection as a cubic funtion
}
\usage{
func_logit_cubic(data_set, ODnTh, t_since_inf, percentile)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_set}{
  The dataset to be evaluated. It should include the following variables eddi (for time since infection), viral_load (for viral load) and ODn (recency test reuslt)
%%     ~~Describe \code{x} here~~
}
\item{ODnTh}{
  The ODn resolution you need to evaluate the probability of being recent at time t (P_r(t))
%%     ~~Describe \code{x} here~~
}
\item{t_since_inf}{
  Timepoints within the intertest interval at which the function should be evaluated. in this case we use daily time steps from 1 to 1000 for this function. These timsteps are prefered because our aim with this is to create a dataset of descrete timesteps for use in the likelihood function.
%%     ~~Describe \code{x} here~~
}
\item{percentile}{
  Sequence of percentils from 0.1 to .9 in steps of 0.1
%%     ~~Describe \code{x} here~~
}
}
\details{
function contains uniroot, a function that estimates roots (time points) at where the percentile is true for the func_logit_squared
%%  ~~ If necessary, more details than the description above ~~
}
\value{
funtcion returnes a dataframe the includes timepoints when the percentile if true for the func_logit_squared.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Joseph B. Sempa
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{uniroot}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
ODnth <- seq(0.01, 5, 0.01)#c(.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)#
time_var <- seq(1, 1000, 1)
x <- data.frame(pct_vs_t = pct_vs_t_fun_cloglog_cubic(
  data_set = data_generate_pr_t,
  ODnTh = ODnth,
  t_since_inf = time_var,
  percentile = seq(0.1, 0.9, 0.1)
))
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
