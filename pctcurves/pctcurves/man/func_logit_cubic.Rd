\name{func_logit_cubic}
\alias{func_logit_cubic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A logit cubic function in time since infection.
%%  ~~function to do ... ~~
}
\description{
Evaluated within an integral function to estimate roots of the equation at a specific percentile and ODn threshold.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
func_logit_cubic(t, parameters)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
  time (in days) when the function is evaluated up to a maximum of 1000 days
%%     ~~Describe \code{x} here~~
}
\item{parameters}{
  A set of 5 glm model parameters generated for every ODn threshold.
%%     ~~Describe \code{x} here~~
}
}
\details{
This function is evaluated for every time point (in days) with in the uniroot function
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

%## The function is currently defined as
%function (x)
%{
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
