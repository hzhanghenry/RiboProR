\name{plotCorrelationMatrix}
\alias{plotCorrelationMatrix}
\title{
Plot Correlation Matrix Data
}
\description{
Plot correlation matrix with corrplot() provided by corrplot package.
}
\usage{
plotCorrelationMatrix(plot_data, cor_method="rcorr",
	cor_type="spearman", shape_type="ellipse", p_threshold=0.01)
}
\arguments{
  \item{plot_data}{
  A numeric matrix, the data with which correlation coefficients will be calculated.
}
  \item{cor_method}{
  Character vector, function name for correlation calculation, either "rcorr" (default) or "cor".
}
  \item{cor_type}{
  Character vector, correlation type, either "spearman" or "pearson".
}
  \item{shape_type}{
  Character vector, shape to represent the correlation coefficients on the plot image, one of "circle", "square", "ellipse", "number", "shade",  "color", or "pie".
}
  \item{p_threshold}{
  Numeric, the threshold for significant level to filter the correlation coefficients.
}
}
\keyword{ methods }

