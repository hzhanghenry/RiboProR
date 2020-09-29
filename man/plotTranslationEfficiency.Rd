\name{plotTranslationEfficiency}
\alias{plotTranslationEfficiency}
\title{
Scatter Plot with Translational Efficiency Data
}
\description{
Scatterplot with too translational efficiency (TE) values.
}
\usage{
plotTranslationEfficiency(dds, ratio_level=1, 
	x_label="TE of Wild Type", y_label="TE of Mutant", 
	title_text="Translation Efficiency Of Mutant", 
	x_pos=0, y_pos=12)
}
\arguments{
  \item{dds}{
  A DESeqDataSet object.
}
  \item{ratio_level}{
  Positive numeric, threshold to change point colors.
}
  \item{x_label}{
  Character vector, text for x axis label.
}
  \item{y_label}{
  Character vector, text for y axis label.
}
  \item{title_text}{
  Character vector, text for title of the plot.
}
  \item{x_pos}{
  Integer, x coordinate for text showing pearson's r
}
  \item{y_pos}{
  Integer, y coordinate for text showing pearson's r
}
}
\examples{
\dontrun{
data("dds_TE")
plotTranslationEfficiency(dds_TE)
}
}
\keyword{ methods }
