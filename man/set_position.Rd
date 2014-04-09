\name{set_position}
\alias{reset_position}
\alias{set_position}
\title{functions to set position and size for insets}
\usage{
set_position(x, y, size = 0.1, bg = NULL)

reset_position(mfcol = TRUE)
}
\arguments{
  \item{x,y}{user coordinates for inset}

  \item{size}{relative width and height of inset (compared
  to main plot)}

  \item{bg}{colour of background for inset}

  \item{mfcol}{logical, should support for multiple panels
  be added?}
}
\description{
\code{set_position} sets up a new plot at the specified
position (in user coordinates), \code{reset_position}
resets the plot to the original setup.
}
\keyword{utilities}

