\name{inset}
\alias{inset}
\alias{reset_pos}
\title{functions to set position and size for insets}
\usage{
inset(inset.x, inset.y, fun, inset.size = 0.1, inset.bg = NULL,
  inset.mfcol = TRUE, ...)

reset_pos(oldpar, mfcol = TRUE)
}
\arguments{
  \item{inset.x,inset.y}{user coordinates for inset}

  \item{fun}{Function to be evaluated for inset (e.g.
  \code{plot})}

  \item{inset.size}{relative width and height of inset
  (compared to main plot)}

  \item{inset.bg}{colour of background for inset}

  \item{inset.mfcol}{logical, is original layout with
  multiple panels?}

  \item{...}{additional arguments passed to \code{fun}}

  \item{oldpar}{plotting parameters to which region is
  reset to}

  \item{mfcol}{logical, should support for multiple panels
  be added?}
}
\description{
\code{set_position} sets up a new plot at the specified
position (in user coordinates), \code{reset_position}
resets the plot to the original setup.
}
\keyword{utilities}

