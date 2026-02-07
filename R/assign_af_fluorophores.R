# assign_af_fluorophores.r

#' @title Assign AF Spectrum By Fluorophore Projection
#'
#' @description
#' Projects the autofluorescence spectral variants into fluorophore unmixed
#' space to determine which best fits each cell (event). A fast approximation for
#' brute force sequential unmixing method in early versions of AutoSpectral.
#' Provides essentially identical results to minimization of fluorophore signal
#' (dist0 method), with minor differences being primarily due to the use of L2
#' (squared) error, rather than L1 (absolute) error. Substantially faster.
#'
#' @param raw.data Expression data from raw FCS files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns.
#'
#' @return Row indices for best-fitting AF spectra (from `af.spectra`)
#'
#' @export

assign.af.fluorophores <- function(
    raw.data,
    spectra,
    af.spectra
) {
  # calculate pseudoinverse
  S <- t( spectra )
  P <- solve( t( S ) %*% S ) %*% t( S )

  # how much each AF variant 'looks like' each fluorophore
  v.library <- P %*% t( af.spectra )
  # squared norms
  v.norms.sq <- colSums( v.library^2 )

  # calculate the 'residual AF' (the part of AF fluorophores can't explain)
  r.library <- t( af.spectra ) - ( S %*% v.library )

  # predicted AF intensity
  numerator <- raw.data %*% r.library

  # denominator (vector of length af.n)
  denominator <- colSums( r.library^2 )

  # k_matrix (cell.n x af.n): estimated AF intensity per cell/variant
  k.matrix <- sweep( numerator, 2, denominator, "/" )

  # initial unmix (no AF)
  unmixed <- raw.data %*% t( P )

  # minimize L2 error via dot product of unmixed signal x AFs
  dot.product <- unmixed %*% v.library

  # score the AF variants
  scores <- ( 2 * k.matrix * dot.product ) - sweep( k.matrix^2, 2, v.norms.sq, "*" )

  # maximum score corresponds to the minimum L2 error
  return( max.col( scores ) )

}
