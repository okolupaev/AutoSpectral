# read_channel.r

#' @title Read Channel Information
#'
#' @description
#' This function reads channel information from control files and corrects
#' channel names based on specified forbidden characters.
#'
#' @importFrom utils read.csv read.table write.table
#' @importFrom flowCore read.FCS exprs
#' @importFrom dplyr filter
#'
#' @param control.dir Directory containing control files.
#' @param control.def.file File containing control definitions.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#'
#' @return A data frame containing the original and corrected channel names.
#'
#' @export

read.channel <- function( control.dir, control.def.file, asp )
{
    # read markers from file if available
    if ( ! is.null( asp$marker.file.name ) && file.exists( asp$marker.file.name ) )
      return(
        read.table(
          asp$marker.file.name, sep = ",",
          stringsAsFactors = FALSE
        )
      )

    # read definition of controls
    control <- read.csv( control.def.file, stringsAsFactors = FALSE )

    # get used channels from controls
    control <- dplyr::filter( control, filename != "" )

    if ( anyDuplicated( control$file.name ) != 0 )
      stop( "duplicated filenames in fcs data", call. = FALSE )

   flow.set.channel <- colnames(
     suppressWarnings(
       flowCore::exprs(
         flowCore::read.FCS(
           file.path( control.dir, control$filename[ 1 ] ),
           truncate_max_range = FALSE,
           emptyValue = FALSE
         )
       )
     )
   )

    # correct channel names
    flow.set.channel.corrected <- flow.set.channel

    for ( fmfc.idx in 1 : nchar( asp$marker.forbidden.char ) )
    {
        fmfc <- substr( asp$marker.forbidden.char, fmfc.idx, fmfc.idx )

        flow.set.channel.corrected <- gsub(
          fmfc, asp$marker.substitution.char,
          flow.set.channel.corrected,
          fixed = TRUE
          )
    }

    # save list of markers
    flow.set.channel.table <- data.frame(
      flow.set.channel,
      flow.set.channel.corrected,
      stringsAsFactors = FALSE
      )

    colnames( flow.set.channel.table ) <- c(
      "flow.set.channel", "flow.set.channel.corrected"
      )

    if ( ! is.null( asp$marker.file.name ) )
        write.table(
          flow.set.channel.table,
          file = asp$marker.file.name,
          row.names = FALSE,
          col.names = FALSE,
          sep = ","
          )

    return( flow.set.channel.table )
}

