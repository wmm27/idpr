#' Sequence Plot
#'
#' This is a graphical function used to visualize numeric
#'   data along an amino acid sequence.
#' @importFrom ggplot2 aes theme
#' @param position numeric vector of residue positions.
#'   Typically c(1 : sequenceLength). ie a sequence with 215 amino acids
#'   has a vector of values 1 to 215. This is the X axis
#' @param property vector of values, typically numeric. Equal in length
#'   to position. This is the Y axis.
#' @param hline,propertyLimits optional, numeric values or numeric vectors.
#'   Prints horizontal lines. Set to NA to skip (default).
#'    *hline* specifies the location for a dashed,
#'   grey line to be printed underneath the plot's data line.
#'   Good for separating cutoff values. *propertyLimits* specifies the location
#'   for a solid, black line to be printed. Good for showing maximum and
#'   minimum values.
#' @param dynamicColor optional vector. Typically numeric. Equal in length
#'   to position. Can be used to set colors based on values. Can be
#'   categorical (discrete) or continuous. Set to NA to skip (default).
#' @param customColors optional vector of colors as character strings.
#'   Used to support custom plot colors. If property is a discrete scale, a
#'   character vector of colors with length = number of unique discrete
#'   observations is required. If property is a continuous scale, a character
#'   vector of the colors for c("highColor","lowColor","midColor").
#'   Set NA to skip custom colors (default). Ignored if
#'   \code{dynamicColor = NA}.
#' @param midpoint needed for proper scales of customColors. The default
#'   value is equal to hline (if provided). If there is no hline, the average
#'   of propertyLimits is the midpoint (if provided). If neither is provided,
#'   the value will be NA. The user can explicitly assign the midpoint to
#'   avoid this or to overwrite the defaults.
#' @param customTitle  optional, character string. Allows adding custom title.
#'   Set to NA to skip (default).
#' @return a ggplot
#' @export
#' @examples
#' #Get a data frame returned from another function
#' aaVector <- c("A", "C", "D", "E", "F",
#'            "G", "H", "I", "K", "L",
#'            "M", "N", "P", "Q", "R",
#'            "S", "T", "V", "W", "Y")
#' exampleDF <- chargeCalculationGlobal(sequence = aaVector)
#' head(exampleDF)
#' #Making a sequence plot
#' sequencePlot(
#'   position = exampleDF$Position,
#'   property = exampleDF$Charge)
#'
#' #Change the horizontal lines
#' sequencePlot(
#' position = exampleDF$Position,
#' property = exampleDF$Charge,
#' hline = 0.0,
#' propertyLimits = c(-1.0, 1.0))
#'
#' #Adding a dynamic colors based on the property values
#' sequencePlot(
#' position = exampleDF$Position,
#'   property = exampleDF$Charge,
#' hline = 0.0,
#' propertyLimits = c(-1.0, 1.0),
#' dynamicColor = exampleDF$Charge,
#' customColors = c("red", "blue", "grey50"),
#' customTitle = "Charge of Each Residue / Terminus")

sequencePlot <- function(position, property,
    hline = NA,
    propertyLimits = NA,
    dynamicColor = NA,
    customColors = NA,
    midpoint = hline,
    customTitle = NA) {

    dataDF <- data.frame(position = position, property = property)
    gg <- ggplot2::ggplot(dataDF, aes(x = position))

    if (!is.na(hline)) {
        gg <- gg + ggplot2::geom_hline(yintercept = hline,
                                        linetype = "dotdash", color = "gray13",
                                        size = 1, alpha = 0.5)
    }
    if (!is.na(dynamicColor[1])) {
        gg <- gg + ggplot2::geom_line(aes(y = property, color = property),
                                    linetype = "solid")
        if (!is.na(customColors[1])) {
            if (plyr::is.discrete(property)) {
                gg <- gg + ggplot2::scale_color_manual(values = customColors)
            } else {
                if (is.na(midpoint)) {
                    midpoint <- sum(propertyLimits) / length(propertyLimits)
                }
            gg <- gg + ggplot2::scale_color_gradient2(high = customColors[1],
                                                        low = customColors[2],
                                                        mid = customColors[3],
                                                        midpoint = midpoint)
            }
        }
    } else {
        gg <- gg + ggplot2::geom_line(aes(y = property),  linetype = "solid")
    }
    if (!is.na(customTitle)) {
        gg <- gg + ggplot2::labs(title = customTitle,
                                x = "Residue", y = "Score")
    } else {
        gg <- gg + ggplot2::labs(x = "Residue", y = "Score")
    }
    gg <- gg + ggplot2::theme_minimal() +
                ggplot2::theme(legend.position = "none")
    if (!is.na(propertyLimits[1])) {
        gg <- gg + ggplot2::geom_hline(yintercept = propertyLimits,
                                        color = "gray2")
    }
    return(gg)
}
