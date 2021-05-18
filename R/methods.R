print.robvarsel <- function(x, digits = getOption("digits"), ...) 
{
  header1 <- c("Variable selection for Robust Gaussian model-based classification")
  header2 <- "Stepwise greedy-forward approach via TBIC"
  sep <- paste0(rep("-", max(nchar(header1)) + 1), collapse = "")
  cat(sep, "\n")
  cat(header1, sep = "\n")
  cat(header2, sep = "\n")
  cat(sep, "\n\n")
  #
  print(x$steps.info, na.print = "", row.names = FALSE, digits = digits, ...)
  cat("\n")
  footer <- strwrap(paste("Selected subset:", 
                          paste(names(x$subset), collapse = ", ")),
                    width = getOption("width"), simplify = TRUE, 
                    indent = 0, exdent = 2)
  cat(footer, sep = "\n")
  invisible()
}
