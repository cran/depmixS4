
formatperc <-
function (x, digits) {
    paste(format(100 * x, trim = TRUE, 
                 scientific = FALSE, digits = digits), "%", sep="")
}
