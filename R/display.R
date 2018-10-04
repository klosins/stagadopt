# print numbers in pretty way
#' @export
print_pretty_number <- function(x){
  as.character(signif(x,3))
}
