#' Convert seconds to human-readable time units 
#'
#' This function converts a number of seconds created by proc.time() - ptm 
#' to days, hours, minutes, seconds. 
#'    
#' @usage h_read(timing)
#' @param timing the number of elapsed seconds, returned by proc.time() - ptm
#' @export

h_read <- function(timing) {
  timing <- as.numeric(timing)
  days <- floor(timing/86400)
  hours <- floor((timing - 86400*days)/3600)
  minutes <- floor((timing - 86400*days - 3600*hours)/60)
  seconds <- floor(timing - 86400*days - 3600*hours - 60*minutes)
  day.in <- ifelse(hours == 1, "1 day, ",ifelse(days == 0, "", paste(days,"days, ")))
  hour.in <- ifelse(hours == 1, "1 hour, ",ifelse(hours == 0, "", paste(hours,"hours, ")))
  minute.in <- ifelse(minutes == 1, "1 minute, ",ifelse(minutes == 0, "", paste(minutes,"minutes, ")))
  second.in <- ifelse(seconds == 1, "1 second", paste(seconds,"seconds"))
  paste(day.in, hour.in, minute.in, second.in, sep = "")
}
