###############################################################################################
## Running CharAnalysis in R to Identify Pulses of Tree Recruitment in Seedling or Stand Age Data
## Loading packages
library(zoo)
library(mixtools)
library(ggplot2)
library(sjPlot)

###############################################################################################
## Reading in data - change this to path on your PC
coreData <- read.csv("/Users/kyle/Downloads/KZ_CharAnalysisLC_11June2020.csv")

###############################################################################################
## Char Analysis Parameters
### NOTE: If your seedlings were dated with annual precision (i.e., root shoot boundary),
  # I would use a binwidth of 1. Breaks should be the maximum range of your data plus one year
  # on either end. Currently I have this based on the maximum and minimum observed establishment
  # years buffered by 1. But, if you think that you could have potentially dated younger
  # and older individuals (e.g., samples taken last summer might mean that no dated individuals
  # after 2013 might actually be accurate - those could be "true zeros"), then you might expand the range
  # a bit. "Width" could be set a little smaller or larger. Larger widths cut more potential peak
  # years off of the ends (so a width of 5 would prevent 1976/1977 and 2012/2013 from being considered
  # peaks, while a width of 7 would prevent 1976-1978 and 2011-2013 from being considered peaks). The
  # trade-off is that a wider width also means that peaks in the middle of the series are more
  # likely to be detected. I would keep the cutoff between 0.9 and 0.99, as this is sort of like setting a
  # p-value in a statistical test, and lower values might give you false positives.

width = 7 ## Specifies the window (# of time periods) for the locally weighted regression and gaussian mixture models
cutoff = 0.95 ## Quantile on tail of gaussian model to be defined as the local cutoff value to detect outliers
min_thresh = 2 ## The minimum number of trees in a bin to be considered a peak
multiplier = 1 ## Normalize counts and analysis by plot area, kept at 1 by default which means that analyses are
               ## performed on raw counts rather than seedlings/ha. This doesn't affect results but
               ## affects the appearance of plots.
method = "median" ## Method to de-trend the time series. Either a running median (robust to outliers), or
                  ## locally weighted regression (lowess). Lowess, in it's R implementation, can be prone to
                  ## outlier effects in some cases (i.e. very narrow peaks in establishment surrounded by "0"
                  ## values are not always identified)
iter = 20 ## Because mixtools randomly selects starting values for mixed models, we need to run things a few times to come up with
          ## a similar solution on each run. Choosing a larger number of iterations makes results more consistent,
          ## but takes longer.

###############################################################################################
## Defining functions
gauss = function(x, cutoff){
  ## Helper function for CharAnalysis. Allows for the rolling calculation of gaussian mixture models
  # and extraction of std. dev of model (sigma; centered at 0) for each time step. Currently there appear
  # to be issues with the model crashing when there are too many zeros in a window, so I wrapped
  # it in a try statement to allow things to continue even if it fails on one iteration. This is what leads
  # to the many red error messages, but doesn't really affect the performance
  if(sd(x) > 1/qnorm(p = c(cutoff))){
    t = try(as.numeric(normalmixEM(x, k = 2, mean.constr = c(0,NA))$sigma[1]))
  }else{t = 0}
  return(t)
}

CharAnalysis = function(inputDF, width, cutoff, min_thresh = 2, 
                        pico = NA, potr = NA, method = c("median", "lowess"), iter = 20){
  # Function output is a data frame that includes center years (time) and counts (counts) for each bin, as
  # well as the fitted lowess values (lowess), residuals from the model (residuals), whether or not
  # the bin was identified as a recruitment peak (peaks), thresholds from the gaussian mixture model
  # (thresholds), as well as times with at least 2 pico or potr, and more than twice as many of that
  # species as the previous one. Variation on that used by Tepley and Veblen (2015) Ecol. Monographs
  ## Make sure the packages are there...
  require(zoo) ## For time series calculations and interpolation of missing values
  require(mixtools) ## For use in creating gaussian mixture models and thresholds

  ## Pulling out counts and time
  counts <- inputDF$counts; time <- inputDF$time
  
  ## Method parameter selects lowess or median smoother (see above for more info)
  if(method == "lowess"){
    model = lowess(counts~time, f = width/length(counts), iter = 20)
    pred = model$y
  }
  if(method == "median"){pred = runmed(counts, width)}
  
  ## Setting lower bounds on model for areas with no trees. Predicted values less than zero make no sense
  pred[pred<0] = 0 
  devs = counts - pred
  
  for(j in 1:iter){
    ## Fitting gaussian mixture models to residuals and identifying local cutoff in each window
    # function "gauss" is defined above and assumes two gaussian distributions at each time step,
    # forcing one to go through zero. Running this several times (in this case 20) makes things slower,
    # but makes results consistent between runs.
    gaussians = rollapply(devs, width = width, FUN = gauss, cutoff = cutoff)
    gaussians[grepl("Error", gaussians)] = NA ## Identify the iterations that crashed in the previous line
    
    ## Mixtools can be a little "buggy", occasionally not reaching conversion and producing lots of red 
    # error messages when zeros are present. For now, use a spline to interpolate the few values in which a mixture 
    # model can not be created. Then reassign negative values to zero for plotting and to prevent weirdness.
    gaussians = as.numeric(gaussians)
    if((!is.na(gaussians) & gaussians>0)){
      sub_gaussians = gaussians[min(which(gaussians > 0)):max(which(gaussians > 0))]
      sub_gaussians = na.spline(sub_gaussians)
      gaussians[min(which(gaussians > 0)):max(which(gaussians > 0))] = sub_gaussians}
    gaussians[is.na(gaussians)] = 0
    gaussians[gaussians < 0] = 0
    if(j == 1){gaussians_list = as.data.frame(gaussians)}
    if(j > 1){gaussians_list = as.data.frame(cbind(gaussians_list, gaussians))}
  }
  gaussians = rowMeans(gaussians_list, na.rm = T)
  
  ## identifies number of SD's to cutoff in standard normal distribution and multiplies the sigmas
  # from the gaussian mixture models by this to define local thresholds.
  cutoffs = qnorm(p = c(cutoff)) * gaussians
  buffer = (width-1)/2
  thresholds = c(rep(NA, buffer), cutoffs, rep(NA, buffer))
  means = rollmean(counts, width - 2, na.pad = T)
  for(j in (width - 2):(length(means)-(width - 4))){
    if(means[j] == 0){
      thresholds[j] = 0}
  } ## Remove extraneous values for thresholds in parts of the graph with zeros
  
  ## Buffering the cutoff values with NAs to make vector lengths equal. Identifying peaks in data, then
  # specifying that they must be larger than the previous bin, may not have a bin that is more than twice
  # as big after them (so a pulse is identified at or near the peak, rather than a low value preceding a peak),
  # Successive peaks may also not exist in successive bins, as these are likely the same cohort. Lastly,
  # there must be at least two trees in a given bin to be called a peak.
  peaks = c(ifelse(((devs > thresholds) & (counts >= min_thresh)), 1, 0)) ## Do the bins lie above the local cutoff?
  peaks[is.na(peaks)] = 0 ## Removes NAs for criteria below, in case there are any
  for(j in 1:(length(counts)-1)){if(counts[j+1] > 1.2 * counts[j]){peaks[j] = 0}} ## Do peaks have a much larger bin after them?
  for(j in 2:length(counts)){if(counts[j-1] >= counts[j]){peaks[j] = 0}} ## Are peaks larger than the previous bin?
  for(j in 2:length(peaks)){if(peaks[j-1] == 1){peaks[j] = 0}} ## Are peaks preceded by another one?
  peaks[is.na(peaks)] = 0 ## And a final error check in case anything weird happened
  
  ## Combine it all for export
  peaks = data.frame("time" = time, "estCounts" = counts, "smooth" = pred, "residuals" = devs, "peaks" = peaks, thresholds)
  return(peaks)
}

plot_pulses = function(peaks, title, dataRange){
  ## Plotting for display. Works fine as is, but can be tweaked as desired.
  # First, creating spline models to make smoother looking curves through same points
  data = peaks[complete.cases(peaks$thresholds),]
  splines = smooth.spline(x = data$time, y = data$smooth, df = nrow(data)-1)
  splines = as.data.frame(predict(splines, seq(min(data$time), max(data$time))))
  
  splines2 = smooth.spline(x = data$time, y = data$thresholds, df = nrow(data)-1)
  splines2 = as.data.frame(predict(splines2, seq(min(data$time), max(data$time))))
  splines2$y = splines$y+splines2$y
  
  ## And creating initial plot
  peaks_plot = ggplot(peaks, aes(x = time, y = counts)) + 
    ggtitle(title) + theme_bw() + theme(axis.title.x=element_blank()) +
    ylab("count") + xlim(min(dataRange), max(dataRange)) +
    ylim(0, max((max(peaks$counts, na.rm = T) + 3*multiplier), (max(splines2$y) + 3 *multiplier))) + 
    geom_line(data = splines2, aes(x = x, y = y), lty = 2, colour = "black") +
    geom_line(data = splines, aes(x = x, y = y), colour = "black") + 
    geom_bar(col = "black", fill = "deepskyblue4", stat= "identity", alpha = 0.9) + 
    theme(plot.title = element_text(hjust = 0.5)) 
  
  ## Adding identifiers for peak locations - black for all trees, green for aspen, orange for lodgepole
  for(j in 1:length(peaks$peaks)){
    if(peaks$peaks[j] == 1){
      peaks_plot = peaks_plot + annotate("point", x = peaks[j,1], y = peaks[j,2] + multiplier, shape = "\u25BC", size = 3)
    }
  }
  
  ## And displaying graph...
  return(peaks_plot)
}

##################################################################################################################
##################################################################################################################
##################################################################################################################
## To run with the code above, the first function argument is a data frame with the
  # columns "count" representing the number or density of establishing trees in a
  # given time period, and "time" representing the time period (middle of the period
  # probably makes the most sense)

## Making up some example data with sine curve and uniform error
inputDF <- data.frame(time = seq(1800, 2000, by = 10),
                      counts = c(25, 20, 10, 10, 50, 45, 10, 37, 12, 11, 65,
                                 12, 44, 34, 54, 44, 20, 15, 10, 20, 21))

## Running it  
peaks <- CharAnalysis(inputDF, width, cutoff = cutoff, method = method, 
                      min_thresh = min_thresh, iter = iter)
plot_pulses(peaks, "Recruitment Pulses - All Plots", dataRange = c(1800, 2000))

