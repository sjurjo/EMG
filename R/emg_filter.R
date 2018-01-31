install.packages("biosignalEMG")
library(biosignalEMG)
library(readxl)
library(tidyverse)

# change graphical parameters to show multiple plots
op <- par(mfrow = c(2, 2))

# Simulate 5 seconds of an EMG
emgx <- syntheticemg(n.length.out = 5000, on.sd = 1, on.duration.mean = 350,
                     on.duration.sd = 10, off.sd = 0.05, off.duration.mean = 300, off.duration.sd = 20,
                     on.mode.pos = 0.75, shape.factor = 0.5, samplingrate = 1000, units = "mV",
                     data.name = "Synthetic EMG")
plot(emgx, main = "Synthetic EMG")
str(emgx)

#Full-wave rectified EMG
emgr <- rectification(emgx, rtype = "fullwave")
plot(emgr, main = "Rectified EMG")

# Integration of the full-wave rectified
#Integration of the full-wave rectified EMG with reset points every
# 200 samples
emgi <- integration(emgr, reset = TRUE, reset.criteria = "samples", vreset = 200)
plot(emgi, main = "Integrated EMG")

#MA-envelope
emgma <- envelope(emgx, method = "RMS", wsize = 60)
plot(emgma)
# Ensemble-averaged EMG
ea <- eaemg(emgma, runs = emgx$on.off, what = 1, timenormalization = "mean",
            scalem = 1, empirical = TRUE, level = 0.9)
plot(ea, lwd = 2, main = "Ensemble-averaged EMG")

# reset graphical parameters
par(op)










###EMG OBJECT, Cycling analysis 2018

emg <- read_excel("data/VLLtest3_20cyc.xlsx") %>% 
  data.frame()%>%
  print()


str(emg)

##Change to EMG signal (function; as.emg)

emg.x <- emg$VLL
plot(emg.x, channels = "all", samples = 0, type = "l",  timeunits = "seconds",
     add = FALSE)

length(emg.x)

emgraw <- as.emg(emg.x, samplingrate = 1500, units = "mV", date.name = "VLL") # coerce to an EMG object # sampling rate 1500Hz to avoid signal loss
is.emg(emgraw) # EMG object check                                             # sampling rate 1500Hz to avoid signal loss
plot(emgraw)

# High pass filter
x <- highpass(emgraw, cutoff = 20) # 20 Hz recommended for sport activities (maybe as much as 30 Hz, avoid artifacts of the signal) Deluca, 2010.
emghpass <- rectification(x)

plot(x, main = "High pass filter")
plot(emghpass, main = "High pass filter rectification")
emghpass #data high pass filter and rectified

## Full-wave rectified EMG
emgrect <- rectification(x, rtype = "fullwave", data.name = "VLL") #fullwave returning all values to positive values
plot(emgrect, main = "Rectified EMG")
emgrect


### Active or passive phase

emgthr <- onoff_singlethres(emghpass, t = 90, data.name = "RMS") #correct

plot(emghpass, main = "Sample EMG")
plot(emgthr, type = "l", main = "Detected phases (single thresholding)")

### EMG RMS (preferred recommendation for smoothing)

op <- par(mfrow = c(2, 1)) # change graphical parameters to show multiple plots

emgrms <- envelope(emghpass, method = "MA", wsize = 100) # wsize; 100 ms, most common for all conditions
emgrms
plot(emgrms, main = "RMS-envelope")

mean(emgrms$values, na.rm=T)
max(emgrms$values, na.rm=T)

par(op) # reset graphical parameters 


# Ensemble-averaged EMG
avgemg <- eaemg(emgrect, runs = emgthr, what = 1, timenormalization = "mean", scalem = 1, empirical = TRUE, level = 0.95)
avgemg
plot(avgemg, lwd = 1.5, main = "Ensemble-averaged EMG")

mean(avgemg$intervals[,2], na.rm=T)
max(avgemg$intervals[,2], na.rm=T)


plot(emgthr, emgrect$values)




# onoffbonato
sigma_n <- sd(tail(emgraw$values, 15))
sigma_n
b <- onoff_bonato(emgraw, sigma_n = sigma_n, Pfa = 0.05, m= 100, minL = 0, data.name = "Threshold")
plot(b, type = "b")



