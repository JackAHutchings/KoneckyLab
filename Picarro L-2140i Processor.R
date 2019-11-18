#Konecky Lab - Picarro Data Processing Version 1.0
rm(list=ls())
{
  ### Uncomment the two lines below to install the required packages for this script. They only need to be installed once.
  # install.packages(c("dplyr","tidyr","ggplot2","readxl","writexl","zoo","lubridate","rstudioapi","BiocManager","bit64"))
  # BiocManager::install("rhdf5") # Required to deal with Picarro's data format  
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  theme_set(theme_classic())
  require(readxl)
  require(writexl)
  require(zoo)
  require(lubridate)
  require(rhdf5)
  require(rstudioapi)
}

# STEP 1: Save this script in the same directory as the Private Data .zip files and the matching Piccaro_Tray_Sheet.
# STEP 2: Provide your name or initials! This will be attached to various output files.
vars$username = "JAH"
# STEP 3: Execute the script in order. User-defined options are confined to the 'vars' object, which is in-line with the code
#         near where the user-defined option is relevant. In general, you should not need to change these, but do pay attention
#         just in case!

# Sets the working directory to the script location
directory = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(directory)

#Extract  all the data from each Picarro h5 file and combine into the 'rawdata' object
{
  sapply(list.files(pattern=".zip"),FUN=unzip,junkpaths=T) # Unzips the h5 files
  files <- list.files(pattern=".h5") # Generate a list of the h5 files
  rawdata <- lapply(files,h5read,name="results",bit64conversion='bit64') # Reads in each file as a data.frame in this list.
  rawdata <- do.call("rbind",rawdata) # collapse the list of data.frames into a single data.frame.
  unlink(files) # Deletes the unzipped h5 files.
  # write.csv(rawdata,"rawdata.csv",row.names=F) # If you want the rawdata as a csv for some reason, uncomment this line. Be warned, the raw CSV will be ~600 mb.
}

# New Data Frame with new variables
# All isotope delta values (and derivative metrics) are based on on the apparent R of a VSMOW2 series of injections made previously.
# These need to be properly scale normalized to SMOW-SLAP using the standards from this run.

data <- rawdata %>% 
  mutate(ref.rD = 0.148605440565915,   #      Observed R of VSMOW2 for the Konecky Lab L2140-i from 2018-12-17
         ref.r17O = 0.589386241299641, #      Observed R of VSMOW2 for the Konecky Lab L2140-i from 2018-12-17
         ref.r18O = 0.966147395314548) %>%  # Observed R of VSMOW2 for the Konecky Lab L2140-i from 2018-12-17
  mutate(rDH = str3_offset / str2_offset, # 2H / 1H
         r1716 = str13_offset / str2_offset, # 17O / 16O
         r1816 = str11_offset / str2_offset, # 18O / 16O
         r1816_alt = str1_offset / str2_offset, #18O / 16O, using the other 18O peak
         # dD_picarro = rDH * 6860.3159 -1028.4299, # dD using the Picarro calibration
         # d17O_picarro = r1716 * 1712.8345 - 1009.0035, # d17O using the Picarro calibration
         # d18O_picarro = r1816 * 1030.6447 - 996.3966, # d18O using the Picarro calibration
         # dxs_picarro = dD_picarro-8*d18O_picarro, # dueterium excess using Picarro calibration,
         dD = (rDH/ref.rD-1)*1000,
         d17O = (r1716/ref.r17O-1)*1000,
         d18O = (r1816/ref.r18O-1)*1000,
         d17O_prime = log(r1716/ref.r17O),
         d18O_prime = log(r1816/ref.r18O),
         D17O = (d17O_prime - 0.528*d18O_prime)*1000,
         dxs = dD-8*d18O,
         unixhours = time/60/60,
         datetime = as.POSIXct(time,origin="1970-01-01",tz='CST6CDT'))

vars <- data.frame(h2o_threshold_lower = 7500, #ppm of H2O to reach before looking for peaks
                   h2o_threshold_upper = 27500, #ppm of H2O to reach before stopping looking for peaks; this helps remove end-of-pulse jumps
                   h2o_smoothing_n = 10, #number of consecutive scans to calculate a rolling mean for
                   h2o_smoothed_diff_threshold = 125, #differential ppm of H2O to ID the front end of a peak
                   h2o_slope_maximum = 2500, # this is an attempt to remove end-of-injection 'spikes' from the peak set
                   h2o_slope_minimum = 75, # also is an attempt to remove end-of-injection 'spikes' from the peak set
                   peak_start_trim = 30, #number of scans at the start of a peak to drop
                   peak_end_trim = 20, #number of scans at the end of a peak to drop
                   short_integration = 180, #count of scans (aka seconds) to use for short integrations
                   peak_interval_sd_threshold = 5, #How many standard deviations from the mean peak-to-peak interval to use to flag false positives?
                   residual_threshold = 2, #Flag samples with values greater than this many standard deviations away from standard injections.
                   baseline_shift_threshold = 18, #Flag samples with values greater than this many standard deviations away from standard injections.
                   baseline_curvature_threshold = 3, #Flag samples with values greater than this many standard deviations away from standard injections.
                   dxs_threshold = 0 #Flag samples with deuterium excess values below this value, typically indicative of vial evaporation.
)

data_pulses <- data %>% filter(ValveMask>1) %>% #filter the data for when the vaporizer valve is set introducing an injection
  mutate(scan_num = 1:n()) #create an index of scans

data_pulse_peaks <- data_pulses %>% filter(H2O > vars$h2o_threshold_lower & H2O < vars$h2o_threshold_upper) %>% #remove scans outside of our H2O thresholds
  select(datetime,H2O,scan_num) %>% 
  mutate(raw_diff = abs(c(NA,diff(H2O)))) %>% filter(raw_diff < vars$h2o_slope_maximum) %>% 
  mutate(smoothed_h2o = rollmean(H2O,vars$h2o_smoothing_n,fill=NA),
         smooth_diff = c(NA,diff(smoothed_h2o))) %>% #create a rolling mean of H2O to smooth out any bumps
  mutate(start_peak = ifelse(smooth_diff > vars$h2o_smoothed_diff_threshold,T,F)) %>%  #identify the rising section of each pulse, with a threshold value
  mutate(start_peak = ifelse(abs(raw_diff) < vars$h2o_slope_minimum,F,start_peak), #Apply a minimum slope threshold
         start_peak = ifelse(scan_num<min(scan_num[which(start_peak==T)]),T,start_peak), #cleaning the start of the data set
         start_peak = ifelse(scan_num>max(scan_num[which(start_peak==F)]),F,start_peak)) #cleaning the end of the data set

data_pulse_peaks_2 <- data_pulse_peaks %>% filter(start_peak==T) %>% 
  mutate(start_peak_index = c(2,diff(scan_num))) %>% #uses the scan index and the diff function to find the start of each peak
  filter(start_peak_index > 1) %>% #filter for just the start of each peak
  mutate(start_peak_index = 1:n()) %>% #create a peak index - this should match the visible number of pulses in the dataset
  select(datetime,start_peak_index) %>% 
  group_by(start_peak_index) %>% 
  mutate(interval = ifelse(start_peak_index == min(.$start_peak_index),NA,as.numeric(datetime-.$datetime[which(.$start_peak_index == start_peak_index-1)]))) %>% #Time between this peak and the previous
  ungroup() %>% 
  mutate(mean_interval = mean(interval,na.rm=T), #Average peak-to-peak interval
         sd_interval = sd(interval,na.rm=T), #SD of peak-to-peak interval
         flag = ifelse(interval < mean_interval-(vars$peak_interval_sd_threshold*sd_interval) | 
                         interval > mean_interval+(vars$peak_interval_sd_threshold*sd_interval),
                       "Bad","Good"), #Flags peaks with >2 SD peak-to-peak intervals
         flag = ifelse(start_peak_index == 1,"Good",flag)) #The first peak often has strange timing, so it is always flagged 'Good' 

peak_flag_override = c() #Enter the 'start_peak_index' from data_pulse_peaks_2 of any peak flagged as Bad that is actually good.
# This should only be necessary when there has been instrument failure, but other things (Daylight Savings) can sometimes be a problem.

data_pulse_peaks_3 <- data_pulse_peaks_2 %>%
  mutate(flag = ifelse(start_peak_index %in% peak_flag_override,"Good",flag)) %>% 
  filter(flag == "Good" | is.na(flag)) %>%   #Select only the first peak (NA) and 'Good' peaks using the 2 SD criterion
  ungroup() %>% 
  mutate(start_peak_index = 1:n(), # Recreate the start peak index after having removed false positives
         mean_minutes = floor(mean_interval-1), # Subtract one minute to account for the variability in peak-to-peak duration
         mean_seconds = (mean_interval - (mean_minutes+1))*60,
         end_peak_time = datetime + minutes(mean_minutes) + seconds(mean_seconds)) %>%  # Ending peak time based on average peak-to-peak duration
  select(datetime,end_peak_time,interval,start_peak_index)




data_pulse_peaks_4 <- full_join(data_pulses,data_pulse_peaks_3,by = "datetime") %>% #rejoin with the pulses data
  select(datetime,H2O,spect_duration,end_peak_time,start_peak_index) %>% 
  mutate(peak_number = na.locf(start_peak_index,na.rm=F)) %>% #use the na.locf function to apply each peak starting point to its full peak
  group_by(peak_number) %>% 
  filter(is.na(peak_number)==F) %>% #remove non-peaks
  filter(datetime < end_peak_time[which(!is.na(end_peak_time))]) %>% #trim peaks based on the end_peak_time
  filter(H2O > vars$h2o_threshold_lower) %>% #re-apply our H2O threshold filter
  mutate(peak_scan_num = 1:n()) %>% #create an index of scans within each peak
  filter(peak_scan_num > vars$peak_start_trim & peak_scan_num < max(peak_scan_num)-vars$peak_end_trim) %>% #trim the front/back of each peak
  mutate(usable_peak_length = sum(spect_duration),#calculate the 'usable' data of each peak. For default injection settings, this should be around 200 s
         peak_seconds = cumsum(spect_duration), #re-create the peak scan index with the trimmed peaks
         short_integration = ifelse(peak_seconds <= vars$short_integration,T,F), #flag the short integration time, typically for D/H
         is_peak = T) %>% #is there a peak at this scan?
  full_join(data,by = c("datetime", "H2O","spect_duration")) %>% #re-join with the full data set
  ungroup() %>% arrange(datetime) %>% 
  group_by(peak_number) %>% 
  mutate(mean_datetime = mean(datetime[which(is_peak)])) #find the middle of each peak in time, for plotting

### Important!!!
# The following is a plot of all identified pulses and their assignment.
# The starting_pulse variable should correspond to the first pulse in your run.
# If peaks are missing or there are extra peaks, then my peak-finding parameters are off...
# To resolve this, you'll want to adjust the first six h2o parameters in the 'vars' object above
# to increase or decrease the sentivitity to finding peak tops, depending on your problem.
# You'll then need to re-run the code starting from line 33 (where the 'data' object is created)

# To troubleshoot, the chunk of code after this plot allows you to plot specific peaks, which can help find problems...

vars <- vars %>% 
  mutate(plot_start_adjust = 1, #Number of hours before the first detected injection to plot. Set this to -Inf to plot from the beginning.
         plot_end_adjust = 1, #Number of hours after the last detected injection to plot. Set this to Inf to plot until the end.
         resolution_factor = 8 #Reduces the resolution of the plot below as a multiple of this factor. Makes plotting faster.
  )

ggplot(data_pulse_peaks_4 %>% ungroup() %>% select(peak_number,is_peak,mean_datetime,datetime,H2O,d18O) %>% gather(var,val,H2O,d18O) %>% 
         filter(datetime > (min(mean_datetime,na.rm=T)-(3600*vars$plot_start_adjust)) & 
                  datetime < (max(mean_datetime,na.rm=T)+(3600*vars$plot_end_adjust))) %>% 
         mutate(remove_point = rep(seq(from = 1, to = vars$resolution_factor, by=1),length.out=n())) %>% 
         filter(remove_point == vars$resolution_factor),
       aes(x=datetime,y=val)) +
  geom_path(aes(color=is_peak,group=1)) +
  geom_text(data=data_pulse_peaks_4  %>% select(mean_datetime,is_peak,peak_number,H2O,d18O) %>%
              gather(var,val,H2O,d18O) %>% group_by(var) %>%
              mutate(val = max(val[which(is_peak==T)])+abs(max(val[which(is_peak==T)]))*0.2) %>% select(is_peak,val,var,mean_datetime,peak_number) %>%
              distinct() %>% na.omit() %>% 
              filter(mean_datetime > (min(mean_datetime,na.rm=T)-(3600*vars$plot_start_adjust)) & 
                       mean_datetime < (max(mean_datetime,na.rm=T)+(3600*vars$plot_end_adjust))),
            aes(x=mean_datetime,y=val,label=peak_number),size=3,angle=90) +
  facet_wrap(~var,ncol=1,scales="free_y") +
  scale_x_datetime(breaks = "4 hour", date_labels = "%b %d\n%I %p") +
  theme(axis.text.x = element_text(angle=0))

# Uncomment the following code and enter a peak number in vars$peak_to_plot. This will plot an hour of measurements centered on that peak.
# Useful for troubleshooting, but not necessary under normal conditions.
{
  # vars$peak_to_plot = 10
  # 
  # ggplot(data_pulse_peaks_4 %>% ungroup() %>% select(peak_number,is_peak,mean_datetime,datetime,H2O,d18O) %>% gather(var,val,H2O,d18O) %>% 
  #          filter(datetime > (mean_datetime[which(peak_number == vars$peak_to_plot)][1] - (30*60)) & 
  #                 datetime < (mean_datetime[which(peak_number == vars$peak_to_plot)][1] + (30*60))) %>% 
  #          mutate(remove_point = rep(seq(from = 1, to = vars$resolution_factor, by=1),length.out=n())) %>% 
  #          filter(remove_point == vars$resolution_factor),
  #        aes(x=datetime,y=val)) +
  #   geom_path(aes(color=is_peak,group=1)) +
  #   geom_text(data=data_pulse_peaks_4  %>% select(datetime,mean_datetime,is_peak,peak_number,H2O,d18O) %>%
  #               gather(var,val,H2O,d18O) %>% group_by(var) %>%
  #               filter(datetime > (mean_datetime[which(peak_number == vars$peak_to_plot)][1] - (30*60)) & 
  #                        datetime < (mean_datetime[which(peak_number == vars$peak_to_plot)][1] + (30*60))) %>% 
  #               mutate(val = max(val[which(is_peak==T)])+abs(max(val[which(is_peak==T)]))*0.2) %>% select(is_peak,val,var,mean_datetime,peak_number) %>%
  #               distinct() %>% na.omit() %>% 
  #               filter(mean_datetime > (min(mean_datetime,na.rm=T)-(3600*vars$plot_start_adjust)) & 
  #                        mean_datetime < (max(mean_datetime,na.rm=T)+(3600*vars$plot_end_adjust))),
  #             aes(x=mean_datetime,y=val,label=peak_number),size=3,angle=90) +
  #   facet_wrap(~var,ncol=1,scales="free_y") +
  #   scale_x_datetime(breaks = "4 hour", date_labels = "%b %d\n%I %p") +
  #   theme(axis.text.x = element_text(angle=0))
}

vars$start_pulse = 1 #Identify the first pulse in your run here!


# Read in and re-format your tray data
tray_info <- read_excel("Picarro_Tray_Sheet.xlsx",sheet=1)[,1:4] %>% na.omit() %>% 
  arrange(position) %>% 
  group_by(position,sample) %>% 
  mutate(injection = list(c(1:injections))) %>% 
  unnest() %>% 
  ungroup() %>% 
  mutate(peak_number = (1:n())-1+vars$start_pulse)

# Join tray info and pulses by peak_number
# This output file should be considered 'raw'
pulse_level_results <- left_join(tray_info,data_pulse_peaks_4,by = "peak_number") %>% 
  group_by(position,sample) %>% 
  mutate(vial_seconds = cumsum(spect_duration))  # Calculate an index of consecutive scans from a vial through all of its injections

# Date the run was started
rundate = as.Date(min(pulse_level_results$datetime),format="YYYY-MM-DD")

# Summarize data by injection.
# Term    Explanation
# short   First 180 seconds of usable peak data
# mean    Arithmetic mean of the peak data
# sd      Standard deviation of the peak data
# se      Standard error of the peak data
# slope   OLS-Slope during the peak: variable ~ peak_seconds

injection_level_results_raw <- pulse_level_results %>%
  group_by(position,sample,sample_type,injection) %>% 
  summarize(n = n(), # Count of data points in the peak top
            n_short = length(short_integration[which(short_integration==T)]), # Count of data points for the 'short' integration
            datetime_start = min(datetime), # Calendar date and time at the start of the peak, the start of an injection is about 3 minutes prior to this
            duration = max(peak_seconds), # Duration of the pulse peak
            duration_short = max(peak_seconds[which(short_integration==T)]), # Duration of the short integration
            h2o_mean = mean(H2O), # H2O, ppm
            h2o_sd = sd(H2O), # H2O, ppm
            h2o_slope = lm(H2O~peak_seconds)$coefficients[2], # H2O, ppm/s
            h2o_mean_short = mean(H2O[which(short_integration==T)]), # H2O, ppm
            h2o_sd_short = sd(H2O[which(short_integration==T)]), # H2O, ppm
            cavity_torr_mean = mean(CavityPressure), # Cavity pressure, torr
            cavity_torr_sd = sd(CavityPressure), # Cavity pressure, torr
            cavity_temp_mean = mean(CavityTemp), # Cavity temperature, torr
            cavity_temp_sd = sd(CavityTemp), # Cavity temperature, torr
            cavity_torr_mean_short = mean(CavityPressure[which(short_integration==T)]), # Cavity pressure, torr
            cavity_torr_sd_short = sd(CavityPressure[which(short_integration==T)]), # Cavity pressure, torr
            cavity_temp_mean_short = mean(CavityTemp[which(short_integration==T)]), # Cavity temperature, torr
            cavity_temp_sd_short = sd(CavityTemp[which(short_integration==T)]), # Cavity temperature, torr
            dD_mean_short = mean(dD[which(short_integration==T)]), # dD relative to VSMOW2, STILL REQUIRES CALIBRATION
            dD_sd_short = sd(dD[which(short_integration==T)]), # dD relative to VSMOW2, STILL REQUIRES CALIBRATION
            d17O_mean = mean(d17O), # d17O relative to VSMOW2, STILL REQUIRES CALIBRATION
            d17O_sd = sd(d17O), # d17O relative to VSMOW2, STILL REQUIRES CALIBRATION
            d18O_mean = mean(d18O), # d18O relative to VSMOW2, STILL REQUIRES CALIBRATION
            d18O_sd = sd(d18O), # d18O relative to VSMOW2, STILL REQUIRES CALIBRATION
            d18O_slope = lm(d18O~peak_seconds)$coefficients[2], # d18O change during pulse
            rawD17O_mean = mean(D17O), # D17O based on uncalibrated data
            dxs_mean_short = mean(dxs[which(short_integration==T)]), # deuterium excess for potential problem flagging
            residual_mean = mean(residuals), # RMS-residual of the spectra fit, used for organics flagging
            baselineshift_mean = mean(baseline_shift), # Spectral baseline shift, used for organics flagging
            baselinecurvative_mean = mean(baseline_curvature), # Spectral baseline curvature based on a quadratic fit, used for organics flagging
            slopeshift_mean = mean(slope_shift) # Slope of the spectral baseline, used for organics flagging (unclear on proper threshold currently)
  ) %>% 
  ungroup() %>% mutate(run_n = seq(1:n()))

# This calculates 'fmem' or the contribution of the current injection to the observed signal. 1-fmem is equal to the contribution of the previous vial.
fmem_terms <- injection_level_results_raw %>% group_by(position,sample) %>% filter(n() > 9) %>% # Find only the 12-injection standards
  filter(position!=1) %>% 
  select(position,sample,sample_type,run_n,injection,dD_mean_short,d17O_mean,d18O_mean) %>%
  gather(variable,measured,dD_mean_short:d18O_mean) %>% 
  group_by(position,sample,variable) %>% 
  mutate(memory_run = ifelse(sample_type=="memory",T,F),
         bar = ifelse(memory_run,
                      mean(measured[which(injection>(max(injection)-7))]), # For memory runs, use the last 8 injections as the 'true' value.
                      mean(measured[which(injection>(max(injection)-3))]))) %>% # For normal runs, use the last 3 injections as the 'true' value.
  mutate(prev_bar = ifelse(position==1,NA,.$bar[which(.$variable == variable & .$position == position-1)])) %>%  # Determine the previous vial's bar
  na.omit() %>% # Removes the first standard from each series as they do not have previous known 'true' values
  mutate(fmem = 1-(bar - measured)/(bar-prev_bar)) %>% # Calculates the contribution of the current injection to the overall signal
  filter(fmem < 1) %>% # fmem is a fraction that can only vary between 0 and 1; values above 1 typically come from the system being memory-free and the numerator in fmem is recording noise
  mutate(memterm_limit = ifelse(unique(sample_type)=="memory",17,9)) %>% 
  filter(injection <= memterm_limit) %>% # We only calculate fmem for the injections assumed to carry memory (we assume the last 3 are memory-less)
  group_by(memory_run,variable,injection) %>% 
  summarize(fmem_mean = mean(fmem), # Mean fmem for all standard transitions (n=3 typically)
            fmem_sd = sd(fmem)) # SD of fmem, these should be small.

# ggplot(injection_level_results_raw %>% filter(position > 1 & position < 5),aes(x=datetime_start,y=dD_mean_short)) +
#   geom_point() + geom_errorbar(aes(ymin = dD_mean_short - dD_sd_short, ymax = dD_mean_short + dD_sd_short))

# Write the current memory terms and their run-type, and then use the latest memory-run terms for correcting the data.
# If this is a memory coefficient run (i.e., all samples are of type "memory"), then the script can stop after running this code.
{
  setwd("C:/Box Sync/Konecky Lab/Data/Picarro L2140/Memory Coefficients")
  fmem_history <- read_excel("memory_coefficients.xlsx")
  fmem_data <- rbind(fmem_history,fmem_terms %>% mutate(rundate = rundate) %>% select(rundate,memory_run:fmem_sd) %>% as.data.frame(.)) %>% 
    distinct()
  fmem_terms <- fmem_data %>% filter(memory_run) %>% filter(rundate == max(rundate))
  write_xlsx(fmem_data,"memory_coefficients.xlsx")
  rm(fmem_history,fmem_data)
  setwd(directory)
  }

# This sets up a data frame containing all the previous 'bars' (a weighted mean of last 3 injections) for every vial.
# These 'prev_bar' value are dissimilar from above. We do not believe they are the true values of the vials, but, instead 
# characterize the state of the instrument's memory. This takes the fmem with injection=1 assigns a decreasing influence
# to the last 3 injections of the previous vial. First, we calculate fmem_m1, or (1-fmem). This is now the contribution of
# the preceding injection. The last injection of the previous vial gets the full weight of this. The second-to-last injection of the
# previous vial has a weight by a factor of fmem_m1*fmem_m1. The third-to-last is weighted fmem_m1*fmem_m1*fmem_m1. In practice,
# this weights the 'prev_bar' to the isotopic composition of the last injection while still accounting for the previous two.
# This approach implicitly accounts for the changing memory during the last 3 injections of the previous vial as it is still
# experiencing memory from its preceding vial (position-2). Based on Groning et al., 2011.

memory_prev_bars <- injection_level_results_raw %>% group_by(position,sample) %>% 
  select(position,sample,run_n,injection,dD_mean_short,d17O_mean,d18O_mean) %>% 
  gather(variable,value,dD_mean_short:d18O_mean) %>% 
  group_by(position,sample,variable) %>% filter(injection>max(injection)-3) %>% 
  full_join(fmem_terms %>% filter(injection==1) %>% select(variable,fmem_mean),by = "variable") %>% 
  mutate(fmem_m1 = (1-fmem_mean)^rev(1:n()), 
         weight = fmem_m1/sum(fmem_m1)) %>% 
  summarize(bar = sum(weight*value)) %>% 
  mutate(prev_bar = ifelse(position==1,NA,.$bar[which(.$position == position-1)])) %>% 
  select(-bar)

injection_level_results_memory_corrected_raw <- injection_level_results_raw %>%
  select(position,sample,sample_type,run_n,injection,dD_mean_short,d17O_mean,d18O_mean) %>%
  gather(variable,value,dD_mean_short:d18O_mean) %>%
  full_join(fmem_terms,by=c("variable","injection")) %>%
  full_join(memory_prev_bars,by = c("position", "sample", "variable")) %>%
  group_by(position,sample,variable) %>%
  mutate(memcorr = ifelse(is.na(prev_bar),value,(value - (1-fmem_mean)*prev_bar)/fmem_mean))
memory_correction_plot_data <- injection_level_results_memory_corrected_raw %>% ungroup() %>% na.omit() %>% 
  rename(uncorrected = value, corrected = memcorr) %>%
  gather(memory,value,corrected,uncorrected) %>%
  filter(position %in% sample(c(1:max(position)),6))
memory_correction_plot <- ggplot(memory_correction_plot_data %>% filter(variable=="dD_mean_short"),aes(x=run_n,y=value)) +
  geom_point(aes(color=memory)) +
  facet_wrap(~position,scales="free") +
  labs(x="Injection in Sequence (#)",y="dD")
memory_correction_plot
ggsave("sample_plot_of_memory_correction.png",memory_correction_plot)


injection_level_results <- injection_level_results_memory_corrected_raw %>% 
  select(-c(value:prev_bar)) %>% 
  spread(variable,memcorr) %>% 
  full_join(injection_level_results_raw %>% 
              rename(dD_mean_short_raw = dD_mean_short,
                     d17O_mean_raw = d17O_mean,
                     d18O_mean_raw = d18O_mean),
            by = c("position", "sample", "sample_type", "run_n", "injection")) %>% 
  na.omit()

vars$n_injections = 3 #The number of injections from a vial to use for calculating vial-level summary, counting backwards from the last injection.

vial_level_results <- injection_level_results %>% 
  group_by(position,sample,sample_type) %>% 
  filter(injection > max(injection)-vars$n_injections) %>% # Use only the last X injections of a vial, where X is vars$n_injections
  summarize(vial_start = min(datetime_start),
            vial_end = max(datetime_start) + duration[which(datetime_start==max(datetime_start))],
            vial_dxs_flag = ifelse(mean(dxs_mean_short)<vars$dxs_threshold,"[Deuterium Excess is Low] ",""),
            vial_dD_mean_short = mean(dD_mean_short),
            vial_dD_mean_short_raw = mean(dD_mean_short_raw),
            vial_dD_sd_short = sd(dD_mean_short),
            vial_d17O_mean = mean(d17O_mean),
            vial_d17O_mean_raw = mean(d17O_mean_raw),
            vial_d17O_sd = sd(d17O_mean),
            vial_d18O_mean = mean(d18O_mean),
            vial_d18O_mean_raw = mean(d18O_mean_raw),
            vial_d18O_sd = sd(d18O_mean),
            vial_rawD17O_mean = mean(rawD17O_mean),
            vial_dxs = vial_dD_mean_short-8*vial_d18O_mean,
            vial_injections = length(position),
            vial_meaninject = mean(run_n),
            vial_duration = sum(duration),
            vial_h2o_mean = mean(h2o_mean),
            vial_h2o_sd = sd(h2o_mean),
            vial_h2o_mean_short = mean(h2o_mean_short),
            vial_h2o_sd_short = sd(h2o_mean_short),
            vial_cavity_torr_mean = mean(cavity_torr_mean),
            vial_cavity_torr_sd = sd(cavity_torr_mean),
            vial_cavity_temp_mean = mean(cavity_temp_mean),
            vial_cavity_temp_sd = sd(cavity_temp_mean),
            vial_residuals = mean(residual_mean),
            vial_baseshift = mean(baselineshift_mean),
            vial_basecurve = mean(baselinecurvative_mean),
            vial_slopeshift = mean(slopeshift_mean)) %>% 
  ungroup() %>% 
  mutate(standards_residuals_mean = mean(vial_residuals[which(sample_type=="standard")]),
         standards_residuals_sd = sd(vial_residuals[which(sample_type=="standard")]),
         residuals_upper = standards_residuals_mean + vars$residual_threshold*standards_residuals_sd,
         residuals_lower = standards_residuals_mean - vars$residual_threshold*standards_residuals_sd,
         standards_baseshift_mean = mean(vial_baseshift[which(sample_type=="standard")]),
         standards_baseshift_sd = sd(vial_baseshift[which(sample_type=="standard")]),
         baseshift_upper = standards_baseshift_mean + vars$baseline_shift_threshold*standards_baseshift_sd,
         baseshift_lower = standards_baseshift_mean - vars$baseline_shift_threshold*standards_baseshift_sd,
         standards_basecurve_mean = mean(vial_basecurve[which(sample_type=="standard")]),
         standards_basecurve_sd = sd(vial_basecurve[which(sample_type=="standard")]),
         basecurve_upper = standards_basecurve_mean + vars$baseline_curvature_threshold*standards_basecurve_sd,
         basecurve_lower = standards_basecurve_mean - vars$baseline_curvature_threshold*standards_basecurve_sd,
         standards_slopeshift_mean = mean(vial_slopeshift[which(sample_type=="standard")]),
         standards_slopeshift_sd = sd(vial_slopeshift[which(sample_type=="standard")]),
         residuals_flag = ifelse(vial_residuals < residuals_lower & sample_type=="sample" | vial_residuals > residuals_upper ,T,F),
         baseshift_flag = ifelse(vial_baseshift < baseshift_lower | vial_baseshift > baseshift_upper,T,F),
         basecurve_flag = ifelse(vial_basecurve < basecurve_lower | vial_basecurve > basecurve_upper,T,F),
         vial_organics_flag = ifelse(residuals_flag|baseshift_flag|basecurve_flag,"[Organic Spectral Contamination]","")) %>% 
  select(-c(vial_residuals:basecurve_flag)) %>% 
  select(position,sample,sample_type,vial_start,vial_end,vial_dxs_flag,vial_organics_flag,vial_meaninject,
         vial_dD_mean_short:vial_cavity_temp_sd) %>% 
  gather(variable,value,vial_dD_mean_short:vial_cavity_temp_sd) %>% 
  mutate(variable = substring(variable,first=6)) %>% 
  spread(variable,value) # %>% filter(vial_meaninject !=281)

# ggplot(vial_level_results %>% ungroup() %>% mutate(rel_h2o = h2o_mean/mean(h2o_mean)*100),aes(x=vial_meaninject,y=rel_h2o)) +
#   geom_line()

# Corrected Results
drift_functions <- vial_level_results %>% 
  filter(sample_type=="drift" & vial_meaninject > 15) %>% 
  select(vial_meaninject,dD_mean_short,d18O_mean,d17O_mean,dD_mean_short_raw,d18O_mean_raw,d17O_mean_raw) %>% 
  gather(variable,value,-vial_meaninject) %>% 
  mutate(memcorr = ifelse(substr(variable,nchar(variable)-2,nchar(variable))=="raw","raw","memcorr")) %>% 
  separate(variable,"variable",sep="_",extra="drop") %>% 
  spread(memcorr,value) %>% 
  group_by(variable) %>% 
  mutate(slope = lm(memcorr~vial_meaninject)$coefficients[2],
         slope_raw = lm(raw~vial_meaninject)$coefficients[2])

# # # # Drift Injections - Run this to help suss out any weirdness in the drift_function_plot
# ggplot(injection_level_results %>% filter(sample_type =="drift") %>% group_by(position),aes(x=run_n,y=dxs_mean_short)) +
#   geom_point()

# Use this plot to determine if the final values should be drift corrected. Typically, the answer is 'yes'.
drift_function_plot <- ggplot(drift_functions,aes(x=vial_meaninject,y=memcorr)) + geom_point() + geom_smooth(se=F,method='lm') + 
  facet_wrap(~variable,scales="free") + labs(x="Sequence Injection #",y="\u2030")
drift_function_plot

# Set the 'drift_correct' variable to T/F (i.e., TRUE/FALSE) to determine if it is to be drift-corrected.
drift_correction_table <- data.frame(variable =      c("d17O","d18O","dD"),
                                     drift_correct = c(F     ,F     ,F)) %>%   mutate(variable = as.character(variable))

ggsave("drift_functions.pdf",drift_function_plot)

# Calibration Plot (should be a very straight line)
cal_plot <- vial_level_results %>%
  filter(sample_type=="standard"|sample_type=="control") %>%
  select(vial_meaninject,sample,d18O_mean,d17O_mean,dD_mean_short) %>%
  left_join(read_excel("Picarro_Tray_Sheet.xlsx",sheet=2),by="sample") %>%
  select(-c(D17O_known,dxs_known)) %>%
  gather(variable,value,d18O_mean,dD_mean_short,d17O_mean,d17O_known,d18O_known,dD_known) %>%
  separate(variable,c("variable","type"),sep="_",extra="drop") %>%
  spread(type,value)


ggplot(cal_plot,aes(x=mean,y=known)) +
  geom_point() + 
  facet_wrap(~variable,scales="free") +
  geom_smooth(se=F,method='lm')

calibration_functions <- vial_level_results %>% 
  filter(sample_type=="standard" & sample!="USGS45") %>% 
  left_join(read_excel("Picarro_Tray_Sheet.xlsx",sheet=2),by="sample") %>% 
  select(vial_meaninject,sample,d18O_mean,d18O_known,d17O_mean,d17O_known,dD_mean_short,dD_known) %>% 
  gather(variable,value,d18O_mean:dD_known) %>% 
  separate(variable,c("variable","class"),sep="_",extra="drop") %>% 
  left_join(drift_correction_table,by = "variable") %>% 
  spread(class,value) %>% 
  full_join(drift_functions %>% select(variable,slope) %>% distinct(),by="variable") %>% 
  mutate(drift_corr =  mean - (slope*vial_meaninject)) %>% 
  group_by(variable,drift_correct) %>% 
  summarize(drift_slope = unique(slope),
            driftcorr_normal_slope = lm(known~drift_corr)$coefficients[2],
            driftcorr_normal_intercept = lm(known~drift_corr)$coefficients[1],
            normal_slope = lm(known~mean)$coefficients[2],
            normal_intercept = lm(known~mean)$coefficients[1])

calibrated_control <- vial_level_results %>% 
  filter(sample_type=="control" | (sample=="USGS45")) %>% 
  group_by(vial_meaninject,sample) %>% 
  select(vial_meaninject,sample,sample_type,dD_mean_short,d18O_mean,d17O_mean) %>% 
  gather(variable,value,dD_mean_short:d17O_mean) %>% 
  separate(variable,"variable",sep="_",extra="drop") %>% 
  full_join(calibration_functions,by="variable") %>% 
  mutate(drift_corrected = value - (drift_slope*vial_meaninject),
         normalized = drift_corrected * driftcorr_normal_slope + driftcorr_normal_intercept,
         normalized_nodriftcorr = value * normal_slope + normal_intercept) %>% 
  ungroup() %>% 
  mutate(calibrated_value = ifelse(drift_correct,normalized,normalized_nodriftcorr)) %>%
  select(vial_meaninject,sample,sample_type,variable,calibrated_value) %>% 
  spread(variable,calibrated_value) %>% 
  mutate(dxs = dD - 8*d18O,
         D17O = (log(d17O/1000+1)-0.528*log(d18O/1000+1))*1000*1000) %>% 
  gather(variable,value,d17O:D17O) %>% 
  left_join(read_excel("Picarro_Tray_Sheet.xlsx",sheet=2) %>% 
              gather(variable,known,-sample) %>% 
              separate(variable,'variable',sep="_",extra="drop"),by=c("sample","variable")) %>% 
  group_by(sample,sample_type,variable) %>% 
  summarize(Observed_Mean = mean(value),
            Observed_Standard_Deviation = sd(value),
            Accepted_Value = unique(known),
            Root_Mean_Square_Error = sqrt(sum((known-value)^2)/(n()-1)),
            Mean_Signed_Difference = sum(value-known)/n()) %>% 
  gather(metric,value,Observed_Mean:Mean_Signed_Difference) %>% 
  mutate(metric = factor(metric,levels=c("Observed_Mean","Observed_Standard_Deviation","Accepted_Value",
                                         "Root_Mean_Square_Error","Mean_Signed_Difference"),ordered=T,
                         labels=c("Observed Mean","Observed Standard Deviation","Accepted Value",
                                  "Root Mean Square Error","Mean Signed Difference"))) %>% 
  spread(variable,value) %>% 
  mutate(d17O = round(d17O,4),
         d18O = round(d18O,4),
         dD = round(dD,2),
         dxs = round(dxs,2),
         D17O = round(D17O,0)) %>% 
  gather(variable,value,d17O:dxs) %>% 
  mutate(variable = factor(variable,levels=c("d18O","d17O","dD","dxs","D17O"),ordered=T,
                           labels=c("\u03B4\u00B9\u2078O (\u2030)",
                                    "\u03B4\u00B9\u2077O (\u2030)",
                                    "\u03B4D (\u2030)",
                                    "D-Excess (\u2030)",
                                    "\u00B9\u2077O-Excess (10\u2076 x \u03B4)"))) %>% 
  spread(variable,value) %>% 
  rename(`Control Standard` = sample,
         `Metric` = metric)

calibrated_samples <- vial_level_results %>% 
  filter(sample_type=="sample") %>% 
  group_by(vial_meaninject,sample) %>% 
  select(vial_meaninject,sample,vial_dxs_flag,vial_organics_flag,dD_mean_short,d18O_mean,d17O_mean) %>% 
  gather(variable,value,dD_mean_short:d17O_mean) %>% 
  separate(variable,"variable",sep="_",extra="drop") %>% 
  full_join(calibration_functions,by="variable") %>% 
  mutate(raw = value,
         drift_corrected = value - (drift_slope*vial_meaninject),
         normalized = drift_corrected * driftcorr_normal_slope + driftcorr_normal_intercept,
         normalized_nodriftcorr = value * normal_slope + normal_intercept) %>% 
  ungroup() %>% 
  mutate(calibrated_value = ifelse(drift_correct,normalized,normalized_nodriftcorr)) %>%
  select(sample,vial_organics_flag,variable,calibrated_value) %>% 
  spread(variable,calibrated_value) %>% 
  mutate(dxs = dD - 8*d18O,
         D17O = (log(d17O/1000+1)-0.528*log(d18O/1000+1))*1000*1000) %>% 
  mutate(d17O = round(d17O,2),
         d18O = round(d18O,2),
         dD = round(dD,2),
         dxs = round(dxs,2),
         D17O = round(D17O,0),
         vial_dxs_flag = ifelse(dxs < vars$dxs_threshold,"[Deuterium Excess is Low] ","")) %>% 
  unite(vial_flags,c(vial_dxs_flag,vial_organics_flag),sep="") %>% 
  gather(variable,value,d17O:D17O) %>% 
  mutate(variable = factor(variable,levels=c("d18O","d17O","dD","dxs","D17O"),ordered=T,
                           labels=c("\u03B4\u00B9\u2078O (\u2030)",
                                    "\u03B4\u00B9\u2077O (\u2030)",
                                    "\u03B4D (\u2030)",
                                    "D-Excess (\u2030)",
                                    "\u00B9\u2077O-Excess (10\u2076 x \u03B4)"))) %>% 
  spread(variable,value)

sample_report <- calibrated_samples %>% 
  rename(`Sample ID` = sample, `Vial Flags` = vial_flags)

report_readme <- data.frame(Sheet = c(rep("Sheet 1",7),rep("Sheet 2",7)),
                            Column = c("Sample ID","Vial Flags","\u03B4\u00B9\u2078O (\u2030)","\u03B4\u00B9\u2077O (\u2030)","\u03B4D (\u2030)",
                                       "D-Excess (\u2030)","\u00B9\u2077O-Excess (10\u2076 x \u03B4)","Control Standard","Metric",
                                       "\u03B4\u00B9\u2078O (\u2030)","\u03B4\u00B9\u2077O (\u2030)","\u03B4D (\u2030)",
                                       "D-Excess (\u2030)","\u00B9\u2077O-Excess (10\u2076 x \u03B4)"),
                            Explanation = c("Sample identifier",
                                            "Flags for two potential problems: low D-excess indicating evaporation or spectral interference from dissolved organic carbon",
                                            "All isotope measurements are reported on the VSMOW-SLAP scale with two-point normalization.","","","",
                                            "Note the 'per-meg' scaling here","Control standard identifier",
                                            "Precision and accuracy metrics. Typically, the most appropriate to report is the Root Mean Square Error",
                                            "","","","",""))


qaqc_tracking <- vial_level_results %>% 
  filter(sample_type=="drift"|sample_type=="control"|sample_type=="control_std") %>% 
  group_by(vial_meaninject,sample) %>% 
  select(vial_meaninject,sample,sample_type,vial_dxs_flag,vial_organics_flag,dD_mean_short,d18O_mean,d17O_mean) %>% 
  gather(variable,value,dD_mean_short:d17O_mean) %>% 
  separate(variable,"variable",sep="_",extra="drop") %>% 
  full_join(calibration_functions,by="variable") %>% 
  mutate(raw = value,
         drift_corrected = value - (drift_slope*vial_meaninject),
         normalized = drift_corrected * driftcorr_normal_slope + driftcorr_normal_intercept,
         normalized_nodriftcorr = value * normal_slope + normal_intercept) %>% 
  ungroup() %>% 
  select(sample,sample_type,vial_meaninject,variable,normalized) %>% 
  spread(variable,normalized) %>% 
  mutate(dxs = dD - 8*d18O,
         D17O = (log(d17O/1000+1)-0.528*log(d18O/1000+1))*1000*1000) %>% 
  mutate(d17O = round(d17O,4),
         d18O = round(d18O,4),
         dD = round(dD,2),
         dxs = round(dxs,2),
         D17O = round(D17O,0),
         rundate = as.Date(min(vial_level_results$vial_start),format="YYYY-MM-DD")) %>% 
  select(rundate,sample:D17O)


# Data Export

# Uncalibrated Raw pulse peaks
write.csv(pulse_level_results,paste(rundate,"_pulses.csv",sep=""),row.names=F)
# Uncalibrated Injection level data
write.csv(injection_level_results,paste(rundate,"_injections.csv",sep=""),row.names=F)
# Uncalibrated Vial level data
write.csv(vial_level_results,paste(rundate,"_vials.csv",sep=""),row.names=F)
# Calibrated Sample Data
write_xlsx(calibrated_samples,paste(rundate,"_CorrectedSamples.xlsx",sep=""))
# Calibrated Sample Report
write_xlsx(list(sample_report,calibrated_control,report_readme),paste(rundate,"_Report.xlsx",sep=""))

# QAQC Tracking
setwd("C:/Box Sync/Konecky Lab/Data/Picarro L2140/")
qaqc_master <- read_excel("QAQC Tracking.xlsx") %>% 
  rbind(qaqc_tracking) %>% distinct()
write_xlsx(qaqc_master,"QAQC Tracking.xlsx")


#####                                                                         #####
##### Below are option/diagnostic plots of your data. Feel free to stop here! #####
#####                                                                         #####

# GMWL Plot of Samples
{
  ggplot(calibrated_samples %>% select(sample,variable,'Corrections: Drift & Calibrated') %>%
           spread(variable,'Corrections: Drift & Calibrated') %>% .[,c(1,2,4)] %>%
           setnames(.,c("sample","d18O","dD")),aes(x=d18O,y=dD)) +
    geom_point() +
    geom_abline(aes(slope=8.2,intercept=11.27))

  ggplot(calibrated_samples %>% select(sample,variable,'Corrections: Drift & Calibrated') %>%
           spread(variable,'Corrections: Drift & Calibrated') %>% .[,c(1,2,4,5,6)] %>%
           setnames(.,c("sample","d18O","dD","D_Excess","D17O")),
         aes(x=D_Excess,y=D17O)) +
    geom_point()
}
# Plot to look at a single pulse for demo purposes...
{
  single_pulse_plot <- ggplot(data_pulse_peaks_4 %>% ungroup() %>% select(peak_number,is_peak,mean_datetime,datetime,H2O,d18O,dD,D17O) %>%
                                mutate(D17O = D17O*1000) %>%
                                gather(var,val,H2O:D17O) %>%
                                filter(datetime > (min(mean_datetime,na.rm=T)+(3600*1.1)) &
                                         datetime < (min(mean_datetime,na.rm=T)+(3600*1.3))) %>%
                                mutate(var_label = factor(var,levels=c("H2O","d18O","dD","D17O"),ordered=T,
                                                          labels=c("H2O (ppm)","\u03B4\u00B9\u2078O (permil)","\u03B4D (permil)","\u0394\u00B9\u2077O (permeg)"))),
                              aes(x=datetime,y=val)) +
    geom_path(aes(color=is_peak,group=1)) +
    facet_wrap(~var_label,ncol=1,scales="free_y",strip.position="left") +
    scale_x_datetime(breaks = "5 min",date_labels = "%I:%M") +
    theme(strip.placement = "outside",
          strip.text.y = element_text(angle=180),
          legend.position = "none") +
    labs(y=NULL,x="Time")
  ggsave("single_pulse_plot.png",single_pulse_plot,width=150,height=125,units="mm")
}

# Code to produce a memory coefficient plot
{
  require(broom)
  fmem_curve <- fmem_terms %>% group_by(variable) %>%
    mutate(n = n()) %>% group_by(variable,n) %>%
    mutate(fmem_mean = log(fmem_mean),
           injection = log(injection)) %>%
    do(lm = lm(.$fmem_mean ~ poly(.$injection,4,raw=T))) %>%
    tidy(lm) %>% select(variable:estimate) %>%
    spread(term,estimate) %>%
    setNames(.,c("variable","n","int","a","b","c","d")) %>%
    mutate(injection = list(seq(1,n[1],0.01))) %>% unnest() %>%
    mutate(x = log(injection),
           fmem_mean = exp(int + a*x + b*x^2 + c*x^3 + d*x^4)) %>%
    ungroup() %>% select(variable,injection,fmem_mean) %>%
    mutate(variable=factor(variable,levels=c("d18O_mean","d17O_mean","dD_mean_short"),ordered=T,
                           labels=c("\u03B4\u00B9\u2078O (\u2030)","\u03B4\u00B9\u2077O (\u2030)","\u03B4D (\u2030)")))
  
  fmem_plot <- ggplot(fmem_terms %>% ungroup() %>% mutate(variable=factor(variable,levels=c("d18O_mean","d17O_mean","dD_mean_short"),ordered=T,
                                                                          labels=c("\u03B4\u00B9\u2078O (\u2030)","\u03B4\u00B9\u2077O (\u2030)","\u03B4D (\u2030)"))),
                      aes(x=(injection))) +
    geom_point(aes(y=fmem_mean*100),size=0.25) +
    geom_errorbar(aes(ymin = (fmem_mean-fmem_sd)*100, ymax = (fmem_mean+fmem_sd)*100),width=0.25,size=0.25) +
    geom_line(data=fmem_curve,aes(y=fmem_mean*100),size=0.25) +
    facet_wrap(~variable) +
    labs(x="Injection #",
         y="Contribution of Current Vial to Injection (%)")
  ggsave("memory_terms_plot.png",fmem_plot)
}
