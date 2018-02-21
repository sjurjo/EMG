
# EMG function 

myEMG <- function(filepath, file, data.name = "VLL"){
  

  emg <- read.csv(filepath, sep =";")
  
  # Change emg to emg object
  emg.x <- emg$VLL 
  
  
  emgraw <- as.emg(emg.x, samplingrate = 1500, units = "mV", date.name = data.name) # sampling rate 1500Hz to avoid signal loss
  
  
  if(is.emg(emgraw) == FALSE) stop(paste0("The file (", as.character(file), ") does not turn out to be an EMG object")) # EMG object check
  
  # High pass filter
  x <- highpass(emgraw, cutoff = 20) # 20 Hz recommended for sport activities (maybe as much as 30 Hz)
  xlh <- lowpass(x, cutoff = 400, n = 5) # 500 Hz most common for high frequencies, n = 5 for signal delay
  emghpass <- rectification(xlh)
  
  #EMG Moving average
  emgma <- envelope(emghpass, method = "MA", wsize = 100) # wsize; 100 ms, most common for all conditions
  
  detected<-onoff_singlethres(emgma, 1, emgma, t = 20)
  
  #Compute the length of each phase detected
  res<-phasestats(emghpass,detected,length)
  
  #Compute the mean of each phase detected
  res<-phasestats(emghpass,detected,mean)

  
  # Only show the means of the active phases
  x <- res$stats[names(res$stats)==1]
  
  # ID variables
  
  id.var <- data.frame(str_split_fixed(file, "_", 3))
  
  id.var$X3 <- gsub(".csv", "", id.var$X3)
  
  colnames(id.var) <- c("subject", "condition", "drag")
  
  # Data frame
  m <- as.data.frame(x)  
  colnames(m)[colnames(m) == 'x'] <- 'mV' 
  m$cycles <- as.numeric(rownames(m))
  m$drag <- id.var[1,3]
  m$subj <- id.var[1,1]
  m$cond <- id.var[1,2]
  m <- m %>% 
    filter(!(mV <0)) %>% 
    filter(!(cycles > 40))

  
  return(m)

} 

fp <- "./data/complete_test_na/"

## List files to be used
file.list <- list.files("./data/complete_test_na/") 

## Store results
results <- list()

# Loop over files
for(i in 1:length(file.list)) {
  

  fp <-  paste0("./data/complete_test_na/",file.list[i])
  print(file.list[i])
  results[[i]] <- myEMG(filepath = fp, file = file.list[i])
 
  }

results <- bind_rows(results)
# Convert to one dataframe and store df's
results <- as.data.frame(results)
data.frame(results)


### mean per subject per drag

meanperfp <- results %>% 
  group_by(drag, subj, cond) %>%
  mutate(m.subj = mean(mV)) %>%
  summarise(m.perdrag = mean(m.subj)) %>%
  filter(!(drag %in% c("n1", "n2", "n2u", "n1u", "n1v", "n2v"))) %>% #removing rows for tidyverse assumptions
  ungroup() %>%
  print()
  
norm <- results %>% 
  group_by(drag, subj, cond) %>%
  mutate(m.subj = mean(mV)) %>%
  summarise(m.perdrag = mean(m.subj)) %>%
  ungroup() %>%
  filter(drag %in% c("n1", "n2", "n2u", "n1u")) %>% 
  mutate(norm = if_else(drag %in% c("n1", "n1u"), "norm200", "norm325"))%>%
  dplyr::select(norm,cond, subj, m.perdrag) %>%
  spread(norm, m.perdrag) %>%
  inner_join(meanperfp) %>%
  data.frame() %>% 
  print()

# normalized to 200 WITHOUT vibration for both conditions (nwov)
nwov <- norm %>%
  spread(drag, m.perdrag) %>% 
  mutate(no.d1 = (d1-norm200)/norm200*100,
         no.d2 = (d2-norm200)/norm200*100,
         no.d3 = (d3-norm200)/norm200*100,
         no.d4 = (d4-norm200)/norm200*100,
         no.d5 = (d5-norm200)/norm200*100,
         no.d6 = (d6-norm200)/norm200*100) %>% 
  gather(drag, norm.values, no.d1:no.d6) 

# ggplot
library(devtools)
install_github("dhammarstrom/publR", build_vignettes = TRUE)


ggplot(nwov,aes(drag, norm.values, color = cond)) + geom_boxplot()


nwovplot <- nwov %>%
  group_by(drag, cond) %>% 
  summarise(m = mean(norm.values),
            s = sd(norm.values)) %>% 
  print()

  ggplot(nwovplot,aes(drag, m, group = cond)) + geom_point(aes(shape = cond), position = position_dodge(0.4), size = 2) + geom_line(aes(linetype=cond), position = position_dodge(0.4)) +
  geom_errorbar(aes(ymax=m+s, ymin=m-s), width = 0.3, position = position_dodge(0.4))+ 
    geom_segment(aes(x=0.7,xend=6.3,y=67.7, yend=67.7)) +
    geom_segment(aes(x = 0.7, y = 67, xend = 0.7, yend = 67.7)) +
    geom_segment(aes(x = 6.3, y = 67, xend = 6.3, yend = 67.7)) +
    scale_y_continuous(breaks = seq(0, 65, by = 10)) +
    labs(x="Intervals",y="mV Change (%)")+
    scale_x_discrete(labels= c("no.d1" = "I.1", "no.d2" = "I.2", "no.d3" ="I.3", "no.d4" ="I.4", "no.d5" ="I.5", "no.d6" ="I.6")) + 
    publr_theme(font_family = "Times", axis_line_width = 0.2, font_size_labels = 10, font_size_axis = 11, text_color = "black") +
    annotate("text", x= 3.5, y = 68, label = "*", size = 5)

    
    
    # analysis nwov


nwov$norm.values <- scale(nwov$norm.values)
  

library(lme4); library(lmerTest)
 
norm200wov <- lmerTest::lmer(norm.values ~  drag * cond + (1|subj), data = nwov)
anova(norm200wov)

temp<-nwov %>%
  group_by(subj, cond)%>%
  summarise(m = mean(norm.values, na.rm=T))%>%
  spread(cond, m) %>%
  print()
  

t.test(temp$uvib, temp$vib, paired = TRUE)





anova(norm200wov)


aov.mod <- aov(norm.values ~ drag * cond + Error(subj), data = nwov)

step()
summary(aov.mod)


?step


lmerTest::anova(norm200wov)

summary(norm200wov)
plot(norm200wov)
qqnorm(residuals(norm200wov))
qqline(residuals(norm200wov))



## normalized to 325 WITHOUT vibration for both conditions
     
n2 <- norm %>%
  spread(drag, m.perdrag) %>% 
  mutate(no.d1 = (d1-norm325)/norm325*100,
         no.d2 = (d2-norm325)/norm325*100,
         no.d3 = (d3-norm325)/norm325*100,
         no.d4 = (d4-norm325)/norm325*100,
         no.d5 = (d5-norm325)/norm325*100,
         no.d6 = (d6-norm325)/norm325*100) %>% 
  gather(drag, norm.values, no.d1:no.d6) %>% 
  group_by(drag, cond) %>% 
  summarise(m = mean(norm.values),
            s = sd(norm.values)) %>% 
  print()

ggplot(n,aes(drag, norm.values, color = cond)) + geom_boxplot()

ggplot(n2,aes(drag, m, group = cond, color = cond)) + geom_point(position = position_dodge(0.4)) + geom_line(position = position_dodge(0.4)) +
  geom_errorbar(aes(ymax=m+s, ymin=m-s), width = 0.3, position = position_dodge(0.4))


norm325wov <- lmerTest::lmer(norm.values ~ drag * cond + (drag|subj), data = n2)
summary(norm325wov)
plot(norm325wov)
qqnorm(residuals(norm325wov))
qqline(residuals(norm325wov))

## Normalize to 200 w with vibration for condition vib and 200 w without vibration for condition without vib
#(325 w with vibration missing in FP 8)

normw <- results %>% 
  group_by(drag, subj, cond) %>%
  mutate(m.subj = mean(mV)) %>%
  summarise(m.perdrag = mean(m.subj)) %>%
  ungroup() %>%
  filter(drag %in% c("n1", "n2", "n2v", "n1v")) %>%
  mutate(norm = if_else(drag %in% c("n1", "n1v"), "norm200", "norm325"))%>%
  dplyr::select(norm,cond, subj, m.perdrag) %>%
  spread(norm, m.perdrag) %>%
  inner_join(meanperfp) %>%
  data.frame() %>% 
  print()

nw <- normw %>%
  spread(drag, m.perdrag) %>% 
  mutate(no.d1 = (d1-norm200)/norm200*100,
         no.d2 = (d2-norm200)/norm200*100,
         no.d3 = (d3-norm200)/norm200*100,
         no.d4 = (d4-norm200)/norm200*100,
         no.d5 = (d5-norm200)/norm200*100,
         no.d6 = (d6-norm200)/norm200*100) %>% 
  gather(drag, norm.values, no.d1:no.d6) %>% 
  group_by(drag, cond) %>% 
  summarise(m = mean(norm.values),
            s = sd(norm.values)) %>% 
  print()


ggplot(nw,aes(drag, m, group = cond, color = cond)) + geom_point(position = position_dodge(0.4)) + geom_line(position = position_dodge(0.4)) +
  geom_errorbar(aes(ymax=m+s, ymin=m-s), width = 0.3, position = position_dodge(0.4))



ggplot(nw,aes(drag, norm.values, color = cond)) + geom_boxplot()


norm200wv <- lmerTest::lmer(norm.values ~ drag * cond + (1|subj), data = nw)
summary(norm200wv)
plot(norm200wv)
qqnorm(residuals(norm200wv))
qqline(residuals(norm200wv))



## Muscle activation increase MAI

MAI <- norm %>%
  spread(drag, m.perdrag) %>%
  gather(drag, norm.values, d1:d6) %>%
  select(subj,drag, cond, norm.values) %>% 
  spread(cond, norm.values) %>% 
  mutate(mai = (vib-uvib)/uvib*100) %>%
  select(subj, drag, mai) %>% 
  spread(drag, mai) %>% 
  gather(drag, mai.in, d1:d6) %>% 
  group_by(drag) %>%
  summarise(m = mean(mai.in),
            s = sd(mai.in)) %>% 
  print()
  
 
ggplot(MAI,aes(drag, m)) + geom_point(position = position_dodge(0.4)) + geom_line(position = position_dodge(0.4)) +
  geom_errorbar(aes(ymax=m+s, ymin=m-s), width = 0.3, position = position_dodge(0.4))

mai <- lmerTest::lmer(mai ~ drag + (1|subj), data = MAI)
summary(mai)

plot(mai)
qqnorm(residuals(mai))
qqline(residuals(mai))   
   
   
   
   


