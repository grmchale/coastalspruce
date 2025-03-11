source("Functions/lecospectR.R")

#Read in data
spec_library<-read.csv("./R_outputs/speclib_dendrometers/dendrometer_canopy_bytreeid_5nm.csv")
dendro_attributes<-read.csv("./R_outputs/speclib_dendrometers/dendro_attributetable_hyp.csv")
#Join attributes to spectral library
spec_library <- spec_library %>%
  left_join(select(dendro_attributes, TreeID, DBH, Dndrmtr, Site), by = "TreeID")

colnames(spec_library)
#Set seed for stable cal/val split
set.seed(1234)

#List of band names
band_names<-colnames(spec_library[,29:ncol(spec_library)])

#List columns names available for modeling 
  str(spec_library[,11:28])
#Set the response variable for modeling
className = "DBH"
#Filter data to include rows with the response variable

    #slice_sample(n =15, replace = F)
spec_library_n25<-spec_library[is.na(spec_library[className])==F,] %>% 
    subset(DBH != "not sampled") %>% #dplyr::filter(GenusSpecies == "White_Spruce") %>%  select(Site, TreeID) %>% unique() %>% dim#tally
    mutate(DBH = as.numeric(DBH), Vigor_class = as.factor(Vigor_class), GenusSpecies = as.factor(GenusSpecies)) %>%
    group_by(Site, TreeID, GenusSpecies) %>% #tally %>% dplyr::select(n) %>% ungroup() %>% summarise(min_pix= min(n), max_pix = max(n), median_pix = median(n))
    slice_sample(n =80, replace = F)

  spec_library_n25 %>% group_by(Site, TreeID, GenusSpecies, eval(className)) %>% tally() %>% print(n=100)
  unique(spec_library_n25$DBH)
  #Create a test and train split
  inTrain <- caret::createDataPartition(
    y = spec_library_n25[[className]],
    p = 0.7,
    list = FALSE#,
    #na.rm = TRUE
  )
  
training <- spec_library_n25[inTrain,]  %>% #tally() %>% print(n=200)
    ungroup %>% dplyr::select(className, band_names)
testing <- spec_library_n25[-inTrain,]  %>% 
    ungroup %>% dplyr::select(className, band_names)


#Alternatively, use cal and val polygons  
#cal<-spec_library_n25 %>%
#    dplyr::filter(Canopy_Type == "cal")  %>% ungroup %>% dplyr::select(className, band_names) 
#val<-spec_library_n25 %>%
#    dplyr::filter(Canopy_Type == "val")  %>% ungroup %>% dplyr::select(className, band_names)

#dim(cal)
#dim(val) 
head(training[,1])

################# Ranger models aka RF
n=1000
     rf_mod <- ranger::ranger(as.formula(paste(className, "~.")),
    data = training,num.trees = n)
    rf_mod_pred<-predict(rf_mod, testing)
    rf_mod$confusion.matrix
    caret::confusionMatrix(data = rf_mod_pred$predictions, reference = testing$Vigor_class)
    windows()
    graphics::plot(rf_mod_pred$predictions, testing$DBH, type="p")
    #plot(hexbin::hexbin(rf_mod_pred$predictions, testing$Vigor_class))
    #abline(lm(rf_mod_pred$predictions~ testing$dbh_cm))
    #abline(0,1)
    R2(rf_mod_pred$predictions, testing$dbh_cm, formula = "corr")
 
 #Determine outliers based on obs vs pred by variable and identify tree and site for QAQC




 
#################Partial least squares regression 
  #tune model: 10-fold cross-validation repeated 3 times
  ctrl <- caret::trainControl(
    method = "repeatedcv",
    number = 10,
    #sampling = "down",
    repeats = 3)
  ncomp = 30
  #Fit model. Note max iterations set to 100000 to allow model convergence
  plsFit <- caret::train(
    as.formula(paste(className, "~.")),
    data = training,
    maxit = 10000,
    method = "pls",
    trControl = ctrl,
    tuneLength = ncomp)

   plsFit_pred<- predict(plsFit, newdata = testing)

   plsFit_pred
   windows()
   plot(hexbin::hexbin(plsFit_pred,testing[[className]]))
    #corrplot::corrplot(plsFit_pred)    
   #windows()
   caret::confusionMatrix(plsFit)
   accuracy = sum(ifelse(plsFit_pred==val$GenusSpecies, 1,0))/length(val$GenusSpecies)
    accuracy    
    plot(plsFit_pred ~ testing$Average_Condensed_TANNINS_percent_DW)
    summary(lm(plsFit_pred ~ testing$Average_Condensed_TANNINS_percent_DW))
    R2(plsFit_pred,testing$Fv_Fm)

caret::R2(plsFit_pred, val$GenusSpecies)
caret::RMSE()
xyplot(plsFit)
median(plsFit$resample$Rsquared)













