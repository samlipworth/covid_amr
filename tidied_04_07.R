

##############define functions#################


prepare_to_plot<- function(model){
  model_data<- smooth_samples(model,n = 100,n_vals = (365*7))
  model_data$time<-min(data$specimen_date_sgss,na.rm = T) + days(round(model_data$timeline,0))
  model_data_credible <- model_data %>%
    group_by(timeline) %>%
    summarise(lower_credible = quantile(.value, probs=c(0.025)),
              upper_credible = quantile(.value, probs=c(0.975)))
  model_data_credible$time<-min(data$specimen_date_sgss,na.rm = T) + days(round(model_data_credible$timeline,0))

  sm_est <- smooth_estimates(model,nvals=(365*7)) |>
    add_confint()
  sm_est$time<-min(data$specimen_date_sgss,na.rm = T) + days(round(sm_est$timeline,0))
  return(list(sm_est = sm_est, model_data = model_data, model_data_credible = model_data_credible))

}

annotate_derivatives<-function(object){
object$first_sig<-ifelse((object$low <0 & object$high <0) | (object$low >0 & object$high >0),1,0)
object$first_breakpoint <- 0
object$first_breakpoint[object$first_sig == 1 &
                               lag(object$first_sig == 0, 1)] <- 1


object$first_change<-0
object$first_change[(object$est <0 &
                            lag(object$est > 0, 1)) | (object$est >0 &
                                                                           lag(object$est <0,1))] <- 1


# calculate direction of change-point
object <- object %>%
  dplyr::mutate(lag_upper = dplyr::lag(high, n = 1, default = NA),
                lag_lower = dplyr::lag(low, n = 1, default = NA))
object$fdirection <- NA

object$fdirection[object$est == 1 &
                         object$low > 0 &
                         object$lag_lower < 0] <- "increasing"

object$fdirection[object$est == 1 &
                         object$high < 0 &
                         object$lag_upper > 0] <- "decreasing"
return(object)
}


############load libraries#############

library(sf)
library(spdep)
library(tidyverse)
library(mgcv)
library(gratia)
library(patchwork)
###########read in working dataset#############

ecoli_data <- read_delim('Z:/Y080_UID_PID/SamLipworth/working_dataset_new_12_06.tsv', delim = ",") %>%
 filter(spec_dt > ymd(20120101)) %>% filter(spec_dt < ymd(20240101))

kleb_data<-read_delim('Z:/Y080_UID_PID/SamLipworth/working_kleb_dataset_new_13_6.tsv',delim=",") %>%
  filter(spec_dt > ymd(20120101)) %>% filter(spec_dt < ymd(20240101))
kleb_data$age<-round(as.numeric((kleb_data$spec_dt - kleb_data$dob)/365),1)
ecoli_data$age<-round(ecoli_data$age,1)
kleb_data$monthno<-month(kleb_data$spec_dt)

ecoli_data$species<-"ecoli"
kleb_data$species<-"kleb"

data<-bind_rows(ecoli_data,kleb_data)

nuts<-read_tsv('Z:/Y080_UID_PID/SamLipworth/rng_lsoa21_nuts.tsv')
nuts<-filter(nuts,!is.na(rgn) & !is.na(lsoa21) & !is.na(nuts))
nuts<-distinct(nuts,lsoa21,.keep_all = T)
data<-left_join(data,nuts,by=c("lsoa21"))

data$nuts_name <- case_when(
  data$nuts == 'UKC' ~ 'North East',
  data$nuts == 'UKD' ~ 'North West',
  data$nuts == 'UKE' ~ 'Yorkshire and The Humber',
  data$nuts == 'UKF' ~ 'East Midlands',
  data$nuts == 'UKG' ~ 'West Midlands',
  data$nuts == 'UKH' ~ 'East of England',
  data$nuts == 'UKI' ~ 'London',
  data$nuts == 'UKJ' ~ 'South East',
  data$nuts == 'UKK' ~ 'South West',
  TRUE ~ 'Other'
)

###################basic tidying of dataset#################
data$year<-year(data$spec_dt)

#we collapse Ampicillin and Amoxicillin to be the same thing and tidy up

data$AMPAMOX <- case_when(
  data$e_AMPAMOX == 'R' ~ 'R',
  data$e_AMPAMOX == 'S' ~ 'S',
  data$e_AMP == 'R' ~ 'R',
  data$e_AMP == 'S' ~ 'S',
  data$e_AMOX == 'R' ~ 'R',
  data$e_AMOX == 'S' ~ 'S',
  TRUE ~ NA_character_
)
data<-select(data,-e_AMP,-e_AMOX,-e_AMPAMOX)


#recode so that 1 is R, 0 is S
data <- data %>%
  mutate(across(.cols = 50:69, .fns = ~case_when(
    . == "R" ~ 1,
    . == "S" ~ 0,
    TRUE ~ NA_real_  # Use NA_real_ for numeric NA
  )))

data$threegc<-case_when(data$e_CTX == 1 ~1,
                     data$e_CTZ ==1 ~ 1,
                     data$e_CEF ==1 ~ 1,
                     data$e_CTRI ==1 ~ 1,
                     data$e_CTX ==0 ~0,
                     data$e_CTZ ==0 ~0,
                     data$e_CEF == 0 ~ 0,
                     data$e_CTRI ==0 ~ 0,
  TRUE ~ NA)

#a timeline from the beginning of data collection
data$timeline<-as.numeric(data$spec_dt - min(data$spec_dt,na.rm = T))

data$AMPAMOX<-ifelse(
  data$AMPAMOX == 'S',0,
  ifelse(data$AMPAMOX =='R',1,
         NA
))

data$ym<-zoo::as.yearmon(data$spec_dt)
data$ym<-my(data$ym)
data$day<-yday(data$spec_dt)

data$healthcare_exposure<-as.factor(data$healthcare_exposure)


data$sex<- case_when(
  data$sex.x == "Female" | data$sex =="Female" ~ "Female",
  data$sex.x == "Male" | data$sex =="Male" ~ "Male",
  TRUE ~ NA
)

#############table 1################
explanatory <- c("age","imd_score","elixhauser","sex","e_AMOXCLAV","AMPAMOX","e_GENT","e_CIP","e_GENT","e_PIPTAZ","e_CXM","threegcS","healthcare_exposure")
data$year2<-as.factor(data$year)
data<-rename(data,rururb=`Rural Urban Classification 2011 (10 fold)`)

dependent="year2"

data %>% finalfit::summary_factorlist(
  dependent,
  explanatory,
  cont="median"
) -> t

t %>% kableExtra::kable(format = 'html')


##########test for seasonality####################

#fit seasonal and non seasonal models and look at lrt/dAIC
seasonality <- function(drug,bug) {
  cat(paste0("\n",drug))
  data_<-filter(data,species==bug)
  formula_seasonal <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25) + s(day, bs='cc', k=12)"))
  formula_non_seasonal <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25)"))


  seasonal_model <- bam(formula_seasonal, data=data_, family=binomial(link="logit"), method="REML")
  non_seasonal_model <- bam(formula_non_seasonal, data=data_, family=binomial(link="logit"), method="REML")


  p <- anova(seasonal_model, non_seasonal_model, test="LRT")$`Pr(>Chi)`[2]
  aic_seasonal <- AIC(seasonal_model)
  aic_non_seasonal <- AIC(non_seasonal_model)


  return(list(abx=drug,bug=bug,p_value = p, AIC_seasonal = aic_seasonal, AIC_non_seasonal = aic_non_seasonal))
}
antibiotics <- c("e_AMOXCLAV", "e_CIP", "e_PIPTAZ", "e_GENT","AMPAMOX","threegc","e_MER")
bugs<-c("ecoli","kleb")
results_list <- lapply(antibiotics, function(drug) {
  lapply(bugs, function(bug) seasonality(drug, bug))
})

results<-bind_rows(results_list)
results$psig<-ifelse(results$p_value <0.05,1,0)
results$AIC_seasonal <- as.numeric(as.character(results$AIC_seasonal))
results$AIC_non_seasonal <- as.numeric(as.character(results$AIC_non_seasonal))

results$aic_diff<-results$AIC_non_seasonal - results$AIC_seasonal
results$aic_sig<-ifelse((results$AIC_non_seasonal - results$AIC_seasonal) >10,1,0)
results$p_value<-round(results$p_value,3)
results$AIC_seasonal<-round(results$AIC_seasonal,0)
results$AIC_non_seasonal<-round(results$AIC_non_seasonal,0)
results$aic_diff<-round(results$aic_diff,0)

results %>% arrange(bug) %>%  knitr::kable() %>% clipr::write_clip()

#########load shape data#################

##############read in shape data##################
d<-data
shp<-st_read('Z:/Y080_UID_PID/SamLipworth/NUTS_1_UK/')
shp<-filter(shp,nuts118cd != 'UKL' & nuts118cd != 'UKM' & nuts118cd != 'UKN')
shp$nuts118nm<-str_replace_all(shp$nuts118nm," \\(England\\)","")
#shp<-filter(shp,nuts118nm %in% d$nuts)

d<-filter(d,nuts_name %in% shp$nuts118nm)

df <- st_drop_geometry(shp)
df <- droplevels(df)


nb <- poly2nb(shp, row.names = df$nuts118nm)
names(nb) <- attr(nb, "region.id")

#################fit national data #################

#this function fits national trends and returns first derivative changepoints

function_1<-function(drug,bug,seasonal){
  if(seasonal==TRUE){
    cat("\nfitting a seasonal model")
    formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25) + s(day, bs='cc', k=12)"))
  }else{
    cat("\nfitting a non-seasonal model")
    formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25)"))

    }
  drug_col <- sym(drug)
data_<-filter(data,species==bug)
model<-bam(formula,
                  data=data_,
                  family = binomial(link = "logit"),
                  method = "REML")


cat("\ngenerating data to predict with")
time<-seq(0,max(data_$timeline,14))
predict_data<-expand_grid(timeline=time)
predict_data$date<-min(data_$spec_dt,na.rm = T) + days(predict_data$timeline)
predict_data$day<-yday(predict_data$date)

cat("\ntaking samples from posterior distribution of the expected value of responses")
fitted<-fitted_samples(model,data = predict_data,n = 100)
cat("\ndownsampling sample to 100 for plotting")
all_fitted<-fitted #%>% group_by(.row) %>% sample_n(100)
dates<-data.frame(.row=row_number(predict_data),date=predict_data$date)
all_fitted<-left_join(all_fitted,dates,by=c(".row"))
all_fitted$bug<-bug
cat("\nGenerating credibility intervals")
fitted<- fitted %>%
  group_by(.row) %>%
  summarise(low=quantile(.fitted,probs = c(0.025)),
            high=quantile(.fitted,probs=c(0.975)),
            est=quantile(.fitted,probs=c(0.5)))

fitted$date<-predict_data$date
fitted$bug<-bug

cat("\nsummarising observed data")
observed<-data_ %>% group_by(ym) %>%
  summarise(drug = mean(!!drug_col,na.rm = T))
observed$bug<-bug

cat("\nmaking derivative samples")
deriv<-derivative_samples(model,order = 1, focal = "timeline",scale = "response",draws = 1000,data = predict_data,select="s(timeline)")
deriv<-deriv %>% group_by(.row) %>%
  summarise(low=quantile(.derivative,probs = c(0.025)),
            high=quantile(.derivative,probs=c(0.975)),
            est=quantile(.derivative,probs=c(0.5)))
deriv$date<-predict_data$date
deriv<-annotate_derivatives(deriv)
deriv$bug<-bug
deriv$drug<-drug

return(list(fitted=fitted,observed=observed,abx=drug,bug=bug,fitted_draws=all_fitted,deriv=deriv))
}



ecoli_coamox<-function_1("e_AMOXCLAV","ecoli",TRUE)
ecoli_gent<-function_1("e_GENT","ecoli",TRUE)
kleb_coamox<-function_1("e_AMOXCLAV","kleb",TRUE)
kleb_gent<-function_1("e_GENT","kleb",FALSE)
ecoli_piptaz<-function_1("e_PIPTAZ","ecoli",TRUE)
kleb_piptaz<-function_1("e_PIPTAZ","kleb",FALSE)
ecoli_cip<-function_1("e_CIP","ecoli",TRUE)
kleb_cip<-function_1("e_CIP","kleb",FALSE)
ecoli_tgc<-function_1("threegc","ecoli",TRUE)
kleb_tgc<-function_1("threegc","kleb",FALSE)
ecoli_cotrim<-function_1("e_CXM","ecoli",TRUE)
kleb_cotrim<-function_1("e_CXM","kleb",FALSE)
ecoli_mero<-function_1("e_MER","ecoli",TRUE)
kleb_mero<-function_1("e_MER","kleb",FALSE)
ecoli_amox<-function_1("AMPAMOX","ecoli",TRUE)

fitted_draws_coamox<-rbind(ecoli_coamox$fitted_draws,kleb_coamox$fitted_draws)
fitted_coamox<-rbind(ecoli_coamox$fitted,kleb_coamox$fitted)
observed_coamox<-rbind(ecoli_coamox$observed,kleb_coamox$observed)
deriv_coamox<-rbind(ecoli_coamox$deriv,kleb_coamox$deriv)

fitted_draws_amox<-rbind(ecoli_amox$fitted_draws)
fitted_amox<-rbind(ecoli_amox$fitted)
observed_amox<-rbind(ecoli_amox$observed)
deriv_amox<-rbind(ecoli_amox$deriv)

fitted_draws_gent<-rbind(ecoli_gent$fitted_draws,kleb_gent$fitted_draws)
fitted_gent<-rbind(ecoli_gent$fitted,kleb_gent$fitted)
observed_gent<-rbind(ecoli_gent$observed,kleb_gent$observed)
deriv_gent<-rbind(ecoli_gent$deriv,kleb_gent$deriv)

fitted_draws_piptaz<-rbind(ecoli_piptaz$fitted_draws,kleb_piptaz$fitted_draws)
fitted_piptaz<-rbind(ecoli_piptaz$fitted,kleb_piptaz$fitted)
observed_piptaz<-rbind(ecoli_piptaz$observed,kleb_piptaz$observed)
deriv_piptaz<-rbind(ecoli_piptaz$deriv,kleb_piptaz$deriv)

fitted_draws_cip<-rbind(ecoli_cip$fitted_draws,kleb_cip$fitted_draws)
fitted_cip<-rbind(ecoli_cip$fitted,kleb_cip$fitted)
observed_cip<-rbind(ecoli_cip$observed,kleb_cip$observed)
deriv_cip<-rbind(ecoli_cip$deriv,kleb_cip$deriv)

fitted_draws_tgc<-rbind(ecoli_tgc$fitted_draws,kleb_tgc$fitted_draws)
fitted_tgc<-rbind(ecoli_tgc$fitted,kleb_tgc$fitted)
observed_tgc<-rbind(ecoli_tgc$observed,kleb_tgc$observed)
deriv_tgc<-rbind(ecoli_tgc$deriv,kleb_tgc$deriv)

fitted_draws_cotrim<-rbind(ecoli_cotrim$fitted_draws,kleb_cotrim$fitted_draws)
fitted_cotrim<-rbind(ecoli_cotrim$fitted,kleb_cotrim$fitted)
observed_cotrim<-rbind(ecoli_cotrim$observed,kleb_cotrim$observed)
deriv_cotrim<-rbind(ecoli_cotrim$deriv,kleb_cotrim$deriv)

fitted_draws_mero<-rbind(ecoli_mero$fitted_draws,kleb_mero$fitted_draws)
fitted_mero<-rbind(ecoli_mero$fitted,kleb_mero$fitted)
observed_mero<-rbind(ecoli_mero$observed,kleb_mero$observed)
deriv_mero<-rbind(ecoli_mero$deriv,kleb_mero$deriv)

#deriv_coamox<-rbind(ecoli_coamox$deriv,kleb_coamox$deriv) %>% filter(first_breakpoint==1)

plot_p1<-function(drug){
  fitted_draws <- get(paste0("fitted_draws_", drug))
  fitted_draws$group<-paste0(fitted_draws$bug,'_',fitted_draws$.draw)
  fitted<- get(paste0("fitted_",drug))
  deriv<-get(paste0("deriv_",drug))
  deriv<-left_join(deriv,fitted,by=c("date","bug"))
  observed_data <- get(paste0("observed_", drug))
  #deriv<-get(paste0("deriv_",drug))
  p<-ggplot(fitted) +
  aes(x=date,y=est,group=bug) +
  geom_line(color="blue") +
  geom_line(data=fitted_draws,aes(x=date,y=.fitted,group=group),alpha=0.05) +
  #geom_ribbon(data=fitted_data,aes(x=date,ymin=low,ymax=high),alpha=0.1) +
  ylab("Proportion of isolates resistant") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle=-90),legend.position =ifelse(drug!="amox", "right","none")) +
  scale_x_date(breaks = seq(as.Date("2012-01-01"), as.Date("2024-01-01"), by="1 year"),
               date_labels = "%Y") +
  ggtitle("") +
  xlab("") +
  geom_line(data=observed_data,aes(x=ym,y=drug,color=bug)) +
  labs(color="Species",linetype="Changepoint") +
  scale_color_brewer(palette = "Dark2") 
  
  params<-ggplot_build(p)
  params<-params$layout$panel_params[[1]]$y$breaks
  params<-abs(params[3]-params[2])
  #geom_vline(data=filter(deriv,first_breakpoint==1),aes(xintercept=date,color=bug),linetype="dashed",ymax=0.05)
  p+ geom_segment(data=filter(deriv,first_breakpoint==1),aes(x=date,xend=(date),y=est.y-(params/2),yend=est.y+(params/2),linetype=bug),linewidth=1,color="darkred") +
    scale_linetype_manual(values=c("dashed","dotted"))
    
}

cip_p1<-plot_p1("cip") + ggtitle("Ciprofloxacin")
coamox_p1<-plot_p1("coamox") + ggtitle("Co-amoxiclav")
gent_p1<-plot_p1("gent") + ggtitle("Gentamicin")
piptaz_p1<-plot_p1("piptaz") + ggtitle("Piperacillin-Tazobactam")
tgc<-plot_p1("tgc") + ggtitle("Third Generation Cephalosporin")
cotrim_p1<-plot_p1("cotrim") + ggtitle("Co-trimoxazole")
mero_p1<-plot_p1("mero") + ggtitle("Meropenem")
amox_p1<-plot_p1("amox") + ggtitle("Amoxicillin") + theme(legened.position="none")

coamox_p1 + cip_p1 + gent_p1 + piptaz_p1 + tgc + cotrim_p1 + mero_p1 + amox_p1 +
  plot_layout(guides="collect",nrow = 2)

#this function returns proportions by year for plotting

function_2<-function(code,bug){
  code <- sym(code)
  data_<-data#filter(data,species==bug)
prop <- data_ %>%
  group_by(species,year, monthno) %>%
  summarise(S=sum(!!code == 0, na.rm = TRUE),
  R=sum(!!code == 1, na.rm = TRUE))

prop$date<-ym(paste0(prop$year,prop$monthno))
prop<-prop %>% pivot_longer(cols = 4:5,names_to = "which",values_to = "n")

p2<-ggplot(prop) +
  aes(x=date,y=n,fill=which) +
  geom_bar(position="stack", stat="identity") +
  scale_x_date(breaks = seq(as.Date("2012-01-01"), as.Date("2024-01-01"), by="2 years"),
               date_labels = "%Y") +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() + labs(fill="") +
  ylab("Number of isolates") +
  facet_wrap(~species)
return(p2)
}

coamox_p2<-function_2("e_AMOXCLAV","ecoli")
cipro_p2<-function_2("e_CIP","ecoli")
gent_p2<-function_2("e_GENT","ecoli")
piptaz_p2<-function_2("e_PIPTAZ","ecoli")
tgc_p2<-function_2("threegc","ecoli")
cotrim_p2<-function_2("e_CXM","ecoli")
mero_p2<-function_2("e_MER","ecoli")


#this function is similar to function 1 but returns values on the linear scale

function_3<-function(drug,bug,seasonal){

  if(seasonal==TRUE){
    formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25) + s(day, bs='cc', k=12)"))
  }else{
    formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25)"))

  }
data_<-filter(data,species==bug)
cat("\nfitting model")
model<-bam(formula,
           data=data_,
           family = binomial(link = "logit"),
           method = "REML")

time<-seq(0,max(data_$timeline,14))
predict_data<-expand_grid(timeline=time)
predict_data$date<-min(data_$spec_dt,na.rm = T) + days(predict_data$timeline)
predict_data$day<-yday(predict_data$date)

cat("\nposterior fitted samples from model")
fitted<-smooth_samples(model,data = predict_data,n = 1000,select="s(timeline)")
fitted<- fitted %>%
  group_by(.row) %>%
  summarise(low=quantile(.value,probs = c(0.025)),
            high=quantile(.value,probs=c(0.975)),
            est=quantile(.value,probs=c(0.5)))
fitted$date<-predict_data$date
fitted$bug<-bug

cat("\nmaking derivative samples")
deriv<-derivative_samples(model,order = 1, focal = "timeline",scale = "linear_predictor",draws = 1000,data = predict_data)
deriv<-deriv %>% group_by(.row) %>%
  summarise(low=quantile(.derivative,probs = c(0.025)),
            high=quantile(.derivative,probs=c(0.975)),
            est=quantile(.derivative,probs=c(0.5)))
deriv$date<-predict_data$date
deriv<-annotate_derivatives(deriv)
deriv$bug<-bug


return(list(deriv=deriv,fitted=fitted,drug=drug,bug=bug))
}

ecoli_coamox<-function_3("e_AMOXCLAV","ecoli",TRUE)
ecoli_gent<-function_3("e_GENT","ecoli",TRUE)
kleb_coamox<-function_3("e_AMOXCLAV","kleb",TRUE)
kleb_gent<-function_3("e_GENT","kleb",FALSE)
ecoli_piptaz<-function_3("e_PIPTAZ","ecoli",TRUE)
kleb_piptaz<-function_3("e_PIPTAZ","kleb",FALSE)
ecoli_cip<-function_3("e_CIP","ecoli",TRUE)
kleb_cip<-function_3("e_CIP","kleb",FALSE)
ecoli_tgc<-function_3("threegc","ecoli",TRUE)
kleb_tgc<-function_3("threegc","kleb",FALSE)
ecoli_cotrim<-function_3("e_CXM","ecoli",TRUE)
kleb_cotrim<-function_3("e_CXM","kleb",FALSE)
ecoli_mero<-function_3("e_MER","ecoli",TRUE)
kleb_mero<-function_3("e_MER","kleb",FALSE)



plot_p3<-function(drug){
  fitted<-rbind(get(paste0("ecoli_",drug))$fitted,get(paste0("kleb_",drug))$fitted)
  deriv<-rbind(get(paste0("ecoli_",drug))$deriv,get(paste0("kleb_",drug))$deriv)
p3<- ggplot(fitted) +
  aes(x=date,y=est,color=bug) +
  geom_line() +
  geom_ribbon(data=fitted,aes(x=date,ymin=low,ymax=high,fill=healthcare_exposure),alpha=0.1) +
  ylab("s(timeline)") +
  theme_minimal() +
  scale_x_date(breaks = seq(as.Date("2012-01-01"), as.Date("2024-01-01"), by="2 years"),
               date_labels = "%Y") +
  ggtitle("") +
  xlab("") +
  geom_vline(data=filter(deriv,first_breakpoint==1),aes(xintercept=date,color=bug),linetype="dashed") +
  facet_wrap()

return(p3)
}

coamox_p3<-plot_p3("coamox")
cip_p3<-plot_p3("cip")
mero_p3<-plot_p3("mero")
gent_p3<-plot_p3("gent")
piptaz_p3<-plot_p3("piptaz")
cotrim_p3<-plot_p3("cotrim")
tgc_p3<-plot_p3("tgc")



coamox_p1 + (coamox_p3 / coamox_p2)
cip_p1 + (cip_p3 / cipro_p2)
mero_p1 + (mero_p3 / mero_p2)
gent_p1 + (gent_p3 / gent_p2)
piptaz_p1 + (piptaz_p3 / piptaz_p2)
cotrim_p1 + (cotrim_p3 / cotrim_p2)
tgc + (tgc_p3 / tgc_p2)


####################stratify by healthcare exposure#################

#this function fits models by healthcare exposure and compares smooths to the community smooth as well as returning changepoints by healthcare exposure

data$healthcare_exposure<-as.factor(data$healthcare_exposure)
he<-function(drug,bug,seasonal){


  if(seasonal==TRUE){
    formula <- as.formula(paste(drug, "~ s(timeline,by=healthcare_exposure, bs='ad', m=3, k=25) + s(day, bs='cc', k=12) + healthcare_exposure"))
  }else{
    formula <- as.formula(paste(drug, "~ s(timeline,by=healthcare_exposure, bs='ad', m=3, k=25) + healthcare_exposure"))

  }
  data_<-filter(data,species==bug)
  cat("\nfitting model")
  model<-bam(formula,
             data=data_,
             family = binomial(link = "logit"),
             method = "REML")

  cat("\n Difference smooths")
dif<-difference_smooths(model,select="s(timeline)",n = length(unique(data$timeline)))
dif$sig<-ifelse(dif$.lower_ci < 0 & dif$.upper_ci < 0 ,"Significant", ifelse(dif$.lower_ci >0 & dif$.upper_ci >0,"Significant","Not Significant"))


healthcare<-unique(data$healthcare_exposure)
healthcare<-healthcare[!is.na(healthcare)]
time<-seq(min(data$timeline,na.rm = T),max(data$timeline,na.rm = T),1)
predict_data<-expand_grid(healthcare_exposure=healthcare,timeline=time)
predict_data$date<-min(data$spec_dt,na.rm = T) + days(predict_data$timeline)
predict_data$day<-yday(predict_data$date)

pred_n<-gratia::smooth_estimates(model,data = predict_data,select="s(timeline):healthcare_exposurenosocomial",n=1000) %>% add_confint()
pred_qn<-gratia::smooth_estimates(model,data = predict_data,select="s(timeline):healthcare_exposurequasi-nosocomial",n=1000)%>% add_confint()
pred_qc<-gratia::smooth_estimates(model,data = predict_data,select="s(timeline):healthcare_exposurequasi-community",n=1000)%>% add_confint()
pred_c<-gratia::smooth_estimates(model,data = predict_data,select="s(timeline):healthcare_exposurecommunity",n=1000)%>% add_confint()

pred<-rbind(pred_n,pred_qn,pred_qc,pred_c)

pred$date<-min(data$spec_dt,na.rm = T) + days(pred$timeline)


dif2<-filter(dif,.level_1 =="community" | .level_2 == "community") %>% filter(sig=="Significant")
dif2$time<-min(data$spec_dt,na.rm = T) + days(round(dif2$timeline,0))

dif2$healthcare_exposure<-ifelse(dif2$.level_1!="community",dif2$.level_1,dif2$.level_2)



dif2 <- dif2 %>%
  rowwise() %>%
  mutate(comparison = paste(sort(c(.level_1, .level_2)), collapse = " vs. ")) %>%
  ungroup()

cat("\nFitted samples")
fitted<-fitted_samples(model,data=predict_data,n=1000)
fitted<-fitted %>% group_by(.row) %>%
  summarise(low=quantile(.fitted,probs = c(0.025)),
            high=quantile(.fitted,probs=c(0.975)),
            est=quantile(.fitted,probs=c(0.5)))


pred<-left_join(pred,dif2,by=c("date"="time","healthcare_exposure"="healthcare_exposure"))
pred$sig2 <- ifelse(is.na(pred$sig), 0, ifelse(pred$sig == "Significant", 1, 0))
pred$bug<-bug
fitted$bug<-bug
fitted$drug<-drug
fitted$date<-predict_data$date
fitted$healthcare_exposure<-predict_data$healthcare_exposure

cat("\nmaking derivative samples")
deriv<-derivative_samples(model,order = 1, focal = "timeline",scale = "response",draws = 1000,data = predict_data)
deriv<-deriv %>% group_by(.row) %>%
  summarise(low=quantile(.derivative,probs = c(0.025)),
            high=quantile(.derivative,probs=c(0.975)),
            est=quantile(.derivative,probs=c(0.5)))
deriv$date<-predict_data$date
deriv<-annotate_derivatives(deriv)
deriv$bug<-bug
deriv$healthcare_exposure<-predict_data$healthcare_exposure

return(list(pred=pred,fitted=fitted,deriv=deriv))
}

he_ecoli_coamox<-he("e_AMOXCLAV","ecoli",TRUE)
he_kleb_coamox<-he("e_AMOXCLAV","kleb",TRUE)
he_ecoli_piptaz<-he("e_PIPTAZ","ecoli",TRUE)
he_kleb_piptaz<-he("e_PIPTAZ","kleb",FALSE)
he_ecoli_gent<-he("e_GENT","ecoli",TRUE)
he_kleb_gent<-he("e_GENT","kleb",FALSE)
he_ecoli_tgc<-he("threegc","ecoli",TRUE)
he_kleb_tgc<-he("threegc","kleb",FALSE)
he_ecoli_cotrim<-he("e_CXM","ecoli",TRUE)
he_kleb_cotrim<-he("e_CXM","kleb",FALSE)
he_ecoli_cipro<-he("e_CIP","ecoli",TRUE)
he_kleb_cipro<-he("e_CIP","kleb",FALSE)
he_ecoli_mero<-he("e_MER","ecoli",TRUE)
he_kleb_mero<-he("e_MER","kleb",FALSE)

plot_p3<-function(drug){
  fitted<-rbind(get(paste0("he_ecoli_",drug))$fitted,get(paste0("he_kleb_",drug))$fitted)
  deriv<-rbind(get(paste0("he_ecoli_",drug))$deriv,get(paste0("he_kleb_",drug))$deriv)
  deriv<-left_join(deriv,fitted,by=c("date","bug","healthcare_exposure"))
  pred<-rbind(get(paste0("he_ecoli_",drug))$pred,get(paste0("he_kleb_",drug))$pred)
  pred<-left_join(pred,fitted,by=c("date","bug","healthcare_exposure"))
  pred$direction<-ifelse(pred$.diff <0, "Positive","Negative")
  p3<- ggplot(fitted) +
    aes(x=date,y=est,color=healthcare_exposure) +
    geom_line() +
    geom_ribbon(data=fitted,aes(x=date,ymin=low,ymax=high),alpha=0.1) +
    ylab("Proportion resistant") +
    theme_minimal() + theme(axis.text.x = element_text(angle=-90)) +
    scale_x_date(breaks = seq(as.Date("2012-01-01"), as.Date("2024-01-01"), by="2 years"),
                 date_labels = "%Y") +
    ggtitle("") +
    xlab("") + labs(color="Healthcare Exposure",fill="Healthcare Exposure") +
    facet_wrap(~bug)
  p3<-p3 +new_scale_color()
   p3<- p3 + geom_segment(data=filter(pred,sig2==1),aes(x=date,xend=date + days(10),y=est,color=direction)) +
      scale_color_manual(values=c("black","darkred")) + labs(color="Comparison with community smooth",linetype="Species")
   
  
  params<-ggplot_build(p3)
  params<-params$layout$panel_params[[1]]$y$breaks
  params<-abs(params[3]-params[2])
  #geom_vline(data=filter(deriv,first_breakpoint==1),aes(xintercept=date,color=bug),linetype="dashed",ymax=0.05)
  p3<-p3+ geom_segment(data=filter(deriv,first_breakpoint==1),aes(x=date,xend=(date),y=est.y-(params/4),yend=est.y+(params/4),linetype=bug),linewidth=1,color="darkred") +
    scale_linetype_manual(values=c("dashed","dotted"))
  
    p4<-ggplot(pred) + aes(x=date,y=.estimate,group=healthcare_exposure,color=healthcare_exposure) +geom_line() + theme_minimal() +
    scale_x_date(breaks = seq(as.Date("2012-01-01"), as.Date("2024-01-01"), by="2 years"),
                 date_labels = "%Y") + theme(legend.position = "None",axis.text.x = element_text(angle=-90)) + labs(linetype="Species") + ylab("s(time)")+
    facet_wrap(~bug) 
  
  return(list(p3=p3))
}

(plot_p3("coamox")$p3 + plot_p3("tgc")$p3)  /
  (plot_p3("gent")$p3 + plot_p3("piptaz")$p3)+  plot_annotation(tag_levels  = "A")  + plot_layout(guides="collect")

(plot_p3("cipro")$p3 + plot_p3("cotrim")$p3) +  plot_annotation(tag_levels  = "A")  + plot_layout(guides="collect")

# / (plot_p3("gent")$p3 |plot_p3("gent")$p4) /
  (plot_p3("tgc")$p3 |plot_p3("tgc")$p4) / (plot_p3("cipro")$p3 | plot_p3("cipro")$p4) +
 
  

coamox_pred<-rbind(he_ecoli_coamox,he_kleb_coamox)
plot_p3<-function(drug){
  fitted<-get(paste0("fitted_",drug))
  pred<-get(paste0(drug,"_pred"))
  deriv<-get(paste0("deriv_",drug))
  p3<- ggplot(fitted) +
    aes(x=date,y=est) +
    geom_line() +
    geom_ribbon(data=fitted,aes(x=date,ymin=low,ymax=high),alpha=0.1) +
    ylab("s(timeline)") +
    theme_minimal() +
    scale_x_date(breaks = seq(as.Date("2012-01-01"), as.Date("2024-01-01"), by="1 year"),
                 date_labels = "%Y") +
    ggtitle("") +
    xlab("") +
    geom_vline(data=filter(deriv,first_breakpoint==1),aes(xintercept=date),linetype="dashed",alpha=0.2) +
    geom_line(data=pred,aes(x=date,y=.estimate,color=healthcare_exposure),alpha=0.5) +
    geom_segment(data=filter(pred,sig2==1),aes(x=date,xend=(date + days(1)),y=.estimate,group=healthcare_exposure,color=healthcare_exposure),linewidth=2,alpha=0.5) +
    labs(color="Healthcare Exposure") +
    facet_wrap(~bug)
  return(p3)
}

plot_p3("coamox")


plot_p4<-function(code,bug){
  code <- sym(code)
  data_<-filter(data,species==bug)
prop <- data %>%
  group_by(year, monthno,healthcare_exposure) %>%
  summarise(S=sum(!!code == 0, na.rm = TRUE),
            R=sum(!!code == 1, na.rm = TRUE))

prop$date<-ym(paste0(prop$year,prop$monthno))
prop<-prop %>% pivot_longer(cols = 4:5,names_to = "which",values_to = "n")

prop$healthcare_exposure<-factor(
  prop$healthcare_exposure,levels=c("community","quasi-community","quasi-nosocomial","nosocomial")
)

col<-c("#470000",   "#700404FD", "#AB0707FD", "#FA0808FD" ,"#001F04"  , "#00520B" ,  "#009E18",   "#00E827")
col<-rev(col)
prop$which<-factor(prop$which,levels=c("S","R"))
p4<-ggplot(prop) +
  aes(x=date,y=n,fill=interaction(healthcare_exposure,which)) +
  geom_bar(position="stack", stat="identity") +
  scale_x_date(breaks = seq(as.Date("2012-01-01"), as.Date("2024-01-01"), by="1 year"),
               date_labels = "%Y") +
  scale_fill_manual(values=col) +
  theme_minimal() + labs(fill="") +
  ylab("Number of isolates")
return(p4)
}

plot_p4("e_AMOXCLAV","ecoli")

(p1/p4 + plot_layout(heights=c(3,1),guides = "collect")) | p3

##########ethnos###################


data$ETHNOS<-as.factor(data$ETHNOS)
data$ETHNOS_short<-case_when(
  data$ETHNOS == 'British (White)' ~ 'White British',
  data$ETHNOS ==  'Any other White background' | data$ETHNOS == 'Irish (White)' ~ 'Other White background',
  data$ETHNOS == 'African (Black or Black British)' | data$ETHNOS == 'Caribbean (Black or Black British)' | data$ETHNOS == 'Any other Black background' ~ 'Black British',
  data$ETHNOS == "Indian (Asian or Asian British)" | data$ETHNOS == "Pakistani (Asian or Asian British)" | data$ETHNOS == "Bangladeshi (Asian or Asian British)" | data$ETHNOS == "Any other Asian background" | data$ETHNOS == "Chinese (other ethnic group)" ~ 'Asian',
  data$ETHNOS == "Any other Mixed background" | data$ETHNOS ==  "White and Asian (Mixed)" | data$ETHNOS == "White and Black African (Mixed)"| data$ETHNOS =="White and Black Caribbean (Mixed)" | data$ETHNOS == "Any other ethnic group" ~ "Other",
  data$ETHNOS == "Unknown"  ~ 'Unknown'
)
data$ETHNOS_short<-factor(data$ETHNOS_short,levels=c('White British','Other White background','Black British','Asian','Other','Unknown'))



i_ethnos<-function(drug,bug,seasonal){

  if(seasonal==TRUE){
    formula <- as.formula(paste(drug, "~ s(timeline, by=ETHNOS_short, bs='ad', m=3, k=25) + s(day, bs='cc', k=12) +
                                ETHNOS_short"))
  }else{
    formula <- as.formula(paste(drug, "~ s(timeline,by=healthcare_exposure) + healthcare_exposure, bs='ad', m=3, k=25)"))

  }
  data_<-filter(data,species==bug)
  cat("\nfitting model")
  model<-bam(formula,
             data=data_,
             family = binomial(link = "logit"),
             method = "REML")


cat("\nDifference smooths")

dif<-difference_smooths(model,select="s(timeline)",n = length(unique(data$timeline)))
dif$sig<-ifelse(dif$.lower_ci < 0 & dif$.upper_ci < 0 ,"Significant", ifelse(dif$.lower_ci >0 & dif$.upper_ci >0,"Significant","Not Significant"))

cat("\nMaking test dataset")
ethnos<-unique(data$ETHNOS_short)
time<-seq(min(data$timeline,na.rm = T),max(data$timeline,na.rm = T),1)
predict_data<-expand_grid(ETHNOS_short=ethnos,timeline=time)
predict_data$date<-min(data$spec_dt,na.rm = T) + days(predict_data$timeline)
predict_data$day<-yday(predict_data$date)

cat("\nTaking smooth estimates")
pred_wb<-gratia::smooth_estimates(model,data = predict_data,select="s(timeline):ETHNOS_shortWhite British",n=1000) %>% add_confint()
pred_ow<-gratia::smooth_estimates(model,data = predict_data,select="s(timeline):ETHNOS_shortOther White background",n=1000) %>% add_confint()
pred_b<-gratia::smooth_estimates(model,data = predict_data,select="s(timeline):ETHNOS_shortBlack British",n=1000)%>% add_confint()
pred_a<-gratia::smooth_estimates(model,data = predict_data,select="s(timeline):ETHNOS_shortAsian",n=1000)%>% add_confint()
pred_m<-gratia::smooth_estimates(model,data = predict_data,select="s(timeline):ETHNOS_shortMixed",n=1000)%>% add_confint()
pred_o<-gratia::smooth_estimates(model,data = predict_data,select="s(timeline):ETHNOS_shortOther ethnic groups",n=1000)%>% add_confint()
pred_u<-gratia::smooth_estimates(model,data = predict_data,select="s(timeline):ETHNOS_shortUnknown",n=1000)%>% add_confint()

pred<-rbind(pred_wb,pred_ow,pred_b,pred_a,pred_m,pred_o,pred_u)

pred$date<-min(data$spec_dt,na.rm = T) + days(pred$timeline)


dif2<-filter(dif,.level_1 =="White British" | .level_2 == "White British") %>% filter(sig=="Significant")
dif2$time<-min(data$spec_dt,na.rm = T) + days(round(dif2$timeline,0))

dif2$ETHNOS_short<-ifelse(dif2$.level_1!="White British",dif2$.level_1,dif2$.level_2)



dif2 <- dif2 %>%
  rowwise() %>%
  mutate(comparison = paste(sort(c(.level_1, .level_2)), collapse = " vs. ")) %>%
  ungroup()


pred<-left_join(pred,dif2,by=c("date"="time","ETHNOS_short"="ETHNOS_short"))
pred$sig2 <- ifelse(is.na(pred$sig), 0, ifelse(pred$sig == "Significant", 1, 0))
pred$bug<-bug
return(pred)
}

ecoli_coamox_ethnos<-i_ethnos("e_AMOXCLAV","ecoli",TRUE)
kleb_coamox_ethnos<-i_ethnos("e_AMOXCLAV","kleb",TRUE)

coamox_ethnos<-rbind(ecoli_coamox_ethnos,kleb_coamox_ethnos)

plot_p4<-function(drug){
  pred<-get(paste0(drug,'_ethnos'))
p4<-ggplot(pred) +
  aes(x=date,y=.estimate,color=ETHNOS_short) +
  geom_line(aes(size=sig2))  +
  scale_size_continuous(range = c(0.5, 1.5)) + theme_minimal() +
  guides(size="none") +
  facet_wrap(~bug)
return(p4)
}

plot_p4("coamox")

ethnos_fitted<-function(drug,bug,seasonal){
  if(seasonal==TRUE){
    formula <- as.formula(paste(drug, "~ s(timeline, by=ETHNOS_short, bs='ad', m=3, k=25) + s(day, bs='cc', k=12) +
                                ETHNOS_short"))
  }else{
    formula <- as.formula(paste(drug, "~ s(timeline,by=healthcare_exposure) + healthcare_exposure, bs='ad', m=3, k=25)"))

  }
data_<-filter(data,species==bug)
cat("\nfitting model")
model<-bam(formula,
           data=data_,
           family = binomial(link = "logit"),
           method = "REML")

cat("\nMaking test dataset")
ethnos<-unique(data$ETHNOS_short)
time<-seq(min(data$timeline,na.rm = T),max(data$timeline,na.rm = T),1)
predict_data<-expand_grid(ETHNOS_short=ethnos,timeline=time)
predict_data$date<-min(data$spec_dt,na.rm = T) + days(predict_data$timeline)
predict_data$day<-yday(predict_data$date)

cat("\nFitted Samples")
fitted<-fitted_samples(model,data = predict_data,n = 1000)
fitted<-fitted %>%
  group_by(.row) %>%
  summarise(low=quantile(.fitted,probs=c(0.025)),
            est=quantile(.fitted,probs=c(0.5)),
            high=quantile(.fitted,probs=c(0.975)))

fitted$date<-predict_data$date
fitted$ETHNOS_short<-predict_data$ETHNOS_short
fitted$bug<-bug
return(fitted)
}
ecoli_coamox_ethnos_fitted<-ethnos_fitted("e_AMOXCLAV","ecoli",TRUE)
kleb_coamox_ethnos_fitted<-ethnos_fitted("e_AMOXCLAV","kleb",TRUE)

coamox_fitted<-rbind(ecoli_coamox_ethnos_fitted,kleb_coamox_ethnos_fitted)

plot_ethnos_fitted<-function(drug){

  fitted<-get(paste0(drug,'_fitted'))

p3<- ggplot(fitted) +
  aes(x=date,y=est,color=ETHNOS_short) +
  geom_line() +
  geom_ribbon(data=fitted,aes(x=date,ymin=low,ymax=high),alpha=0.1) +
  ylab("Proportion of resistance isolates") +
  theme_minimal() + facet_wrap(~ETHNOS_short~bug) +
  scale_x_date(breaks = seq(as.Date("2012-01-01"), as.Date("2024-01-01"), by="2 years"),
               date_labels = "%Y") +
  ggtitle(drug) +
  xlab("")

return(p3)
}

plot_ethnos_fitted("coamox")

(p1/p2 + plot_layout(heights=c(3,1),guides = "collect")) | p3

prop <- data %>%
  group_by(year, monthno,healthcare_exposure) %>%
  summarise(S=sum(e_AMOXCLAV == 0, na.rm = TRUE),
            R=sum(e_AMOXCLAV == 1, na.rm = TRUE))

prop$date<-ym(paste0(prop$year,prop$monthno))
prop<-prop %>% pivot_longer(cols = 4:5,names_to = "which",values_to = "n")

prop$healthcare_exposure<-factor(
  prop$healthcare_exposure,levels=c("community","quasi-community","quasi-nosocomial","nosocomial")
)

col<-c("#470000",   "#700404FD", "#AB0707FD", "#FA0808FD" ,"#001F04"  , "#00520B" ,  "#009E18",   "#00E827")
col<-rev(col)
prop$which<-factor(prop$which,levels=c("S","R"))
p4<-ggplot(prop) +
  aes(x=date,y=n,fill=interaction(healthcare_exposure,which)) +
  geom_bar(position="stack", stat="identity") +
  scale_x_date(breaks = seq(as.Date("2012-01-01"), as.Date("2024-01-01"), by="1 year"),
               date_labels = "%Y") +
  scale_fill_manual(values=col) +
  theme_minimal() + labs(fill="") +
  ylab("Number of isolates")

(p1/p4 + plot_layout(heights=c(3,1),guides = "collect")) | p3

######tempero-spatial###############

#this fits national data with an interaction between space (NUTS region) and time

d$nuts_name<-factor(d$nuts_name,levels=names(nb))
d$day<-yday(d$spec_dt)

m1 <- bam(e_AMOXCLAV ~ s(nuts_name, bs = 'mrf', xt = list(nb = nb)) +
            s(timeline,k=25,bs="tp")  + s(day,k=12,bs="cc") +
            ti(nuts_name, timeline, bs = c("mrf","tp"),
               xt = list(nb = nb)),
          data = d,
          method = 'REML',
          family = binomial(link="logit"))

nuts<-unique(d$nuts_name)
time<-seq(0,max(d$timeline),28)
new_data<-expand_grid(nuts_name=nuts,timeline=time)
new_data$date<-min(d$spec_dt,na.rm = T) + days(new_data$timeline)
new_data$day<-yday(new_data$date)

m<-fitted_values(m1,data = new_data)


p1a<- ggplot(fitted_national_estimate) +
  aes(x=date,y=est) +
  geom_line(color="blue") +
  geom_ribbon(data=fitted_national_estimate,aes(x=date,ymin=low,ymax=high),alpha=0.1) +
  ylab("Proportion of isolates co-amoxiclav resistant") +
  theme_minimal() +
  scale_x_date(breaks = seq(as.Date("2012-01-01"), as.Date("2024-01-01"), by="1 year"),
               date_labels = "%Y") +
  ggtitle("") +
  xlab("") +
  geom_line(data=m,aes(x=date,y=.fitted,color=nuts_name),alpha=0.5) + scale_color_viridis_d() +
  labs(color="Region")

coamox_plot1<-(p1a/p4 + plot_layout(heights=c(3,1),guides = "collect")) | p3
coamox_plot1 + plot_annotation(title = "Co-amoxiclav - E. coli")

samples<-smooth_samples(m1,select = "ti(nuts_name,timeline)",n = 1000,data = new_data)
samples<-samples %>% group_by(.row,nuts_name) %>%
  summarise(low=quantile(.value,probs = c(0.025)),
            high=quantile(.value,probs=c(0.975)),
            est=quantile(.value,probs=c(0.5)))

samples$date<-new_data$date

spatial_interactions<-ggplot(samples) +
  aes(x=date,y=est) +
  geom_line() +
  geom_ribbon(data=samples,aes(x=date,ymin=low,ymax=high),alpha=0.1) +
  facet_wrap(~nuts_name) + theme_minimal() +
  scale_x_date(breaks = seq(as.Date("2012-01-01"), as.Date("2024-01-01"), by="1 year"),
               date_labels = "%Y") +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("ti(region,time)") +
  xlab("Date")

spatial_fitted<-ggplot(data=m) +
  aes(x=date,y=.fitted) +
  geom_line() +
  geom_ribbon(data=m,aes(x=date,ymin=.lower_ci,ymax=.upper_ci),alpha=0.2) +
  facet_wrap(~nuts_name) +
  theme_minimal() +
  scale_x_date(breaks = seq(as.Date("2012-01-01"), as.Date("2024-01-01"), by="2 year"),
               date_labels = "%Y") +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Estimated proportion resistant") +
  xlab("Date")

coamox_space_time<-spatial_interactions + spatial_fitted
coamox_space_time + plot_annotation("E. coli - co-amoxiclav")

############variable-time interactions#########
########continuous########
#this is a general function for fitting modules for continuous variables with and without an interaction for time 

data$age<-ifelse(data$age > 100,100,data$age)

#first fit a model with no interaction with time
main_continuous<-function(drug,bug,seasonal,variable,min,max,step,freedom){
  if(seasonal==TRUE){
    formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25) + s(day, bs='cc', k=12) +
                                s(",variable,",bs='tp',k=",freedom,")"))
  }else{
    formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25) +
                                s(",variable,",bs='tp',k=",freedom,")"))

  }
  data_<-filter(data,species==bug)
  cat("\nfitting model")
  model<-bam(formula,
             data=data_,
             family = binomial(link = "logit"),
             method = "REML")

  cat("\nMaking dummy data")
  times<-seq(0,max(data$timeline,na.rm=T),14)
  variable_<-seq(min,max,step)
  predict_data<-expand_grid(timeline=times,variable=variable_)
  predict_data$date<-min(data$spec_dt,na.rm=T) + days(predict_data$timeline)
  predict_data$day<-yday(predict_data$date)
  predict_data<-filter(predict_data,timeline==4018) #no int with
  y=sym(variable)
  predict_data<-rename(predict_data, !!y := variable)


  x=paste0("s(",variable,")")
  cat("\nTaking smooth samples")
  #samples<-smooth_samples(model,select = x,n = 1000,data = predict_data)
  samples<-fitted_samples(model,data = predict_data,n=1000)
  cat("\nGenerating posterior credibility intervals")
samples<-samples %>%
  group_by(.row) %>%
  summarise(low=quantile(.fitted,probs=c(0.025)),
            est=quantile(.fitted,probs=c(0.5)),
            high=quantile(.fitted,probs=c(0.975)))

samples[[variable]]<-predict_data[[variable]]
samples$drug<-drug
samples$bug<-bug

return(samples)
}

ecoli_coamox_age<-main_continuous("e_AMOXCLAV","ecoli",TRUE,"age",0,100,1,25)
kleb_coamox_age<-main_continuous("e_AMOXCLAV","kleb",TRUE,"age",0,100,1,25)
ecoli_gent_age<-main_continuous("e_GENT","ecoli",TRUE,"age",0,100,1,25)
kleb_gent_age<-main_continuous("e_GENT","kleb",FALSE,"age",0,100,1,25)
ecoli_piptaz_age<-main_continuous("e_PIPTAZ","ecoli",TRUE,"age",0,100,1,25)
kleb_piptaz_age<-main_continuous("e_PIPTAZ","kleb",FALSE,"age",0,100,1,25)
ecoli_cip_age<-main_continuous("e_CIP","ecoli",TRUE,"age",0,100,1,25)
kleb_cip_age<-main_continuous("e_CIP","kleb",FALSE,"age",0,100,1,25)
ecoli_mero_age<-main_continuous("e_MER","ecoli",FALSE,"age",0,100,1,25)
kleb_mero_age<-main_continuous("e_MER","kleb",FALSE,"age",0,100,1,25)
ecoli_tgc_age<-main_continuous("threegc","ecoli",TRUE,"age",0,100,1,25)
kleb_tgc_age<-main_continuous("threegc","kleb",FALSE,"age",0,100,1,25)
ecoli_cotrim_age<-main_continuous("e_CXM","ecoli",TRUE,"age",0,100,1,25)
kleb_cotrim_age<-main_continuous("e_CXM","kleb",FALSE,"age",0,100,1,25)

ecoli_coamox_imd<-main_continuous("e_AMOXCLAV","ecoli",TRUE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),1,25)
kleb_coamox_imd<-main_continuous("e_AMOXCLAV","kleb",TRUE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),1,25)
ecoli_gent_imd<-main_continuous("e_GENT","ecoli",TRUE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),1,25)
kleb_gent_imd<-main_continuous("e_GENT","kleb",FALSE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),1,25)
ecoli_piptaz_imd<-main_continuous("e_PIPTAZ","ecoli",TRUE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),1,25)
kleb_piptaz_imd<-main_continuous("e_PIPTAZ","kleb",FALSE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),1,25)
ecoli_cip_imd<-main_continuous("e_CIP","ecoli",TRUE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),1,25)
kleb_cip_imd<-main_continuous("e_CIP","kleb",FALSE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),1,25)
#ecoli_mero_imd<-main_continuous("e_MER","ecoli",FALSE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),1)
#kleb_mero_imd<-main_continuous("e_MER","kleb",FALSE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),1)
ecoli_tgc_imd<-main_continuous("threegc","ecoli",TRUE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),1,25)
kleb_tgc_imd<-main_continuous("threegc","kleb",FALSE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),1,25)
ecoli_cotrim_imd<-main_continuous("e_CXM","ecoli",TRUE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),1,25)
kleb_cotrim_imd<-main_continuous("e_CXM","kleb",FALSE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),1,25)

data$n_admissions<-ifelse(data$n_admissions>10,10,data$n_admissions)
ecoli_coamox_n_admissions<-main_continuous("e_AMOXCLAV","ecoli",TRUE,"n_admissions",0,10,1,6)
kleb_coamox_n_admissions<-main_continuous("e_AMOXCLAV","kleb",TRUE,"n_admissions",0,10,1,6)
ecoli_gent_n_admissions<-main_continuous("e_GENT","ecoli",TRUE,"n_admissions",0,10,1,6)
kleb_gent_n_admissions<-main_continuous("e_GENT","kleb",FALSE,"n_admissions",0,10,1,6)
ecoli_piptaz_n_admissions<-main_continuous("e_PIPTAZ","ecoli",TRUE,"n_admissions",0,10,1,6)
kleb_piptaz_n_admissions<-main_continuous("e_PIPTAZ","kleb",FALSE,"n_admissions",0,10,1,6)
ecoli_cip_n_admissions<-main_continuous("e_CIP","ecoli",TRUE,"n_admissions",0,10,1,6)
kleb_cip_n_admissions<-main_continuous("e_CIP","kleb",FALSE,"n_admissions",0,10,1,6)
#ecoli_mero_n_admissions<-main_continuous("e_MER","ecoli",FALSE,"n_admissions",0,10,1,5)
#kleb_mero_n_admissions<-main_continuous("e_MER","kleb",FALSE,"n_admissions",0,10,1,5)
ecoli_tgc_n_admissions<-main_continuous("threegc","ecoli",TRUE,"n_admissions",0,10,1,6)
kleb_tgc_n_admissions<-main_continuous("threegc","kleb",FALSE,"n_admissions",0,10,1,6)
ecoli_cotrim_n_admissions<-main_continuous("e_CXM","ecoli",TRUE,"n_admissions",0,10,1,6)
kleb_cotrim_n_admissions<-main_continuous("e_CXM","kleb",FALSE,"n_admissions",0,10,1,6)

ecoli_coamox_elixhauser<-main_continuous("e_AMOXCLAV","ecoli",TRUE,"elixhauser",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.975),na.rm = T),1,10)
kleb_coamox_elixhauser<-main_continuous("e_AMOXCLAV","kleb",TRUE,"elixhauser",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.975),na.rm = T),1,10)
ecoli_gent_elixhauser<-main_continuous("e_GENT","ecoli",TRUE,"elixhauser",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.975),na.rm = T),1,10)
kleb_gent_elixhauser<-main_continuous("e_GENT","kleb",FALSE,"elixhauser",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.975),na.rm = T),1,10)
ecoli_piptaz_elixhauser<-main_continuous("e_PIPTAZ","ecoli",TRUE,"elixhauser",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.975),na.rm = T),1,10)
kleb_piptaz_elixhauser<-main_continuous("e_PIPTAZ","kleb",FALSE,"elixhauser",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.975),na.rm = T),1,10)
ecoli_cip_elixhauser<-main_continuous("e_CIP","ecoli",TRUE,"elixhauser",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.975),na.rm = T),1,10)
kleb_cip_elixhauser<-main_continuous("e_CIP","kleb",FALSE,"elixhauser",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.975),na.rm = T),1,10)
#ecoli_mero_elixhauser<-main_continuous("e_MER","ecoli",FALSE,"elixhauser",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),1,10)
#kleb_mero_elixhauser<-main_continuous("e_MER","kleb",FALSE,"elixhauser",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),1,10)
ecoli_tgc_elixhauser<-main_continuous("threegc","ecoli",TRUE,"elixhauser",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.975),na.rm = T),1,10)
kleb_tgc_elixhauser<-main_continuous("threegc","kleb",FALSE,"elixhauser",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.975),na.rm = T),1,10)
ecoli_cotrim_elixhauser<-main_continuous("e_CXM","ecoli",TRUE,"elixhauser",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.975),na.rm = T),1,10)
kleb_cotrim_elixhauser<-main_continuous("e_CXM","kleb",FALSE,"elixhauser",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.975),na.rm = T),1,10)


age_data_frames <- ls(pattern = "*_*_age")
age <- map_dfr(age_data_frames, get)

imd_data_frames <- ls(pattern = "*_*_imd")
imd_score <- map_dfr(imd_data_frames, get)

admi_data_frames<-ls(pattern = "*_*_n_admissions")
n_admissions<-map_dfr(admi_data_frames,get)

elix_data_frames<-ls(pattern = "*_elixhauser")
elixhauser<-map_dfr(elix_data_frames,get)

plot_main_continuous<-function(variable,labs,variable_clean){
  fitted<-get(paste0(variable))
  fitted$drug<- case_when(
    fitted$drug == "e_AMOXCLAV" ~ "co-amoxiclav",
    fitted$drug =="e_CIP" ~ "ciprofloxacin",
    fitted$drug == "e_CXM" ~ "co-trimoxazole",
    fitted$drug == "e_GENT" ~ "gentamicin",
    fitted$drug == "e_MER" ~ "meropenem",
    fitted$drug == "e_PIPTAZ" ~ "piperacillin-tazobactam",
    fitted$drug == "threegc" ~ "third generation cephalosporin"
  )
  p<-ggplot(fitted) +
    aes_string(x=variable,y="est",color="bug") +
    geom_line() +
    geom_ribbon(data=fitted,aes_string(x=variable,ymin="high",ymax="low",fill="bug"),alpha=0.1) +
    theme_minimal() + theme(legend.position="None") +
    ylab(ifelse(labs==TRUE,"Estimated proportion resistant","")) +
    xlab(variable_clean) +
    guides(fill="none") +
    facet_wrap(~drug,ncol = 1) +
    labs(color="Species",fill="Species")
  return(p)
}

plot_main_continuous("age")+ plot_main_continuous("imd_score") +plot_main_continuous("n_admissions") + plot_main_continuous("elixhauser") + 
  plot_layout(nrow = 1)

time_cut<-function(drug,bug,seasonal,variable,min,max,step){
  if(seasonal==TRUE){
    formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25) + s(day, bs='cc', k=12) +
                                s(",variable,",bs='tp') +
                                ti(",variable,",timeline)"))
  }else{
    formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25)  +
                                s(",variable,",bs='tp') +
                                    ti(",variable,",timeline)"))
    
  }
  data_<-filter(data,species==bug)
  cat("\nfitting model")
  model<-bam(formula,
             data=data_,
             family = binomial(link = "logit"),
             method = "REML")
  
  cat("\nMaking dummy data")
  times<-seq(0,max(data$timeline,na.rm=T),1)
  variable_<-seq(min,max,step)
  predict_data<-expand_grid(timeline=times,variable=variable_)
  predict_data$date<-min(data$spec_dt,na.rm=T) + days(predict_data$timeline)
  predict_data$day<-yday(predict_data$date)
  predict_data<-filter(predict_data,timeline==0 | timeline==2920 | timeline ==4300) #no int with
  y=sym(variable)
  predict_data<-rename(predict_data, !!y := variable)
  
  
  x=paste0("s(",variable,")")
  cat("\nTaking smooth samples")
  #samples<-smooth_samples(model,select = x,n = 1000,data = predict_data)
  samples<-fitted_samples(model,data = predict_data,n=1000)
  cat("\nGenerating posterior credibility intervals")
  samples<-samples %>%
    group_by(.row) %>%
    summarise(low=quantile(.fitted,probs=c(0.025)),
              est=quantile(.fitted,probs=c(0.5)),
              high=quantile(.fitted,probs=c(0.975)))
  
  samples[[variable]]<-predict_data[[variable]]
  samples$drug<-drug
  samples$bug<-bug
  samples$date<-predict_data$date
  
  return(samples)
}

age_cut<-time_cut("e_AMOXCLAV","ecoli",TRUE,"imd_score",0,40,1)

age_cut$date<-as.character(age_cut$date)
a<-ggplot(age_cut) +
  aes(x=imd_score,y=est,group=date,color=date,fill=date) +
  geom_line() +
  geom_ribbon(aes(x=imd_score,ymin=low,ymax=high),alpha=0.1) + theme_minimal() +
  ylab("Estimated proportion of resistance") +
  xlab("IMD Score") +
  labs(color="Date",fill="Date") +
  ggtitle("Co-amox / IMD")

age_cut<-time_cut("e_AMOXCLAV","ecoli",TRUE,"elixhauser",-6,46,1)

age_cut$date<-as.character(age_cut$date)
b<-ggplot(age_cut) +
  aes(x=elixhauser,y=est,group=date,color=date,fill=date) +
  geom_line() +
  geom_ribbon(aes(x=elixhauser,ymin=low,ymax=high),alpha=0.1) + theme_minimal() +
  ylab("Estimated proportion of resistance") +
  xlab("Elixhauser Score") +
  labs(color="Date",fill="Date") +
  ggtitle("Co-amox / Elixhauser")

age_cut<-time_cut("threegc","ecoli",TRUE,"imd_score",0,40,1)

age_cut$date<-as.character(age_cut$date)
c<-ggplot(age_cut) +
  aes(x=imd_score,y=est,group=date,color=date,fill=date) +
  geom_line() +
  geom_ribbon(aes(x=imd_score,ymin=low,ymax=high),alpha=0.1) + theme_minimal() +
  ylab("Estimated proportion of resistance") +
  xlab("IMD Score") +
  labs(color="Date",fill="Date") +
  ggtitle("Third-gen ceph / IMD")

age_cut<-time_cut("e_CIP","ecoli",TRUE,"age",0,100,1)

age_cut$date<-as.character(age_cut$date)
d<-ggplot(age_cut) +
  aes(x=age,y=est,group=date,color=date,fill=date) +
  geom_line() +
  geom_ribbon(aes(x=age,ymin=low,ymax=high),alpha=0.1) + theme_minimal() +
  ylab("Estimated proportion of resistance") +
  xlab("Age") +
  labs(color="Date",fill="Date") +
  ggtitle("Cipro / Age")

a + b + c + d + plot_layout(guides="collect")

#here we fit an interaction between time and a continuous variable

i_continuous<-function(drug,bug,seasonal,variable,min,max,step,freedom){
  if(seasonal==TRUE){
    formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25) + s(day, bs='cc', k=12) +
                                s(",variable,",bs='tp') +
                                ti(",variable,",timeline)"))
  }else{
    formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25)  +
                                s(",variable,",bs='tp') +
                                    ti(",variable,",timeline)"))

  }
  data_<-filter(data,species==bug)
  cat("\nfitting model")
  model<-bam(formula,
             data=data_,
             family = binomial(link = "logit"),
             method = "REML")

cat("\nmaking dummy data")
times<-seq(0,max(data$timeline,na.rm=T),14)

variable_<-seq(min,max,step)
predict_data<-expand_grid(timeline=times,variable=variable_)
predict_data$date<-min(data$spec_dt,na.rm=T) + days(predict_data$timeline)
predict_data$day<-yday(predict_data$date)
y=sym(variable)
predict_data<-rename(predict_data, !!y := variable)
x=paste0("ti(",variable,",timeline)")
cat("\nTaking smooth samples")
samples<-smooth_samples(model,select = x,n = 1000,data = predict_data)

cat("\nCredibility intervals")
samples<-samples %>%
  group_by(.row) %>%
  summarise(low=quantile(.value,probs=c(0.025)),
            est=quantile(.value,probs=c(0.5)),
            high=quantile(.value,probs=c(0.975)))
samples$date<-predict_data$date
samples[[variable]]<-predict_data[[variable]]
samples$drug<-drug
samples$bug<-bug

cat("\nFitted samples")
fitted<-fitted_samples(model = model,n = 1000,data = predict_data)
fitted<-fitted %>% 
  group_by(.row) %>%
  summarise(low=quantile(.fitted,probs=c(0.025)),
            est=quantile(.fitted,probs=c(0.5)),
            high=quantile(.fitted,probs=c(0.975)))
fitted$date<-predict_data$date
fitted[[variable]]<-predict_data[[variable]]
fitted$bug<-bug
fitted$drug<-drug
cat("\nDone!")
return(list(partial_effects=samples,fitted=fitted))
}

ecoli_coamox_age_i<-i_continuous("e_AMOXCLAV","ecoli",TRUE,"age",quantile(data$age,probs = c(0.005),na.rm = T),quantile(data$age,probs = c(0.995),na.rm = T),15)
kleb_coamox_age_i<-i_continuous("e_AMOXCLAV","kleb",TRUE,"age",quantile(data$age,probs = c(0.005),na.rm = T),quantile(data$age,probs = c(0.995),na.rm = T),15)
ecoli_gent_age_i<-i_continuous("e_GENT","ecoli",TRUE,"age",quantile(data$age,probs = c(0.005),na.rm = T),quantile(data$age,probs = c(0.995),na.rm = T),15)
kleb_gent_age_i<-i_continuous("e_GENT","kleb",FALSE,"age",quantile(data$age,probs = c(0.005),na.rm = T),quantile(data$age,probs = c(0.995),na.rm = T),15)
ecoli_piptaz_age_i<-i_continuous("e_PIPTAZ","ecoli",TRUE,"age",quantile(data$age,probs = c(0.005),na.rm = T),quantile(data$age,probs = c(0.995),na.rm = T),15)
kleb_piptaz_age_i<-i_continuous("e_PIPTAZ","kleb",FALSE,"age",quantile(data$age,probs = c(0.005),na.rm = T),quantile(data$age,probs = c(0.995),na.rm = T),15)
ecoli_cip_age_i<-i_continuous("e_CIP","ecoli",TRUE,"age",quantile(data$age,probs = c(0.005),na.rm = T),quantile(data$age,probs = c(0.995),na.rm = T),15)
kleb_cip_age_i<-i_continuous("e_CIP","kleb",FALSE,"age",quantile(data$age,probs = c(0.005),na.rm = T),quantile(data$age,probs = c(0.995),na.rm = T),15)
#ecoli_mero_age_i<-i_continuous("e_MER","ecoli",FALSE,"age",quantile(data$age,probs = c(0.005),na.rm = T),quantile(data$age,probs = c(0.995),na.rm = T),10)
#kleb_mero_age_i<-i_continuous("e_MER","kleb",FALSE,"age",quantile(data$age,probs = c(0.005),na.rm = T),quantile(data$age,probs = c(0.995),na.rm = T),10)
ecoli_tgc_age_i<-i_continuous("threegc","ecoli",TRUE,"age",quantile(data$age,probs = c(0.005),na.rm = T),quantile(data$age,probs = c(0.995),na.rm = T),15)
kleb_tgc_age_i<-i_continuous("threegc","kleb",FALSE,"age",quantile(data$age,probs = c(0.005),na.rm = T),quantile(data$age,probs = c(0.995),na.rm = T),15)
ecoli_cotrim_age_i<-i_continuous("e_CXM","ecoli",TRUE,"age",quantile(data$age,probs = c(0.005),na.rm = T),quantile(data$age,probs = c(0.995),na.rm = T),15)
kleb_cotrim_age_i<-i_continuous("e_CXM","kleb",FALSE,"age",quantile(data$age,probs = c(0.005),na.rm = T),quantile(data$age,probs = c(0.995),na.rm = T),15)

ecoli_coamox_imd_i<-i_continuous("e_AMOXCLAV","ecoli",TRUE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),15)
kleb_coamox_imd_i<-i_continuous("e_AMOXCLAV","kleb",TRUE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),15)
ecoli_gent_imd_i<-i_continuous("e_GENT","ecoli",TRUE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),15)
kleb_gent_imd_i<-i_continuous("e_GENT","kleb",FALSE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),15)
ecoli_piptaz_imd_i<-i_continuous("e_PIPTAZ","ecoli",TRUE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),15)
kleb_piptaz_imd_i<-i_continuous("e_PIPTAZ","kleb",FALSE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),15)
ecoli_cip_imd_i<-i_continuous("e_CIP","ecoli",TRUE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),15)
kleb_cip_imd_i<-i_continuous("e_CIP","kleb",FALSE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),15)
#ecoli_mero_imd_i<-i_continuous("e_MER","ecoli",FALSE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),10)
#kleb_mero_imd_i<-i_continuous("e_MER","kleb",FALSE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),10)
ecoli_tgc_imd_i<-i_continuous("threegc","ecoli",TRUE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),15)
kleb_tgc_imd_i<-i_continuous("threegc","kleb",FALSE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),15)
ecoli_cotrim_imd_i<-i_continuous("e_CXM","ecoli",TRUE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),15)
kleb_cotrim_imd_i<-i_continuous("e_CXM","kleb",FALSE,"imd_score",quantile(data$imd_score,probs = c(0.005),na.rm = T),quantile(data$imd_score,probs = c(0.995),na.rm = T),15)

ecoli_coamox_n_admissions_i<-i_continuous("e_AMOXCLAV","ecoli",TRUE,"n_admissions",quantile(data$n_admissions,probs = c(0.005),na.rm = T),quantile(data$n_admissions,probs = c(0.995),na.rm = T),1)
kleb_coamox_n_admissions_i<-i_continuous("e_AMOXCLAV","kleb",TRUE,"n_admissions",quantile(data$n_admissions,probs = c(0.005),na.rm = T),quantile(data$n_admissions,probs = c(0.995),na.rm = T),1)
ecoli_gent_n_admissions_i<-i_continuous("e_GENT","ecoli",TRUE,"n_admissions",quantile(data$n_admissions,probs = c(0.005),na.rm = T),quantile(data$n_admissions,probs = c(0.995),na.rm = T),1)
kleb_gent_n_admissions_i<-i_continuous("e_GENT","kleb",FALSE,"n_admissions",quantile(data$n_admissions,probs = c(0.005),na.rm = T),quantile(data$n_admissions,probs = c(0.995),na.rm = T),1)
ecoli_piptaz_n_admissions_i<-i_continuous("e_PIPTAZ","ecoli",TRUE,"n_admissions",quantile(data$n_admissions,probs = c(0.005),na.rm = T),quantile(data$n_admissions,probs = c(0.995),na.rm = T),1)
kleb_piptaz_n_admissions_i<-i_continuous("e_PIPTAZ","kleb",FALSE,"n_admissions",quantile(data$n_admissions,probs = c(0.005),na.rm = T),quantile(data$n_admissions,probs = c(0.995),na.rm = T),1)
ecoli_cip_n_admissions_i<-i_continuous("e_CIP","ecoli",TRUE,"n_admissions",quantile(data$n_admissions,probs = c(0.005),na.rm = T),quantile(data$n_admissions,probs = c(0.995),na.rm = T),1)
kleb_cip_n_admissions_i<-i_continuous("e_CIP","kleb",FALSE,"n_admissions",quantile(data$n_admissions,probs = c(0.005),na.rm = T),quantile(data$n_admissions,probs = c(0.995),na.rm = T),1)
#ecoli_mero_n_admissions_i<-i_continuous("e_MER","ecoli",FALSE,"n_admissions",quantile(data$n_admissions,probs = c(0.005),na.rm = T),quantile(data$n_admissions,probs = c(0.995),na.rm = T),10)
#kleb_mero_n_admissions_i<-i_continuous("e_MER","kleb",FALSE,"n_admissions",quantile(data$n_admissions,probs = c(0.005),na.rm = T),quantile(data$n_admissions,probs = c(0.995),na.rm = T),10)
ecoli_tgc_n_admissions_i<-i_continuous("threegc","ecoli",TRUE,"n_admissions",quantile(data$n_admissions,probs = c(0.005),na.rm = T),quantile(data$n_admissions,probs = c(0.995),na.rm = T),1)
kleb_tgc_n_admissions_i<-i_continuous("threegc","kleb",FALSE,"n_admissions",quantile(data$n_admissions,probs = c(0.005),na.rm = T),quantile(data$n_admissions,probs = c(0.995),na.rm = T),1)
ecoli_cotrim_n_admissions_i<-i_continuous("e_CXM","ecoli",TRUE,"n_admissions",quantile(data$n_admissions,probs = c(0.005),na.rm = T),quantile(data$n_admissions,probs = c(0.995),na.rm = T),1)
kleb_cotrim_n_admissions_i<-i_continuous("e_CXM","kleb",FALSE,"n_admissions",quantile(data$n_admissions,probs = c(0.005),na.rm = T),quantile(data$n_admissions,probs = c(0.995),na.rm = T),1)

ecoli_coamox_elixhauser_i<-i_continuous("e_AMOXCLAV","ecoli",TRUE,"elixhauser",quantile(data$elixhauser,probs = c(0.005),na.rm = T),quantile(data$elixhauser,probs = c(0.995),na.rm = T),15)
kleb_coamox_elixhauser_i<-i_continuous("e_AMOXCLAV","kleb",TRUE,"elixhauser",quantile(data$elixhauser,probs = c(0.005),na.rm = T),quantile(data$elixhauser,probs = c(0.995),na.rm = T),15)
ecoli_gent_elixhauser_i<-i_continuous("e_GENT","ecoli",TRUE,"elixhauser",quantile(data$elixhauser,probs = c(0.005),na.rm = T),quantile(data$elixhauser,probs = c(0.995),na.rm = T),15)
kleb_gent_elixhauser_i<-i_continuous("e_GENT","kleb",FALSE,"elixhauser",quantile(data$elixhauser,probs = c(0.005),na.rm = T),quantile(data$elixhauser,probs = c(0.995),na.rm = T),15)
ecoli_piptaz_elixhauser_i<-i_continuous("e_PIPTAZ","ecoli",TRUE,"elixhauser",quantile(data$elixhauser,probs = c(0.005),na.rm = T),quantile(data$elixhauser,probs = c(0.995),na.rm = T),15)
kleb_piptaz_elixhauser_i<-i_continuous("e_PIPTAZ","kleb",FALSE,"elixhauser",quantile(data$elixhauser,probs = c(0.005),na.rm = T),quantile(data$elixhauser,probs = c(0.995),na.rm = T),15)
ecoli_cip_elixhauser_i<-i_continuous("e_CIP","ecoli",TRUE,"elixhauser",quantile(data$elixhauser,probs = c(0.005),na.rm = T),quantile(data$elixhauser,probs = c(0.995),na.rm = T),15)
kleb_cip_elixhauser_i<-i_continuous("e_CIP","kleb",FALSE,"elixhauser",quantile(data$elixhauser,probs = c(0.005),na.rm = T),quantile(data$elixhauser,probs = c(0.995),na.rm = T),15)
#ecoli_mero_elixhauser_i<-i_continuous("e_MER","ecoli",FALSE,"elixhauser",quantile(data$elixhauser,probs = c(0.005),na.rm = T),quantile(data$elixhauser,probs = c(0.995),na.rm = T),10)
#kleb_mero_elixhauser_i<-i_continuous("e_MER","kleb",FALSE,"elixhauser",quantile(data$elixhauser,probs = c(0.005),na.rm = T),quantile(data$elixhauser,probs = c(0.995),na.rm = T),10)
ecoli_tgc_elixhauser_i<-i_continuous("threegc","ecoli",TRUE,"elixhauser",quantile(data$elixhauser,probs = c(0.005),na.rm = T),quantile(data$elixhauser,probs = c(0.995),na.rm = T),15)
kleb_tgc_elixhauser_i<-i_continuous("threegc","kleb",FALSE,"elixhauser",quantile(data$elixhauser,probs = c(0.005),na.rm = T),quantile(data$elixhauser,probs = c(0.995),na.rm = T),15)
ecoli_cotrim_elixhauser_i<-i_continuous("e_CXM","ecoli",TRUE,"elixhauser",quantile(data$elixhauser,probs = c(0.005),na.rm = T),quantile(data$elixhauser,probs = c(0.995),na.rm = T),15)
kleb_cotrim_elixhauser_i<-i_continuous("e_CXM","kleb",FALSE,"elixhauser",quantile(data$elixhauser,probs = c(0.005),na.rm = T),quantile(data$elixhauser,probs = c(0.995),na.rm = T),15)


plot_interaction_continuous<-function(variable,drug_,bug_,variable_clean,drug_clean,labs){
  fitted<-get(paste0(variable))$fitted
  
  fitted<-filter(fitted,drug==drug_ & bug ==bug_)
  #fitted$age<-as.character(fitted$age)
  #pe<-get(paste0(variable))$partial_effects %>% filter(drug==drug_)
  #pe<-filter(pe,drug==drug_ & bug==bug_)
  #pe$age<-as.character(pe$age)
  #fitted<-left_join(fitted,pe,by=c("bug","age","date","drug"))
  #fitted$sig<-ifelse(fitted$low.y <0 & fitted$high.y <0,"1",ifelse(fitted$low.y >0 & fitted$high.y >0,"1","0"))
  #fitted$age<-as.numeric(fitted$age)
  fitted$drug<- case_when(
    fitted$drug == "e_AMOXCLAV" ~ "co-amoxiclav",
    fitted$drug =="e_CIP" ~ "ciprofloxacin",
    fitted$drug == "e_CXM" ~ "co-trimoxazole",
    fitted$drug == "e_GENT" ~ "gentamicin",
    fitted$drug == "e_MER" ~ "meropenem",
    fitted$drug == "e_PIPTAZ" ~ "piperacillin-tazobactam",
    fitted$drug == "threegc" ~ "third generation cephalosporin"
  )
  pe$drug<- case_when(
    pe$drug == "e_AMOXCLAV" ~ "co-amoxiclav",
    pe$drug =="e_CIP" ~ "ciprofloxacin",
    pe$drug == "e_CXM" ~ "co-trimoxazole",
    pe$drug == "e_GENT" ~ "gentamicin",
    pe$drug == "e_MER" ~ "meropenem",
    pe$drug == "e_PIPTAZ" ~ "piperacillin-tazobactam",
    pe$drug == "threegc" ~ "third generation cephalosporin"
  )
  

  p<-ggplot(fitted) +
    aes_string(x="date",y="est",group=variable,color=variable) +
    geom_line() +
    #geom_segment(data=filter(fitted,sig==1),aes(x=date,xend=date+days(14),y=est.x),size=1.1) +
    theme_minimal() + theme(axis.text.x = element_text(angle=45),legend.position = "None") +
    ylab(ifelse(labs==TRUE,"Proportion resistant","")) + xlab("") +
    guides(fill="none") +
    labs(color=variable_clean) +
      scale_color_viridis_c() +
    ggtitle(ifelse(labs==TRUE,drug_clean,""))
  
  
  pp<-ggplot(pe)+
  aes(x=date,y=est,group=bug)  +
  geom_line() +
  geom_ribbon(data=pe,aes(x=date,ymin=low,ymax=high),alpha=0.1)+
    geom_hline(yintercept = 0,linetype="dashed")+
  theme_minimal() + theme(axis.text.x = element_text(angle=45)) +
    ylab(paste0("s(",variable_clean,")")) + xlab("") +
 facet_wrap(as.formula(paste("~", variable)))
 
  return(list(p=p,pp=pp))
}

age_data_frames<-ls(pattern = "*_age_i")
age<-map_dfr(age_data_frames,get)


imd_data_frames<-ls(pattern="*_imd_i")
imd_score<-map_dfr(imd_data_frames,get)

admi_data_frames<-ls(pattern="*_n_admissions_i")
n_admissions<-map_dfr(admi_data_frames,get)

elix_data_frames<-ls(pattern="*_elixhauser_i")
elixhauser<-map_dfr(elix_data_frames,get)


age_names <- ls(pattern = c("*_age_i"))
dfs<-mget(age_names)
fitted_age<-lapply(dfs,function(x) x$fitted)
fitted_age<-bind_rows(fitted_age)
fitted_age$group<-"age"
fitted_age<-rename(fitted_age,variable=age)

imd_names <- ls(pattern = c("*_imd_i"))
dfs<-mget(imd_names)
fitted_imd<-lapply(dfs,function(x) x$fitted)
fitted_imd<-bind_rows(fitted_imd)
fitted_imd$group<-"imd"
fitted_imd<-rename(fitted_imd,variable=imd_score)

admi_names <- ls(pattern = c("*_n_admissions_i"))
dfs<-mget(admi_names)
fitted_admi<-lapply(dfs,function(x) x$fitted)
fitted_admi<-bind_rows(fitted_admi)
fitted_admi$group<-"admissions"
fitted_admi<-rename(fitted_admi,variable=n_admissions)

elix_names <- ls(pattern = c("*_elixhauser_i"))
dfs<-mget(elix_names)
fitted_elix<-lapply(dfs,function(x) x$fitted)
fitted_elix<-bind_rows(fitted_elix)
fitted_elix$group<-"elixhauser"
fitted_elix<-rename(fitted_elix,variable=elixhauser)

fitted<-rbind(fitted_age,fitted_imd,fitted_admi,fitted_elix)

fitted$drug <- case_when(
  fitted$drug == "e_AMOXCLAV" ~ "co-amoxiclav",
  fitted$drug == "e_CIP" ~ "ciprofloxacin",
  fitted$drug == "e_CXM" ~ "co-trimoxazole",
  fitted$drug == "e_GENT" ~ "gentamicin",
  fitted$drug == "e_MER" ~ "meropenem",
  fitted$drug == "e_PIPTAZ" ~ "piperacillin-tazobactam",
  fitted$drug == "threegc" ~ "third generation cephalosporin",
  TRUE ~ fitted$drug 
)

plots <- fitted %>% #filter(bug=="ecoli") %>%
split(.$group) %>%
lapply(function(df) {
ggplot(df, aes(x = date, y = est, group = variable, color = variable)) +
geom_line() +
facet_grid(drug ~ bug, scales = "free_y") +
scale_color_viridis_c() +
theme_minimal() +
theme(legend.position = "bottom",strip.text = element_text(size=6),strip.background = element_rect(fill="white")) + 
ggtitle(unique(df$group)) +
    labs(color=df$group)
})

combined_plot <- wrap_plots(plots, ncol = length(plots))


dfs<-mget(ecol_names)
derivs_ecol<-lapply(dfs,function(x) x$deriv)
derivs_ecol<-bind_rows(derivs_ecol)
fitted_ecol<-lapply(dfs,function(x) x$fitted)
fitted_ecol<-bind_rows(fitted_ecol)
fitted_kleb<-bind_rows(fitted_kleb)













int<-ggplot(samples) +
  aes(x=date,y=est) +
  geom_line() +
  facet_wrap(~age) +
  scale_x_date(breaks = seq(as.Date("2012-01-01"), as.Date("2024-01-01"), by="2 years"),
               date_labels = "%Y")  +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90))+
  ylab("ti(time,age)")

data$age_group <- cut(data$age,
                    breaks = c(0, 1, 18, 40, 60, 80, Inf),
                    include.lowest = TRUE,
                    right = FALSE,
                    labels = c("<1", "1-17", "18-39", "40-59", "60-79", "80+"))

data$monthno
prop <- data %>%
  group_by(year, monthno,age_group) %>%
  summarise(S=sum(e_AMOXCLAV == 0, na.rm = TRUE),
            R=sum(e_AMOXCLAV == 1, na.rm = TRUE))

prop$date<-ym(paste0(prop$year,prop$monthno))
prop<-prop %>% pivot_longer(cols = 4:5,names_to = "which",values_to = "n")

cols<-viridis(6)
cols2<-magma(6)
cols<-c(cols2,cols)

prop<-ggplot(prop) +
  aes(x=date,y=n,fill=interaction(age_group,which)) +
  geom_bar(position="stack", stat="identity") +
  scale_x_date(breaks = seq(as.Date("2012-01-01"), as.Date("2024-01-01"), by="1 year"),
               date_labels = "%Y") +
  theme_minimal() + labs(fill="") +
  scale_fill_manual(values=cols)+
  ylab("Number of isolates")



pe<-ggplot(samples) +
  aes(x=age,y=est) +
  geom_line()  +
  theme_minimal() +
  ylab("s(age)")

times<-seq(0,max(data$timeline,na.rm=T),14)
ages<-seq(0,100,10)
predict_data<-expand_grid(timeline=times,age=ages)
predict_data$date<-min(data$spec_dt,na.rm=T) + days(predict_data$timeline)
predict_data$day<-yday(predict_data$date)

samples<-fitted_samples(coamox_model_age,data=predict_data,n=1000)
samples<-samples %>%
  group_by(.row) %>%
  summarise(low=quantile(.fitted,probs=c(0.025)),
            est=quantile(.fitted,probs=c(0.5)),
            high=quantile(.fitted,probs=c(0.975)))
samples$age<-predict_data$age
samples$date<-predict_data$date

fe<-ggplot(samples) +
  aes(x=date,y=est,color=age,group=age) +
  geom_line()  +
  theme_minimal() +
  ylab("Estimated proportion of co-amoxiclav resistance") +
  scale_color_viridis_c()

coamox_age<-(int/prop + plot_layout(heights=c(3,1), guides="collect")) | (pe/fe)
coamox_age + plot_annotation("Coamoxiclav:age")









#######################categorical########################

#these functions do the same as above but for categorical variables

main_discrete<-function(drug,bug,seasonal,variable){
  if(seasonal==TRUE){
    formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25) + s(day, bs='cc', k=12) +
                                " ,variable))
  }else{
    formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25) +
                                " ,variable))
    
  }
  data_<-filter(data,species==bug)
  cat("\nfitting model")
  model<-bam(formula,
             data=data_,
             family = binomial(link = "logit"),
             method = "REML")
  
  cat("\nMaking dummy data")
  times<-seq(0,max(data$timeline,na.rm=T),14)
  variable_<-unique(data[[variable]])
  variable_<-variable_[!is.na(variable_)]
  predict_data<-expand_grid(timeline=times,variable=variable_)
  predict_data$date<-min(data$spec_dt,na.rm=T) + days(predict_data$timeline)
  predict_data$day<-yday(predict_data$date)
  predict_data<-filter(predict_data,timeline==4018) #no int with
  y=sym(variable)
  predict_data<-rename(predict_data, !!y := variable)
  
  
  x=paste0("s(",variable,")")
  cat("\nTaking smooth samples")
  #samples<-smooth_samples(model,select = x,n = 1000,data = predict_data)
  samples<-fitted_samples(model,data = predict_data,n=1000)
  cat("\nGenerating posterior credibility intervals")
  samples<-samples %>%
    group_by(.row) %>%
    summarise(low=quantile(.fitted,probs=c(0.025)),
              est=quantile(.fitted,probs=c(0.5)),
              high=quantile(.fitted,probs=c(0.975)))
  
  samples[[variable]]<-predict_data[[variable]]
  samples$drug<-drug
  samples$bug<-bug
  
  return(samples)
}

ecoli_coamox_he<-main_discrete("e_AMOXCLAV","ecoli",TRUE,"healthcare_exposure")
ecoli_gent_he<-main_discrete("e_GENT","ecoli",TRUE,"healthcare_exposure")
ecoli_cip_he<-main_discrete("e_CIP","ecoli",TRUE,"healthcare_exposure")
ecoli_piptaz_he<-main_discrete("e_PIPTAZ","ecoli",TRUE,"healthcare_exposure")
ecoli_cotrim_he<-main_discrete("e_CXM","ecoli",TRUE,"healthcare_exposure")
ecoli_tgc_he<-main_discrete("threegc","ecoli",TRUE,"healthcare_exposure")

kleb_coamox_he<-main_discrete("e_AMOXCLAV","kleb",TRUE,"healthcare_exposure")
kleb_gent_he<-main_discrete("e_GENT","kleb",FALSE,"healthcare_exposure")
kleb_cip_he<-main_discrete("e_CIP","kleb",FALSE,"healthcare_exposure")
kleb_piptaz_he<-main_discrete("e_PIPTAZ","kleb",FALSE,"healthcare_exposure")
kleb_cotrim_he<-main_discrete("e_CXM","kleb",FALSE,"healthcare_exposure")
kleb_tgc_he<-main_discrete("threegc","kleb",FALSE,"healthcare_exposure")

ecoli_coamox_ethnos<-main_discrete("e_AMOXCLAV","ecoli",TRUE,"ETHNOS_short")
ecoli_gent_ethnos<-main_discrete("e_GENT","ecoli",TRUE,"ETHNOS_short")
ecoli_cip_ethnos<-main_discrete("e_CIP","ecoli",TRUE,"ETHNOS_short")
ecoli_piptaz_ethnos<-main_discrete("e_PIPTAZ","ecoli",TRUE,"ETHNOS_short")
ecoli_cotrim_ethnos<-main_discrete("e_CXM","ecoli",TRUE,"ETHNOS_short")
ecoli_tgc_ethnos<-main_discrete("threegc","ecoli",TRUE,"ETHNOS_short")

kleb_coamox_ethnos<-main_discrete("e_AMOXCLAV","kleb",TRUE,"ETHNOS_short")
kleb_gent_ethnos<-main_discrete("e_GENT","kleb",FALSE,"ETHNOS_short")
kleb_cip_ethnos<-main_discrete("e_CIP","kleb",FALSE,"ETHNOS_short")
kleb_piptaz_ethnos<-main_discrete("e_PIPTAZ","kleb",FALSE,"ETHNOS_short")
kleb_cotrim_ethnos<-main_discrete("e_CXM","kleb",FALSE,"ETHNOS_short")
kleb_tgc_ethnos<-main_discrete("threegc","kleb",FALSE,"ETHNOS_short")

ecoli_coamox_rururb<-main_discrete("e_AMOXCLAV","ecoli",TRUE,"rururb")
ecoli_gent_rururb<-main_discrete("e_GENT","ecoli",TRUE,"rururb")
ecoli_cip_rururb<-main_discrete("e_CIP","ecoli",TRUE,"rururb")
ecoli_piptaz_rururb<-main_discrete("e_PIPTAZ","ecoli",TRUE,"rururb")
ecoli_cotrim_rururb<-main_discrete("e_CXM","ecoli",TRUE,"rururb")
ecoli_tgc_rururb<-main_discrete("threegc","ecoli",TRUE,"rururb")

kleb_coamox_rururb<-main_discrete("e_AMOXCLAV","kleb",TRUE,"rururb")
kleb_gent_rururb<-main_discrete("e_GENT","kleb",FALSE,"rururb")
kleb_cip_rururb<-main_discrete("e_CIP","kleb",FALSE,"rururb")
kleb_piptaz_rururb<-main_discrete("e_PIPTAZ","kleb",FALSE,"rururb")
kleb_cotrim_rururb<-main_discrete("e_CXM","kleb",FALSE,"rururb")
kleb_tgc_rururb<-main_discrete("threegc","kleb",FALSE,"rururb")

ecoli_coamox_sex<-main_discrete("e_AMOXCLAV","ecoli",TRUE,"sex")
ecoli_gent_sex<-main_discrete("e_GENT","ecoli",TRUE,"sex")
ecoli_cip_sex<-main_discrete("e_CIP","ecoli",TRUE,"sex")
ecoli_piptaz_sex<-main_discrete("e_PIPTAZ","ecoli",TRUE,"sex")
ecoli_cotrim_sex<-main_discrete("e_CXM","ecoli",TRUE,"sex")
ecoli_tgc_sex<-main_discrete("threegc","ecoli",TRUE,"sex")

kleb_coamox_sex<-main_discrete("e_AMOXCLAV","kleb",TRUE,"sex")
kleb_gent_sex<-main_discrete("e_GENT","kleb",FALSE,"sex")
kleb_cip_sex<-main_discrete("e_CIP","kleb",FALSE,"sex")
kleb_piptaz_sex<-main_discrete("e_PIPTAZ","kleb",FALSE,"sex")
kleb_cotrim_sex<-main_discrete("e_CXM","kleb",FALSE,"sex")
kleb_tgc_sex<-main_discrete("threegc","kleb",FALSE,"sex")


he_data_frames<-ls(pattern = "*_he")
healthcare_exposure<-map_dfr(he_data_frames,get)
healthcare_exposure$healthcare_exposure<-factor(healthcare_exposure$healthcare_exposure,levels=c("nosocomial","quasi-nosocomial","quasi-community","community"))

ethnos_data_frames<-ls(pattern="*_ethnos")
ETHNOS_short<-map_dfr(ethnos_data_frames,get)

rururb_data_frames<-ls(pattern="*_rururb")
rururb<-map_dfr(rururb_data_frames,get)
rururb$rururb<-factor(rururb$rururb,levels=c("Urban conurbation","Urban","Town","Village"))

sex_data_frames<-ls(pattern = "*_sex")
sex<-map_dfr(sex_data_frames,get)

library(ggplot2)
library(dplyr)

plot_main_discrete <- function(variable,labs,variable_clean) {
  
  fitted <- get(paste0(variable))
  
  fitted$drug <- case_when(
    fitted$drug == "e_AMOXCLAV" ~ "co-amoxiclav",
    fitted$drug == "e_CIP" ~ "ciprofloxacin",
    fitted$drug == "e_CXM" ~ "co-trimoxazole",
    fitted$drug == "e_GENT" ~ "gentamicin",
    fitted$drug == "e_MER" ~ "meropenem",
    fitted$drug == "e_PIPTAZ" ~ "piperacillin-tazobactam",
    fitted$drug == "threegc" ~ "third generation cephalosporin",
    TRUE ~ fitted$drug # default case to handle any unmatched values
  )
  
 
  p <- ggplot(fitted, aes_string(x = variable, y = "est", color = "bug")) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes_string(ymax = "high", ymin = "low", fill = "bug"), 
                  position = position_dodge(width = 0.5), width = 0.3, alpha = 0.4, size = 1) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle=-45),
      legend.position = "none"
    ) +
    ylab(ifelse(labs==TRUE,"Estimated Proportion Resistant","")) +
    xlab(variable_clean) +
    guides(fill = guide_legend(title = "Bug")) +
    facet_wrap(~drug, ncol = 1, scales = "free_y") +
    labs(color = "Antibiotic", fill = "Antibiotic")
  
  return(p)
}



plot_main_discrete("healthcare_exposure",TRUE)
plot_main_discrete("ETHNOS_short")
plot_main_discrete("rururb")
plot_main_discrete("sex")

plot_main_continuous("age",TRUE,"Age")+ plot_main_continuous("imd_score",FALSE,"IMD Score") +plot_main_continuous("n_admissions",FALSE,"N. admissions") + plot_main_continuous("elixhauser",FALSE,"Elixhauser") + 
  plot_main_discrete("healthcare_exposure",FALSE,"Healthcare Exposure") + plot_main_discrete("ETHNOS_short",FALSE,"Ethnicity") + plot_main_discrete("rururb",FALSE,"Rural/Urban") + plot_main_discrete("sex",FALSE,"Sex") +
  plot_layout(nrow = 1)

data$rururb<- case_when(
  data$`Rural Urban Classification 2011 (10 fold)` == 'Rural town and fringe' ~ 'Town',
  data$`Rural Urban Classification 2011 (10 fold)` == 'Rural town and fringein a sparse setting' ~ 'Town',
  data$`Rural Urban Classification 2011 (10 fold)` == 'Rural village and dispersed' ~ 'Village',
  data$`Rural Urban Classification 2011 (10 fold)` == 'Rural village and dispersed in a sparse setting' ~ 'Village',
  data$`Rural Urban Classification 2011 (10 fold)` == 'Urban city and town' ~ 'Urban',
  data$`Rural Urban Classification 2011 (10 fold)` == 'Urban city and town in a sparse setting' ~ 'Urban',
  data$`Rural Urban Classification 2011 (10 fold)` == 'Urban major conurbation' ~ 'Urban conurbation',
  data$`Rural Urban Classification 2011 (10 fold)` == 'Urban minor conurbation ' ~ 'Urban conurbation'
    
)

plot_main_discrete("rurub")


i_discrete<-function(drug,bug,seasonal,variable){
  if(seasonal==TRUE){
    formula <- as.formula(paste(drug, "~ s(timeline, by=",variable,", bs='ad', m=3) + s(day, bs='cc', k=12) +
                                ",variable))
  }else{
    formula <- as.formula(paste(drug, "~ s(timeline, by=",variable,", bs='ad', m=3) +
                                ",variable))
    
  }
  data_<-filter(data,species==bug)
  cat("\nfitting model")
  model<-bam(formula,
             data=data_,
             family = binomial(link = "logit"),
             method = "REML")
  
  cat("\nmaking dummy data")
  times<-seq(0,max(data$timeline,na.rm=T),14)
  
  variable_<-unique(data_[[variable]])
  variable_<-variable_[!is.na(variable_)]
  predict_data<-expand_grid(timeline=times,variable=variable_)
  predict_data$date<-min(data$spec_dt,na.rm=T) + days(predict_data$timeline)
  predict_data$day<-yday(predict_data$date)
  y=sym(variable)
  predict_data<-rename(predict_data, !!y := variable)
  

  cat("\nFitted samples")
  fitted<-fitted_samples(model = model,n = 1000,data = predict_data)
  fitted<-fitted %>% 
    group_by(.row) %>%
    summarise(low=quantile(.fitted,probs=c(0.025)),
              est=quantile(.fitted,probs=c(0.5)),
              high=quantile(.fitted,probs=c(0.975)))
  fitted$date<-predict_data$date
  fitted[[variable]]<-predict_data[[variable]]
  fitted$bug<-bug
  fitted$drug<-drug
  cat("\nDone!")
  return(list(fitted=fitted))
}
ecoli_coamox_ethnos_i<-i_discrete("e_AMOXCLAV","ecoli",seasonal = TRUE,variable = "ETHNOS_short")
ecoli_gent_ethnos_i<-i_discrete("e_GENT","ecoli",seasonal = TRUE,variable = "ETHNOS_short")
ecoli_piptaz_ethnos_i<-i_discrete("e_PIPTAZ","ecoli",seasonal = TRUE,variable = "ETHNOS_short")
ecoli_cipro_ethnos_i<-i_discrete("e_CIP","ecoli",seasonal = TRUE,variable = "ETHNOS_short")
ecoli_cotrim_ethnos_i<-i_discrete("e_CXM","ecoli",seasonal = TRUE,variable = "ETHNOS_short")
ecoli_tgc_ethnos_i<-i_discrete("threegc","ecoli",seasonal = TRUE,variable = "ETHNOS_short")

kleb_coamox_ethnos_i<-i_discrete("e_AMOXCLAV","kleb",seasonal = TRUE,variable = "ETHNOS_short")
kleb_gent_ethnos_i<-i_discrete("e_GENT","kleb",seasonal = TRUE,variable = "ETHNOS_short")
kleb_piptaz_ethnos_i<-i_discrete("e_PIPTAZ","kleb",seasonal = TRUE,variable = "ETHNOS_short")
kleb_cipro_ethnos_i<-i_discrete("e_CIP","kleb",seasonal = TRUE,variable = "ETHNOS_short")
kleb_cotrim_ethnos_i<-i_discrete("e_CXM","kleb",seasonal = TRUE,variable = "ETHNOS_short")
kleb_tgc_ethnos_i<-i_discrete("threegc","kleb",seasonal = TRUE,variable = "ETHNOS_short")

data$sex<-as.factor(data$sex)
ecoli_coamox_sex_i<-i_discrete("e_AMOXCLAV","ecoli",seasonal = TRUE,variable = "sex")
ecoli_gent_sex_i<-i_discrete("e_GENT","ecoli",seasonal = TRUE,variable = "sex")
ecoli_piptaz_sex_i<-i_discrete("e_PIPTAZ","ecoli",seasonal = TRUE,variable = "sex")
ecoli_cipro_sex_i<-i_discrete("e_CIP","ecoli",seasonal = TRUE,variable = "sex")
ecoli_cotrim_sex_i<-i_discrete("e_CXM","ecoli",seasonal = TRUE,variable = "sex")
ecoli_tgc_sex_i<-i_discrete("threegc","ecoli",seasonal = TRUE,variable = "sex")

kleb_coamox_sex_i<-i_discrete("e_AMOXCLAV","kleb",seasonal = TRUE,variable = "sex")
kleb_gent_sex_i<-i_discrete("e_GENT","kleb",seasonal = FALSE,variable = "sex")
kleb_piptaz_sex_i<-i_discrete("e_PIPTAZ","kleb",seasonal = FALSE,variable = "sex")
kleb_cipro_sex_i<-i_discrete("e_CIP","kleb",seasonal = FALSE,variable = "sex")
kleb_cotrim_sex_i<-i_discrete("e_CXM","kleb",seasonal = FALSE,variable = "sex")
kleb_tgc_sex_i<-i_discrete("threegc","kleb",seasonal = FALSE,variable = "sex")

data$rururb<-factor(data$rururb,levels=c("Urban conurbation","Urban","Town","Village"))
ecoli_coamox_rururb_i<-i_discrete("e_AMOXCLAV","ecoli",seasonal = TRUE,variable = "rururb")
ecoli_gent_rururb_i<-i_discrete("e_GENT","ecoli",seasonal = TRUE,variable = "rururb")
ecoli_piptaz_rururb_i<-i_discrete("e_PIPTAZ","ecoli",seasonal = TRUE,variable = "rururb")
ecoli_cipro_rururb_i<-i_discrete("e_CIP","ecoli",seasonal = TRUE,variable = "rururb")
ecoli_cotrim_rururb_i<-i_discrete("e_CXM","ecoli",seasonal = TRUE,variable = "rururb")
ecoli_tgc_rururb_i<-i_discrete("threegc","ecoli",seasonal = TRUE,variable = "rururb")

kleb_coamox_rururb_i<-i_discrete("e_AMOXCLAV","kleb",seasonal = TRUE,variable = "rururb")
kleb_gent_rururb_i<-i_discrete("e_GENT","kleb",seasonal = FALSE,variable = "rururb")
kleb_piptaz_rururb_i<-i_discrete("e_PIPTAZ","kleb",seasonal = FALSE,variable = "rururb")
kleb_cipro_rururb_i<-i_discrete("e_CIP","kleb",seasonal = FALSE,variable = "rururb")
kleb_cotrim_rururb_i<-i_discrete("e_CXM","kleb",seasonal = FALSE,variable = "rururb")
kleb_tgc_rururb_i<-i_discrete("threegc","kleb",seasonal = FALSE,variable = "rururb")

sex_names <- ls(pattern = c("*_sex_i"))
dfs<-mget(sex_names)
fitted_sex<-lapply(dfs,function(x) x$fitted)
fitted_sex<-bind_rows(fitted_sex)
fitted_sex$group<-"sex"
fitted_sex<-rename(fitted_sex,variable=sex)

rururb_names <- ls(pattern = c("*_rururb_i"))
dfs<-mget(rururb_names)
fitted_rururb<-lapply(dfs,function(x) x$fitted)
fitted_rururb<-bind_rows(fitted_rururb)
fitted_rururb$group<-"rururb"
fitted_rururb<-rename(fitted_rururb,variable=rururb)

ethnos_names<-ls(pattern=c("*_ethnos_i"))
dfs<-mget(ethnos_names)
fitted_ethnos<-lapply(dfs,function(x) x$fitted)
fitted_ethnos<-bind_rows(fitted_ethnos)
fitted_ethnos$group<-"ethnos"
fitted_ethnos<-rename(fitted_ethnos,variable=ETHNOS_short)

fitted<-rbind(fitted_sex,fitted_rururb,fitted_ethnos)

fitted$drug <- case_when(
  fitted$drug == "e_AMOXCLAV" ~ "co-amoxiclav",
  fitted$drug == "e_CIP" ~ "ciprofloxacin",
  fitted$drug == "e_CXM" ~ "co-trimoxazole",
  fitted$drug == "e_GENT" ~ "gentamicin",
  fitted$drug == "e_MER" ~ "meropenem",
  fitted$drug == "e_PIPTAZ" ~ "piperacillin-tazobactam",
  fitted$drug == "threegc" ~ "third generation cephalosporin",
  TRUE ~ fitted$drug 
)

plots <- fitted %>% #filter(bug=="ecoli") %>%
  split(.$group) %>%
  lapply(function(df) {
    ggplot(df, aes(x = date, y = est, group = variable, color = variable)) +
      geom_line() +
      facet_grid(drug ~ bug, scales = "free_y") +
      scale_color_viridis_d() +
      theme_minimal() +
      theme(legend.position = "bottom",strip.text = element_text(size=6),strip.background = element_rect(fill="white")) + 
      ggtitle(unique(df$group)) +
      labs(color=df$group)
  })

combined_plot <- wrap_plots(plots, ncol = length(plots))


##########multivariable###################

main_continuous<-function(drug,bug,seasonal,variable,slices){
  if(seasonal==TRUE){
    formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25) + s(day, bs='cc', k=12) +
                                s(age,bs='tp', k=25) + s(imd_score, bs='tp', k=25) + s(n_admissions, bs='tp',k=6) +
                                s(elixhauser, bs='tp', k=12) + healthcare_exposure + ETHNOS_short + rururb + sex"))
                                
  }else{
    formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25)  +
                                s(age,bs='tp', k=25) + s(imd_score, bs='tp', k=25) + s(n_admissions, bs='tp',k=6) +
                                s(elixhauser, bs='tp', k=12) + healthcare_exposure + ETHNOS_short + rururb + sex"))
    
  }
  data_<-filter(data,species==bug)
  cat("\nfitting model")
  model<-bam(formula,
             data=data_,
             family = binomial(link = "logit"),
             method = "REML")
  
  cat("\nMaking dummy data")
  predict_data<-select(data_,timeline,age,imd_score,n_admissions,elixhauser,healthcare_exposure,ETHNOS_short,rururb,sex)
  predict_data$ETHNOS_short<-as.factor(predict_data$ETHNOS_short)
  predict_data$rururb<-as.factor(predict_data$rururb)
  predict_data$healthcare_exposure<-as.factor(predict_data$healthcare_exposure)
  predict_data$sex<-as.factor(predict_data$sex)
  predict_data<-data_slice(predict_data,!!sym(variable) := evenly(!!sym(variable), n =slices),
                           healthcare_exposure = level(healthcare_exposure,"community"))
  predict_data$date<-min(data$spec_dt,na.rm=T) + days(predict_data$timeline)
  predict_data$day<-yday(predict_data$date)
  

  
  
  x=paste0("s(",variable,")")
  cat("\nTaking smooth samples")
  #samples<-smooth_samples(model,select = x,n = 1000,data = predict_data)
  samples<-fitted_samples(model,data = predict_data,n=1000)
  cat("\nGenerating posterior credibility intervals")
  samples<-samples %>%
    group_by(.row) %>%
    summarise(low=quantile(.fitted,probs=c(0.025)),
              est=quantile(.fitted,probs=c(0.5)),
              high=quantile(.fitted,probs=c(0.975)))
  
  samples[[variable]]<-predict_data[[variable]]
  samples$drug<-drug
  samples$bug<-bug
  
  return(samples)
}


i_continuous<-function(drug,bug,seasonal,variable,slices){
  if(seasonal==TRUE){
    formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25) + s(day, bs='cc', k=12) +
                                s(age,bs='tp', k=25) + s(imd_score, bs='tp', k=25) + s(n_admissions, bs='tp',k=6) +
                                s(elixhauser, bs='tp', k=12) + healthcare_exposure + ETHNOS_short + rururb + sex +
                                ti(",variable,",timeline)"))
  }else{
    formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25) +
                                  s(age,bs='tp', k=25) + s(imd_score, bs='tp', k=25) + s(n_admissions, bs='tp',k=6) +
                                  s(elixhauser, bs='tp', k=12) + healthcare_exposure + ETHNOS_short + rururb + sex +
                                ti(",variable,",timeline)"))
    
  }
  data_<-filter(data,species==bug)
  cat("\nfitting model")
  model<-bam(formula,
             data=data_,
             family = binomial(link = "logit"),
             method = "REML")
  
  cat("\nmaking dummy data")
  predict_data<-select(data_,timeline,age,imd_score,n_admissions,elixhauser,healthcare_exposure,ETHNOS_short,rururb,sex)
  predict_data$ETHNOS_short<-as.factor(predict_data$ETHNOS_short)
  predict_data$rururb<-as.factor(predict_data$rururb)
  predict_data$healthcare_exposure<-as.factor(predict_data$healthcare_exposure)
  predict_data$sex<-as.factor(predict_data$sex)
  predict_data<-data_slice(predict_data,!!sym(variable) := evenly(!!sym(variable), n =slices),
                           timeline=evenly(timeline,n=1000),
                           healthcare_exposure = level(healthcare_exposure,"community"))
  predict_data$date<-min(data$spec_dt,na.rm=T) + days(round(predict_data$timeline,0))
  predict_data$day<-yday(predict_data$date)
  
  x=paste0("ti(",variable,",timeline)")
  cat("\nTaking smooth samples")
  samples<-smooth_samples(model,select = x,n = 1000,data = predict_data)
  
  cat("\nCredibility intervals")
  samples<-samples %>%
    group_by(.row) %>%
    summarise(low=quantile(.value,probs=c(0.025)),
              est=quantile(.value,probs=c(0.5)),
              high=quantile(.value,probs=c(0.975)))
  samples$date<-predict_data$date
  samples[[variable]]<-predict_data[[variable]]
  samples$drug<-drug
  samples$bug<-bug
  
  cat("\nFitted samples")
  fitted<-fitted_samples(model = model,n = 1000,data = predict_data)
  fitted<-fitted %>% 
    group_by(.row) %>%
    summarise(low=quantile(.fitted,probs=c(0.025)),
              est=quantile(.fitted,probs=c(0.5)),
              high=quantile(.fitted,probs=c(0.975)))
  fitted$date<-predict_data$date
  fitted[[variable]]<-predict_data[[variable]]
  fitted$bug<-bug
  fitted$drug<-drug
  cat("\nDone!")
  return(list(partial_effects=samples,fitted=fitted))
}


main_discrete<-function(drug,bug,seasonal,variable){
  if(seasonal==TRUE){
    formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25) + s(day, bs='cc', k=12) +
                                s(age,bs='tp', k=25) + s(imd_score, bs='tp', k=25) + s(n_admissions, bs='tp',k=6) +
                                s(elixhauser, bs='tp', k=12) + healthcare_exposure + ETHNOS_short + rururb + sex"))
  }else{
    formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25) +
                                s(age,bs='tp', k=25) + s(imd_score, bs='tp', k=25) + s(n_admissions, bs='tp',k=6) +
                                s(elixhauser, bs='tp', k=12) + healthcare_exposure + ETHNOS_short + rururb + sex"))
    
  }
  data_<-filter(data,species==bug)
  cat("\nfitting model")
  model<-bam(formula,
             data=data_,
             family = binomial(link = "logit"),
             method = "REML")
  
  cat("\nMaking dummy data")
  predict_data<-select(data_,timeline,age,imd_score,n_admissions,elixhauser,healthcare_exposure,ETHNOS_short,rururb,sex)
  predict_data$ETHNOS_short<-as.factor(predict_data$ETHNOS_short)
  predict_data$rururb<-as.factor(predict_data$rururb)
  predict_data$healthcare_exposure<-as.factor(predict_data$healthcare_exposure)
  predict_data$sex<-as.factor(predict_data$sex)
  predict_data<-data_slice(predict_data,!!sym(variable) := factor_combos(model, vars = !!sym(variable)))
  predict_data[[variable]]<-predict_data[[variable]][[variable]]
  predict_data$date<-min(data$spec_dt,na.rm=T) + days(round(predict_data$timeline,0))
  predict_data$day<-yday(predict_data$date)
  
  x=paste0("s(",variable,")")
  cat("\nTaking smooth samples")
  #samples<-smooth_samples(model,select = x,n = 1000,data = predict_data)
  samples<-fitted_samples(model,data = predict_data,n=1000)
  cat("\nGenerating posterior credibility intervals")
  samples<-samples %>%
    group_by(.row) %>%
    summarise(low=quantile(.fitted,probs=c(0.025)),
              est=quantile(.fitted,probs=c(0.5)),
              high=quantile(.fitted,probs=c(0.975)))
  
  samples[[variable]]<-predict_data[[variable]]
  samples$drug<-drug
  samples$bug<-bug
  
  return(samples)
}


i_discrete<-function(drug,bug,seasonal,variable){
  if(seasonal==TRUE){
    formula <- as.formula(paste(drug, "~ s(timeline, by=",varible," bs='ad', m=3, k=25) + s(day, bs='cc', k=12) +
                                s(age,bs='tp', k=25) + s(imd_score, bs='tp', k=25) + s(n_admissions, bs='tp',k=6) +
                                s(elixhauser, bs='tp', k=12) + healthcare_exposure + ETHNOS_short + rururb + sex"))
  }else{
    formula <- as.formula(paste(drug, "~ s(timeline, by=",varible," bs='ad', m=3, k=25) +
                                s(age,bs='tp', k=25) + s(imd_score, bs='tp', k=25) + s(n_admissions, bs='tp',k=6) +
                                s(elixhauser, bs='tp', k=12) + healthcare_exposure + ETHNOS_short + rururb + sex"))
    
  }
  data_<-filter(data,species==bug)
  cat("\nfitting model")
  model<-bam(formula,
             data=data_,
             family = binomial(link = "logit"),
             method = "REML")
  
  cat("\nmaking dummy data")
  predict_data<-select(data_,timeline,age,imd_score,n_admissions,elixhauser,healthcare_exposure,ETHNOS_short,rururb,sex)
  predict_data$ETHNOS_short<-as.factor(predict_data$ETHNOS_short)
  predict_data$rururb<-as.factor(predict_data$rururb)
  predict_data$healthcare_exposure<-as.factor(predict_data$healthcare_exposure)
  predict_data$sex<-as.factor(predict_data$sex)
  predict_data<-data_slice(predict_data,!!sym(variable) := factor_combos(model, vars = !!sym(variable)),
                           timeline=evenly(timeline,n=1000))
  predict_data[[variable]]<-predict_data[[variable]][[variable]]
  predict_data$date<-min(data$spec_dt,na.rm=T) + days(round(predict_data$timeline,0))
  predict_data$day<-yday(predict_data$date)
  
  
  cat("\nFitted samples")
  fitted<-fitted_samples(model = model,n = 1000,data = predict_data)
  fitted<-fitted %>% 
    group_by(.row) %>%
    summarise(low=quantile(.fitted,probs=c(0.025)),
              est=quantile(.fitted,probs=c(0.5)),
              high=quantile(.fitted,probs=c(0.975)))
  fitted$date<-predict_data$date
  fitted[[variable]]<-predict_data[[variable]]
  fitted$bug<-bug
  fitted$drug<-drug
  cat("\nDone!")
  return(list(fitted=fitted))
}


####other random bits of exploratory code##############

#############plot overall seasonality##################

seasonality<-function(drug,bug){

  formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25) + s(day, bs='cc', k=12)"))

  data_<-filter(data,species==bug)
  cat("\nfitting model")
  model<-bam(formula,
             data=data_,
             family = binomial(link = "logit"),
             method = "REML")

  times<-seq(0,max(data$timeline,na.rm=T),1)
  predict_data<-expand_grid(timeline=times)
  predict_data$date<-min(data$spec_dt,na.rm=T) + days(predict_data$timeline)
  predict_data$day<-yday(predict_data$date)

  estimates<-smooth_samples(model,data = predict_data,n = 1000,select = "s(day)")
  estimates<-estimates %>%
    group_by(.row) %>%
    summarise(low=quantile(.value,probs=c(0.025)),
              est=quantile(.value,probs=c(0.5)),
              high=quantile(.value,probs=c(0.975)))
  estimates$date<-predict_data$date
  estimates$day<-yday(predict_data$date)
  estimates<-distinct(estimates,day,.keep_all = T)
  estimates$drug<-drug
  estimates$bug<-bug
  return(estimates)
}

ecoamox<-seasonality("e_AMOXCLAV","ecoli")
egent<-seasonality("e_GENT","ecoli")
ecip<-seasonality("e_CIP","ecoli")
etaz<-seasonality("e_PIPTAZ","ecoli")
etgc<-seasonality("threegc","ecoli")
emer<-seasonality("e_MER","ecoli")
ecotrim<-seasonality("e_CXM","ecoli")
eamp<-seasonality("AMPAMOX","ecoli")

all<-rbind(ecoamox,egent,ecip,etaz,etgc,emer,ecotrim)

ecoli<-ggplot(all) +
  aes(x=day,y=est,group=drug) +
  geom_line() +
  theme_minimal() +
  geom_ribbon(data=all,aes(x=day,ymin=low,ymax=high),alpha=0.1) +
  facet_wrap(~drug) +
  scale_x_continuous(
    breaks=seq(15,346,60),
    labels=month.abb[seq(1,12,2)]
  ) +
  xlab("Month") + ylab("s(Day)") +
  ggtitle("E. coli")

kcoamox<-seasonality("e_AMOXCLAV","kleb")
kgent<-seasonality("e_GENT","kleb")
kcip<-seasonality("e_CIP","kleb")
ktaz<-seasonality("e_PIPTAZ","kleb")
ktgc<-seasonality("threegc","kleb")
kmer<-seasonality("e_MER","kleb")
kcotrim<-seasonality("e_CXM","kleb")


all<-rbind(kcoamox,kgent,kcip,ktaz,ktgc,kmer,kcotrim)
kleb<-ggplot(all) +
  aes(x=day,y=est,group=drug) +
  geom_line() +
  theme_minimal() +
  geom_ribbon(data=all,aes(x=day,ymin=low,ymax=high),alpha=0.1) +
  facet_wrap(~drug) +
  scale_x_continuous(
    breaks=seq(15,346,60),
    labels=month.abb[seq(1,12,2)]
  ) + theme(legend.position = "none") + xlab("Month") + ylab("s(Day)") +
  ggtitle("Klebsiella pneumoniae")


ecoli+ kleb + plot_layout(guides="collect")
###########seasonality by HE############################

seasonal_he<-function(drug,bug){

formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25) + s(day,by=healthcare_exposure, bs='cc', k=12) + healthcare_exposure"))

data_<-filter(data,species==bug)
cat("\nfitting model")
model<-bam(formula,
           data=data_,
           family = binomial(link = "logit"),
           method = "REML")



cat("\nSmooth samples")
checkc<-smooth_samples(model,select="s(day):healthcare_exposurecommunity",n = 1000,n_vals = 365)
checkqc<-smooth_samples(model,select="s(day):healthcare_exposurequasi-community",n=1000,n_vals = 365)
checkn<-smooth_samples(model,select="s(day):healthcare_exposurenosocomial",n=1000,n_vals = 365)
checkqn<-smooth_samples(model,select="s(day):healthcare_exposurequasi-nosocomial",n=1000,n_vals = 365)
check<-rbind(checkc,checkqc,checkn,checkqn)

cat("\nCredibility intervals")
check<-check %>%
  group_by(healthcare_exposure,.row) %>%
  summarise(low=quantile(.value,probs = c(0.025)),
            high=quantile(.value,probs=c(0.975)),
            est=quantile(.value,probs=c(0.5)))

check$drug<-drug
check$bug<-bug
return(check)
}

gent_k<-seasonal_he("e_GENT","kleb")
gent_e<-seasonal_he("e_GENT","ecoli")
coamox_k<-seasonal_he("e_AMOXCLAV","kleb")
coamox_e<-seasonal_he("e_AMOXCLAV","ecoli")
piptaz_k<-seasonal_he("e_PIPTAZ","kleb")
piptaz_e<-seasonal_he("e_PIPTAZ","ecoli")
tgc_k<-seasonal_he("threegc","kleb")
tgc_e<-seasonal_he("threegc","ecoli")
cip_k<-seasonal_he("e_CIP","kleb")
cip_e<-seasonal_he("e_CIP","ecoli")
cotrim_e<-seasonal_he("e_CXM","ecoli")
cotrim_k<-seasonal_he("e_CXM","kleb")
amox_e<-seasonal_he("AMPAMOX","ecoli")

plot<-function(data){
  data_e <- get(paste0(data, "_e"))
  if(data!='amox'){
  data_k <- get(paste0(data, "_k"))
  }
  p1<-ggplot(data_e) +
    aes(x=.row,y=est,color=healthcare_exposure,group=healthcare_exposure) +
    geom_line() +
    geom_ribbon(data=data_e,aes(x=.row,ymin=low,ymax=high,fill=healthcare_exposure),alpha=0.1) +
    theme_minimal() +
    facet_wrap(~healthcare_exposure)+ labs(color="Healthcare exposure",fill="Healthcare exposure") +
    xlab("Day")
  if(data!="amox"){
  p2<-ggplot(data_k) +
    aes(x=.row,y=est,color=healthcare_exposure,group=healthcare_exposure) +
    geom_line() +
    geom_ribbon(data=data_k,aes(x=.row,ymin=low,ymax=high,fill=healthcare_exposure),alpha=0.1) +
    theme_minimal() +
    facet_wrap(~healthcare_exposure) + labs(color="Healthcare exposure",fill="Healthcare exposure") +
    xlab("Day")
  return(p1 + p2)
  }else{
    return(p1)
  }
}



plot("cip") + plot_annotation("Ciprofloxacin") + plot_layout(guides="collect")
plot("gent") +   plot_annotation("Gentamicin") + plot_layout(guides="collect")
plot("piptaz") +   plot_annotation("Piptaz") + plot_layout(guides="collect")
plot("coamox") +   plot_annotation("Co-amoxiclav") + plot_layout(guides="collect")
plot("tgc") +   plot_annotation("Third-gen Cephalosporin") + plot_layout(guides="collect")
plot("cotrim") + plot_annotation("Co-trimoxazole") + plot_layout(guides="collect")
plot("amox")

plot_p4_R<-function(code,bug){
  code <- sym(code)
  data_<-filter(data,species==bug)
  prop <- data %>%
    group_by(year, monthno,healthcare_exposure) %>%
    summarise(S=sum(!!code == 0, na.rm = TRUE),
              R=sum(!!code == 1, na.rm = TRUE))

  prop$date<-ym(paste0(prop$year,prop$monthno))
  prop<-prop %>% pivot_longer(cols = 4:5,names_to = "which",values_to = "n")

  prop$healthcare_exposure<-factor(
    prop$healthcare_exposure,levels=c("community","quasi-community","quasi-nosocomial","nosocomial")
  )

  col<-c("#470000",   "#700404FD", "#AB0707FD", "#FA0808FD" ,"#001F04"  , "#00520B" ,  "#009E18",   "#00E827")
  col<-rev(col)
  prop$which<-factor(prop$which,levels=c("S","R"))
  p4<-ggplot(filter(prop,which=='R' & !is.na(healthcare_exposure))) +
    aes(x=date,y=n,fill=interaction(healthcare_exposure,which)) +
    geom_bar(position="stack", stat="identity") +
    scale_x_date(breaks = seq(as.Date("2012-01-01"), as.Date("2024-01-01"), by="1 year"),
                 date_labels = "%Y") +
    scale_fill_manual(values=col) +
    theme_minimal() + labs(fill="") +
    ylab("Number of isolates")
  return(p4)
}

plot_p4_R("e_CIP","ecoli")
plot_p4_R("e_PIPTAZ","ecoli")
plot_p4_R("e_GENT","ecoli")

i_continuous_day<-function(drug,bug,variable,min,max,step){

    formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25) + s(day, bs='cc', k=12) +
                                s(",variable,",bs='tp',k=25) +
                                ti(",variable,",day,bs=c('tp','cc'))"))

  data_<-filter(data,species==bug)
  cat("\nfitting model")
  model<-bam(formula,
             data=data_,
             family = binomial(link = "logit"),
             method = "REML")

  cat("\nmaking dummy data")
  times<-seq(0,max(data$timeline,na.rm=T),14)

  variable_<-seq(min,max,step)
  predict_data<-expand_grid(timeline=times,variable=variable_)
  predict_data$date<-min(data$spec_dt,na.rm=T) + days(predict_data$timeline)
  predict_data$day<-yday(predict_data$date)
  y=sym(variable)
  predict_data<-rename(predict_data, !!y := variable)
  x=paste0("ti(",variable,",day)")
  cat("\nTaking smooth samples")
  samples<-smooth_samples(model,select = x,n = 1000,data = predict_data)

  cat("\nCredibility intervals")
  samples<-samples %>%
    group_by(.row) %>%
    summarise(low=quantile(.value,probs=c(0.025)),
              est=quantile(.value,probs=c(0.5)),
              high=quantile(.value,probs=c(0.975)))
  samples$date<-predict_data$date
  samples[[variable]]<-predict_data[[variable]]
  cat("\nDone!")
  return(samples)
}
test<-i_continuous_day("e_PIPTAZ","ecoli","age",18,90,5)
test$day<-yday(test$date)
ggplot(test) +
   aes(x=day,y=est,color=age,group=age) +
   geom_line() +
   geom_ribbon(data=test,aes(ymin=low,ymax=high,fill=age),alpha=0.1) + theme_minimal() + facet_wrap(~age)


###################seasonality over time#################

i_continuous_day<-function(drug,bug){

  formula <- as.formula(paste(drug, "~ s(timeline, bs='ad', m=3, k=25) + s(day, bs='cc', k=12) +
                                ti(timeline,day,bs=c('tp','cc'))"))

  data_<-filter(data,species==bug)
  cat("\nfitting model")
  model<-bam(formula,
             data=data_,
             family = binomial(link = "logit"),
             method = "REML")

  cat("\nmaking dummy data")
  times<-seq(0,max(data$timeline,na.rm=T),1)


  predict_data<-expand_grid(timeline=times)
  predict_data$date<-min(data$spec_dt,na.rm=T) + days(predict_data$timeline)
  predict_data$day<-yday(predict_data$date)

  cat("\nTaking smooth samples")
  samples<-smooth_samples(model,select = x,n = 1000,data = predict_data)

  cat("\nCredibility intervals")
  samples<-samples %>%
    group_by(.row) %>%
    summarise(low=quantile(.value,probs=c(0.025)),
              est=quantile(.value,probs=c(0.5)),
              high=quantile(.value,probs=c(0.975)))
  samples$date<-predict_data$date
  samples$drug<-drug
  samples$bug<-bug

  cat("\nDone!")
  return(samples)
}
epiptaz<-i_continuous_day("e_PIPTAZ","ecoli")
egent<-i_continuous_day("e_GENT","ecoli")

all<-rbind(epiptaz,egent)
all$day<-yday(all$date)
all$year<-year(all$date)
all$year<-as.character(all$year)
ggplot(all) +
  aes(x=day,y=est,color=drug,group=drug) +
  geom_line()+ facet_wrap(~year) +
  geom_ribbon(data=all,aes(ymin=low,ymax=high,fill=drug,group=drug),alpha=0.1) + theme_minimal()  +
  scale_x_continuous(
    breaks=seq(15,346,60),
    labels=month.abb[seq(1,12,2)]
  )

coamox + piptaz


#########space time################

function_st_he<-function(drug,bug,seasonal){
  if(seasonal==TRUE){
    cat("\nfitting a seasonal model")
    formula <- as.formula(paste(drug, "~ s(timeline, bs='tp', m=3, k=25) + s(day, bs='cc', k=12) +
                                s(nuts_name, bs = 'mrf', xt = list(nb = nb)) + healthcare_exposure +
                                ti(nuts_name, timeline,healthcare_exposure, bs = c('mrf','tp','re'),
                                xt = list(nb = nb))"))
  }else{
    cat("\nfitting a non-seasonal model")
    formula <- as.formula(paste(drug, "~ s(timeline, bs='tp', m=3, k=25)  +
                                s(nuts_name, bs = 'mrf', xt = list(nb = nb)) + healthcare_exposure +
                                ti(nuts_name, timeline,healthcare_exposure, bs = c('mrf','tp','re'),
                                xt = list(nb = nb))"))

  }

  data_<-filter(d,species==bug )
  cat("\nfitting model")
  model<-bam(formula,
             data=data_,
             family = binomial(link = "logit"),
             method = "REML")


  nuts<-unique(d$nuts_name)
  time<-seq(0,max(d$timeline),28)
  he<-unique(data_$healthcare_exposure)
  he<-he[!is.na(he)]
  new_data<-expand_grid(nuts_name=nuts,timeline=time,healthcare_exposure=he)
  new_data$date<-min(d$spec_dt,na.rm = T) + days(new_data$timeline)
  new_data$day<-yday(new_data$date)

  cat("\nPredicting from model")
  m<-fitted_values(model,data = new_data)
  m$bug<-bug
  m$drug<-drug


  cat("\nSamples of marginal effects from posterior")
  samples<-smooth_samples(model,select = "ti(nuts_name,timeline,healthcare_exposure)",n = 1000,data = new_data)

  samples<-samples %>%
    group_by(.row) %>%
    summarise(low=quantile(.value,probs=c(0.025)),
              est=quantile(.value,probs=c(0.5)),
              high=quantile(.value,probs=c(0.975)))

  cat("\nTaking derivative samples")
  deriv<-derivative_samples(model,focal="timeline",order=1,scale="response",data=new_data,draws = 1000)
  cat("\nCreating credibility intervals")
  deriv<-deriv %>% group_by(.row) %>%
    summarise(low=quantile(.derivative,probs = c(0.025)),
              high=quantile(.derivative,probs=c(0.975)),
              est=quantile(.derivative,probs=c(0.5)))

  
  deriv$bug<-bug
  deriv$drug<-drug
  deriv$date<-new_data$date
  deriv$nuts_name<-new_data$nuts_name
  deriv$healthcare_exposure<-new_data$healthcare_exposure
  deriv_n<-filter(deriv,healthcare_exposure=='nosocomial')
  deriv_n<-annotate_derivatives(deriv_n)
  deriv_qn<-filter(deriv,healthcare_exposure=='quasi-nosocomial')
  deriv_qn<-annotate_derivatives(deriv_qn)
  deriv_qc<-filter(deriv,healthcare_exposure=='quasi-community')
  deriv_qc<-annotate_derivatives(deriv_qc)
  deriv_c<-filter(deriv,healthcare_exposure=='community')
  deriv_c<-annotate_derivatives(deriv_c)
  deriv<-rbind(deriv_n,deriv_qn,deriv_qc,deriv_c)
  
  samples$drug<-drug
  samples$bug<-bug
  samples$nuts_name<-new_data$nuts_name
  samples$date<-new_data$date
  samples$healthcare_exposure<-new_data$healthcare_exposure
  return(list(fitted=m,deriv=deriv,pe=samples))
}

coamox_ecoli_sthe<-function_st_he("e_AMOXCLAV","ecoli",TRUE)
coamox_kleb_sthe<-function_st_he("e_AMOXCLAV","kleb",TRUE)
gent_ecoli_sthe<-function_st_he("e_GENT","ecoli",TRUE)
gent_kleb_sthe<-function_st_he("e_GENT","kleb",FALSE)
piptaz_ecoli_sthe<-function_st_he("e_PIPTAZ","ecoli",TRUE)
piptaz_kleb_sthe<-function_st_he("e_PIPTAZ","kleb",FALSE)
cotrim_ecoli_sthe<-function_st_he("e_CXM","ecoli",TRUE)
cotrim_kleb_sthe<-function_st_he("e_CXM","kleb",FALSE)
tgc_ecoli_sthe<-function_st_he("threegc","ecoli",TRUE)
tgc_kleb_sthe<-function_st_he("threegc","kleb",FALSE)
mero_ecoli_sthe<-function_st_he("e_MER","ecoli",TRUE)
mero_kleb_sthe<-function_st_he("e_MER","kleb",FALSE)
cip_ecoli_sthe<-function_st_he("e_CIP","ecoli",TRUE)
cip_kleb_sthe<-function_st_he("e_CIP","kleb",FALSE)

sthe_data_frames <- ls(pattern = "*_*_sthe")
age <- map_dfr(age_data_frames, get)


plot<-function(drug,bug){
  fitted<-get(paste0(drug,'_',bug))$fitted
  deriv<-get(paste0(drug,'_',bug))$deriv
  pe<-get(paste0(drug,'_',bug))$pe
  deriv<-left_join(deriv,fitted,by=c("date","healthcare_exposure","nuts_name"))
  deriv$direction<-ifelse(deriv$est < 0,'reducing','increasing')
  p1a<- ggplot(fitted) +
    aes(x=date,y=.fitted,color=healthcare_exposure) +
    geom_line() +
    #geom_ribbon(data=fitted,aes(x=date,ymin=.lower_ci,ymax=.upper_ci,fill=healthcare_exposure),alpha=0.1) +
    geom_segment(data=filter(deriv,first_breakpoint==1),aes(x=date,xend=(date + days(20)),y=.fitted,group=healthcare_exposure,color=direction),linewidth=5,alpha=1) +
    ylab("Proportion of isolates resistant") +
    theme_minimal() +
    theme(legend.position = "none",axis.text.x=element_text(angle=45)) +

    scale_x_date(breaks = seq(as.Date("2012-01-01"), as.Date("2024-01-01"), by="2 years"),
                 date_labels = "%Y") +
    xlab("") +
    labs(color="Region") +
    facet_wrap(~nuts_name)
  p1b<- ggplot(pe) +
    aes(x=date,y=est,color=healthcare_exposure) +
    geom_line() +
    #geom_ribbon(data=pe,aes(x=date,ymin=low,ymax=high,fill=healthcare_exposure),alpha=0.1) +
    facet_wrap(~nuts_name) +
    theme_minimal()  +
    theme(axis.text.x = element_text(angle=45))  +
    scale_x_date(breaks = seq(as.Date("2012-01-01"), as.Date("2024-01-01"), by="2 years"),
                 date_labels = "%Y")


  return(p1a + p1b )
}
plot("coamox","kleb")

ggplot()

function_he2<-function(drug,bug,seasonal){
  if(seasonal==TRUE){
    cat("\nfitting a seasonal model")
    formula <- as.formula(paste(drug, "~ s(timeline, by =healthcare_exposure, bs='tp', m=3, k=25) + s(day, bs='cc', k=12) +
                                 healthcare_exposure"))
  }else{
    cat("\nfitting a non-seasonal model")
    formula <- as.formula(paste(drug, "~ s(timeline, bs='tp', m=3, k=25)  +
                                s(nuts_name, bs = 'mrf', xt = list(nb = nb)) +
                                ti(nuts_name, timeline, bs = c('mrf','tp'))"))
    
  }
  
  data_<-filter(d,species==bug )
  cat("\nfitting model")
  model<-bam(formula,
             data=data_,
             family = binomial(link = "logit"),
             method = "REML")
  
  
  
  time<-seq(0,max(d$timeline),28)
  he<-unique(data_$healthcare_exposure)
  
  new_data<-expand_grid(timeline=time,healthcare_exposure=he)
  new_data$date<-min(d$spec_dt,na.rm = T) + days(new_data$timeline)
  new_data$day<-yday(new_data$date)
  
  cat("\nPredicting from model")
  m<-fitted_values(model,data = new_data)
  m$bug<-bug
  m$drug<-drug
  

  
  cat("\nTaking derivative samples")
  deriv<-derivative_samples(model,focal="timeline",order=1,scale="response",data=new_data,draws = 1000)
  cat("\nCreating credibility intervals")
  deriv<-deriv %>% group_by(.row) %>%
    summarise(low=quantile(.derivative,probs = c(0.025)),
              high=quantile(.derivative,probs=c(0.975)),
              est=quantile(.derivative,probs=c(0.5)))
  
  
  
  deriv$bug<-bug
  deriv$drug<-drug
  deriv$date<-new_data$date
  deriv$healthcare_exposure<-new_data$healthcare_exposure
  deriv_n<-filter(deriv,healthcare_exposure=='nosocomial')
  deriv_n<-annotate_derivatives(deriv_n)
  deriv_qn<-filter(deriv,healthcare_exposure=='quasi-nosocomial')
  deriv_qn<-annotate_derivatives(deriv_qn)
  deriv_qc<-filter(deriv,healthcare_exposure=='quasi-community')
  deriv_qc<-annotate_derivatives(deriv_qc)
  deriv_c<-filter(deriv,healthcare_exposure=='community')
  deriv_c<-annotate_derivatives(deriv_c)
  deriv<-rbind(deriv_n,deriv_qn,deriv_qc,deriv_c)
  
  
  return(list(fitted=m,deriv=deriv))
}

coamox_ecoli<-function_he2("e_AMOXCLAV","ecoli",TRUE)


plot<-function(drug,bug){
  fitted<-get(paste0(drug,'_',bug))$fitted
  deriv<-get(paste0(drug,'_',bug))$deriv
  pe<-get(paste0(drug,'_',bug))$pe
  deriv<-left_join(deriv,fitted,by=c("date","healthcare_exposure"))
  deriv$direction<-ifelse(deriv$est < 0,'reducing','increasing')
  p1a<- ggplot(fitted) +
    aes(x=date,y=.fitted,color=healthcare_exposure) +
    geom_line() +
    geom_ribbon(data=fitted,aes(x=date,ymin=.lower_ci,ymax=.upper_ci,fill=healthcare_exposure),alpha=0.1) +
    geom_segment(data=filter(deriv,first_breakpoint==1),aes(x=date,xend=(date + days(20)),y=.fitted,group=healthcare_exposure,color=direction),linewidth=5,alpha=1) +
    ylab("Proportion of isolates resistant") +
    theme_minimal() +
    #theme(legend.position = "none",axis.text.x=element_text(angle=45)) +
    
    scale_x_date(breaks = seq(as.Date("2012-01-01"), as.Date("2024-01-01"), by="2 years"),
                 date_labels = "%Y") +
    xlab("") +
    labs(color="Region")
  p1b<- ggplot(pe) +
    aes(x=date,y=est,color=healthcare_exposure) +
    geom_line() +
    #geom_ribbon(data=pe,aes(x=date,ymin=low,ymax=high,fill=healthcare_exposure),alpha=0.1) +
    theme_minimal()  +
    #theme(legend.position = "none",axis.text.x = element_text(angle=45))  +
    scale_x_date(breaks = seq(as.Date("2012-01-01"), as.Date("2024-01-01"), by="2 years"),
                 date_labels = "%Y")
  
  
  return(p1a + p1b)
}
plot("coamox","ecoli")
