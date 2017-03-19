
# last modification: sep2016 ----------------------------------------------------------------------------------------------------------------------------------

```{r descript0, echo=TRUE}
#library("haven")

#PQ_main<-read_stata("Primaquine base 05apr2014-Meddly_Julio.dta") #----> CONSERVE THE LABELS / AFTER write.csv THEY ARE GONE!
SV_main<-read_stata("severe_malaria_vivax SEP2016.dta")             #----> WITHOUT LABELS
dim(SV_main)


#
SV_main$codigo <- gsub('LIM0','LIM',SV_main$codigo) #---> http://stackoverflow.com/questions/7641666/changing-text-in-a-data-frame
#
#write.csv(PQ_main, file="PQ_05apr2014-Meddly_Julio.csv")
write.csv(SV_main, file="SEVERE_sep2016.csv")
#
tabAll <- read.table("SEVERE_sep2016.csv", row.names=2, header=TRUE, sep=",")
##

#sevData<-read.table("SEVERA_covariates_224samples_EXTENDED_18ago2016.csv",row.names=1,header=TRUE, sep=",", 
#                        colClasses = c(rep("factor",5),"numeric","factor","numeric","numeric",rep("factor",5),"numeric","numeric",
#                                       rep("factor",5),rep("numeric",9),"factor","numeric",rep("factor",2),rep("numeric",2)))
# ABOVE replaced by BELOW "SELECTEDcovarariates"

SELECTEDcovarariates <- c("edad", 
                          "sexo",
                          
                          "admin_uci_corr", "cns_malaria", 
                          "dias_uci", 
                          "postracion_corr", "confusion_corr", "coma_corr", 
                          "glas_corr", 
                          "convulsion_corr", #"num_convulsion", 
                          "sdra", "lung_injury", 
                          "disnea_corr", "taquipnea_corr", "infiltracion_corr", "edema_pulm_corr", 
                          "ictericia_corr", "hiperbil", 
                          "bilirrubina", 
                          "shock", "bleeding", "gingivorragia_corr", "epistaxis_corr", "sang_intest_corr", 
                          "trombocitopenia", 
                          "cont_plaquetas", 
                          "glic_normal", #"hipo_glic1", "hipo_glic2", 
                          "hiperparasitemia", #"erit_asex", "cien_paras", 
                          "hiperpirexia", 
                          "hemoglobinuria_corr", "anemia_severa", 
                          "oliguria_corr", 
                          "temperatura", 
                          "fcardiaca", "frespiratoria", "parterial", 
                          "talla", "peso",
                          
                          "malestar_us_corr", "fiebre_us_corr", "escalosfrios_us_corr", 
                          "mialgia_us_corr", "artralgia_us_corr", "cefalea_us_corr", 
                          "diarrea_us_corr", "vomitos_us_corr", "dif_resp_us_corr", "incons_us_corr",
                          
                          "malestar_hoy_corr", "fiebre_hoy_corr", "escalosfrios_hoy_corr", 
                          "mialgia_hoy_corr", "artralgia_hoy_corr", "cefalea_hoy_corr", 
                          "diarrea_hoy_corr", "vomitos_hoy_corr", "dif_resp_hoy_corr", "inconciencia_hoy_corr",
                          
                          "malaria_pasado",
                          "malaria_pasado_num", 
                          
                          "presenta_comorb", 
                          "diabetes", "enf_snc", "enf_pulm",
                          "enf_cardiac", "enf_renal", "enf_otra"
                          ) # NOW able to ADD whatever you want!

tabSELECT <- tabAll[,SELECTEDcovarariates]
#tabSELECT <- tabAll[,colnames(sevData)]
#tabSELECT <- SV_main[,colnames(sevData)] #--------> On dplyr format
#sapply(tabSELECT, class)

#subset(sevData,select = c(cns_malaria,postracion_corr:convulsion_corr))
#subset(sevData,select = c(sdra:edema_pulm_corr,frespiratoria))     #LIM1076 not SDRA?????
#subset(sevData,select = c(ictericia_corr:bilirrubina))
#subset(sevData,select = c(shock:sang_intest_corr,fcardiaca,parterial))
#subset(sevData,select = c(hiperpirexia,temperatura,fiebre_us_corr,fiebre_hoy_corr))    # DIFF us or hoy?????
#subset(sevData,select = c(trombocitopenia,cont_plaquetas))
#subset(sevData,select = c(hemoglobinuria_corr,anemia_severa))


###################
write.csv(tabSELECT, file="SEVERE_sep2016_SELECTEDcovariates.csv")

sevData<-read.table("SEVERE_sep2016_SELECTEDcovariates.csv",row.names=1,header=TRUE, sep=",", 
                    colClasses = c("factor", #toma en cuenta los row.names
                                   "numeric",
                                   rep("factor",3),
                                   "numeric",
                                   rep("factor",3),
                                   "numeric",
                                   rep("factor",9),
                                   "numeric",
                                   rep("factor",6),
                                   "numeric",
                                   rep("factor",6),
                                   rep("numeric",3),
                                   "factor",
                                   rep("numeric",2),
                                   rep("factor",21),
                                   "numeric",
                                   rep("factor",7)))
#sapply(sevData, class)
###################

#subset to SEVERE Sample Names
eset<-normal.log2.raw.eSet.VIVAX.SEV.FILTER
sampleSEV <- sevData[sampleNames(eset),]
#
#sapply(sampleSEV, class)
```

