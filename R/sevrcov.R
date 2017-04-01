
# 30mar2017 -----------------------------------------------------------------------------------------------
##essential
library(Hmisc)
library("knitr")      # To display nice tables
library("tidyverse")
#ggplot2, tibble, tidyr, readr, purr, dplyyr
##extra
library("Rmisc")      # multiplot ggplots
library("haven")      # import STATA df
library("forcats")    # factor vectors
library("stringr")    # for strings
library("viridis")    # color visualization

# INPUT ----------------------------------------------------------------------------------------

sevdb <- haven::read_dta("data-raw/severe_malaria_vivax SEP2016.dta") #%>% haven::as_factor()
#glimpse(sevdb)
#glimpse(sevdb %>% select(-ends_with("_corr")))
#glimpse(sevdb %>% select(ends_with("_us_corr")))
#glimpse(sevdb %>% select(ends_with("_hoy_corr")))
#Hmisc::contents(sevdb)

sevdb %>% select(contains("fiebre"),temperatura) # inconsistent

# NO creatinina -> ???
# NO parasitemia -> hiperparasitemia
sevdb_1 <- sevdb %>% 
  select(codigo, edad, sexo,
         eritro_parasit, admin_uci_corr:temperatura, 
         ends_with("_corr"), -ends_with("_us_corr"), -ends_with("_hoy_corr"), 
         presenta_comorb, ends_with("_riesgo"), starts_with("malaria_"),
         -glas_ojos, -glas_verbal, -glas_motor, -glas_total, -num_convulsion, 
         -bilirrubina, -cont_plaquetas,
         -starts_with("droga"), -starts_with("tx_")
         ) %>% 
  mutate(codigo= str_replace(codigo,"LIM0","LIM"),
         glas_11= if_else(glas_corr<11 & glas_corr>0,1,0),
         shock= str_replace(shock,"99",replacement = NA_character_),
         hemoglobinuria_corr= str_replace(hemoglobinuria_corr,"99",replacement = NA_character_),
         sdra= str_replace(hemoglobinuria_corr,"NaN",replacement = NA_character_),
         hiperbil= str_replace(hemoglobinuria_corr,"NaN",replacement = NA_character_)) %>% 
  mutate(shock= as.double(shock),
         hemoglobinuria_corr= as.double(hemoglobinuria_corr),
         sdra= as.double(sdra),
         hiperbil= as.double(hiperbil),
         malaria_pasado= as.factor(malaria_pasado)
         ) %>% 
  mutate(malaria_pasado= fct_recode(malaria_pasado, "expuesto"="1", "no-expuesto"="0"))
glimpse(sevdb_1)

sevdb_1 %>% dplyr::count(malaria_pasado)

#Hmisc::contents(sevdb_1)

# CNS malaria = postracion, coma, convulcion, confusion
sevdb_1 %>% group_by(cns_malaria,postracion_corr, coma_corr,convulsion_corr,confusion_corr) %>% dplyr::count()
sevdb_1 %>% group_by(cns_malaria,glas_corr, glas_11) %>% dplyr::count()

# SDRA != LUNG_INJURY = at least 2 criteria of lung related-injuries
sevdb_1 %>% group_by(lung_injury,disnea_corr,taquipnea_corr,infiltracion_corr,edema_pulm_corr) %>%  dplyr::count()

# DIFFERENT VARIABLES
sevdb_1 %>% group_by(ictericia_corr,hiperbil) %>%  dplyr::count()             #ictericia!=hiperbili
sevdb_1 %>% group_by(anemia_severa, hemoglobinuria_corr) %>%  dplyr::count()  #anemia!=hemoglobinuria
sevdb_1 %>% group_by(glic_normal,hipo_glic1,hipo_glic2) %>%  dplyr::count()   #
#sevdb_1 %>% group_by(admin_uci_corr,dias_uci) %>% dplyr::count()

# SEVERE CRITERIA ----------------------------------------------------------------------------------------------------------
sevdb_2 <- sevdb_1 %>% 
  select(1:4,contains("uci"),
         cns_malaria,
         
         glas_11,
         sdra,
         lung_injury,
         ictericia_corr,hiperbil,
         shock,
         bleeding,
         glic_normal, #contains("glic") = hipo_glic1 + hipo_glic2
         hiperparasitemia,
         anemia_severa, hemoglobinuria_corr,
         
         presenta_comorb,
         malaria_pasado, malaria_pasado_num
         ) %>% 
  mutate(sev_WHO_num= rowSums(.[8:18],na.rm = T),
         sev_WHO= if_else(sev_WHO_num>0,1,0)) %>% 
  mutate(sev_WHO= as.factor(sev_WHO)) %>% 
  mutate(sev_WHO= fct_recode(sev_WHO, "no-severo"="0","severo"="1"))
glimpse(sevdb_2)

Hmisc::contents(sevdb_2)

sevdb_2 %>% dplyr::count(malaria_pasado)
sevdb_2 %>% dplyr::count(sev_WHO)
sevdb_2 %>% group_by(malaria_pasado, sev_WHO) %>% dplyr::count()

#sevdb_2 %>% group_by(cns_malaria,glas_11) %>% dplyr::count()

# abnomic SUBSET ---------------------------------------------------------------------------------------------------------
subsev <- data_frame(
  codigo= c( "LIM2048", "LIM1025", "LIM1029", "LIM2035", "LIM1043", "LIM2017", "LIM1054", 
             "LIM2102", "LIM2010", "LIM1038", "LIM2050", "LIM2115", "LIM2043", "LIM1004", 
             "LIM2030", "LIM2071", "LIM1035", "LIM2104", "LIM1051", "LIM2107", "LIM1062", 
             "LIM2072", "LIM2074", "LIM1056", "LIM1050", "LIM2110", "LIM2026", "LIM1047", 
             "LIM1057", "LIM1016", "LIM1063", "LIM1044", "LIM2044", "LIM1037", "LIM2007", 
             "LIM2086", "LIM2042", "LIM2101", "LIM1045", "LIM2085", "LIM1048", "LIM2066", 
             "LIM2067", "LIM1021", "LIM2077", "LIM1061", "LIM2056", "LIM2034", "LIM2016", 
             "LIM2111", "LIM1030", "LIM2096", "LIM1064", "LIM2093", "LIM1041", "LIM2040", 
             "LIM1052", "LIM2013", "LIM1059", "LIM2089"),
  abnomic= TRUE
)

sevdb_3 <- dplyr::inner_join(sevdb_2,subsev)
Hmisc::contents(sevdb_3)
sevdb_3 %>% dplyr::count(malaria_pasado)
sevdb_3 %>% dplyr::count(sev_WHO)
sevdb_3 %>% group_by(malaria_pasado, sev_WHO) %>% dplyr::count()


#sevdb_3 %>% group_by(admin_uci_corr,dias_uci) %>% dplyr::count()
#sevdb_3 %>% group_by(sdra,lung_injury) %>%  dplyr::count()
#sevdb_3 %>% group_by(ictericia_corr,hiperbil) %>%  dplyr::count()
#sevdb_3 %>% group_by(anemia_severa, hemoglobinuria_corr) %>%  dplyr::count()

# SUM SEVERITY CRITERIAS
sevdb_3 %>% 
  group_by(#cns_malaria,
           glas_11,
           sdra,
           lung_injury,
           ictericia_corr,hiperbil,
           shock,
           bleeding,
           glic_normal,
           hiperparasitemia,
           anemia_severa, hemoglobinuria_corr
           ) %>% 
  dplyr::count() #%>% View()


## HMISC summaries ------------------------------------------------------------------------------------------------------------
sevdb_4 <- sevdb_3 %>% haven::as_factor()
sevdb_4 %>% Hmisc::contents()

s1 <- Hmisc::summaryM(edad + sexo + eritro_parasit + presenta_comorb ~ malaria_pasado,
               data=sevdb_4,
               overall=FALSE, test=TRUE)
plot(s1, which='categorical')
Key(0,1)
plot(s1, which='continuous')

latex(s1, caption='Clinical data',
      exclude1=TRUE, npct='both', 
      digits=3,
      prmsd=TRUE, brmsd=TRUE, msdsize=mu$smaller2, #NOT-EVALUATE if PDF
      middle.bold=TRUE, long = TRUE,
      #legend.bottom = TRUE, #insert.bottom = TRUE, 
      what="%", html = TRUE
      ) #change here for LaTeX PDF

# -----------------------------------------------------------------------------------------------

