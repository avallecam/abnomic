library(tidyverse)

# importar ag reactivos ---------------------------------------------------

td <- read_rds("data/06-reactive_compare-all_reactives.rds")
td %>% count(reactivity)

pvx_sal1_snp_2017 <- read_tsv("data/04-listgen.tsv") %>% 
  janitor::clean_names() %>% 
  select(gene_id,product_description,
         total_sn_ps_all_strains,
         non_syn_syn_snp_ratio_all_strains,
         transcript_length) %>% 
  filter(str_detect(gene_id,"PVX"))

pvx_sal1_snp_2017

# importar referencia -----------------------------------------------------

## _p01 -----
pvp01 <- read_tsv("data/20211223-pvivax-p01-ortholog_group-more.tsv") %>% 
  janitor::clean_names() %>% 
  # glimpse()
  mutate(number_tm_domains=as.numeric(number_tm_domains)) %>% 
  select(-organism,-x13) %>% 
  select(-signal_p_scores,-number_tm_domains) %>% 
  select(-ec_numbers_from_ortho_mcl) %>% 
  select(-ortholog_count,-paralog_count) %>%
  # select(-ortholog_group) %>% 
  rename(gene_id_p01=gene_id,
         source_id_p01=source_id) 

## _sal1 ----
pvsal <- read_tsv("data/20211223-pvivax-sal1-ortholog_group-more.tsv") %>% 
  janitor::clean_names() %>% 
  mutate(number_tm_domains=as.numeric(number_tm_domains)) %>% 
  select(-organism,-x13) %>% 
  select(-signal_p_scores,-number_tm_domains) %>% 
  select(-ec_numbers_from_ortho_mcl) %>% 
  select(-total_sn_ps_all_strains,
         -non_syn_syn_snp_ratio_all_strains) %>% 
  select(-ortholog_count,-paralog_count) %>%
  # select(-ortholog_group) %>% 
  rename(gene_id_sal=gene_id,
         source_id_sal=source_id)

pvp01 %>% glimpse()
pvsal %>% glimpse()

pvp01 %>% 
  filter(str_detect(ortholog_group,"OG6_149563"))

# unir referencia a reactivos ---------------------------------------------

unir_sal <- td %>% 
  # solo ver antigenos y estatus de reactividad
  select(id=ID,reactivity) %>% 
  # solo matchear a los antigenos reactivos
  filter(reactivity=="reactive") %>% 
  # extraer codigos de genes del nombre del antigeno
  mutate(gene_id_sal = str_replace(id,"(PVX_.{6}).+","\\1")) %>% 
  # reemplazar por aliases conocidos
  mutate(gene_id_sal = case_when(
    # https://plasmodb.org/plasmo/app/record/gene/PVX_080690#Alias
    gene_id_sal == "PVX_080695" ~ "PVX_080690",
    # https://plasmodb.org/plasmo/app/record/gene/PVX_001660#Alias
    gene_id_sal == "PVX_001665" ~ "PVX_001660",
    # https://plasmodb.org/plasmo/app/record/gene/PVX_091434#Alias
    str_detect(id, "PVX_091435_6o7") ~ "PVX_091434",
    # https://plasmodb.org/plasmo/app/record/gene/PVX_091436#Alias
    str_detect(id,"PVX_091435_1o7_S2") ~ "PVX_091436",
    TRUE ~ gene_id_sal
  )) %>% 
  # unir con base de sal1 para obtener los codigos de grupos ortologos
  left_join(pvsal)

# verificar missings
# unir_sal %>%
  # naniar::miss_var_summary()
  # filter(is.na(ortholog_group)) %>% 
  # select(1:3)
  # epihelper::print_inf()

# verificar que solo haya un ortholog_group en base left

#33
sal_ortholog_vario <- unir_sal %>%   
  count(ortholog_group,sort = T) %>% 
  filter(n>1)

#260
sal_ortholog_unico <- unir_sal %>%   
  count(ortholog_group,sort = T) %>% 
  filter(n==1)

#25
sal_product_vario <- unir_sal %>%   
  count(product_description,sort = T) %>% 
  filter(n>1)

#139
sal_product_unico <- unir_sal %>%   
  count(product_description,sort = T) %>% 
  filter(n==1)

sal_ortholog_vario
sal_ortholog_unico
sal_product_vario
sal_product_unico

#' # _test one ortho ---------------------------------------------------------
#' 
#' unir_sal %>%   
#'   count(ortholog_group,sort = T)
#' unir_sal %>%   
#'   count(ortholog_group,product_description,sort = T)
#' unir_sal %>%
#'   filter(ortholog_group=="OG6_100908")
#' 
#' pvp01 %>% 
#'   filter(ortholog_group=="OG6_100908") %>%
#'   select(-source_id_p01) %>% 
#'   # rename(product_description_p01 = product_description) %>%
#'   # rename(ortholog_group_p01 = ortholog_group) %>%
#'   left_join(
#'     unir_sal %>%
#'       filter(ortholog_group=="OG6_100908") %>%
#'       select(-id)
#'   ) %>% 
#'   distinct(gene_id_p01,.keep_all = T) %>%
#'   epihelper::print_inf()
#'   # identity() %>% 
#'   # count(reactivity)
#' 
#' #' input n=17
#' #' input n=05
#' #' 
#' #' resultados de join + distinct
#' #' by product -> 03/17 reactivos
#' #' by ortholog -> 17/17 reactivos
#' #' by product + ortholog -> 03/17  reactivos
#' #' 
#' #' resultados de join (sin distinct)
#' #' by product + ortholog -> 07/21  reactivos

# _pretest ----------------------------------------------------------------

#' plan
#' 
#' primero,
#' unir genes p01 y sal1 
#' a traves de ortholog_group y product_name
#' debido a replicas en sal1
#' retirar genes replicados en p01
#' 
#' segundo,
#' del complemento de genes que no se hayan unido en ambas bases
#' unir genes de p01 y sal1
#' a traves de ortholog_group (unicamente)
#' debido a replicas en sal1
#' retirar genes replicados en p01
#' 
#' tercero,
#' el complemento de ambos conjuntos anteriores
#' al intentar unirlos, solo identifica
#' 02 genes hipoteticos de sal1
#' por lo tanto
#' decidimos no tomar en cuenta esta union
#' y solo usar los valores de snp del 2017
#' 

# __ortho+prod ------------------------------------------------------------

sal_test <- unir_sal %>% 
  select(gene_id_sal,ortholog_group,product_description,reactivity)

p01_test <- pvp01 %>% 
  select(gene_id_p01,ortholog_group,product_description,
         non_syn_syn_snp_ratio_all_strains,total_sn_ps_all_strains)

p01_sal_by_ortho_prod_unique <- 
  p01_test %>% 
  inner_join(sal_test) %>% 
  distinct(gene_id_p01,.keep_all = T)

# n=164 por ortholog + product + distinct
# p01_sal_by_ortho_prod_unique #%>% 
#   # epihelper::print_inf()

p01_sal_by_ortho_prod_unique %>% 
  # count(gene_id_sal,sort = T)
  arrange(gene_id_p01,desc(non_syn_syn_snp_ratio_all_strains)) %>% 
  distinct(gene_id_sal,.keep_all = T) %>%
  # count(gene_id_sal,sort = T) %>% 
  identity()

p01_sal_by_ortho_prod_unique

# __ortho only ------------------------------------------------------------

p01_sal_by_ortho_prod_unique_complement_ortho <- p01_test %>% 
  filter(!magrittr::is_in(gene_id_p01,
                         p01_sal_by_ortho_prod_unique %>% 
                           pull(gene_id_p01))) %>% 
  rename(product_description_p01 = product_description) %>%
  # rename(ortholog_group_p01 = ortholog_group) %>%
  inner_join(
    # 220 subset
    sal_test %>% 
      filter(!magrittr::is_in(gene_id_sal,
                             p01_sal_by_ortho_prod_unique %>% 
                               distinct(gene_id_sal,.keep_all = T) %>%
                               pull(gene_id_sal)))
  ) %>% 
  distinct(gene_id_p01,.keep_all = T) %>%
  # select(-non_syn_syn_snp_ratio_all_strains,-total_sn_ps_all_strains) %>% 
  # epihelper::print_inf()
  # view()
  identity()

p01_sal_by_ortho_prod_unique_complement_ortho %>% 
  distinct(gene_id_sal,.keep_all = T)

p01_sal_by_ortho_prod_unique
p01_sal_by_ortho_prod_unique_complement_ortho

# __prod only -------------------------------------------------------------

#' esta lista de genes no será incluida
#' por que solo pegan por nombre 
#' 02 genes de sal1 a
#' 58 genes de p01

p01_test %>% 
  filter(!magrittr::is_in(gene_id_p01,
                          p01_sal_by_ortho_prod_unique %>% 
                            pull(gene_id_p01))) %>% 
  filter(!magrittr::is_in(gene_id_p01,
                          p01_sal_by_ortho_prod_unique_complement_ortho %>% 
                            pull(gene_id_p01))) %>% 
  # rename(product_description_p01 = product_description) %>%
  rename(ortholog_group_p01 = ortholog_group) %>%
  inner_join(
    sal_test %>% 
      filter(!magrittr::is_in(gene_id_sal,
                              p01_sal_by_ortho_prod_unique %>% 
                                distinct(gene_id_sal,.keep_all = T) %>%
                                pull(gene_id_sal))) %>% 
      filter(!magrittr::is_in(gene_id_sal,
                              p01_sal_by_ortho_prod_unique_complement_ortho %>% 
                                distinct(gene_id_sal,.keep_all = T) %>%
                                pull(gene_id_sal))) %>% 
      # epihelper::print_inf() %>% 
      identity()
  ) %>% 
  distinct(gene_id_p01,.keep_all = T) %>%
  distinct(gene_id_sal,.keep_all = T) %>%
  select(-non_syn_syn_snp_ratio_all_strains,-total_sn_ps_all_strains) %>% 
  epihelper::print_inf()

p01_sal_by_ortho_prod_unique_complement_ortho_complement_sal1 <- 
  sal_test %>% 
  filter(!magrittr::is_in(gene_id_sal,
                          p01_sal_by_ortho_prod_unique %>% 
                            distinct(gene_id_sal,.keep_all = T) %>%
                            pull(gene_id_sal))) %>% 
  filter(!magrittr::is_in(gene_id_sal,
                          p01_sal_by_ortho_prod_unique_complement_ortho %>% 
                            distinct(gene_id_sal,.keep_all = T) %>%
                            pull(gene_id_sal))) %>% 
  # epihelper::print_inf() %>% 
  # identity()
  left_join(
    pvx_sal1_snp_2017,
    by = c("gene_id_sal" = "gene_id",
           "product_description" = "product_description")
  )

p01_sal_by_ortho_prod_unique
p01_sal_by_ortho_prod_unique_complement_ortho
p01_sal_by_ortho_prod_unique_complement_ortho_complement_sal1

# __unite all -------------------------------------------------------------

p01_sal_union_strategy_p01 <- p01_sal_by_ortho_prod_unique %>% 
  union_all(p01_sal_by_ortho_prod_unique_complement_ortho)

p01_sal_union_strategy_sal1 <- 
  p01_sal_by_ortho_prod_unique_complement_ortho_complement_sal1 %>% 
  mutate(non_synonymous_sn_ps_all_strains = 
           total_sn_ps_all_strains - 
           (total_sn_ps_all_strains / (1+non_syn_syn_snp_ratio_all_strains))) %>% 
  mutate(synonymous_sn_ps_all_strains = 
           total_sn_ps_all_strains -non_synonymous_sn_ps_all_strains) %>% 
  # make unique names
  mutate(gene_id_sal=janitor::make_clean_names(gene_id_sal))

p01_sal_union_strategy_p01_sal1 <- p01_sal_union_strategy_p01 %>% 
  union_all(p01_sal_union_strategy_sal1)

p01_sal_union_strategy_p01_sal1 %>% 
  naniar::vis_miss()

p01_sal_union_strategy_p01_sal1 %>% 
  count(reactivity)
td %>% 
  count(reactivity)

p01_sal_union_strategy_sal1 %>% 
  # mutate(gene_id_sal=janitor::make_clean_names(gene_id_sal)) %>% 
  count(gene_id_sal,sort = T)

# # _test ---
# 
#' wrongwayround <- pvp01 %>% 
#'   # left join only by ortholog
#'   # rename(product_description_p01 = product_description) %>%
#'   # rename(ortholog_group_p01 = ortholog_group) %>%
#'   left_join(unir_sal) %>% 
#'   # count(reactivity)
#'   mutate(reactivity=if_else(is.na(reactivity),
#'                             "not_reactive",
#'                             reactivity)) 
#' 
#' wrongwayround %>% 
#'   count(reactivity)
#'   # naniar::vis_miss()
#'   # naniar::miss_var_summary()
#' 
#' wrongwayround %>% 
#'   group_by(gene_id_p01) %>% 
#'   filter(n()>1)
#' 
#' wrongwayround_unique <- wrongwayround %>% 
#'   # arrange(gene_id_p01,desc(non_syn_syn_snp_ratio_all_strains)) %>% 
#'   distinct(gene_id_p01,.keep_all = T) %>%
#'   identity()
#'  
#' wrongwayround_unique %>% 
#'   count(reactivity)
#' td %>% 
#'   count(reactivity)
#' 
#' #' 340/5655
#' #' 164/6861
#' 
#' wrongwayround_final <- wrongwayround_unique %>% 
#'   select(reactivity,id,gene_id_sal,gene_id_p01) 

# wrongwayround_final
# wrongwayround_unique %>% 
#   select(gene_id_p01,product_description_p01,ortholog_group,
#          id,gene_id_sal,product_description) %>% 
#   epihelper::print_inf()

# # _con ortholog ---
# 
# unir_conortho <- unir_sal %>% 
#   # 260
#   # filter(magrittr::is_in(ortholog_group,
#   #                        sal_ortholog_unico %>% pull(ortholog_group))) %>% 
#   # filter(!str_starts(gene_id_sal,"PVAD")) %>% 
#   left_join(
#     pvp01 %>% 
#       select(-product_description) %>%
#       filter(!str_detect(ortholog_group,"N/A")) %>%
#       arrange(gene_id_p01,
#               # product_description,
#               ortholog_group,
#               non_syn_syn_snp_ratio_all_strains) %>% 
#       # paso clave:
#       # aqui decido conservar
#       # la informacion de un solo antigeno por grupo ortologos
#       # de 6861 observaciones
#       # pasamos a 5363 observaciones
#       # distinct(product_description,.keep_all = T)
#       distinct(ortholog_group,.keep_all = T)
#   ) %>% 
#   # naniar::miss_var_summary()
#   # epihelper::print_inf()
#   select(-starts_with("source")) %>% 
#   arrange(gene_id_sal,product_description,
#           #ortholog_group,
#           non_syn_syn_snp_ratio_all_strains,
#           #total_sn_ps_all_strains
#   )
# 
# unir_conortho %>% glimpse()
# unir_conortho %>% naniar::miss_var_summary()
# 
# unir_conortho %>% 
#   filter(is.na(gene_id_p01)) #%>% 
#   # count(ortholog_group,sort = T)
# 
# 
# # _sin ortholog ---
# 
# unir_sinortho <- unir_conortho %>% 
#   filter(is.na(gene_id_p01)) %>% 
#   select(-ortholog_group,
#          -gene_id_p01,
#          -non_syn_syn_snp_ratio_all_strains,
#          -total_sn_ps_all_strains) %>% 
#   # 80
#   # filter(magrittr::is_in(ortholog_group,
#   #                        sal_ortholog_vario %>% pull(ortholog_group))) %>% 
#   # 78
#   # filter(!str_detect(ortholog_group,"N/A"))
#   # filter(!str_starts(gene_id_sal,"PVAD")) %>% 
#   left_join(
#     pvp01 %>% 
#       filter(!str_detect(ortholog_group,"N/A")) %>%
#       # select(-product_description) %>%
#       select(-ortholog_group) %>%
#       arrange(gene_id_p01,
#               product_description,
#               # ortholog_group,
#               non_syn_syn_snp_ratio_all_strains) %>% 
#       # paso clave:
#       # aqui decido conservar
#       # la informacion de un solo antigeno por grupo ortologos
#       # de 6861 observaciones
#       # pasamos a 5363 observaciones
#       distinct(product_description,.keep_all = T)
#       # distinct(ortholog_group,.keep_all = T)
#   ) %>% 
#   # naniar::miss_var_summary()
#   # epihelper::print_inf()
#   select(-starts_with("source")) %>% 
#   arrange(gene_id_sal,product_description,
#           #ortholog_group,
#           non_syn_syn_snp_ratio_all_strains,
#           #total_sn_ps_all_strains
#   )
# 
# # unir_sinortho %>% glimpse()
# unir_sinortho %>% naniar::miss_var_summary()
# 
# # unir_sinortho %>% 
# #   filter(is.na(gene_id_p01)) %>% 
# # # count(ortholog_group,sort = T)
# #   select(1:5) %>% 
# #   epihelper::print_inf()
# 
# 
# # final ---
# 
# unir_conortho %>% naniar::miss_var_summary()
# 
# unir_final <- unir_conortho %>% 
#   select(reactivity,id,gene_id_sal,gene_id_p01) %>% 
#   filter(!is.na(gene_id_p01))

# usando pv p01 -----------------------------------------------------------

# warning
# replace p01 dataset
pvp01 <- read_tsv("data/20220103-pvivax-p01-variability-metadata-v2.tsv") %>% 
  janitor::clean_names() %>% 
  # glimpse()
  select(gene_id,
         total_sn_ps_all_strains:transcript_length) %>% 
  rename(gene_id_p01=gene_id)

pvp01 %>% glimpse()

all_listgen <- pvp01 %>% 
  left_join(p01_sal_union_strategy_p01 %>% 
              select(-total_sn_ps_all_strains,
                     -non_syn_syn_snp_ratio_all_strains)) %>% 
  union_all(p01_sal_union_strategy_sal1) %>%
  # count(reactivity)
  # left_join(unir_final) %>% 
  mutate(reactivity=if_else(is.na(reactivity),
                            "not_reactive",
                            reactivity))

all_listgen %>% glimpse()
all_listgen %>% naniar::miss_var_summary()
all_listgen %>% 
  count(reactivity)

# all_listgen %>% 
#   filter(reactivity=="reactive") %>% 
#   epihelper::print_inf()

# code from 06-...Rmd ------------------------------------------------------

# exp_t <- td %>% 
#   # arrange(desc(AveExpr)) %>% 
#   dplyr::arrange(desc(sample)) %>% 
#   select(-feature_specie,-t,-B,-p.order,-significance,-AveExpr) %>% 
#   # select(ID,#Gene.ID,Description,
#   #        AveExpr,p.order,logFC) %>% 
#   slice(1:16) #%>% 
# # inner_join(all, by="Gene.ID") %>% 
# # inner_join(lg %>% dplyr::rename(seg.num=n), by="Gene.ID")
# 
# # exp_t

# _compare diversity - selection pressure -----------

all_listgen_snp <- 
  all_listgen %>%
  # count(gene_id_p01,sort = T)
  # glimpse()
  mutate(gene_id=case_when(
    is.na(gene_id_p01) ~ gene_id_sal,
    TRUE ~ gene_id_p01 
  )) %>% 
  # count(gene_id,sort = T)
  select(reactivity,
         gene_id,
         # non_syn_syn_snp_ratio_all_strains,
         total_sn_ps_all_strains,
         synonymous_sn_ps_all_strains,
         non_synonymous_sn_ps_all_strains,
         transcript_length) %>%
  mutate(total_sn_ps_all_strains_norm=
           (total_sn_ps_all_strains/transcript_length)*(10^3),
         synonymous_sn_ps_all_strains_norm=
           (synonymous_sn_ps_all_strains/transcript_length)*(10^3),
         non_synonymous_sn_ps_all_strains_norm=
           (non_synonymous_sn_ps_all_strains/transcript_length)*(10^3))

all_listgen_snp %>% glimpse()
all_listgen_snp %>% naniar::miss_var_summary()

# __distribution ---------------------------------------------------------


all_listgen_snp %>%
  pivot_longer(
    cols = -c(reactivity,gene_id)) %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~name,scales = "free")
ggsave("figure/08-fig01-reactive_compare-histogram-all.png",
       height = 6,width = 6,dpi = "retina")

all_listgen_snp %>%
  pivot_longer(
    cols = -c(reactivity,gene_id)) %>%
  ggplot(aes(x = value, fill = reactivity)) +
  geom_histogram(position = position_identity()) +
  facet_wrap(~name,scales = "free")
ggsave("figure/08-fig02-reactive_compare-histogram-strata.png",
       height = 6,width = 6,dpi = "retina")

# __table -----------------------------------------------------------------


all_listgen_snp_table <- all_listgen_snp %>% 
  compareGroups::compareGroups(
    formula = 
      reactivity ~ 
      total_sn_ps_all_strains +
      synonymous_sn_ps_all_strains +
      non_synonymous_sn_ps_all_strains +
      total_sn_ps_all_strains_norm +
      synonymous_sn_ps_all_strains_norm +
      non_synonymous_sn_ps_all_strains_norm,
    data = .,
    method = list(total_sn_ps_all_strains = 2,
                  synonymous_sn_ps_all_strains = 2,
                  non_synonymous_sn_ps_all_strains = 2,
                  total_sn_ps_all_strains_norm = 2,
                  synonymous_sn_ps_all_strains_norm = 2,
                  non_synonymous_sn_ps_all_strains_norm = 2)) %>% 
  compareGroups::createTable()

all_listgen_snp_table

all_listgen_snp_table %>% 
  compareGroups::export2xls("table/08-table01-reactive_compare-snp_variability.xls")

# __raw -------------------------------------------------------------------


all_listgen_snp %>% 
  select(gene_id,reactivity,non_synonymous_sn_ps_all_strains) %>% 
  mutate(reactivity = as.factor(reactivity)) %>% 
  wilcox.test(non_synonymous_sn_ps_all_strains ~ reactivity, 
              data = .) %>% 
  broom::tidy()

all_listgen_snp %>% 
  select(gene_id,reactivity,non_synonymous_sn_ps_all_strains_norm) %>% 
  mutate(reactivity = as.factor(reactivity)) %>% 
  wilcox.test(non_synonymous_sn_ps_all_strains_norm ~ reactivity, 
              data = .) %>% 
  broom::tidy()

# [pendiente] ---------------------------------------------------------------
#' 
#' hecho: usar solo pvp01 como referencia
#' posibilidad: conseguir lista completa de genes y snp values

