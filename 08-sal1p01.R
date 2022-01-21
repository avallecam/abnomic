library(tidyverse)

# importar ag reactivos ---------------------------------------------------

td <- read_rds("data/06-reactive_compare-all_reactives.rds")
td %>% count(reactivity)

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

# _test -------------------------------------------------------------------

wrongwayround <- pvp01 %>% 
  # left join only by ortholog
  rename(product_description_p01 = product_description) %>%
  # rename(ortholog_group_p01 = ortholog_group) %>%
  left_join(unir_sal) %>% 
  mutate(reactivity=if_else(is.na(reactivity),
                            "not_reactive",
                            reactivity)) 

wrongwayround %>% 
  count(reactivity)
  # naniar::vis_miss()
  # naniar::miss_var_summary()

wrongwayround_unique <- wrongwayround %>% 
  arrange(gene_id_p01,desc(non_syn_syn_snp_ratio_all_strains)) %>% 
  distinct(gene_id_p01,.keep_all = T)
 
wrongwayround_unique %>% 
  count(reactivity)

wrongwayround_final <- wrongwayround_unique %>% 
  select(reactivity,id,gene_id_sal,gene_id_p01) 

# # _con ortholog ---------------------------------------------------------
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
# # _sin ortholog ---------------------------------------------------------
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
# # final -------------------------------------------------------------------
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

all_listgen <- pvp01 %>% 
  left_join(wrongwayround_final)
  # left_join(unir_final) %>% 
  # mutate(reactivity=if_else(is.na(reactivity),
  #                           "not_reactive",
  #                           reactivity))

all_listgen %>% glimpse()
all_listgen %>% 
  count(reactivity)

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
  # glimpse()
  select(reactivity,
         gene_id=gene_id_p01,
         # non_syn_syn_snp_ratio_all_strains,
         total_sn_ps_all_strains,
         synonymous_sn_ps_all_strains,
         non_synonymous_sn_ps_all_strains,
         transcript_length) %>%
  mutate(total_sn_ps_all_strains_norm=(total_sn_ps_all_strains/transcript_length)*(10^3),
         synonymous_sn_ps_all_strains_norm=(synonymous_sn_ps_all_strains/transcript_length)*(10^3),
         non_synonymous_sn_ps_all_strains_norm=(non_synonymous_sn_ps_all_strains/transcript_length)*(10^3))

all_listgen_snp %>% glimpse()

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

# [pendiente] ---------------------------------------------------------------
#' 
#' hecho: usar solo pvp01 como referencia
#' posibilidad: reemplazar sal1-2017 values de reactive antigens [no]

