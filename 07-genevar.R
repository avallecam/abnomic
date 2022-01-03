library(tidyverse)

# importar ----------------------------------------------------------------

pvp01 <- read_tsv("data/20211223-pvivax-p01-ortholog_group-more.tsv") %>% 
  janitor::clean_names() %>% 
  # glimpse()
  mutate(number_tm_domains=as.numeric(number_tm_domains)) %>% 
  select(-organism,-x13) %>% 
  select(-signal_p_scores,-number_tm_domains) %>% 
  select(-ec_numbers_from_ortho_mcl) %>% 
  # select(-ortholog_group,-ortholog_count,-paralog_count) %>% 
  rename(gene_id_p01=gene_id,
         source_id_p01=source_id) 
pvsal <- read_tsv("data/20211223-pvivax-sal1-ortholog_group-more.tsv") %>% 
  janitor::clean_names() %>% 
  mutate(number_tm_domains=as.numeric(number_tm_domains)) %>% 
  select(-organism,-x13) %>% 
  select(-signal_p_scores,-number_tm_domains) %>% 
  select(-ec_numbers_from_ortho_mcl) %>% 
  select(-total_sn_ps_all_strains,
         -non_syn_syn_snp_ratio_all_strains) %>% 
  # select(-ortholog_group,-ortholog_count,-paralog_count) %>% 
  rename(gene_id_sal=gene_id,
         source_id_sal=source_id)

pvp01 %>% glimpse()
pvsal %>% glimpse()

# revisar replicas --------------------------------------------------------

pvsal %>% 
  count(#product_description,
        ortholog_group,
        # ortholog_count,paralog_count,
        # ec_numbers_from_ortho_mcl,
        # signal_p_scores,number_tm_domains,
        sort = T) %>% 
  filter(n>1)

pvp01 %>% 
  count(#product_description,
        ortholog_group,
        # ortholog_count,paralog_count,
        # ec_numbers_from_ortho_mcl,
        # signal_p_scores,number_tm_domains,
        sort = T) %>% 
  filter(n>1)

# como unir ---------------------------------------------------------------

# pvsal %>% count(ortholog_group,ortholog_count,paralog_count,
#                 sort = T)

# pvp01 %>% 
#   filter(ortholog_group=="OG6_145873") %>% 
#   glimpse()

pvp01 %>% 
  count(ortholog_group, ortholog_count, paralog_count,sort = T) %>% 
  filter(n>1) #%>% 
  # epihelper::print_inf()

pvp01 %>% 
  filter(str_detect(ortholog_group,"N/A")) %>% 
  count(product_description, sort = T)

# pvp01 %>% 
#   filter(product_description=="tryptophan-rich protein")
# pvsal %>% 
#   filter(product_description=="tryptophan-rich antigen (Pv-fam-a)")

# variacion snp por ortolog group -----------------------------------------

pvp01 %>% 
  select(ortholog_group,non_syn_syn_snp_ratio_all_strains,
         total_sn_ps_all_strains) %>% 
  group_by(ortholog_group) %>% 
  filter(n()>1) %>% 
  ungroup() %>% 
  ggplot(aes(x = non_syn_syn_snp_ratio_all_strains, y = ortholog_group)) +
  geom_point()

pvp01 %>% 
  select(ortholog_group,non_syn_syn_snp_ratio_all_strains,
         total_sn_ps_all_strains) %>% 
  group_by(ortholog_group) %>% 
  filter(n()>1) %>% 
  ungroup() %>% 
  ggplot(aes(x = total_sn_ps_all_strains, y = ortholog_group)) +
  geom_point()

pvp01 %>% 
  select(ortholog_group,non_syn_syn_snp_ratio_all_strains,
         total_sn_ps_all_strains) %>% 
  group_by(ortholog_group) %>% 
  filter(n()>1) %>% 
  skimr::skim() %>% 
  as_tibble() %>% 
  arrange(desc(numeric.sd)) %>% 
  ggplot(aes(x = numeric.sd)) +
  geom_histogram() +
  facet_grid(~skim_variable,scales = "free")

# unir --------------------------------------------------------------------

# _ con ortholog ----------------------------------------------------------

unir_pre <- pvsal %>% 
  # filter(!str_starts(gene_id_sal,"PVAD")) %>% 
  left_join(
    pvp01 %>% 
      select(-product_description) %>% 
      filter(!str_detect(ortholog_group,"N/A")) %>%
      arrange(gene_id_p01,
              # product_description,
              ortholog_group,
              non_syn_syn_snp_ratio_all_strains) %>% 
      # paso clave:
      # aqui decido conservar
      # la informacion de un solo antigeno por grupo ortologos
      # de 6861 observaciones
      # pasamos a 5363 observaciones
      # distinct(product_description,.keep_all = T)
      distinct(ortholog_group,.keep_all = T)
            ) %>% 
  select(-starts_with("source")) %>% 
  arrange(gene_id_sal,product_description,
          #ortholog_group,
          non_syn_syn_snp_ratio_all_strains,
          #total_sn_ps_all_strains
          )

unir_pre %>% glimpse()
unir_pre %>% naniar::vis_miss()
unir_pre %>% naniar::miss_var_summary()

unir_pre %>% 
  filter(str_detect(ortholog_group,"N/A")) %>%
  count(product_description, sort = T)

unir_pre %>% 
  filter(is.na(gene_id_p01)) %>% 
  count(ortholog_group,sort = T)

# _ sin ortholog ----------------------------------------------------------

unir_pre_con_ortholog <- unir_pre %>% 
  filter(!str_detect(ortholog_group,"N/A"))

unir_pre_sin_ortholog <- unir_pre %>% 
  filter(str_detect(ortholog_group,"N/A")) %>% 
  select(#-ortholog_group,-ortholog_count,-paralog_count,
         -gene_id_p01,
         -non_syn_syn_snp_ratio_all_strains,-total_sn_ps_all_strains) %>% 
  left_join(
    pvp01 %>% 
      # select(-product_description) %>% 
      filter(str_detect(ortholog_group,"N/A")) %>%
      arrange(gene_id_p01,
              product_description,
              #ortholog_group,
              non_syn_syn_snp_ratio_all_strains) %>% 
      # paso clave:
      # aqui decido conservar
      # la informacion de un solo antigeno por grupo ortologos
      # de 6861 observaciones
      # pasamos a 5363 observaciones
      distinct(product_description,.keep_all = T) %>% 
      # distinct(ortholog_group,.keep_all = T)
      select(-ortholog_group,-ortholog_count,-paralog_count)
  ) %>% 
  select(-starts_with("source")) %>% 
  arrange(gene_id_sal,product_description,
          #ortholog_group,
          non_syn_syn_snp_ratio_all_strains,
          #total_sn_ps_all_strains
  )

unir_pre # 5665
unir_pre_con_ortholog # 5530
unir_pre_sin_ortholog # 135

unir_pos <- 
  unir_pre_con_ortholog %>% 
  union_all(unir_pre_sin_ortholog)

unir_pos %>% glimpse()
unir_pos %>% naniar::vis_miss()
unir_pos %>% naniar::miss_var_summary()

unir_pos %>%
  filter(is.na(gene_id_p01)) %>%
  arrange(product_description)
#   epihelper::print_inf()

# unir: evaluar 

# unir_pre %>% 
#   group_by(gene_id_sal,product_description,ortholog_group) %>%
#   filter(n()>1) %>%
#   # count(gene_id_sal,sort = T) %>% 
#   # epihelper::print_inf()
#   ungroup() %>% 
#   select(-ortholog_group,-ortholog_count,
#          -paralog_count,-ec_numbers_from_ortho_mcl,
#          -signal_p_scores,-number_tm_domains) %>%
#   # filter(!str_detect(product_description,"tRNA")) %>% 
#   epihelper::print_inf() #%>%
#   # filter(gene_id_sal=="PVX_000000") %>% 
#   # filter(is.na(gene_id_p01)) %>% 
#   # filter(!str_detect(product_description,"unknown")) %>% 
#   # filter(!str_detect(product_description,"VIR")) %>% 
#   # naniar::miss_var_summary()

# unir: final 


# comparar null vs reactives ----------------------------------------------

unir_pos %>% 
  select(-ortholog_count,-paralog_count) %>% 
  write_rds("data/20211223-pvivax-sal1-p01-snp_variant_solution.rds")

#' pendiente:
#' - usar esta lista como distribución nula para comparación de snp por gene en 06-uscntrl.Rmd!
#' - solución preliminar a respuesta de plasmoDB

