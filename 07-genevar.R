library(tidyverse)

# importar ----------------------------------------------------------------

pvp01 <- read_tsv("data/20211223-pvivax-p01-ortholog_group-more.tsv") %>% 
  janitor::clean_names() %>% 
  # glimpse()
  select(-organism,-x13,
         # -signal_p_scores,-number_tm_domains,
         # -ec_numbers_from_ortho_mcl
         ) %>% 
  rename(gene_id_p01=gene_id,
         source_id_p01=source_id) %>% 
  mutate(number_tm_domains=as.numeric(number_tm_domains))
pvsal <- read_tsv("data/20211223-pvivax-sal1-ortholog_group-more.tsv") %>% 
  janitor::clean_names() %>% 
  select(-organism,-x13,
         # -signal_p_scores,-number_tm_domains,
         # -ec_numbers_from_ortho_mcl,
         -total_sn_ps_all_strains,
         -non_syn_syn_snp_ratio_all_strains) %>% 
  rename(gene_id_sal=gene_id,
         source_id_sal=source_id) %>% 
  mutate(number_tm_domains=as.numeric(number_tm_domains))

pvp01 %>% glimpse()
pvsal %>% glimpse()

# revisar replicas --------------------------------------------------------

pvsal %>% 
  count(product_description,ortholog_group,
        ortholog_count,paralog_count,
        ec_numbers_from_ortho_mcl,
        signal_p_scores,number_tm_domains,
        sort = T) %>% 
  filter(n>1)

pvp01 %>% 
  count(product_description,ortholog_group,
        ortholog_count,paralog_count,
        ec_numbers_from_ortho_mcl,
        signal_p_scores,number_tm_domains,
        sort = T) %>% 
  filter(n>1)

# unir --------------------------------------------------------------------

unir_pre <- pvsal %>% 
  # filter(!str_starts(gene_id_sal,"PVAD")) %>% 
  left_join(pvp01) %>% 
  select(-starts_with("source")) %>% 
  arrange(gene_id_sal,product_description,ortholog_group,
          non_syn_syn_snp_ratio_all_strains,
          #total_sn_ps_all_strains
          )

# unir: evaluar -----------------------------------------------------------

unir_pre %>% 
  group_by(gene_id_sal,product_description,ortholog_group) %>%
  filter(n()>1) %>%
  # count(gene_id_sal,sort = T) %>% 
  # epihelper::print_inf()
  ungroup() %>% 
  select(-ortholog_group,-ortholog_count,
         -paralog_count,-ec_numbers_from_ortho_mcl,
         -signal_p_scores,-number_tm_domains) %>%
  # filter(!str_detect(product_description,"tRNA")) %>% 
  epihelper::print_inf() #%>%
  # filter(gene_id_sal=="PVX_000000") %>% 
  # filter(is.na(gene_id_p01)) %>% 
  # filter(!str_detect(product_description,"unknown")) %>% 
  # filter(!str_detect(product_description,"VIR")) %>% 
  # naniar::miss_var_summary()

# unir: final -------------------------------------------------------------

#' pendiente:
#' - conservar un valor de snp por  gene, priorizar a lo más conservador, el menor ratio.
#' - usar esta lista como distribución nula para comparación de snp por gene en 06-uscntrl.Rmd!
#' - solución preliminar a respuesta de plasmoDB

pvsal

unir_pre %>% 
  distinct(gene_id_sal,product_description,ortholog_group,
           non_syn_syn_snp_ratio_all_strains)

