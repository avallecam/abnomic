library(tidyverse)

# importar ag reactivos ---------------------------------------------------

td <- write_rds("data/06-reactive_compare-all_reactives.rds")
td %>% 
  select(id=ID,reactivity) %>% 
  mutate(gene_id_sal = str_replace(id,"(PVX_.{6}).+","\\1"))


# importar referencia -----------------------------------------------------

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

