library(tidyverse)


# awesome gert output -----------------------------------------------------

gert::git_log() %>% 
  select(time,files,message) %>% 
  avallecam::print_inf()
  # view()
  # identity()

# year frequency ----------------------------------------------------------


gert::git_log() %>% 
  mutate(time=lubridate::as_date(time)) %>% 
  ggplot(aes(x = time)) +
  geom_histogram()

# weekly frequency by month -----------------------------------------------


gert::git_log() %>% 
  mutate(time=lubridate::as_date(time)) %>% 
  mutate(SEMANA_EPI = aweek::date2week(time,
                                       floor_day = TRUE,
                                       week_start = "Sunday")) %>% 
  count(SEMANA_EPI) %>% 
  mutate(year=lubridate::year(SEMANA_EPI),
         week=lubridate::epiweek(SEMANA_EPI)) %>% 
  cdcper::cdc_yearweek_to_date(year_integer = year,week_integer = week) %>% 
  # mutate(month=lubridate::month(epi_date,label = TRUE)) %>% 
  # select(month,n) %>% 
  # # #expand dataframe for weeks witout reports
  # # complete(month = full_seq(month,1),
  # #          fill = list(var_event_count=0))
  # # pull(month)
  # # glimpse()
  # identity()
  # # mutate(SEMANA_EPI = as.factor(SEMANA_EPI)) %>% 
  # ggplot(aes(x = month,y = n)) +
  ggplot(aes(x = epi_date,y = n)) +
  geom_col() +
  # NULL
  scale_x_date(date_breaks = "1 month",date_labels = "%Y-%m") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # facet_grid(~year)
  NULL

# hour frequency by year --------------------------------------------------


gert::git_log() %>% 
  mutate(hour=lubridate::hour(time),
         year=lubridate::year(time)) %>% 
  count(year,hour) %>% 
  mutate(n=as.integer(n)) %>% 
  ggplot(aes(x = hour,y = year,fill = n)) +
  # geom_histogram() +
  geom_tile() +
  scale_x_continuous(breaks = 0:23,labels = 0:23) +
  # NULL
  colorspace::scale_fill_continuous_sequential() +
  # scale_fill_continuous()
  # colorspace::scale_fill_binned_sequential()+
  theme_bw()
  # guides(fill = guide_colorsteps(barheight = unit(11, "cm")))

# volcano_long <- data.frame(
#   x = as.vector(col(volcano)),
#   y = as.vector(row(volcano)),
#   z = as.vector(volcano)
# )
# ggplot(volcano_long, aes(x, y, z = z)) + 
#   geom_contour_filled(aes(fill = stat(level))) + 
#   guides(fill = guide_colorsteps(barheight = unit(10, "cm")))

