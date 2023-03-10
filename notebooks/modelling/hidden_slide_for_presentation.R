```{r}
library(readr)
library(dplyr)
library(knitr)
data_strain <- read_csv("data/features_strain.csv", na = c("NA", "-Inf"))

# TWO SCHEMES IN TWO FILES
hr_testing12 <- read_csv("data/heteroresistance_testing.csv", col_select = c(strain, Gr12)) %>% 
  filter(!is.na(strain)) %>% 
  rename("resistance" = Gr12 )

data_strain <- data_strain %>% 
  left_join(hr_testing12, by = "strain")

data_strain <- data_strain %>% 
  mutate(n.beta.lac.3 = factor(ifelse(n.beta.lac > 3, "yes", "no"))) %>% 
  mutate(n.beta.lac.4 = factor(ifelse(n.beta.lac > 4, "yes", "no"))) %>% 
  relocate(n.beta.lac.3, n.beta.lac.4, .before = "n.plasmids") %>% 
  filter(resistance != "R") %>% 
  mutate(resistance = factor(resistance, levels = c("HR", "nonHR")),
         chrom.status = factor(chrom.status))

# hr_testing12 %>% 
#   group_by(resistance) %>% 
#   count() %>% 
#   kable()
```


```{r}
library(ggplot2)
data_strain %>% 
  group_by(resistance, chrom.status) %>% 
  count() %>% 
  ggplot(aes(chrom.status, n)) +
  geom_col(position = "dodge") + 
  facet_grid(cols = vars(resistance))
```