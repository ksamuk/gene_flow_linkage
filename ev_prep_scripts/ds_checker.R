# omg testing

library(dplyr)
library(ggplot2)

ds <- read.table(file = "evs/ds.txt", stringsAsFactors = FALSE, header = TRUE)
ds_old <- read.table(file = "evs/ds_old.txt", stringsAsFactors = FALSE, header = TRUE)

ds_old <- ds_old %>%
  group_by(lg, pos1) %>%
  summarise(pos2 = mean(pos2), ds = mean(ds)) %>%
  ungroup

ds.lg1 <- ds %>%
  filter(lg == 1) %>%
  filter(ds > 10)
ds_old.lg1 <- ds_old %>%
  filter(lg == 1)%>%
  filter(ds < 10)

ds_old.lg1 <- ds_old.lg1[ds_old.lg1$pos1 %in% ds.lg1$pos1,]

ds.lg1 <- ds.lg1[ds.lg1$pos1 %in% ds_old.lg1$pos1,]

ds_old.lg1$ds %>% length
ds.lg1$ds %>% length


# plot

ds %>%
  ggplot(aes(x = pos1, y = ds)) +
  geom_smooth() +
  geom_smooth(data = ds_old, aes(x = pos1, y = ds)) +
  facet_wrap(~lg)



