# prepare glazer recombination estimates for reanalysis

library("dplyr")

glz_recomb <- read.table("evs/recomb_rate_glazer.txt", h = T)

glz_recomb <- glz_recomb %>%
	rename(lg = chr , glazer_rate = recomb_rate) %>%
	select(lg, pos1, pos2, glazer_rate)

write.table(glz_recomb, "evs/glazer_rate.txt", quote = FALSE, row.names = FALSE)
