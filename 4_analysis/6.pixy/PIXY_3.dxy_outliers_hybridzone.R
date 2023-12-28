library(tidyverse)

#list of scaffolds to remove
badscafs <- c("scaffold_2531","scaffold_2446","scaffold_2151","scaffold_1085")


#read in data
dxy_10kb <- read_delim("fullgenome_fullspecies_10kb_dxy.txt")
dxy_50kb <- read_delim("fullgenome_fullspecies_50kb_dxy.txt")

#filter data to make plots
filtered_10kb <- dxy_10kb %>%
  filter(!chromosome %in% badscafs) %>%
  na.omit() %>%
  mutate(pop1 = as.character(str_sub_all(pop1, 3, 5))) %>%
  mutate(pop2 = as.character(str_sub_all(pop2, 3, 5))) %>%
  mutate(midpoint = ceiling((window_pos_1-1 + window_pos_2)/2)/1000000) %>%
  filter(pop1 %in% c("dav", "new", "rup") & pop2 %in% c("dav", "new", "rup")) %>%
  mutate(pop_comparison = paste(pop1, pop2, sep = "_")) %>%
  filter(pop_comparison != "new_rup") %>%
  mutate(chromosome = gsub("fold", "", chromosome))

filtered_50kb <- dxy_50kb %>%
  filter(!chromosome %in% badscafs) %>%
  na.omit() %>%
  mutate(pop1 = as.character(str_sub_all(pop1, 3, 5))) %>%
  mutate(pop2 = as.character(str_sub_all(pop2, 3, 5))) %>%
  mutate(midpoint = ceiling((window_pos_1-1 + window_pos_2)/2)/1000000) %>%
  filter(pop1 %in% c("dav", "new", "rup") & pop2 %in% c("dav", "new", "rup")) %>%
  mutate(pop_comparison = paste(pop1, pop2, sep = "_")) %>%
  filter(pop_comparison != "new_rup") %>%
  mutate(chromosome = gsub("fold", "", chromosome))


#identify outlier dxy values based on z score
outliers_10kb <- filtered_10kb %>%
  mutate(zscore = ((avg_dxy - mean(avg_dxy))/sd(avg_dxy))) %>%
  filter(abs(zscore) > 3)

outliers_50kb <- filtered_50kb %>%
  mutate(zscore = ((avg_dxy - mean(avg_dxy))/sd(avg_dxy))) %>%
  filter(abs(zscore) > 3)


#plot these results
ggplot(filtered_10kb, aes(x = midpoint, y = avg_dxy)) +
  geom_point(size = 0.2) +
  geom_point(data = outliers_10kb, colour = "red", size = 0.4) +
  facet_grid(pop_comparison~chromosome, scales = 'free_x', space = 'free_x')

ggplot(filtered_50kb, aes(x = midpoint, y = avg_dxy)) +
  geom_point(size = 0.2) +
  geom_point(data = outliers_50kb, colour = "red", size = 0.4) +
  facet_grid(pop_comparison~chromosome, scales = 'free_x', space = 'free_x')



#check if any outlier windows are shared between the two comparisons
shared_outliers_10kb <- outliers_10kb %>%
  group_by(chromosome, midpoint) %>%
  filter(n() > 1) %>%
  ungroup()

shared_outliers_50kb <- outliers_50kb %>%
  group_by(chromosome, midpoint) %>%
  filter(n() > 1) %>%
  ungroup()


#plot shared outliers
a <- ggplot(filtered_10kb, aes(x = midpoint, y = avg_dxy)) +
  geom_point(size = 0.2) +
  geom_point(data = shared_outliers_10kb, colour = "red", size = 0.4) +
  facet_grid(pop_comparison~chromosome, scales = 'free_x', space = 'free_x') +
  theme_bw() +
  labs(x = "Position on chromosome (Mbp)",
       y = expression(paste("Average ", italic(d[xy])))) +
  theme(panel.spacing.x = unit(0.01, "in")) +
  ggtitle("Shared dxy outliers 10kb windows")
write.csv(shared_outliers_10kb, "shared_dxy_outliers_10kb.csv", row.names = F)
png(filename = "shared_dxy_outliers_10kb.png",
    height = 2.5, width = 9, units = "in", res = 400)
a
dev.off()

b <- ggplot(filtered_50kb, aes(x = midpoint, y = avg_dxy)) +
  geom_point(size = 0.2) +
  geom_point(data = shared_outliers_50kb, colour = "red", size = 0.4) +
  facet_grid(pop_comparison~chromosome, scales = 'free_x', space = 'free_x') +
  theme_bw() +
  labs(x = "Position on chromosome (Mbp)",
       y = expression(paste("Average ", italic(d[xy])))) +
  theme(panel.spacing.x = unit(0.01, "in")) +
  ggtitle("Shared dxy outliers 50kb windows")
write.csv(shared_outliers_50kb, "shared_dxy_outliers_50kb.csv", row.names = F)
png(filename = "shared_dxy_outliers_50kb.png",
    height = 2.5, width = 9, units = "in", res = 400)
b
dev.off()


#write bed file for shared outliers to ID genes in these regions later
write_delim(shared_outliers_10kb %>%
              filter(pop_comparison == "dav_new") %>%
              select(c(chromosome, window_pos_1, window_pos_2, zscore)) %>%
              mutate(chromosome = gsub("scaf", "scaffold", chromosome)),
            "shared_dxy_outliers_10kb.bed",
            delim = '\t', col_names = F)


write_delim(shared_outliers_50kb %>%
              filter(pop_comparison == "dav_new") %>%
              select(c(chromosome, window_pos_1, window_pos_2, zscore)) %>%
              mutate(chromosome = gsub("scaf", "scaffold", chromosome)),
            "shared_dxy_outliers_50kb.bed",
            delim = '\t', col_names = F)


