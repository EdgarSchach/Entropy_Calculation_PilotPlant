# Loading packages, preparing data and 
library(tidyverse)
library(RColorBrewer)
library(compositions)
load("./MLA_data/MLA_Particles.RData")
source(file="./funtions.paper.2.R")
source("./mlaXMLSession.R")
source("./Entropy_Functions-toSubmit.R")

# Colors and Leves for Graphs####

cols_component.dist <- brewer.pal(n = 5, name = "Blues")[2:5]
names(cols_component.dist) <- c("c0j0","cij0","c0jk","cijk")

cols_particle.dist <- brewer.pal(n = 4, name = "Greens")
names(cols_particle.dist) <- c("p00k","p0jk","pi0k","pijk")

cols_partition.vec <- brewer.pal(n = 4, name = "Reds")
names(cols_partition.vec) <- c("ti00","ti0k","tij0","tijk")

entropy_colors <- c(cols_component.dist,
                    cols_particle.dist,
                    cols_partition.vec)

type.levels <- c("ti00","ti0k","tij0","tijk", # partitions
                 "p00k","p0jk","pi0k","pijk" ,# distributions
                 "c0j0","cij0","c0jk","cijk") # compositions

sample_lvl <- c("Feed","Recal Feed","Magnetic","Non-magnetic","Product","Coarse","Fine")
machine_lvl <- c("Sieve","Magnetic separator","Mill")

#Grouping MLA Data ####

mlaSetXMLGroupsAndPalette("./Mineral_grouping.xml") # pay attention to the right path
aux = options()$mlaGrouping
rownames(aux$groupingMatrix)[25]  <-"Tourmaline (Schorl-Dravite_Na-Mg)"
options(mlaGrouping = aux)

gm <- aux$groupingMatrix
gm <- as.data.frame(t(gm))
gm[,"Mn-Silicate"] <- c(1,rep(0,times = 11))
gm <- as.data.frame(t(gm))

mineral_masses <-   lapply(SamplesMLA, function(x){
  res <- as.data.frame(as.matrix(x$minMassComp) %*% as.matrix(gm[colnames(x$minMassComp),]))
  return(res)
})

mineral_masses_clean <- lapply(mineral_masses,function(x)select(x,-contains("Tourmaline")))
part_dat <- mineral_masses_clean

#Calculating grouped particle data####



of.interest <- c("Cassiterite","Sphalerite","Chalcopyrite","Arsenopyrite","Iron Oxides")

group_minerals <- function(x){
  rest <- x %>% 
    select(-any_of(of.interest)) %>% 
    rowSums()
  x <- x %>%
    select(any_of(of.interest)) %>% 
    mutate(Rest = rest)
  return(x)
}

part_dat_grouped <- map(part_dat,group_minerals)

#loading machine data####

load("./machine_data.RData")

machines$sieve$input <- c("S01","S06")
machines$sieve$weights <- c(S01 = 810.05,S06 = 477.63, S02 = 776.89, S03 = 510.79)
machines$magnetic_separator$weights = c(S02 = 776.89, machines$magnetic_separator$weights)

machines$mill1$weights <- c(S04 = 316.98, S06 = 316.98)
machines$mill2$weights <- c(S16 = 134, S18 = 134)

machine <- machines$magnetic_separator
n = 10
th.mass = 1
machines$sieve$weights <- c(S01 = 810.05,S06 = 477.63, S02 = 776.89, S03 = 510.79)
machines$magnetic_separator$weights = c(S02 = 776.89, machines$magnetic_separator$weights)

machines$hydrocyclone <- list(input = c("S12"),
                              output = c("S13","S14"),
                              weights = c("S12" = 87,
                                          "S13" = 70.8,
                                          "S14" = 16.2)
)

####calculations for the sieve####

#1) Calculation of batch entropies for each sample - mass split already incorporated

res.sieve_batch <- batch_entropies(part_dat = part_dat,machine = machines$sieve, n = 10, th.mass = 1)


res.sieve.agg <- map(res.sieve_batch,bind_rows,.id = "Nr.") %>% 
  bind_rows(.id = "sample") %>% 
  group_by(sample,aggregation,entropy_type) %>% 
  summarise(mean_entropy = mean(entropy), ent_err = sd(entropy))
res.sieve.agg.batch <- group_by(res.sieve.agg,sample,aggregation) %>% mutate(new_ent_mean= cumsum(mean_entropy))
save(res.sieve.agg.batch,file = "sieve_agg_batch.RData")
res.sieve.agg <- filter(res.sieve.agg, !(sample %in% c("feed","S06","S01")))

res.sieve.agg$sample <- str_replace_all(res.sieve.agg$sample,c("S02" = "Fine", "S03" = "Coarse","recal_feed" = "Feed"))
save(res.sieve.agg,file = "res.sieve.agg.batches.RData")
type.levels <- c("ti00","ti0k","tij0","tijk", # partitions
                 "p00k","p0jk","pi0k","pijk" ,# distributions
                 "c0j0","cij0","c0jk","cijk") # compositions

g <- ggplot(res.sieve.agg,aes(x = aggregation, y = mean_entropy)) +
  geom_bar(stat = "identity",aes(fill = factor(entropy_type,levels = type.levels))) +
  geom_errorbar(aes(x = aggregation, ymin = new_ent_mean-ent_err, ymax = new_ent_mean+ent_err), width = 0.3) +
  facet_grid(.~factor(sample,levels = sample_lvl), drop = T, scales = "free", space = "free") +
  scale_fill_manual(values = entropy_colors) +
  theme_bw() +
  labs(x = "aggregation", y = "entropy value", fill = "contribution")
plot(g)

ggsave(filename = "sieve.entropies.batches.png",device = "png",width = 15, height = 10, units = "cm")

#2) calculating entropies stage wise:

res.sieve_stage <- pp_machine_entropy(machine = machines$sieve,part_dat = part_dat,n = 10,th.mass = 1)

res_entropies.sieve <- res.sieve_stage$entropies
res.sieve_stage.agg <- bind_rows(res_entropies.sieve,.id = "Nr.") %>% 
  group_by(stage,aggregation,entropy_type) %>% 
  summarise(ent_mean = mean(entropy),ent_err = sd(entropy))
res.sieve_stage.agg <- group_by(res.sieve_stage.agg,aggregation) %>% mutate(new_ent_mean= cumsum(ent_mean))

type.levels <- c("ti00","ti0k","tij0","tijk", # partitions
                 "p00k","p0jk","pi0k","pijk" ,# distributions
                 "c0j0","cij0","c0jk","cijk") # compositions

g <- ggplot(res_agg,aes(x = aggregation, y = ent_mean)) +
  geom_bar(stat = "identity",aes(fill = factor(entropy_type,levels = type.levels)))+
  geom_errorbar(aes(x = aggregation, ymin = new_ent_mean-ent_err, ymax = new_ent_mean+ent_err), width = 0.3) +
  facet_grid(.~stage, drop = T, scales = "free", space = "free") +
  scale_fill_manual(values = entropy_colors) +
  theme_bw() +
  labs(x = "aggregation", y = "entropy value", fill = "contribution")
plot(g)

res.sieve_stage.grouped <- pp_machine_entropy(machine = machines$sieve,part_dat = part_dat_grouped,n = 10,th.mass = 1)

res_entropies.sieve <- res.sieve_stage.grouped$entropies
res.sieve_stage.grouped.agg <- bind_rows(res_entropies.sieve,.id = "Nr.") %>% 
  group_by(stage,aggregation,entropy_type) %>% 
  summarise(ent_mean = mean(entropy),ent_err = sd(entropy))
res.sieve_stage.grouped.agg <- group_by(res.sieve_stage.grouped.agg,aggregation) %>% mutate(new_ent_mean= cumsum(ent_mean))

g <- ggplot(res.sieve_stage.grouped.agg,aes(x = aggregation, y = ent_mean)) +
  geom_bar(stat = "identity",aes(fill = factor(entropy_type,levels = type.levels)))+
  geom_errorbar(aes(x = aggregation, ymin = new_ent_mean-ent_err, ymax = new_ent_mean+ent_err), width = 0.3) +
  facet_grid(.~stage, drop = T, scales = "free", space = "free") +
  scale_fill_manual(values = entropy_colors) +
  theme_bw() +
  labs(x = "aggregation", y = "entropy value", fill = "contribution")
plot(g)
ggsave(filename = "sieve.entropies.grouped.paper2.png",device = "png",width = 20, height = 10, units = "cm")


####calculations for the magnetic separator####

#1) Calculation of batch entropies for each sample - mass split already incorporated

res.magsep_batch <- batch_entropies(part_dat = part_dat,machine = machines$magnetic_separator, n = 10, th.mass = 1)

res.magsep.agg <- map(res.magsep_batch,bind_rows,.id = "Nr.") %>% 
  bind_rows(.id = "sample") %>% 
  group_by(sample,aggregation,entropy_type) %>% 
  summarise(mean_entropy = mean(entropy), ent_err = sd(entropy))
res.magsep.agg <- group_by(res.magsep.agg,sample,aggregation) %>% mutate(new_ent_mean= cumsum(mean_entropy))


res.magsep.agg <- filter(res.magsep.agg, !(sample %in% c("feed","S02")))

res.magsep.agg$sample <- str_replace_all(res.magsep.agg$sample,c("S07" = "Non-magnetic", "S08" = "Magnetic","recal_feed" = "Feed"))

type.levels <- c("ti00","ti0k","tij0","tijk", # partitions
                 "p00k","p0jk","pi0k","pijk" ,# distributions
                 "c0j0","cij0","c0jk","cijk") # compositions

h <- ggplot(res.magsep.agg,aes(x = aggregation, y = mean_entropy)) +
  geom_bar(stat = "identity",aes(fill = factor(entropy_type,levels = type.levels))) +
  geom_errorbar(aes(x = aggregation, ymin = new_ent_mean-ent_err, ymax = new_ent_mean+ent_err), width = 0.3) +
  facet_grid(.~factor(sample,levels = sample_lvl), drop = T, scales = "free", space = "free") +
  scale_fill_manual(values = entropy_colors) +
  theme_bw() +
  labs(x = "aggregation", y = "entropy value", fill = "contribution")
plot(h)
save(res.magsep.agg,file = "res.magsep.batches.agg.RData")
ggsave(filename = "magsep.entropies.batches.png",device = "png",width = 15, height = 10, units = "cm")

#2) calculating entropies stage wise:

#for "normaly" grouped minerals:

res.magsep_stage <- pp_machine_entropy(machine = machines$magnetic_separator,part_dat = part_dat,n = 10,th.mass = 1)

res_entropies.magsep <- res.magsep_stage$entropies
res.magsep_stage.agg <- bind_rows(res_entropies.magsep,.id = "Nr.") %>% 
  group_by(stage,aggregation,entropy_type) %>% 
  summarise(ent_mean = mean(entropy),ent_err = sd(entropy))
res.magsep_stage.agg <- group_by(res.magsep_stage.agg,aggregation) %>% mutate(new_ent_mean= cumsum(ent_mean))

type.levels <- c("ti00","ti0k","tij0","tijk", # partitions
                 "p00k","p0jk","pi0k","pijk" ,# distributions
                 "c0j0","cij0","c0jk","cijk") # compositions

g <- ggplot(res.magsep_stage.agg,aes(x = aggregation, y = ent_mean)) +
  geom_bar(stat = "identity",aes(fill = factor(entropy_type,levels = type.levels)))+
  geom_errorbar(aes(x = aggregation, ymin = new_ent_mean-ent_err, ymax = new_ent_mean+ent_err), width = 0.3) +
  facet_grid(.~stage, drop = T, scales = "free", space = "free") +
  scale_fill_manual(values = entropy_colors) +
  theme_bw() +
  labs(x = "aggregation", y = "entropy value", fill = "contribution")
plot(g)

res.magsep_stage.grouped <- pp_machine_entropy(machine = machines$magnetic_separator,part_dat = part_dat_grouped,n = 10,th.mass = 1)

res_entropies.magsep <- res.magsep_stage.grouped$entropies
res.magsep_stage.grouped.agg <- bind_rows(res_entropies.magsep,.id = "Nr.") %>% 
  group_by(stage,aggregation,entropy_type) %>% 
  summarise(ent_mean = mean(entropy),ent_err = sd(entropy))
res.magsep_stage.grouped.agg <- group_by(res.magsep_stage.grouped.agg,aggregation) %>% mutate(new_ent_mean= cumsum(ent_mean))

g <- ggplot(res.magsep_stage.grouped.agg,aes(x = aggregation, y = ent_mean)) +
  geom_bar(stat = "identity",aes(fill = factor(entropy_type,levels = type.levels)))+
  geom_errorbar(aes(x = aggregation, ymin = new_ent_mean-ent_err, ymax = new_ent_mean+ent_err), width = 0.3) +
  facet_grid(.~stage, drop = T, scales = "free", space = "free") +
  scale_fill_manual(values = entropy_colors) +
  theme_bw() +
  labs(x = "aggregation", y = "entropy value", fill = "contribution")
plot(g)

ggsave(filename = "magsep.entropies.grouped.paper2.png",device = "png",width = 20, height = 10, units = "cm")

####calculations for the mill####

res.mill <- batch_entropies(part_dat = part_dat_grouped,machine = machines$mill1,n = 10, th.mass = 1)


res.mill.agg <- map(res.mill,bind_rows,.id = "Nr.") %>% 
  bind_rows(.id = "sample") %>% 
  group_by(sample,aggregation,entropy_type) %>% 
  summarise(mean_entropy = mean(entropy), ent_err = sd(entropy))
res.mill.agg <- group_by(res.mill.agg,sample,aggregation) %>% mutate(new_ent_mean= cumsum(mean_entropy))

res.mill.agg <- filter(res.mill.agg, !(sample %in% c("recal_feed")))

res.mill.agg$sample <- str_replace_all(res.mill.agg$sample,c("S04" = "Feed", "S06" = "Product"))


i <- ggplot(res.mill.agg,aes(x = aggregation, y = mean_entropy)) +
  geom_bar(stat = "identity",aes(fill = factor(entropy_type,levels = type.levels)))+
  geom_errorbar(aes(x = aggregation, ymin = new_ent_mean-ent_err, ymax = new_ent_mean+ent_err), width = 0.3) +
  facet_grid(.~sample, drop = T, scales = "free", space = "free") +
  scale_fill_manual(values = entropy_colors) +
  theme_bw() +
  labs(x = "aggregation", y = "entropy value", fill = "contribution")
plot(i)

ggsave(filename = "mill.entropies.batches.png",device = "png",width = 15, height = 10, units = "cm")


#2)for the respective stage:

res.mill_stage <- pp_machine_entropy(machine = machines$mill1,part_dat = part_dat,n = 10,th.mass = 1)

res_entropies.mill <- res.mill_stage$entropies
res.mill_stage.agg <- bind_rows(res_entropies.mill,.id = "Nr.") %>% 
  group_by(stage,aggregation,entropy_type) %>% 
  summarise(ent_mean = mean(entropy),ent_err = sd(entropy))
res.mill_stage.agg <- group_by(res.mill_stage.agg,aggregation) %>% mutate(new_ent_mean= cumsum(ent_mean))

type.levels <- c("ti00","ti0k","tij0","tijk", # partitions
                 "p00k","p0jk","pi0k","pijk" ,# distributions
                 "c0j0","cij0","c0jk","cijk") # compositions

g <- ggplot(res.mill_stage.agg,aes(x = aggregation, y = ent_mean)) +
  geom_bar(stat = "identity",aes(fill = factor(entropy_type,levels = type.levels)))+
  geom_errorbar(aes(x = aggregation, ymin = new_ent_mean-ent_err, ymax = new_ent_mean+ent_err), width = 0.3) +
  facet_grid(.~stage, drop = T, scales = "free", space = "free") +
  scale_fill_manual(values = entropy_colors) +
  theme_bw() +
  labs(x = "aggregation", y = "entropy value", fill = "contribution")
plot(g)

res.mill_stage.grouped <- pp_machine_entropy(machine = machines$mill1,part_dat = part_dat_grouped,n = 10,th.mass = 1)

res_entropies.mill <- res.mill_stage.grouped$entropies
res.mill_stage.grouped.agg <- bind_rows(res_entropies.mill,.id = "Nr.") %>% 
  group_by(stage,aggregation,entropy_type) %>% 
  summarise(ent_mean = mean(entropy),ent_err = sd(entropy))
res.mill_stage.grouped.agg <- group_by(res.mill_stage.grouped.agg,aggregation) %>% mutate(new_ent_mean= cumsum(ent_mean))

type.levels <- c("ti00","ti0k","tij0","tijk", # partitions
                 "p00k","p0jk","pi0k","pijk" ,# distributions
                 "c0j0","cij0","c0jk","cijk") # compositions

g <- ggplot(res.mill_stage.grouped.agg,aes(x = aggregation, y = ent_mean)) +
  geom_bar(stat = "identity",aes(fill = factor(entropy_type,levels = type.levels)))+
  geom_errorbar(aes(x = aggregation, ymin = new_ent_mean-ent_err, ymax = new_ent_mean+ent_err), width = 0.3) +
  facet_grid(.~stage, drop = T, scales = "free", space = "free") +
  scale_fill_manual(values = entropy_colors) +
  theme_bw() +
  labs(x = "aggregation", y = "entropy value", fill = "contribution")
plot(g)

ggsave(filename = "mill.entropies.grouped.paper2.png",device = "png",width = 20, height = 10, units = "cm")

# calculations for the Spiral separator ####
#Berechnung für gruppierten Datensatz
res.spiralsep_stage <- pp_machine_entropy(machine = machines$spiral_separator,part_dat = part_dat_grouped,n = 10,th.mass = 1)

res_entropies.spiralsep <- res.spiralsep_stage$entropies
res.spiralsep_stage.agg <- bind_rows(res_entropies.spiralsep,.id = "Nr.") %>% 
  group_by(stage,aggregation,entropy_type) %>% 
  summarise(ent_mean = mean(entropy),ent_err = sd(entropy))
res.spiralsep_stage.agg <- group_by(res.spiralsep_stage.agg,aggregation) %>% mutate(new_ent_mean= cumsum(ent_mean))
save(res.spiralsep_stage.agg,file = "res.spiralsep.grouped.RData")

type.levels <- c("ti00","ti0k","tij0","tijk", # partitions
                 "p00k","p0jk","pi0k","pijk" ,# distributions
                 "c0j0","cij0","c0jk","cijk") # compositions

g <- ggplot(res.spiralsep_stage.agg,aes(x = aggregation, y = ent_mean)) +
  geom_bar(stat = "identity",aes(fill = factor(entropy_type,levels = type.levels)))+
  geom_errorbar(aes(x = aggregation, ymin = new_ent_mean-ent_err, ymax = new_ent_mean+ent_err), width = 0.3) +
  facet_grid(.~stage, drop = T, scales = "free", space = "free") +
  scale_fill_manual(values = entropy_colors) +
  theme_bw() +
  labs(x = "aggregation", y = "entropy value", fill = "contribution")
plot(g)

#Calculation for sulfide Flotation ####

res.sulfideflot_stage <- pp_machine_entropy(machine = machines$sulphide_flotation,part_dat = part_dat_grouped,n = 10,th.mass = 1)

res_entropies.sulfideflot <- res.sulfideflot_stage$entropies
res.sulfideflot_stage.agg <- bind_rows(res_entropies.sulfideflot,.id = "Nr.") %>% 
  group_by(stage,aggregation,entropy_type) %>% 
  summarise(ent_mean = mean(entropy),ent_err = sd(entropy))
res.sulfideflot_stage.agg <- group_by(res.sulfideflot_stage.agg,aggregation) %>% mutate(new_ent_mean= cumsum(ent_mean))
save(res.sulfideflot_stage.agg,file = "res.sulfideflot.grouped.RData")


g <- ggplot(res.sulfideflot_stage.agg,aes(x = aggregation, y = ent_mean)) +
  geom_bar(stat = "identity",aes(fill = factor(entropy_type,levels = type.levels)))+
  geom_errorbar(aes(x = aggregation, ymin = new_ent_mean-ent_err, ymax = new_ent_mean+ent_err), width = 0.3) +
  facet_grid(.~stage, drop = T, scales = "free", space = "free") +
  scale_fill_manual(values = entropy_colors) +
  theme_bw() +
  labs(x = "aggregation", y = "entropy value", fill = "contribution")
plot(g)

#calculations for screw classifier ####

res.screwclass_stage <- pp_machine_entropy(machine = machines$screw_classifier,part_dat = part_dat_grouped,n = 10,th.mass = 1)

res_entropies.screwclass <- res.screwclass_stage$entropies
res.screwclass_stage.agg <- bind_rows(res_entropies.screwclass,.id = "Nr.") %>% 
  group_by(stage,aggregation,entropy_type) %>% 
  summarise(ent_mean = mean(entropy),ent_err = sd(entropy))
res.screwclass_stage.agg <- group_by(res.screwclass_stage.agg,aggregation) %>% mutate(new_ent_mean= cumsum(ent_mean))
save(res.screwclass_stage.agg,file = "res.screwclass.grouped.RData")


g <- ggplot(res.screwclass_stage.agg,aes(x = aggregation, y = ent_mean)) +
  geom_bar(stat = "identity",aes(fill = factor(entropy_type,levels = type.levels)))+
  geom_errorbar(aes(x = aggregation, ymin = new_ent_mean-ent_err, ymax = new_ent_mean+ent_err), width = 0.3) +
  facet_grid(.~stage, drop = T, scales = "free", space = "free") +
  scale_fill_manual(values = entropy_colors) +
  theme_bw() +
  labs(x = "aggregation", y = "entropy value", fill = "contribution")
plot(g)

#Calculations for shaking table ####

res.shaking_table_stage <- pp_machine_entropy(machine = machines$shaking_table,part_dat = part_dat_grouped,n = 10,th.mass = 1)

res_entropies.shaking_table <- res.shaking_table_stage$entropies
res.shaking_table_stage.agg <- bind_rows(res_entropies.shaking_table,.id = "Nr.") %>% 
  group_by(stage,aggregation,entropy_type) %>% 
  summarise(ent_mean = mean(entropy),ent_err = sd(entropy))
res.shaking_table_stage.agg <- group_by(res.shaking_table_stage.agg,aggregation) %>% mutate(new_ent_mean= cumsum(ent_mean))
save(res.shaking_table_stage.agg,file = "res.shaking_table.grouped.RData")


g <- ggplot(res.shaking_table_stage.agg,aes(x = aggregation, y = ent_mean)) +
  geom_bar(stat = "identity",aes(fill = factor(entropy_type,levels = type.levels)))+
  geom_errorbar(aes(x = aggregation, ymin = new_ent_mean-ent_err, ymax = new_ent_mean+ent_err), width = 0.3) +
  facet_grid(.~stage, drop = T, scales = "free", space = "free") +
  scale_fill_manual(values = entropy_colors) +
  theme_bw() +
  labs(x = "aggregation", y = "entropy value", fill = "contribution")
plot(g)

#calculations mill2 ####

res.mill2_stage <- pp_machine_entropy(machine = machines$mill2,part_dat = part_dat_grouped,n = 10,th.mass = 1)

res_entropies.mill2 <- res.mill2_stage$entropies
res.mill2_stage.agg <- bind_rows(res_entropies.mill2,.id = "Nr.") %>% 
  group_by(stage,aggregation,entropy_type) %>% 
  summarise(ent_mean = mean(entropy),ent_err = sd(entropy))
res.mill2_stage.agg <- group_by(res.mill2_stage.agg,aggregation) %>% mutate(new_ent_mean= cumsum(ent_mean))
save(res.mill2_stage.agg,file = "res.mill2.grouped.RData")


g <- ggplot(res.mill2_stage.agg,aes(x = aggregation, y = ent_mean)) +
  geom_bar(stat = "identity",aes(fill = factor(entropy_type,levels = type.levels)))+
  geom_errorbar(aes(x = aggregation, ymin = new_ent_mean-ent_err, ymax = new_ent_mean+ent_err), width = 0.3) +
  facet_grid(.~stage, drop = T, scales = "free", space = "free") +
  scale_fill_manual(values = entropy_colors) +
  theme_bw() +
  labs(x = "aggregation", y = "entropy value", fill = "contribution")
plot(g)

#calculations hydrocyclone ####

res.hydrocyclone_stage <- pp_machine_entropy(machine = machines$hydrocyclone,part_dat = part_dat_grouped,n = 10,th.mass = 1)

res_entropies.hydrocyclone <- res.hydrocyclone_stage$entropies
res.hydrocyclone_stage.agg <- bind_rows(res_entropies.hydrocyclone,.id = "Nr.") %>% 
  group_by(stage,aggregation,entropy_type) %>% 
  summarise(ent_mean = mean(entropy),ent_err = sd(entropy))
res.hydrocyclone_stage.agg <- group_by(res.hydrocyclone_stage.agg,aggregation) %>% mutate(new_ent_mean= cumsum(ent_mean))
save(res.hydrocyclone_stage.agg,file = "res.hydrocyclone.grouped.RData")


g <- ggplot(res.hydrocyclone_stage.agg,aes(x = aggregation, y = ent_mean)) +
  geom_bar(stat = "identity",aes(fill = factor(entropy_type,levels = type.levels)))+
  geom_errorbar(aes(x = aggregation, ymin = new_ent_mean-ent_err, ymax = new_ent_mean+ent_err), width = 0.3) +
  facet_grid(.~stage, drop = T, scales = "free", space = "free") +
  scale_fill_manual(values = entropy_colors) +
  theme_bw() +
  labs(x = "aggregation", y = "entropy value", fill = "contribution")
plot(g)

#darstellung Entropie ####

entropy.res <- list("Sieve" = res.sieve_stage.grouped.agg,
                    "Magnetic separator" = res.magsep_stage.grouped.agg,
                    "Mill" = res.mill_stage.grouped.agg)
entropy.res <- bind_rows(entropy.res,.id = "machine")

entropy.res <- entropy.res %>% 
  filter(aggregation %in% c("ti00.cij0.pijk","ti00.pi0k.cijk","c0j0.p0jk","p00k.c0jk"))

entropy.res$entropy_type <- str_replace_all(entropy.res$entropy_type,c("c0j0" = "cij0","c0jk" = "cijk","p00k" = "pi0k","p0jk"="pijk"))

feed.res <- entropy.res %>% filter(stage == "feed")
prod.res <- entropy.res %>% filter(stage == "product")
by = join_by(machine,entropy_type)
entropy.diff <- right_join(feed.res,prod.res,by = by, suffix = c(".feed",".prod"))
entropy.diff <- entropy.diff %>% mutate(entropy.diff = (ent_mean.prod - ent_mean.feed)/ent_mean.feed*100)
entropy.diff <- na.omit(entropy.diff)

g <- ggplot(entropy.res,aes(x = entropy_type, y = entropy.diff)) +
  geom_bar(stat = "identity")+
 # geom_errorbar(aes(x = stage, ymin = ent_mean-ent_err, ymax = ent_mean+ent_err), width = 0.3) +
  facet_grid(~machine, drop = T, scales = "free", space = "free") +
  theme_bw() +
  labs(x = "contribution", y = "entropy value", fill = "contribution")
plot(g)
ggsave(filename = "entropy.differences.paper2.png",device = "png",width = 20, height = 10, units = "cm")

g <- ggplot(entropy.res,aes(x = stage, y = ent_mean)) +
  geom_bar(stat = "identity")+
  geom_errorbar(aes(x = stage, ymin = ent_mean-ent_err, ymax = ent_mean+ent_err), width = 0.3) +
  facet_grid(factor(machine,levels = machine_lvl)~entropy_type, drop = T, scales = "free", space = "free") +
  theme_bw() +
  labs(x = "stage", y = "entropy value", fill = "contribution")
plot(g)
ggsave(filename = "entropy.absolutevalues.paper2.png",device = "png",width = 20, height = 10, units = "cm")


# 
# save(res.sieve_stage,file = "res.sieve_stage.RData")
# save(res.magsep_stage,file = "res.magsep_stage.RData")
# save(res.mill_stage,file = "res.mill_stage.RData")
# save(res.sieve_stage.grouped,file = "res.sieve_stage.grouped.RData")
# save(res.magsep_stage.grouped,file = "res.magsep_stage.grouped.RData")
# save(res.mill_stage.grouped,file = "res.mill_stage.grouped.RData")

#Loading entropy data
library(tidyverse)
load(file = "./res.sieve_stage.grouped.RData")
load(file = "./res.spiralsep.grouped.RData")
load(file = "./res.mill_stage.grouped.RData")
load(file = "./res.magsep_stage.grouped.RData")
load(file = "./res.sulfideflot.grouped.RData")
load(file = "./res.spiralsep.grouped.RData")
load(file = "./res.mill2.grouped.RData")

res.entropies.pp <- list(sieve = res.sieve_stage.grouped,
                         spiral_separator = res.spiralsep_stage.agg,
                         mill = res_mill_stage.grouped,
                         magnetic_separator = res.magsep_stage.grouped,
                         sulfide_flotation = res.sulfideflot_stage.agg,
                         screw_classifier = res.screwclass_stage.agg,
                         shaking_table = )


## calcultaing particle size distributions ####

#1) for the sieve: 

samples <- c(SamplesMLA[c(machines$sieve$input,machines$sieve$output)])
idx <- res.sieve_stage$idx

psd.sieve <- psd(idx = idx, machine = machines$sieve)
psd.sieve <- bind_rows(psd.sieve, .id = "sample")
psd.sieve <- filter(psd.sieve, dist == "Q3" & sample != "feed")

psd.sieve$sample <- str_replace_all(psd.sieve$sample,c("S02" = "Fine", "S03" = "Coarse","recal_feed" = "Feed"))

ggplot(psd.sieve,aes(log10(s.classes),Qi,col = sample)) +
  geom_line()

idx <- res.magsep_stage$idx

#2) for the magnetic separation: 

idx <- res.magsep_stage$idx

psd.magsep <- psd(idx = idx, machine = machines$magnetic_separator)
psd.magsep <- bind_rows(psd.magsep, .id = "sample")
psd.magsep <- filter(psd.magsep, dist == "Q3" & sample != "feed")

psd.magsep$sample <- str_replace_all(psd.magsep$sample,c("S07" = "Non-magnetic", "S08" = "Magnetic","recal_feed" = "Feed"))

ggplot(psd.magsep,aes(log10(s.classes),Qi,col = sample)) +
  geom_line()

#3) for the mill: 

idx <- res.mill_stage$idx

psd.mill <- psd(idx = idx, machine = machines$mill1)
psd.mill <- bind_rows(psd.mill, .id = "sample")
psd.mill <- filter(psd.mill, dist == "Q3" & sample != "recal_feed")

psd.mill$sample <- str_replace_all(psd.mill$sample,c("S06" = "Product","recal_feed" = "Feed", "feed" = "Feed"))

ggplot(psd.mill,aes(log10(s.classes),Qi,col = sample)) +
  geom_line()

psd.res <- list("Sieve" = psd.sieve,"Magnetic separator" = psd.magsep, "Mill" = psd.mill )
psd.res <- bind_rows(psd.res, .id = "machine")

g <- ggplot(psd.res, aes(x = log10(s.classes), y = Qi, col = factor(sample,levels = sample_lvl))) +
  geom_point(size = 1.5) +
  geom_line() +
  facet_grid(.~factor(machine,levels = machine_lvl))+
  theme_bw() +
  scale_color_brewer(type = "qual",palette = "Set1") +
  labs(x = "log(x) in µm", y = "Q3 in %",col = "Sample") +
  theme(text = element_text(size=14))
  
print(g)

ggsave(filename = "psd.paper2.png",device = "png",width = 20, height = 10, units = "cm")

# calculating component ####

#1) for sieve:

idx <- res.sieve_stage$idx

comp.sieve <- comp_function(idx = idx, machine = machines$sieve, part_dat = part_dat_grouped)
comp.sieve <- bind_rows(comp.sieve, .id = "sample")
comp.sieve <- filter(comp.sieve, sample != "feed")

comp.sieve$sample <- str_replace_all(comp.sieve$sample,c("S02" = "Fine", "S03" = "Coarse","recal_feed" = "Feed"))

ggplot(comp.sieve,aes(sample,mean_comp,fill = Mineral)) +
  geom_bar(stat = "identity")

#2) for magnetic separation:

idx <- res.magsep_stage$idx

comp.magsep <- comp_function(idx = idx, machine = machines$magnetic_separator, part_dat = part_dat_grouped)
comp.magsep <- bind_rows(comp.magsep, .id = "sample")
comp.magsep <- filter(comp.magsep, sample != "feed")

comp.magsep$sample <- str_replace_all(comp.magsep$sample,c("S07" = "Non-magnetic", "S08" = "Magnetic","recal_feed" = "Feed"))

ggplot(comp.magsep,aes(sample,mean_comp,fill = Mineral)) +
  geom_bar(stat = "identity")

idx <- res.mill_stage$idx

comp.mill <- comp_function(idx = idx, machine = machines$mill1, part_dat = part_dat_grouped)
comp.mill <- bind_rows(comp.mill, .id = "sample")
comp.mill <- filter(comp.mill, sample != "recal_feed")

comp.mill$sample <- str_replace_all(comp.mill$sample,c("S06" = "Product","recal_feed" = "Feed", "feed" = "Feed"))

ggplot(comp.mill,aes(sample,mean_comp,fill = Mineral)) +
  geom_bar(stat = "identity")

comp.res <- list("Sieve" = comp.sieve,"Magnetic separator" = comp.magsep, "Mill" = comp.mill )
comp.res <- bind_rows(comp.res, .id = "machine")

g <- ggplot(comp.res,aes(factor(sample,levels = sample_vec),mean_comp,fill = Mineral)) +
  geom_bar(stat = "identity") +
  facet_grid(.~factor(machine, levels = machine_lvl),scales = "free_x") +
  theme_classic() +
  labs(x = "",y = "Composition in %") +
  scale_fill_brewer(type = "qual",palette = "Set2") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size=14)) 
print(g)
ggsave(filename = "comps.paper2.png",device = "png",width = 20, height = 10, units = "cm")

# plotting batch entropies together: 

plots = list("sieve" = g,
             "magnetic separator" = h,
             "mill" = i)

ggarrange(plotlist = plots,ncol = 2, nrow = 2,
          common.legend = T,
          legend = "right",
          hjust = c(-2,-0.4,-1),
          vjust = 0.5,
          heights = c(0.7,0.7,0.7),
          legend.grob = get_legend(plots),
          widths = c(1,1,1))


ggsave(filename= "flowsheet.png",device = "png", units = "cm",dpi = 300,width = 20,height = 10)

batch.entropies <- list(
  "Sieve" = res.sieve.agg,
  "Magnetic Separator" = res.magsep.agg,
  "Mill" = res.mill.agg
)

batch.entropies <- batch.entropies %>% bind_rows(.id = "machine")

g <- ggplot(batch.entropies,aes(x = aggregation, y = mean_entropy)) +
  geom_bar(stat = "identity",aes(fill = factor(entropy_type,levels = type.levels))) +
  geom_errorbar(aes(x = aggregation, ymin = new_ent_mean-ent_err, ymax = new_ent_mean+ent_err), width = 0.3) +
  facet_grid(factor(machine,levels = machine_lvl)~factor(sample,levels = sample_lvl), drop = T, scales = "free", space = "free") +
  scale_fill_manual(values = entropy_colors) +
  theme_bw() +
  labs(x = "aggregation", y = "entropy value", fill = "contribution")
plot(g)

#this is just a test for github - lets see if it works
