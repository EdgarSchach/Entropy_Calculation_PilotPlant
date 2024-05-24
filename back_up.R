machine <- machines$sieve

a <- particle.numbers[machine$output]
b <- clo(machine$weights[machine$output])
clo(a*b)*100000

machine <- machines$sieve
entropies <- Entropies.agg
machines$screw_classifier$input <- c("S09","S018")
machines$screw_classifier$weights <- c(machines$screw_classifier$weights,"S09" = 372.66,"S018" = 134)

machine.entropy <- function(machine,entropies){
  entropies_input <- entropies %>% filter(Sample %in% c(machine$input))
  weights_input <- tibble(Sample = machine$input,
                          weights = machine$weights[machine$input])
  if (length(machine$input == 1)) {
    weights_input$weights = 1
  }
  input <- right_join(entropies_input,weights_input,by = "Sample") %>% group_by(type,aggregation) %>% 
    mutate(weighted_entropy = weights/sum(weights)*mean_entropy) %>% 
    summarise(entropy = sum(weighted_entropy)) %>% 
    ungroup() %>% 
    add_row( type = "ti00", aggregation = "c0j0.p0jk", entropy = .entropy(clo(machine$weights[machine$input]))) %>% 
    add_row( type = "ti00", aggregation = "p00k.c0jk", entropy = .entropy(clo(machine$weights[machine$input]))) 
  input$type <-  str_replace_all(input$type,c("c0j0"="cij0","c0jk"="cijk","p00k"="pi0k","p0jk"="pijk"))
  input$aggregation <-  str_replace_all(input$aggregation,c("c0j0"="ti00.cij0","c0jk"="cijk","p00k"="ti00.pi0k","p0jk"="pijk")) 
  
  entropies_output <- entropies %>% filter(Sample %in% c(machine$output))
  weights_output <- tibble(Sample = machine$output,
                           weights = machine$weights[machine$output])
  if (length(machine$output == 1)) {
    weights_output$weights = 1
  }
  output <- right_join(entropies_output,weights_output) %>% group_by(type,aggregation) %>% 
    mutate(weighted_entropy = weights/sum(weights)*mean_entropy) %>% 
    summarise(entropy = sum(weighted_entropy)) %>% 
    ungroup() %>% 
    add_row( type = "ti00", aggregation = "c0j0.p0jk", entropy = .entropy(clo(machine$weights[machine$output]))) %>% 
    add_row( type = "ti00", aggregation = "p00k.c0jk", entropy = .entropy(clo(machine$weights[machine$output]))) 
  output$type <-  str_replace_all(output$type,c("c0j0"="cij0","c0jk"="cijk","p00k"="pi0k","p0jk"="pijk"))
  output$aggregation <-  str_replace_all(output$aggregation,c("c0j0"="ti00.cij0","c0jk"="cijk","p00k"="ti00.pi0k","p0jk"="pijk")) 
  res <- list(feed = input,
              product = output)
  res <- bind_rows(res,.id = "stage")
  
  
  return(res)
}

entropies_pp_ungrouped <- map(machines,~machine.entropy(machine = .x,entropies = Entropies.agg))
entropies_pp_ungrouped <- bind_rows(entropies_pp_ungrouped,.id = "machine") 
entropies_pp_ungrouped$entropy[is.na(entropies_pp_ungrouped$entropy)] <- 0


g <- ggplot(entropies_pp_ungrouped,aes(x = aggregation,y= entropy, fill = type)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  facet_grid(factor(machine,levels = machine.lvl.extended)~stage, scales = "free_x") +
  scale_fill_manual(values = entropy_colors)
print(g)  
ggsave(filename = "entropy.machines.normalized.png",device = "png",width = 20, height = 30, units = "cm")

# calculating absolute differences for machines

feed.res <- entropies_pp_ungrouped %>% filter(stage == "feed")
prod.res <- entropies_pp_ungrouped %>% filter(stage == "product")
by = join_by(machine,type,aggregation)
entropy.diff <- right_join(feed.res,prod.res,by = by, suffix = c(".feed",".prod"))
entropy.diff <- entropy.diff %>% mutate(entropy.diff = (entropy.prod - entropy.feed))
entropy.diff <- na.omit(entropy.diff)

entropy.diff$machine <- str_replace_all(entropy.diff$machine,machine.nms)

g <- ggplot(entropy.diff,aes(x = type, y = entropy.diff,fill = type)) +
  geom_bar(stat = "identity")+
  # geom_errorbar(aes(x = stage, ymin = ent_mean-ent_err, ymax = ent_mean+ent_err), width = 0.3) +
  facet_grid(~factor(machine,levels = machine.lvl.extended), drop = T, scales = "free", space = "free") +
  scale_fill_manual(values = entropy_colors)+
  theme_bw() +
  labs(x = "contribution", y = "entropy value", fill = "contribution")
plot(g)
ggsave(filename = "entropy.differences.paper2.png",device = "png",width = 20, height = 10, units = "cm")


machine <- machines$mill1
part.dat <- part_dat
part.number <- 100000
entropies <- Entropies.agg

machine_entropy <- function(machine,part.dat,entropies,part.number = 100000){
  
  #1) entropy of recalculated feed:####
  
  # getting particle data for output streams:
  if (length(machine$output) > 1){
    output <- part.dat[machine$output]
    part.number = sum(particle.numbers[machine$output]*clo(machine$weights[machine$output]))
  }else{
    output = part.dat[machine$input]
    part.number = particle.numbers[machine$input]
  }
  
  # particle weights: mass fraction product * (1/number of particles)
  pw <- map2(clo(machine$weights[machine$output]),sapply(output,nrow),~.x/.y)
  #pw <- clo(clo(machine$weights[machine$output])/mean_masses[machine$output])
  pw <- map(output,~1/nrow(.x))
  # particle numbers:
  pn <- lapply(output,nrow)
  # weights as vector of pn*pw for the different samples:
  w <- unlist(map2(pw,pn,~c(rep(.x,.y))))
  # gluing the product particles together: 
  output <- bind_rows(output)
  
  #function for parallized calculation of entropies for bootstraps - 100000 particles seem to be a reliant number:
  
  batch_entropy.parallel <- function(MLA.sample,part.numbers,w = 1){
    idx <- sample.int(nrow(MLA.sample),size = part.numbers,replace = T,prob = w)
    pd <- MLA.sample[idx,]
    entropy = batch.entropy.norm(pd)
    return(entropy)
  }
  
  # preparing parallel programming: 
  
  availableCores()
  plan(multisession,workers = availableCores() -1)
  # calculating entropies : 
  recal.feed.entropy <- future_replicate(20,batch_entropy.parallel(MLA.sample = output,part.numbers = 100000,w = w),future.seed = T,simplify = F)
  
  # aggregation of entropies for boots:
  recal.feed.agg <- recal.feed.entropy %>% bind_rows(.id = "N") %>% 
    group_by(entropy_type,aggregation) %>% 
    summarize(mean_entropy = mean(entropy), sd_entropy = sd(entropy))
  
  #rescaling the particle entropy with the number of particles
  if(length(machine$input == 1)) {
    weight_input = 1
  }else{
    weight_input = clo(machine$weights[machine$input])
  }
  
  
  recal.part.entropy <- filter(recal.feed.agg,entropy_type %in% c("p00k","p0jk")) %>% 
    mutate(mean_entropy = mean_entropy + log(part.number))
  recal.feed.agg <- rows_update(recal.feed.agg,recal.part.entropy,by = "entropy_type") %>% 
    ungroup() %>% 
    add_row( entropy_type = "ti00", aggregation = "c0j0.p0jk", mean_entropy = .entropy(weight_input)) %>% 
    add_row( entropy_type = "ti00", aggregation = "p00k.c0jk", mean_entropy = .entropy(weight_input)) 
  recal.feed.agg$entropy_type <-  str_replace_all(recal.feed.agg$entropy_type,c("c0j0"="cij0","c0jk"="cijk","p00k"="pi0k","p0jk"="pijk"))
  recal.feed.agg$aggregation <-  str_replace_all(recal.feed.agg$aggregation,c("c0j0"="ti00.cij0","c0jk"="cijk","p00k"="ti00.pi0k","p0jk"="pijk"))
  recal.feed.agg <-  rename(recal.feed.agg,c("type" = "entropy_type","entropy" = "mean_entropy")) %>% select(-sd_entropy)
  
  #calculating entropies for the product stage weigthed with ti00:
  entropies_output <- entropies %>% filter(Sample %in% c(machine$output))
  weights_output <- tibble(Sample = machine$output,
                           weights = machine$weights[machine$output])
  if (length(machine$output == 1)) {
    weights_output$weights = 1
  }
  output <- right_join(entropies_output,weights_output) %>% group_by(type,aggregation) %>% 
    mutate(weighted_entropy = weights/sum(weights)*mean_entropy) %>% 
    summarise(entropy = sum(weighted_entropy)) %>% 
    ungroup() %>% 
    add_row( type = "ti00", aggregation = "c0j0.p0jk", entropy = .entropy(clo(machine$weights[machine$output]))) %>% 
    add_row( type = "ti00", aggregation = "p00k.c0jk", entropy = .entropy(clo(machine$weights[machine$output]))) 
  output$type <-  str_replace_all(output$type,c("c0j0"="cij0","c0jk"="cijk","p00k"="pi0k","p0jk"="pijk"))
  output$aggregation <-  str_replace_all(output$aggregation,c("c0j0"="ti00.cij0","c0jk"="cijk","p00k"="ti00.pi0k","p0jk"="pijk"))
  
  res <- output %>% rename("entropy.product" = "entropy") %>% 
    mutate(entropy.feed = recal.feed.agg$entropy,
           entropy.diff = entropy.product - entropy.feed)
  return(res)
  
}
output %>% group_by(aggregation) %>% summarise(total_entropy = sum(entropy))
recal.feed.agg %>% group_by(aggregation) %>% summarise(total_entropy = sum(entropy))


machine.entropies <- map(machines,~machine_entropy(.x,entropies = Entropies.agg,part.dat = part_dat))
### plotting log N relations: 

load("./error.results.RData")

machine.entropies.agg <- machine.entropies %>% bind_rows(.id = "machine")

g <- ggplot(machine.entropies.agg, aes(x = type, y = entropy.diff,fill = type)) +
  geom_bar(stat = "identity") +
  facet_grid(.~machine) +
  scale_fill_manual(values = entropy_colors)+
  theme_bw() +
  labs(x = "contribution", y = "entropy value", fill = "contribution")
plot(g)