sample_point <- "S01"
part.dat <- part_dat_grouped$S01
batch.entropy.norm <- function(part.dat,sample_point){

  # get particle data from name
  con <- grouped_concentrations[sample_point,]/100
  
  con_weighting <- 1/(1-con$Rest)
  con <- select(con,-Rest)
  # Calculating Extensive Properties for Feed Batch:
  
  N000 <- sum(part.dat)
  N0j0 <- colSums(part.dat)
  N00k <- rowSums(part.dat)
  N0jk <- part.dat
  
  #Calculating Extensive Properties for Product Batch(es)
  
  #Intensive Größen für Feed:   
  
  #c0j0 <- colSums(part.dat)/sum(part.dat) # composition of batchstage
  c0j0 <- unlist(con)
  
  
  c0jk <- N0jk/N00k # integrowth of batchstage
  p00k <- N00k/N000 # mass distribution in batchstage
  p0jk <- sweep(N0jk,2,N0j0,"/") #component distribution in batchstage
  
  #Entropy Calculation for Feed####
  #c0j0.p0jk
  
  ent_c0j0 <- .entropy(c0j0)*con_weighting #- log(length(c0j0))
  ent_p0jk <- apply(p0jk,2,.entropy)
  ent_p0jk <- ent_p0jk - log(length(p00k))
  ent_p0jk <- ent_p0jk[names(ent_p0jk) %in% of.interest]  
  ent_p0jk <- c0j0 %*% ent_p0jk
  
  #p00k.cojk
  
  ent_p00k <- .entropy(p00k)
  ent_p00k <- ent_p00k - log(length(p00k))
  ent_c0jk <- apply(c0jk,1,.entropy)# - log(length(c0j0))
  ent_c0jk <- p00k %*% ent_c0jk
  
  c0j0.p0jk <- c(c0j0 = ent_c0j0, p0jk = ent_p0jk)
  p00k.c0jk <- c(p00k = ent_p00k, c0jk = ent_c0jk)
  
  res_batch <- list(c0j0.p0jk = c0j0.p0jk,
                    p00k.c0jk = p00k.c0jk)
  
  res_batch <- bind_rows(res_batch,.id = "aggregation") %>% 
    gather(key = "entropy_type",value = "entropy",-aggregation) %>% na.omit()
  
  
  return(res_batch)
}
particle_data <- part_dat_grouped
Sample <- names(part_dat)[1]
sample_points <- names(part_dat)
sample_point <- sample_points[1]
names(sample_points) <- names(part_dat)
part.numbers <- 1000



batch_entropy.parallel <- function(sample_point,part.numbers,particle_data){
  pd <- particle_data[[sample_point]]
  idx <- sample(1:nrow(pd),size = part.numbers,replace = T)
  pd <- pd[idx,]
  entropy = batch.entropy.norm(pd,sample_point = sample_point)
  return(entropy)
}

test <- map(sample_points,~batch_entropy.parallel(sample_point = .x,part.numbers = 10000,particle_data = part_dat))


entropyforboots.particles <- function(feed,products,weights){
  pl <- length(products)
  # Calculating Extensive Properties for Feed Batch:
  feed <- feed[[1]]
  F000 <- sum(feed)
  F0j0 <- colSums(feed)
  F00k <- rowSums(feed)
  F0jk <- feed
  
  #Calculating Extensive Properties for Product Batch(es)
  
  if (pl > 1){
    Pi00 <- lapply(products,sum)
    Pi0k <- lapply(products,rowSums)
    Pij0 <- lapply(products,colSums) 
    Pijk <- products
  }else{
    products <- products[[1]]
    P000 <- sum(products)
    P0j0 <- colSums(products)
    P00k <- rowSums(products)
    P0jk <- products
  }
  #Intensive Größen für Feed:   
  
  c0j0 <- colSums(feed)/sum(feed) # compositioF of batchstage
  c0jk <- F0jk/F00k # integrowth of batchstage
  p00k <- F00k/F000 # mass distribution in batchstage
  p0jk <- sweep(F0jk,2,F0j0,"/") #component distribution in batchstage
  
  #Entropy Calculation for Feed####
  #c0j0.p0jk
  
  ent_c0j0 <- .entropy(c0j0)
  ent_p0jk <- apply(p0jk,2,.entropy)
  ent_p0jk <- c0j0 %*% ent_p0jk
  ent_p0jk.norm <- ent_p0jk - log(nrow(feed))
  
  
  #p00k.cojk
  
  ent_p00k <- .entropy(p00k) 
  ent_p00k.norm <- .entropy(p00k) - log(length(p00k))
  ent_c0jk <- apply(c0jk,1,.entropy)
  ent_c0jk <- p00k %*% ent_c0jk
  
  ti00.cij0.pijk <- c(ti00=  0, cij0 = ent_c0j0, pijk = ent_p0jk, p0jk.norm = ent_p0jk.norm)
  ti00.pi0k.cijk <- c(ti00 = 0, pi0k = ent_p00k, cijk = ent_c0jk, p00k.norm = ent_p00k.norm)
  
  res_feedstage <- list(ti00.cij0.pijk = ti00.cij0.pijk,
                        ti00.pi0k.cijk = ti00.pi0k.cijk)
  res_feedstage <- bind_rows(res_feedstage,.id = "aggregation") %>% 
    gather(key = "entropy_type",value = "entropy",-aggregation) %>% na.omit()
  
  # Calculation of Product Stage 
  
  # Calculation of intensive properties for product stages ####
  
  if (pl == 1){
    cij0 <- colSums(feed)/sum(feed) # composition of productstage
    cijk <- P0jk/P00k # integrowth of productstage
    pi0k <- P00k/P000 # mass distribution in productstage
    pijk <- sweep(P0jk,2,P0j0,"/") #component distribution in productstage
    
    #c0j0.p0jk
    
    ent_cij0 <- .entropy(cij0)
    ent_pijk <- apply(pijk,2,.entropy)
    ent_pijk <- cij0 %*% ent_pijk
    ent_pijk.norm <- ent_pijk - log(length(pijk))
    
    #p00k.c0jk
    
    ent_pi0k <- .entropy(pi0k) 
    ent_pi0k.norm <- .entropy(pi0k) - log(length(pi0k))
    ent_cijk <- apply(cijk,1,.entropy)
    ent_cijk <- pi0k %*% ent_cijk
    
    ti00.cij0.pijk <- c(ti00 = 0,cij0 = ent_cij0, pijk = ent_pijk, pijk.norm = ent_pijk.norm)
    ti00.pi0k.cijk <- c(ti00 = 0,pi0k = ent_pi0k, cijk = ent_cijk, pi0k.norm = ent_pi0k.norm)
    
    res_productstage <- list(ti00.cij0.pijk = ti00.cij0.pijk,
                             ti00.pi0k.cijk = ti00.pi0k.cijk)
    res_productstage <- bind_rows(res_productstage,.id = "aggregation") %>% 
      gather(key = "entropy_type",value = "entropy",-aggregation) %>% na.omit()
  }else{
    P0j0 <- Reduce("+",Pij0)   # recalculated mass of minerals from the products
    cij0 <- map(Pij0,~.x/sum(.x)) #composition of product stage 
    cijk <- map2(Pijk,Pi0k,~sweep(.x,1,.y,"/")) #integrowth of product stage
    ti00 <- unlist(map_dfr(products,~sum(.x)/sum(feed))) #mass partition into products
    
    tij0 <- map_dfc(Pij0,~.x/P0j0)
    pi0k <- map(Pi0k,~.x/sum(.x))  # mass distribution in product stage
    pijk <- map2(Pijk, Pij0, ~ sweep(.x,2,.y,"/")) #component distribution in product stage
    pijk <- map(pijk,~mutate_all(.x,~replace(., is.na(.), 0)))
    
    #ti00.cij0.pijk
    
    ent_ti00 <- .entropy(ti00)
    ent_cij0 <- ti00 %*% sapply(cij0,.entropy)
    ent_pijk <- lapply(pijk,function(x)apply(x,2,.entropy))
    ent_pijk <- map2(cij0,ent_pijk,~.x%*%.y,)
    #ent_pijk <- map2(w,ent_pijk,~.x*.y)
    ent_pijk <- ti00 %*% unlist(ent_pijk)
    
    ti00.cij0.pijk <- c(ti00 = ent_ti00, cij0 = ent_cij0, pijk = ent_pijk)
    
    #c0j0.tij0.pijk
    
    ent_c0j0 <- .entropy(c0j0)
    ent_tij0 <- apply(tij0,1,.entropy)
    ent_tij0 <- c0j0 %*% ent_tij0
    
    c0j0.tij0.pijk <- c(c0j0 = ent_c0j0,tij0 = ent_tij0,pijk = ent_pijk)
    
    #ti00.pi0k.cijk
    
    ent_ti00 <- .entropy(ti00)
    ent_pi0k <- sapply(pi0k,.entropy)
    #ent_pi0k <- w*sapply(pi0k,.entropy)
    ent_pi0k <- ti00 %*% ent_pi0k
    ent_cijk <- lapply(cijk,function(x)apply(x,1,.entropy))
    ent_cijk <- map2(pi0k,ent_cijk,~.x%*%.y)
    ent_cijk <- ti00 %*% unlist(ent_cijk)
    ti00.pi0k.cijk <- c(ti00 = ent_ti00,pi0k = ent_pi0k,cijk = ent_cijk)
    
    
    
    
    res_productstage <- list(c0j0.tij0.pijk = c0j0.tij0.pijk,
                             ti00.cij0.pijk = ti00.cij0.pijk,
                             ti00.pi0k.cijk = ti00.pi0k.cijk)
    res_productstage <- bind_rows(res_productstage,.id = "aggregation") %>% 
      gather(key = "entropy_type",value = "entropy",-aggregation) %>% na.omit()
    
    
  }
  
  res <-  bind_rows(feed= res_feedstage,product = res_productstage,.id = "stage")
  return(res)
}


pp_machine_entropy <-  function(machine, part_dat, n, th.mass = NA,SamplesMLA){
  
  samples <- SamplesMLA[c(machine$input,machine$output)]
  
  
  relevants <- function(x) data.frame(ECD = sqrt(4*x$geomProps$ParticleArea/pi),
                                      mass = x$geomProps$Mass,
                                      dens = x$geomProps$Mass/x$geomProps$ParticleArea/100)
  rel <- lapply(samples,relevants)
  
  #calculating mean Volume based on ECD: 
  
  
  cal_mean_d <- function(x){
    volume <- pi/6*x$ECD^3
    mean_d <- sum(volume*x$ECD)/sum(volume)
    mean_dens <- sum(volume*x$dens)/sum(volume)
    mean_mass <- 1/6*pi*mean_d^3 * mean_dens * 10^(-12)
    res <- data.frame(ecd = mean_d,
                      dens = mean_dens,
                      mass = mean_mass)
    return(res)
  }
  mean_d <- sapply(rel,cal_mean_d)
  mean_m <- mean_d["mass",]
  
  
  
  # calcuating mass for bootstraps
  if(is.na(th.mass)){
    w <- c(Feed = sum(machine$weights[machine$output]),machine$weights[machine$output])
    names(w) <- names(samples)
    m <- w/1000
  }else{
    w <- c(rep(th.mass,length(c(machine$input,machine$output))))
    names(w) <- names(samples)
    if(length(machine$input) == 1){
      m <- w*c(1,machine$weights[machine$output]/sum(machine$weights[machine$output]))
    }else{
      m <- w*c(machine$weights[machine$input]/sum(machine$weights[machine$input]),machine$weights[machine$output]/sum(machine$weights[machine$output]))
    }
  }
  N <- m/unlist(mean_m[names(m)])
  
  
  #create indexes for each sample of the product - gives a matrix where each column contain indexes of one sample
  
  idx <- map2(samples,N,~replicate(n,sample(seq(1,nrow(.x)),size=.y, replace = T)))
  samples <- part_dat[c(machine$input,machine$output)]
  boots <- map2(samples,idx,~ apply(.y,2,function(y).x[y,])) 

  
  boots_p <-  transpose(boots[machine$output])
  boots_f <-  transpose(boots[machine$input])
  
  remove(boots)


  boots_f <-  lapply(boots_f,bind_rows)
  # getting the indexes for the feed sample: 
  boots_f_recal <- lapply(boots_p,bind_rows)
  
  
  #calculating the entropies using the entropy for boots function:
  
  if(length(machine$output)>1){
    booted_entropies <- map2(boots_f_recal,boots_p,entropyforboots.particles,weights = machine$weights) 
  }else{
    booted_entropies <- map2(boots_f,boots_p,entropyforboots.particles,weights = machine$weights) 
  }
  
  #booted_entropies <- mapply(entropyforboots, boots_f,boots_w,MoreArgs = list(weights = machine$weights) , SIMPLIFY = F)
  res = list(
    idx = idx,
    entropies = booted_entropies,
    num_cal = mean_d,
    numbers = N)
  return(res)
}



batch_entropies <- function(part_dat,machine,n,th.mass = NA,SamplesMLA){
  

  samples <- SamplesMLA[c(machine$input,machine$output)]
  
  
  relevants <- function(x) data.frame(ECD = sqrt(4*x$geomProps$ParticleArea/pi),
                                      mass = x$geomProps$Mass,
                                      dens = x$geomProps$Mass/x$geomProps$ParticleArea/100)
  rel <- lapply(samples,relevants)
  
  #calculating mean Volume based on ECD: 
  
  
  cal_mean_d <- function(x){
    volume <- pi/6*x$ECD^3
    mean_d <- sum(volume*x$ECD)/sum(volume)
    mean_dens <- sum(volume*x$dens)/sum(volume)
    mean_mass <- 1/6*pi*mean_d^3 * mean_dens * 10^(-12)
    res <- data.frame(ecd = mean_d,
                      dens = mean_dens,
                      mass = mean_mass)
    return(res)
  }
  mean_d <- sapply(rel,cal_mean_d)
  mean_m <- mean_d["mass",]
  
  
  
  # calcuating mass for bootstraps
  if(is.na(th.mass)){
    w <- c(Feed = sum(machine$weights[machine$output]),machine$weights[machine$output])
    names(w) <- names(samples)
    m <- w/1000
  }else{
    w <- c(rep(th.mass,length(c(machine$input,machine$output))))
    names(w) <- names(samples)
    m <- w*c(machine$weights[machine$input]/sum(machine$weights[machine$input]),machine$weights[machine$output]/sum(machine$weights[machine$output]))
  }
  N <- m/unlist(mean_m[names(m)])
  
  idx <- map2(samples,N,~replicate(n,sample(seq(1,nrow(.x)),size=.y, replace = T)))
  samples <- part_dat[c(machine$input,machine$output)]
  boots <- map2(samples,idx,~ apply(.y,2,function(y).x[y,])) 
  
  boots_p <-  transpose(boots[machine$output])
  # getting the indexes for the feed sample: 
  boots_f_recal <- lapply(boots_p,bind_rows)

  entropy_p <- lapply(boots,function(x)lapply(x,batch.entropy.cal))  
  entropy_f <- lapply(boots_f_recal,batch.entropy.cal)
  entropy_p$recal_feed <- entropy_f

  return(entropy_p)

  }

batch.entropy.cal <- function(part.dat){
    # Calculating Extensive Properties for Feed Batch:
  
  N000 <- sum(part.dat)
  N0j0 <- colSums(part.dat)
  N00k <- rowSums(part.dat)
  N0jk <- part.dat
  
  #Calculating Extensive Properties for Product Batch(es)
  
  #Intensive Größen für Feed:   
  
  c0j0 <- colSums(part.dat)/sum(part.dat) # compositioF of batchstage
  c0jk <- N0jk/N00k # integrowth of batchstage
  p00k <- N00k/N000 # mass distribution in batchstage
  p0jk <- sweep(N0jk,2,N0j0,"/") #component distribution in batchstage
  
  #Entropy Calculation for Feed####
  #c0j0.p0jk
  
  ent_c0j0 <- .entropy(c0j0)
  ent_p0jk <- apply(p0jk,2,.entropy)
  ent_p0jk <- c0j0 %*% ent_p0jk
  
  
  #p00k.cojk
  
  ent_p00k <- .entropy(p00k)
  ent_c0jk <- apply(c0jk,1,.entropy)
  ent_c0jk <- p00k %*% ent_c0jk
  
  c0j0.p0jk <- c(c0j0 = ent_c0j0, p0jk = ent_p0jk)
  p00k.c0jk <- c(p00k = ent_p00k, c0jk = ent_c0jk)
  
  res_batch <- list(c0j0.p0jk = c0j0.p0jk,
                        p00k.c0jk = p00k.c0jk)
  res_batch <- bind_rows(res_batch,.id = "aggregation") %>% 
    gather(key = "entropy_type",value = "entropy",-aggregation) %>% na.omit()
}


psd <- function(idx,machine,SamplesMLA){
  
  samples <- SamplesMLA[c(machine$input,machine$output)]
  
  relevants <- function(x) data.frame(ECD = sqrt(4*x$geomProps$ParticleArea/pi),
                                      mass = x$geomProps$Mass,
                                      dens = x$geomProps$Mass/x$geomProps$ParticleArea/100)
  
  ll <- lapply(samples,relevants)
  
  parts <-map2(idx,ll,~apply(.x,MARGIN = 2,function(x).y[x,]))
  
  parts_p <- parts[machine$output]
  parts_p <- transpose(parts_p)
  parts$recal_feed = map(parts_p,bind_rows)
  
  if(length(machine$input) > 1){
    parts_f <- parts[machine$input]
    parts_f <- transpose(parts_f)
    parts$feed = map(parts_f,bind_rows)
    parts[machine$input] <- NULL
  }else{
    parts$feed = parts[[machine$input]]
    parts[machine$input] <- NULL
  }
  
  psd2d <- function(x){
    res <- x %>% 
      mutate(volume = pi/6*x$ECD^3,
             mass_q3 = volume*dens) %>% 
      mutate(s.classes = cut(x$ECD,breaks = c(0,10,20,40,80,160,320,640,800,1000,1300,1800,2600),
                             labels = c("10","20","40","80","160","320","640","800","1000","1300","1800","2600"))) %>% 
      group_by(s.classes) %>% 
      summarise(m_q2 = sum(mass),m_q3 = sum(mass_q3)) %>% 
      mutate( q2= cumsum(clo(m_q2))*100, q3 = cumsum(clo(m_q3))*100)
    
    return(res)
  }
  res <- lapply(parts,function(x)lapply(x,psd2d))
  
  agg.func <- function(psd.dat){
    res <- bind_rows(psd.dat,.id = "N")
    res$s.classes <- as.numeric(as.character(res$s.classes))
    res_q2 <- res %>% 
      group_by(s.classes) %>% 
      summarize(Qi = mean(q2), sdi = sd(q2))
    res_q3 <- res %>% 
      group_by(s.classes) %>% 
      summarize(Qi = mean(q3), sdi = sd(q3))
    res <- list(Q2 = res_q2, Q3 = res_q3)
    res <- bind_rows(res,.id = "dist")
    return(res)
  }
  res <- lapply(res,agg.func)
  return(res)
}


comp_function <- function(idx,machine,part_dat){
  
  part_dat <- part_dat[c(machine$input,machine$output)]
  
  parts <-map2(idx,part_dat,~apply(.x,MARGIN = 2,function(x).y[x,]))
  
  parts_p <- parts[machine$output]
  parts_p <- transpose(parts_p)
  parts$recal_feed = map(parts_p,bind_rows)
  
  if(length(machine$input) > 1){
    parts_f <- parts[machine$input]
    parts_f <- transpose(parts_f)
    parts$feed = map(parts_f,bind_rows)
    parts[machine$input] <- NULL
  }else{
    parts$feed = parts[[machine$input]]
    parts[machine$input] <- NULL
  }

  comp <- function(x){
    res <- colSums(x)/sum(x)*100
    return(res)
  }
  res <- lapply(parts,function(x)t(sapply(x,comp)))
  
  agg.func <- function(res){
    res <- gather(as.data.frame(res),key = "Mineral",value = "value") %>% 
    group_by(Mineral) %>% 
    summarise(mean_comp = mean(value), sd_comp = sd(value))}
  res <- lapply(res,agg.func)
  return(res)
}



#calculating mean Volume based on ECD:




cal_mean_d <- function(x){
  volume <- pi/6*x$ECD^3
  mean_d <- sum(volume*x$ECD)/sum(volume)
  mean_dens <- sum(volume*x$dens)/sum(volume)
  mean_mass <- 1/6*pi*mean_d^3 * mean_dens * 10^(-12)
  res <- data.frame(ecd = mean_d,
                    dens = mean_dens,
                    mass = mean_mass)
  return(res)
}


cal_mean_d <- function(x){
  volume <- pi/6*x$ECD^3
  mass <- volume * x$dens * 10^(-12)
  mean_mass <- mean(mass)
  res <- data.frame(mass = mean_mass)
  return(res)
}

