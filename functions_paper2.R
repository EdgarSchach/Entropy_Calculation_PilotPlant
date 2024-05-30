
batch.entropy.corr <- function(part.dat,sample_point,corr = T){
  # Calculating Extensive Properties for Feed Batch:
  
  N000 <- sum(part.dat)
  N0j0 <- colSums(part.dat)
  N00k <- rowSums(part.dat)
  N0jk <- part.dat
  
  if(corr == T){
    con <- corrected_concentrations[sample_point,]/100
    c0j0 <- unlist(con)
  }else{
    c0j0 <- colSums(part.dat)/sum(part.dat)
  }
  c0jk <- N0jk/N00k # integrowth of batchstage
  p00k <- N00k/N000 # mass distribution in batchstage
  p0jk <- sweep(N0jk,2,N0j0,"/") #component distribution in batchstage
  
  #Entropy Calculation for Feed####
  #c0j0.p0jk
  
  ent_c0j0 <- .entropy(c0j0)  #*con_weighting #- log(length(c0j0))
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


entropyforboots.corrected <- function(feed,products,weights,recal = T, con.corr = T,idx = idx){
  #CALCULATIONS FOR FEEDSTAGE
  
  # extensive properties for feed stage:
  
  Fi00 <- lapply(feed,sum)
  Fi0k <- lapply(feed,rowSums)
  Fij0 <- lapply(feed,colSums) # this might be corrected later on
  Fijk <- feed
  
  # calculation of intensive properties:
  if (con.corr == F){
    cij0 <- map(Fij0,~.x/sum(.x)) 
  }else{
    if(recal == T){
      cij0 <- t(grouped_concentrations[names(products),])
      cij0 <- cij0 %*% weights[names(products)]/100
      cij0 <- list(recal_feed = as.data.frame(t(cij0))) 
    }else{
      cij0 <- map(names(feed),~grouped_concentrations[.x,]/100)
    }
    con_weighting <- map(cij0,~1/(1-.x$Rest)) 
    cij0 <- map(cij0,~select(.x,-Rest))
    cij0 <- map(cij0,unlist) 
    
  }
  
  cijk <- map2(Fijk,Fi0k,~sweep(.x,1,.y,"/")) #integrowth of product stage
  
  if (recal == T){
    ti00 <-  1
  }else{
    vm <- vol_masses[names(feed)]
    vm <- map2(vm,idx[names(feed)],~.x[.y,])
    vm <- map(vm,~sum(.x$mass_vol))
    ti00 <- clo(unlist(vm)) 
  }
  
  ti001 <- weights[names(feed)]
  
  pi0k <- map(Fi0k,~.x/sum(.x))  # mass distribution in product stage
  pijk <- map2(Fijk, Fij0, ~ sweep(.x,2,.y,"/")) #component distribution in product stage
  pijk <- map(pijk,~mutate_all(.x,~replace(., is.na(.), 0)))
  # pijk <- map(pijk,~select(.x,-Rest))
  
  ent_ti00 <- .entropy(ti00)
  if(con.corr == T){
    ent_cij0 <- sapply(cij0,.entropy) *  unlist(con_weighting)
    ent_cij0 <- ti00 %*% ent_cij0
    pijk <- map(pijk,~select(.x,-Rest))
  }else{
    ent_cij0 <- ti00 %*% sapply(cij0,.entropy) 
  }
  
  ent_pijk <- lapply(pijk,function(x)apply(x,2,.entropy))
  ent_pijk <- map2(cij0,ent_pijk,~.x%*%.y,)
  ent_pijk <- ti00 %*% unlist(ent_pijk)
  
  ti00.cij0.pijk <- c(ti00 = ent_ti00, cij0 = ent_cij0, pijk = ent_pijk)
  
  #ti00.pi0k.cijk
  
  ent_ti00 <- .entropy(ti00)
  ent_pi0k <- sapply(pi0k,.entropy)
  ent_pi0k <- ti00 %*% ent_pi0k
  ent_cijk <- lapply(cijk,function(x)apply(x,1,.entropy))
  ent_cijk <- map2(pi0k,ent_cijk,~.x%*%.y)
  ent_cijk <- ti00 %*% unlist(ent_cijk)
  
  ti00.pi0k.cijk <- c(ti00 = ent_ti00,pi0k = ent_pi0k,cijk = ent_cijk)
  
  res_feedstage <- list(ti00.cij0.pijk = ti00.cij0.pijk,
                        ti00.pi0k.cijk = ti00.pi0k.cijk)
  res_feedstage <- bind_rows(res_feedstage,.id = "aggregation") %>% 
    gather(key = "entropy_type",value = "entropy",-aggregation) %>% na.omit()
  
  # CALCULATION FOR PRODUCT STAGE####
  
  Pi00 <- lapply(products,sum)
  Pi0k <- lapply(products,rowSums)
  Pij0 <- lapply(products,colSums) 
  Pijk <- products
  
  if (con.corr == F){
    cij0 <- map(Pij0,~.x/sum(.x)) 
  }else{
    cij0 <- map(names(products),~grouped_concentrations[.x,]/100)
    con_weighting <- map(cij0,~1/(1-.x$Rest)) 
    cij0 <- map(cij0,~select(.x,-Rest))
    cij0 <- map(cij0,unlist) 
  }
  
  cijk <- map2(Pijk,Pi0k,~sweep(.x,1,.y,"/")) #integrowth of product stage
  
  vm <- vol_masses[names(products)]
  vm <- map2(vm,idx[names(products)],~.x[.y,])
  vm <- map(vm,~sum(.x$mass_vol))
  
  ti00 <- clo(unlist(vm))
  ti001 <- weights[names(products)]
  
  pi0k <- map(Pi0k,~.x/sum(.x))  # mass distribution in product stage
  pijk <- map2(Pijk, Pij0, ~ sweep(.x,2,.y,"/")) #component distribution in product stage
  pijk <- map(pijk,~mutate_all(.x,~replace(., is.na(.), 0)))
  
  
  #ti00.cij0.pijk
  
  ent_ti00 <- .entropy(ti00)
  if(con.corr == T){
    ent_cij0 <- sapply(cij0,.entropy) *  unlist(con_weighting)
    ent_cij0 <- ti00 %*% ent_cij0
    pijk <- map(pijk,~select(.x,-Rest))
  }else{
    ent_cij0 <- ti00 %*% sapply(cij0,.entropy)
  }
  
  ent_pijk <- lapply(pijk,function(x)apply(x,2,.entropy))
  ent_pijk <- map2(cij0,ent_pijk,~.x%*%.y,)
  ent_pijk <- ti00 %*% unlist(ent_pijk)
  
  ti00.cij0.pijk <- c(ti00 = ent_ti00, cij0 = ent_cij0, pijk = ent_pijk)
  
  #ti00.pi0k.cijk
  
  ent_ti00 <- .entropy(ti00)
  ent_pi0k <- sapply(pi0k,.entropy)
  ent_pi0k <- ti00 %*% ent_pi0k
  ent_cijk <- lapply(cijk,function(x)apply(x,1,.entropy))
  ent_cijk <- map2(pi0k,ent_cijk,~.x%*%.y)
  ent_cijk <- ti00 %*% unlist(ent_cijk)
  ti00.pi0k.cijk <- c(ti00 = ent_ti00,pi0k = ent_pi0k,cijk = ent_cijk)
  
  res_productstage <- list(ti00.cij0.pijk = ti00.cij0.pijk,
                           ti00.pi0k.cijk = ti00.pi0k.cijk)
  
  res_productstage <- bind_rows(res_productstage,.id = "aggregation") %>% 
    gather(key = "entropy_type",value = "entropy",-aggregation) %>% na.omit()
  
  res <-  bind_rows(feed= res_feedstage,product = res_productstage,.id = "stage")
  
  
  res <- list(entropy = res,
              t_val = ti00)
  return(res)
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

