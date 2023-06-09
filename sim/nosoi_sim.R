## Initialize stable working environment and store time of intiation
rm(list=ls())

#setwd("/Users/macbook/Dropbox (UFL)/DYNAMITE/HIVdynamite/nosoi_simulations")
# List of packages for session
.packages <-  c("phytools", "ape", "parallel", "ggplot2", "viridis", "igraph", "ggnetwork", 
                "ggpubr", "ggtree", "treeio", "ape", "remotes", "dplyr", "plyr", "phyclust", "phangorn") 
.github_packages <- c("slequime/nosoi", "emillykkejensen/familyR") 

# Install CRAN packages (if not already installed)
#.inst <- .packages %in% installed.packages()
#if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
#.inst_github <- .github_packages %in% installed.packages()
## Install GitHub packages(if not already installed)
#if(length(.github_packages[!.inst_github]) > 0) try(remotes::install_github(.github_packages[!.inst_github]))
#if(length(.github_packages[!.inst_github]) > 0) try(devtools::install_github(.github_packages[!.inst_github]))

# Load packages into session 
lapply(.packages, require, character.only=TRUE)
lapply(gsub(".+\\/(.+)", "\\1", .github_packages), require, character.only=TRUE)

### Set seeds #############################################################################
# get simulation parameters

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0) {
  sim_index=1
  seeds = 111
} else{
  sim_index = as.numeric(args[1]) # arg 1 fraction new
  seeds = readLines('seeds.txt')
}
#sim_index="test"
set.seed(seeds[sim_index])
numCores = detectCores()
`%notin%` <- Negate(`%in%`)



## Matrix generation #######################################################################
traits <- data.frame(location=rbind('A', 'B', 'C', 'D', 'E', 'F', 'G'))
Q <-list()
for (column in 1:ncol(traits)) {
  suppressWarnings({ ## We know it fills diagonals with NAs
    Q[[column]] <- diag(unique(traits[,column]), nrow = length(unique(traits[,column])))
  })
  #    diag(Q[[column]]) = 1-nrow(Q[[column]])
  diag(Q[[column]]) = 0
  non.diag <- 1/(nrow(Q[[column]])-1)
  Q[[column]][lower.tri(Q[[column]])] <- non.diag
  Q[[column]][upper.tri(Q[[column]])] <- non.diag
  colnames(Q[[column]]) <- rownames(Q[[column]]) <- unique(traits[,column])
  ## Modify Q matrix to allow for birth of cluster I from cluster H
  #  Q[[column]][nrow(traits)-1,] <- c(rep(0.0, nrow(traits)-1),1.0) # H only gives rise to I
  #  Q[[column]][1,] <- c(0.0,rep(1/(nrow(traits)-2),nrow(traits)-2),0.0) # A cannot give rise to I; the remaining clusters can stay at 0.125 because they have a 0 probability of leaving (see below)
}

Q <- plyr::compact(Q)

## Increase relative probability of growing and shrinking clusters
Q[[1]][1,] = c(0, rep(0.1, 4), rep(0.3,2))


# OR
#create matrix layout
# df.matrix = df %>% spread(To, value=N, fill= 0)
# df.matrix = as.matrix(df.matrix[-1])
# rownames(df.matrix) = colnames(df.matrix)
# 
# #Get probabilities (beware, rows should sum up to 1)
# df_transpose = t(df.matrix)
# probabilities <- apply(df_transpose, 1, function(i) i/sum(i))
# transition.matrix = t(probabilities)



# #pExit daily probability for a host to leave the simulation (either cured, died, etc.).
p_Exit_fct  <- function(t, t_sero){ 
  if(t <= t_sero){p=0}
  if(t > t_sero){p=0.95} 
  return(p)
}
t_sero_fct <- function(x){rnorm(x,14,2)} #Approximately 14 days

#pMove probability (per unit of time) for a host to do move, i.e. to leave its current state (for example, leaving state “A”). 
#It should not be confused with the probabilities extracted from the structure.matrix, which represent the probability to go 
#to a specific location once a movement is ongoing (for example, going to “B” or “C” while coming from “A”).
p_Move_fct  <- function(t, current.in, host.count){
  
  if(t>=1) {
    if(current.in=="A" & host.count > 2){return(0.0018)}
    if(current.in=="A" & host.count <= 2){return(0)} 
    if(current.in=="B"){return(0)}
    if(current.in=="C"){return(0)}
    if(current.in=="D"){return(0)}
    if(current.in=="E"){return(0)}
    if(current.in=="F"){return(0)}
    if(current.in=="G"){return(0)}
  } else{
    return(0)
  }
}
#  t_clust_fct <- function(x){rnorm(x,mean = 3,sd=1)}

n_contact_fct = function(t, current.in, host.count){
  
  
  n0_low <- abs(round(rnorm(1, 4, 1)))
  assign("n0_low", n0_low, envir = globalenv())
  
  growth <- list()
  r1 <- abs(rnorm(1, 0.0018, 0.0009))
  r2 <- abs(rnorm(1, 0.018, 0.009))
  # max_r1 = length(growth$log_growth[which(growth$log_growth<=20)])
  # max_r2 = length(growth$lin_growth[which(growth$lin_growth<=20)])

  assign("r1", r1, envir = globalenv())
  assign("r2", r2, envir = globalenv())
  
  growth$log_growth <- round(20*n0_low/(n0_low+(20-n0_low)*exp(-r1*host.count)))
  growth$lin_growth = r2*host.count+n0_low
  growth$poly_growth = 0.1*t^3 - 0.5 * t^2 - t + 5 + rnorm(length(t),0,10) #polyomial regression curve

  growth <- lapply(growth, function(x) {
    x[x<0] = 0
    return(x)
  })
  
  n0_hi <- round(rnorm(1, 6, 1))
  assign("n0_hi", n0_hi, envir = globalenv())
  
  decay <- list()
  r3 <- abs(rnorm(1, 0.018, 0.009))
  r4 <- abs(rnorm(1, 0.18, 0.09))
  assign("r3", r3, envir = globalenv())
  assign("r4", r4, envir = globalenv())
  
  
  decay$log_decay <- round(10*n0_hi/(n0_hi+(10-n0_hi)*exp(r3*host.count)))
  decay$lin_decay = -r4*host.count+n0_hi
  decay$poly_decay = n0_hi*sin(t^-40)# See sim_parametric.R!
  decay <- lapply(decay, function(x) {
    x[x<0] = 0
    return(x)
  })
  
  
  if(current.in=="A") {p=abs(round(rnorm(1, 16, 1), 0))}
  if(current.in=="B") {p=abs(round(rnorm(1, 4, 1), 0))}
  if(current.in=="C") {p=abs(round(rnorm(1, 4, 1), 0))}
  if(current.in=="D") {p=abs(round(rnorm(1, 6, 1), 0))}
  if(current.in=="E") {p=abs(round(rnorm(1, 4, 1), 0))}
  # if(current.in=="F" & host.count <= max_r2 ) {p=growth$lin_growth}
  # if(current.in=="F" & host.count > max_r2) {p=0}
  if(current.in=="F") {p=growth$lin_growth}
  if(current.in=="G") {p=decay$lin_decay}
  return(p)
}

#pTrans represents the probability of transmission over time (when a contact occurs).
# in the form of a threshold function: before a certain amount of time since initial infection, the host does not transmit (incubation time, which we call t_incub), and after that time, it will transmit with a certain (constant) probability (which we call p_max). This function is dependent of the time since the host’s infection t.
p_Trans_fct <- function(t, current.in, t_incub){
  if(t < t_incub){p=0}
  if(t >= t_incub & current.in=="A"){p=0.015}
  if(t >= t_incub & current.in=="B"){p=0.09}
  if(t >= t_incub & current.in=="C"){p=0.11}
  if(t >= t_incub & current.in=="D"){p=0.09}
  if(t >= t_incub & current.in=="E"){p=0.12}
  if(t >= t_incub & current.in=="F"){p=0.11}
  if(t >= t_incub & current.in=="G"){p=0.11}
  return(p)
}

t_incub_fct <- function(x){rnorm(x,mean = 5,sd=2)} #Approximately 4.2 days
#p_max_fct <- function(x){rbeta(x,shape1 = 1,shape2=8)}# Mean of roughly 0.10

# Starting the simulation ------------------------------------

#set.seed(111)

num_sims = 100
i = 1
while(i<=num_sims){
  SimulationSingle <- nosoiSim(type="single", # Number of hosts
                               popStructure="discrete", #discrete or continuous
                               structure.matrix = Q[[1]], # prob matrix defined above (row sums to 1, with diags as zero)
                               length.sim = 365, # Max number of time units (can be days, months, weeks, etc.)
                               max.infected = 10000, #maximum number of individuals that can be infected during the simulation.
                               init.individuals = 1, #number of individuals (an integer above 1) that will start a transmission chain. Keep in mind that you will have as many transmission chains as initial individuals, which is equivalent as launching a number of independent nosoi simulations.
                               init.structure = "A",
                               
                               pExit = p_Exit_fct,
                               param.pExit=list(t_sero=t_sero_fct),
                               timeDep.pExit=FALSE,
                               diff.pExit=FALSE,
                               
                               pMove = p_Move_fct,
                               hostCount.pMove=TRUE,
                               param.pMove=NA,
                               timeDep.pMove=FALSE,
                               diff.pMove=TRUE,
                               
                               nContact=n_contact_fct,
                               hostCount.nContact=TRUE,
                               param.nContact = NA,
                               timeDep.nContact=FALSE,
                               diff.nContact=TRUE,
                               
                               pTrans = p_Trans_fct,
                               hostCount.pTrans=FALSE,
                               param.pTrans = list(t_incub=t_incub_fct),
                               timeDep.pTrans=FALSE,
                               diff.pTrans=TRUE,
                               
                               prefix.host="S",
                               print.progress=TRUE,
                               print.step=100)
  
  if (SimulationSingle[["host.info.A"]]$N.infected > 50 ) {
    sum_sim <- summary(SimulationSingle)
    cumulative.table <- getCumulative(SimulationSingle)
    dynamics.table <- getDynamic(SimulationSingle)
    
    # cum.p <- ggplot(data=cumulative.table, aes(x=t, y=Count)) + geom_line() + theme_minimal() +
    #   labs(x="Time (t)",y="Cumulative count of infected hosts") + scale_y_log10()
    # cum.p.c <- ggplot(data=dynamics.table, aes(x=t, y=Count, color=state)) + geom_line() + theme_minimal() +
    #   labs(x="Time (t)",y="Number of active infected hosts") + scale_y_log10()
    # ggpubr::ggarrange(cum.p, cum.p.c, widths = 2, heights = 1, legend="right")
    # # ggsave("simtest_max10000.png", plot=last_plot())
    
    
    ## Grab tree #########################################################################################
    save.data <- function(){
      sim.tree <- getTransmissionTree(SimulationSingle)
      assign("sim.tree", sim.tree, envir = globalenv())
      # ggtree(test.nosoiA.tree, color = "gray30") + geom_nodepoint(aes(color=state)) + geom_tippoint(aes(color=state)) + 
      #   theme_tree2() + xlab("Time (t)") + theme(legend.position = c(0,0.8), 
      #                                            legend.title = element_blank(),
      #                                            legend.key = element_blank()) 
      
      # Get sampled tree
      getSample <- function(SimulationSingle) {
        int.nodes <- sample((Ntip(sim.tree@phylo)+2):(Ntip(sim.tree@phylo)+sim.tree@phylo$Nnode)) #randomize order of internal nodes
        n <- unique(sum_sim$dynamics$state[sum_sim$dynamics$state != "A"]) 
        s <- seq(20,200,1)
        min_s=min(s)
        
        state.list <- list()
        for (i in 1:length(n)) {
          s_n=sim.tree@data$node[sim.tree@data$state==n[i]]
          if(length(s_n)>=min_s) {
            nodes <- sample(s_n, replace=F)
            state <- n[i]
            state.list[[i]] <- data.frame(nodes=nodes, state=state)
            state.list[[i]] <- dplyr::filter(state.list[[i]],
                                             nodes %in% int.nodes)
          } else {
            state.list[[i]]=NULL
          }
          
        } 
        
        state.list = compact(state.list)
        
        #End state.list generation      
        
        sampleState <- function(state) {
          max_length <- length(state$nodes)
          tcs <- data.frame(taxa=NA, state=NA)
          i=1
          while (isTRUE(is.na(tcs$taxa) & i <= max_length)) {
            n=state$nodes[i]    
            tcs.phylo <- extract.clade(sim.tree@phylo, n)
            
            ## Above might sample only a part of the clade, so need full clade.
            repeat {
              root_tcs=findMRCA(sim.tree@phylo, tcs.phylo$tip.label)
              state_root=sim.tree@data$state[sim.tree@data$node==root_tcs]
              parent=Ancestors(sim.tree@phylo, root_tcs, type="parent")
              state_parent = sim.tree@data$state[sim.tree@data$node==parent]
              if(isTRUE(state_parent==state_root)) {
                tcs.phylo <- extract.clade(sim.tree@phylo, parent)
              } else {
                break
              }
            }
            
            tcs.taxa <- subset(sim.tree@data, sim.tree@data$host %in% tcs.phylo$tip.label)
            unique_hosts <- unique(tcs.taxa$host)
            
            lo <- 0.05 ; hi <- 1
            r <- 0.002
            C <- exp(-r*hi); D <- exp(-r*lo)
            n <- 100
            U <- runif(n,min=C,max=D)
            X <- (1/r)*log(1/U)
            #hist(X,breaks=10,xlim=c(0,1))
            sf <- sample(X, 1)
            
            sampled_hosts <- sample(unique_hosts, round(sf*length(unique_hosts)), replace=F)
            
            sampled_states <- dplyr::select(SimulationSingle$host.info.A$table.hosts, hosts.ID, current.in) %>%
              filter(hosts.ID %in% sampled_hosts) %>%
              distinct() %>%
              dplyr::select(current.in)
            sampled_states <- as.vector(sampled_states$current.in)
            
            # Remember that background individuals infected by high-risk individuals are still considered part of the cluster, so need the following:        
            sampled_state_origins <- dplyr::select(SimulationSingle$host.info.A$table.hosts, hosts.ID, inf.in) %>%
              filter(hosts.ID %in% sampled_hosts) %>%
              distinct() %>%
              dplyr::select(inf.in)
            sampled_state_origins <- as.vector(sampled_state_origins$inf.in)
            sampled_state_origins[1]=state$state[1]
            
            
            if(isTRUE(length(grep(state$state[1], tcs.taxa$state)) >= 0.51*length(tcs.taxa$state) &
                      length(sampled_hosts) %in% s)) {
              tcs <- data.frame(taxa=sampled_hosts, 
                                state=sampled_states, 
                                sampling_fraction=sf, 
                                cluster_size=length(unique_hosts),
                                #                 fraction_state=length(sampled_states[sampled_states==state$state[1]])/length(sampled_states)) 
                                fraction_state=length(sampled_state_origins[sampled_state_origins==state$state[1]])/length(sampled_state_origins))
              i=i
            } else {
              tcs <- data.frame(taxa=NA, state=NA, sampling_fraction=NA, cluster_size=NA, perc_cluster=NA)
              i=i+1}
          } # End while loop
          return(tcs)
        } # End function
        
        print("Sampling high-risk groups")
        tcs.list <- mclapply(state.list, sampleState, mc.cores=numCores)
        tcs.list <- mclapply(tcs.list, function(x) { # Removes small clusters that have not passed the test and resulted in NA taxa
          if(nrow(x)==1 & is.na(x$taxa[1])) {
            x=NULL
          } else{
            x=x
          }
          return(x)
        }, mc.cores=numCores)
        tcs.list=compact(tcs.list)
        
        #         s_rand <- round(rnorm(10,50,10))
        #         group_A <- sample(sim.tree@data$node[sim.tree@data$state=="A"])
        #         group_A <- group_A[group_A %in% int.nodes]
        #         
        #         sampleA <- function(group_A) {
        #           cluster_A <- NULL
        #           while(is.null(cluster_A)){
        #             true.cluster.phylo <- extract.clade(sim.tree@phylo, sample(group_A, 1))
        #             true.cluster.taxa <- unique(subset(sim.tree@data, sim.tree@data$host %in% true.cluster.phylo$tip.label))
        #             unique_hosts <- unique(true.cluster.taxa$host)
        #             
        #             #       lo <- 0.05 ; hi <- 1
        #             #        r <- 0.002
        #             #        C <- exp(-r*hi); D <- exp(-r*lo)
        #             #        n <- 100
        #             #        U <- runif(n,min=C,max=D)
        #             #        X <- (1/r)*log(1/U)
        #             #        #hist(X,breaks=10,xlim=c(0,1))
        #             #        sf <- sample(X, 1)
        #             
        #             #        sampled_hosts <- sample(unique_hosts, round(sf*length(unique_hosts)), replace=F)
        #             sf=1
        #             sampled_hosts=unique_hosts
        #             
        #             if(isTRUE(length(sampled_hosts) %in% s_rand &
        #                       length(sampled_hosts) >= 5 &
        #                       length(grep("A", true.cluster.taxa$state)) >= 0.95*length(true.cluster.taxa$state))) {
        #               cluster_A <- data.frame(taxa=sampled_hosts, state="A", sampling_fraction=sf, cluster_size=length(unique_hosts)) 
        #             }
        #           }
        #           return(cluster_A)
        #         }
        #         
        #         print("Extracting complete transmission chain from background population (A)")
        #         group_A_cluster <- sampleA(group_A)
        #         tcs.list <- append(list(group_A_cluster), tcs.list) %>%
        #           Filter(function(a) any(!is.na(a)), .)
        
        print("Sampling background population (A)") 
        table.hosts <- getTableHosts(SimulationSingle, pop="A")
        sf = 0.16*(SimulationSingle$host.info.A$N.infected - length(unlist(lapply(tcs.list, function(x) x$taxa))))
        sampled.hosts <- sample(table.hosts$hosts.ID[table.hosts$current.in=="A"], sf, replace=F)
        # Determne which of samlped hosts not in tcs.list in order to calculate sampling fraction of the background population
        clustered <- unname(unlist(lapply(tcs.list, "[", 'taxa')))
        bg_sample <- sampled.hosts[sampled.hosts %notin% clustered]
        bg_size <- length(bg_sample)
        assign("bg_size", bg_size, envir = globalenv())
        bg_sampling_frac <- bg_size/length(table.hosts$hosts.ID)
        assign("bg_sampling_frac", bg_sampling_frac, envir = globalenv())
        ## Add these individuals to list of randomly sampled individuals
        sampled.hosts <- c(bg_sample, clustered)
        ## Extract tree for list of individuals from the full simulation tree
        sampled.tree <- sampleTransmissionTreeFromExiting(sim.tree, sampled.hosts)        
        #sampled.tree <- keep.tip(sim.tree@phylo, sampled.hosts)
        
        
        
        list_nodes <- mclapply(tcs.list, function(x) {
          n <- findMRCA(sampled.tree@phylo, x$taxa)
          x2 <- as_tibble(extract.clade(sampled.tree@phylo, n))
          return(x2)
        }, mc.cores=numCores)
        
        # Give names to clusters (c1..n)
        for (i in seq_along(list_nodes)) {
          names(list_nodes)[i] <- paste0("c", i)
          names(tcs.list)[i] <- paste0("c", i)
        }
        
        assign("list_nodes", list_nodes, envir = globalenv()) # Remember now a tibble
        
        trans_clusters <- do.call(rbind, tcs.list)
        ## Create background for export 
        
        sampled_hosts <- sampled.hosts[sampled.hosts %notin% trans_clusters$taxa]
        
        background <- list(data.frame(taxa=sampled_hosts, 
                                      state="A", 
                                      sampling_fraction=bg_sampling_frac) %>%
                             mutate(cluster_size=nrow(.),
                                    fraction_state = 1))
        names(background)="Background"
        
        tcs.list <- append(background, tcs.list) %>%
          Filter(function(a) any(!is.na(a)), .)
        assign("tcs.list", tcs.list, envir = globalenv()) # Remember now a tibble
        
        return(as_tibble(sampled.tree))
      }
      sampled.tree <- getSample(SimulationSingle)      
      
      # Transform sampled tree into a tbl object and assign cluster IDs to internal and external nodes found in list of clusters
      print("Annotating tree according to cluster/group")
      sampled.tree$cluster_id <- "Background"
      sim.tree@data$cluster_id <- "Background"
      for (i in seq_along(list_nodes)) {
        for (j in 1:length(sampled.tree$label)) {
          if (isTRUE(sampled.tree$label[j] %in% list_nodes[[i]]$label)) {
            sampled.tree$cluster_id[j] = names(list_nodes)[i]
          } # End if-else statement
        } #End for loop
        for (k in 1:length(sim.tree@data$host)) {
          if (isTRUE(sim.tree@data$host[k] %in% list_nodes[[i]]$label)) {
            sim.tree@data$cluster_id[k] = names(list_nodes)[i]
          } # End if-else statement
        } # End for loop
      } # End for loop
      
      # If background individual found in cluster and is a part of a new clade of background infections, 
      # keep the originating individual as clustered, but assign "Background" to all descendants
      for (i in 1:nrow(sampled.tree)) {
        if(isTRUE(sampled.tree$cluster_id[i] %notin% "Background")) {
          #     if(isTRUE(sampled.tree$state[i]=="A" & sampled.tree$state[sampled.tree$node==sampled.tree$parent[i]]=="A")) {
          if(isTRUE(sampled.tree$state[i]=="A" & 
                    SimulationSingle$host.info.A$table.hosts$inf.in[SimulationSingle$host.info.A$table.hosts$hosts.ID==sampled.tree$label[i]]=="A")) {
            sampled.tree$cluster_id[i]="Background"
          }
        }
      }
      
      
      print("Assigning arbitrary dates to sequences and including in sequence names")
      sampled.tree$date <- as.Date(sampled.tree$time, origin="2019-12-31")
      sampled.tree$label <- paste(sampled.tree$label, sampled.tree$state, sampled.tree$date, sep="|")
      sampled.tree <- dplyr::select(sampled.tree, parent, node, branch.length, label, cluster_id, state, date) %>%
        as_tibble()
      class(sampled.tree) = c("tbl_tree", class(sampled.tree))
      
      t2 <- as.treedata(sampled.tree)
      write.beast(t2, paste0('sim_', sim_index, "_sampled.tree")) # NEED THIS OUTPUT
      
      print("Saving sequence assignment information in csv")
      
      write.csv(sampled.tree, paste0("sim_", sim_index, "_data.csv"), quote=F, row.names=F)
      
      true.cluster.dyn <- function(conf.level=0.95){
        cluster_dynamics <- data.frame(sim=sim_index,
                                       state=rbind('A', 'B', 'C', 'D', 'E', 'F', 'G'),
                                       dynamic = rbind('static', 'static', 'static', 'static', 'static', 
                                                       'growth', 'decay'),
                                       dynamic_rparam = rbind(0, 0, 0, 0, 0, r2, r4) )
        
        cluster_dynamics$cluster_id <- NA
        cluster_dynamics$v_sampling_fraction <- NA
        cluster_dynamics$v_cluster_size <- NA
        for (i in seq_along(tcs.list)) {
          for (j in 1:nrow(cluster_dynamics)) {
            if(isTRUE(cluster_dynamics$state[j] == names(which.max(table(tcs.list[[i]]$state))))) {
              cluster_dynamics$cluster_id[j] <- names(tcs.list)[i]
              cluster_dynamics$v_sampling_fraction[j] <- tcs.list[[i]]$sampling_fraction[1]
              cluster_dynamics$v_cluster_size[j] <- tcs.list[[i]]$cluster_size[1]
            } 
          }
        }
        cluster_dynamics <- cluster_dynamics %>%
          filter(!is.na(cluster_id))
        
        # host.data = dplyr::select(SimulationSingle$host.info.A$table.hosts,
        #                           hosts.ID, inf.time) %>%
        #   dplyr::rename(host=hosts.ID) %>%
        #   dplyr::left_join(., dplyr::select(sim.tree@data, host, state, cluster_id), by="host") %>%
        #   dplyr::filter(cluster_id != "Background") %>%
        #   dplyr::distinct() %>%
        #   dplyr::group_by(cluster_id) %>%
        #   dplyr::group_split()
        # 
        # coef.fcn.decay = function(y, x) {
        #   if(max(x)<max(dynamics.table$t)-1) {
        #     x=c(x, max(dynamics.table$t))
        #     y=c(y,0)
        #   }
        #   m = which(y==max(y))[1]
        #   x=x[m:length(x)]
        #   y=y[m:length(y)]
        #   p=summary(lm(y~x))$coefficients
        #   if(nrow(p)>1) {
        #     result=data.frame(Intercept=p[1,1], 
        #                       slope=p[2,1],
        #                       p_val=p[2,4])
        #   } else {
        #     result = data.frame(Intercept=NA, slope=NA, p_val=NA)
        #   }
        #   return(result)
        # }
        # 
        # coef.fcn.growth = function(y, x) {
        #   df = data.frame(x=x, y=y)
        #             dplyr::mutate(ratio = y/lag(y)) %>%
        #              filter(ratio>=1)
        #   p=summary(lm(df$ratio~df$x.in))$coefficients
        #   if(nrow(p)>1) {
        #     result=data.frame(Intercept=p[1,1], 
        #                       slope=p[2,1],
        #                       p_val=p[2,4])
        #   } else {
        #     result = data.frame(Intercept=NA, slope=NA, p_val=NA)
        #   }
        #   return(result)
        # }
        
        # host_rates = suppressMessages(
        #   do.call(rbind, lapply(host.data, function(x) {
        #     result=data.frame(cluster_id=x$cluster_id[1],
        #                       state=x$state[1])
        #     df = dplyr::group_by(x, cluster_id, inf.time) %>%
        #                   dplyr::summarize(n = length(unique(host)))
        # 
        #     if(x$state[1]=="G") {
        #       eq = coef.fcn.decay(df$n, df$inf.time)}
        #     else {
        #       eq = coef.fcn.growth(df$n, df$inf.time)
        #     }
        #     result = cbind(result, eq)
        #     return(result)
        #   }) ))
        # 
        # cluster_dynamics = dplyr::left_join(cluster_dynamics, host_rates, by=c("state", "cluster_id")) %>%
        #   dplyr::mutate(dynamic = if_else(state == 'G' &
        #                                     is.na(Intercept) &
        #                                     ( slope > 0 | is.na(slope) ) & 
        #                                     ( p_val < 0.05 | is.na(p_val) ), 
        #                                   "static", dynamic))  %>%
        #   dplyr::mutate(dynamic = if_else(state == 'F' &
        #                                     p_val >= 0.10, 
        #                                   "static", dynamic)) 
        
        
        ## Estimate R0 for each cluster based on paramaters used in simulation
        # s <- 100 # sample size of 100
        # n_contacts_E <- rnorm(s, 4, 1)
        # 
        # t_incub <- rnorm(s, 5, 2)
        # t_exit <- rnorm(s, 14, 2)
        # p_trans_E <- 0.2
        # 
        # Z=qnorm(0.5*(1 + conf.level))
        # 
        # R0_E = sample(n_contacts_E*p_trans_E*(t_exit-t_incub),
        #               size=s, replace=F)
        # 
        # logR0_E = log(R0_E)
        # 
        # SElogR0_E = sd(logR0_E, na.rm=T)/sqrt(s) # standard deviation of 2, sample size of 10
        # lower_E = exp(mean(logR0_E, na.rm=T) - Z*SElogR0_E) 
        # upper_E = exp(mean(logR0_E, na.rm=T) + Z*SElogR0_E) 
        # 
        # #R0_A <- read.table(paste0("R0_", sim_index, ".tab"), header=T)
        # #cluster_dynamics$mean_R0[1] = as.numeric(R0_A[1]) # Can't use exact R0 because standard deviation includes the high R0 of E and 0.
        # cluster_dynamics$mean_R0[1] = sum_sim$R0$R0.mean
        # cluster_dynamics$upper_R0[1] = quantile(sum_sim$R0$R0.dist, 0.95)
        # cluster_dynamics$lower_R0[1] = quantile(sum_sim$R0$R0.dist, 0.05)
        # 
        # cluster_dynamics$mean_R0[5] = mean(R0_E)
        # cluster_dynamics$upper_R0[5] = upper_E
        # cluster_dynamics$lower_R0[5] = lower_E
        # 
        
        return(cluster_dynamics)
        
      }# End function
      cluster_dynamics <- true.cluster.dyn()
      
      print("Saving dynamics summary information in csv")
      write.csv(cluster_dynamics, paste0("sim_", sim_index, "_summary.csv"), 
                quote=F, row.names = F)
    }  
    
    
    ## If any cluster exceeds background population in number, do not output results; however, if not, proceed with tree extraction.
    max.t <- max(dynamics.table$t)
    if(isTRUE(length(dynamics.table$Count[dynamics.table$state == 'A' & 
                                          dynamics.table$t == max.t])==0) |
       any(dynamics.table$Count[dynamics.table$state %in% c('B', 'C', 'D', 'E', 'F', 'G') & 
                                dynamics.table$t == max.t] > 
           dynamics.table$Count[dynamics.table$state == 'A' & 
                                dynamics.table$t == max.t] )) {
      
      i=i+1
    } else {
      print(paste0("Simulating tree for seed number ", seeds[sim_index], " run number ", sim_index))
      save.data()
      i=1000} 
  } else {
    i=i+1
  }
}# End while loop


## End ###################################################################################
