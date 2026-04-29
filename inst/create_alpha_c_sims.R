t_duration <- 700
t_steadied <- 100
t_blended <- 200

species_vec <- species_params(params)$species

p <- expand.grid(
    species = species_vec,
    c = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10),
    alpha = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
)

p$status_Y <- NA
p$BH_Y <- NA
p$status_B <- NA
p$BH_B <- NA
p$status_SSB <- NA
p$BH_SSB <- NA
p$status_N <- NA
p$BH_N <- NA

#Yield_list <-list()
#sim_list <- list()

for(c_val in unique(p$c)){

    for(alpha_val in unique(p$alpha)){

        my_blender <- make_blended_ssBH_FMort(
            t_max_blend = t_blended,
            target_c = c_val,
            t_steady = t_steadied,
            alpha_max = alpha_val
        )

        sim_BH_start <- setRateFunction(ps, "FMort", "my_blender")
        sim_BH <- project(sim_BH_start, t_max = t_duration, effort = 1)

        sim_list[[paste0("c_",c_val,"_a_",alpha_val)]] <- sim_BH

        Yield <- getYield(sim_BH)
        Yield_list[[paste0("c_",c_val,"_a_",alpha_val)]] <- Yield

        Biomass <- getBiomass(sim_BH)
        SSB <- getSSB(sim_BH)
        N <- getN(sim_BH)

        for(species in unique(species_params(ps)$species)){
            rows <- p$species == species & p$c == c_val & p$alpha == alpha_val
            p$status_Y[rows] <- Yield[99, species]
            p$BH_Y[rows] <- Yield[699, species]
            p$status_B[rows] <- Biomass[99, species]
            p$BH_B[rows] <- Biomass[699, species]

            p$status_SSB[rows] <- SSB[99, species]
            p$BH_SSB[rows] <- SSB[699, species]

            p$status_N[rows] <- N[99, species]
            p$BH_N[rows] <- N[699, species]
        }
    }

}

#save simlist
##saveRDS(sim_list, "/Users/jessicawestworth/Desktop/BH simulations/sim_list_extended/sim_list.rds")
#save p
##saveRDS(p, "/Users/jessicawestworth/Desktop/BH simulations/p_extended_time/p_extended_time.rds")
#save Yield
##saveRDS(Yield_list, "/Users/jessicawestworth/Desktop/BH simulations/yield_extended/Yield_extended_time.rds")


p<-readRDS("/Users/jessicawestworth/Desktop/BH simulations/p_extended_time/p_extended_time.rds")
#Full BH plots (where: alpha=1)
p_full<-p%>%
    filter(alpha==1)

ggplot(p_full, aes(x=c,y=log(BH_Y), group=species))+
    geom_line(aes(color=species))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(y= "Log Yield (g/year/m^2)", x= "c", color="Species")

ggplot(p_full, aes(x=c,y=log(BH_B), group=species))+
    geom_line(aes(color=species))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(y= "Log Biomass (g/m^2)", x= "c", color="Species")

ggplot(p_full, aes(x=c,y=log(BH_SSB), group=species))+
    geom_line(aes(color=species))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(y= "Log Spawning Stock Biomass (g/m^2)", x= "c", color="Species")

ggplot(p_full, aes(x=c,y=log(BH_N), group=species))+
    geom_line(aes(color=species))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(y= "Log Number of Individuals (per m^2)", x= "c", color="Species")

p_grouped<-p%>%
    group_by(c,alpha)%>%
    summarise(status_Y=sum(status_Y),
              BH_Y=sum(BH_Y),
              status_B=sum(status_B),
              BH_B=sum(BH_B),
              status_SSB=sum(status_SSB),
              BH_SSB=sum(BH_SSB),
              status_N=sum(status_N),
              BH_N=sum(BH_N)
    )

#LFI Catches: fish larger than 50 cm

for(c_val in unique(p_grouped$c)){
    for(alpha_val in unique(p_grouped$alpha)){

        sim<-sim_list[[paste0("c_",c_val,"_a_",alpha_val)]]
        my_blender <- make_blended_ssBH_FMort(
            t_max_blend = t_blended,
            target_c = c_val,
            t_steady = t_steadied,
            alpha_max = alpha_val
        )

        sim_BH_start <- setRateFunction(ps, "FMort", "my_blender")

        sim<-sim_list[[paste0("c_",c_val,"_a_",alpha_val)]]
        f<-getFMort(sim, drop = FALSE)
        n<-sim@n

        biomass <- sweep(sim@n, 3, sim@params@w * sim@params@dw, "*")
        yield<-f * biomass
        catch<-f*n
        number<-n

        status_LFY_s<-0
        status_Y_s<-0
        status_LFY<-0
        status_Y<-0

        BH_LFY_s<-0
        BH_Y_s<-0
        BH_LFY<-0
        BH_Y<-0

        status_LFB_s<-0
        status_B_s<-0
        status_LFB<-0
        status_B<-0

        BH_LFB_s<-0
        BH_B_s<-0
        BH_LFB<-0
        BH_B<-0

        for(s in 1:nrow(sim@params@species_params)){
            l<-50
            a<-sim@params@species_params$a[s]
            b<-sim@params@species_params$b[s]
            w_cut<-a*(l^b)

            #yields
            L_yield<-yield[,s,w>w_cut]

            status_LFY_s<-sum(L_yield[99,])
            status_Y_s<-sum(yield[99,s,])
            status_LFY<-status_LFY+status_LFY_s
            status_Y<-status_Y+status_Y_s

            BH_LFY_s<-sum(L_yield[699,])
            BH_Y_s<-sum(yield[699,s,])
            BH_LFY<-BH_LFY+BH_LFY_s
            BH_Y<-BH_Y+BH_Y_s

            #Biomass
            L_biomass<-biomass[,s,w>w_cut]

            status_LFB_s<-sum(L_biomass[99,])
            status_B_s<-sum(biomass[99,s,])
            status_LFB<-status_LFB+status_LFB_s
            status_B<-status_B+status_B_s

            BH_LFB_s<-sum(L_biomass[699,])
            BH_B_s<-sum(biomass[699,s,])
            BH_LFB<-BH_LFB+BH_LFB_s
            BH_B<-BH_B+BH_B_s
    }

    status_LFY_t<-status_LFY
    status_LFY<-status_LFY/status_Y
    BH_LFY_t<-BH_LFY
    BH_LFY<-BH_LFY/BH_Y

    status_LFB_t<-status_LFB
    status_LFB<-status_LFB/status_B
    BH_LFB_t<-BH_LFB
    BH_LFB<-BH_LFB/BH_B

    #insert values into p_grouped
    rows <- p_grouped$c == c_val & p_grouped$alpha == alpha_val

    #large fish proportion in yield
    p_grouped$status_LFY_t[rows]<-status_LFY_t
    p_grouped$status_LFY_t[rows]<-BH_LFY_t
    p_grouped$status_LFY[rows]<-status_LFY
    p_grouped$BH_LFY[rows]<-BH_LFY

    #large fish proportion in system biomass
    p_grouped$status_LFB_t[rows]<-status_LFB_t
    p_grouped$BH_LFB_t[rows]<-BH_LFB_t
    p_grouped$status_LFB[rows]<-status_LFB
    p_grouped$BH_LFB[rows]<-BH_LFB
    }
}

#saveRDS(p_grouped, "/Users/jessicawestworth/Desktop/BH simulations/p_grouped.rds")
p_grouped<-readRDS("/Users/jessicawestworth/Desktop/BH simulations/p_grouped.rds")
#to extract specific simulations: sim <- sim_list[["c_0.5_a_0.3"]]
#example extract c=0.5 and alpha=0.3

plot_ly(data=p_grouped, x= ~c,y= ~alpha,z= ~BH_Y, type= "heatmap")
plot_ly(data=p_grouped, x= ~c,y= ~alpha,z= ~BH_B, type= "heatmap")
plot_ly(data=p_grouped, x= ~c,y= ~alpha,z= ~BH_SSB, type= "heatmap")
plot_ly(data=p_grouped, x= ~c,y= ~alpha,z= ~BH_N, type= "heatmap")
plot_ly(data=p_grouped, x= ~c,y= ~alpha,z= ~BH_LFY_t, type= "heatmap")
plot_ly(data=p_grouped, x= ~c,y= ~alpha,z= ~BH_LFB_t, type= "heatmap")
plot_ly(data=p_grouped, x= ~c,y= ~alpha,z= ~BH_LFY, type= "heatmap")
plot_ly(data=p_grouped, x= ~c,y= ~alpha,z= ~BH_LFB, type= "heatmap")


















#Increasing and decreasing effort for current fishing
t_duration <- 700
t_steadied <- 100
t_blended <- 200

species_vec <- species_params(params)$species

current_fish <- expand.grid(
    species = species_vec,
    alpha = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10)
)

current_fish$status_Y <- NA
current_fish$proj_Y <- NA
current_fish$status_B <- NA
current_fish$proj_B <- NA
current_fish$status_SSB <- NA
current_fish$proj_SSB <- NA
current_fish$status_N <- NA
current_fish$proj_N <- NA

current_fish_Yield_list <-list()
current_fish_sim_list <- list()

for(alpha_val in unique(current_fish$alpha)){
    my_blender <- make_blended_increase_effort(
        t_max_blend = t_blended,
        t_steady = t_steadied,
        alpha_max = alpha_val
    )

    sim_start <- setRateFunction(ps, "FMort", "my_blender")
    sim <- project(sim_start, t_max = t_duration, effort = 1)

    current_fish_sim_list[[paste0("effort_",alpha_val)]] <- sim

    Yield <- getYield(sim)
    current_fish_Yield_list[[paste0("effort_",alpha_val)]] <- Yield

    Biomass <- getBiomass(sim)
    SSB <- getSSB(sim)
    N <- getN(sim)

    for(species in unique(species_params(ps)$species)){
        rows <- current_fish$species == species & current_fish$alpha == alpha_val
        current_fish$status_Y[rows] <- Yield[99, species]
        current_fish$proj_Y[rows] <- Yield[699, species]
        current_fish$status_B[rows] <- Biomass[99, species]
        current_fish$proj_B[rows] <- Biomass[699, species]

        current_fish$status_SSB[rows] <- SSB[99, species]
        current_fish$proj_SSB[rows] <- SSB[699, species]

        current_fish$status_N[rows] <- N[99, species]
        current_fish$proj_N[rows] <- N[699, species]
    }
}
#saveRDS(current_fish_sim_list,"/Users/jessicawestworth/Desktop/BH simulations/current_fish_sim_list.rds")
#saveRDS(current_fish_Yield_list,"/Users/jessicawestworth/Desktop/BH simulations/current_fish_Yield_list.rds")
#saveRDS(current_fish,"/Users/jessicawestworth/Desktop/BH simulations/current_fish.rds")
current_fish<-readRDS("/Users/jessicawestworth/Desktop/BH simulations/current_fish.rds")

ggplot(current_fish, aes(x=alpha,y=log(proj_Y), group=species))+
    geom_line(aes(color=species))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(y= "Log Yield (g/year/m^2)", x= "Effort", color="Species")+
    coord_cartesian(ylim = c(-10, 5))

ggplot(current_fish, aes(x=alpha,y=log(proj_B), group=species))+
    geom_line(aes(color=species))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(y= "Log Biomass (g/m^2)", x= "Effort", color="Species")+
    coord_cartesian(ylim = c(-8, 5))

ggplot(current_fish, aes(x=alpha,y=log(proj_SSB), group=species))+
    geom_line(aes(color=species))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(y= "Log Spawning Stock Biomass (g/m^2)", x= "Effort", color="Species")+
    coord_cartesian(ylim = c(-10, 5))

ggplot(current_fish, aes(x=alpha,y=log(proj_N), group=species))+
    geom_line(aes(color=species))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(y= "Log Number of Individuals (per m^2)", x= "Effort", color="Species")+
    coord_cartesian(ylim = c(-10, 5))

#group and get system LFI values
current_fish_grouped<-current_fish%>%
    group_by(alpha)%>%
    summarise(status_Y=sum(status_Y),
              proj_Y=sum(proj_Y),
              status_B=sum(status_B),
              proj_B=sum(proj_B),
              status_SSB=sum(status_SSB),
              proj_SSB=sum(proj_SSB),
              status_N=sum(status_N),
              proj_N=sum(proj_N)
    )

#current_fish_sim_list<-readRDS("/Users/jessicawestworth/Desktop/BH simulations/current_fish_sim_list.rds")

for(alpha_val in unique(current_fish_grouped$alpha)){
        sim<-current_fish_sim_list[[paste0("effort_",alpha_val)]]
        my_blender <- make_blended_increase_effort(
            t_max_blend = t_blended,
            t_steady = t_steadied,
            alpha_max = alpha_val
        )

        sim_start <- setRateFunction(ps, "FMort", "my_blender")

        sim<-current_fish_sim_list[[paste0("effort_",alpha_val)]]
        f<-getFMort(sim, drop = FALSE)
        n<-sim@n

        biomass <- sweep(sim@n, 3, sim@params@w * sim@params@dw, "*")
        yield<-f * biomass
        catch<-f*n
        number<-n

        status_LFY_s<-0
        status_Y_s<-0
        status_LFY<-0
        status_Y<-0

        proj_LFY_s<-0
        proj_Y_s<-0
        proj_LFY<-0
        proj_Y<-0

        status_LFB_s<-0
        status_B_s<-0
        status_LFB<-0
        status_B<-0

        proj_LFB_s<-0
        proj_B_s<-0
        proj_LFB<-0
        proj_B<-0

        for(s in 1:nrow(sim@params@species_params)){
            l<-50
            a<-sim@params@species_params$a[s]
            b<-sim@params@species_params$b[s]
            w_cut<-a*(l^b)

            #yields
            L_yield<-yield[,s,w>w_cut]

            status_LFY_s<-sum(L_yield[99,])
            status_Y_s<-sum(yield[99,s,])
            status_LFY<-status_LFY+status_LFY_s
            status_Y<-status_Y+status_Y_s

            proj_LFY_s<-sum(L_yield[699,])
            proj_Y_s<-sum(yield[699,s,])
            proj_LFY<-proj_LFY+proj_LFY_s
            proj_Y<-proj_Y+proj_Y_s

            #Biomass
            L_biomass<-biomass[,s,w>w_cut]

            status_LFB_s<-sum(L_biomass[99,])
            status_B_s<-sum(biomass[99,s,])
            status_LFB<-status_LFB+status_LFB_s
            status_B<-status_B+status_B_s

            proj_LFB_s<-sum(L_biomass[699,])
            proj_B_s<-sum(biomass[699,s,])
            proj_LFB<-proj_LFB+proj_LFB_s
            proj_B<-proj_B+proj_B_s
        }

    status_LFY<-status_LFY/status_Y
    proj_LFY<-proj_LFY/proj_Y

    status_LFB<-status_LFB/status_B
    proj_LFB<-proj_LFB/proj_B

    #insert values into current_fish_grouped
    rows <- current_fish_grouped$alpha == alpha_val

    #large fish proportion in yield
    current_fish_grouped$status_LFY[rows]<-status_LFY
    current_fish_grouped$proj_LFY[rows]<-proj_LFY

    #large fish proportion in system biomass
    current_fish_grouped$status_LFB[rows]<-status_LFB
    current_fish_grouped$proj_LFB[rows]<-proj_LFB
}

#saveRDS(current_fish_grouped, "/Users/jessicawestworth/Desktop/BH simulations/current_fish_grouped.rds")
current_fish_grouped<-readRDS("/Users/jessicawestworth/Desktop/BH simulations/current_fish_grouped.rds")

plot_ly(data=current_fish_grouped, x= ~alpha,y= ~proj_Y, mode='lines')%>%
    layout(yaxis = list(title = 'Yield'), xaxis = list(title = 'Effort'),
           shapes = list(vline(1)))%>%
    add_annotations(text="status quo", x =1.9, y=18,showarrow = FALSE)

plot_ly(data=current_fish_grouped, x= ~alpha,y= ~proj_B, mode='lines')%>%
    layout(yaxis = list(title = 'Biomass'), xaxis = list(title = 'Effort'),
           shapes = list(vline(1)))%>%
    add_annotations(text="status quo", x =1.9, y=45,showarrow = FALSE)

plot_ly(data=current_fish_grouped, x= ~alpha,y= ~proj_SSB, mode='lines')%>%
    layout(yaxis = list(title = 'Spawning Stock Biomass'), xaxis = list(title = 'Effort'),
           shapes = list(vline(1)))%>%
    add_annotations(text="status quo", x =1.9, y=19,showarrow = FALSE)

plot_ly(data=current_fish_grouped, x= ~alpha,y= ~proj_N, mode='lines')%>%
    layout(yaxis = list(title = 'Number of Individuals'), xaxis = list(title = 'Effort'),
           shapes = list(vline(1)))%>%
    add_annotations(text="status quo", x =1.9, y=71,showarrow = FALSE)

plot_ly(data=current_fish_grouped, x= ~alpha,y= ~proj_LFY, mode='lines')%>%
    layout(yaxis = list(title = 'Proportion Large Fish in Yield'), xaxis = list(title = 'Effort'),
           shapes = list(vline(1)))%>%
    add_annotations(text="status quo", x =1.9, y=0.6,showarrow = FALSE)

plot_ly(data=current_fish_grouped, x= ~alpha,y= ~proj_LFB, mode='lines')%>%
    layout(yaxis = list(title = 'Proportion Large Fish in Biomass'), xaxis = list(title = 'Effort'),
           shapes = list(vline(1)))%>%
    add_annotations(text="status quo", x =1.9, y=0.12,showarrow = FALSE)

#Cod goes extinct at effort = 2
#Haddock goes extinct effort = 6
#Boarfish goes extinct effort = 7
#Monkfish begins to go extinct effort = 10

#determine important full BH states
pg<-p_grouped%>%filter(alpha==1)
colSums(pg)
pg$Total_BH_Y<-249.094755
pg$Total_BH_B<-352.493550
pg$Total_BH_SSB<-153.522653
pg$Total_BH_N<-920.604007
pg$Total_BH_LFB<-4.122696

pb<-na.omit(pg)
colSums(pb)

pg$Total_BH_LFY<-0.7527794

pg$Y<-pg$BH_Y/pg$Total_BH_Y
pg$N<-pg$BH_N/pg$Total_BH_N
pg$B<-pg$BH_B/pg$Total_BH_B
pg$SSB<-pg$BH_SSB/pg$Total_BH_SSB
pg$LFB<-pg$BH_LFB/pg$Total_BH_LFB
pg$LFY<-pg$BH_LFY/pg$Total_BH_LFY

pg[1,26]<-NA
colSums(pg)

plot(pg$c, pg$N, type="l",col="red", ylim=c(0,0.15))
points(pg$c, pg$B, type="l", col="orange")
points(pg$c, pg$Y, type="l", col="yellow3")
points(pg$c, pg$SSB, type="l", col="green3")
points(pg$c, pg$LFY, type="l", col="blue")
points(pg$c, pg$LFB, type="l", col="lightblue")

pg[1,29]<-0

pg$points<-pg$N+pg$B+pg$SSB+pg$LFY+pg$LFB
max(pg$points)




