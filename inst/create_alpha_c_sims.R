#Create Model ssBH Simulations
t_duration <- 700
t_steadied <- 100
t_blended <- 200

#species from the model
species_vec <- species_params(params)$species

#define fishing intensity (c) and weightings (alpha)
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

#Run the simulation code and create the model and save it for each c and alpha
#combination
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

#Sim list with the 220 models contained within:
sim_list<-readRDS("/Users/jessicawestworth/Desktop/BH simulations/sim_list_extended/sim_list.rds")

#Species specific Metrics values from the Simulation
p<-readRDS("/Users/jessicawestworth/Desktop/BH simulations/p_extended_time/p_extended_time.rds")

#Plotting how in full ssBH regimes how the steady state is different across
#simulations with different fishing intensities
p_full<-p%>%
    filter(c>0,alpha==1)

g1<-ggplot(p_full, aes(x=c,y=log(BH_Y), group=species))+
    geom_line(aes(color=species))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(y= "Log Yield (g/year/m^2)", x= "c", color="Species")

g2<-ggplot(p_full, aes(x=c,y=log(BH_B), group=species))+
    geom_line(aes(color=species))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(y= "Log Biomass (g/m^2)", x= "c", color="Species")

g3<-ggplot(p_full, aes(x=c,y=log(BH_SSB), group=species))+
    geom_line(aes(color=species))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(y= "Log Spawning Stock Biomass (g/m^2)", x= "c", color="Species")

g4<-ggplot(p_full, aes(x=c,y=log(BH_N), group=species))+
    geom_line(aes(color=species))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(y= "Log Number of Individuals (per m^2)", x= "c", color="Species")

legend <- get_legend(
    g1 +
        theme(legend.position = "right") +
        guides(color = guide_legend(ncol = 1))  # force vertical
)

g1 <- g1 + theme(legend.position = "none")
g2 <- g2 + theme(legend.position = "none")
g3 <- g3 + theme(legend.position = "none")
g4 <- g4 + theme(legend.position = "none")

plot_grid(
    plot_grid(g1, g2, g3, g4, ncol = 2),
    legend,
    ncol = 2,
    rel_widths = c(1, 0.25)
)

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
    p_grouped$BH_LFY_t[rows]<-BH_LFY_t
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
p_grouped<-p_grouped%>%filter(c>0)

plot_ly(data=p_grouped, x= ~c,y= ~alpha,z= ~BH_Y, type= "heatmap", colorbar = list(title = "Y"))%>%
    layout(yaxis = list(title = ' '), xaxis = list(title = ' '))
plot_ly(data=p_grouped, x= ~c,y= ~alpha,z= ~BH_B, type= "heatmap", colorbar = list(title = "B"))%>%
    layout(yaxis = list(title = ' '), xaxis = list(title = ' '))
plot_ly(data=p_grouped, x= ~c,y= ~alpha,z= ~BH_SSB, type= "heatmap", colorbar = list(title = "SSB"))%>%
    layout(yaxis = list(title = ' '), xaxis = list(title = ' '))
plot_ly(data=p_grouped, x= ~c,y= ~alpha,z= ~BH_N, type= "heatmap", colorbar = list(title = "N"))%>%
    layout(yaxis = list(title = ' '), xaxis = list(title = ' '))
plot_ly(data=p_grouped, x= ~c,y= ~alpha,z= ~BH_LFY_t, type= "heatmap", colorbar = list(title = "BLFY"))%>%
    layout(yaxis = list(title = ' '), xaxis = list(title = ' '))
plot_ly(data=p_grouped, x= ~c,y= ~alpha,z= ~BH_LFB_t, type= "heatmap", colorbar = list(title = "BLFS"))%>%
    layout(yaxis = list(title = ' '), xaxis = list(title = ' '))
plot_ly(data=p_grouped, x= ~c,y= ~alpha,z= ~BH_LFY, type= "heatmap", colorbar = list(title = "PBLFY"))%>%
    layout(yaxis = list(title = ' '), xaxis = list(title = ' '))
plot_ly(data=p_grouped, x= ~c,y= ~alpha,z= ~BH_LFB, type= "heatmap", colorbar = list(title = "PBLFS"))%>%
    layout(yaxis = list(title = ' '), xaxis = list(title = ' '))

whole_BH<-p_grouped%>%filter(c>0, alpha==1)

g1<-ggplot(whole_BH, aes(x=c,y=BH_Y))+
    geom_line()+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(y= "Y", x= "c")
g2<-ggplot(whole_BH, aes(x=c,y=BH_B))+
    geom_line()+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(y= "B", x= "c")
g3<-ggplot(whole_BH, aes(x=c,y=BH_SSB))+
    geom_line()+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(y= "SSB", x= "c")
g4<-ggplot(whole_BH, aes(x=c,y=BH_N))+
    geom_line()+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(y= "N", x= "c")
g5<-ggplot(whole_BH, aes(x=c,y=BH_LFB))+
    geom_line()+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(y= "PBLFS", x= "c")
g6<-ggplot(whole_BH, aes(x=c,y=BH_LFB_t))+
    geom_line()+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(y= "BLFS", x= "c")
g7<-ggplot(whole_BH, aes(x=c,y=BH_LFY))+
    geom_line()+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(y= "PBLFY", x= "c")
g8<-ggplot(whole_BH, aes(x=c,y=BH_LFY_t))+
    geom_line()+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(y= "BLFY", x= "c")


plot<-plot_grid(g1,g2, g3, g4,g5,g6,g7,g8, ncol = 2, label_y = 1.1, labels = letters[1:8], rel_heights = c(1, 1.1))

ggdraw() +
    draw_plot(plot, y = 0, height = 0.95)

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

#determine important whole BH states
pg<-p_grouped%>%filter(alpha==1)
pg[1,13]<-0
pg$Y<-(pg$BH_Y-min(pg$BH_Y))/(max(pg$BH_Y)-min(pg$BH_Y))
pg$N<-(pg$BH_N-min(pg$BH_N))/(max(pg$BH_N)-min(pg$BH_N))
pg$B<-(pg$BH_B-min(pg$BH_B))/(max(pg$BH_B)-min(pg$BH_B))
pg$SSB<-(pg$BH_SSB-min(pg$BH_SSB))/(max(pg$BH_SSB)-min(pg$BH_SSB))
pg$LFB<-(pg$BH_LFB-min(pg$BH_LFB))/(max(pg$BH_LFB)-min(pg$BH_LFB))
pg$LFB_t<-(pg$BH_LFB_t-min(pg$BH_LFB_t))/(max(pg$BH_LFB_t)-min(pg$BH_LFB_t))
pg$LFY<-(pg$BH_LFY-min(pg$BH_LFY))/(max(pg$BH_LFY)-min(pg$BH_LFY))
pg$LFY_t<-(pg$BH_LFY_t-min(pg$BH_LFY_t))/(max(pg$BH_LFY_t)-min(pg$BH_LFY_t))

pg$conservation<-(pg$N+pg$B+pg$SSB+pg$LFB_t+pg$LFB)/5
pg$economic<-(pg$Y+pg$LFY+pg$LFY_t)/3

pg$points<-pg$conservation+pg$economic

max(pg$points)

#determine best state
pg<-p_grouped
pg[11,13]<-0
pg$Y<-(pg$BH_Y-min(pg$BH_Y))/(max(pg$BH_Y)-min(pg$BH_Y))
pg$N<-(pg$BH_N-min(pg$BH_N))/(max(pg$BH_N)-min(pg$BH_N))
pg$B<-(pg$BH_B-min(pg$BH_B))/(max(pg$BH_B)-min(pg$BH_B))
pg$SSB<-(pg$BH_SSB-min(pg$BH_SSB))/(max(pg$BH_SSB)-min(pg$BH_SSB))
pg$LFB<-(pg$BH_LFB-min(pg$BH_LFB))/(max(pg$BH_LFB)-min(pg$BH_LFB))
pg$LFB_t<-(pg$BH_LFB_t-min(pg$BH_LFB_t))/(max(pg$BH_LFB_t)-min(pg$BH_LFB_t))
pg$LFY<-(pg$BH_LFY-min(pg$BH_LFY))/(max(pg$BH_LFY)-min(pg$BH_LFY))
pg$LFY_t<-(pg$BH_LFY_t-min(pg$BH_LFY_t))/(max(pg$BH_LFY_t)-min(pg$BH_LFY_t))

pg$conservation<-(pg$N+pg$B+pg$SSB+pg$LFB_t+pg$LFB)/5
pg$economic<-(pg$Y+pg$LFY+pg$LFY_t)/3

pg$points<-pg$conservation+pg$economic

max(pg$points)

pg<-pg%>%filter(c>0)
plot_ly(data=pg, x= ~c,y= ~alpha,z= ~points, type= "heatmap", colorbar = list(title = "MCCS"))%>%
    layout(yaxis = list(title = 'A (weighting)'), xaxis = list(title = 'c (fishing intensity)'))



#plot the Yields
my_blender <- make_blended_ssBH_FMort(
    t_max_blend = t_blended,
    target_c = 0.4,
    t_steady = t_steadied,
    alpha_max = 1
)

sim_BH_start <- setRateFunction(ps, "FMort", "my_blender")

#f_whole<-getFMort(sim_list[["c_0.4_a_1"]], drop = FALSE)

my_blender <- make_blended_ssBH_FMort(
    t_max_blend = t_blended,
    target_c = 0.4,
    t_steady = t_steadied,
    alpha_max = 0.5
)

sim_BH_start <- setRateFunction(ps, "FMort", "my_blender")


#f_hybrid<-getFMort(sim_list[["c_0.4_a_0.5"]], drop = FALSE)

f_status<-f_whole[99,,]
f_status<-as.data.frame.table(f_status)
f_status$w<-as.character(f_status$w)
f_status$w<-as.numeric(f_status$w)

#size dependent fishing mortality
g1<-ggplot(f_status, aes(x=w,y=Freq, group=sp))+
    geom_line(aes(color=sp))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_log10()+
    labs(title="Status Quo",y= "F", x= "w (g)", color="Species")

f_whole_data<-f_whole[699,,]
f_whole_data<-as.data.frame.table(f_whole_data)
f_whole_data$w<-as.character(f_whole_data$w)
f_whole_data$w<-as.numeric(f_whole_data$w)
g2<-ggplot(f_whole_data, aes(x=w,y=Freq, group=sp))+
    geom_line(aes(color=sp))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_log10()+
    labs(title="Whole ssBH",y= "F", x= "w (g)", color="Species")


f_hybrid_data<-f_hybrid[699,,]
f_hybrid_data<-as.data.frame.table(f_hybrid_data)
f_hybrid_data$w<-as.character(f_hybrid_data$w)
f_hybrid_data$w<-as.numeric(f_hybrid_data$w)
g3<-ggplot(f_hybrid_data, aes(x=w,y=Freq, group=sp))+
    geom_line(aes(color=sp))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_log10()+
    labs(title="Hybrid",y= "F", x= "w (g)", color="Species")

g1 <- g1 + theme(legend.position = "none")
g2 <- g2 + theme(legend.position = "none")
g3 <- g3 + theme(legend.position = "none")

plot_grid(
    plot_grid(g1, g2, g3, ncol = 1),
    legend,
    ncol = 2,
    rel_widths = c(1, 0.25)
)


#Fishing Mortaltiy (1/m^2)
#Weight (g)
biomass_whole <- sweep(sim_list[["c_0.4_a_1"]]@n, 3, sim_list[["c_0.4_a_1"]]@params@w * sim_list[["c_0.4_a_1"]]@params@dw, "*")
biomass_hybrid <- sweep(sim_list[["c_0.4_a_0.5"]]@n, 3, sim_list[["c_0.4_a_0.5"]]@params@w * sim_list[["c_0.4_a_0.5"]]@params@dw, "*")

yield_hybrid<-f_whole[699,,]*biomass_hybrid[699,,]
yield_hybrid<-as.data.frame.table(yield_hybrid)
yield_hybrid$w<-as.character(yield_hybrid$w)
yield_hybrid$w<-as.numeric(yield_hybrid$w)

yield_whole<-f_whole[699,,]*biomass_whole[699,,]
yield_whole<-as.data.frame.table(yield_whole)
yield_whole$w<-as.character(yield_whole$w)
yield_whole$w<-as.numeric(yield_whole$w)

yield_status<-f_whole[99,,]*biomass_whole[99,,]
yield_status<-as.data.frame.table(yield_status)
yield_status$w<-as.character(yield_status$w)
yield_status$w<-as.numeric(yield_status$w)

g4<-ggplot(yield_status, aes(x=w,y=Freq, group=sp))+
    geom_line(aes(color=sp))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_log10()+
    labs(title= "Status Quo", y= "Yield", x= "w (g)", color="Species")+
    theme(legend.position = "none")

g5<-ggplot(yield_whole, aes(x=w,y=Freq, group=sp))+
    geom_line(aes(color=sp))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_log10()+
    labs(title= "Whole ssBH",y= "Yield", x= "w (g)", color="Species")+
    theme(legend.position = "none")

g6<-ggplot(yield_hybrid, aes(x=w,y=Freq, group=sp))+
    geom_line(aes(color=sp))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_log10()+
    labs(title= "Hybrid", y= "Yield", x= "w (g)", color="Species")+
    theme(legend.position = "none")

legend <- get_legend(
    g1 +
        theme(legend.position = "right") +
        guides(color = guide_legend(ncol = 1))  # force vertical
)

g1 <- g1 + theme(legend.position = "none")
g2 <- g2 + theme(legend.position = "none")
g3 <- g3 + theme(legend.position = "none")

plot_grid(
    plot_grid(g1,g4,g2,g5,g3,g6, ncol = 2),
    NULL,
    legend,
    ncol = 3,
    rel_widths = c(1, 0.05,0.25)
)


biomass_hybrid<-biomass_hybrid[699,,]
biomass_hybrid<-as.data.frame.table(biomass_hybrid)
biomass_hybrid$w<-as.character(biomass_hybrid$w)
biomass_hybrid$w<-as.numeric(biomass_hybrid$w)

biomass_status<-biomass_whole[99,,]
biomass_status<-as.data.frame.table(biomass_status)
biomass_status$w<-as.character(biomass_status$w)
biomass_status$w<-as.numeric(biomass_status$w)

biomass_whole<-biomass_whole[699,,]
biomass_whole<-as.data.frame.table(biomass_whole)
biomass_whole$w<-as.character(biomass_whole$w)
biomass_whole$w<-as.numeric(biomass_whole$w)



g7<-ggplot(biomass_status, aes(x=w,y=Freq, group=sp))+
    geom_line(aes(color=sp))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_log10()+
    labs(title= "Status Quo", y= "Biomass", x= "w (g)", color="Species")+
    theme(legend.position = "none")

g8<-ggplot(biomass_whole, aes(x=w,y=Freq, group=sp))+
    geom_line(aes(color=sp))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_log10()+
    labs(title= "Whole ssBH",y= "Biomass", x= "w (g)", color="Species")+
    theme(legend.position = "none")

g9<-ggplot(biomass_hybrid, aes(x=w,y=Freq, group=sp))+
    geom_line(aes(color=sp))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_log10()+
    labs(title= "Hybrid", y= "Biomass", x= "w (g)", color="Species")+
    theme(legend.position = "none")

plot_grid(
    plot_grid(g1,g4,g7,g2,g5,g8,g3,g6,g9, ncol = 3),
    NULL,
    legend,
    ncol = 3,
    rel_widths = c(1, 0.05,0.20)
)

#this weight dependency allows for these yields to be so so high. Because mackerel
#is being fished between

#each species plot with the yield from each regime
yield_hybrid$type<-"hybrid"
yield_whole$type<-"whole"
yield_status$type<-"status"

yield<-rbind(yield_hybrid,yield_whole,yield_status)

horse_mackerel_yield<-subset(yield, sp=="Horse mackerel")
g1<-ggplot(horse_mackerel_yield, aes(x=w,y=Freq, group=type))+
    geom_line(aes(color=type))+
    theme_cowplot(12)+
    labs(title= "Horse mackerel")+
    theme(legend.position = "none")

mackerel_yield<-subset(yield, sp=="Mackerel")
g2<-ggplot(mackerel_yield, aes(x=w,y=Freq, group=type))+
    geom_line(aes(color=type))+
    theme_cowplot(12)+
    labs(title= "Mackerel")+
    theme(legend.position = "none")

blue_whiting_yield<-subset(yield, sp=="Blue whiting")
g3<-ggplot(blue_whiting_yield, aes(x=w,y=Freq, group=type))+
    geom_line(aes(color=type))+
    theme_cowplot(12)+
    labs(title= "Blue whiting")+
    theme(legend.position = "none")

boarfish_yield<-subset(yield, sp=="Boarfish")
g4<-ggplot(boarfish_yield, aes(x=w,y=Freq, group=type))+
    geom_line(aes(color=type))+
    theme_cowplot(12)+
    labs(title= "Boarfish")+
    theme(legend.position = "none")

cod_yield<-subset(yield, sp=="Cod")
g5<-ggplot(cod_yield, aes(x=w,y=Freq, group=type))+
    geom_line(aes(color=type))+
    theme_cowplot(12)+
    labs(title= "Cod")+
    theme(legend.position = "none")

haddock_yield<-subset(yield, sp=="Haddock")
g6<-ggplot(haddock_yield, aes(x=w,y=Freq, group=type))+
    geom_line(aes(color=type))+
    theme_cowplot(12)+
    labs(title= "Haddock")+
    theme(legend.position = "none")

hake_yield<-subset(yield, sp=="Hake")
g7<-ggplot(hake_yield, aes(x=w,y=Freq, group=type))+
    geom_line(aes(color=type))+
    theme_cowplot(12)+
    labs(title= "Hake")+
    theme(legend.position = "none")

herring_yield<-subset(yield, sp=="Herring")
g8<-ggplot(herring_yield, aes(x=w,y=Freq, group=type))+
    geom_line(aes(color=type))+
    theme_cowplot(12)+
    labs(title= "Herring")+
    theme(legend.position = "none")

megrim_yield<-subset(yield, sp=="Megrim")
g9<-ggplot(megrim_yield, aes(x=w,y=Freq, group=type))+
    geom_line(aes(color=type))+
    theme_cowplot(12)+
    labs(title= "Megrim")+
    theme(legend.position = "none")

monkfish_yield<-subset(yield, sp=="Monkfish")
g10<-ggplot(monkfish_yield, aes(x=w,y=Freq, group=type))+
    geom_line(aes(color=type))+
    theme_cowplot(12)+
    labs(title= "Monkfish")+
    theme(legend.position = "none")

plaice_yield<-subset(yield, sp=="Plaice")
g11<-ggplot(plaice_yield, aes(x=w,y=Freq, group=type))+
    geom_line(aes(color=type))+
    theme_cowplot(12)+
    labs(title= "Plaice")+
    theme(legend.position = "none")

red_gurnard_yield<-subset(yield, sp=="Red gurnard")
g12<-ggplot(red_gurnard_yield, aes(x=w,y=Freq, group=type))+
    geom_line(aes(color=type))+
    theme_cowplot(12)+
    labs(title= "Red gurnard")+
    theme(legend.position = "none")

sole_yield<-subset(yield, sp=="Sole")
g13<-ggplot(sole_yield, aes(x=w,y=Freq, group=type))+
    geom_line(aes(color=type))+
    theme_cowplot(12)+
    labs(title= "Sole")+
    theme(legend.position = "none")

whiting_yield<-subset(yield, sp=="Whiting")
g14<-ggplot(whiting_yield, aes(x=w,y=Freq, group=type))+
    geom_line(aes(color=type))+
    theme_cowplot(12)+
    labs(title= "Whiting")+
    theme(legend.position = "none")

legend <- get_legend(
    g1 +
        theme(legend.position = "right") +
        guides(color = guide_legend(ncol = 1))  # force vertical
)


plot_grid(
    plot_grid(g1, g2, g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14, ncol = 3),
    legend,
    ncol = 2,
    rel_widths = c(2, 0.25)
)



#plot the Biomass
biomass_whole <- sweep(sim_list[["c_0.4_a_1"]]@n, 3, sim_list[["c_0.4_a_1"]]@params@w * sim_list[["c_0.4_a_1"]]@params@dw, "*")
biomass_hybrid <- sweep(sim_list[["c_0.4_a_0.5"]]@n, 3, sim_list[["c_0.4_a_0.5"]]@params@w * sim_list[["c_0.4_a_0.5"]]@params@dw, "*")


yield_whole<-biomass_whole*f_whole
yield_whole<-as.data.frame.table(yield_whole)
yield_whole$time<-as.character(yield_whole$time)
yield_whole$time<-as.numeric(yield_whole$time)
yield_whole<-yield_whole%>%group_by(time, sp)%>%
    summarise(Freq=sum(Freq))

yield_hybrid<-biomass_hybrid*f_hybrid
yield_hybrid<-as.data.frame.table(yield_hybrid)
yield_hybrid$time<-as.character(yield_hybrid$time)
yield_hybrid$time<-as.numeric(yield_hybrid$time)
yield_hybrid<-yield_hybrid%>%group_by(time, sp)%>%
    summarise(Freq=sum(Freq))

g1<-ggplot(yield_whole, aes(x=time,y=Freq, group=sp))+
    geom_line(linewidth = 0.8,aes(color=sp))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_log10(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(title="Whole ssBH", y= "Yield", x= "Year", color="Species")+
    theme(legend.position = "none")

g2<-ggplot(yield_hybrid, aes(x=time,y=Freq, group=sp))+
    geom_line(linewidth = 0.8,aes(color=sp))+
    scale_colour_manual(values = params@linecolour)+
    theme_cowplot(12)+
    scale_y_log10(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(title="Hybrid", y= "Yield", x= "Year", color="Species")+
    theme(legend.position = "none")

g3<-plotBiomass(sim_list[["c_0.4_a_1"]])+theme_cowplot(12)+
    scale_y_log10(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(title="Whole ssBH")+
    theme(legend.position = "none")
g4<-plotBiomass(sim_list[["c_0.4_a_0.5"]])+theme_cowplot(12)+
    scale_y_log10(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    labs(title="Hybrid")+
    theme(legend.position = "none")
legend <- get_legend(
    g1 +
        theme(legend.position = "right") +
        guides(color = guide_legend(ncol = 1))  # force vertical
)
plot_grid(
    plot_grid(g1, g2, g3,g4, ncol = 2),
    legend,
    ncol = 2,
    rel_widths = c(2, 0.39)
)

