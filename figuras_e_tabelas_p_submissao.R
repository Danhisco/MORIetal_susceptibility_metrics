# primeiro rodar: source.R
# Instalar e carregar o pacote extrafont
# install.packages("extrafont")
library(extrafont)

# Importar fontes do Windows para o R
font_import()  # Pode demorar alguns minutos

# Registrar as fontes
loadfonts(device = "win")

#---------------------------------------------------------------------
#----------------------------- FIGURAS -------------------------------
#---------------------------------------------------------------------

#--------------------------- matriz de contato -----------------------
p <- df_contato %>%
    mutate(y = factor(y,
                      levels=unique(df_contato$y)[5:1])
    ) %>%
    ggplot(aes(x=x,y=y,fill=contact_rate)) +
    geom_tile() +
    geom_label(aes(label=round(contact_rate,2)),
               fill="white",
               size = 11/.pt) +
    scale_fill_viridis_c("") +
    labs(x="",y="") +
    scale_y_discrete(expand=FALSE,position = "right") +
    scale_x_discrete(expand=FALSE) +
    theme(legend.position = "none",
          text = element_text(family = "Times New Roman", size = 12))
ggsave(filename = "8_pastas_auxiliares/artigo_suscetiveis/figuras/figSI1.jpeg",
       plot=p,
       width = 6.37,
       height = 7.15)
# f_img("8_pastas_auxiliares/artigo_suscetiveis/figuras/mat_contato.jpeg")


#--------------------------- raw values ------------------------------
f_plot_raw <- \(dfi,
                fillcol="p_sus",
                vrange=c(0,100),
                fillname="Prop. Suscetíveis",
                facetlab="TESTE",
                stripsize=11,
                textsize=5,
                vtrans="identity"){
    if(is.null(attr(dfi,"sf_column"))) dfi <- inner_join(x=dfnb,y=dfi)
    mutate(dfi,
           facetlab=facetlab) %>%
        ggplot() +
        geom_sf(aes(fill=.data[[fillcol]]),
                linetype = "dashed",
                linewidth=0.25,
                color="darkgray") +
        scale_fill_distiller(fillname,
                             palette="Spectral",
                             limits=vrange,
                             transform=vtrans) +
        geom_sf(data=mutate(shp_gve,
                            facetlab=facetlab),
                col="black",
                alpha=0.1,
                linewidth=0.5) +
        coord_sf(expand = FALSE) +
        facet_wrap(~facetlab) +
        theme_bw() +
        annotation_north_arrow(location = "bl",
                               which_north = "true",
                               pad_x = unit(0.75, "in"),
                               pad_y = unit(0.5, "in"),
                               style = north_arrow_fancy_orienteering) +
        theme(legend.position = "bottom",
              text = element_text(family = "Times New Roman", size = 12),
              legend.background = element_rect(fill = "white", color = "black"),
              legend.box.background = element_rect(fill = "white", color = "black"),
              legend.margin = margin(10, 10, 10, 15),
              legend.key = element_rect(fill = "white"),
              strip.text = element_text(size=stripsize,
                                        hjust=0,
                                        family = "Times New Roman",
                                        face = "bold"),
              strip.background = element_rect(fill="white",color="white")
              )
}
f_wrap_plots <- \(lp){
    Reduce("+",lp) + guide_area() +
        plot_layout(guides = "collect",ncol=3) &
        theme(legend.position = "right")
}

##@ overall population
vrange <- range(df_ad$POP)
dfp <- filter(df_ad,fxet!="total")
lp <- split(dfp,dfp$fxet)
names(lp) <- paste0(LETTERS[1:length(lp)],") ",names(lp))
lp <- lapply(names(lp),\(li){
    dfi <- mutate(lp[[li]],
                  logPOP=log(POP))
    f_plot_raw(dfi,
               facetlab=li,
               fillcol="POP",
               vtrans = "log10",
               vrange = vrange,
               fillname="Population size\n(number of ind.)"
    )

})
pfinal_rawpop <- f_wrap_plots(lp)
ggsave(filename="8_pastas_auxiliares/artigo_suscetiveis/figuras/map_raw_values/overall_population.jpeg",
       plot=pfinal_rawpop,
       dpi=200,
       width=14,
       height=10)

##@ proportion of susceptibles
dfp <- filter(df_ad,fxet!="total")
lp <- split(dfp,dfp$fxet)
names(lp) <- paste0(LETTERS[1:length(lp)],") ",names(lp))
lp <- lapply(names(lp),\(li){
    dfi <- lp[[li]]
    f_plot_raw(dfi,
               facetlab = li,
               fillname="proportion of\nsusceptibles")
})
pfinal_rawprop <- f_wrap_plots(lp)
ggsave(filename="8_pastas_auxiliares/artigo_suscetiveis/figuras/map_raw_values/prop_susceptibles.jpeg",
       plot=pfinal_rawprop,
       dpi=200,
       width=14,
       height=10)

##@ pocket of susceptibility
vrange <- range(df_Bi$Bi)
lp <- split(df_Bi,df_Bi$fxet)
names(lp) <- paste0(LETTERS[1:length(lp)],") ",names(lp))
lp <- lapply(names(lp),\(li){
    dfi <- lp[[li]]
    f_plot_raw(dfi,
               facetlab=li,
               vrange = vrange,
               fillcol = "Bi",
               vtrans="log10",
               fillname="Pockets of\nsusceptibility")
})
pfinal_rawpockets <- f_wrap_plots(lp)
ggsave(filename="8_pastas_auxiliares/artigo_suscetiveis/figuras/map_raw_values/pockets_susceptibility.jpeg",
       plot=pfinal_rawpockets,
       dpi=200,
       width=14,
       height=10)


##@ Marginal impact of susceptibility
vrange <- range(df_Bi$Ii)
lp <- split(df_Bi,df_Bi$fxet)
names(lp) <- paste0(LETTERS[1:length(lp)],") ",names(lp))
lp <- lapply(names(lp),\(li){
    dfi <- lp[[li]]
    f_plot_raw(dfi,
               facetlab=li,
               vrange = vrange,
               fillcol = "Ii",
               fillname="Marginal\nsusceptible\ncontribution")
})
pfinal_rawpockets <- f_wrap_plots(lp)
ggsave(filename="8_pastas_auxiliares/artigo_suscetiveis/figuras/map_raw_values/marginal_susceptibility_impact.jpeg",
       plot=pfinal_rawpockets,
       dpi=200,
       width=14,
       height=10)


#--------------------------- heatmaps ------------------------------
#@ a ideia é plotar todos os heatmaps de global uni e bivariate Moran heatmap

vrangeI <- range(
    c(
        df_global_bsus$I,
        filter(df_global_psus,x!="total") %>% pull(I),
        filter(df_global_pop,x!="total") %>% pull(I),
        df_global_Isus$I
    )
)
f_heatmap <- \(dfi,rangeI=vrangeI){
    dfi %>%
        ggplot(aes(x=x,y=y,fill=I)) +
        geom_tile() +
        scale_x_discrete(expand = FALSE,position = "top") +
        scale_y_discrete(expand = FALSE) +
        geom_label(aes(label=label),
                   fill = "white",
                   size = 11/.pt) +
        scale_fill_viridis_c("Moran's I",option="magma",
                             limits = rangeI) +
        theme_classic() +
        labs(x="",y="") +
        theme(legend.position = "none",
              text = element_text(family="Times New Roman",
                                  size=12)
        )
}
f_add_title <- \(ggobj,vtitle){
    ggobj +
        labs(title = vtitle) +
        theme(plot.title = element_text(hjust = 0,
                                        face="bold",
                                        size=12,
                                        margin = margin(b = -5))
              )
}
lp <- list()
#@ overall population
lp$pop <- f_add_title(
    f_heatmap(filter(df_global_pop,x!="total")),
    vtitle="A) Overall population"
)

#@ overall population
lp$Psus <- f_add_title(
    f_heatmap(filter(df_global_psus,x!="total")),
    vtitle="B) Proportion of susceptibles"
)

#@ global uni e bivariate Moran heatmap
lp$Bsus <- f_add_title(
    f_heatmap(df_global_bsus),
    vtitle="C) Pockets of susceptibility"
)

#@ global uni e bivariate Moran heatmap
lp$Isus <- f_add_title(
    f_heatmap(df_global_Isus),
    vtitle = "D) Marginal susceptible contribution"
)

#@ gráfico composto
heatmap_final <- wrap_plots(lp,ncol=2)
ggsave(filename="8_pastas_auxiliares/artigo_suscetiveis/figuras/fig3.jpeg",
       plot = heatmap_final,
       dpi=120,
       width = 8,
       height = 8)

lf <- ls(pattern = "df_global")
df_global <- lapply(lf,\(li){
    dfi <- get(li)
    mutate(dfi,
           base=gsub("df_global_","",li) %>%
               gsub("bsus","bolsao",.) %>%
               gsub("Isus","MSI",.) %>%
               gsub("pop","pop",.) %>%
               gsub("psus","prop. suscep.",.))
}) %>% do.call("rbind",.) %>%
    filter(x!="total")

df_global %>%
    arrange(desc(I)) %>%
    group_by(base) %>%
    filter(I==max(I))

#----------------------------- LISA --------------------------------
f_plotLISA <- \(dfi,
                dfref=df_global_psus,
                legpos = c(0.825,0.85),
                textpos=c(Inf,-Inf),
                stripsize=11,
                vscalefill="LISA clusters",
                textsize=12){
    vanno1 <- "Global Moran’s I = "
    vanno2 <- "\ntwo-sided p-value < 0.05"
    vlabfac <- c('High-High','Low-Low','Low-High','High-Low','Non-significant')
    vfilt <- dfi$fxet[1]
    if(!grepl(":",vfilt)) vfilt <- paste(vfilt,vfilt,sep=" : ")
    vgmoran <- mutate(dfref,id=paste(x,y,sep = " : ")) %>%
        filter(id==vfilt) %>% pull(I)
    fxet <- dfi$fxet[1]
    vpal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "darkgray")
    names(vpal) <- vlabfac
    pret <- dfi %>%
        mutate(facet_label=fxet,
               class_agregado = factor(
                   class,
                   levels=c('High-High','Low-Low','Low-High','High-Low','Non-significant'),
                   labels=vlabfac)
        ) %>%
        ggplot() +
        geom_sf(aes(fill=class_agregado),
                linetype = "dashed",
                linewidth=0.25,
                color="white") +
        scale_fill_manual(vscalefill,values=vpal) +
        geom_sf(data=mutate(shp_gve,
                            facet_label=fxet),
                col="black",
                alpha=0.1,
                linewidth=1.2) +
        coord_sf(expand = FALSE) +
        theme_bw() +
        annotation_scale(location = "bl", width_hint = 0.5) +
        annotation_north_arrow(location = "bl",
                               which_north = "true",
                               pad_x = unit(0.75, "in"),
                               pad_y = unit(0.5, "in"),
                               style = north_arrow_fancy_orienteering) +
        annotate(
            geom = "label",
            x = textpos[1],
            y = textpos[2],
            label = paste0(vanno1,
                           round(vgmoran,3),
                           vanno2
            ),
            hjust = 1.05,
            vjust = -0.2,
            size = textsize,
            family = "Times New Roman",
            fill = alpha("white", 0.7),
            color = "black",
            label.padding = unit(0.2, "lines"),
            label.r = unit(0.1, "lines")
        ) +
        theme(legend.position = legpos,
              text=element_text(family="Times New Roman",size=12),
              legend.background = element_rect(fill = "white", color = "black"),
              legend.box.background = element_rect(fill = "white", color = "black"),
              legend.key = element_rect(fill = "white"),
              strip.text = element_text(size=stripsize,
                                        hjust=0,
                                        family="Times New Roman",
                                        face = "bold"),
              strip.background = element_rect(fill="white",color="white")) +
        facet_wrap(~facet_label)
    return(pret)
}

##@ overall pop
lp <- lapply(lLISA_pop,f_plotLISA)
names(lp) <- gsub(" : ","_",names(lp))
lapply(names(lp),\(li){
    vpath <- paste0("8_pastas_auxiliares/artigo_suscetiveis/figuras/LISA/pop/",li,".jpeg")
    ggsave(filename = vpath,
           plot=lp[[li]],
           dpi=120,
           width=12,
           height=8)
})

##@ psus
lp <- lapply(lLISA_psus,f_plotLISA)
names(lp) <- gsub(" : ","_",names(lp))
lapply(names(lp),\(li){
    vpath <- paste0("8_pastas_auxiliares/artigo_suscetiveis/figuras/LISA/psus/",li,".jpeg")
    ggsave(filename = vpath,
           plot=lp[[li]],
           dpi=120,
           width=12,
           height=8)
})

##@ bsus
lp <- lapply(lLISA_bsus,f_plotLISA,dfref=df_global_bsus)
names(lp) <- gsub(" : ","_",names(lp))
lapply(names(lp),\(li){
    vpath <- paste0("8_pastas_auxiliares/artigo_suscetiveis/figuras/LISA/bsus/",li,".jpeg")
    ggsave(filename = vpath,
           plot=lp[[li]],
           dpi=120,
           width=12,
           height=8)
})

##@ Isus
lp <- lapply(lLISA_Isus,f_plotLISA,dfref=df_global_Isus)
names(lp) <- gsub(" : ","_",names(lp))
lapply(names(lp),\(li){
    vpath <- paste0("8_pastas_auxiliares/artigo_suscetiveis/figuras/LISA/Isus/",li,".jpeg")
    ggsave(filename = vpath,
           plot=lp[[li]],
           dpi=120,
           width=12,
           height=8)
})


#------------------------- summary LISA PLOT -------------------------------
f_plot_summaryLISA <- \(dff,
                        vposleg=c(0.9,0.8),
                        stripsize=12){
    facetlab <- gsub("bsus","C) Pockets of susceptibility",dff$base[1]) %>%
        gsub("Isus","D) Marginal susceptible contribution",.) %>%
        gsub("psus","B) Proportion of susceptibles",.) %>%
        gsub("pop","A) Overall population",.)
    legendlab <- paste0("High–High\nLISA (%)","\nn = ",dff$nK[1])
    mutate(dff,
           facetlab=facetlab) %>%
        ggplot() +
        geom_sf(aes(fill=prop),
                linetype = "dashed",
                linewidth=0.25,
                color="darkgray") +
        scale_fill_distiller(legendlab,
                             # limits=c(0,100),
                             palette="Spectral") +
        geom_sf(data=mutate(shp_gve,
                            facetlab=facetlab),
                col="black",
                alpha=0.1,
                linewidth=0.5) +
        coord_sf(expand = FALSE) +
        facet_wrap(~facetlab) +
        theme_bw() +
        annotation_north_arrow(location = "bl",
                               which_north = "true",
                               pad_x = unit(0.75, "in"),
                               pad_y = unit(0.5, "in"),
                               style = north_arrow_fancy_orienteering) +
        theme(legend.position = "inside",
              text=element_text(family="Times New Roman",size=12),
              legend.position.inside = vposleg,
              legend.background = element_rect(fill = "white", color = "black"),
              legend.box.background = element_rect(fill = "white", color = "black"),
              legend.margin = margin(10, 10, 10, 15),
              legend.key = element_rect(fill = "white"),
              strip.text = element_text(face="bold",
                                        size=15,
                                        hjust=0,
                                        color="black",
                                        family = "Times New Roman"),
              strip.background = element_rect(fill="white",color="white")
        )
}
lp_summaryLISA <- lapply(l_summaryLISA[c("pop","psus","bsus","Isus")],f_plot_summaryLISA)
psummaryLISA <- wrap_plots(lp_summaryLISA,nrow=2)
ggsave(filename = "8_pastas_auxiliares/artigo_suscetiveis/figuras/fig4.jpeg",
       plot=psummaryLISA,
       dpi=200,
       width = 16,
       height = 14)
# f_img("8_pastas_auxiliares/artigo_suscetiveis/figuras/summaryLISA.jpeg")


#------------------- contiguidade ------------#
# dff <- l_graphHH[["psus"]]
f_plotAgregados <- \(dff,
                     vposleg=c(0.9,0.8),
                     stripsize=9){
    facetlab <- gsub("bsus","C) Pockets of susceptibility",dff$base[1]) %>%
        gsub("Isus","D) Marginal susceptible contribution",.) %>%
        gsub("psus","B) Proportion of susceptibles",.) %>%
        gsub("pop","A) Overall population",.)
    dfad <- st_drop_geometry(dff) %>% count(comp_conectado) %>% filter(!is.na(comp_conectado))
    dfp <- left_join(dff,dfad)
    dfp %>%
        mutate(facetlab=facetlab) %>%
        ggplot() +
        geom_sf(aes(fill=n),
                linetype = "dashed",
                linewidth=0.25,
                color="gray") +
        scale_fill_gradient("component size\n(number of mun.)",
                            low = "blue", high = "red") +
        geom_sf(data=mutate(shp_gve,
                            facetlab=facetlab),
                col="black",
                alpha=0.1,
                linewidth=0.5) +
        coord_sf(expand = FALSE) +
        facet_wrap(~facetlab) +
        theme_bw() +
        annotation_north_arrow(location = "bl",
                               which_north = "true",
                               pad_x = unit(0.75, "in"),
                               pad_y = unit(0.5, "in"),
                               style = north_arrow_fancy_orienteering) +
        theme(legend.position = "inside",
              text = element_text(family="Times New Roman",size=9),
              legend.position.inside = vposleg,
              legend.background = element_rect(fill = "white", color = "black"),
              legend.box.background = element_rect(fill = "white", color = "black"),
              legend.margin = margin(10, 10, 10, 15),
              legend.key = element_rect(fill = "white"),
              strip.text = element_text(size=15,
                                        hjust=0,
                                        family="Times New Roman",
                                        face = "bold"),
              strip.background = element_rect(fill="white",color="white")
              )
}
lp_cont <- lapply(l_graphHH[c("pop","psus","bsus","Isus")],f_plotAgregados,vposleg=c(0.85,0.8))
p_compsize <- wrap_plots(lp_cont,ncol=2)
ggsave(filename = "8_pastas_auxiliares/artigo_suscetiveis/figuras/figSI2.jpeg",
       plot=p_compsize,
       dpi=200,
       width = 16,
       height = 14)
# f_img("8_pastas_auxiliares/artigo_suscetiveis/figuras/componenteSizeLISA.jpeg")


#----------------------------- Graf. Exp --------------------------------
#@ boxplot das 4 variáveis de interesse:
f_boxplot <- \(dfi,ycol="POP",vtrans="log10",facetlab="A) Overall population (log scale)"){
    dfi %>%
        mutate(flab=facetlab) %>%
        ggplot(aes(x=fxet,y=.data[[ycol]])) +
        geom_boxplot(color="red") +
        geom_jitter(alpha=0.25) +
        scale_y_continuous(trans = vtrans) +
        facet_wrap(~flab) +
        labs(x="age group (years)",y="") +
        theme_bw() +
        theme(strip.text = element_text(face="bold",
                                        size=15,hjust=0,color="black",
                                        family="Times New Roman"),
              text=element_text(family="Times New Roman",size=14),
              strip.background = element_rect(fill="white",color="white")
              )
}
lp <- list()
lp$overallpop <- f_boxplot(filter(df_ad,fxet!="total"))
lp$propsus <- f_boxplot(filter(df_ad,fxet!="total"),
                        ycol="p_sus",
                        vtrans="identity",
                        facetlab="B) Proportion of susceptibles")
lp$pockets <- f_boxplot(df_Bi,
                        ycol="Bi",
                        vtrans="log10",
                        facetlab="C) Pocket of susceptibility (log scale)")
lp$MSI <- f_boxplot(df_Bi,
                    ycol="Ii",
                    vtrans="identity",
                    facetlab="D) Marginal susceptible contribution")
boxplot_final <- wrap_plots(lp,ncol=2) +
    plot_layout(axes = "collect", guides = "collect") &
    theme(
        plot.margin = margin(5, 5, 5, 5),
        axis.title = element_text(size = 9)
    )
ggsave(filename="8_pastas_auxiliares/artigo_suscetiveis/figuras/fig1.jpeg",
       plot=boxplot_final,
       dpi=120,
       width=10,
       height=10)

#@ heatmap triangular da correlaçaõ entre fxet por variável de interesse:
df_fxet
v2fx <- df_fxet[1,]
f_corr_byvar <- \(v2fx,dfref=df_ad,cname="pop_z"){
    #
    dfret <- v2fx
    v2fx <- unlist(v2fx)
    vcols <- c("ID_MUN",cname)
    #
    fxet1 <- filter(dfref,fxet==v2fx[1]) %>%
        st_drop_geometry() %>%
        select(all_of(vcols))
    fxet2 <- filter(dfref,fxet==v2fx[2]) %>%
        st_drop_geometry() %>%
        select(all_of(vcols))
    #
    dfr <- inner_join(x=fxet1,y=fxet2,by=c("ID_MUN"))
    #
    dfret$cor <- cor(x=dfr[[paste0(cname,".x")]],
                     y=dfr[[paste0(cname,".y")]],
                     method="pearson")
    dfret$var <- cname
    #
    return(dfret)
}

vlevfx1 <- vlevfx1[vlevfx1!="total"]
vlevfx <- vlevfx1[length(vlevfx1):1]
df_cor <- list(
    split(df_fxet,seq(nrow(df_fxet))) %>%
        lapply(f_corr_byvar) %>%
        do.call("rbind",.),
    split(df_fxet,seq(nrow(df_fxet))) %>%
        lapply(f_corr_byvar,cname="p_sus_z") %>%
        do.call("rbind",.),
    split(df_fxet,seq(nrow(df_fxet))) %>%
        lapply(f_corr_byvar,dfref=df_Bi,cname="Bi_z") %>%
        do.call("rbind",.),
    split(df_fxet,seq(nrow(df_fxet))) %>%
        lapply(f_corr_byvar,dfref=df_Bi,cname="Ii_z") %>%
        do.call("rbind",.)
) %>% do.call("rbind",.) %>%
    rename("x"="V1","y"="V2") %>%
    mutate(x=factor(x,levels=vlevfx1),
           y=factor(y,levels=vlevfx),
           label=round(cor,3),
           var=factor(var,
                      levels=c("pop_z",
                               "p_sus_z",
                               "Bi_z",
                               "Ii_z"),
                      labels=c("A) Overall population",
                               "B) Proportion of susceptibles",
                               "C) Pocket of susceptibility",
                               "D) Marginal susceptible contribution"))
    )

p <- df_cor %>%
    ggplot(aes(x=x,y=y,fill=cor)) +
    geom_tile() +
    scale_x_discrete(expand = FALSE,position = "top") +
    scale_y_discrete(expand = FALSE) +
    geom_label(aes(label=label),fill = "white",
               size = 10/.pt) +
    scale_fill_viridis_c("Pearson\ncorrelation") +
    theme_classic() +
    labs(x="",y="") +
    facet_wrap(~var,ncol=2) +
    theme(legend.position = "none",
          text=element_text(family="Times New Roman",size=10),
          strip.text = element_text(face="bold",size=12,hjust=0,
                                    family="Times New Roman"),
          strip.placement = "outside",
          strip.background = element_blank(),
          plot.title = element_text(margin = margin(b = -10))
    )
ggsave(filename="8_pastas_auxiliares/artigo_suscetiveis/figuras/fig2.jpeg",
       plot = p,
       width = 6.5,
       height = 7.5,
       dpi=120)

# df_Bi %>%
#     pivot_longer(cols=starts_with("Bi")) %>%
#     filter(name!="Bi") %>%
#     ggplot(aes(y=value,x=name)) +
#     geom_boxplot() +
#     geom_jitter() +
#     facet_wrap(~fxet,nrow=2)
# df_Bi %>%
#     ggplot(aes(x=Bi,y=Ii)) +
#     geom_point() +
#     scale_x_log10() +
#     facet_wrap(~fxet,nrow=2,scales="free")



#---------------------- tabelas de interesse ----------------------------
#@ quais são os casos em que o Moran Global foi classificado como relevante (pvalor<0.05 e I>=0.1)?
ldf <- lapply(ls(pattern = "sig"),get)
names(ldf) <- gsub("df_sig_","",ls(pattern = "sig"))
ldf <- ldf[names(ldf)!="pop"]
dftofx <- lapply(names(ldf),\(li){
    dfi <- ldf[[li]]
    dfi$label <- NULL
    dfi$base <- gsub("bsus","Pockets of susceptibility",li) %>%
        gsub("Isus","Marginal susceptible contribution",.) %>%
        gsub("psus","Proportion of susceptibles",.)
    return(dfi)
}) %>% do.call("rbind",.)
rownames(dftofx) <- NULL
names(dftofx) <- gsub("base","metric",names(dftofx))
dftofx$metric <- factor(dftofx$metric,
                        levels=c("Proportion of susceptibles",
                                 "Pockets of susceptibility",
                                 "Marginal susceptible contribution"))
dftofx <- arrange(dftofx,metric)

library(flextable)
library(officer)

vcap <- as_paragraph(
    as_chunk("Table 1. ", props = fp_text(bold = TRUE)), # Apenas o Table 1. em negrito
    as_chunk("Significant instances where Moran's I ≥ 0.100 and p-value < 0.05. P-values were derived from a bootstrap procedure with 10,001 simulations.")
)

f_flextable <- \(dfi=dftofx){
    #
    df_grouped <- as_grouped_data(dfi, groups = "metric")
    #
    as_flextable(df_grouped) %>%
        # Estética geral
        bg(bg = "white", part = "all") %>%
        theme_vanilla() %>%
        #
        flextable::set_header_labels(
            values = list(
                x = "age group 1",
                y = "age group 2",
                I = "Moran's I",
                pvalue = "p-value"
            )
        ) %>%
        # Destacar as linhas que são os grupos
        bold(i = ~ !is.na(metric), bold = TRUE, part = "body") %>%
        bg(i = ~ !is.na(metric), bg = "#F7F7F7", part = "body") %>% #
        # Formatação de números
        colformat_double(digits = 3) %>%
        #colformat_double(j = "weight", digits = 2) %>% #

        # Alinhamento e legenda
        align(align = "center", part = "all") %>%
        align(j = 1, align = "left", part = "body") %>% #
        # set_caption(caption = vcap) %>%
        add_header_lines(values = as_paragraph(
            as_chunk("Table SI1. ", props = fp_text(bold = TRUE)),
            as_chunk("Significant global spatial autocorrelation for susceptibility metrics (10,001 bootstrap simulations).")
        )) %>%
        #
        set_table_properties(layout = "fixed") %>%
        autofit() # Ajusta as colunas ao conteúdo
    # Define a largura apenas para a coluna 'modelo'
    # flextable::width(j = "modelo", width = 9.5/5, unit = "in") %>%
    # # Define a largura para todas as OUTRAS colunas
    # flextable::width(j = c("dAICc", "df", "weight", "dev.expl"), width = 6.5/5, unit = "in")
}
dffx <- f_flextable()
save_as_image(dffx, path = "8_pastas_auxiliares/artigo_suscetiveis/figuras/tabela1.png")
save_as_docx(
    dffx,
    path = "8_pastas_auxiliares/artigo_suscetiveis/figuras/tabela1.docx"
)


#------------------ municípios classificados como HH -------------------#
dff <- l_summaryLISA[[1]]
f_tabHH <- \(dff){
    dfr <- st_drop_geometry(dff)
    dfr$base <- gsub("bsus","Pockets of susceptibility",dfr$base) %>%
        gsub("Isus","Marginal susceptible contribution",.) %>%
        gsub("psus","Proportion of susceptibles",.) %>%
        gsub("pop","Overall population",.)
    filter(dfr,nHH>0) %>%select(-prop)
}
df_tabHH <- lapply(l_summaryLISA,f_tabHH) %>% do.call("rbind",.)
rownames(df_tabHH) <- NULL
df_tabHH$base <- factor(df_tabHH$base,
                        levels=c("Overall population",
                                 "Proportion of susceptibles",
                                 "Pockets of susceptibility",
                                 "Marginal susceptible contribution")
)
df_tabHH %>%
    group_by(base) %>%
    count()

#---------------- quais os municipios em relação aos agregados estruturais ----#
l_count_graphHH <- lapply(l_graphHH,\(li) li %>%
                              st_drop_geometry %>%
                              count(comp_conectado) %>%
                              filter(!is.na(comp_conectado))
)

df_ret <- lapply(l_graphHH,\(li) li %>%
                     st_drop_geometry %>%
                     filter(!is.na(comp_conectado))) %>%
    do.call("rbind",.) %>%
    rename(id_component = comp_conectado) %>%
    select(-prop) %>%
    mutate(base = gsub("bsus","Pockets of susceptibility",base) %>%
               gsub("Isus","Marginal susceptible contribution",.) %>%
               gsub("psus","Proportion of susceptibles",.) %>%
               gsub("pop","Overall population",.),
           base = factor(base,
                         levels=c("Overall population",
                                  "Proportion of susceptibles",
                                  "Pockets of susceptibility",
                                  "Marginal susceptible contribution"))
    )


#@ os componentes estão contidos em um único GVE?
f_gveComp <- \(dff){
    #
    vbase <- gsub("bsus","Pockets of susceptibility",dff$base[1]) %>%
        gsub("Isus","Marginal susceptible contribution",.) %>%
        gsub("psus","Proportion of susceptibles",.) %>%
        gsub("pop","Overall population",.)
    dfres <- inner_join(dff,dfgeo)
    dfres <- st_drop_geometry(dfres) %>%
        filter(!is.na(comp_conectado)) %>%
        select(comp_conectado,ID_MUN,GVE) %>%
        distinct() %>%
        arrange(comp_conectado,ID_MUN) %>%
        rename(id_comp=comp_conectado)
    dfsum <- group_by(dfres,id_comp) %>%
        summarise(nMUN=length(unique(ID_MUN)),
                  nGVE=length(unique(GVE)),
                  namesGVE=paste(unique(GVE),collapse=", "),
                  idsMUN=paste(unique(ID_MUN),collapse=", "),
                  base=vbase
        )
    #list("idMUN"=dfres,"namesGVE"=dfsum)
    return(dfsum)
}
df_comp <- lapply(l_graphHH,f_gveComp) %>% do.call("rbind",.) %>%
    mutate(base = factor(base,
                         levels=c("Overall population",
                                  "Proportion of susceptibles",
                                  "Pockets of susceptibility",
                                  "Marginal susceptible contribution"))
    )

f_flextable <- \(dfi=df_comp){
    #
    df_grouped <- rename(dfi,metric=base) %>%
        arrange(metric) %>%
        as_grouped_data(groups = "metric")
    #
    as_flextable(df_grouped) %>%
        # Estética geral
        bg(bg = "white", part = "all") %>%
        theme_vanilla() %>%
        #
        flextable::set_header_labels(
            values = list(
                id_comp = "comp. ID",
                nMUN = "n mun.",
                nGVE = "n GVE",
                namesGVE = "GVE",
                idsMUN = "id mun"
            )
        ) %>%
        # Destacar as linhas que são os grupos
        bold(i = ~ !is.na(metric), bold = TRUE, part = "body") %>%
        bg(i = ~ !is.na(metric), bg = "#F7F7F7", part = "body") %>% #
        #colformat_double(j = "weight", digits = 2) %>% #

        # Alinhamento e legenda
        align(align = "center", part = "all") %>%
        align(j = 1, align = "left", part = "body") %>% #
        #set_caption(caption = vcap) %>%
        add_header_lines(values = as_paragraph(
            as_chunk("Table SI2. ", props = fp_text(bold = TRUE)),
            as_chunk("Summary of spatial components by health surveillance groups. The table presents the component ID (Comp. ID), the number of municipalities per component (n Mun.), and the number of affected Epidemiological surveillance regions (n GVE) with their respective names (GVE), categorized by susceptibility metric.")
        )) %>%
        #
        set_table_properties(layout = "fixed") %>%
        autofit() %>%
        #
        width(j = "namesGVE", width = 4) %>%
        width(j = "idsMUN", width = 3)
}
dffx <- f_flextable()
save_as_image(dffx, path = "8_pastas_auxiliares/artigo_suscetiveis/figuras/tabelaSI2.png")
save_as_docx(
    dffx,
    path = "8_pastas_auxiliares/artigo_suscetiveis/figuras/tabelaSI2.docx",
    pr_section = prop_section(
        page_size = page_size(orient = "landscape"),
        type = "continuous"
    )
)


#-------------------------- arrumar todas as figuras
lf <- list.files(path="8_pastas_auxiliares/artigo_suscetiveis/figuras",
                 pattern=".jpeg",
                 full.names=TRUE,
                 recursive=TRUE)
lapply(lf,f_img)
