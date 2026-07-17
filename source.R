##################################################################################
# pacotes
library(patchwork)
library(scales)
library(future.apply)
library(ggspatial)
library(spdep)
library(tmap)
library(magick)
library(tidyverse)
library(geobr)
library(sf)
library(rsatscan)
library(ggthemes)
library(ggpmisc)
library(igraph)
f_z <- \(vx){
    (vx - mean(vx)) / sd(vx)

}
f_img <- \(imgpath){
    vimg <- image_read(imgpath) %>%
        image_trim()
    image_write(vimg,imgpath)
}
plan(multisession, workers = 10)
#---------------------------------------------------------------------
#---------------------------- Dados ----------------------------------
#---------------------------------------------------------------------

#------------- dados espaciais ------------#
# shapefile GVE
oname <- load("raw_data/shp_gve.RData")
shp_gve <- get(oname);rm(list=oname);gc()
shp_gve <- shp_gve %>%
    mutate(z=0)
## dados dos municípios por GVE ##
oname <- load("raw_data/geo_full.RData")
dfgeo <- get(oname);rm(list=oname);gc()
dfgeo <- select(dfgeo,ID_MUN,MUN,GVE)
## shapefiles dos municípios ##
dfmun0 <- read_municipality(code="SP",year=2024)
dfmun0 <- dfmun0 %>%
    mutate(ID_MUN = substr(code_muni, 1, 6) %>% as.character()) %>%
    select(ID_MUN,name_muni)
dfmun <- select(dfmun0,-name_muni)
##@ objetos espaciais @##
dfnb <- right_join(x=dfmun,y=dfgeo)

#------------- dados suscetibilidade ------------#
##@ população total e proporção de suscetíveis ##@
df_ad <- readRDS(file="7_produtos/agregados_espaciais/df_suscetiveis.rds") %>%
    mutate(fxet = gsub(" anos","",fx_et) %>% gsub(" a ","-",.)) %>%
    select(-c(fx_et:n_D2))
df_ad <- rbind(
    df_ad,
    group_by(df_ad,id_mn) %>%
        summarise(n_imn=sum(n_imn),POP=sum(POP)) %>%
        mutate(p_sus=100*(POP-n_imn)/POP,
               p_sus=ifelse(p_sus<0,0,p_sus),
               fxet="total")
)
df_ad <- group_by(df_ad,fxet) %>%
    mutate(p_sus_z = f_z(p_sus),
           pop_z = f_z(log(POP))) %>%
    right_join(x=dfnb,y=.,by=c("ID_MUN"="id_mn"))

##@ bolsão de suscetibilidade @##
f_arrumar_nome <- \(vname){
    #
    ln <- str_split_1(vname,"-") %>% as.numeric
    vmax <- max(ln)-1
    #
    lret <- c(min(ln),vmax) %>% as.character
    lret <- sapply(lret, \(vi){
        if(nchar(vi)==1){
            paste0("0",vi)
        }else{
            vi
        }
    })
    #
    paste(lret,collapse="-")
}
#@1) matriz de contatos
df_contact <- read.csv("4_bases/finalsize/contact_matrix.csv")
vnames <- gsub("_","-",df_contact$X) %>%
    gsub("00","01",.) %>%
    sapply(.,f_arrumar_nome)
names(vnames) <- df_contact$X
rownames(df_contact) <- vnames
df_contact <- select(df_contact,-X)
#- seleção dos elementos de interesse:
vfxet <- unique(df_ad$fxet)
vfxet <- vfxet[vfxet!="total"]
vnames <- vnames[vnames %in% vfxet]
df_contato <- df_contact[rownames(df_contact)%in%vnames,
                         paste0("X",names(vnames))]
names(df_contato) <- vnames
df_contato <- df_contato %>%
    tibble::rownames_to_column(var = "x") %>%
    pivot_longer(
        cols = -x,
        names_to = "y",
        values_to = "contact_rate"
    )

#@2) calculo de Bij
f_df_p_bolsao <- \(dfi){
    # setup
    vfxet <- dfi[,c("x","y")] %>% unlist
    dfx1 <- filter(df_ad,fxet==vfxet[1]) %>%
        st_drop_geometry() %>%
        select(ID_MUN,n_imn,POP) %>%
        rename(Si=n_imn,Ni=POP) %>%
        mutate(Si = ifelse(Si>Ni,Ni,Si))
    dfx2 <- filter(df_ad,fxet==vfxet[2]) %>%
        st_drop_geometry() %>%
        select(ID_MUN,n_imn,POP) %>%
        rename(Sj=n_imn,Nj=POP) %>%
        mutate(Sj = ifelse(Sj>Nj,Nj,Sj))
    # juntar
    dfr <- inner_join(dfx1,dfx2)
    # padronizar e retornar
    mutate(dfr,
           Cij = dfi$contact_rate,
           Bij = Cij * Si * Sj / Ni,
           Bij_z = f_z(log(Bij)),
           fxet = paste(unique(vfxet),collapse = ":")) %>%
        relocate(fxet)
}
df_Bij <- split(df_contato,seq(nrow(df_contato))) %>%
    future_lapply(.,f_df_p_bolsao) %>%
    do.call("rbind",.)
write_csv(df_Bij,file = "8_pastas_auxiliares/artigo_suscetiveis/csv/df_Bij.csv")

#@3) Calculo de Bi
df_Bi <- lapply(vnames,\(li){
    dftosum <- filter(df_Bij,grepl(paste0("^",li),fxet))
    dfret <- group_by(dftosum,ID_MUN) %>%
        summarise(Bi=sum(Bij),
                  Si=unique(Si)) %>%
        mutate(fxet=li,
               Ii=Bi/Si,
               Ii_z=f_z(Ii),
               Ii_logz=f_z(log(Ii)),
               Bi_z=f_z(log(Bi))) # por se tratar de uma variável multiplicativa se aplica o log
    relocate(dfret,fxet)
}) %>% do.call("rbind",.)
df_Bi <- right_join(dfnb,df_Bi)
write_csv(df_Bi,file = "8_pastas_auxiliares/artigo_suscetiveis/csv/df_Bi.csv")

#--------- tabela guia para as análises -------#
f_fxet_comb <- \(vfxet=unique(df_ad$fxet)){
    allcomb <- combn(vfxet,2)
    as.data.frame(t(allcomb))
}
vfxet0 <- unique(df_ad$fxet)
vfxet0 <- vfxet0[vfxet0!="total"]
df_fxet <- f_fxet_comb(vfxet = vfxet0)

#---------------------------------------------------------------------
#--------------------- Global autocorrelation ------------------------
#---------------------------------------------------------------------

#@ f_unimoran: Moran univariado
f_unimoran <- \(v1fx,dfref=df_ad,colname="p_sus_z"){
    #@ v1fx: 1 unica fxet
    df_moran <- filter(dfref,fxet==v1fx)
    # vizinhança espacial
    nb <- poly2nb(df_moran, queen=TRUE)
    l_MNIsolado <- which(card(nb) == 0)
    for(i in l_MNIsolado){
        dist <- st_distance(dfnb[i,],dfnb)
        vizinho <- order(dist)[2] #
        #
        nb[[i]] <- vizinho
        #
        nb[[vizinho]] <- c(nb[[vizinho]],i)
    }
    nbw <- nb2listw(nb, style = "W")
    #
    gmoran <- moran.mc(df_moran[[colname]],
                       nbw,
                       nsim = 10000,
                       alternative = "two.sided")
    #
    dfr <- data.frame(x=v1fx,
                      y=v1fx,
                      I=gmoran$statistic,
                      pvalue=gmoran$p.value)
    #
    list(
        "data"=dfr,
        "uniMoran"=gmoran
    )
}
f_pseudo_pvalue <- \(obj_moranbv){
    I_obs <- obj_moranbv$t0
    I_sim <- obj_moranbv$t
    (sum(abs(I_sim) >= abs(I_obs)) + 1) / (length(I_sim) + 1)
}
f_bvmoran <- \(v2fx,dfref=df_ad,colname="p_sus_z"){
    #
    v2fx <- unlist(v2fx)
    vfxet <- paste0(c("x","y"),v2fx)
    # padronização
    vcols <- c("ID_MUN",colname)
    fxet1 <- filter(dfref,fxet==v2fx[1]) %>% select(all_of(vcols))
    fxet2 <- filter(dfref,fxet==v2fx[2]) %>% select(all_of(vcols)) %>% st_drop_geometry()
    dfr <- inner_join(fxet1,fxet2,by=c("ID_MUN"))
    names(dfr)[grepl("\\.x|\\.y",names(dfr))] <- vfxet
    #
    nb <- poly2nb(dfr, queen=TRUE)
    l_MNIsolado <- which(card(nb) == 0)
    for(i in l_MNIsolado){
        dist <- st_distance(dfnb[i,],dfnb)
        vizinho <- order(dist)[2] #
        #
        nb[[i]] <- vizinho
        #
        nb[[vizinho]] <- c(nb[[vizinho]],i)
    }
    nbw <- nb2listw(nb, style = "W")
    # padronização
    fxet1 <- dfr[[vfxet[[1]]]]
    fxet2 <- dfr[[vfxet[[2]]]]
    #
    l_bvm <- list()
    # Moran bivariado 12
    l_bvm$bv12 <- moran_bv(x=fxet1,y=fxet2,listw=nbw,nsim = 10000)
    l_bvm$bv12$pseudo_pvalue <- f_pseudo_pvalue(l_bvm$bv12)
    # Moran bivariado 21
    l_bvm$bv21 <- moran_bv(x=fxet2,y=fxet1,listw=nbw,nsim = 10000)
    l_bvm$bv21$pseudo_pvalue <- f_pseudo_pvalue(l_bvm$bv21)
    # Correlation Pearson
    vPearsonCor <- cor(x=fxet1,y=fxet2)
    #
    dfr <- data.frame(x=v2fx,
                      y=v2fx[2:1],
                      I=sapply(l_bvm,"[[","t0"),
                      pvalue=sapply(l_bvm,"[[","pseudo_pvalue"),
                      row.names = NULL)
    #
    list(
        "data"=dfr,
        "Pcorrelation"=vPearsonCor,
        "biMoran"=l_bvm
    )
}

#-------- população por faixa etária ----------------
##@ univariado
l_unimoran <- future_lapply(unique(df_ad$fxet),f_unimoran,colname="pop_z",future.seed = TRUE)
##@ bivariado
l_bvmoran <- split(df_fxet, seq(nrow(df_fxet))) %>%
    future_lapply(., f_bvmoran, colname="pop_z", future.seed = TRUE)
names(l_bvmoran) <- mutate(df_fxet,
                           vname = paste(V1,V2,sep=" : ")) %>%
    pull(vname)
##@ Sumário: matriz do I de Moran
vlevfx1 <- unique(df_ad$fxet)
vlevfx <- vlevfx1[length(vlevfx1):1]
df_global_pop <- rbind(
    lapply(l_bvmoran, "[[","data") %>% do.call("rbind",.),
    lapply(l_unimoran, "[[","data") %>% do.call("rbind",.)
) %>%
    mutate(x=factor(x,levels=vlevfx1),
           y=factor(y,levels=vlevfx),
           label=paste0(
               round(I,3),
               ifelse(pvalue<0.05,"*","")
           )
    )
write_csv(df_global_pop,file = "8_pastas_auxiliares/artigo_suscetiveis/csv/globalMoran_pop.csv")


#-------- proporção de suscetiveis ----------------
##@ univariado
l_unimoran <- future_lapply(unique(df_ad$fxet),f_unimoran,future.seed = TRUE)
##@ bivariado
l_bvmoran <- split(df_fxet, seq(nrow(df_fxet))) %>%
    future_lapply(., f_bvmoran, future.seed = TRUE)
names(l_bvmoran) <- mutate(df_fxet,
                           vname = paste(V1,V2,sep=" : ")) %>%
    pull(vname)
##@ Sumário: matriz do I de Moran
vlevfx1 <- unique(df_ad$fxet)
vlevfx <- vlevfx1[length(vlevfx1):1]
df_global_psus <- rbind(
    lapply(l_bvmoran, "[[","data") %>% do.call("rbind",.),
    lapply(l_unimoran, "[[","data") %>% do.call("rbind",.)
) %>%
    mutate(x=factor(x,levels=vlevfx1),
           y=factor(y,levels=vlevfx),
           label=paste0(
               round(I,3),
               ifelse(pvalue<0.05,"*","")
           )
    )
write_csv(df_global_psus,file = "8_pastas_auxiliares/artigo_suscetiveis/csv/globalMoran_psus.csv")

#-------- bolsão de suscetibilidade ----------------
##@ univariado
l_unimoran <- future_lapply(unique(df_Bi$fxet),
                            f_unimoran,
                            dfref=df_Bi,
                            colname="Bi_z",
                            future.seed=TRUE)
##@ bivariado
l_bvmoran <- split(df_fxet, seq(nrow(df_fxet))) %>%
    future_lapply(.,f_bvmoran,
                  dfref=df_Bi,
                  colname="Bi_z",
                  future.seed = TRUE)
names(l_bvmoran) <- mutate(df_fxet,
                           vname = paste(V1,V2,sep=" : ")) %>%
    pull(vname)
##@ Sumário
df_global_bsus <- rbind(
    lapply(l_bvmoran, "[[","data") %>% do.call("rbind",.),
    lapply(l_unimoran, "[[","data") %>% do.call("rbind",.)
) %>%
    mutate(x=factor(x,levels=vlevfx1),
           y=factor(y,levels=vlevfx),
           label=paste0(
               round(I,3),
               ifelse(pvalue<0.05,"*","")
           )
    )
write_csv(df_global_bsus,file = "8_pastas_auxiliares/artigo_suscetiveis/csv/globalMoran_bsus.csv")


#-------- Impacto marginal na suscetibilidade ----------------
##@ univariado
l_unimoran <- future_lapply(unique(df_Bi$fxet),
                            f_unimoran,
                            dfref=df_Bi,
                            colname="Ii_z",
                            future.seed=TRUE)
##@ bivariado
l_bvmoran <- split(df_fxet, seq(nrow(df_fxet))) %>%
    future_lapply(.,f_bvmoran,
                  dfref=df_Bi,
                  colname="Ii_z",
                  future.seed = TRUE)
names(l_bvmoran) <- mutate(df_fxet,
                           vname = paste(V1,V2,sep=" : ")) %>%
    pull(vname)
##@ Sumário
df_global_Isus <- rbind(
    lapply(l_bvmoran, "[[","data") %>% do.call("rbind",.),
    lapply(l_unimoran, "[[","data") %>% do.call("rbind",.)
) %>%
    mutate(x=factor(x,levels=vlevfx1),
           y=factor(y,levels=vlevfx),
           label=paste0(
               round(I,3),
               ifelse(pvalue<0.05,"*","")
           )
    )
write_csv(df_global_Isus,file = "8_pastas_auxiliares/artigo_suscetiveis/csv/globalMoran_Isus.csv")


#---------------------------------------------------------------------
#----------------------------- LISA  --------------------------------
#---------------------------------------------------------------------
#@
f_class_lbvm <- \(vfac){
    vret <- vfac
    factor(vret,levels=c(levels(vret),"Non-significant"))
}
#@ f_unimoran: Moran univariado
f_unimoranLISA <- \(v1fx,dfref=df_ad,colname="p_sus_z"){
    #@ v1fx: 1 unica fxet
    df_moran <- filter(dfref,fxet==v1fx)
    # vizinhança espacial
    nb <- poly2nb(df_moran, queen=TRUE)
    l_MNIsolado <- which(card(nb) == 0)
    for(i in l_MNIsolado){
        dist <- st_distance(dfnb[i,],dfnb)
        vizinho <- order(dist)[2] #
        #
        nb[[i]] <- vizinho
        #
        nb[[vizinho]] <- c(nb[[vizinho]],i)
    }
    nbw <- nb2listw(nb, style = "W")
    #
    lmoran <- localmoran_perm(df_moran[[colname]],
                              nbw,
                              nsim = 10000,
                              alternative = "two.sided")
    #
    dfr <- data.frame(fxet = v1fx,
                      I = lmoran[,1],
                      pvalue = lmoran[,5],
                      class = f_class_lbvm(attr(lmoran,"quadr")$mean))
    dfr$class[dfr$pvalue>=0.05] <- "Non-significant"
    #
    dfr <- cbind(dfr,df_moran[,c("ID_MUN","geometry")])
    st_as_sf(dfr)
}
f_bvmoranLISA <- \(v2fx,dfref=df_ad,colname="p_sus_z"){
    #
    v2fx <- unlist(v2fx)
    vfxet <- paste0(c("x","y"),v2fx)
    # padronização
    vcols <- c("ID_MUN",colname)
    fxet1 <- filter(dfref,fxet==v2fx[1]) %>% select(all_of(vcols))
    fxet2 <- filter(dfref,fxet==v2fx[2]) %>% select(all_of(vcols)) %>% st_drop_geometry()
    dfr <- inner_join(fxet1,fxet2,by=c("ID_MUN"))
    names(dfr)[grepl("\\.x|\\.y",names(dfr))] <- vfxet
    #
    nb <- poly2nb(dfr, queen=TRUE)
    l_MNIsolado <- which(card(nb) == 0)
    for(i in l_MNIsolado){
        dist <- st_distance(dfnb[i,],dfnb)
        vizinho <- order(dist)[2] #
        #
        nb[[i]] <- vizinho
        #
        nb[[vizinho]] <- c(nb[[vizinho]],i)
    }
    nbw <- nb2listw(nb, style = "W")
    # padronização
    fxet1 <- dfr[[vfxet[[1]]]]
    fxet2 <- dfr[[vfxet[[2]]]]
    #
    l_bvm <- list()
    # Moran bivariado 12
    l_bvm$bv12 <- localmoran_bv(x=dfr[[vfxet[1]]],
                                y=dfr[[vfxet[2]]],
                                listw=nbw,
                                nsim = 10000)
    dfr$Ixy <- l_bvm$bv12[,1]
    dfr$pvalue_xy <- l_bvm$bv12[,5]
    dfr$class_xy <- f_class_lbvm(attr(l_bvm$bv12,"quadr")$mean)
    dfr$class_xy[dfr$pvalue_xy>=0.05] <- "Non-significant"
    # Moran bivariado 21
    l_bvm$bv21 <- localmoran_bv(x=dfr[[vfxet[2]]],
                                y=dfr[[vfxet[1]]],
                                listw=nbw,
                                nsim = 10000)
    dfr$Iyx <- l_bvm$bv21[,1]
    dfr$pvalue_yx <- l_bvm$bv21[,5]
    dfr$class_yx <- f_class_lbvm(attr(l_bvm$bv21,"quadr")$mean)
    dfr$class_yx[dfr$pvalue_yx>=0.05] <- "Non-significant"
    #
    # relocate(dfr,geom,.after = last_col())
    return(dfr)
}
f_apply_localM <- \(dfim,dfref=df_ad,colname="p_sus_z"){
    vfxet <- dfim[,c("x","y")] %>% unlist %>% unique
    if(length(vfxet)==1){
        dfret <- f_unimoranLISA(vfxet,dfref,colname)
    }else{
        dfret <- f_bvmoranLISA(vfxet,dfref,colname)
        vnames <- grep(paste(vfxet,collapse="|"),
                       names(dfret),value = TRUE) %>%
            str_extract(.,"[a-zA-Z]+") %>%
            paste(collapse="")
        dfret <- select(dfret,ID_MUN,contains(vnames)) %>%
            mutate(fxet = paste(vfxet,collapse=" : ")) %>%
            relocate(fxet)
        names(dfret) <-
            gsub("_xy","",names(dfret)) %>%
            gsub("xy","",.)
    }
    return(dfret)
}

#------------- população por faixa etária ---------
##@ apenas os casos significativos
df_sig_pop <- filter(df_global_pop,I>=0.100 & pvalue<0.05 & x!="total")
# rotina de LISA
lLISA_pop <- split(df_sig_pop,seq(nrow(df_sig_pop))) %>%
    future_lapply(.,f_apply_localM,colname="pop_z", future.seed = TRUE)
names(lLISA_pop) <- mutate(df_sig_pop,vname = paste(x,y,sep=" : ")) %>%
    pull(vname)
dftow <- lapply(lLISA_pop,st_drop_geometry) %>% do.call("rbind",.)
rownames(dftow) <- NULL
write_csv(dftow,file = "8_pastas_auxiliares/artigo_suscetiveis/csv/LISA_pop.csv")


#------------- proporção de suscetíveis ---------
##@ apenas os casos significativos
df_sig_psus <- filter(df_global_psus,I>=0.10 & pvalue<0.05 & x!="total")
# rotina de LISA
lLISA_psus <- split(df_sig_psus,seq(nrow(df_sig_psus))) %>%
    future_lapply(.,f_apply_localM, future.seed = TRUE)
names(lLISA_psus) <- mutate(df_sig_psus,vname = paste(x,y,sep=" : ")) %>%
    pull(vname)
dftow <- lapply(lLISA_psus,st_drop_geometry) %>% do.call("rbind",.)
rownames(dftow) <- NULL
write_csv(dftow,file = "8_pastas_auxiliares/artigo_suscetiveis/csv/LISA_psus.csv")


#------------- bolsão de suscetíveis ---------
##@ apenas os casos significativos
df_sig_bsus <- filter(df_global_bsus,I>=0.10 & pvalue<0.05)
# rotina de LISA
lLISA_bsus <- split(df_sig_bsus,seq(nrow(df_sig_bsus))) %>%
    future_lapply(.,f_apply_localM, dfref=df_Bi, colname="Bi_z", future.seed = TRUE)
names(lLISA_bsus) <- mutate(df_sig_bsus,vname = paste(x,y,sep=" : ")) %>%
    pull(vname)
dftow <- lapply(lLISA_bsus,st_drop_geometry) %>% do.call("rbind",.)
rownames(dftow) <- NULL
write_csv(dftow,file = "8_pastas_auxiliares/artigo_suscetiveis/csv/LISA_bsus.csv")


#------------- Impacto marginal na sustibilidade ---------
##@ apenas os casos significativos
df_sig_Isus <- filter(df_global_Isus,I>=0.10 & pvalue<0.05)
# rotina de LISA
lLISA_Isus <- split(df_sig_Isus,seq(nrow(df_sig_Isus))) %>%
    future_lapply(.,f_apply_localM, dfref=df_Bi, colname="Ii_z", future.seed = TRUE)
names(lLISA_Isus) <- mutate(df_sig_Isus,vname = paste(x,y,sep=" : ")) %>%
    pull(vname)
dftow <- lapply(lLISA_Isus,st_drop_geometry) %>% do.call("rbind",.)
rownames(dftow) <- NULL
write_csv(dftow,file = "8_pastas_auxiliares/artigo_suscetiveis/csv/LISA_Isus.csv")


#---------------------------------------------------------------------
#--------------------------- summario LISA ---------------------------
#---------------------------------------------------------------------
#@
library(data.table)
f_summarioLISA <- \(vpath,LISAclass="High-High"){
    #
    dflisa <- fread(vpath)
    #
    Nfxet <- length(unique(dflisa$fxet))
    vbase <- str_extract(vpath, "(?<=LISA_).*?(?=\\.csv)")
    #
    dfret <- dflisa[, .(nHH = sum(class == LISAclass, na.rm = TRUE)), by = ID_MUN]
    dfret$prop <- dfret$nHH*100/Nfxet
    dfret$nK <- Nfxet
    dfret$base <- vbase
    #
    return(dfret)
}
f_ij <- \(dff){
    dff$ID_MUN <- as.character(dff$ID_MUN)
    inner_join(x=dfmun,dff)
}
#
lf <- list.files(path="8_pastas_auxiliares/artigo_suscetiveis/csv",
                 pattern="LISA_",
                 full.names = TRUE)
names(lf) <- str_extract(lf, "(?<=LISA_).*?(?=\\.csv)")
l_summaryLISA <- lapply(lf,f_summarioLISA) %>% lapply(.,f_ij)

# criação dos outros sumários
v_todasLISA <- c('High-High','Low-High','High-Low')
l_summary_todasLISA <- lapply(v_todasLISA,\(li){
    ldf <- lapply(lf,f_summarioLISA,LISAclass=li)
    dfret <- do.call("rbind",ldf)
    dfret$LISAclass=li
    return(dfret)
})
df_summary_todasLISA <- rbindlist(l_summary_todasLISA)
df_summary_todasLISA$base <- gsub("bsus","Pockets of susceptibility",df_summary_todasLISA$base) %>%
    gsub("Isus","Marginal susceptible impact",.) %>%
    gsub("psus","Proportion of susceptibles",.) %>%
    gsub("pop","Overall population",.)

#---------------------------------------------------------------------
#----------------------- contiguidade entre HH -----------------------
#---------------------------------------------------------------------
library(igraph)
# dff <- l_summaryLISA[[1]]
f_continuidade <- \(dff){
    # vizinhança de contiguidade e inclusão do vizinho mais próximo da municipio Ilha
    nb <- poly2nb(dff, queen=TRUE)
    l_MNIsolado <- which(card(nb) == 0)
    for(i in l_MNIsolado){
        dist <- st_distance(dff[i,],dff)
        vizinho <- order(dist)[2] #
        #
        nb[[i]] <- vizinho
        #
        nb[[vizinho]] <- c(nb[[vizinho]],i)
    }
    # quais os municípios que foram HH em pelo menos 1 caso
    vi <- which(dff$nHH!=0)
    # matriz de adj
    adj_matrix <- nb2mat(nb, style = "B", zero.policy = TRUE)
    # apenas os mun de interesse
    adj_matrix_filtrada <- adj_matrix[vi,vi]
    # grafos de conexão estrutural
    g <- graph_from_adjacency_matrix(adj_matrix_filtrada,mode = "undirected")
    # componentes
    componentes <- components(g)
    #
    # atribuição de componente
    dff$comp_conectado <- NA
    dff$comp_conectado[vi] <- componentes$membership
    # retorno
    return(dff)
}
l_graphHH <- lapply(l_summaryLISA,f_continuidade)
