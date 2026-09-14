rm(list=ls())
cat("\014")

library(terra)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(hydroTSM)
"%>%" = magrittr::`%>%`

year <- 2026
month <- 8

## load data ----
SPI <- paste0("chirps-v3.0.spi3_pearson_category.1981-", year,"-", sprintf("%02d", month),".colombia.nc") %>% 
  rast(., subds = "precip")
AH <- vect("../shp_Area_hid_col/Zonificacion_hidrografica_2013_AH.shp")

nom_AH <- AH$NOM_AH # c("Caribe", "Magdalena Cauca", "Orinoco", "Pacifico", "Amazonas")

dates <- seq(as.Date("1981-01-01"), 
             paste0(terra::nlyr(SPI) %/% 12 + 1981,"-",terra::nlyr(SPI) %% 12,"-01") %>% as.Date(), 
             by = "month")
n.dates <- length(dates)

## preprocess cells ----
data <- terra::extract(SPI[[3]], AH, cells = TRUE, xy =T)

coords <- data[,c("cell","ID","x","y")]
coord.area <- split(coords, coords$ID %>% as.factor())
# n.cells <- list(3370, 8871, 11261, 2516, 11079); names(n.cells) <- nom_AH

areas.value <- sapply(coord.area, nrow) 

## classification by hydrographic area ----
# save(aux, "aux", file= "SPI-CHIRPS_monitoring.RData")

class.area.f <- function(aux.f, n.dates, dates, nom_AH){
  
  class.area_f <- list()
  
  for( i in 1:length(nom_AH)){
    print(nom_AH[i])
    class.area_f[[i]] <- matrix(NA, nrow = n.dates, ncol= 7)
    colnames(class.area_f[[i]]) <- c("extremely wet","severely wet","moderately wet",
                                   "Normal",
                                   "moderately dry","severely dry","extremely dry")
    cells.area <- aux.f$ID == i
    aux.f.date <- aux.f[ ,2:(n.dates + 1)] %>% data.frame(date= .)
    
    n.cell.cat <- apply(aux.f.date, 2, # number of cells in each SPI category in each month
                        function(x, c.area){ 
                          y <- x[ c.area ] %>% factor() %>% summary()
                          return(y)
                        },
                        cells.area) 
    
    if(n.cell.cat %>% is.array()){
      exist.na <- which(n.cell.cat %>% row.names() == "NA's") 
      if ( length(exist.na) > 0) n.cell.cat <- n.cell.cat[ -exist.na ,1]
      n.cell.cat <- n.cell.cat %>% as.data.frame() %>% as.list()
      bool.all.cat <- TRUE
      } else bool.all.cat <- FALSE
    
    # check internally how many time steps have non-NAN categories and which
    cat.month <- sapply(n.cell.cat, function(x) x %>% names() %>% as.numeric() %>% sum() %>% is.na()) 
    
    for( j in seq(1,n.dates)[ !cat.month ]  ){
      
      if (!bool.all.cat) col.aux <- n.cell.cat[[j]] %>% names() %>% as.numeric() else col.aux <- 1:7
      class.area_f[[i]] [ j , col.aux] <- n.cell.cat[[j]]
      class.area_f[[i]] [ j , c(1:7)[-col.aux]] <- 0
    }
    
    rm(bool.all.cat)
    class.area_f[[i]] <- (class.area_f[[i]] / areas.value[i]) * 100
  }
  
  class.area_f <- lapply(class.area_f, function(y){z <- as.data.frame(y);  row.names(z) <- dates;  return(z)})
  names(class.area_f) <- nom_AH
  
  return(class.area_f)
}

if (file.exists("SPI-CHIRPS_monitoring.RData")){
  
  load("SPI-CHIRPS_monitoring.RData")
  load("SPI-CHIRPS_monitoring_AreaClass.RData")
  
  last.dates <- row.names(class.area[[1]] ) %>% tail(1)
  if ( last.dates != tail(dates, 1)){
    
    comp.dates <- dates %in% row.names(class.area[[1]])
    aux.new <- terra::extract(SPI[[  !comp.dates ]], AH, cells = TRUE, xy =T)
    
    class.area_new <- class.area.f(aux.new, sum(!comp.dates), dates[ ! dates %in% row.names(class.area[[1]])], nom_AH)
    
    class.area <- mapply(function(x,y) rbind(x,y), class.area, class.area_new, SIMPLIFY = FALSE)
    
    aux <- cbind(aux[,1:(n.dates+1)], aux.new[,2:(sum(!comp.dates)+1)], aux[,(n.dates+2):ncol(aux)])
  }
  
  
  
} else {
  aux <- terra::extract(SPI, terra::vect(AH), cells = TRUE, xy =T)
  
  class.area <- class.area.f(aux, n.dates, dates, nom_AH)
  
  # save(class.area, "class.area", file="SPI-CHIRPS_monitoring_AreaClass.RData")
}

save(aux, "aux", file= "SPI-CHIRPS_monitoring.RData")
save(class.area, "class.area", file="SPI-CHIRPS_monitoring_AreaClass.RData")

## joining areas for area analysis ----

sel.col.areas <- c(1,2,4)  # c(1,2,4)

SPI.class <- class.area[sel.col.areas] # [c(1,2,4)]
tot.area <- areas.value[sel.col.areas] %>% sum() # [c(1,2,4)]

SPI.avg <- SPI.class[[1]]

for ( i in 1:ncol(SPI.avg)){
  SPI.avg[,i] <- mapply( 
    function(x,y){
      x[,i] * y
    },
    SPI.class, # [c(1,2,4)]
    areas.value[sel.col.areas] / tot.area) %>% as.data.frame() %>%
    apply(., MARGIN = 1, FUN = sum)
  # SPI.class[["Caribe"]][,i]* n.cells$Caribe / tot.area + ...
}

data.g <- SPI.avg %>%
  zoo(, order.by = row.names(.)) %>%
  fortify(melt = TRUE) %>% within(Index <- as.Date(Index))

## continuous area plot ----
paleta <- brewer.pal(9,name = "BrBG")[c(1,2,3,5,7,8,9)] %>% rev(); paleta[4] <- "00"

p.title <- paste("SPI-3 Classification of affected area in", paste(nom_AH[sel.col.areas], collapse = ", ") ,"regions in Colombia")

p <- ggplot()+
  geom_area(data=data.g, aes(Index, Value, fill=Series)) + 
  scale_fill_manual(values = paleta, name="Class.") + 
  scale_y_continuous(expand = c(0.01,0.01)) +
  xlab("Date") + ylab("[%] Area affected by precipitation deficit (SPI-3)")+ 
  ggtitle(p.title)+
  theme_bw()+theme(legend.position="bottom", legend.title = element_text(size = 14), legend.text = element_text(size = 12), legend.background = element_rect(linetype="solid", colour ="black"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed"),
                   axis.title = element_text(size=14, face = "bold"), axis.text = element_text(size=12),
                   axis.title.x = element_blank())


ggsave("SPI3-CHIRPS_monitoring_Col.png", plot= p + 
         geom_text(aes(x=as.Date("2006-01-01"),y=50,label="@duque-gardeazabal\nData: CHIRPS\nNo Orinoco, No Amazon")) +
         scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = c(0.01,0.01)),
       dpi=300,width = 1400*3/300,height = 600*3/300, units = "in")

library(ggforce)

p1 <- p + 
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0.01,0.01)) +
  facet_zoom(xlim=c("2013-01-01","2016-12-01") %>% as.Date()) + 
  theme(legend.position = c(0.7,0.5), legend.direction = "horizontal",
          axis.text.x = element_text(margin=margin(1,5,5,5),vjust = -1, size=8))
ggsave("SPI3-CHIRPS_monitoring_Col_zoom2015-16.png", plot= p1,
       dpi=300,width = 1400*3/300,height = 600*3/300, units = "in")


p2 <- ggplot()+
  geom_area(data=data.g %>%  subset(., Index >= as.Date("2025-07-01")), 
            aes(Index, Value, fill=Series)) + 
  scale_fill_manual(values = paleta, name="Class.") + 
  scale_x_date(date_breaks = "3 months", date_labels = "%Y-%m", expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  geom_text(data= data.frame(x=as.Date("2025-12-01"),y=20,texto="@duque-gardeazabal\ndata: CHIRPS"),
            aes(x, y, label= texto))+
  xlab("Date") + ylab("[%] Area affected by precipitation deficit (SPI-3)")+ 
  ggtitle(p.title)+
  theme_bw()+theme(legend.position="bottom", legend.title = element_text(size = 14), legend.text = element_text(size = 12), legend.background = element_rect(linetype="solid", colour ="black"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed"),
                   axis.title = element_text(size=14, face = "bold"), axis.text = element_text(size=12),
                   axis.title.x = element_blank()) 
ggsave("SPI3-CHIRPS_monitoring_Col_zoom2026.png", plot= p2,
       dpi=300,width = 1400*3/300,height = 600*3/300, units = "in")


## East and central Pacific  ----

# organise the intensity of droughts
data.g2 <- subset(data.g, 
                  Series %in% c("extremely dry", "severely dry", "moderately dry") & 
                    !is.na(Value)) 

f.intensity <- function(x, y){
  pre <- x %>% dplyr::slice_max(order_by = Value, n = 1) %>% dplyr::select(Series) %>% unlist()
  if( length(pre) > 1) pre <- pre[1]
  
  if (pre == "moderately dry"){
    sev <- with(x, {Value[Series == "severely dry"] + Value[Series == "extremely dry"] })
    if ( sev > with(x, 1.3* Value[Series == "moderately dry"])){
      pre <- "severely dry"
    }
  }
  
  return(as.data.frame(pre))
}

Intensity <- data.g2 %>%
  within(., Series <- as.character(Series)) %>%
  dplyr::group_by(Index) %>%
  dplyr::group_modify(f.intensity)

data.g2 <- data.g2 %>%
  dplyr::group_by(Index) %>%
  dplyr::summarise(, drought = sum(Value)) %>%
  merge(., Intensity, by= "Index") %>%
  dplyr::mutate(, Month = format(Index, format="%m") %>% as.numeric(),
                Year = format(Index, format="%Y") %>% as.numeric())

# organise in ENSO diversity
ENSO.e <- c(1982, 1997, 2015)
ENSO.c <- c(1986, 1991, 1994, 2002, 2009,  2023)

ELI.e <- data.g2 %>%
  subset(., Year %in% ENSO.e | Year %in% (ENSO.e+1)) %>% 
  split(., rep( 1:length(ENSO.e), each=24)) %>% 
  lapply(., function(x){
    y <- subset(x, select= -Year)
    y <- subset(y, select= -Index)
    y[13:24,"Month"] <- 13:24
    return(y)} )
names(ELI.e) <- paste(ENSO.e, ENSO.e+1, sep="-")

ELI.c <- data.g2 %>%
  subset(Year %in% ENSO.c | Year %in% (ENSO.c+1)) %>% 
  split(., rep( 1:length(ENSO.c), each=24)) %>% 
  lapply(., function(x){
    y <- subset(x, select= -Year)
    y <- subset(y, select= -Index)
    y[13:24,"Month"] <- 13:24
    return(y)} )
names(ELI.c) <- paste(ENSO.c, ENSO.c+1, sep="-")

ELI.2026 <- data.g2 %>% subset(Year>=2026, select = -Year) %>% subset(select = -Index) %>% 
  within(., L1 <- "2026") %>% magrittr::set_colnames(c("drought","Series","Month","Event"))

# organise for plotting

div <- "East Pacific events"
ELI.e.g <- reshape2::melt(ELI.e, id=c("Month","drought")) %>% 
  subset(, select = -variable) %>%
  magrittr::set_colnames(c("Month","drought","Series","Event")) %>%
  rbind(., ELI.2026[,c(3,1,2,4)]) %>% 
  cbind(., div) %>% 
  within(., Series <- factor(Series, levels = c("extremely dry", "severely dry", "moderately dry")))
div <- "Central Pacific events"
ELI.c.g <- reshape2::melt(ELI.c, id=c("Month","drought")) %>% subset(,select = -variable) %>%
  magrittr::set_colnames(c("Month","drought","Series","Event")) %>% 
  cbind(., div) %>% 
  within(., Series <- factor(Series, levels = c("extremely dry", "severely dry", "moderately dry")))


# plotting east  vs central pacific ----
lab.sign <- data.frame(x=2.5,y=10,texto="@duque-gardeazabal\ndata: CHIRPS", div) %>% 
  rbind(.,.); lab.sign[2,4] <- "East Pacific events"
paleta <- brewer.pal(9, "BrBG")[1:3] %>% magrittr::set_names(c("extremely dry", "severely dry", "moderately dry"))
p.title <- paste("Evolution of area under Drought El Niño events in the", paste(nom_AH[sel.col.areas], collapse = ", ") ,"regions - Colombia")


t.last.month <- rbind( ELI.c.g %>% subset(., Month == 12),
                       ELI.e.g %>% subset(., Month == 12))
library(ggrepel)
ggplot() + facet_wrap(.~ div, ncol = 1) +
  geom_line(data= ELI.c.g, aes(Month, drought, group = Event, colour = Series), linewidth = 2)+
  
  geom_line(data= ELI.e.g, aes(Month, drought, group = Event, colour = Series), linewidth = 2)+ 
  scale_color_manual(values=paleta, breaks =c("extremely dry", "severely dry", "moderately dry"),
                     name="Intensity")+
  
  scale_x_continuous(labels= rep(month.abb,2), breaks=1:24, expand = c(0.01,0.01))+
  geom_text(data = lab.sign, aes(x, y, label= texto))+
  geom_text_repel(data = t.last.month, aes(x = Month , y = drought, label= Event), nudge_x = -0.5, nudge_y = 5)+
  labs(y="[%] Area affected by precipitation deficit (SPI-3)", 
       x="Month since January of year 0 and year 1",
       title=p.title)+
  theme_bw(  )+ theme(legend.key.width = unit(1, "cm"), legend.position = c(0.25,0.95), 
                      legend.background = element_rect(color="black"), legend.direction = "horizontal",
                      axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed", color="grey85"),
                      axis.title = element_text(size=14, face = "bold"), axis.text = element_text(size=12),
                      strip.text = element_text(size=16, face = "bold"),
                      title = element_text(size=16))+
  guides(color = guide_legend(override.aes = list(linewidth = 3), theme = theme(legend.text = element_text(size=12))),
         linewidth = "none")

ggsave("SPI-3_Evolution_central_East_Pacific_ENSO.png",
       dpi=300,width = 1400*3/300,height = 900*3/300, units = "in")


last.month <- substr(last.dates,6,7) %>% as.numeric()
t.last.month <- rbind( ELI.c.g %>% subset(., Month == last.month +1),
                       ELI.e.g %>% subset(., Month == last.month +1))

ggplot() + facet_wrap(.~ div, ncol = 2) +
  geom_line(data= ELI.c.g %>% subset(., Month <= last.month +1),
            aes(Month, drought, group = Event, colour = Series), linewidth = 2)+
  
  geom_line(data= ELI.e.g %>% subset(., Month <= last.month +1),
            aes(Month, drought, group = Event, colour = Series), linewidth = 3)+ 
  scale_color_manual(values=paleta, breaks =c("extremely dry", "severely dry", "moderately dry"),
                     name="Intensity")+
  
  scale_x_continuous(labels= rep(month.abb,2), breaks=1:24, expand = c(0.01,0.01))+
  geom_text(data = lab.sign, aes(x, y, label= texto))+
  geom_text_repel(data = t.last.month, aes(x = Month , y = drought, label= Event), nudge_x = -0.5, nudge_y = 5)+
  labs(y="[%] Area affected by precipitation deficit (SPI-3)", 
       x="Month since January of year 0 and year 1",
       title=p.title)+
  theme_bw(  )+ theme(legend.key.width = unit(1, "cm"), legend.position = c(0.25,0.95), 
                      legend.background = element_rect(color="black"), legend.direction = "horizontal",
                      axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed", color="grey85"),
                      axis.title = element_text(size=17, face = "bold"), axis.text = element_text(size=16),
                      strip.text = element_text(size=16, face = "bold"),
                      title = element_text(size=16),
                      panel.spacing = unit(1, "cm"))+
  guides(color = guide_legend(override.aes = list(linewidth = 3), theme = theme(legend.text = element_text(size=12))),
         linewidth = "none")

ggsave("SPI-3_Monitoring_Evolution_ENSO.png",
       dpi=300,width = 1400*3/300,height = 700*3/300, units = "in")
