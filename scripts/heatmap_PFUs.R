library(tidyverse)
library(scales)
library(ggpubr)
library(data.table)
library(vegan)

# with the final dataset for our region

library(tidyverse)

rm(list=ls())

aci <- read_tsv("aci_in_lab_20230628.tsv")

aci_all <- read_csv("Acinetobacter baumannii strains_Phages_20230802.csv")

# merge without adding space
aci_all$Strain <- paste(aci_all$Name, aci_all$`Identifier of the strain`, sep="")

#Left  join in R:  merge() function takes df1 and df2 as argument along with all.x=TRUE  
#there by returns all rows from the left table, and any rows with matching keys from the right table.

Aci_comb = merge(x=aci_all,y=aci,by="Strain",all.x=TRUE) 

#Aci_comb %>% distinct(country)

# új oszlopok létrehozása 2 összevonásából

Aci_comb$strain_city <- paste(Aci_comb$Strain, "-", Aci_comb$city)

Aci_comb$ST_KL <- paste(Aci_comb$MLST, "-", Aci_comb$KL)

# Aci_comb %>% distinct(ST_KL)

# filter out those for which we don't know the ST_KL

Aci_comb <- Aci_comb %>% filter (ST_KL != "NA - NA")

# make the table shorter

ACI  <- Aci_comb %>% select(c(1, 5:21, 31, 33, 85, 86))

###### filter for the 5 countries from the region we got samples

ctry3 <- ACI %>% filter(country == "hungary" | country == "romania" | country == "serbia" | country == "bosnia_and_herzegovina" | country == "montenegro") %>%
  filter (ST_KL == "ST2 - KL3"| ST_KL == "ST636 - KL40" | ST_KL == "ST492 - KL104" | ST_KL == "ST2 - KL2" |
            ST_KL == "ST1 - KL1" | ST_KL == "ST2 - KL9" | ST_KL == "ST2 - KL12" | ST_KL == "ST1 - KL17" |
            ST_KL == "ST2 - KL7" | ST_KL == "ST2 - KL32" | ST_KL == "ST2 - KL77")

# ctry3_F <- ctry3 %>% filter (!is.na (Highwayman) & !is.na (Silvergun) & !is.na (Fanak) & !is.na (PhT2_Hun) & !is.na (Porter) & !is.na (Dino) & !is.na (Fishpie) & !is.na (Tama) & !is.na (Margaret) & !is.na (ABW132) &  
#!is.na (ABW311) & !is.na (Navy4_Hun) & !is.na (Rocket) & !is.na (Konradin_Hun) & !is.na (KissB))

# kiszűrni minden sort, ahol bármelyik fág NA, illetve fertőz, de nincs PFU (YES, Yes, yes) 
#összevonvan felírva függvénnyel

ctry3_F <- ctry3 %>% filter (if_all(c(Highwayman,Silvergun,Fanak,PhT2_Hun,Porter,Dino,Fishpie,Tama,Margaret,ABW132,ABW311,Navy4_Hun,Rocket,Konradin_Hun,KissB), function(x) !is.na(x))) %>%
  filter (if_all(c(Highwayman,Silvergun,Fanak,PhT2_Hun,Porter,Dino,Fishpie,Tama,Margaret,ABW132,ABW311,Navy4_Hun,Rocket,Konradin_Hun,KissB), function(x) !`%in%`(x, c("YES", "Yes", "yes"))))

# I could do it like this for each phage separately
# filter(! Silvergun %in% c("YES", "Yes", "yes"), ! Highwayman %in% c("YES", "Yes", "yes"))


ctry3_F %>% distinct(country) #5
ctry3_F %>% distinct(ST_KL) #11


ACI_long <- pivot_longer(data = ctry3_F, cols = 4:18, values_to = "value")

ACI_long$value[ACI_long$value==1] <- "NA"
ACI_long$value[ACI_long$value=="NO"] <- 0

ACI_long <- transform(ACI_long, value=as.numeric(value))


# ST-KL és fágsorrend meghatározása

ACI_long$ST_KL <- factor(ACI_long$ST_KL, levels=c("ST2 - KL3", "ST636 - KL40", "ST492 - KL104", "ST2 - KL2", "ST1 - KL1", "ST2 - KL9" , "ST2 - KL12" , "ST1 - KL17" ,
                                                  "ST2 - KL7" , "ST2 - KL32" , "ST2 - KL77"))

ACI_long$name <- factor(ACI_long$name, levels=c("Highwayman", "Silvergun", "Fanak", "PhT2_Hun", "Porter", "Dino", "Fishpie", "Tama", "Margaret" ,"ABW132",  
                                                "ABW311", "Navy4_Hun", "Rocket", "Konradin_Hun", "KissB"))

# separate MLST and KL

plot <- ggplot(ACI_long)+
  geom_tile(aes(strain_city,name,fill=value,color="Lysis from without"))+
  scale_fill_gradientn(name="log10 titer (PFU/mL)", colours = c("#ffefd8","#feea8a","#cd5d5c"), 
                       values = rescale(c(0,6,12)),
                       guide = "colorbar", limits=c(0,12),
                       na.value = "gray")+
  scale_colour_manual(values=NA)+
  guides(colour=guide_legend("", override.aes=list(fill="gray", colour="transparent")))+
  ggh4x::facet_nested(.~MLST + KL, scales="free_x", space="free")+
  ylab("")+
  xlab("")+
  theme(axis.text.x = element_text(angle=90, hjust=0, vjust=0.5),
        strip.background = element_rect(color="gray70", fill="white"),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_text(face="bold"),
        strip.text.x = element_text(face="bold", size=6),
        legend.position = "top",
        panel.background = element_blank(),
        panel.border = element_rect(color="transparent", fill="transparent"))+
  scale_x_discrete(position="bottom")+
  scale_y_discrete(limits=rev)

# ST-KL 

plot <- ggplot(ACI_long)+
geom_tile(aes(strain_city,name,fill=value,color="Lysis from without"))+
  scale_fill_gradientn(name="log10 titer (PFU/mL)", colours = c("#ffefd8","#feea8a","#cd5d5c"), 
                       values = rescale(c(0,6,12)),
                       guide = "colorbar", limits=c(0,12),
                       na.value = "gray")+
  scale_colour_manual(values=NA)+
  guides(colour=guide_legend("", override.aes=list(fill="gray", colour="transparent")))+
  ggh4x::facet_nested(.~ST_KL, scales="free_x", space="free")+
  ylab("")+
  xlab("")+
  theme(axis.text.x = element_text(angle=90, hjust=0, vjust=0.5),
        strip.background = element_rect(color="gray70", fill="white"),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_text(face="bold"),
        strip.text.x = element_text(face="bold", size=6),
        legend.position = "top",
        panel.background = element_blank(),
        panel.border = element_rect(color="transparent", fill="transparent"))+
  scale_x_discrete(position="bottom")+
  scale_y_discrete(limits=rev)

# blue coloring

plot <- ggplot(ACI_long)+
  geom_tile(aes(strain_city,name,fill=value,color="Lysis from without"))+
  scale_fill_gradientn(name="log10 titer (PFU/mL)", colours = c("white","cadetblue1","steelblue"), 
                       values = rescale(c(0,6,12)),
                       guide = "colorbar", limits=c(0,12),
                       na.value = "azure3")+
  scale_colour_manual(values=NA)+
  guides(colour=guide_legend("", override.aes=list(fill="azure3", colour="transparent")))+
  ggh4x::facet_nested(.~ST_KL, scales="free_x", space="free")+
  ylab("")+
  xlab("")+
  theme(axis.text.x = element_text(angle=90, hjust=0, vjust=0.5),
        strip.background = element_rect(color="gray70", fill="white"),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_text(face="bold"),
        strip.text.x = element_text(face="bold", size=6),
        legend.position = "top",
        panel.background = element_blank(),
        panel.border = element_rect(color="transparent", fill="transparent"))+
  scale_x_discrete(position="bottom")+
  scale_y_discrete(limits=rev)

ggsave("heatmap_phages_PFU_0802.pdf", device = "pdf", dpi = 800, units = "mm", width = 400, height = 120, 
       bg = "transparent")
