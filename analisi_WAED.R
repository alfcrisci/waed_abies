# install.packages("Compositional")
# install.packages("compositions")
# install.packages("ToolsForCoDa")
# install.packages("easyCODA")

##########################################################################
library(Compositional)
library(easyCODA)
library(vegan)
library(openxlsx)
library(lubridate)
library(ggplot2)
library(viridis)
library(dplyr)
library(purrr)
library(ggfortify)
library(cluster)
library(robCompositions)
library(highmean)
cat("\014") 
setwd("/home/alf/Scrivania/lav_waed_coda/")

# cat("\014") 
##########################################################################
mat_on_tot=read.xlsx("matwork_A_alba_nebrodensis.xlsx" ,3)
mat_on_mono=read.xlsx("matwork_A_alba_nebrodensis.xlsx" ,4)
##########################################################################
# cleaning names

names(mat_on_tot)=gsub(".Results","",names(mat_on_tot))
names(mat_on_mono)=gsub(".Results","",names(mat_on_mono))
names(mat_on_tot)=gsub(" ","",names(mat_on_tot))
names(mat_on_mono)=gsub(" ","",names(mat_on_mono))

##########################################################################
# matrix 

mat_on_tot_ord=as.matrix(mat_on_tot[,7:30])
mat_on_tot_ord.rel <-decostand(mat_on_tot_ord, method = "total") # relA

mat_on_mono_ord=as.matrix(mat_on_mono[,7:23])
mat_on_mono_ord.rel <-decostand(mat_on_mono_ord, method = "total") # relA
##########################################################################
# https://mb3is.megx.net/gustame/dissimilarity-based-methods/nmds
mat_on_tot_ord.rel_NMS <-
  metaMDS(mat_on_tot_ord.rel,
          distance = "bray",
          k = 3,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)

goodness(mat_on_tot_ord.rel_NMS)
pdf(file = "Full_RelAbundance_stressplot.pdf",width=8, height=6)
stressplot(mat_on_tot_ord.rel_NMS) # Produces a Shepards diagram
dev.off()

pdf(file = "NMS_monoontotal.pdf",width=8, height=6)

plot(mat_on_tot_ord.rel_NMS, "sites")   # Produces distance 
orditorp(mat_on_tot_ord.rel_NMS, "species",col="red",air=0.01,cex = 0.65)   
ordiellipse(
  mat_on_tot_ord.rel_NMS,
  mat_on_tot$flag,
  display = "sites",
  draw = c("polygon"),
  col = NULL,
  border = c("green","red","orange"),
  lty = c(1, 2, 2),
  lwd = 1.5
)
title("Ordination ( nMDS)  MONO ON TOTAL ")
legend('topleft', mat_on_mono$flag, col=c("green","red","orange"),lty = c(1, 2, 2),cex=0.75)
dev.off()

#################################################################################
# Aalb_vuoto Aneb_pieno Aneb_vuoto
##########################################################################
mat_on_mono_ord.rel_NMS <-
  metaMDS(mat_on_mono_ord.rel,
          distance = "bray",
          k = 3,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)

goodness(mat_on_mono_ord.rel_NMS)
pdf(file = "monoontotmono_stressplot.pdf",width=8, height=6)
stressplot(mat_on_mono_ord.rel_NMS) # Produces a Shepards diagram
dev.off()

pdf(file = "NMS_monoontotmono.pdf",width=8, height=6)
plot(mat_on_mono_ord.rel_NMS, "sites")   # Produces distance 
orditorp(mat_on_mono_ord.rel_NMS, "species",col="red",air=0.01,cex = 0.65)   
ordiellipse(
  mat_on_mono_ord.rel_NMS,
  mat_on_mono$flag,
  display = "sites",
  draw = c("polygon"),
  col = NULL,
  border = c("green","red","orange"),
  lty = c(1, 2, 2),
  lwd = 1.5
)
title("Ordination ( nMDS)  MONO ON TOT MONO ")
legend('topleft', mat_on_mono$flag, col=c("green","red","orange"),lty = c(1, 2, 2),cex=0.75)
dev.off()

#############################################################################################
autoplot(pam(mat_on_tot[7:30], 2),data = mat_on_tot[7:30],label = TRUE, label.size = 5,frame = TRUE, frame.type = 'norm')
autoplot(pam(mat_on_mono[7:23], 2),data = mat_on_mono[7:23],label = TRUE, label.size = 5,frame = TRUE, frame.type = 'norm')

data_mat=cbind(mat_on_tot[,c(2,4)],mat_on_tot[7:30])

for (i in 3:26) {
message(names(data_mat)[i])
print(summary(glm(data_mat[,i]~SP_ID+seed_status,data=data_mat)))
print(anova(glm(data_mat[,i]~SP_ID+seed_status,data=data_mat)))

}



mat_delta_on_tot=psych::describeBy(mat_on_tot[7:30], mat_on_tot$SP_ID, mat = TRUE)
mat_delta_on_tot=psych::describeBy(mat_on_tot[7:30], mat_on_tot$flag, mat = TRUE)


maov(mat_on_tot_ord[,c("Germacrene.D","α-humulene","β-caryophillene","Limonene","α-pinene")],as.factor(mat_on_tot$SP_ID))

# https://rstudio-pubs-static.s3.amazonaws.com/533895_305b8af2a0104ecdb6d30ad9d8b9fffa.html NMDS Ordination
# https://www.rpubs.com/RGrieger/545184
# https://academic.oup.com/biomet/article/105/1/115/4591648 # highmean