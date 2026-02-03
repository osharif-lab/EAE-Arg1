library(readxl)
library(pheatmap)

################
# Read in data #
################

metabo_df <- read_excel("metabolomics/Report_M085_Arg_metabolism_20251121.xlsx",na = "N/F")
metabo_dat<-as.data.frame(metabo_df[,c("Argininosuccinate","Citrulline","Creatine","L-Arginine","L-Glutamic Acid","Ornithine",'Urea')])
colnames(metabo_dat)<-c("Argininosuccinate","Citrulline","Creatine","Arginine","Glutamate","Ornithine",'Urea')
rownames(metabo_dat)<-metabo_df$`Sample ID`

############
# Heatmaps #
############

# Extended Data Figure 9b
pdf("metabolomics/Heatmap_metabo_healthy_EAE_Cells.pdf", width=6,height=2.5)
p<-pheatmap(log2(t(metabo_dat)),
            color = rev(colorRampPalette(c("#B2182B", "#FDDBC7", "#2166AC"))(n = 256)),
            border_color = "black",
            scale = "row",
            angle_col = 45,
            fontsize_row = 9,
            cluster_cols = F,
            cluster_rows = T,
            cellwidth = 15,
            annotation_names_col = F
)
print(p)
dev.off()

for (i in 1:ncol(metabo_dat)){
  res<-t.test(log2(metabo_dat)[1:5,i],log2(metabo_dat)[6:10,i])
  print(colnames(log2(metabo_dat))[i])
  print(res$p.value)
}



##############
# Polyamines
##############
#arginine, orthinine, putrescine, spermidine, sperminine

# SC
metabo_df <- read_excel("metabolomics/AVMKe33_AVMKe33a_SC HPR polyamines.xls",sheet = "data")
colnames(metabo_df)
head(metabo_df)

metabo_df$Group
metabo_df$Condition<-rep('Peak',nrow(metabo_df))
metabo_df$Condition[grepl('ealth',metabo_df$Group)]<-'Healthy'
metabo_df$Condition[grepl('emiss',metabo_df$Group)]<-'Remission'
metabo_df$Condition
cbind(metabo_df$Group,metabo_df$Condition)
metabo_df$Experiment<-rep(1,nrow(metabo_df))
metabo_df$Experiment[grepl('BL6',metabo_df$Group)]<-2

ids<-c(1,2,3,4,5,13,14,15,16,17,18,6,7,19,20,21,22,23,24)
metabo_df2<-metabo_df[ids,]
metabo_dat<-as.data.frame(metabo_df2[,c("Ornithine","Putrescine","Spermidine","N1-Acetylspermine","N8-Acetylspermidine or isomer", "1,3 Diaminopropane")])
metabo_dat<-as.data.frame(metabo_df2[,c("Ornithine","Putrescine","Spermidine","Spermine","N1-Acetylspermine","N8-Acetylspermidine or isomer", "1,3 Diaminopropane")])

metabo_dat_log<-log2(metabo_dat)

# batch correction

exp1_means<-colMeans(metabo_dat_log[metabo_df2$Experiment==1,])
exp1_std<-apply(metabo_dat_log[metabo_df2$Experiment==1,],MARGIN = 2,FUN = sd)

exp2_means<-colMeans(metabo_dat_log[metabo_df2$Experiment==2,])
exp2_std<-apply(metabo_dat_log[metabo_df2$Experiment==2,],MARGIN = 2,FUN = sd)

metabo_dat_log_scaled<-matrix(rep(NA, nrow(metabo_dat)*ncol(metabo_dat)),ncol=ncol(metabo_dat))
for (i in 1:nrow(metabo_df2)){
  if(metabo_df2$Experiment[i]==1){
    metabo_dat_log_scaled[i,]<-as.matrix(((metabo_dat_log[i,] - exp1_means)/exp1_std),nrow=1)
  } else {
    metabo_dat_log_scaled[i,]<-as.matrix(((metabo_dat_log[i,] - exp2_means)/exp2_std),nrow=1)
  }
  
}
dim(metabo_dat_log_scaled)
rownames(metabo_dat_log_scaled)<-metabo_df2$Group
colnames(metabo_dat_log_scaled)<-colnames(metabo_dat)

rownames(metabo_dat)<-metabo_df2$Group

# Figure 2j
pdf("metabolomics/Heatmap_metabo_polyamines_SC.pdf", width=10,height=5)
p<-pheatmap(t(metabo_dat_log_scaled),#log2(t(metabo_dat)),
            color = rev(colorRampPalette(c("#B2182B", "#FDDBC7", "#2166AC"))(n = 256)),
            border_color = "black",
            #scale = "row",
            column_names_rot = 45,
            fontsize_row = 9,
            cluster_cols = F,
            cluster_rows = F,
            cellwidth = 15,
            cellheight = 20,
            # annotation_col =  colnames(log_lip_sum),
            # annotation_row = rownames(log_lip_sum),
            annotation_names_col = F
)
dev.off()

for (i in 1:ncol(metabo_dat)){
  res<-t.test(metabo_dat_log_scaled[1:11,i],metabo_dat_log_scaled[12:19,i])
  print(colnames(metabo_dat_log_scaled)[i])
  print(res$p.value)
}


###################
#CSF
metabo_df <- read_excel("metabolomics/AVMKe33_AVMKe33a_CSF HPR polyamines.xls",sheet = "data")
colnames(metabo_df)
head(metabo_df)
metabo_dat<-as.data.frame(metabo_df[,c("Ornithine","Putrescine","Spermidine","Spermine","N1-Acetylspermine","N8-Acetylspermidine or isomer", "1,3 Diaminopropane")])
#metabo_dat<-as.data.frame(metabo_df[,c("Arginine","Ornithine","Putrescine","Spermidine","Spermine")])
metabo_dat
metabo_df$Group
rownames(metabo_dat)<-metabo_df$Group

library(stringr)
metabo_df$Condition<-str_sub(metabo_df$Group,1,-3)
take_ids<-which(metabo_df$Condition!='Remission')


# Figure 2j
pdf("metabolomics/Heatmap_metabo_polyamines_CSF_Putrescine.pdf", width=8,height=5)
p<-pheatmap(log2(t(metabo_dat[take_ids,2])),
            color = rev(colorRampPalette(c("#B2182B", "#FDDBC7", "#2166AC"))(n = 256)),
            border_color = "black",
            scale = "row",
            angle_col = 45,
            fontsize_row = 9,
            cluster_cols = F,
            cluster_rows = F,
            cellwidth = 15,
            cellheight = 20,
            # annotation_col =  colnames(log_lip_sum),
            # annotation_row = rownames(log_lip_sum),
            annotation_names_col = F
)
print(p)
dev.off()

