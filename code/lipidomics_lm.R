library(readxl)
library(stringr)
library(limma)
library(xlsx)
library(EnhancedVolcano)

################
# Read in data #
################

info<-read_excel("lipidomics/EAECM Lipidomics Sample List.xlsx")
lipids_df<-read_excel("lipidomics/Report_M085_LIPID01_20251117.xlsx",na = "<LOD")
lipids_df<-(as.data.frame(lipids_df))

lip_data<-lipids_df[,-c(1:4)]

# Replace <LODs with values

rat<-1 #coefficient nect to <LOD replacement (min. measured) value
verbose<-T
dim(lip_data)
rownames(lip_data)
lip_data[lip_data==0]<-NA
lip_data_repl<-lip_data
for (i in 1:ncol(lip_data)){
  zero.ids<-which(is.na(lip_data[,i]))
  if (length(zero.ids)>0 & length(zero.ids)<nrow(lip_data)){
    zero.ids
    replacement_value<-rat*min(lip_data[-c(zero.ids),i],na.rm=T)
    print(replacement_value)
    lip_data_repl[zero.ids,i]<-replacement_value
    if (verbose){
      print(paste0(length(zero.ids),' zero values replaced in column ',i, ': ',colnames(lip_data)[i]))
    }
  }
}

which(colSums(is.na(lip_data_repl))!=0)
lip_data<-lip_data[,-which(colnames(lip_data) %in% (c('PA 36:4','PA 37:4')))]
lip_data_repl<-lip_data_repl[,-which(colnames(lip_data_repl) %in% (c('PA 36:4','PA 37:4')))]
lipids_df<-lipids_df[,-which(colnames(lipids_df) %in% (c('PA 36:4','PA 37:4')))]

log_df<-log2(lip_data_repl)
df<-data.frame(cbind(info,log_df))

# data log2 and replaced <LOD values
write.xlsx2(df, paste0("lipidomics/log_data_replaced_NAs.xlsx"), row.names= F, quote = F)

# Replace with formula-friendly symbols
info$replicate<-str_sub(info$CONDITION, start= -1)
info$condition<-str_sub(info$CONDITION, 1,-3)
info$condition<-str_replace(info$condition,' Replicate','')
info$condition<-str_replace(info$condition,'-','.')
info$condition<-str_replace(info$condition,' ','_')
info$condition<-str_replace(info$condition,'-','.')
info$condition<-str_replace(info$condition,' ','_')
condition<-info$condition
table(condition)

#Lipid class colors for plotting
lipid_colors<-c(
  Cer = 'darkblue',   
  SM  = "red3",    
  DAG = "darkorchid1",  
  LPA = "goldenrod1",   
  LPC = "slategray",   
  LPE = "darkgoldenrod4",   
  PA  = "steelblue2",   
  PC  = "mediumpurple4",   
  PE  = "forestgreen",   
  PG  = "dodgerblue4",   
  PI  = "greenyellow",   
  PS  = "pink1",   
  TAG = "coral1"   
)

lipid_class<-word(colnames(lip_data), 1)

all_lipid_colors<-rep('ligtgrey',ncol(lip_data))
names(all_lipid_colors)<-colnames(lip_data)
for (cl in names(lipid_colors)){
  ids<-which(lipid_class==cl)
  all_lipid_colors[ids]<-lipid_colors[cl]
  names(all_lipid_colors)[ids]<-cl
}

###########################
# Differential expression #
###########################

# "Healthy_CM","EAE.CM"
levs = c("Healthy_CM","EAE.CM")
take_ids<-which(info$condition %in% levs)
log_df<-log2(t(lip_data)[,take_ids])
group<-factor(condition[take_ids])
design<-model.matrix(~ 0 + group)
colnames(design)<-levels(group)
fit<-lmFit(log_df, design)
contrast.matrix<-makeContrasts(EAE.CM - Healthy_CM, levels = design)
fit2<-contrasts.fit(fit, contrast.matrix)
fit2<-eBayes(fit2)
results<-topTable(fit2, adjust = "fdr", number=Inf,sort.by = 'none')
res_healthy_eae_CM<-results

write.xlsx2(results, paste0("lipidomics/de_",paste0(levs,collapse='_'),".xlsx"), row.names= T, quote = F)

# Volcano

# lipids to label in volcano plot
ids<-which(res_healthy_eae_CM$adj.P.Val<0.05 & res_healthy_eae_CM$logFC>1)
select_lips<-rownames(results)[ids]

pdf(width=6,height=9,file=paste0('lipidomics/volcano_',paste0(levs,collapse='_'),'_colored.pdf'))
v<-EnhancedVolcano(results,
                   lab = rownames(results),
                   selectLab = select_lips,
                   labSize = 3,
                   drawConnectors = TRUE,
                   maxoverlapsConnectors=20,
                   colCustom = all_lipid_colors,
                   colConnectors = all_lipid_colors[ids],
                   labCol = all_lipid_colors[ids],
                   ylim=c(-0.5,7.5),
                   xlim=c(-2.5,4),
                   border='full',
                   gridlines.major = FALSE,
                   gridlines.minor = FALSE,
                   cutoffLineType = 'longdash',
                   cutoffLineCol = 'black',
                   cutoffLineWidth = 0.8,
                   typeConnectors = 'open',
                   arrowheads = FALSE,
                   x = 'logFC',
                   y = 'P.Value', pCutoff=0.05)
print(v)
dev.off()


# "EAE.CM_Cells","EAE.CM_L.ARG_Cells"
levs = c("EAE.CM_Cells","EAE.CM_L.ARG_Cells")
take_ids<-which(info$condition %in% levs)
log_df<-log2(t(lip_data)[,take_ids])
group<-factor(condition[take_ids])
design<-model.matrix(~ group)
colnames(design)<-levels(group)
fit<-lmFit(log_df, design)
contrast.matrix<-makeContrasts(EAE.CM_L.ARG_Cells - EAE.CM_Cells, levels = design)
fit2<-contrasts.fit(fit, contrast.matrix)
fit2<-eBayes(fit2)
results<-topTable(fit2, adjust = "fdr", number=Inf,sort.by = 'none')
res_argeae_eae_cells<-results

write.xlsx2(results, paste0("lipidomics/de_",paste0(levs,collapse='_'),".xlsx"), row.names= T, quote = F)

# Volcano

# lipids to label in volcano plot
ids<-which(res_healthy_eae_CM$adj.P.Val<0.05 & res_healthy_eae_CM$logFC>1 & res_argeae_eae_cells$adj.P.Val<0.05 & res_argeae_eae_cells$logFC<(-1) )
select_lips<-rownames(results)[ids]

pdf(width=6,height=9,file=paste0('lipidomics/volcano_',paste0(levs,collapse='_'),'_colored_labels_sig_healthyvEAE.pdf'))
v<-EnhancedVolcano(results,
                   lab = rownames(results),
                   selectLab = select_lips,
                   labSize = 3,
                   drawConnectors = TRUE,
                   maxoverlapsConnectors=20,
                   colCustom = all_lipid_colors,
                   colConnectors = all_lipid_colors[ids],
                   labCol = all_lipid_colors[ids],
                   ylim=c(-0.5,15),
                   xlim=c(-11,5),
                   border='full',
                   gridlines.major = FALSE,
                   gridlines.minor = FALSE,
                   cutoffLineType = 'longdash',
                   cutoffLineCol = 'black',
                   cutoffLineWidth = 0.8,
                   typeConnectors = 'open',
                   arrowheads = FALSE,
                   x = 'logFC',
                   y = 'P.Value', pCutoff=0.05)
print(v)
dev.off()


######################
# Heatmap & analysis #
#  per lipid class  #
######################

# Cells in three conditions
levs = c("Untreated_Cells","EAE.CM_Cells","EAE.CM_L.ARG_Cells")

# sum by lipid classes

take_ids<-which(condition %in% levs)
take_ids
unique_lipid_class<-unique(lipid_class)
unique_lipid_class

lip_sum<-matrix(NA,nrow=length(unique_lipid_class), ncol=length(take_ids))
colnames(lip_sum)<-info$`EPPI LID LABEL`[take_ids]
rownames(lip_sum)<-unique_lipid_class

for (i in 1:length(unique_lipid_class)){
  row_ids<-which(lipid_class==unique_lipid_class[i])
  lip_sum[i,]<-colSums(t(lip_data[take_ids,row_ids]), na.rm=T)
}

lip_sum
write.xlsx2(lip_sum, paste0("lipidomics/Cells_sum_by_lipid_class.xlsx"), row.names= T, quote = F)

# then log
log_lip_sum<-log2(lip_sum)
write.xlsx2(log_lip_sum, paste0("lipidomics/Cells_sum_by_lipid_class_log.xlsx"), row.names= T, quote = F)

# then z-score
library(pheatmap)
library(grid)

pdf("lipidomics/Heatmap_classSums_EAE_Cells.pdf", width=6,height=4)
p<-pheatmap(log_lip_sum,
            color = rev(colorRampPalette(c("#B2182B", "#FDDBC7", "#2166AC"))(n = 256)),
            border_color = "black",
            scale = "row",
            angle_col = 45,
            fontsize_row = 9,
            cluster_cols = F,
            cluster_rows = T,
            cellwidth = 15,
            # annotation_col =  colnames(log_lip_sum),
            # annotation_row = rownames(log_lip_sum),
            annotation_names_col = F
)
print(p)
dev.off()

# T-test
for (i in 1:ncol(log_lip_sum)){
  res<-t.test(log_lip_sum[i,5:8],log_lip_sum[i,9:12])
  print(rownames(log_lip_sum)[i])
  print(res$p.value)
}

# ANOVA + Tukey's HSD
group_labels <- factor(c(rep("Healthy", 4),rep("EAE", 4),rep("EAE_ARG", 4)))
tukey_results <- list()

for (i in 1:nrow(log_lip_sum)) {
  lipid_values <- as.numeric(log_lip_sum[i, 1:12])
  lipid_name <- rownames(log_lip_sum)[i]
  df_temp <- data.frame(
    Value = lipid_values,
    Group = group_labels
  )
  anova_res <- aov(Value ~ Group, data = df_temp)
  tukey_res <- TukeyHSD(anova_res)
  tukey_results[[lipid_name]] <- tukey_res$Group[,"p adj"]
}

final_stats <- do.call(rbind,tukey_results)
head(final_stats)
write.xlsx2(final_stats,paste0("lipidomics/Cells_sum_by_lipid_class_ANOVA+TukeyHSD.xlsx"),row.names= T,quote = F)


