info <- read_excel("lipidomics/EAECM Lipidomics Sample List.xlsx")
lipids_df <- read_excel("lipidomics/Report_M085_LIPID01_20251117.xlsx",na = "<LOD")
colnames(lipids_df)
class(lipids_df)
lipids_df<-(as.data.frame(lipids_df))

lipids_df[1:5,1:10]

lip_data<-lipids_df[,-c(1:4)]
lip_data[1:5,1:5]

table(lipids_df$`Reported Unit`)
lipids_df[,c('Reported Unit','Sample Description')]

lipid_class <- word(colnames(lip_data), 1)
table(lipid_class)

lipids_df$`Sample ID`
info #(NO. from info ~ rownames(lip_data))

rat<-1
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
colSums(is.na(lip_data_repl))
which(colSums(is.na(lip_data_repl))!=0)
lip_data<-lip_data[,-which(colnames(lip_data) %in% (c('PA 36:4','PA 37:4')))]
lip_data_repl<-lip_data_repl[,-which(colnames(lip_data_repl) %in% (c('PA 36:4','PA 37:4')))]
lipids_df<-lipids_df[,-which(colnames(lipids_df) %in% (c('PA 36:4','PA 37:4')))]

saveRDS(lip_data_repl,'lip_data_repl1.Rds')

log_df<-log2(lip_data_repl)
head(log_df)





