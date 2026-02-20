setwd('~/GCH_Hsiesh/Mp2/github')
#----------------------- data loading -----------------------
meta <- read.csv('data/meta.csv',sep=',')
rownames(meta) <- meta$Sample_ID

vcf <- read.table('data/outgroup_ref_no_pre_post_M129.summary_of_snp_distribution.vcf', sep = '\t',comment.char = '#')
vcf <- vcf[,(1:ncol(vcf)-1)]
hh <- readLines('data/outgroup_ref_no_pre_post_M129.summary_of_snp_distribution.vcf', n = 4)[4]
colnames(vcf) <- hh <- strsplit(hh,"\\s+")[[1]]
snp_info <- vcf[,1:5]
snp <- vcf[,10:ncol(vcf)]
snp <- snp[,meta$Sample_ID]

#----------------------- Table S1 -----------------------
submeta <- meta[!is.na(meta$Subclude),]
submeta <- submeta[submeta$Subclude == 'subclade_pre',]
submeta$condition <- ifelse(submeta$Year > 2020,'post','pre')
subsnp <- snp[,rownames(submeta)]
table(submeta$condition)

stat_df <- data.frame()
for(i in 1:nrow(subsnp)){
  snpvec1 <- as.character(subsnp[i,submeta$condition == 'post'])
  snpvec0 <- as.character(subsnp[i,submeta$condition != 'post'])
  
  all_bases <- as.character(subsnp[i,])
  all_bases_clean <- all_bases[!all_bases %in% c("-", "N", "")]
  if(length(all_bases_clean) == 0){next}
  
  snpvec0_clean <- snpvec0[!snpvec0 %in% c("-", "N", "")]
  if(length(snpvec0_clean) == 0){next}
  
  snpvec1_clean <- snpvec1[!snpvec1 %in% c("-", "N", "")]
  if(length(snpvec1_clean) == 0){next}
  
  ref_base <- names(sort(table(all_bases_clean), decreasing = T))[1]
  rest_base <- all_bases[!all_bases %in% c(ref_base, "N", "")]
  rest_base <- paste0(unique(rest_base), collapse = ',')
  success_post <- sum(snpvec1 == ref_base, na.rm = TRUE)
  n_post       <- sum(snpvec1 %in% c("A","T","C","G","-"))
  
  success_pre  <- sum(snpvec0 == ref_base, na.rm = TRUE)
  n_pre        <- sum(snpvec0 %in% c("A","T","C","G","-"))
  
  test <- prop.test(
    x = c(success_post, success_pre),
    n = c(n_post,     n_pre),
    alternative = "two.sided",
  )
  
  rest_base <- paste0(strsplit(rest_base,',')[[1]], collapse = '&')
  test <- data.frame('site'=snp_info[i,2],'REF'=ref_base,'ALT'=rest_base,'G_prop'=test$estimate[1],'nG_prop'=test$estimate[2],'p'=test$p.value)
  
  stat_df <- rbind(stat_df, test)
}

stat_df <- stat_df[!is.na(stat_df$p),]
stat_df$padj <- p.adjust(stat_df$p, method = 'bonferroni')
stat_df <- stat_df[stat_df$padj < 1e-5,]
nrow(stat_df)

plot(stat_df$site, log(stat_df$padj)*-1)
plot(stat_df$G_prop, stat_df$nG_prop, xlim = c(0,1), ylim = c(0,1))
write.csv(stat_df,'subclade_pre.csv', row.names = F, quote = F)
#----------------------- Table S2 -----------------------
submeta <- meta[meta$ST == '3',]
submeta$Subclude[is.na(submeta$Subclude)] <- "other ST3"
subsnp <- snp[,rownames(submeta)]
table(submeta$Subclude)

stat_df <- data.frame()
for(cl in c('subclade_post_a','subclade_post_b','subclade_post_c','subclade_post_d','subclade_pre')){
  for(i in 1:nrow(subsnp)){
    snpvec1 <- as.character(subsnp[i,submeta$subclude == cl])
    snpvec0 <- as.character(subsnp[i,submeta$subclude != cl])
    
    all_bases <- as.character(subsnp[i,])
    all_bases_clean <- all_bases[!all_bases %in% c("-", "N", "")]
    if(length(all_bases_clean) == 0){next}
    
    snpvec0_clean <- snpvec0[!snpvec0 %in% c("-", "N", "")]
    if(length(snpvec0_clean) == 0){next}
    
    snpvec1_clean <- snpvec1[!snpvec1 %in% c("-", "N", "")]
    if(length(snpvec1_clean) == 0){next}
    
    ref_base <- names(sort(table(all_bases_clean), decreasing = T))[1]
    rest_base <- all_bases[!all_bases %in% c(ref_base, "N", "")]
    rest_base <- paste0(unique(rest_base), collapse = ',')
    success_1 <- sum(snpvec1 == ref_base, na.rm = TRUE)
    n_1       <- sum(snpvec1 %in% c("A","T","C","G","-"))
    
    success_0  <- sum(snpvec0 == ref_base, na.rm = TRUE)
    n_0        <- sum(snpvec0 %in% c("A","T","C","G","-"))
    
    test <- prop.test(
      x = c(success_1, success_0),
      n = c(n_1,     n_0),
      alternative = "two.sided",
    )
    
    rest_base <- paste0(strsplit(rest_base,',')[[1]], collapse = '&')
    test <- data.frame('clade'=cl,'site'=snp_info[i,2],'REF'=ref_base,'ALT'=rest_base,'clade_prop'=test$estimate[1],'others_prop'=test$estimate[2],'p'=test$p.value)
    stat_df <- rbind(stat_df, test)
  }
}
stat_df <- stat_df[!is.na(stat_df$p),]
stat_df$padj <- p.adjust(stat_df$p, method = 'bonferroni')
stat_df <- stat_df[stat_df$padj < 1e-5,]
nrow(stat_df)
plot(stat_df$site, log(stat_df$padj)*-1)
plot(stat_df$G_prop, stat_df$nG_prop, xlim = c(0,1), ylim = c(0,1))

write.csv(stat_df,'clade_specific_snp.csv', row.names = F, quote = F)