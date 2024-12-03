setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(readr)
library(lattice)
library(colorspace)
library(quantmod)
library(stats)

# Rearrange TRGT output into tables by region
tmp.trgt = read.delim('_trgt_output_pathogenic/198_pathogenic_trgt.vcf', stringsAsFactors = F, skip = 211, header = T)
tmp.trgt$gt = sapply(strsplit(tmp.trgt$X198.GRCh38.deepvariant.haplotagged, ':'), '[', 5)
tmp.trgt$STRUCTURE = sapply(strsplit(tmp.trgt$INFO, '='), '[', 5)

tmp.trgt$motif = sapply(strsplit(tmp.trgt$INFO, '='), '[', 4)
tmp.trgt$motif = sapply(strsplit(tmp.trgt$motif, ';'), '[', 1)
# tmp.trgt$motif.len = nchar(sapply(strsplit(tmp.trgt$motif, ','), '[', 1))

tmp.trgt$gt.unit = paste0(round(as.numeric(sapply(strsplit(tmp.trgt$gt, ','), '[', 1)), 0), '|', 
                          round(as.numeric(sapply(strsplit(tmp.trgt$gt, ','), '[', 2)), 0))
tmp.trgt$gt.unit[tmp.trgt$gt == '.'] = '.|.'

tmp.trgt$gene = sapply(strsplit(tmp.trgt$INFO, '='), '[', 2)
tmp.trgt$gene = sapply(strsplit(tmp.trgt$gene, ';'), '[', 1)

tmp.trgt$end = sapply(strsplit(tmp.trgt$INFO, '='), '[', 3)
tmp.trgt$end = sapply(strsplit(tmp.trgt$end, ';'), '[', 1)

tmp.trgt$id = paste(tmp.trgt$X.CHROM, tmp.trgt$POS, tmp.trgt$end, tmp.trgt$motif,sep = '.')
tmp.trgt$gene.id = paste(tmp.trgt$gene, tmp.trgt$motif,sep = '.')

for (gene.id in tmp.trgt$gene.id){
  tmp.table = data.frame()
  for (file in list.files('_trgt_output_custom/')){
    tmp.file = read.delim(paste0('_trgt_output_custom/', file), stringsAsFactors = F, skip = 211, header = T)
    colnames(tmp.file)[10] = 'deepvariant.haplotagged'
    # tmp.file$gt = sapply(strsplit(tmp.file$deepvariant.haplotagged, ':'), '[', 2)
    tmp.file$gt = sapply(strsplit(tmp.file$deepvariant.haplotagged, ':'), '[', 5)
    
    tmp.file$gene = sapply(strsplit(tmp.file$INFO, '='), '[', 2)
    tmp.file$gene = sapply(strsplit(tmp.file$gene, ';'), '[', 1)
    
    tmp.file$motif = sapply(strsplit(tmp.file$INFO, '='), '[', 4)
    tmp.file$motif = sapply(strsplit(tmp.file$motif, ';'), '[', 1)
    tmp.file$motif.len = nchar(tmp.file$motif)
    
    # tmp.file$gt.1 = round(as.numeric(sapply(strsplit(tmp.file$gt, ','), '[', 1))/tmp.file$motif.len, 0)
    # tmp.file$gt.2 = round(as.numeric(sapply(strsplit(tmp.file$gt, ','), '[', 2))/tmp.file$motif.len, 0)
    
    tmp.file$gt.1 = round(as.numeric(sapply(strsplit(tmp.file$gt, ','), '[', 1)), 0)
    tmp.file$gt.2 = round(as.numeric(sapply(strsplit(tmp.file$gt, ','), '[', 2)), 0)
    
    for (i in 1:nrow(tmp.file)){
      if (!is.na(tmp.file$gt.1[i]) & !is.na(tmp.file$gt.2[i])){
        if (tmp.file$gt.1[i] != tmp.file$gt.2[i]){
          tmp.file$allele.1[i] = min(tmp.file$gt.1[i], tmp.file$gt.2[i])
          tmp.file$allele.2[i] = max(tmp.file$gt.1[i], tmp.file$gt.2[i])
        }else{
          tmp.file$allele.1[i] = tmp.file$gt.1[i]
          tmp.file$allele.2[i] = tmp.file$gt.2[i]
        }
      }
    }
    tmp.file$gt.unit = paste0(tmp.file$allele.1, '|', tmp.file$allele.2)
  
    tmp.file$gt.unit[tmp.file$gt == '.'] = '.|.'
    # now rearrange them into tables by region 
    tmp.file$end = sapply(strsplit(tmp.file$INFO, '='), '[', 3)
    tmp.file$end = sapply(strsplit(tmp.file$end, ';'), '[', 1)
    
    tmp.file$id = paste(tmp.file$X.CHROM, tmp.file$POS, tmp.file$end, tmp.file$motif,sep = '.')
    tmp.file$gene.id = paste(tmp.file$gene, tmp.file$motif,sep = '.')
    
    tmp.file$sample = strsplit(file, '_')[[1]][1]
    tmp.table = rbind(tmp.table, tmp.file[tmp.file$gene.id == gene.id,])
    colnames(tmp.table)[c(21,12,13,23,19)]
    colnames(tmp.file[tmp.file$gene.id == gene.id,])
  }
  write.table(tmp.table[c(21,12,13,23,19)], paste0('_trgt_custom_named_UPD/', tmp.table$id[1], '.txt'), col.names = T, row.names = F, quote = F, sep = '\t')
  
}
colnames(tmp.table)[c(17,13,14,19,12)]
# combine two tables (Eh and TRGT) - plot distribution
# then check the alternative, if needed - calculate them separately

sum.table = data.frame()

for (i in 1:length(list.files('_trgt_custom_named_UPD/'))){ # 69
  tmp = read.delim(paste0('_trgt_custom_named_UPD/', list.files('_trgt_custom_named_UPD/')[i]), stringsAsFactors = F)
  tmp = tmp[tmp$sample %in% tof.meta$WGS.ID,]# 48 -> 36
  region.tmp = sub(".[^.]+$", "", sub(".[^.]+$", "", list.files('_trgt_custom_named_UPD//')[i]))
  gene.tmp = tmp$gene[1]
  motif.tmp = tmp$motif[1]
  
  tmp$all1 = sapply(strsplit(tmp$gt.unit, '\\|'), '[', 1)
  tmp$all2 = sapply(strsplit(tmp$gt.unit, '\\|'), '[', 2)
  tmp = tmp[ ! (tmp$all1 == '.' & tmp$all2 == '.'), ]
  tmp$all1 = as.numeric(tmp$all1)
  tmp$all2 = as.numeric(tmp$all2)
  
  for (j in 1:nrow(tmp)){
    tmp$sex[j] = tof.meta.samples$Gender[tof.meta.samples$WGS.ID == tmp$sample[j]]
  }
  for (j in 1:nrow(tmp)){
    if (strsplit(region.tmp, '\\.')[[1]][1] == 'chrX'){
      tmp$allele1[j] = ifelse(tmp$sex[j] == 'M', mean(tmp$all1[j], tmp$all2[j]), tmp$all1)
      tmp$allele2[j] = ifelse(tmp$sex[j] == 'M', '.', tmp$all2)
    }else{
      tmp$allele1[j] = tmp$all1[j]
      tmp$allele2[j] = tmp$all2[j]
    }
    
  }
  tmp$all1 = tmp$allele1
  tmp$all2 = tmp$allele2
  
  tmp$allele2 = ifelse( tmp$all2 == '.', tmp$all1, tmp$all2)# homozyg
  # legend('topleft', legend = category, bty = "n", inset = c(-0.05, -0.1), xpd = TRUE, col = 'darkgrey')
  # legend('topleft', legend = gene.tmp, bty = "n", inset = c(-0.1, -0.1), xpd = TRUE, col = 'darkgrey')
  
  eh.tof.tmp = eh.tof[,c(sort(colnames(eh.tof))[i], 'Sample')]
  tof.region = colnames(eh.tof.tmp)[1]
  tof.region.tmp = sub(".[^.]+$", "", tof.region)
  
  eh.tof.tmp = eh.tof.tmp[eh.tof.tmp$Sample %in% tmp$sample,]# 36 -> 30
  
  eh.tof.tmp$all1 = sapply(strsplit(eh.tof.tmp[,1], '\\|'), '[', 1)
  eh.tof.tmp$all2 = sapply(strsplit(eh.tof.tmp[,1], '\\|'), '[', 2)
  for (j in 1:nrow(eh.tof.tmp)){
    eh.tof.tmp$sex[j] = tof.meta.samples$Gender[tof.meta.samples$WGS.ID == eh.tof.tmp$Sample[j]]
  }
  eh.tof.tmp = eh.tof.tmp[ ! (eh.tof.tmp$all1 == '.' & eh.tof.tmp$all2 == '.'), ]# 29
  # here go back to TRGT and keep only these samples
  tmp = tmp[tmp$sample %in% eh.tof.tmp$Sample,]
  tmp.freq = as.data.frame(table(tmp[,9:10]), stringsAsFactors = F)#162
  tmp.freq.nonzero = tmp.freq[tmp.freq$Freq != 0,]#35
  
  tmp.freq.nonzero$size = 1
  tmp.freq.nonzero$size[tmp.freq.nonzero$Freq <= 10] <- 1.75
  tmp.freq.nonzero$size[tmp.freq.nonzero$Freq <= 5] <- 1.5
  tmp.freq.nonzero$size[tmp.freq.nonzero$Freq <= 2] <- 1.25
  tmp.freq.nonzero$size[tmp.freq.nonzero$Freq == 1] <- 1
  
  eh.tof.tmp$allele1 = eh.tof.tmp$all1
  eh.tof.tmp$allele2 = ifelse( is.na(eh.tof.tmp$all2), eh.tof.tmp$all1, eh.tof.tmp$all2)# homozyg
  
  eh.tof.tmp$allele1 = as.numeric(eh.tof.tmp$allele1)
  eh.tof.tmp$allele2 = as.numeric(eh.tof.tmp$allele2)
    
  eh.tmp.freq = as.data.frame(table(eh.tof.tmp[,6:7]), stringsAsFactors = F)#162
  eh.tmp.freq.nonzero = eh.tmp.freq[eh.tmp.freq$Freq != 0,]#35
  
  eh.tmp.freq.nonzero$size = 1
  eh.tmp.freq.nonzero$size[eh.tmp.freq.nonzero$Freq <= 10] <- 1.75
  eh.tmp.freq.nonzero$size[eh.tmp.freq.nonzero$Freq <= 5] <- 1.5
  eh.tmp.freq.nonzero$size[eh.tmp.freq.nonzero$Freq <= 2] <- 1.25
  eh.tmp.freq.nonzero$size[eh.tmp.freq.nonzero$Freq == 1] <- 1
  
  
  pdf(paste0('~/Documents/_Ryan/_Seema_Myocardium/_PacBio/output_plot_eh_trgt_v7/', gene.tmp, '_', motif.tmp, '_', '_v7.pdf'), w = 10, h = 5)
  layout(matrix(c(1,2,3,4,5,6,7,8,9), nrow = 3, ncol = 3),#  byrow = TRUE 
         widths=c(rep(1,9)), heights=c(rep(c(1.5,1,1),3)))
  par(mar = c(4,4,3,6))
  
  xmin = min(as.numeric(tmp.freq.nonzero[,1]), as.numeric(eh.tmp.freq.nonzero[,1]))
  xmax = max(as.numeric(tmp.freq.nonzero[,1]), as.numeric(eh.tmp.freq.nonzero[,1]))
  
  ymin = min(as.numeric(tmp.freq.nonzero[,2]), as.numeric(eh.tmp.freq.nonzero[,2]))
  ymax = max(as.numeric(tmp.freq.nonzero[,2]), as.numeric(eh.tmp.freq.nonzero[,2]))

  # Plot #1 - allele size distribution for LONG read sequencing
  plot(tmp.freq.nonzero[,1], tmp.freq.nonzero[,2], col = 'black', bg = alpha('maroon3', 0.7), pch = 21, cex = tmp.freq.nonzero$size,
       xlab = 'Smaller allele size', ylab = 'Larger allele size', main = 'LONG', # paste0('LONG: ', gene.tmp)
       xlim = c(xmin, xmax), ylim = c(ymin, ymax), yaxt = 'n')
  axis(2, las = 2)
  legend('topright', legend = region.tmp, bty = "n", inset = c(0, -0.15), xpd = TRUE)
  legend('topleft', legend = strsplit(gene.tmp, ',')[[1]][1], bty = "n", inset = c(-0.07, -0.15), xpd = TRUE)
  legend('topright', legend = paste0(c('1', '2','5','10')),col='black',pt.bg = alpha('maroon3', 0.7),
         pch=21,pt.cex=c(1,1.25,1.5,1.75),cex = 0.8,title = 'Frequency', xpd = TRUE, inset = c(-0.25,0), bty = 'n')
  # legend('bottomright', legend = paste0(c('control', 'GENCOV')),col='black',
  #        pt.bg = c('grey', 'orange'),pch=21,pt.cex=1,
  #        cex = 0.8,title = '', xpd = TRUE, inset = c(-0.04,0.1), bty = 'n')
  
  # Plot #2 - smaller allele size histogram for LONG read sequencing
  hist(as.numeric(tmp$all1), breaks = 30, col = alpha('maroon3', 0.8), main = 'Smaller allele size', xlab = '', # breaks = 30, 
       xlim = c(min(as.numeric(tmp$allele1), as.numeric(tmp$allele2)), max(as.numeric(tmp$allele1), as.numeric(tmp$allele2))), yaxt = 'n')
  axis(2, las = 2)
  legend('topright', legend = paste0('n = ', length(tmp$all1[tmp$all1 != '.'])), bty = "n", inset = c(0, -0.1), xpd = TRUE)

  # Plot #3 - larger allele size histogram for LONG read sequencing
  hist(as.numeric(tmp$all2[tmp$all2 != '.']), breaks = 30, col = alpha('maroon3', 0.8), main = 'Larger allele size', xlab = '',# breaks = 30, 
       xlim = c(min(as.numeric(tmp$allele1), as.numeric(tmp$allele2)), max(as.numeric(tmp$allele1), as.numeric(tmp$allele2))), yaxt = 'n')
  axis(2, las = 2)
  legend('topright', legend = paste0('n = ',length(tmp$all2[tmp$all2 != '.'])), bty = "n", inset = c(0, -0.1), xpd = TRUE)
  
  # Plot #4 - allele size distribution for SHORT read sequencing
  plot(eh.tmp.freq.nonzero[,1], eh.tmp.freq.nonzero[,2], col = 'black', bg = alpha('black', 0.5), pch = 21, cex = tmp.freq.nonzero$size,
       xlab = 'Smaller allele size', ylab = 'Larger allele size', main = 'SHORT',# paste0('SHORT: ', gene.tmp)
       xlim = c(xmin, xmax), ylim = c(ymin, ymax), yaxt = 'n')
  axis(2, las = 2)
  legend('topright', legend = tof.region.tmp, bty = "n", inset = c(0, -0.15), xpd = TRUE)
  legend('topleft', legend = strsplit(gene.tmp, ',')[[1]][1], bty = "n", inset = c(-0.07, -0.15), xpd = TRUE)
  legend('topright', legend = paste0(c('1', '2','5','10')),col='black',pt.bg = alpha('black', 0.5),
         pch=21,pt.cex=c(1,1.25,1.5,1.75),cex = 0.8,title = 'Frequency', xpd = TRUE, inset = c(-0.25,0), bty = 'n')

  # Plot #5 - smaller allele size histogram for SHORT read sequencing
  hist(as.numeric(eh.tof.tmp$all1), breaks = 30, col = alpha('black', 0.5), main = 'Smaller allele size', xlab = '',
       xlim = c(min(as.numeric(eh.tof.tmp$allele1), as.numeric(eh.tof.tmp$allele2)), max(as.numeric(eh.tof.tmp$allele1), as.numeric(eh.tof.tmp$allele2))), yaxt = 'n')
  axis(2, las = 2)
  legend('topright', legend = paste0('n = ', length(eh.tof.tmp$all1)), bty = "n", inset = c(0, -0.1), xpd = TRUE)
  
  # Plot #6 - larger allele size histogram for SHORT read sequencing
  hist(as.numeric(eh.tof.tmp$all2), breaks = 30, col = alpha('black', 0.5), main = 'Larger allele size', xlab = '',
       xlim = c(min(as.numeric(eh.tof.tmp$allele1), as.numeric(eh.tof.tmp$allele2)), max(as.numeric(eh.tof.tmp$allele1), as.numeric(eh.tof.tmp$allele2))), yaxt = 'n')
  axis(2, las = 2)
  legend('topright', legend = paste0('n = ', length(na.omit(eh.tof.tmp$all2))), bty = "n", inset = c(0, -0.1), xpd = TRUE)

  # Plot #7 - allele size ditribution for both LONG and SHORT read sequencing
  plot(eh.tmp.freq.nonzero[,1], eh.tmp.freq.nonzero[,2], col = 'black', bg = alpha('black', 0.5), pch = 21, cex = tmp.freq.nonzero$size,
       xlab = 'Smaller allele size', ylab = 'Larger allele size', main = 'LONG and SHORT',# paste0('LONG and SHORT: ', gene.tmp)
       xlim = c(xmin, xmax), ylim = c(ymin, ymax), yaxt = 'n')
  axis(2, las = 2)
  legend('topleft', legend = strsplit(gene.tmp, ',')[[1]][1], bty = "n", inset = c(-0.07, -0.15), xpd = TRUE)
  points(tmp.freq.nonzero[,1], tmp.freq.nonzero[,2], col = 'black', bg = alpha('maroon3', 0.8), pch = 21, cex = tmp.freq.nonzero$size)
  legend('topright', legend = paste0(c('1', '2','5','10')),col='black',pt.bg = alpha('black', 0.5),
         pch=21,pt.cex=c(1,1.25,1.5,1.75),cex = 0.8,title = 'Frequency', xpd = TRUE, inset = c(-0.25,0), bty = 'n')
  legend('bottomright', legend = paste0(c('SHORT', 'LONG')),col='black',
         pt.bg = c(alpha('black', 0.5), alpha('maroon3', 0.8)),pch=21,pt.cex=1.75,
         cex = 0.8,title = '', xpd = TRUE, inset = c(-0.25,0.2), bty = 'n')
  
  tmp.plot = tmp[tmp$sample %in% eh.tof.tmp$Sample,]
  eh.tof.tmp.plot = eh.tof.tmp[eh.tof.tmp$Sample %in% tmp$sample,]
  colnames(eh.tof.tmp.plot)[2] = 'sample'
  
  both.plot = merge(tmp.plot, eh.tof.tmp.plot, by = 'sample')
  
  # Plot #8 - correlation of a smaller allele size for LONG and SHORT read sequencing
  par(mar = c(4,4,1,15), mgp = c(2.2,1,0))
  plot(both.plot$all1.x, both.plot$all1.y, col = 'black', bg = alpha('darkgrey', 0.5), pch = 21, yaxt = 'n',
       xlab = 'LONG', ylab = 'SHORT', main = '')
  axis(2, las = 2)
  if (sd(as.numeric(both.plot$all1.x)) != 0){
    if (sd(as.numeric(both.plot$all1.y)) != 0){
      # sm.test = cor.test(as.numeric(both.plot$allele1.x), as.numeric(both.plot$allele1.y))
      sm.test = cor.test(as.numeric(both.plot$all1.x), as.numeric(both.plot$all1.y))
      legend('topright', legend = paste0('Comparison of smaller allele size\nin short read vs long read data:'), 
             bty = "n", inset = c(-2.2, -0.1), xpd = TRUE)
      legend('topright', legend = paste0('cor = ', format(sm.test$estimate, digits = 2)), 
             bty = "n", inset = c(-1, 0.3), xpd = TRUE, col = 'darkgrey')
      legend('topright', legend = paste0('p-val = ', format(sm.test$p.value, scientific = T, digits = 2)), 
             bty = "n", inset = c(-1.3, 0.5), xpd = TRUE, col = 'darkgrey')
      legend('topright', legend = paste0('n = ', nrow(both.plot)), 
             bty = "n", inset = c(-1, 0.7), xpd = TRUE, col = 'darkgrey')
      if(!is.na(sm.test$p.value)){
        # abline(lm(both.plot$allele1.y ~ as.numeric(both.plot$allele1.x)), lty = 2)
        abline(lm(both.plot$all1.y ~ as.numeric(both.plot$all1.x)), lty = 2)
        
      }
    }
  }
  
  # Plot #9 - correlation of a larger allele size for LONG and SHORT read sequencing
  plot(both.plot$all2.x, both.plot$all2.y, col = 'black', bg = alpha('darkgrey', 0.5), pch = 21, yaxt = 'n',
       xlab = 'LONG', ylab = 'SHORT', main = '')
  axis(2,las = 2)
  # cor.test(as.numeric(both.plot$allele1), as.numeric(both.plot$med_1))
  if (sd(as.numeric(both.plot$all2.x[both.plot$all2.x != '.'])) != 0){
    if (sd(as.numeric(both.plot$all2.y[!is.na(both.plot$all2.y)])) != 0){
      lg.test = cor.test(as.numeric(both.plot$allele2.x), as.numeric(both.plot$allele2.y))
      legend('topright', legend = paste0('Comparison of larger allele size\nin short read vs long read data:'), 
             bty = "n", inset = c(-2.1, -0.1), xpd = TRUE)
      legend('topright', legend = paste0('cor = ', format(lg.test$estimate, digits = 2)), 
             bty = "n", inset = c(-1, 0.3), xpd = TRUE, col = 'darkgrey')
      legend('topright', legend = paste0('p-val = ', format(lg.test$p.value, scientific = T, digits = 2)), 
             bty = "n", inset = c(-1.3, 0.5), xpd = TRUE, col = 'darkgrey')
      legend('topright', legend = paste0('n = ', nrow(both.plot)), 
             bty = "n", inset = c(-1, 0.7), xpd = TRUE, col = 'darkgrey')
      if(!is.na(lg.test$p.value)){
        abline(lm(as.numeric(both.plot$all2.y[!is.na(both.plot$all2.y)]) ~ as.numeric(both.plot$all2.x[both.plot$all2.x != '.'])), lty = 2)
      }
    }
  }
  
  dev.off()
  sum.table = rbind(sum.table, data.frame('long.varid' = tmp$id[1], 
                                          'short.varid' = colnames(eh.tof.tmp)[1],
                                          'gene' = tmp$gene[1], 'motif' = tmp$motif[1],
                                          'long.all.1' = paste0(both.plot$all1.x, collapse = ','),
                                          'long.all.2' = paste0(both.plot$all2.x, collapse = ','),
                                          'sh.all.1' = paste0(both.plot$all1.y, collapse = ','),
                                          'sh.all.2' = paste0(both.plot$all2.y, collapse = ','),
                                          'n.samples' = nrow(both.plot),
                                          'all.1.cor' = sm.test$estimate,
                                          'all.1.pval' = sm.test$p.value,
                                          'all.2.cor' = lg.test$estimate,
                                          'all.2.pval' =  lg.test$p.value))
}
write.table(sum.table.ref, "summary.table.txt", sep="\t", row.names = F, col.names=T, quote=F)

### Plot for all regions
sum.table = read.delim('summary.table.txt', stringsAsFactors = F)
sum.table$start = sapply(strsplit(sum.table$short.varid, '\\.'), '[', 2)
sum.table$end = sapply(strsplit(sum.table$short.varid, '\\.'), '[', 3)
sum.table$length = as.numeric(sum.table$end) - as.numeric(sum.table$start)
sum.table$chr = sub('chr', '', sapply(strsplit(sum.table$short.varid, '\\.'), '[', 1))

for (i in 1:nrow(sum.table)){
  sum.table$num[i] = length(sum.table$gene[sum.table$gene == sum.table$gene[i]])
}

sum.table.alt = sum.table[sum.table$num == 2,]
sum.table.alt[c(1,6,7,9,10,11,13,14,18,19),c(3,4)]
sum.table.ref = rbind(sum.table[sum.table$num == 1,], sum.table.alt[c(1,6,7,9,10,11,13,14,18,19),])
# sum.table.ref = sum.table[sum.table$num == 1,]

sum.table.ref$nmotif = nchar(sum.table.ref$motif)
sum.table.ref$len.units = sum.table.ref$length/sum.table.ref$nmotif
sum.table.ref$motif.plot = sum.table.ref$motif
sum.table.ref$motif.plot[sum.table.ref$motif == 'CGG'] = 'CCG'
sum.table.ref$motif.plot[sum.table.ref$motif == 'CGC'] = 'CCG'
sum.table.ref$motif.plot[sum.table.ref$motif == 'GCG'] = 'CCG'
sum.table.ref$motif.plot[sum.table.ref$motif == 'GCC'] = 'CCG'
sum.table.ref$motif.plot[sum.table.ref$motif == 'GGC'] = 'CCG'

sum.table.ref$motif.plot[sum.table.ref$motif == 'CTG'] = 'CAG'
sum.table.ref$motif.plot[sum.table.ref$motif == 'GCA'] = 'CAG'
sum.table.ref$motif.plot[sum.table.ref$motif == 'GTC'] = 'CAG'

sum.table.ref$motif.plot[sum.table.ref$motif == 'TAAAA'] = "AAAAT"
sum.table.ref$motif.plot[sum.table.ref$motif == 'AAAAT'] = "AAAAT"
sum.table.ref$motif.plot[sum.table.ref$motif == 'TTTTA'] = "AAAAT"

genes.1 = c('SAMD12', 'XYLT', 'TNRC6A', 'AFF2', 'MARCH6', 'BEAN1')# 1 down
genes.2 = c('CSTB', 'RAPGEF2', 'HOXA13', 'FOXL2', 'EIF4A3', 'TBP', 'DAB1', 'STARD7', 'RUNX2', 'KCNN3')# 2 left
genes.3 = c('HRAS', 'PHOX2B', 'RFC1', 'RELN', 'YEATS2', 'NOTCH2NLC', 'JPH3')# 3 up
genes.4 = c('CNBP')# 4 right

motifs = as.data.frame(table(sum.table.ref$motif.plot))
colnames(motifs) = c('motif', 'frequency')
motifs$motif = as.character(motifs$motif)
motifs$gc = str_count(motifs$motif, 'G|C')/nchar(motifs$motif)

motifs = motifs[order(motifs$gc, decreasing = F),]
motifs$color <- diverge_hcl(nrow(motifs), palette = "Blue-Red 2")

for (i in 1:nrow(sum.table.ref)){
  sum.table.ref$col.plot[i] = motifs$color[motifs$motif == sum.table.ref$motif.plot[i]]
}

pdf(paste0('cor_larger_smaller_59_trgt_eh_COLOR_LEGEND_B-R2_GC.pdf'), w = 8, h = 4.2)

par(mfrow=c(1,1), mar = c(4,4,1,17), mgp = c(3,1,0)) 
plot(sum.table.ref$all.1.cor, sum.table.ref$all.2.cor, pch = 21, col = 'black', bg = alpha(sum.table.ref$col.plot, 1), yaxt = 'n', # bg = sum.table.ref$color,
     xlab = 'Correlation on a smaller allele', ylab = 'Correlation on a larger allele', cex = sum.table.ref$length*0.01+1,# cex = sum.table.ref$nmotif*0.2+1,
     xlim = c(-0.15,1.1))# , ylim = c(0,1)
abline(h = 0.5, lty = 2)
abline(v = 0.5, lty = 2)
axis(2, las = 2)

text(sum.table.ref$all.1.cor, sum.table.ref$all.2.cor, #
     labels = ifelse(sum.table.ref$gene %in% genes.1, sum.table.ref$gene, ''), pos = 1, cex = 0.8)
text(sum.table.ref$all.1.cor, sum.table.ref$all.2.cor, #
     labels = ifelse(sum.table.ref$gene %in% genes.2, sum.table.ref$gene, ''), pos = 2, cex = 0.8)
text(sum.table.ref$all.1.cor, sum.table.ref$all.2.cor, #
     labels = ifelse(sum.table.ref$gene %in% genes.3, sum.table.ref$gene, ''), pos = 3, cex = 0.8)
text(sum.table.ref$all.1.cor, sum.table.ref$all.2.cor, #
     labels = ifelse(sum.table.ref$gene %in% genes.4, sum.table.ref$gene, ''), pos = 4, cex = 0.8)

legend('topright', legend = paste0(c('10', '50','100', ' 150')),col='black',pt.bg = alpha('grey', 0.5),
       pch=21,pt.cex=c(10,50,100,150)*0.01+1,cex = 0.8,title = 'Repeat size', xpd = TRUE, inset = c(-0.25,0), bty = 'n')
legend('bottomright', legend = motifs$motif,col='black',
       pt.bg = motifs$color,pch=21,pt.cex=1.75, cex = 0.8,title = '', xpd = TRUE, inset = c(-0.85,0), bty = 'n')

dev.off()
