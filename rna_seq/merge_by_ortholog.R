#merge_by_ortholog.R
rm(list=ls())
library(DESeq2)

# CHOOSE DATASETS ---------------------------------------------------------


dat1_file = 'rna_seq/data/cladeC_ortholog_counts.Rdata'
dat2_file = 'rna_seq/data/cladeD_ortholog_counts.Rdata'


# LOAD AND COMBINE --------------------------------------------------------------------

#LOAD
load_clade = function(dat_file){
  ll=load(dat_file)
  colnames(counts) = paste(colnames(counts), clade, sep='_')
  colnames(counts)[1]<-'orthogroup'
  coldata$Run = paste(coldata$Run, clade, sep='_')
  res = list(counts,
             coldata,
             clade)
  names(res)=c('counts',
               'coldata',
               'clade')
  return(res)
}

#load data1
dat1 = load_clade(dat1_file)
counts1 = dat1[['counts']]
coldata1 = dat1[['coldata']]
clade1 = dat1[['clade']]

#load data2
dat2 = load_clade(dat2_file)
counts2 = dat2[['counts']]
coldata2 = dat2[['coldata']]
clade2 = dat2[['clade']]



#COMBINE INTO SINGLE DFS
counts0 = inner_join(counts1, counts2, by = 'orthogroup') %>% 
  column_to_rownames('orthogroup') %>% 
  data.frame()
coldata = rbind(coldata1, coldata2)


# FILTER THE COMBINED SET -------------------------------------------------

#FILTER LOW COVERAGE SAMPLES
sample_cutoff = 1e5
sample_totals = apply(counts0, 2, sum)
plot(density(sample_totals))
low_samples = colnames(counts0)[sample_totals < sample_cutoff]
length(low_samples)
filt1_counts = counts0[!colnames(counts0) %in% low_samples]
dim(sample_filt)
new_totals = apply(filt1_counts, 2, sum)
lines(density(new_totals), col='blue')



#FILTER LOW COVERAGE GENES
gene_cutoff=3
means=apply(filt1_counts,1,mean)
table(means>gene_cutoff)
filt2_counts = filt1_counts[means > gene_cutoff,]
dim(filt2_counts)
filt_coldata = data.frame(coldata)
rownames(filt_coldata) = coldata$Run
filt_coldata = filt_coldata[colnames(filt2_counts),]


#GET RLD

ddsHTSeq<-DESeqDataSetFromMatrix(filt2_counts,
                                 colData = filt_coldata,
                                 design = formula(~1))
dds = DESeq(ddsHTSeq)
rld = vst(dds)
rld_df=assay(rld)
colnames(rld_df) = colnames(filt2_counts)


# CHECK CORRELATION -------------------------------------------------------

#GET RLD FOR EACH
prep_for_cor = function(rld_df, clade){
  rld_sub = rld_df[,grep(clade, colnames(rld_df))] %>% 
    as_tibble() %>% 
    mutate(orthogroup = rownames(rld_df))
  colnames(rld_sub) = sub(paste('_', clade, sep=''), '', colnames(rld_sub))
  rld_long = rld_sub %>% 
    pivot_longer(-orthogroup,
                 names_to = 'sample',
                 values_to = 'count')
  return(rld_long)
}

rld1 = prep_for_cor(rld_df, clade1)
rld2 = prep_for_cor(rld_df, clade2)
test_df = rld1 %>% 
  inner_join(rld2,
             by=c('orthogroup','sample'))




get_long = function(counts, clade){
  long = counts %>% 
    pivot_longer(-orthogroup,
                 names_to = 'sample',
                 values_to = 'count') %>% 
    mutate(sample = sub(paste('_', clade, sep=''), '', sample))
  return(long)
}

long1 = get_long(counts1, clade1)
long2 = get_long(counts2, clade2)
test_dat = inner_join(long1, long2, by = c('orthogroup', 'sample'))


test_dat %>% 
  ggplot(aes(x=count.x,
             y=count.y)) +
  geom_point()

