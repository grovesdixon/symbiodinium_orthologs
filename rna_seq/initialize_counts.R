#initialize_counts.R
rm(list=ls())
source('project_functions.R')


# SELECT COUNTS SET -------------------------------------------------------

source('rna_seq/data/cladeC_count_source.R') #set variables for cladeC
source('rna_seq/data/cladeD_count_source.R') #set variables for cladeD

# READ IN -----------------------------------------------------------------

#load counts based on the selected species above
counts1 = read_tsv(multicov_input,
                 trim_ws = TRUE) %>% 
  select(-start, -end) %>% 
  dplyr::rename(gene=chr)
colnames(counts1) = sub(to_remove, '', colnames(counts1))
colnames(counts1) = sub(paste('_', clade, sep=''), '', colnames(counts1))


# ORGANIZE COLDATA --------------------------------------------------------

#load trait data and ensure sample names match
coldata = read_csv('metadata/allSymbiodinium_coldata.csv')
samples = counts1 %>% 
  dplyr::select(-gene) %>% 
  colnames()
sum(samples %in% coldata$Run)==nrow(coldata)




# MERGE WITH ORTHOLOG GROUPS ----------------------------------------------

#read in the single copy orthologs
odat = read_tsv('ortholog_results/singleCopyOrthos.txt',
                col_names = c('orthogroup', 'gene')) %>% 
  mutate(gene = sapply(gene, function(x) strsplit(x, '.p')[[1]][1])) %>% 
  filter(gene %in% counts1$gene)

#read in the set of collapsed contigs (these are output from paraPrune.py as of 8-9-20 update)
cdat = read_tsv('ortholog_results/collapsed_contigs.tsv') %>% 
  mutate(kept = sapply(kept, function(x) strsplit(x, '.p')[[1]][1]),
         collapsed = sapply(collapsed, function(x) strsplit(x, '.p')[[1]][1]),
         keep = kept %in% counts1$gene & collapsed %in% counts1$gene) %>% 
  filter(keep,
         kept %in% odat$gene)
cdat

#add the collapsed contigs onto the single copy
collapsed_dat = cdat %>% 
  dplyr::select(orthogroup, collapsed) %>% 
  dplyr::rename(gene = collapsed)
odat2 = rbind(odat, collapsed_dat)


#do the merge
counts2 = counts1 %>% 
  inner_join(odat2, by = 'gene')
dim(counts2)

#check 
counts2 %>% 
  group_by(orthogroup) %>% 
  summarize(N=n()) %>% 
  filter(N==3)

#sum over collapsed
counts3 = counts2 %>% 
  select(-gene) %>% 
  pivot_longer(-orthogroup,
               names_to = 'sample',
               values_to = 'count') %>% 
  group_by(orthogroup, sample) %>% 
  summarize(orthoSum = sum(count)) %>% 
  pivot_wider(names_from = sample,
              values_from = orthoSum)



# OUTPUT ------------------------------------------------------------------

dim(counts3)
out_name = paste('rna_seq/data/', clade, '_ortholog_counts.Rdata', sep='')
counts = counts3
save(coldata, counts, clade, file=out_name)
