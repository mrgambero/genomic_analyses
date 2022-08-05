#load dependencies or install 
library("tidyverse") #for various
library("taxizedb") #for taxonomic infos
library("ggplot2")#for graphs
library("ggtree")#for phylogenetic trees
library("doParallel") #for parallel computing 
library("foreach") #for improved compiler
library("DECIPHER") #For seqeunces allignments
library("phangorn") #For phylo tree calculations
library("seqinr") #for fasta operations
library("Biostrings") #for fasta operations 
library("treeio") # for phylogenetic trees
library("seqRFLP")#for DNA operations
library("rtracklayer") #for reading gff
library("gggenes")#For arrowplots

#### ### ### ### ### ### ### ### ### ### ### ### ##
###Code used for the data analysis of the genomic data presented in the paper ### insert title
###the code Performs the following things:
###1)Reconstruct the whole taxonomy base on species names
###2)Elaborate the hmm results from Dunivin et al 2019
###3)
#### ### ### ### ### ### ### ### ### ### ### ### ##

#Import the names of the isolates ---------------------------------------------------------------------------------------------------
#First obtain a list of all the names of the strains scanned, I used the list from the hmm-search results (first method I actually run).
#read in hmnsearch results (results have to be all in the same folder)
setwd(dir = "/Users/gabri/OneDrive - University of Arizona/Arsenic Project/Second_part_genomes/January_results/total_scanned/")
names=list.files(pattern="*.tbl.txt")
data <- do.call(rbind, lapply(names, function(X) { 
data.frame(id = basename(X), read.table(X, fill = TRUE)[,1:22])})) #due to name differences
#remove unnecessary columns
data <- data %>%
  mutate(id = gsub(".tbl.txt", "", id))  %>%
  select(-c(V2, V5))
#add column names
#known from http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf
colnames(data) <- c("Sample", "Gene",  "t.length", "lo.name", "q.length", "e.value", "score1", "bias1", "#", "of", "c.evalue", "i.evalue", "score2", "bias2", "from.hmm", "to.hmm", "from.ali", "to.ali", "from.env", "to.env", "acc")
















#### ### ### ### ### ### ### ### ### ### ### ### ##
#taxonomic classification-------------------------------------------------------------------------------------------------------------------------------------
#### ### ### ### ### ### ### ### ### ### ### ### ##

####given that my isolates are all there, I take the chance to reconstruct the total taxonomy by looking at the species name. 
###First_step: get clean names#
names_un = unique(data$Sample)
names = substr(names_un,1,nchar(names_un)-6) #remove anvi-o
names = gsub("_[0-9]*$","",names) #remove ID number Number 
names = gsub("_"," ",names) #remove underscore
names = gsub("sp$","sp.",names) #remove underscore
#given that I do not need the full name, I will only use until the genus and add sp.(to avoid 0 matches)
names_genuses = sub(" .*", "", names)
#names_genuses = paste0(names_genuses, " sp.")
names_genuses[names_genuses == "Proteus"] = "Proteus mirabilis" # I keep this like this because there are also samalandreae colled proteus
names_genuses[names_genuses == "Bacillus"] = "Bacillus alkalisoli"  # Same for baciullus (a walking stick is called like that, so I pick up a random bacillus species)
###look for taxonomy
taxonomy = data.frame()
for (i in 1:length(names_genuses))
{
  ids <- name2taxid(names_genuses[i], out_type="summary")
  classification = taxizedb::classification(ids)
  classification = as.data.frame(classification[ids$name])
  classification = classification[!duplicated(classification[,2]),]
  phylum = classification[classification[,2] == "phylum", 1]
  class = classification[classification[,2] == "class", 1]
  order = classification[classification[,2] == "order", 1]
  family = classification[classification[,2] == "family", 1]
  genus = classification[classification[,2] == "genus", 1]
  taxonomy1 = as.vector(c(phylum,class,order,family, genus))
  taxonomy = rbind(taxonomy, taxonomy1)
}
###cbind
taxonomy = cbind(taxonomy, names)
colnames(taxonomy) = c("phylum", "class","order" ,"family" , "genus", "species")
###Add their original accession  
rownames(taxonomy) = names_un
###Phylum overview of my samples
tax_aggregate = as.data.frame(table(taxonomy$phylum))
tax_aggregate$total = "total"
ggplot(tax_aggregate, aes(fill=Var1, y=Freq, x = total)) + 
  geom_bar(position="stack", stat = "identity") + theme_bw()  + ylab("Count") +  scale_fill_brewer(palette = "Dark2", name = "Phylum")+
  xlab(NULL) 

#save taxonomy file 
#saveRDS(taxonomy,"taxonomy.RDS")




##Compare species composition found here with murine community found in Schiro et al paper------------------------------------------
bac.tax.table<-read.table('/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/preliminary study/data/Mice/All_Seqs_silva.txt', sep='\t', header=T, row.names=1)
species = bac.tax.table[!is.na(bac.tax.table$Species),]
####Be careful, blanks were not removed
species = paste(species$Genus, species$Species) #List of species found in the amplicon sequencing data
matches = names[names %in% species]


###Read the data in (genomes data)
setwd("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/Second_part_genomes/January_results/Genomes")
filesInfolder <- list.files()
File_names = substr(filesInfolder,1,nchar(filesInfolder)-4)
#they are the same so I can load all them and append them in a list
list_files = list()
for (i in 1:length(File_names))#load all the genomes
{
  x <- readDNAStringSet(filesInfolder[i], format='fasta')
  assign(File_names[i], x)
  list_files[[i]] = get(File_names[i])
}
rm(list = File_names)#Remove files not in a list
names(list_files) = File_names #Rename the list objects
## read the sequences in from the amplicon sequencing project
seq = read.csv("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/preliminary study/data/sequences.csv",  row.names = 1 )
seq_tot = cbind(seq, bac.tax.table[match( seq$sequences,rownames(bac.tax.table)),])

#Results_16s = data.frame()
#####takes a while to run (10 minutes more or less)
#for (i in 1:length(list_files))
#  {
#    for (y  in 1:nrow(seq_tot))
#      {
#        results = vmatchPattern(seq_tot[y,2],list_files[[i]])
#        results = as.data.frame(results)
#        if (nrow(results) >0)
#          {
#              results$organism = File_names[i]
#              results[7:13] = seq_tot[y,3:9]
#              results$ASV = seq_tot[y,1]
#              Results_16s  = rbind(Results_16s, results)
#          }
#      }
#  }

## Save file ##
#write.table(Results_16s, file = "/Users/gabri/OneDrive - University of Arizona/Arsenic Project/Second_part_genomes/results/16_match/forward_match.txt")
Results_16s = read.table("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/Second_part_genomes/January_results/16_match/forward_match.txt")

#####Now reverse complement#
seq_tot$rf = NA
for ( i in 1:nrow(seq_tot))
{
  seq_tot$rf[i] = revComp(seq_tot$sequences[i])  
}


#Results_16s_rf = data.frame()
#####takes a while to run (10 minuntes more or less)
#for (i in 1:length(list_files))
#{
#  for (y  in 1:nrow(seq_tot))
#  {
#    results = vmatchPattern(seq_tot[y,10],list_files[[i]])
#    results = as.data.frame(results)
#    if (nrow(results)>0)
#    {
#      results$organism = File_names[i]
#      results[7:13] = seq_tot[y,3:9]
#      results$ASV = seq_tot[y,1]
#      Results_16s_rf  = rbind(Results_16s_rf, results)
#    }
#  }
# }
#write.table(Results_16s_rf, file = "/Users/gabri/OneDrive - University of Arizona/Arsenic Project/Second_part_genomes/results/16_match/reverse_match.txt")
Results_16s_rf = read.table("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/Second_part_genomes/January_results/16_match/reverse_match.txt")
####Match the tables###
results = rbind(Results_16s,Results_16s_rf)
results = results[!duplicated(results$organism),]
rm(x)




















####################                     ###
#HMM scan ---------------------------------------------------------------------- 
##Filter the results------------------------------------------------------------
####################                     ###

#Calculate the length of the alignment
data <- data %>%
  mutate(length = to.ali - from.ali) %>%
  mutate(perc.ali = length / q.length)

#plot data quality distribution
(quality <- ggplot(data, aes(x = perc.ali, y = score1)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~Gene, scales = "free_y") +
    ylab("Score") +
    xlab("Percent alignment") +
    theme_bw(base_size = 10))
quality
#remove rows that do not have at least 80% 
#of hmm length (std)
data.90 <- data[which(data$perc.ali > 0.90 & data$score1 > 100),]
#examine if any HMM hits apply to two genes
#of the duplicates, all are arrA/aioA mixed hits
#we will accept the one with a higher score
#arrange data by score
data.90 <- data.90[order(data.90$lo.name, abs(data.90$score1), decreasing = TRUE), ] 
#remove duplicates that have the lower score
data.90 <- data.90[!duplicated(data.90$lo.name),]
#dissim must have a score > 1000 so it doesnt
#pick up thiosulfate reductases ( do not use this method to subset because if integer(0)
#) is produces it cancels everything
data.90 <- subset(data.90, !(data.90$Gene == "arxA" & data.90$score1 < 1000))
data.90 <- subset(data.90, !(data.90$Gene == "aioA" & data.90$score1 < 1000)) 
#plot data quality distribution
ggplot(data.90, aes(x = perc.ali, y = score1)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~Gene, scales = "free_y")
#save table to output
#write.table(data.90, "/Users/gabri/OneDrive - University of Arizona/Arsenic Project/Second_part_genomes/results/results_hmmsearch/results_hmmsearch.txt", sep="\t") 



























####################                     ###
##HMM scan Statistics         ####
####################                     ###
length(unique(data.90$Sample)) # 34 / 73
#Ok now let's do a barplot with phylum
taxonomy_ars = taxonomy[rownames(taxonomy) %in% unique(data.90$Sample),]
tax_aggregate_ars = as.data.frame(table(taxonomy_ars$phylum))
tax_aggregate_ars$total = "Ars genes"
graph = rbind(tax_aggregate, tax_aggregate_ars)
ggplot(graph, aes(fill=Var1, y=Freq, x=total)) + 
  geom_bar(position="stack", stat="identity") + theme_bw()  + ylab("Count") +
  xlab(NULL) 
sort(table(data.90$Sample))
genes_count = as.data.frame(sort(table(data.90$Gene), decreasing =TRUE))
genes_count
#arsC      ACR3 arsC_thio      arsB      arsD arsC_glut      arsM      arsA      arrA 
#44        15        13         8         6         5         5         2         1 
p<-ggplot(data=genes_count, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity") + 
  theme_classic() +
  ylab("Count") +
  xlab("Gene") 
p
###create a table with counts
data_for_table = data.90[,1:2]
casted = data.frame(t(reshape2::dcast(data_for_table,  Gene ~ Sample)))
colnames(casted) = casted[1,]
casted = casted[-1,]
rownames(casted) =substr(rownames(casted),1,nchar(rownames(casted))-6)  
###export the table
#write.table(casted, "/Users/gabri/OneDrive - University of Arizona/Arsenic Project/Second_part_genomes/results/results_hmmsearch/casted_hmmsearch.txt", sep="\t") 








































































#KO orthologs ----------------------------------------------------------------
agg_file = read.csv("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/Second_part_genomes/Bio_informatics/KOfam/ko2level_Jan2021_right.csv")
arsenic = agg_file[grep("arse",agg_file$KO_description,fixed=TRUE),]
#Arsenic = grep("Arse",agg_file$KO_description,fixed=TRUE) #does not find anything
arsenic = arsenic[!duplicated(arsenic$KO),]
setwd(dir = "/Users/gabri/OneDrive - University of Arizona/Arsenic Project/Second_part_genomes/Bio_informatics/KOfam/KOresults")
names=list.files(pattern="*.txt")
data1 <- do.call(rbind, lapply(names, function(X) { 
  data.frame(id = basename(X), read.table(X, sep = '\t', header = FALSE, fill =TRUE, col.names=c("locus","KO_ortho")))})) #due to name differences
data1 = data1[which(data1$KO_ortho %in% arsenic$KO), ]





length(unique(data1$id)) # 59/73
#Ok now let's do a barplot with phylum
#change ID name 
rownames(taxonomy) = substr(rownames(taxonomy),1,nchar(rownames(taxonomy))-6) #remove anvi-o
data1$id = substr(data1$id,1,nchar(data1$id)-4)
taxonomy_ars = taxonomy[rownames(taxonomy) %in% unique(data1$id),]
tax_aggregate_ars = as.data.frame(table(taxonomy_ars$phylum))
tax_aggregate_ars$total = "Ars genes"
graph = rbind(tax_aggregate, tax_aggregate_ars)
ggplot(graph, aes(fill=Var1, y=Freq, x=total)) + 
  geom_bar(position="stack", stat="identity") + theme_bw()  + ylab("Count") +
  xlab(NULL) 
sort(table(data.90$Sample))
genes_count = as.data.frame(sort(table(data1$KO_ortho), decreasing =TRUE))
genes_count$name = arsenic[match(genes_count$Var1, arsenic$KO), 5]
genes_count$name_short = gsub(",.*$", "", genes_count$name)
genes_count$name_short = gsub(";.*$", "", genes_count$name_short)
genes_count$name_short[genes_count$name == "arsC; arsenate reductase (glutaredoxin) [EC:1.20.4.1]"] = "arsC_g"
genes_count$name_short[genes_count$name == "arsC; arsenate reductase (thioredoxin) [EC:1.20.4.4]"] = "arsC_t"
genes_count$name_short = factor(genes_count$name_short, levels = genes_count$name_short)
genes_count1 = genes_count[,c(2,4)]
genes_count1
#Freq name_short
#1   79       arsR
#2   49     arsC_g
#3   20     arsC_t
#4   18       ACR3
#5   12       arsA
#6   11       arsB
#7    5      AS3MT
p<-ggplot(data=genes_count1, aes(x=name_short, y=Freq)) +
  geom_bar(stat="identity") + 
  theme_classic() +
  ylab("Count") +
  xlab("Gene") 
p



###create a table with counts
data_for_table_KO = data1[,c(1,3)]
casted = data.frame(t(reshape2::dcast(data_for_table_KO,  KO_ortho ~ id)))
colnames(casted) = casted[1,]
casted = casted[-1,]
colnames(casted) = genes_count[match(genes_count$Var1, colnames(casted)),4 ]
###export the table
#write.table(casted, "/Users/gabri/OneDrive - University of Arizona/Arsenic Project/Second_part_genomes/results/results_hmmsearch/casted_hmmsearch.txt", sep="\t") 





















#Arsenic from Prokka annotation-----------------------------------------------------------------------
ars_ann = data.frame()
setwd(dir = "/Users/gabri/OneDrive - University of Arizona/Arsenic Project/Second_part_genomes/January_results/gff_files")
names=list.files(pattern="*.gff")
for (i in 1:length(names))
{
  genbank = as.data.frame(import.gff(names[i]))
  gen_ars = genbank[grep("arse",genbank$product,fixed=TRUE),]
  Gen_ars = genbank[grep("Arse",genbank$product,fixed=TRUE),]
  if (nrow(Gen_ars)>0) gen_ars = rbind(gen_ars, Gen_ars)
  if (nrow(gen_ars)>0) 
    {
      for (a in 1:nrow(gen_ars))gen_ars$inference1[a] = unlist(gen_ars$inference[[a]])[2]
      gen_ars$isolate = basename(names[i])
      ars_ann = rbind(ars_ann, gen_ars)
  }
}





###Fill the names of NAs
ars_annNA = ars_ann[is.na(ars_ann$gene),]
###names manually searched from the file "/Users/gabri/OneDrive - University of Arizona/Arsenic Project/Second_part_genomes/Bio_informatics/hmm_PGAP-2
data_frame_names = data.frame(inference = unique(ars_annNA$inference1), names = c("arsenic_eff", "acr3", "arsM", "arsN2", "perox_w_seleSAM", "tranport_ArsG", "glyco_like_mftF","arsC", "arsS", "glyco_like_cofC")) 
ars_ann$Name = data_frame_names[match(ars_ann$inference1, data_frame_names$inference),2]
###
ars_ann$Name <- ifelse(is.na(ars_ann$Name), ars_ann$gene, ars_ann$Name)
length(unique(ars_ann$isolate)) #55
ars_ann$annotation = "Prokka"


#Put results togheter---------------------------------------------------------------
####FUrther annotation from KEGGS
####Loci Keggs
KEGGS_loc = data.frame() 
for (i in 1:length(names))
{
  genbank = as.data.frame(import.gff(names[i]))
  gen_ars = genbank[genbank$ID %in% data1$locus ,]
  if (nrow(gen_ars)>0) 
  {
    for (a in 1:nrow(gen_ars))gen_ars$inference1[a] = unlist(gen_ars$inference[[a]])[2]
    gen_ars$isolate = basename(names[i])
    KEGGS_loc = rbind(KEGGS_loc, gen_ars)
  }
}
KEGGS_loc$annotation = "KOfam"
ars_ann = ars_ann[!( ars_ann$ID %in% KEGGS_loc$ID),]


####Loci DUVININ
DUV_loc = data.frame() 
for (i in 1:length(names))
{
  genbank = as.data.frame(import.gff(names[i]))
  gen_ars = genbank[genbank$ID %in% data.90$lo.name ,]
  if (nrow(gen_ars)>0) 
  {
    for (a in 1:nrow(gen_ars))gen_ars$inference1[a] = unlist(gen_ars$inference[[a]])[2]
    gen_ars$isolate = basename(names[i])
    DUV_loc = rbind(DUV_loc, gen_ars)
  }
}
DUV_loc$annotation = "Dunivin"
ars_ann = ars_ann[!( ars_ann$ID %in% DUV_loc$ID),]
DUV_loc = DUV_loc[!(DUV_loc$ID %in% KEGGS_loc$ID),]
final_table = rbind(DUV_loc, KEGGS_loc)
final_table = rbind(final_table, ars_ann)





###Insert_Consensus (high trustable source KOfam, second Dunivin, last Prokka)
final_table$KOfam = data1[match(final_table$locus_tag, data1$locus),3]
final_table$KOfam_gene = genes_count[match(final_table$KOfam, genes_count$Var1),4]
final_table$Dunivin = data.90[match(final_table$locus_tag, data.90$lo.name),2]

for (i in 1:nrow(final_table))
{
  if (final_table$annotation[i] == "Dunivin") final_table$final_gene[i] = final_table$Dunivin[i]
  else if  (final_table$annotation[i] == "KOfam") final_table$final_gene[i] = toString(final_table$KOfam_gene[i])
  else final_table$final_gene[i]  = final_table$Name[i]
}
sort(table(final_table$final_gene))

final_table$final_gene = str_remove(final_table$final_gene, "_[0-9]")
sort(table(final_table$final_gene))
final_table[final_table$final_gene == 'acr3',25] = "ACR3"
sort(table(final_table$final_gene))
final_table[final_table$final_gene == 'arsC_thio',25] = "arsC_t"
sort(table(final_table$final_gene))


to_aggregate = final_table[,c(25,20)]
casted = data.frame(t(reshape2::dcast(to_aggregate,  final_gene ~ isolate)))
colnames(casted) = casted[1,]
casted=  casted[-1,]
rownames(casted) = substr(rownames(casted),1,nchar(rownames(casted))-4) 
#write.table(casted, "/Users/gabri/OneDrive - University of Arizona/Arsenic Project/Second_part_genomes/results/final_table.txt", sep="\t") 
ncol(casted) # too many genes for a nice plot, lets reduce it to ten 9 + others



###Make a list for output
output_table = data.frame(isolate = final_table$isolate, contig = final_table$seqnames,locus = final_table$locus_tag ,start = final_table$start, end = final_table$end, strand = final_table$strand, 
                        gene = final_table$final_gene, annotation = final_table$annotation, Kofam_name = final_table$KOfam_gene, KEGG_ID = final_table$KOfam, Dunuvin_gene = final_table$Dunivin, 
                        Prokka_gene = final_table$Name, Prokka_inference = final_table$inference1, Prokka_description = final_table$product)

#write.table(output_table, "/Users/gabri/OneDrive - University of Arizona/Arsenic Project/Second_part_genomes/results/final_table_details.txt", sep="\t") 








#Make a nice graph----------------------------------------------------------
###first I need to remove some genes (too many, I will only use 10)
`%nin%` = Negate(`%in%`)
to_aggregate$final_gene = ifelse(to_aggregate$final_gene %nin% names(sort(table(final_table$final_gene),decreasing = TRUE))[1:9],  "others", to_aggregate$final_gene)
casted = data.frame(t(reshape2::dcast(to_aggregate,  final_gene ~ isolate)))
colnames(casted) = casted[1,]
casted=  casted[-1,]
rownames(casted) = substr(rownames(casted),1,nchar(rownames(casted))-4) 


ncol(casted) # too 
#read the tree
tree <- read.tree("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/Second_part_genomes/January_results/RAxML_bestTree.NT_seqs_refined.tre")
### Short the names (better visualization)
name_tips = tree$tip.label
#name_tips = gsub("_[0-9]*$","",name_tips) #remove ID number Number 
name_tips = gsub("_"," ",name_tips) #remove underscore
#name_tips = gsub("sp$","sp.",name_tips) #remove underscore
name_tips[ name_tips == "Bifidobacterium longum subsp animalis 26074"] = "Bifidobacterium longum 26074"
name_tips[ name_tips == "Bacteroidetes xylanisolvens 100502"] = "Bacteroides xylanisolvens 100015" #There was a mistake with the name
tree$tip.label  = name_tips
p <- ggtree(tree, branch.length='none', layout='circular')  + 
  geom_tiplab(align= T, linetype=NA, 
              size=2.5, offset=0.5) 
p
##add genes
m <- matrix(0, ncol =  10, nrow = (length(name_tips) - nrow(casted)))
colnames(m) = colnames(casted)
rownames(casted) = gsub("_"," ",rownames(casted)) #remove underscore
rownames(casted)[ rownames(casted) == "Bifidobacterium longum subsp animalis 26074"] = "Bifidobacterium longum 26074"
casted1 = rbind(casted, m)
rownames(casted1)[(nrow(casted)+1):nrow(casted1)] = subset(name_tips, !(name_tips %in% rownames(casted)))
rn <- rownames(casted1)
colors = c("white" ,"green" ,"yellow", "orange", "red","dark red" , "purple")
names(colors) <- c(0,1,2,3,4,5,7)
gheatmap(p, casted1, offset = 30, color="grey", 
         colnames_position="top", 
         colnames_angle=90, colnames_offset_y = -0.1, 
         hjust=0, font.size=2.2) +
  scale_fill_manual(values=colors, breaks=c(0,1,2,3,4,5,7), name = "Counts")  



##Final stats
length(unique(final_table$isolate)) # 64 / 73
#Ok now let's do a barplot with phylum
final_table$isolate = substr(final_table$isolate,1,nchar(final_table$isolate)-4) #remove anvi-o


taxonomy_ars = taxonomy[rownames(taxonomy) %in% unique(final_table$isolate),]
tax_aggregate_ars = as.data.frame(table(taxonomy_ars$phylum))
tax_aggregate_ars$total = "Ars genes"
graph = rbind(tax_aggregate, tax_aggregate_ars)
ggplot(graph, aes(fill=Var1, y=Freq, x=total)) + 
  geom_bar(position="stack", stat="identity") + theme_bw()  + ylab("Count") +
  xlab(NULL) 

genes_count = as.data.frame(sort(table(final_table$final_gene), decreasing =TRUE))
genes_count
#arsC      ACR3 arsC_thio      arsB      arsD arsC_glut      arsM      arsA      arrA 
#44        15        13         8         6         5         5         2         1 
p<-ggplot(data=genes_count, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity") + 
  theme_classic() +
  ylab("Count") +
  xlab("Gene") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
###create a table with counts
data_for_table = data.90[,1:2]
casted = data.frame(t(reshape2::dcast(data_for_table,  Gene ~ Sample)))
colnames(casted) = casted[1,]
casted = casted[-1,]
rownames(casted) =substr(rownames(casted),1,nchar(rownames(casted))-6)  
















###GGgene example###
example = gggenes::example_genes
ggplot(example_genes, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3")
df = example[-c(1:nrow(example)),]
output_table <- output_table[order(output_table$locus),] 
row.names(output_table) <- NULL
df[1:5,]= NA
df[1:5,1] = output_table$isolate[21:25]
df[1:5,2] = output_table$gene[21:25]
df[1:5,3] = output_table$start[21:25]
df[1:5,4] = output_table$end[21:25]
df[1:5,5] = output_table$strand[21:25]

df[6:12,1] = output_table$isolate[148:154]
df[6:12,2] = output_table$gene[148:154]
df[6:12,3] = output_table$start[148:154]
df[6:12,4] = output_table$end[148:154]
df[6:12,5] = output_table$strand[148:154]

df[13:17,1] = output_table$isolate[217:221]
df[13:17,2] = output_table$gene[217:221]
df[13:17,3] = output_table$start[217:221]
df[13:17,4] = output_table$end[217:221]
df[13:17,5] = output_table$strand[217:221]

df[18:22,1] = output_table$isolate[303:307]
df[18:22,2] = output_table$gene[303:307]
df[18:22,3] = output_table$start[303:307]
df[18:22,4] = output_table$end[303:307]
df[18:22,5] = output_table$strand[303:307]

df[,1] =  substr(df[,1],1,nchar(df[,1])-4) 


ggplot(df, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3")





















#construct phylogenetic tree
#file containing all the sequences from the loci found, generated with seqtk #see list at #"/Volumes/Gabri_Bio/Arsenic_project/Genomes/newfile1.fasta"
sequences = seqinr::read.fasta("/Volumes/Gabri_Bio/Arsenic_project/Genomes/newfile1.fasta")
sequences <- sequences[!duplicated(names(sequences))] #remove duplicates
###first_do arsC
subset = output_table[output_table$gene=="arsenic_eff",]
sequences_arsC = sequences[names(sequences) %in% subset$locus]
FUN = function(x) paste(getSequence(x), collapse = "") # Function to make it a DNA_String_set
seq = as(vapply(sequences_arsC, FUN, character(1)), "DNAStringSet") # Function to make it a DNA_String_set
alignment <- AlignSeqs(DNAStringSet(seq), anchor=NA,verbose=FALSE) #align sequences
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA") #translate it into a phangorn file
dm <- dist.ml(phangAlign) #calculate distance matrix
treeNJ <- NJ(dm) # NJ tree estimation
fit = pml(treeNJ, data=phangAlign) #computes the likelihood of a tree
fitJC <- optim.pml(fit, model="JC", rearrangement = "stocastic")
bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
fit = plotBS(fitJC$tree, bs, p = 50, type="p")
more = output_table$isolate[match(fit$tip.label,output_table$locus,)]
substrRight <- function(x, n){#Define a function to isolate last 5 characters
  substr(x, nchar(x)-n+1, nchar(x)) 
}
locuses = substrRight(fit$tip.label, 5) #Isolate last 5 characters
locuses <- paste0("(", locuses, ")") #Add parentheses
more = substr(more,1,nchar(more)-4) #remove anvi-o
more = gsub("_[0-9]*$","",more) #remove ID number Number 
more = gsub("_"," ",more) #remove underscore
more[more == "Bifidobacterium longum subsp animalis"] = "Bifidobacterium longum"
complete = paste(more, locuses)
fit$tip.label = complete
f <- ggtree(fit)  + theme(plot.margin = margin(10, 100, 10, 10))+ geom_tiplab() +geom_nodelab(size = 2.5, nudge_x = 0.01) 
f
ggsave("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/Second_part_genomes/images/plot.pdf", device = "pdf")







