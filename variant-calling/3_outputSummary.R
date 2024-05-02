#### Prior to this script, run the Flumina pipeline to obtain variant calls in VCF files
#### Next, run 2_convertVCFtoTable.R to get an amino acid table, all_sample_amino_acids.txt


#Takes seconds to run on 500 samples

#Main working directory
work.directory = "/path/to/working/directory"

# output.directory used from 1_convertVCFtoTable.R
output.directory = "variant_analysis"

# curated csv database with the columns "Gene", "Amino_Acid", and "Type" (category of site) of interest
# Other columns can be added and joined with the variants
aa.table.path = "curated_database.csv"

#Grouping category column name joined in step 2 with the variant sample data
#Set to NULL if there are no groupings to use
group.names = "discrete_host"


#############################################
#### Should not need to modify below here
#############################################

#sets working directory 
setwd(work.directory)

#read in previous database
sample.data = read.table(paste0(output.directory, "/all_sample_amino_acids.txt"), sep = "\t", header = T)

#Reads in curated database
best.aa = read.csv(aa.table.path, header = TRUE, sep = ",")
best.aa[is.na(best.aa) == TRUE] = "NA"

#Obtains gene names
gene.names = unique(sample.data$locus)

#Combines all the data into 1 data frame
save.sample = c()
for (i in 1:length(gene.names)){
  
  #Subsets data
  temp.sample = sample.data[sample.data$locus %in% gene.names[i],]
  temp.fun = best.aa[best.aa$Gene %in% gene.names[i],]
  
  #Merges with curated database
  merge.sample = merge(temp.sample, temp.fun, by.x = "aa_position", by.y = "Amino_Acid")  
  save.sample = rbind(save.sample, merge.sample)

} #end i 

#Save large tab delimited table of all the amino acids
write.table(save.sample, paste0(output.directory, "/curated_amino_acids.txt"),
            row.names = F, quote = F, sep = "\t")


################################
### Make summary table
################################

#Only amino acid changing sites
red.samples = save.sample[save.sample$aa_changing == "YES",]

#Removes duplicate stuff so they aren't counted twice
red.samples$sample = gsub("-r_", "_", red.samples$sample)
red.samples$sample = gsub("-v_", "_", red.samples$sample)

#minimum allele frequency 
red.samples = red.samples[red.samples$allele_frequency >= 0.005,]

#Obtains different animal groups, if desired
if (is.null(group.names) != TRUE){
  group.values = unique(red.samples[[group.names]])
} else{ group.values = "all" }

all.data = c()
for (i in 1:length(group.values)){

  #Obtains different animal groups
  if (is.null(group.names) != TRUE){
    group.data = red.samples[red.samples[[group.names]] %in% group.values[i],]
  } else{ group.data = red.samples }  
  
  #Obtains gene names
  gene.names = unique(group.data$locus)
  
  #Creates table
  sample.table = c()
  for (j in 1:length(gene.names)){
    
    #Subsets data
    gene.data = group.data[grep(gene.names[j], group.data$locus),]
    
    #Gets unique positions
    aa.pos = unique(gene.data$aa_position)
    
    #empty table
    temp.table = data.frame(animal = as.character(),
                            gene = as.character(),
                            aa_position = as.numeric(),
                            mutation = as.character(),
                            funct = as.character(),
                            no_animal = as.numeric(),
                            ave_alle_freq = as.numeric(),
                            consensus = as.numeric(),
                            low_freq = as.numeric())
    
    #Loops through positions
    for (k in 1:length(aa.pos)){
      
      #empty data for this section
      vector.table = data.frame(animal = as.character(),
                                gene = as.character(),
                                aa_position = as.numeric(),
                                mutation = as.character(),
                                funct = as.character(),
                                no_animal = as.numeric(),
                                ave_alle_freq = as.numeric(),
                                consensus = as.numeric(),
                                low_freq = as.numeric())
      
      #Subsets to amino acid data
      aa.data = gene.data[gene.data$aa_position == aa.pos[k],]
      
      ##### Groupings categories
      vector.table[1,1] = group.values[i]
      vector.table[1,2] = gene.names[j]
      vector.table[1,3] = aa.pos[k]
      
      #Formats mutation names
      ref.chars = unique(aa.data$reference_aa)
      aa.chars = unique(aa.data$alternative_aa)
      new.chars = paste0(ref.chars, aa.pos[k], aa.chars)
      vector.table[1,4] = paste0(new.chars, collapse = ", ")
      
      #saves AA type
      vector.table[1,5] = unique(aa.data$Type)
      
      #Counts samples
      vector.table[1,6] = length(unique(aa.data$sample))
      
      #mean allele frequency
      vector.table[1,7] = round(mean(aa.data$allele_frequency), 3)
      
      #gets those in consensus
      freq.aa = aa.data[aa.data$allele_frequency > 0.50,]
      vector.table[1,8] = length(unique(freq.aa$sample))
      
      #gets those in low frequency
      freq.aa = aa.data[aa.data$allele_frequency < 0.50,]
      vector.table[1,9] = length(unique(freq.aa$sample))
      
      #combines togther
      temp.table = rbind(temp.table, vector.table)
      
    }# end k loop
    
    #saves sample data
    sample.table = rbind(sample.table, temp.table)
    
  }#end j loop
    
  #saves all the data
  all.data = rbind(all.data, sample.table)
  
}# i loop

#Saves final data
write.table(all.data, paste0(output.directory, "/animal_summary_curated_sites.txt"),
            row.names = F, quote = F, sep = "\t")












