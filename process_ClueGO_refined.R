



library(xlsx)


process_ClueGO = function(cluego_data)
{

data = read.xlsx(cluego_data, 1)

#head(data)


########################################
#####	upregulated or cluster 1   #####
########################################


data_upreg = data[which(data$Cluster == 'Specific for Cluster #1'),] ##if you don't have two clusters, substitute this line with "data_upreg = data"

#head(data_upreg)

data_upreg_freq = as.data.frame(table(data_upreg$GOGroups))

data_upreg_freq$perc = (data_upreg_freq$Freq/sum(data_upreg_freq$Freq))*100

data_upreg_freq = data_upreg_freq[order(data_upreg_freq$Freq, decreasing = TRUE),]

colnames(data_upreg_freq)[1] = 'GOGroups'

data_upreg_freq_merge = merge(data_upreg_freq, data_upreg, by = 'GOGroups', sort = FALSE)

write.table(data_upreg_freq_merge, file = paste(unlist(strsplit(as.character(cluego_data), split = '.xls')), 'upreg_cluego_freqs.xls', sep = '_'), sep = '\t', row.names = FALSE)


sig_terms = data.frame()

for (i in 1:length(unique(data_upreg_freq_merge$GOGroups)))
{
	
	sorted_subtable = data_upreg_freq_merge[which(data_upreg_freq_merge$GOGroups == unique(data_upreg_freq_merge$GOGroups)[1]), ]
	sorted_subtable = sorted_subtable[order(sorted_subtable$X..Associated.Genes, decreasing = TRUE),]
	most_sig_term = sorted_subtable[1,]
	sig_terms = rbind(sig_terms, most_sig_term) 
}

write.table(sig_terms, file = paste(unlist(strsplit(as.character(cluego_data), split = '.xls')), 'upreg_compact_cluego_freqs.xls', sep = '_'), sep = '\t', row.names = FALSE)



########################################
#####	downregulated or cluster 2 #####
########################################

###if you dont'have 2 clusters, remove from here... 

data_downreg = data[which(data$Cluster == 'Specific for Cluster #2'),]

#head(data_downreg)

data_downreg_freq = as.data.frame(table(data_downreg$GOGroups))

data_downreg_freq$perc = (data_downreg_freq$Freq/sum(data_downreg_freq$Freq))*100

data_downreg_freq = data_downreg_freq[order(data_downreg_freq$Freq, decreasing = TRUE),]

colnames(data_downreg_freq)[1] = 'GOGroups'

data_downreg_freq_merge = merge(data_downreg_freq, data_downreg, by = 'GOGroups', sort = FALSE)

write.table(data_downreg_freq_merge, file = paste(unlist(strsplit(as.character(cluego_data), split = '.xls')), 'downreg_cluego_freqs.xls', sep = '_'), sep = '\t', row.names = FALSE)


sig_terms = data.frame()


for (i in 1:length(unique(data_downreg_freq_merge$GOGroups)))
{
	
	sorted_subtable = data_downreg_freq_merge[which(data_downreg_freq_merge$GOGroups == unique(data_downreg_freq_merge$GOGroups)[i]), ]
	sorted_subtable = sorted_subtable[order(sorted_subtable$X..Associated.Genes, decreasing = TRUE),]
	most_sig_term = sorted_subtable[1,]
	sig_terms = rbind(sig_terms, most_sig_term) 
}

write.table(sig_terms, file = paste(unlist(strsplit(as.character(cluego_data), split = '.xls')), 'downreg_compact_cluego_freqs.xls', sep = '_'), sep = '\t', row.names = FALSE)


###...to here

}



#process_ClueGO('BRD2_KO_old_BP_by_cluster.xls')

