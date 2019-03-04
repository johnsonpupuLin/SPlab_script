### input ######
colnameList <- c("Motif", "Consensus", "pVal", "ln_pVal", "qVal", "Num.motifs", "Percent.motifs", "Num.bg.motifs", "Percent.bg.motifs")
fileList <- c("HOMER.piRClust_upregulated", 
              "HOMER.piRClust_downregulated",
              "HOMER.piRClust_strict_upregulated",
              "HOMER.piRClust_strict_downregulated")


nameList <- c("1.5_up","1.5_down","2_up","2_down")
clusters <- list()

for(i in 1:length(nameList)) {
	clusters[[i]] <- matrix(scan(paste(fileList[i],"knownResults.txt",sep="/"), "character", sep="\t"),ncol=9,byrow=TRUE)
	colnames(clusters[[i]]) <- colnameList
	clusters[[i]] <- data.frame(clusters[[i]][-1,], stringsAsFactors=FALSE)
	for(idx in 3:9) {
		clusters[[i]][,idx] <- as.numeric(sub("%","",clusters[[i]][,idx]))
	}
	clusters[[i]]$pVal <- exp(clusters[[i]]$ln_pVal)
	clusters[[i]]$qVal <- p.adjust(clusters[[i]]$pVal, method="fdr")
}
names(clusters) <- nameList

####### top 30% lowest qValue + top 30% highest percent motifs #########
qVal_TH <- 0.30
pMotif_TH <- 0.30
sigClusters <- list()
for(i in 1:length(clusters)) {
	lo_qVal <- sort(clusters[[i]]$qVal)[floor( nrow(clusters[[i]]) * qVal_TH )]
	hi_pMotif <- sort(clusters[[i]]$Percent.motifs, decreasing=TRUE)[floor( nrow(clusters[[i]]) * pMotif_TH )]
	sigClusters[[i]] <- clusters[[i]][which(clusters[[i]]$qVal <= lo_qVal & clusters[[i]]$Percent.motifs >= hi_pMotif),]
}
names(sigClusters) <- nameList

output <- c()
for (i in nameList) {
	output <- rbind(output, cbind(sigClusters[[i]][,c(1,2,5,7,9)],sig=i))
}
write.csv(output, "sigClust.motif.csv", row.names=FALSE)
write.table(output, "sigClust.motif.table", sep="\t", quote=FALSE, row.names=FALSE)

#write.csv(rbind(cbind(sigClusters$E11Gup.u9k[,c(1,2,5,7,9)],sig="E11G")
#		   ,cbind(sigClusters$E14Gup.u9k[,c(1,2,5,7,9)],sig="E14G"))
#	   ,"sigClust.u9k.staged.csv")

