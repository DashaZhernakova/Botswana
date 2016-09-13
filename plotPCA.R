library("SNPRelate")
genofile <- snpgdsOpen("/Volumes/rset1/Botswana_chip/cleaned_all_final_no_outliers.gds")
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)
pc.percent <- pca$varprop*100
pdf("/Volumes/rset1/Botswana_chip/cleaned_all_final_no_outliers.PCA.pdf")
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], labels=lbls)
dev.off()

ibd <- snpgdsIBDMoM(genofile, snp.id=snpset.id, maf=0.05, missing.rate=0.05, num.thread=2)
ibd.coeff <- snpgdsIBDSelection(ibd)
pdf("/Volumes/rset1/Botswana_chip/cleaned_all_final_no_outliers.IBD.pdf")
plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1), xlab="k0", ylab="k1", main = "IBD plot")
lines(c(0,1), c(1,0), col="red", lty=2)

dev.off()
