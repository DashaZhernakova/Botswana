library(gdsfmt)
library(SNPRelate)
#
#Filtering
#
dat.f<-openfn.gds("/Volumes/rset1/Botswana_chip/cleaned_all_final.gds")
id<-read.gdsn(index.gdsn(dat.f,"sample.id"))
closefn.gds(dat.f)

to_remove <- c("NG70367","NG50797","TH9049702001506","NG51827","NG51838","BN10714","TH8977302001505","TH8988402001506")
to_keep <- id[which(! id %in% to_remove)]
snpgdsCreateGenoSet("/Volumes/rset1/Botswana_chip/cleaned_all_final.gds", "/Volumes/rset1/Botswana_chip/cleaned_all_final_no_outliers.gds", sample.id=to_keep)

dat.f1<-openfn.gds("/Volumes/rset1/Botswana_chip/default.gds")
id1<-read.gdsn(index.gdsn(dat.f1,"sample.id"))
closefn.gds(dat.f1)

to_keep1 <- id1[which(! id1 %in% to_remove)]
snpgdsCreateGenoSet("/Volumes/rset1/Botswana_chip/default.gds", "/Volumes/rset1/Botswana_chip/default_no_outliers.gds", sample.id=to_keep1)


#
# LD and IBD
#
dat.f<-openfn.gds("/Volumes/rset1/Botswana_chip/cleaned_all_final_no_outliers.gds")

set.seed(1000)
snpset <- snpgdsLDpruning(dat.f, ld.threshold=0.2)
pca <- snpgdsPCA(dat.f, snp.id=snpset.id, num.thread=2)
pdf("/Volumes/rset1/Botswana_chip/cleaned_all_final_no_outliers.PCA.pdf")
plot(pca$eigenvect[,1], pca$eigenvect[,2])
dev.off()

ibd <- snpgdsIBDMoM(dat.f, maf=0.05, missing.rate=0.05, num.thread=2)
ibd.coeff <- snpgdsIBDSelection(ibd)
head(ibd.coeff)
pdf("/Volumes/rset1/Botswana_chip/cleaned_all_final_no_outliers.IBD.pdf")
plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1), xlab="k0", ylab="k1", main="Botswana IBD")
lines(c(0,1), c(1,0), col="red", lty=2)
dev.off()


id<-read.gdsn(index.gdsn(dat.f,"sample.id"))
chr<-read.gdsn(index.gdsn(dat.f,"snp.chromosome"))
pos<-read.gdsn(index.gdsn(dat.f,"snp.position"))
snp<-read.gdsn(index.gdsn(dat.f,"snp.id"))
all<-read.gdsn(index.gdsn(dat.f,"snp.allele"))
gen<-read.gdsn(index.gdsn(dat.f,"genotype"))
closefn.gds(dat.f)
#
id.df<-data.frame(id=id,id1=id)
#
dat.f1<-openfn.gds("/Volumes/rset1/Botswana_chip/default_no_outliers.gds")
id1<-read.gdsn(index.gdsn(dat.f1,"sample.id"))
chr1<-read.gdsn(index.gdsn(dat.f1,"snp.chromosome"))
pos1<-read.gdsn(index.gdsn(dat.f1,"snp.position"))
snp1<-read.gdsn(index.gdsn(dat.f1,"snp.id"))
all1<-read.gdsn(index.gdsn(dat.f1,"snp.allele"))
gen1<-read.gdsn(index.gdsn(dat.f1,"genotype"))
closefn.gds(dat.f1)
#
ftr<-array(FALSE,dim=length(id))
l<-length(id)
for (i in 1:l) {if (sum(as.numeric(id1==id[i]))==1) ftr[i]<-TRUE}
#
id2<-id[ftr]
gen2<-gen[ftr,]
#
ftr1<-array(FALSE,dim=length(id1))
l1<-length(id1)
for (i in 1:l1) {if (sum(as.numeric(id2==id1[i]))==1) ftr1[i]<-TRUE}
# 
id3<-id1[ftr1]
gen3<-gen1[ftr1,]
# 
ord2<-order(id2)
id2<-id2[ord2]
gen2<-gen2[ord2,]
# 
ord3<-order(id3)
id3<-id3[ord3]
gen3<-gen3[ord3,]
dd<-length(id3)
# 
CHR<-NULL
SNP<-NULL
POS<-NULL
ALL<-NULL
SNP1<-NULL
POS1<-NULL
ALL1<-NULL
ZP<-NULL
ZPA<-NULL
ZD<-NULL
ZDD<-NULL
ZR<-NULL
for (ch in levels(as.factor(chr1)))
{
	ff.c<-chr==ch
	ff1.c<-chr1==ch
	gen.c<-gen2[,ff.c]
	gen1.c<-gen3[,ff1.c]
	pos.c<-pos[ff.c]
	pos1.c<-pos1[ff1.c]
	snp.c<-snp[ff.c]
	snp1.c<-snp1[ff1.c]
	all.c<-all[ff.c]
	all1.c<-all1[ff1.c]
 	#
	df.c<-data.frame(pos=pos.c)
	df1.c<-data.frame(pos=pos1.c)
	df2.c<-merge(df.c,df1.c,by="pos")
	ftr.c<-array(FALSE,dim=length(pos.c))
	ftr1.c<-array(FALSE,dim=length(pos1.c))
	pos2.c<-df2.c$pos	
	for (i in 1:length(pos.c))
	{
		if (sum(as.numeric(pos2.c==pos.c[i]))==1) ftr.c[i]<-TRUE
	}
	for (i in 1:length(pos1.c))
	{
	if (sum(as.numeric(pos2.c==pos1.c[i]))==1) ftr1.c[i]<-TRUE
	}
	gen2.c<-gen.c[,ftr.c]
	gen3.c<-gen1.c[,ftr1.c]
	snp2.c<-snp.c[ftr.c]
	snp3.c<-snp1.c[ftr1.c]
	pos2.c<-pos.c[ftr.c]
	pos3.c<-pos1.c[ftr1.c]	
	all2.c<-all.c[ftr.c]
	all3.c<-all1.c[ftr1.c]		
	chr.c<-array(ch,dim=length(pos2.c))
	dm<-dim(gen2.c)
	zp<-colSums(gen2.c!=gen3.c)/dd
	zpa<-colSums((gen2.c!=gen3.c)&(gen2.c!=3 & gen3.c!=3))/colSums(gen2.c!=3 & gen3.c!=3)
	#do not count errors of this sort: rec homo = hetero
	zd<-1-colSums(array(as.numeric((gen2.c==1|gen2.c==2) & (gen3.c==1|gen3.c==2)),dim=dm))/ colSums(array(as.numeric((gen2.c==1|gen2.c==2) | (gen3.c==1|gen3.c==2)),dim=dm))
	zd[(zd=="NaN"|is.na(zd))]<-0
	zdd<-1-colSums(array(as.numeric((gen2.c==0) & (gen3.c==0)),dim=dm))/ colSums(array(as.numeric((gen2.c==0|gen3.c==0)),dim=dm))
	zdd[(zdd=="NaN"|is.na(zdd))]<-0
	zr<-1-colSums(array(as.numeric((gen2.c==2) & (gen3.c==2)),dim=dm))/ colSums(array(as.numeric(gen2.c==2 | gen3.c==2),dim=dm))
	zr[(zr=="NaN"|is.na(zd))]<-0
	CHR<-c(CHR,chr.c)
	SNP<-c(SNP,snp2.c)
	POS<-c(POS,pos2.c)
	ALL<-c(ALL,all2.c)
	SNP1<-c(SNP1,snp3.c)
	POS1<-c(POS1,pos3.c)
	ALL1<-c(ALL1,all3.c)
	ZP<-c(ZP,zp)
	ZPA<-c(ZPA,zpa)
	ZD<-c(ZD,zd)
	ZDD<-c(ZDD,zdd)
	ZR<-c(ZR,zr)
}
d.out<-data.frame(chr=CHR,snp=SNP,pos=POS,all=ALL,snp1=SNP1,pos1=POS1,all1=ALL1,zp=ZP,zpa=ZPA,zd=ZD, zr=ZR,zdd=ZDD)
write.table(d.out,file="/Volumes/rset1/Botswana_chip//qc.csv",row.names=FALSE,sep=",")



# Allele control
d.in<-read.table("/home/smalov/BOTSWANA3/DATA/qc.csv",header=TRUE,sep=",")





ord1<-order(id1)
id1<-id1[ord1]
snp1<-snp1[ord1]
gen1<-gen1[ord1,]
#
snp.df<-data.frame(snp=snp,snp1=snp,ord=c(1:length(snp)))
snp.df1<-data.frame(snp=snp1,snp2=snp1,ord1=c(1:length(snp1)))
snp.dff<-merge(snp.df,snp.df1,by="snp",all=TRUE)
# 
#ftr2<-chr==1
dat.f2<-openfn.gds("/home/smalov/BOTSWANA3/DATA/with-rare-variants.gds")
snp11<-read.gdsn(index.gdsn(dat.f2,"snp.id"))
closefn.gds(dat.f2)
#
l2<-length(snp)
ftr2<-array(FALSE,dim=l2)
for (i in 1:l2) {if (sum(as.numeric((pos1==pos[i])&(chr1==chr[i])))==1) ftr2[i]<-TRUE}
# 




df1<-data.frame(pos=pos)
df2<-data.frame(pos=pos1)
df3<-merge(df1,df2,by="pos")
ftr4<-array(FALSE,dim=length(pos))
ftr5<-array(FALSE,dim=length(pos1))
pos2<-df3$pos
for (i in 1:length(pos))
{
	if (sum(as.numeric(pos2==pos[i]))==1) ftr4[i]<-TRUE
}
for (i in 1:length(pos1))
{
	if (sum(as.numeric(pos2==pos1[i]))==1) ftr5[i]<-TRUE
}
gen4<-gen2[,ftr4]
gen5<-gen3[,ftr5]
# 
z<-colSums(gen4!=gen5)
zn<-z/216
# 
gen6<-array(as.numeric(gen4>=1),dim=dim(gen4))
gen7<-array(as.numeric(gen5>=1),dim=dim(gen5))
zz<-1-colSums(array(as.numeric(gen4>=1 & gen5>=1),dim=dim(gen4)))/colSums(array(as.numeric(gen4>=1 | gen5>=1),dim=dim(gen4)))
zz[zz==NaN]<-0
#
zz1<-1-colSums(array(as.numeric(gen4==2 & gen5==2),dim=dim(gen4)))/colSums(array(as.numeric(gen4==2 | gen5==2),dim=dim(gen4)))
zz1[zz1==NaN]<-0
# 
ftr6<-zn>0.5
znz<-zn[ftr6]
pos6<-pos[ftr4]
pos7<-pos1[ftr5]
snp6<-snp[ftr4]
snp7<-snp1[ftr5]
#
pos8<-pos6[ftr6]
pos9<-pos7[ftr6]
snp8<-snp6[ftr6]
snp9<-snp7[ftr6]
# 

dd<-read.table("/home/smalov/BOTSWANA3/DATA/qc.csv",header=TRUE,sep=",")
read.gdsn


