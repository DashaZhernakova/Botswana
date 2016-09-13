library(gdsfmt)
library(SNPRelate)


dat.f<-openfn.gds("/Volumes/rset1/Botswana_chip/cleaned_all_final_no_outliers.gds")
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
	zd<-1-colSums(array(as.numeric((gen2.c==1|gen2.c==2) & (gen3.c==1|gen3.c==2)),dim=dm))/ colSums(array(as.numeric((gen2.c==1|gen2.c==2) | (gen3.c==1|gen3.c==2)),dim=dm))
	zd[(zd=="NaN"|is.na(zd))]<-0
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
	ZR<-c(ZR,zr)
}
d.out<-data.frame(chr=CHR,snp=SNP,pos=POS,all=ALL,snp1=SNP1,pos1=POS1,all1=ALL1,zp=ZP,zpa=ZPA,zd=ZD,zr=ZR)
write.table(d.out,file="/Volumes/rset1/Botswana_chip//qc.csv",row.names=FALSE,sep=",")

