#### No BioMart & From InterVar to ClinAnnotation ####
FromInterVarNoBioMart<-function(outdir="~/VarDetect/OneSample/Flow1/Flow1.1/Output3/",
sourcedir="/colossus/home/vorthunju/VarDetect/OneSample/Flow1/Flow1.1/Flow1.1.1/PART3/Source/"){

library(stringr)
library(jsonlite)

ENSP=paste0(sourcedir,"ENSP.txt")
RefENSP<-read.table(ENSP,sep="\t",header=T,stringsAsFactors=F)

temp<-read.table(paste0(outdir,"intervar.hg19_multianno.txt"),header=F,sep="\t",fill=T,stringsAsFactors=F)
colnames(temp)<-temp[1,]
temp<-temp[-1,]

temp2<-as.data.frame(str_split_fixed(strsplit(paste0(temp[temp$AAChange.ensGene!=".",]$AAChange.ensGene,collapse=","),",")[[1]], ":", 5))
colnames(temp2)<-c("ENSG","ENST","exon","c","p")
temp2<-data.frame(temp2[,1:2],ENSP="",temp2[3:5],CANO=FALSE,RefSeq.peptide.ID="",len=0,stringsAsFactors=F)

RefENSP<-RefENSP[RefENSP[,2] %in% temp2[,2],]

for(i in RefENSP[,2]){
 	temp2[temp2$ENST==i,3]<-as.character(RefENSP[RefENSP[,2]==i,1])
	temp2[temp2$ENST==i,8]<-as.character(RefENSP[RefENSP[,2]==i,4])
	temp2[temp2$ENST==i,9]<-as.character(RefENSP[RefENSP[,2]==i,6])
	temp2[temp2$ENST==i,7]<-RefENSP[RefENSP[,2]==i,3]
}

temp2<-data.frame(temp2,CDS_position=0)
temp2<-temp2[temp2$ENSP != "",]

idx<-data.frame(Gene=temp2$ENSG,Feature=temp2$ENST,SYMBOL="",HGNC_ID="-",
		CANONICAL=temp2$CANO,CCDS="-",ENSP=temp2$ENSP,Ref_Seq=temp2$RefSeq.peptide.ID,len=temp2$len,stringsAsFactors=F)
idx<-idx[!duplicated(idx$ENSP),]

for(i in idx[!duplicated(idx[,1]),1]) idx[idx[,1]==i,3] <- as.character(temp[temp$Gene.ensGene == i,]$Gene.knownGene[1])


tryCatch(file.copy(paste(sourcedir,idx$ENSP,".json",sep=""),outdir,overwrite=T),warning=function(x){return("")})
tryCatch(file.copy(paste(sourcedir,"PTM_",idx$ENSP,".json",sep=""),outdir,overwrite=T),warning=function(x){return("")})

temp2$c<-gsub("[A-Z]","",gsub("[a-z]\\.[A-Z]","",temp2$c))
temp2$p<-gsub("[A-Z]","",gsub("[a-z]\\.[A-Z]","",temp2$p))

input<-data.frame(ENST=temp2$ENST,ENSP=temp2$ENSP,Protein_position=temp2$p,Consequence="",IMPACT=0.5,stringsAsFactors=F)


for(i in 1:length(input[,1])){
	a=grep(input$ENST[i],temp$AAChange.ensGene)
	b=grep(input$Protein_position[i],temp$AAChange.ensGene)
	a<-a[a %in% b]
	input[i,4] <- temp$ExonicFunc.knownGene[a]
}

input[grep("frameshift",input$Consequence),]$IMPACT<-input[grep("frameshift",input$Consequence),]$IMPACT+5
input[grep("stop",input$Consequence),]$IMPACT<-input[grep("stop",input$Consequence),]$IMPACT+4
input[grep("nonframe",input$Consequence),]$IMPACT<-input[grep("nonframe",input$Consequence),]$IMPACT+3
input[grep("nonsyn",input$Consequence),]$IMPACT<-input[grep("nonsyn",input$Consequence),]$IMPACT+2

for(i in 1:length(idx[,1])){
	tmp<-input[input$ENSP==idx$ENSP[i],]
	tmp<-data.frame(coord=tmp$Protein_position,category=tmp$Consequence,value=tmp$IMPACT,stringsAsFactors=F)
	temp2<-tryCatch(fromJSON(paste(outdir,"PTM_",idx$ENSP[i],".jxson",sep="")),error=function(x){return(tmp[0,])})
	write(toJSON(rbind(tmp,temp2)),paste(outdir,"PTM_",idx$ENSP[i],".json",sep=""))
}

temp<-
data.frame(
Chr=temp$Chr,Start=temp$Start,End=temp$End,Ref=temp$Ref,Alt=temp$Alt,ExonicFunc.refGene=temp$ExonicFunc.refGene,AAChange.ensGene=temp$AAChange.ensGene,avsnp147=temp$avsnp147,gnomAD_genome_ALL=temp$gnomAD_genome_ALL,SIFT_pred=temp$SIFT_pred,Polyphen2_HDIV_pred=temp$Polyphen2_HDIV_pred,Polyphen2_HVAR_pred=temp$Polyphen2_HVAR_pred,LRT_pred=temp$LRT_pred,MutationTaster_pred=temp$MutationTaster_pred,MutationAssessor_pred=temp$MutationAssessor_pred,FATHMM_pred=temp$FATHMM_pred,PROVEAN_pred=temp$PROVEAN_pred,MetaSVM_pred=temp$MetaSVM_pred,MetaLR_pred=temp$MetaLR_pred,stringsAsFactors=F)															


tempx<-cbind(input[0,],temp[0,])
for(i in 1:length(input[,1])){
	a=grep(input$ENST[i],temp$AAChange.ensGene)
	b=grep(input$Protein_position[i],temp$AAChange.ensGene)
	a<-a[a %in% b]
	tempx<-rbind(tempx,cbind(input[i,],temp[a,]))
}

tempx<-tempx[,-c(11,12)]

write.table(tempx,paste(outdir,"out.Part3.vep.txt",sep=""),sep="\t",quote=F,row.names=F)

#return(idx)
return(rbind(tmp,temp2))
}

#################### PART2 Reporting ACMG ################################################

ACMG<-function(indir="/colossus/home/vorthunju/VarDetect/OneSample/Flow1/Flow1.1/Input/",outdir="/colossus/home/vorthunju/VarDetect/OneSample/Flow1/Flow1.1/Output/",sourcedir="/colossus/home/vorthunju/VarDetect/OneSample/Flow1/Flow1.1/Source/"
){
	##LoadFiles
	intervar<-read.table(paste(outdir,"intervar.hg19_multianno.txt.intervar",sep=""),sep="\t")
	#clnid<-read.table(paste(sourcedir,"ClinVar/ClnVarID.txt",sep=""),sep="\t",header=T)
	#clnsbmt<-read.csv(paste(sourcedir,"ClinVar/ClnVarSubmitted.csv",sep=""),header=T)
	inp<-read.table(paste(indir,"Part3.vep.txt",sep=""),sep="\t",header=T)




	########################################### Split REF files ###################################
	 ###Split ClnVarID###
	#temp<-gsub("chr","",gsub("_.*","",clnid[,1]))
	#temp<-temp[!is.na(temp)]
	#temp<-temp[!duplicated(temp)]
	#
	#for(i in temp){
	#	dir.create(paste(sourcedir,"ClinVar/",i,sep=""))
	#	for(j in 1:9) write.table(clnid[gsub("chr","",gsub("_.*","",clnid[,1]))==i & as.numeric(substr(gsub(".*:","",clnid[,2]),1,1)) == j,],paste(sourcedir,"ClinVar/",i,"/",j,".txt",sep=""),sep="\t",quote=F,row.names=F)
	#}
	# ###Split ClnVarSubmitted####
	#	write.table(clnsbmt[nchar(clnsbmt[,1])<3,],paste(sourcedir,"ClinVar/ClnVarSubmitted/l3.txt",sep=""),sep="\t",quote=F,row.names=F)
	#for(i in 3:6){
	#	dir.create(paste(sourcedir,"ClinVar/ClnVarSubmitted/",i,sep=""))
	#	temp<-clnsbmt[nchar(clnsbmt[,1])==i,]
	#	for(j in 0:(10^(i-2)-1)) write.table(temp[as.numeric(substr(temp[,1],nchar(temp[,1])-i+3,nchar(temp[,1])))==j,],paste(sourcedir,"ClinVar/ClnVarSubmitted/",i,"/",j,".txt",sep=""),sep="\t",quote=F,row.names=F)
	#}
	#
	#############################################################################################

	##InterVar

	intervar<-data.frame(id = paste(paste("chr",intervar[,1],sep=""),intervar[,2],
			paste(intervar[,4],intervar[,5],sep="/"),sep="_"),
			InterVar=gsub(" PVS1.*","",gsub(".*rVar: ","",intervar[,14]))
			,PSV1=gsub(" .*","",gsub(".*PVS1=","",intervar[,14]))
			,PS=gsub("] .*","]",gsub(".*PS=","",intervar[,14]))
			,PM=gsub("] .*","]",gsub(".*PM=","",intervar[,14]))
			,PP=gsub("] .*","]",gsub(".*PM=","",intervar[,14]))
			,BA1=gsub(" .*","",gsub(".*BA1=","",intervar[,14]))
			,BS=gsub("] .*","]",gsub(".*BS=","",intervar[,14]))
			,BP=gsub("] .*","]",gsub(".*BS=","",intervar[,14]))
			,ACMG="ACMG:",stringsAsFactors=F)

	for(i in 1:length(intervar[,1])){
		temp<-intervar[i,]
		intervar$ACMG[i]<-paste(if(temp$PSV1==1){"PSV1;"}
					,if(gsub("1.*","",temp$PS)!=temp$PS){
						paste("PS",nchar(gsub("1.*","",temp$PS)) - nchar(gsub(",","",gsub("1.*","",temp$PS))) + 1,";",sep="")}
					,if(gsub("1.*","",temp$PM)!=temp$PM){
						paste("PM",nchar(gsub("1.*","",temp$PM)) - nchar(gsub(",","",gsub("1.*","",temp$PM))) + 1,";",sep="")}
					,if(gsub("1.*","",temp$PP)!=temp$PP){
						paste("PP",nchar(gsub("1.*","",temp$PP)) - nchar(gsub(",","",gsub("1.*","",temp$PP))) + 1,";",sep="")}
					,if(temp$BA1==1){"BA1;"}
					,if(gsub("1.*","",temp$BS)!=temp$BS){
						paste("BS",nchar(gsub("1.*","",temp$BS)) - nchar(gsub(",","",gsub("1.*","",temp$BS))) + 1,";",sep="")}
					,if(gsub("1.*","",temp$BP)!=temp$BP){
						paste("BP",nchar(gsub("1.*","",temp$BP)) - nchar(gsub(",","",gsub("1.*","",temp$BP))) + 1,";",sep="")},sep="")
		}


	write.table(intervar,paste(outdir,"Intervar.txt",sep=""),sep="\t",quote=F,row.names=F)

	##LOVD
	lovd<-data.frame(id=inp[,1],inp$LOVD,stringsAsFactors=F)
	lovd<-lovd[lovd[,2]!="-",]
	lovd<-lovd[!duplicated(lovd[,1]),]
	write.table(lovd,paste(outdir,"LOVD.txt",sep=""),sep="\t",quote=F,row.names=F)

	##ClinVar

	##### Load Ref ClnID #####
	chr<-gsub("chr","",gsub("_.*","",intervar[,1]))
	chr<-chr[!duplicated(chr)]
	clnid<-data.frame(id="",id2="",id3="",id4="",stringsAsFactors=F)[0,]
	for(i in chr){
	temp<-intervar[gsub("chr","",gsub("_.*","",intervar[,1]))==i,]
	temp<-substr(gsub(".*_","",gsub("_[A-Z].*","",temp[,1])),1,1)
	temp<-temp[!duplicated(temp)]
	for(j in temp)	clnid<-rbind(clnid,read.table(paste(sourcedir,"ClinVar/ClnVarID/",i,"/",j,".txt",sep=""),sep="\t",header=T))
	}

	id<-inp[inp[,1] %in% clnid[,1],1]
	id<-as.character(id[!duplicated(id)])
	clnid<-clnid[clnid[,1] %in% id,]

	#### Load Ref clnsbmt ####
	cnt<-nchar(clnid$id3)
	cnt<-cnt[!duplicated(cnt)]

	temp<-data.frame(id3="",sig="",pheno="",descrpt="",stringsAsFactors=F)[0,]
	if(length(cnt[cnt<3])>0) temp<-rbind(temp,read.table(paste(sourcedir,"ClinVar/ClnVarSubmitted/1-2.txt",sep=""),sep="\t",header=T))
	if(length(cnt[cnt>2])>0){
	for(i in cnt[cnt>2]){
	nmb<-as.numeric(substr(clnid[nchar(clnid$id3)==i,3],3,i))
	nmb<-nmb[!duplicated(nmb)]
	for(j in nmb) temp<-rbind(temp,read.table(paste(sourcedir,"ClinVar/ClnVarSubmitted/",i,"/",j,".txt",sep=""),sep="\t",header=T))
		}
	}

	clnsbmt<-temp

	clnsbmt<-clnsbmt[clnsbmt[,1] %in% clnid[,3],]

	write.table(clnid,paste(outdir,"ClnVarID.txt",sep=""),sep="\t",quote=F,row.names=F)
	write.table(clnsbmt,paste(outdir,"ClnVarDetail.txt",sep=""),sep="\t",quote=F,row.names=F)

	###InputIndex
	inp<-inp[inp[,1] %in% clnid[,1],]
	inp<-inp[!duplicated(inp[,1]),c(1,19)]

	write.table(inp,paste(outdir,"Part3.index.txt",sep=""),sep="\t",quote=F,row.names=F)

	#clnid<-data.frame(id=paste(paste("chr",var_sum[,19],sep=""),var_sum[,20],
	#			paste(var_sum[,22],var_sum[,23],sep="/"),sep="_"),
	#	id2=paste(var_sum[,17],": Chr",var_sum[,19],":",var_sum[,20],"-",var_sum[,21],sep=""),
	#	id3=var_sum[,31], id4=var_sum[,3])
	#write.table(clnid,paste(clnpth,"ClnVarID.txt",sep=""),sep="\t",row.names=F,quote=F)
	#
	#clnsbmt<-data.frame(id3=clnsbmt[,1],sig=clnsbmt[,2],pheno=clnsbmt[,6],descrpt=clnsbmt[,4],stringsAsFactors=F)

}

#####GET COORDINATE FROM ENTREZ####
SelectGene<-function(genes="11079;9374;10928",sourcedir="/colossus/home/vorthunju/VarDetect/OneSample/Flow1/Flow1.1/Flow1.1.1/PART3/Source/"){
	if(genes=="ALL") return("")
	if(genes!="ALL"){
		#bcftools view VarDetect/tempfile/Dent_32exomes_Unrelated.raw.g.vcf.gz -r chr1:2323267-2336883,chr6:32121218-32134011,chr18:9475007-9538114
		coord<-read.table(paste(sourcedir,"entrez_local.txt",sep=""),sep="\t",header=T)
		temp<-unlist(strsplit(genes,";"))
		temp<-coord[coord[,5] %in% temp,]
		return(cat(paste("chr",temp[,2],":",temp[,3],"-",temp[,4],sep="",collapse=",")))
	}
}

########################################GETPEPTIDE###########################################################################################################
getpeptide<-function(file="inputVCF.vep.cano.exon.dbnsfp.txt",

postt="/colossus/home/vorthunju/VarDetect/tools/HPRD/FLAT_FILES_072010/POST_TRANSLATIONAL_MODIFICATIONS.txt"
,
hprd="/colossus/home/vorthunju/VarDetect/tools/HPRD/FLAT_FILES_072010/PROTEIN_ARCHITECTURE.txt"
,
pfmmp="/colossus/home/vorthunju/VarDetect/tools/HPRD/pdb_pfam_mapping.txt",outdir="/colossus/home/vorthunju/VarDetect/OneSample/Flow1/Flow1.1/Flow1.1.1/PART3/Output/")
{

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos='http://cran.us.r-project.org')
if (!requireNamespace("biomaRt", quietly = TRUE))
    install.packages("biomaRt",repos='http://cran.us.r-project.org')
if (!requireNamespace("jsonlite", quietly = TRUE))
    install.packages("jsonlite",repos='http://cran.us.r-project.org')


library("biomaRt")
library("jsonlite")

##UploadFiles##
pfm<-read.table(pfmmp,sep="\t",header=T)
ptm<-read.table(postt,sep="\t",header=F)
#hpr<-read.table(hprd,sep="\t",head=F)
inp<-read.table(file,sep="\t",header=T)
#input<-inp$V20

##PrepareBioMartParameters##
input=as.character(inp[!(duplicated(inp$ENSP) | inp$ENSP=="-"),]$ENSP)
can<-inp[,c(27,23)]
can<-can[can[,2] == "YES" & can[,1] != "-",]
can<-can[!duplicated(can[,1]),]
grch37 = useEnsembl(biomart="ensembl",GRCh=37,dataset="hsapiens_gene_ensembl")
temp2<-substr(listAttributes(grch37)[,1],1,4)
interpro<-listAttributes(grch37)[temp2=="inte",1]
pfam<-listAttributes(grch37)[temp2=="pfam",1]

##PullDataFromBioMart##
interprox<-getBM(attributes = c("hgnc_symbol","refseq_peptide","ensembl_peptide_id",interpro), filters = "ensembl_peptide_id" , values = input, mart = grch37)
pfamx<-getBM(attributes = c("hgnc_symbol","refseq_peptide","ensembl_peptide_id",pfam), filters = "ensembl_peptide_id" , values = input, mart = grch37 )
refens<-getBM(attributes = c("hgnc_symbol","refseq_peptide","ensembl_peptide_id","peptide"), filters = "ensembl_peptide_id" , values = input, mart = grch37)

refens<-refens[refens[,4] != "",]
refens$refseq_peptide<- gsub("\\;.*","",refens$refseq_peptide)
#refens$refseq_peptide<-gsub(".*;","",refens$refseq_peptide)

##PrepareProteinLength## 
pl<-data.frame(hgnc_symbol=refens[,3],refseq_peptide=refens[,4],ensembl_peptide_id=refens[,1],source="-",id="-",short="-",long="-",start=rep(1,length(refens[,1])),end=nchar(refens[,2])-1,stringsAsFactors=FALSE)


##PreparePfamData##
pfamx<-data.frame(pfamx[,1:3],source="pfam",id=pfamx[,4],short="-",long="-",start=pfamx[,5],end=pfamx[,6],stringsAsFactors=FALSE)

tmp<-data.frame(pfm[gsub("[.]\\d+","",pfm$PFAM_ACC) %in% pfamx[!duplicated(pfamx$id),]$id,],stringsAsFactors=FALSE)

tmp<-tmp[!duplicated(tmp$PFAM_ACC),5:7]

tmp$PFAM_ACC<-gsub("[.]\\d+","",tmp$PFAM_ACC)

for(i in tmp[,1]){
	pfamx[as.character(pfamx$id) == i,6]<-as.character(tmp[tmp[,1]==i,2])
	pfamx[as.character(pfamx$id) == i,7]<-as.character(tmp[tmp[,1]==i,3])
}

pfamx<-pfamx[pfamx$id %in% tmp$PFAM_ACC,]

##PrepareInterProData##
interprox<-interprox[interprox[,1]!=interprox[,5],]
interprox<-interprox[!duplicated(interprox[,4]),]
interprox<-data.frame(interprox[,1:3],source="interpro",id=interprox[,4],short=interprox[,5],long=interprox[,6],start=interprox[,7],end=interprox[,8],stringsAsFactors=FALSE)

##PreparePostTranslationalModificatonSites##
temp<-ptm[ptm[,4] %in% refens[,4],]
ptm<-data.frame(hgnc_symbol=temp[,2],refseq_peptide=temp[,4],ensembl_peptide_id="",source="PTM",id=gsub("[,]","|",temp[,11]),short=temp[,7],long=temp[,9],start=temp[,5],end=NA,stringsAsFactors=FALSE)

for(i in 1:length(refens[,3]))
	ptm[ptm[,2]==refens[i,2],3]=rep(refens[i,3],length(ptm[ptm[,2]==refens[i,2],3]))

colnames(pfamx)<-colnames(pl)
colnames(interprox)<-colnames(pl)
colnames(ptm)<-colnames(pl)

rtn<-rbind(pfamx,interprox,ptm)

temp<-rtn[!duplicated(rtn[,2:3]) & rtn[,3]!="",2:3]
for(i in 1:length(temp[,1]))
	rtn[rtn[,2]==temp[i,1],3]<-temp[i,2]

##Create Json Coordinate&Mutation##
inp<-cbind(inp,0)
inp[inp$IMPACT=="MODIFIER",65]=1
inp[inp$IMPACT=="LOW",65]=2
inp[inp$IMPACT=="MODERATE",65]=4
inp[inp$IMPACT=="HIGH",65]=8

rtn[rtn[,3] %in% can[,1],3]<-paste("",rtn[rtn[,3] %in% can[,1],3],sep="")
inp$ENSP<-as.character(inp$ENSP)
inp[inp$ENSP %in% can[,1],]$ENSP<-paste("",inp[inp$ENSP %in% can[,1],]$ENSP,sep="")


temp<-rtn[!duplicated(rtn[,2]),3]

idx<-inp[0,]

for(i in 1:length(temp)){
	temp2<-rtn[rtn[,3]==temp[i] & !is.na(rtn$end),]
	temp3<-rtn[rtn[,3]==temp[i] & is.na(rtn$end),]
	temp4<-data.frame(coord=inp[inp$ENSP==temp[i],10],category=inp[inp$ENSP==temp[i],7],value=inp[inp$ENSP==temp[i],65],stringsAsFactors=FALSE)
	if(length(temp3$start)!=0){
		write(toJSON(data.frame(name=temp2$long,coord=paste(temp2$start,"-",temp2$end,sep=""))),
		paste(outdir,temp[i],".json",sep=""))
		write(toJSON(rbind(data.frame(coord=temp3$start,category=temp3$long,value=-1.4,stringsAsFactors=FALSE),temp4)),paste(outdir,"PTM_",temp[i],".json",sep=""))
		idx<-rbind(idx,inp[inp$ENSP==temp[i],])
	}
}

idx<-idx[!duplicated(idx$ENSP),]

idx<-data.frame(idx[,c(4:5,19,21,23,26,27)],Ref_Seq="",len=0,stringsAsFactors=F)

pl[pl[,3] %in% can[,1],3]<-paste("",pl[pl[,3] %in% can[,1],3],sep="")

for(i in idx$ENSP){
	 idx[idx$ENSP==i,]$Ref_Seq <- pl[pl$ensembl_peptide_id==i,]$refseq_peptide
	 idx[idx$ENSP==i,]$len <- pl[pl$ensembl_peptide_id==i,]$end
}


write.table(idx,paste(outdir,"index.txt",sep=""),sep="\t",row.names=F,quote=F)
outp<-inp[,c(1:5,7,10,11,13,14,19,23,26,27,32,33,38:55,57,60)]
outp<-outp[!duplicated(outp[,c(7,14)]),]
outp$PUBMED <- as.character(outp$PUBMED)
outp$PUBMED <- gsub("[,]",":",outp$PUBMED)

write.table(outp,paste(outdir,"out.Part3.vep.txt",sep=""),row.names=F,sep="\t",quote=F)


###BindResultsTogether###
return(rtn)

}
#################################################################################
##########################################   NO BIOMART   ######################################################################################################
##################################################################################

getpeptideNobiomaRt<-function(

ENSP="/colossus/home/vorthunju/VarDetect/OneSample/Flow1/Flow1.1/Flow1.1.1/PART3/Source/ENSP.txt"
,
file="/colossus/home/vorthunju/VarDetect/OneSample/Flow1/Flow1.1/Flow1.1.1/PART3/Output/Part3.vep.txt"
,
outdir="/colossus/home/vorthunju/VarDetect/OneSample/Flow1/Flow1.1/Flow1.1.1/PART3/Output2/"
,
sourcedir="/colossus/home/vorthunju/VarDetect/OneSample/Flow1/Flow1.1/Flow1.1.1/PART3/ENSP/"

)
{

if (!requireNamespace("jsonlite", quietly = TRUE))
    install.packages("jsonlite",repos='http://cran.us.r-project.org')

library("jsonlite")

##UploadFiles##

RefENSP<-read.table(ENSP,sep="\t",header=T,stringsAsFactors=F)

input<-read.table(file,sep="\t",header=T,stringsAsFactors=F)
input<-input[input$CDS_position!="-",]

idx<-data.frame(Gene=input$Gene,Feature=input$Feature,SYMBOL=input$SYMBOL,HGNC_ID=input$HGNC_ID,CANONICAL=input$CANONICAL,CCDS=input$CCDS,ENSP=input$ENSP,Ref_Seq="",len=0,stringsAsFactors=F)


idx<-idx[!duplicated(idx$ENSP),]
idx<-idx[idx$ENSP %in% RefENSP[,1],]


tryCatch(file.copy(paste(sourcedir,idx$ENSP,".json",sep=""),outdir,overwrite=T),warning=function(x){return("")})
tryCatch(file.copy(paste(sourcedir,"PTM_",idx$ENSP,".json",sep=""),outdir,overwrite=T),warning=function(x){return("")})


for(i in 1:length(idx[,1])){
	idx[i,8]<-RefENSP[RefENSP[,1]==idx[i,7],2][1]
	idx[i,9]<-RefENSP[RefENSP[,1]==idx[i,7],4][1]
	temp<-input[input$ENSP==idx$ENSP[i],]
	temp<-data.frame(coord=temp$Protein_position,category=temp$Consequence,value=temp$IMPACT,stringsAsFactors=F)
	temp[temp$value=="MODIFIER",3]=1
	temp[temp$value=="LOW",3]=2
	temp[temp$value=="MODERATE",3]=4
	temp[temp$value=="HIGH",3]=8
	temp2<-tryCatch(fromJSON(paste(outdir,"PTM_",idx$ENSP[i],".jxson",sep="")),error=function(x){return(temp[0,])})
	write(toJSON(rbind(temp,temp2)),paste(outdir,"PTM_",idx$ENSP[i],".json",sep=""))
}

write.table(data.frame(X_Uploaded_variation=input$X_Uploaded_variation,Location=input$Location,Allele=input$Allele,Gene=input$Gene,Feature=input$Feature,Consequence=input$Consequence,Protein_position=input$Protein_position,Amino_acids=input$Amino_acids,Existing_variation=input$Existing_variation,IMPACT=input$IMPACT,SYMBOL=input$SYMBOL,CANONICAL=input$CANONICAL,CCDS=input$CCDS,ENSP=input$ENSP,SIFT=input$SIFT,PolyPhen=input$PolyPhen,AF=input$AF,AFR_AF=input$AFR_AF,AMR_AF=input$AMR_AF,EAS_AF=input$EAS_AF,EUR_AF=input$EUR_AF,SAS_AF=input$SAS_AF,AA_AF=input$AA_AF,EA_AF=input$EA_AF,gnomAD_AF=input$gnomAD_AF,gnomAD_AFR_AF=input$gnomAD_AFR_AF,gnomAD_AMR_AF=input$gnomAD_AMR_AF,gnomAD_ASJ_AF=input$gnomAD_ASJ_AF,gnomAD_EAS_AF=input$gnomAD_EAS_AF,gnomAD_FIN_AF=input$gnomAD_FIN_AF,gnomAD_NFE_AF=input$gnomAD_NFE_AF,gnomAD_OTH_AF=input$gnomAD_OTH_AF,gnomAD_SAS_AF=input$gnomAD_SAS_AF,MAX_AF=input$MAX_AF,CLIN_SIG=input$CLIN_SIG,PUBMED=input$PUBMED),paste(outdir,"out.Part3.vep.txt",sep=""),sep="\t",quote=F,row.names=F)

return(idx)
}

##########################################GETTRANSCRIPT######################################################################################################
gettranscript<-function(wkdir="/colossus/home/vorthunju/VarDetect/AfterVCF/Output/",file="inputVCF.vep.cano.exon.dbnsfp.txt",
postt="/colossus/home/vorthunju/VarDetect/tools/HPRD/FLAT_FILES_072010/POST_TRANSLATIONAL_MODIFICATIONS.txt",
hprd="/colossus/home/vorthunju/VarDetect/tools/HPRD/FLAT_FILES_072010/PROTEIN_ARCHITECTURE.txt",
pfmmp="/colossus/home/vorthunju/VarDetect/tools/HPRD/pdb_pfam_mapping.txt")
{
library("biomaRt")
setwd(wkdir)
input<-read.table(file,sep="\t",header=F,stringsAsFactors=FALSE)[,4:5]
input<-input[!duplicated(input),]

TF = useEnsembl(biomart="ENSEMBL_MART_FUNCGEN",GRCh=37,dataset="hsapiens_motif_feature")
OTH = useEnsembl(biomart="ENSEMBL_MART_FUNCGEN",GRCh=37,dataset="hsapiens_external_feature")
MI = useEnsembl(biomart="ENSEMBL_MART_FUNCGEN",GRCh=37,dataset="hsapiens_mirna_target_feature")
grch37 = useEnsembl(biomart="ensembl",GRCh=37,dataset="hsapiens_gene_ensembl")


####SetAttributes####
ginputAtr<-c(listAttributes(grch37)[c(1,3,9,10,11,59),1],"strand")
atTF<-listAttributes(TF)[,1]
atOTH<-listAttributes(OTH)[,1]
atMI<-listAttributes(MI)[,1]
atEx = c("ensembl_gene_id","ensembl_transcript_id","chromosome_name", "exon_chrom_start", 
            "exon_chrom_end", "external_gene_name", "strand")
######################


#####GetInputPosition#####
input<-getBM(attributes = ginputAtr, filters = "ensembl_transcript_id", values = input[,2], mart = grch37)
input[,4]<-input[,4]-1000
input[,5]<-input[,5]+1000
inv<-list(input[,3],input[,4],input[,5])
##########################


####GetTFdata#############
TF<-getBM(attributes = atTF, filters = c("chromosome_name","start","end"),
values = inv , mart = TF)
TF<-data.frame(ensembl_gene_id="",ensembl_transcript_id="",chromosome_name=TF$chromosome_name,chromosome_start=TF$chromosome_start,
chromosome_end=TF$chromosome_end,source="TF",id=TF$binding_matrix_id,short=TF$so_name,long=TF$display_label,other=""
,stringsAsFactors=FALSE)

####GetOtherdata#########
OTH<-getBM(attributes = atOTH, filters = c("chromosome_name","start","end"),
values = inv, mart = OTH)

OTH<-data.frame(ensembl_gene_id="",ensembl_transcript_id="",chromosome_name=OTH$chromosome_name,chromosome_start=OTH$chromosome_start,
chromosome_end=OTH$chromosome_end,source=OTH$feature_type,id=OTH$display_label,
short=OTH$feature_type_class,long=OTH$feature_type_description,other="",stringsAsFactors=FALSE)
##########################


####GetmiRNAdata#########
MI<-getBM(attributes = atMI, filters = c("chromosome_name","start","end"),
values = inv, mart = MI)

MI<-data.frame(ensembl_gene_id="",ensembl_transcript_id="",chromosome_name=MI$chromosome_name,
chromosome_start=MI$chromosome_start,chromosome_end=MI$chromosome_end,
source=MI$linkage_annotation,id=MI$accession,short=MI$feature_type_description,
long=MI$display_label,other="",stringsAsFactors=FALSE)
#########################

####GetExonPosition#######

EX<-getBM(attributes = atEx, filters = "ensembl_transcript_id",
values = input[,2], mart = grch37)

EX<-data.frame(ensembl_gene_id=EX[,1],ensembl_transcript_id=EX[,2],chromosome_name=EX[,3],
chromosome_start=EX[,4],chromosome_end=EX[,5],
source=EX[,6],id="-",short="EXON",
long=EX$strand,other="",stringsAsFactors=FALSE)

####GetGENEPOSITION######
GN<-data.frame(input[,1:3],chromosome_start=input[,4],chromosome_end=input[,5],
source=input[,6],id="-",short="-",
long=input$strand,other="",stringsAsFactors=FALSE)
#########################


###MERGE&CLEAN###
rtn<-rbind(GN,EX,TF,OTH,MI)
rtn<-rtn[!duplicated(rtn$id) | rtn$id=="-",]

for(i in 1:length(input[,1])){
rtn[rtn$chromosome_name==input[i,3] & rtn$chromosome_start>input[i,4] & 
rtn$chromosome_end<input[i,5],1]<-as.character(input[i,1])
rtn[rtn$chromosome_name==input[i,3] & rtn$chromosome_start>input[i,4] & 
rtn$chromosome_end<input[i,5],2]<-as.character(input[i,2])
}

rtn<-data.frame(rtn,start=0,end=0)
for(i in 1:length(GN[,1])){
rtn[rtn[,1]==GN[i,1],]$start <- rtn[rtn[,1]==GN[i,1],4] - GN[i,4] + 1
rtn[rtn[,1]==GN[i,1],]$end <- rtn[rtn[,1]==GN[i,1],5] - GN[i,4] + 1
}
#################

return(rtn[rtn[,1] %in% input[,1],])
}
