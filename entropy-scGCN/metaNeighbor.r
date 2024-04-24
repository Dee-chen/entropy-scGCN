# **********************************************************
# * Author        : YanYuting
# * Email         : yanyuting@genomics.cn
# * Create time   : 2022-07-31 16:10
# * Filename      : metaNeighbor.r
# * Description   : 
# **********************************************************
library(Rcpp)
Rcpp::sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix asMatrix(NumericVector rp,
                       NumericVector cp,
                       NumericVector z,
                       int nrows,
                       int ncols){

  int k = z.size() ;

  NumericMatrix  mat(nrows, ncols);

  for (int i = 0; i < k; i++){
      mat(rp[i],cp[i]) = z[i];
  }

  return mat;
}
' )
as_matrix <- function(mat){

  row_pos <- mat@i
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])

  tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
                 nrows =  mat@Dim[1], ncols = mat@Dim[2])

  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}
downsample<-function(object,n=5000){
        types<-table(object$cluster)
        Idents(object)<-object$cluster
        if(types[[1]]>=n){
            all_obj<-subset(object,subset=cluster==names(types[1]),downsample = n)
        }else{
            all_obj<-subset(object,subset=cluster==names(types[1]))
        }
        if(length(types)==1){
            return(all_obj)
        }else{
            for(i in (2:length(types))){
                if(types[[i]]>=n){
                    sub_obj<-subset(object,subset=cluster==names(types[i]),downsample = n)
                }else{
                    sub_obj<-subset(object,subset=cluster==names(types[i]))
                }
                all_obj<-merge(all_obj,sub_obj)
            }
            return(all_obj)
        }
    }
args=commandArgs(T)
myMetaNeighbor<-function(refFile,qurFile,refCluster,qurClusterFile,geneTable){
	suppressMessages(library(Seurat))
	suppressMessages(library(entropy))
	suppressMessages(library(Matrix))
	suppressMessages(library(dplyr))
	library(MetaNeighbor)
	library(SummarizedExperiment)
	library("magrittr")
	library(gplots)
	library(RColorBrewer)
	ref<-readRDS(refFile)
	qur<-readRDS(qurFile)
	tmp.cluster<-ref@meta.data[,refCluster]
	tmp.cluster<-as.vector(unlist(tmp.cluster))
	ref$cluster<-tmp.cluster
	ref<-downsample(ref,n=3500)
	ref.counts<-ref$RNA@counts
	qur.counts<-qur$RNA@counts
	ref.counts<-as_matrix(ref.counts)
	qur.counts<-as_matrix(qur.counts)
	ref.genes<-rownames(ref.counts)
	qur.genes<-rownames(qur.counts)
	ref.genes<-toupper(ref.genes)
	qur.genes<-toupper(qur.genes)
	inter.genes<-intersect(ref.genes,qur.genes)
	ref.counts<-ref.counts[rownames(ref.counts)%in%inter.genes,]
	if(geneTable=="F"){
        	if(length(inter.genes)>=2000){
			nf=2000
        	}else{
            		nf=nrow(ref.counts)
        	}
        	data<-ref.counts
        	label<-ref@meta.data[,refCluster]
        	M <- nrow(data);
        	new.label <- label
        	pv1 <- sapply(1:M, function(i){
        	mydataframe <- data.frame(y=as.numeric(data[i,]), ig=new.label)
        	fit <- aov(y ~ ig, data=mydataframe)
        	summary(fit)[[1]][["Pr(>F)"]][1]
        	})
        	names(pv1) <- rownames(data)
        	pv1.sig <- names(pv1)[order(pv1)[1:nf]]
        	sel.features <- unique(pv1.sig)
    	}else{
        	use.genes<-read.table(genesFile)
        	use.genes<-use.genes$V1
        	use.genes<-toupper(use.genes)
        	sel.features<-intersect(inter_genes,use.genes)
        	sel.features<-unique(sel.features)
    	}
	ref.counts<-ref.counts[rownames(ref.counts)%in%sel.features,]
	qur.counts<-qur.counts[rownames(qur.counts)%in%sel.features,]
	all.counts<-cbind(ref.counts,qur.counts)
	print(str(all.counts))
	ref.labels<-ref@meta.data[,refCluster]
	qur.labels<-read.table(qurClusterFile)
	qur.labels<-qur.labels$V2
	all.labels<-c(ref.labels,qur.labels)
	print("all.counts and all.labels has done!")
	ref.study<-rep("ref",ncol(ref))
	qur.study<-rep("qur",ncol(qur))
	all.study<-c(ref.study,qur.study)
	ref.samples<-colnames(ref)
	qur.samples<-colnames(qur)
	all.samples<-c(ref.samples,qur.samples)
	pheno<-data.frame(Sample_ID=all.samples)
	#pheno$Sample_ID<-all.samples
	pheno$Study_ID<-all.study
	pheno$Celltype<-all.labels
	print(str(pheno))
	indatest<-SummarizedExperiment(assays=all.counts,colData= pheno)
	celltype_NV = MetaNeighborUS(var_genes = sel.features,
					dat = indatest,
					study_id = pheno[,2],
					cell_type = pheno[,3],
					fast_version = TRUE)
	print("celltype_NV has done")
	cols=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(100))
	breaks=seq(0,1,length=101)
	pdf("MetaNeighbora.sel.features.pdf",height = 10,width = 12)
	heatmap.2(celltype_NV,trace="none",density.info="none",col=cols,breaks=breaks,cexRow=0.6,cexCol=0.6)
	dev.off()
	print("all has done!")
}
myMetaNeighbor(args[1],args[2],args[3],args[4],args[5])

