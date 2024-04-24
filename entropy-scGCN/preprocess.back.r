suppressMessages(library(Seurat));
suppressMessages(library(entropy));
suppressMessages(library(Matrix));
suppressMessages(library(dplyr));
args=commandArgs(T)
preprocess_data<-function(ref_dir,qur_dir,ref_cluster,qur_cluster,chunks=10000){
    chunks<-chunks
    ref<-readRDS(ref_dir)
    qur<-readRDS(qur_dir)
#测试matrix分块处理解决数据过大无法转换为matrix的问题
# Created by: yanyuting
# Created on: 2022/5/9
    deal_sparse_matrix <- function(x, rmatrix, chunk = chunks){
        passes <- nrow(x) %/% chunk
        remaining <- nrow(x) %% chunk
        if(passes > 0){
            inx <- seq_len(chunk)
            yy <- x[inx, , drop = FALSE]
            yy <- as.matrix(yy)
            passes <- passes - 1L
            for(i in seq_len(passes)){
                inx <- inx + chunk
                y <- x[inx, , drop = FALSE]
                y <- as.matrix(y)
                yy <- rbind(yy,y)
            }
            if(remaining > 0){
                inx <- inx + remaining
                y <- x[inx, , drop = FALSE]
                y <- as.matrix(y)
                yy<-rbind(yy,y)
            }
        } else if(remaining > 0){
            inx <- seq_len(remaining)
            yy <- x[inx, , drop = FALSE]
            yy <- as.matrix(yy)
        }
        return (yy)
    }
#################################################################
#测试ref数据集过大而采取下采样的方式
# Created by: yanyuting
# Created on: 2022/5/20
    downsample<-function(object,n=5000){
        types<-table(object$cluster)
        Idents(object)<-object$cluster
        if(types[[1]]>=6000){
            all_obj<-subset(object,subset=cluster==names(types[1]),downsample = n)
        }else{
            all_obj<-subset(object,subset=cluster==names(types[1]))
        }
        for(i in (2:length(types))){
            if(types[[i]]>=6000){
                sub_obj<-subset(object,subset=cluster==names(types[i]),downsample = n)
            }else{
                sub_obj<-subset(object,subset=cluster==names(types[i]))
            }
            all_obj<-merge(all_obj,sub_obj)
        }
        return(all_obj)
    }
    tmp.cluster<-ref@meta.data[,ref_cluster]
    tmp.cluster<-as.vector(unlist(tmp.cluster))
    ref$cluster<-tmp.cluster
    ref<-downsample(ref,n=3500)
#################################################################
    ref_counts<-ref$RNA@counts
    ref_counts<-as.matrix(ref_counts)
    qur_counts<-qur$RNA@counts
    qur_counts<-as.matrix(qur_counts)
#分块读取matrix的时候部分行会重复！
    ref_genes<-rownames(ref_counts)
    qur_genes<-rownames(qur_counts)
    ref_upper<-toupper(ref_genes)
    qur_upper<-toupper(qur_genes)
    rownames(ref_counts)<-ref_upper
    rownames(qur_counts)<-qur_upper
    inter_genes<-intersect(x=ref_upper,y=qur_upper)
    print(length(inter_genes))
    ref_counts<-ref_counts[rownames(ref_counts)%in%inter_genes,]
    qur_counts<-qur_counts[rownames(qur_counts)%in%inter_genes,]
    count_list<-list(ref_counts,qur_counts)
    print(paste("未拆分的counts.list构建完毕",str(ref_counts),str(qur_counts)))
#切分待测的矩阵matrix
    splitMatrix<-function(qur_matrix,chunk=chunks){
        qur_matrix<-t(qur_matrix)
        n<-nrow(qur_matrix)
        cc<-list()
        if(n<=chunk){
            #cc<-append(cc,t(qur_matrix))
	    cc[[1]]<-t(qur_matrix)
            return(cc)
        }else{
            passes<-n%/%chunk
            remaining<-n%%chunk
            for(i in 1:passes){
                sub_matrix<-qur_matrix[c(((i-1)*chunk+1):(i*chunk)),,drop=FALSE]
                #cc<-append(cc,t(sub_matrix))
		cc[[i]]<-t(sub_matrix)
            }
            if(remaining>0){
                sub_matrix<-qur_matrix[c((passes*chunk+1):(passes*chunk+remaining)),,drop=FALSE]
                cc[[passes+1]]<-t(sub_matrix)
            }
            return(cc)
        }
    }
    qur_labels<-qur@meta.data[,qur_cluster]
    qur_labels<-as.data.frame(qur_labels)
    colnames(qur_labels)<-"type"
    print(paste("qur_labels是:",str(qur_labels)))
    qur_c<-splitMatrix(qur_counts)
    str(qur_c)
#在for循环外面对ref_counts进行操作，只需要执行一次Data1.csv
    ref_label<-ref@meta.data[,ref_cluster]
    ref_label<-as.data.frame(ref_label)
    colnames(ref_label)<-"type"
    print(paste("ref_label是:",str(ref_label)))
    if(length(inter_genes)>=2000){
        nf=2000
    }else{
        nf=nrow(data)
    }
    print("nf是：")
    print(nf)
    data<-count_list[[1]]
    label<-ref_label
    M <- nrow(data);
    new.label <- label[,1]
    pv1 <- sapply(1:M, function(i){
    mydataframe <- data.frame(y=as.numeric(data[i,]), ig=new.label)
    fit <- aov(y ~ ig, data=mydataframe)
    summary(fit)[[1]][["Pr(>F)"]][1]
    })
    names(pv1) <- rownames(data)
    pv1.sig <- names(pv1)[order(pv1)[1:nf]]
    sel.features <- unique(pv1.sig)
    print("sel.features:")
    print(str(sel.features))
    dir.create('input');
    #write.csv(t(count.list[[1]]),file='input/Data1.csv',quote=F,row.names=T)
    print("ref的Data1.csv创建完毕")
    saveRDS(qur_c,"qur_c.rds")
############################################################
    for(i in 1:length(qur_c)){
        qur_counts<-qur_c[[i]]
	count_list<-list(ref_counts,qur_counts)
        if(i!=length(qur_c)){
        	qur_label<-qur_labels[((i-1)*chunks+1):(i*chunks),]
		qur_label<-as.data.frame(qur_label)
		colnames(qur_label)<-"type"
        }else{
		#end_n<-nrow(qur_labels)
                qur_label<-qur_labels[((i-1)*chunks+1):(nrow(qur_labels)),] #nrow..的外侧一定要有()
		qur_label<-as.data.frame(qur_label)
		colnames(qur_label)<-"type"
        }
	print(paste("for循环",i,"中的qur_label:",str(qur_label)))      
        label_list<-list(ref_label,qur_label)
        str(label_list[[1]])
        str(label_list[[2]])
        print(Sys.time())
        print("--------------数据加载成功，list构建完毕---------------------")
        count.list<- list(count_list[[1]][sel.features,],count_list[[2]][sel.features,])
        print("save counts data to certain path: input")
	first_dir<-paste0('input/',as.character(i))
	dir.create(first_dir);
        write.csv(t(count.list[[2]]),file=paste0(first_dir,'/Data2.csv'),quote=F,row.names=T)
        print(Sys.time())
        print("---------------input创建完毕-----------------")
        label.list<-label_list
        new.dat1 <- count.list[[1]]; new.dat2 <- count.list[[2]]
        new.lab1 <- label.list[[1]]; new.lab2 <- label.list[[2]]
        #count.list[[1]]类似于一个matrix
        GenerateGraph <- function(Dat1,Dat2,Lab1,K,check.unknown){
            object1 <- CreateSeuratObject(counts=Dat1,project = "1",assay = "Data1",
                                  min.cells = 0,min.features = 0,
                                  names.field = 1,names.delim = "_")

            object2 <- CreateSeuratObject(counts=Dat2,project = "2",assay = "Data2",
                                  min.cells = 0,min.features =0,names.field = 1,
                                  names.delim = "_")
            print(Sys.time())
            print("--------------创建了seurat对象---------------")
            objects <- list(object1,object2)
            objects1 <- lapply(objects,function(obj){
                obj <- NormalizeData(obj,verbose=F)
                obj <- FindVariableFeatures(obj,
                                       selection.method = "vst",
                                       nfeatures = 2000,verbose=F)
                obj <- ScaleData(obj,features=rownames(obj),verbose=FALSE)
                obj <- RunPCA(obj, features=rownames(obj), verbose = FALSE)
                return(obj)
            })
    #'  Inter-data graph
            object.nn <- FindIntegrationAnchors(object.list = objects1,k.anchor=K,verbose=F)
            arc=object.nn@anchors
            d1.arc1=cbind(arc[arc[,4]==1,1],arc[arc[,4]==1,2],arc[arc[,4]==1,3])
            grp1=d1.arc1[d1.arc1[,3]>0,1:2]-1
            if (check.unknown){
                obj <- objects1[[2]]
                obj <- RunPCA(obj, features = VariableFeatures(object = obj),npcs=30,verbose=F)
                obj <- FindNeighbors(obj,verbose=F)
                obj <- FindClusters(obj, resolution = 0.5,verbose=F)
                hc <- Idents(obj); inter.graph=grp1+1
                scores <- metrics(lab1=Lab1,inter_graph=inter.graph,clusters=hc)
                saveRDS(scores,file=paste0(first_dir,'/statistical_scores.RDS'))
            }
            print(Sys.time())
            print("--------------创建了Inter-data graph---------------")
    #'  Intra-data graph
            d2.list <- list(objects1[[2]],objects1[[2]])
            d2.nn <- FindIntegrationAnchors(object.list =d2.list,k.anchor=K,verbose=F)
            d2.arc=d2.nn@anchors
            d2.arc1=cbind(d2.arc[d2.arc[,4]==1,1],d2.arc[d2.arc[,4]==1,2],d2.arc[d2.arc[,4]==1,3])
            d2.grp=d2.arc1[d2.arc1[,3]>0,1:2]-1
            final <- list(inteG=grp1,intraG=d2.grp)
            return (final)
        } 
        graphs <- suppressWarnings(GenerateGraph(Dat1=new.dat1,Dat2=new.dat2,
                                             Lab1=new.lab1,K=5,
                                             check.unknown=FALSE))
        write.csv(graphs[[1]],file=paste0(first_dir,'/inter_graph.csv'),quote=F,row.names=T)
        write.csv(graphs[[2]],file=paste0(first_dir,'/intra_graph.csv'),quote=F,row.names=T)
        write.csv(new.lab1,file=paste0(first_dir,'/label1.csv'),quote=F,row.names=F)
        write.csv(new.lab2,file=paste0(first_dir,'/label2.csv'),quote=F,row.names=F)
    }
    a1<-list.files('input')
    for(i in a1){
        write.csv(t(count.list[[1]]),file=paste0('input/',i,'/Data1.csv'),quote=F,row.names=T)
    }
}
#args[1]为参考物种的s4对象；args[2]为待测物种的s4对象
#!!!s4对象中物种细胞类型名必须为s4$cell_type
preprocess_data("mice.downsample3500.rds","mice.downsample3500.rds","cluster","cluster",chunk=10000)
