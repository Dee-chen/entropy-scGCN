# **********************************************************
# * Author        : YanYuting
# * Email         : yanyuting@genomics.cn
# * Create time   : 2022-05-24 16:57
# * Filename      : splitQurSeurat.r
# * Description   : 分割待测的rds文件，并保存在二级子目录下 
# **********************************************************
library(Seurat)
args<-commandArgs(T)
splitSeurat<-function(qur_rds,second_dir,chunk=20000){
	dir.create(second_dir)
	scRNA<-readRDS(qur_rds)
	col_num<-ncol(scRNA$RNA@counts)
	passes<-col_num%/%chunk
	remaining<-col_num%%chunk
	if(passes<1){
		save_dir<-paste0(second_dir,"/",qur_rds)
		saveRDS(scRNA,save_dir)
	}else{
		for(i in 1:passes){
			sub_data<-scRNA[,((i-1)*chunk+1):(i*chunk)]
			save_dir<-paste0(second_dir,"/",as.character(i),".rds")
			saveRDS(sub_data,save_dir)
		}
		if(remaining>0){
			sub_data<-scRNA[,(passes*chunk+1):(passes*chunk+remaining)]
			save_dir<-paste0(second_dir,"/",as.character(passes+1),".rds")
			saveRDS(sub_data,save_dir)
		}
	}
}
splitSeurat("./new_homo.RDS","test",20000)
