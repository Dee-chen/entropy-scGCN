source("splitQurSeurat.r")
source("preprocess_utilis.r")
args<-commandArgs(T)
splitSeurat("../new.turtle.rds","test",chunk=20000)
#args[1] refRDS args[2] qurRDS args[3] splitQurSeuratDir args[4] splitQurSeurat cell number
#args[5] cluster1 args[6] cluster2
#splitSeurat(args[2],args[3])
multipleCreateInput("test","../human.rds",cluster1="CellType",cluster2="Celltype")
#multipleCreateInput(args[3],args[1],args[4],args[5])
