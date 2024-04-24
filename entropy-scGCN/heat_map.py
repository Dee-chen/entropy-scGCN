import pandas as pd
# data1=pd.read_csv("./input/Data1.csv",index_col=0)
# data2=pd.read_csv("./input/Data2.csv",index_col=0)
# data3=pd.concat([data1,data2])
# print(data1.shape[0])
#
# print(data2.shape[0])
# print(data3.shape[0])
# dict={'a':1,'b':2,'c':3}
# print(dict["a"])
import os
def createCounts(data_dir,out_dir,type):
    #data_dir: ./hot_map/input/human
    if(type=="query"):
        dataa_dir=data_dir+"/Data2.csv"
        if out_dir=="macaca":
            label_dir=data_dir+"/scGCN_data2_labels_2_prediction.csv.txt"
        else:
            label_dir = data_dir + "/scGCN_data2_labels_2_prediction.csv1.txt"
        data2 = pd.read_csv(dataa_dir, index_col=0)
        label2 = pd.read_csv(label_dir, sep="\t", header=None)
        data = pd.DataFrame()
        cell_type = []
        data["cell_name"] = data2.index.tolist()
        data["cell_type"] = label2[1].tolist()
        cell_type = data["cell_type"].tolist()
        uniq = list(set(cell_type))
        cell_name = []
        for i in range(0, len(uniq)):
            dataa = data[data["cell_type"] == uniq[i]]
            cell_name.append(dataa)
        dict = {}
        for j in range(0, len(cell_name)):
            key = cell_name[j]["cell_type"].tolist()[0]
            names = cell_name[j]["cell_name"].tolist()
            count = data2[data2.index.isin(names)]
            dict[key] = count
    else:
        dataa_dir = data_dir + "/Data1.csv"
        label_dir = data_dir + "/Label1.csv"
        data2 = pd.read_csv(dataa_dir, index_col=0)
        label2 = pd.read_csv(label_dir, delimiter=" ", header=0)
        data = pd.DataFrame()
        cell_type = []
        data["cell_name"] = data2.index.tolist()
        data["cell_type"] = label2["type"].tolist()
        cell_type = data["cell_type"].tolist()
        uniq = list(set(cell_type))
        cell_name = []
        for i in range(0, len(uniq)):
            dataa = data[data["cell_type"] == uniq[i]]
            cell_name.append(dataa)
        dict = {}
        for j in range(0, len(cell_name)):
            key = cell_name[j]["cell_type"].tolist()[0]
            names = cell_name[j]["cell_name"].tolist()
            count = data2[data2.index.isin(names)]
            dict[key] = count
    print(out_file)
    for key, value in dict.items():
        print(key)
        tmp_count=dict[key]
        # tmp_count["cell_type"]=[]
        lens=len(dict[key].index)
        tmp_types=[]
        for n in range(0,lens):
            tmp_types.append(key)
        print(tmp_types)
        tmp_count["cell_type"]=tmp_types
        if key.find("/"):
            keyy=key.split("/")[0]
            file = "./hot_map/output/" + out_dir + "/" + keyy + ".csv"
            tmp_count.to_csv(file)
        else:
            file = "./hot_map/output/"+out_dir+"/"+key + ".csv"
            tmp_count.to_csv(file)
path = "./hot_map/input"
files= os.listdir(path)
txts = []
for file in files:
    position = path+'/'+ file
    out_file=str(file)
    if out_file=="human":
        createCounts(position,out_file,"reference")
    else:
        createCounts(position,out_file,"query")


