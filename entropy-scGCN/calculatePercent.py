#!/user/bin/python
# **********************************************************
# * Author        : YanYuting
# * Email         : yanyuting@genomics.cn
# * Create time   : 2022-06-24 11:19
# * Filename      : test.py
# * Description   : 
# **********************************************************
import pandas as pd
import numpy as np
import openpyxl
import os
def createFinallist():
	path="./results/"
	files=os.listdir(path)
	allfiles=[]
	aa=range(1,len(files)+1)
	for i in aa:
		filespath="./results/"+str(i)+"/scGCN_data2_labels_2_prediction.csv"
		print("./results/"+str(i)+"/scGCN_data2_labels_2_prediction.csv")
		filein=pd.read_table(filespath,header=None,sep=" ")
		allfiles.append(filein)
	if(len(allfiles)==1):
		allfiles[0].columns=["raw_class","pre_class"]
		finalfiles=allfiles[0]
		finalfiles.to_csv("./results/final.list",index=False,header=None,sep="\t")
	elif(len(allfiles)==2):
		finalfiles=pd.concat([allfiles[0],allfiles[1]],axis=0)
		finalfiles.columns=["raw_class","pre_class"]
		finalfiles.to_csv("./results/final.list",index=False,header=None,sep="\t")
	else:
		finalfiles=pd.concat([allfiles[0],allfiles[1]],axis=0)
		for j in range(2,len(allfiles)):
			finalfiles=pd.concat([finalfiles,allfiles[j]],axis=0)
		finalfiles.columns=["raw_class","pre_class"]
		finalfiles.to_csv("./results/final.list",index=False,header=None,sep="\t")
	return finalfiles
def calculatePercent(df):
	data=df	
	#data=pd.read_csv("input.list",sep="\t",header=None)
	data.columns=["original","predicted"]
	aa=[1]*len(data.index)
	data["num"]=aa
	result = pd.pivot_table(data,index='predicted',columns='original',values='num',aggfunc = np.sum)
	result = result.fillna(0)
	num_result=result
	result = round((result / result.sum(axis=0)),4)
	result.index.name=None
	result.columns.name=None
	result=result.T
	num_result=num_result.T
	mosttype = result.idxmax(axis=1)
	mosttype=mosttype.tolist()
	result["mosttype"]=mosttype
	num_result["mosttype"]=mosttype
	result.insert(0,"raw_class",result.index.tolist())
	result.to_csv('cal_percent.tsv',sep="\t")
	num_result.insert(0,"raw_class",result.index.tolist())
	num_result.to_csv('num.cal_percent.tsv',sep="\t")
	lists=pd.DataFrame()
	lists["original"]=result.index.tolist()
	lists["predicted"]=mosttype
	print(lists)
	lists.to_csv("original2mosttype.list",sep="\t",index=False)
	#lists.to_csv("original2mosttype.list",sep="\t",header=None,index=False)
	left_obj=createFinallist()
	right_obj=result
	num_right_obj=num_result
	#raw_class=right_obj.index
	#right_obj.insert(0,"raw_class",raw_class)
	final_result=pd.merge(left_obj,right_obj,left_on="raw_class",right_on="raw_class",how="left")
	final_result.to_csv('final_cal_percent.tsv',index=False)
	num_final_result=pd.merge(left_obj,num_right_obj,left_on="raw_class",right_on="raw_class",how="left")
	num_final_result.to_csv('num_final_cal_percent.tsv',index=False)

#createFinallist()
	


