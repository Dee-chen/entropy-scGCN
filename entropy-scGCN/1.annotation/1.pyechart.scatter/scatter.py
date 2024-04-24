import pandas as pd
from pyecharts.charts import Scatter,Grid
from pyecharts import options as opts
def my_scatter(umap_file,cal_percent_file,title):  #title="tutrle.scatter"
    def my_sort(dict):
        new_dict=sorted(dict.items(),key=lambda d:d[1],reverse=True)
        return new_dict
    #all_data=pd.read_csv("./aa.table",sep="\t")
    all_data = pd.read_csv(umap_file, sep="\t")
    #percent_data=pd.read_csv("./cal_percent.tsv",sep="\t",index_col=0)
    percent_data=pd.read_csv(cal_percent_file,sep="\t",index_col=0)
    del percent_data["raw_class"]
    del percent_data["mosttype"]
    cell_type=list(set(all_data["cell_type"]))
    all_dimensions=percent_data.columns
    dimensions_list=[]
    x_data=[]
    y_data=[]
    labels=[]
    for i in cell_type:
        sub_data=all_data[all_data["cell_type"]==i]
        x_data.append(sub_data["UMAP_1"])
        y_data.append(sub_data["UMAP_2"])
        labels.append(i)
        sub_percent=percent_data[percent_data.index==i]
        cc = {}
        #cc=[]
        all_str=""
        for j in sub_percent.columns:
            cc[j] = sub_percent.loc[i, j]
        #print(cc)
        cc=my_sort(cc)
        #print(cc)
        for i in cc:
            tmp_str=str(i[0])+": "+str(round(i[1]*100,4))+"%</br>"
            all_str=all_str+tmp_str
            #cc.append(tmp_str)
        #cc = json.dumps(cc)
        #dd=[cc for i in sub_data["UMAP_1"]]
        #print(cc)
        #print(all_str)
        dimensions_list.append(all_str)
        #print(dimensions_list)
    c=Scatter(init_opts=opts.InitOpts(
        height="750px",
        width="100%"
    ))
    for i in range(len(x_data)):
        # c.add_dataset(
        #     dimensions=all_dimensions,
        #     source=dimensions_list[i]
        # )
        c.add_xaxis(xaxis_data=x_data[i])
        c.add_yaxis(labels[i],y_axis=y_data[i],symbol_size=3.0,
            label_opts=opts.LabelOpts(is_show=False),
            tooltip_opts=opts.TooltipOpts(
            formatter=dimensions_list[i],
                position=["50%","10%"]
            )) #不在点上刻画value
    c.set_series_opts()
    c.set_global_opts(
        xaxis_opts=opts.AxisOpts(
            type_="value",
            #splitline_opts=opts.SplitLineOpts(is_show=True) #坐标轴中的网格线
            axislabel_opts=opts.LabelOpts(rotate=50, interval=0),
            axisline_opts=opts.AxisLineOpts(
                is_on_zero=False
            )
        ),
        yaxis_opts=opts.AxisOpts(
            type_="value",
            axistick_opts=opts.AxisTickOpts(is_show=True),
            #splitline_opts=opts.SplitLineOpts(is_show=True)
            axisline_opts=opts.AxisLineOpts(
                is_on_zero=False
            )
        ),
        tooltip_opts=opts.TooltipOpts(
            is_show=True,
            # formatter="{a}:{b}"
        ),
        #title_opts=opts.TitleOpts(title="Turtle_seurat.clusters"),
	title_opts=opts.TitleOpts(title=title),
    )
    c.render("scatter.pyechart.html")
my_scatter("aa.table","cal_percent.tsv","Turtle_seurat.clusters")
