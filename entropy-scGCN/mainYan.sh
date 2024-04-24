#!/bin/bash
# **********************************************************
# * Author        : YanYuting
# * Email         : yanyuting@genomics.cn
# * Create time   : 2022-05-23 15:09
# * Filename      : preprocess.sh
# * Description   : 
# **********************************************************
conda activate /jdfssz1/ST_SUPERCELLS/P20Z10200N0059/1.Macaca_shanghai/yanyuting/soft/miniconda3/envs/test
python mainYan.py -r "mice" -q "mice" -s1 "./mice.downsample3500.rds" -s2 "mice.downsample3500.rds" -i "test" -c1 "cluster" -c2 "cluster" -p "./"
