#!/bin/bash
# **********************************************************
# * Author        : YanYuting
# * Email         : yanyuting@genomics.cn
# * Create time   : 2022-07-31 16:55
# * Filename      : metaNeighbor.sh
# * Description   : 
# **********************************************************
conda activate /jdfssz1/ST_SUPERCELLS/P20Z10200N0059/1.Macaca_shanghai/yanyuting/soft/miniconda3/envs/test
Rscript metaNeighbor.r final.mice.0621.rds turtle.cell.type.rds cell.type results/final.list F
