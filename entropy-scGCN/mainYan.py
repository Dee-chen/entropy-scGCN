import argparse
import time 
import subprocess
import sys
import os
import logging
def err_exit():
    sys.exit('\033[1;31;47m!!The program exited abnormally, please check the log file !!\033[0m')
def run_scGCN(ref_species, query_species,ref_seurat,query_seurat,input_unity,cluster_1,cluster_2,currentPath):
    #os.mkdir(os.sep.join([outpath,ref_species+"_"+query_species]))
    logging.info('Start processing type predictions[%s to %s].' % (ref_species, query_species))
    logging.info('Data pre-processing[%s to %s]...' % (ref_species, query_species))
    #os.chdir(os.sep.join([outpath,ref_species+"_"+query_species]))
    #stdout = subprocess.run( ['Rscript', os.sep.join([currentPath, 'preprocess.r']), ref_seurat,query_seurat,input_unity,cluster_1,cluster_2],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print('Rscript preprocess.r '+ref_seurat+" "+query_seurat+" "+input_unity+" "+cluster_1+" "+cluster_2)
    piple1=('Rscript preprocess.r '+ref_seurat+" "+query_seurat+" "+input_unity+" "+cluster_1+" "+cluster_2)
    os.system(piple1)
    #logging.debug('{} stdout:\n'.format('Data pre-processing') + stdout.stdout.decode('utf-8'))
    #logging.debug('{} stderr:\n'.format('Data pre-processing') + stdout.stderr.decode('utf-8'))
    logging.info("Data pre-processing finish[%s to %s]." % (ref_species, query_species))
    logging.info("Cell type predictions[%s to %s]..." % (ref_species, query_species))
    #stdout = subprocess.run(['python', os.sep.join([currentPath, 'train.py'])], stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    piple2=('python train.py')
    os.system(piple2)
    #logging.debug('{} stdout:\n'.format('Data pre-processing') + stdout.stdout.decode('utf-8'))
    #logging.debug('{} stderr:\n'.format('Data pre-processing') + stdout.stderr.decode('utf-8'))
    logging.info("Cell type predictions finish[%s to %s]." % (ref_species, query_species))
def sub_run(args):
    print(args)
    run_scGCN(args.r,args.q,args.s1,args.s2,args.i,args.c1,args.c2,args.p)
def main():
    pwd= os.getcwd()
    parser = argparse.ArgumentParser(description='scGCN pipline')
    parser.add_argument('-r', type=str, required=True, metavar='ref_specie',help='The name of reference specie')
    parser.add_argument('-q', type=str, required=True, metavar='qur_specie',help='The name of query specie')
    parser.add_argument('-s1', type=str, required=True, metavar='ref_seurat',help='The .rds file of ref specie')
    parser.add_argument('-s2', type=str, required=True, metavar='qur_seurat',help='The .rds file of qur specie')
    parser.add_argument('-i', type=str, required=True, metavar='input_unity',help='The directionary of the splited qur')
    parser.add_argument('-c1', type=str,required=True,metavar='cluster1',help='the type name of ref')
    parser.add_argument('-c2', type=str, required=True,metavar='cluster2',help='the type name of qur')
    parser.add_argument('-p', type=str, metavar='current_path',help='the current path')
    args = parser.parse_args()
    sub_run(args)
if __name__ == '__main__':
    main()
