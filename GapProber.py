#!/usr/bin/python3


import subprocess
import argparse
import os
import datetime, random



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='GapProber is a software to fill up gaps in genome assemblies.')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    # wrapper specific arguments
    parser.add_argument('-f', '--filled', required=False, type=str, help="<output file for filled scaffolds>", default=os.path.join(os.path.dirname(__file__), 'filled_contig.fastq'))
    #Reads can be given in either a set of fasta-formatted reads
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-r', '--reads', type=str, help="<fragment reads>", required=True)
    requiredNamed.add_argument('-s', '--scaffolds', type=str, help="<scaffolds with gaps>", required=True)
    
    #parser.print_help()
    args = vars(parser.parse_args())
    args['reads'] = args['reads'].split(',')
    wd = os.getcwd()
    wd_new = os.path.join(wd, 'trash_me_' + str(random.randrange(1000, 9999)))
    subprocess.check_call(['mkdir', wd_new])
    os.chdir(wd_new)
    log = ""
    
    if not args['scaffolds'][0] == '/': args['scaffolds'] = os.path.join(wd,args['scaffolds'])
    if not args['reads'][0][0] == '/':  args['reads'][0] = os.path.join(wd,args['reads'][0])
    if not args['reads'][1][0] == '/':  args['reads'][1] = os.path.join(wd,args['reads'][1])
    
    try:
        log = subprocess.check_call(['bowtie2-build', args['scaffolds'],wd_new+'/indexFile'])
        #print('bowtie2', '-x indexFile -1 '+args['reads'][0]+'  -2 '+args['reads'][1] +'  -S result.sam')
        log = subprocess.check_call(['bowtie2', '-x ',wd_new+'/indexFile', '-1', args['reads'][0],'-2 ', args['reads'][1],'  -S result.sam'])
        log = subprocess.check_call(['g++', '../main.cpp', '-o', 'main.out'])
        log = subprocess.check_call(['g++', '../figbird.cpp',  '-o', 'figbird.out'])
        log = subprocess.check_call(['mkdir', '-p', 'Gaps'])
        log = subprocess.check_call(['./main.out', args['scaffolds'], '250'])
        log = subprocess.check_call(['./figbird.out', args['scaffolds'], ' 0', '250'])
        log = subprocess.check_call(['python', '../reference.py' ,'filledContigs.fa', '../'+args['filled']])
        log = subprocess.check_call(['rm -rf', 'trash_me_*'])
    except:
        print(log)
        exit(1)
