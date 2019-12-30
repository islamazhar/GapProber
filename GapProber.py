#!/usr/bin/python3


from subprocess import call
import argparse
import os
import datetime, random



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='GapProber 1.0')

    # wrapper specific arguments
    parser.add_argument('-f', '--filled', required=False, type=str, help="<output file for filled scaffolds>", default=os.path.join(os.path.dirname(__file__), 'filled_contig.fastq'))
    #Reads can be given in either a set of fasta-formatted reads
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-r', '--reads', type=str, help="< fragment reads >", required=True)
    requiredNamed.add_argument('-s', '--scaffolds', type=str, help="< scaffolds with gaps >", required=True)
    #parser.print_help()
    args = vars(parser.parse_args())
    args['reads'] = args['reads'].split(',')
    call('mkdir -p Gaps')
    wd = os.getcwd()
    wd_new = os.path.join(wd, 'tmp.' + str(random.randrange(1000, 9999)))
    subprocess.check_call(['mkdir', wd_new])
    os.chdir(wd_new)
    try:
        log = subprocess.check_call('bowtie2-build '+args['scaffolds']+' indexFile')
        log = subprocess.check_call('bowtie2 -x indexFile -1 '+ args['reads'][0] +' -2 '+args['reads'][1] +' $read_dir2  -S result.sam')
        log = subprocess.check_call('g++ ./main.cpp  -o main.out')
        log = subprocess.check_call('g++ ./figbird.cpp  -o figbird.out')
        #log = subprocess.check_call('cp '+wd+'/builds/main.out '+wd_new)
        #log = subprocess.check_call('cp '+wd+'/builds/figbird.out '+wd_new)
        log = subprocess.check_call('./main.out '+args['scaffolds'])
        log = subprocess.check_call('./figbird.out '+args['scaffolds']+' 0')
        log = subprocess.check_call('mv filledContigs.fa ./'+args['filled'])
    except:
        print(log)
        exit(1)
