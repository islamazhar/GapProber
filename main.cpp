//
//  main.cpp
//  gap
//
//

#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <cstring>

using namespace std;
#include "cgal.h"
#include "mlga.h"
#include <math.h>

#define MAX_REC_LEN 1024
#define MAX_READLENGTH 200
#define MAX_NAMELENGTH 200
#define HASH_TABLE_SIZE 10001

//===================Defined Structs========================

struct SAM
{
	char qname[MAX_NAMELENGTH];
	int flag;
	char rname[MAX_NAMELENGTH];
	int pos;
	int mapq;
	char cigar[MAX_READLENGTH];
	char rnext[MAX_NAMELENGTH];
	int pnext;
	int tlen;
	char seq[MAX_NAMELENGTH];
	char qual[MAX_NAMELENGTH];
	char md[MAX_NAMELENGTH];
	unsigned long ih;
	int nm;
	long int contigNo;
};

struct Gap
{
    int contigNo;
    long int gapStart;
    int gapLength;
};


struct Contig
{
	char contigName[1000];
	long int contigNo;
};

//==================================================

//==================File variables==============
FILE *outFile;
FILE *linkingFile;
FILE *singletonFile;
FILE *mapFile;
FILE *contigFile;
FILE *unFile;
FILE *statFile;
FILE *scaffoldFile;
FILE **gapFiles;
FILE *gapInfoFile;

//=================Int variables================
long int noGaps=0;
int **readCounts;
long int unCount=0;
long int totalCount=0;
unsigned long maxReadLength=0;
long int MAX_FRAGMENT_SIZE=5000;
long int noContigs;
int maxDistance=200;


//=================Double vars===================

double *contigReadCounts;
double insertSizeMean=0;
double insertSizeVar=0;
double squaredError=0;

//==================All vectors==================
vector <SAM *> reads1;
vector <SAM *> reads2;

vector <SAM *> mixedReads1;
vector <SAM *> mixedReads2;

vector<char*> contigs;
vector<char*> contigNames;
vector<unsigned long> contigLengths;

vector<Contig *> hashTable[HASH_TABLE_SIZE];
vector<Gap *> gaps;

//================Char arrays====================

char line1[MAX_REC_LEN];
char line2[MAX_REC_LEN];


//=================Declaration Done================================


void reverse(char *reverse, char *read)
{
    
	char ch='A';
	unsigned long readLength=strlen(read);
	for(unsigned long i=1; i<=readLength;i++)
	{
		ch=read[i-1];
		if(ch=='A')
			reverse[readLength-i]='T';
		else if(ch=='C')
			reverse[readLength-i]='G';
		else if(ch=='G')
			reverse[readLength-i]='C';
		else if(ch=='T')
			reverse[readLength-i]='A';
		else
			reverse[readLength-i]='N';
	}
	reverse[readLength]='\0';
    
}

unsigned long getHash(char *str)
{
    unsigned long hash = 5381;
    int c;
    
    while (c = *str++)
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
    
    return hash;
}

void initHashTable()
{
    
	unsigned long index;
	for(long int i=0;i<noContigs;i++)
	{
		index=getHash(contigNames[i]) % HASH_TABLE_SIZE;
		Contig *c=new Contig();
		strcpy(c->contigName,contigNames[i]);
		c->contigNo=i;
		hashTable[index].push_back(c);
	}
}

int getContigNo(char *contigName)
{
	unsigned long index=getHash(contigName) % HASH_TABLE_SIZE;
    
	for(int i=0;i<hashTable[index].size();i++)
	{
		if(strcmp(hashTable[index][i]->contigName,contigName)==0)
			return hashTable[index][i]->contigNo;
	}
	return -1;
    
}


void writeSam(SAM* read, FILE *out)
{
    /*	fprintf(out,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t\%s\tIH:i:%d\n",read->qname,read->flag,read->rname,read->pos,
     read->mapq,read->cigar,read->rnext,read->pnext,read->tlen,read->seq,read->qual,read->md,read->ih);
     */
	fprintf(out,"%s\t%d\t%ld\t%d\t%s\t%d\t%s\t\%s\tIH:i:%ld\n",read->qname,read->flag,read->contigNo,read->pos,
            read->cigar,read->tlen,read->seq,read->md,read->ih);
    
}

void printSam(SAM* read)
{
	printf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t\%s\tIH:i:%ld\n",read->qname,read->flag,read->rname,read->pos,
           read->mapq,read->cigar,read->rnext,read->pnext,read->tlen,read->seq,read->qual,read->md,read->ih);
    
    
    /*	printf("%s\t%d\t%d\t%s\t%d\t%s\t\%s\tIH:i:%d\n",read->qname,read->flag,read->pos,
     read->cigar,read->tlen,read->seq,read->md,read->ih);
     */
}


int checkPos(int contigNo, long int pos, int strandNo)
{
    for(int i=0;i<gaps.size();i++)
    {
        if(contigNo==gaps[i]->contigNo && ((strandNo==0 && pos>(gaps[i]->gapStart-maxDistance) && pos<(gaps[i]->gapStart+gaps[i]->gapLength)) || (strandNo==1 && pos>(gaps[i]->gapStart+gaps[i]->gapLength) && pos<(gaps[i]->gapStart+gaps[i]->gapLength+maxDistance))))
        {
            return i;
            
        }
        
        //same two conditions
        //    if(contigNo2==0 && ((strandNo2==0 && pos2>(1000000-200) && pos2<(1000000+200)) || (strandNo2==1 && pos2>(1000000+200) && pos2<(1000000+400))))
        
        if(contigNo==gaps[i]->contigNo && ((strandNo==0 && pos>(gaps[i]->gapStart-maxDistance) && pos<(gaps[i]->gapStart+gaps[i]->gapLength)) || (strandNo==1 && pos>(gaps[i]->gapStart+gaps[i]->gapLength) && pos<(gaps[i]->gapStart+gaps[i]->gapLength+maxDistance))))
        {
            return i;
        }    
    }
        
    return -1;
    
    
}

void printVectors(FILE * out, FILE * u)
{
    
	SAM *read1, *read2;
	unsigned long readLength1,readLength2;
	int pos1, pos2;
    //	int insertSize;
    
	char temp[MAX_READLENGTH];
    
	unsigned long ih=reads1.size();
    
	int ih1_0=0;
	int ih2_0=0;
    
    
	for(int i=0;i<ih;i++)
	{
		if(strcmp(reads1[i]->rname,"*")!=0 && reads1[i]->nm==0)
		{
			ih1_0++;
		}
	}
    
	for(int i=0;i<ih;i++)
	{
		if(strcmp(reads2[i]->rname,"*")!=0  && reads2[i]->nm==0)
		{
			ih2_0++;
		}
	}
	if(ih1_0>0)
	{
        for(int i=0;i<ih;i++)
        {
            if(strcmp(reads1[i]->rname,"*")!=0  && reads1[i]->nm==0)
            {
                contigReadCounts[reads1[i]->contigNo]+=1/(double)ih1_0;
            }
        }
	}
    
	if(ih2_0>0)
	{
        for(int i=0;i<ih;i++)
        {
            if(strcmp(reads2[i]->rname,"*")!=0  && reads2[i]->nm==0)
            {
                contigReadCounts[reads2[i]->contigNo]+=1/(double)ih2_0;
            }
        }
	}
    
	for(unsigned long i=0;i<ih;i++)
	{
		read1=reads1[i];
		read2=reads2[i];

		if(strcmp(read1->rname,"*")==0 || strcmp(read2->rname,"*")==0)
		{
		    //Never here
		    //		int i=0;
            //		int j=0;
			int ncount1=0;
			int ncount2=0;
            
			readLength1=strlen(read1->seq);
			if(readLength1>maxReadLength)
				maxReadLength=readLength1;
            
			readLength2=strlen(read2->seq);
			if(readLength2>maxReadLength)
				maxReadLength=readLength2;
            
            
			for(int i=0;i<readLength1;i++)
			{
				if(read1->seq[i]=='N'||read1->seq[i]=='n')
				{
					ncount1++;
				}
			}
            
			for(int i=0;i<readLength2;i++)
			{
                
				if(read2->seq[i]=='N'||read2->seq[i]=='n')
				{
					ncount2++;
				}
			}
            
			if(ncount1/(double)readLength1 < 0.8 && ncount2/(double)readLength2 < 0.8)
			{
			    //cout<<"PrintVectors*************"<<ncount1<<" $$$$"<<ncount2<<endl;
				unCount++;
				totalCount++;
                
                
				fputs("@",u);
				fputs(read1->qname,u);
				fputs("\n",u);
                
                
				int strandNo=(read1->flag&16)>>4;
				if(strandNo==1)
				{
					reverse(temp,read1->seq);
					fputs(temp,u);
				}
				else
				{
					fputs(read1->seq,u);
				}
				fputs("\n",u);
                
				fputs("+",u);
				fputs(read1->qname,u);
				fputs("\n",u);
                
				fputs(read1->qual,u);
				fputs("\n",u);
                
				fputs("@",u);
				fputs(read2->qname,u);
				fputs("\n",u);
                
				strandNo=(read2->flag&16)>>4;
				if(strandNo==1)
				{
					reverse(temp,read2->seq);
					fputs(temp,u);
				}
				else
				{
					fputs(read2->seq,u);
                    
				}
				fputs("\n",u);
                
				fputs("+",u);
				fputs(read2->qname,u);
				fputs("\n",u);
                
				fputs(read2->qual,u);
				fputs("\n",u);
				
				for(int i=0;i<reads1.size();i++)
					delete reads1[i];
                
				for(int i=0;i<reads2.size();i++)
					delete reads2[i];
                
                
				reads1.clear();
				reads2.clear();
				
				return;
			}
            
		}
		else if(strcmp(read1->rname,read2->rname)!=0)
		{
            //Never here
			readCounts[read1->contigNo][read2->contigNo]++;
			printSam(read1);
			printSam(read2);
			cout<<"Warning: contig names different"<<endl;
            
            //		getchar();
            
		}
		else
		{
			readCounts[read1->contigNo][read2->contigNo]++;
            
			pos1=read1->pos;
			readLength1=strlen(read1->seq);
			if(readLength1>maxReadLength)
				maxReadLength=readLength1;
			read1->ih=ih;
            
            
			pos2=read2->pos;
			readLength2=strlen(read2->seq);
			if(readLength2>maxReadLength)
				maxReadLength=readLength2;
			read2->ih=ih;
            
			writeSam(read1,out);
			writeSam(read2,out);
            //cout<<"In else--->>>>>>>>>>>>>>>>>>>. \n";
		}
	}
	totalCount++;
    
	
    
	for(int i=0;i<reads1.size();i++)
		delete reads1[i];
    
	for(int i=0;i<reads2.size();i++)
		delete reads2[i];
	
    
	reads1.clear();
	reads2.clear();
    
}


void printMixedVectors(FILE * linking, FILE * single, FILE * u)
{
	SAM *read1, *read2;
	unsigned long readLength1,readLength2;
    long int pos1, pos2;
    int contigNo1, contigNo2;
    int strandNo1, strandNo2;
    
    int gapFileNo;
    //	int insertSize;
    
    //	char temp[MAX_READLENGTH];
    
    
    char temp[MAX_READLENGTH];
    
	int ih1=mixedReads1.size();
    int ih2=mixedReads2.size();
    
	int ih1_0=0;
	int ih2_0=0;
    
    
	for(int i=0;i<ih1;i++)
	{
		if(strcmp(mixedReads1[i]->rname,"*")!=0 && mixedReads1[i]->nm==0)
		{
			ih1_0++;
		}
	}
    
	for(int i=0;i<ih2;i++)
	{
		if(strcmp(mixedReads2[i]->rname,"*")!=0  && mixedReads2[i]->nm==0)
		{
			ih2_0++;
		}
	}
	
	if(ih1_0>0)
	{
        for(int i=0;i<ih1;i++)
        {
            if(strcmp(mixedReads1[i]->rname,"*")!=0 && mixedReads1[i]->nm==0)//Unnecessary??? Done before at line 459
            {
                contigReadCounts[mixedReads1[i]->contigNo]+=1/(double)ih1_0;
            }
        }
	}
	
	if(ih2_0>0)
	{
        for(int i=0;i<ih2;i++)
        {
            if(strcmp(mixedReads2[i]->rname,"*")!=0 && mixedReads2[i]->nm==0)
            {
                contigReadCounts[mixedReads2[i]->contigNo]+=1/(double)ih2_0;
            }
        }
	}
	
	//cout<<"ih1 = "<<ih1<<"\tih2 = "<<ih2<<endl;//They are always both 1, everytime called upon
	
    if(ih1>1 || ih2>1)
    {
        for(int i=0;i<mixedReads1.size();i++)
        {
            read1=mixedReads1[i];
        }
        for(int j=0;j<mixedReads2.size();j++)
        {
            read1=mixedReads2[j];
        }
    }
    
    //read1 and read2 are two SAM records
    //cout<<mixedReads1.size()<<" **** "<<mixedReads2.size()<<endl;
	for(int i=0;i<mixedReads1.size();i++)
	{
		read1=mixedReads1[i];
        
		for(int j=0;j<mixedReads2.size();j++)
		{
			read2=mixedReads2[j];
            if(i==0 && j==0)//Why special case?
            {
                readLength1=strlen(read1->seq);
                if(readLength1>maxReadLength)
                    maxReadLength=readLength1;
                
                readLength2=strlen(read2->seq);
                if(readLength2>maxReadLength)
                    maxReadLength=readLength2;
                
                
                int ncount1=0;
                int ncount2=0;
                
                
                for(int i=0;i<readLength1;i++)
                {
                    if(read1->seq[i]=='N'||read1->seq[i]=='n')
                    {
                        ncount1++;
                    }
                }
                
                for(int i=0;i<readLength2;i++)
                {
                    
                    if(read2->seq[i]=='N'||read2->seq[i]=='n')
                    {
                        ncount2++;
                    }
                }

                //IF both the mates in the read pair contains less than 80% Ns, we put them in unmapped.txt??why???
                //cout<<"PrintMixedVectors=>ncount1 = "<<ncount1<<" ncount2 = "<<ncount2<<endl;
                if(ncount1/(double)readLength1 < 0.8 && ncount2/(double)readLength2 < 0.8)
                {
                    //cout<<"PrintMixedVectors###########"<<ncount1<<" $$$$"<<ncount2<<endl;
                    unCount++;
                    totalCount++;
                    fputs("@",u);
                    fputs(read1->qname,u);
                    fputs("\n",u);
                    
                    
                    int strandNo=(read1->flag&16)>>4;
                    if(strandNo==1)// If the seq is reverse complemented, then reverse it again
                    {
                        reverse(temp,read1->seq);
                        fputs(temp,u);
                    }
                    else
                    {
                        fputs(read1->seq,u);
                    }
                    fputs("\n",u);
                    
                    fputs("+",u);
                    fputs(read1->qname,u);
                    fputs("\n",u);
                    
                    fputs(read1->qual,u);
                    fputs("\n",u);
                    
                    fputs("@",u);
                    fputs(read2->qname,u);
                    fputs("\n",u);
                    
                    strandNo=(read2->flag&16)>>4;
                    if(strandNo==1)
                    {
                        reverse(temp,read2->seq);
                        fputs(temp,u);
                    }
                    else
                    {
                        fputs(read2->seq,u);
                        
                    }
                    fputs("\n",u);
                    
                    fputs("+",u);
                    fputs(read2->qname,u);
                    fputs("\n",u);
                    
                    fputs(read2->qual,u);
                    fputs("\n",u);
                }
                else
                {
                    //Otherwise they are cleared from the mixed reads
                    //Never in this block
                    for(int i=0;i<mixedReads1.size();i++)
                        delete mixedReads1[i];
                    
                    for(int i=0;i<mixedReads2.size();i++)
                        delete mixedReads2[i];
                    
                    
                    mixedReads1.clear();
                    mixedReads2.clear();
                    
                    return;
                    
                }
            }
            
            
            //=====================For rest of the i and j combinations including 0 and 0==========================================
            //4 cases
            
            //Case1: If both are unmapped,just clear,Do nothing
            if(((read1->flag) & 4) != 0 && ((read2->flag)& 4) != 0)
			{
                for(int i=0;i<mixedReads1.size();i++)
                    delete mixedReads1[i];
                for(int i=0;i<mixedReads2.size();i++)
                    delete mixedReads2[i];
                
                //cout<<"Case -1 ####################\n";
				mixedReads1.clear();
				mixedReads2.clear();
				
				return;
			}
            //Case2: If one of them read 1 is not unmapped and other one(read 2) is unmapped, write in singleton file
            //Because read2 is unmapped, find which gap it falls in
            else if(((read1->flag) & 4) == 0 && ((read2->flag)& 4) != 0)
            {
            //cout<<"Case -2 ####################\n";
                for(int i=0;i<mixedReads1.size();i++)
                {
                    writeSam(mixedReads1[i], single);
                }
                for(int i=0;i<mixedReads2.size();i++)
                {
                    writeSam(mixedReads2[i], single);
                }
                
                for(int i=0;i<mixedReads1.size();i++)
                {
                    contigNo1=mixedReads1[i]->contigNo;
                    pos1=mixedReads1[i]->pos;
                    strandNo1=(mixedReads1[i]->flag&16)>>4;
                    
             //       if(contigNo1==0 && ((strandNo1==0 && pos1>(1000000-200) && pos1<(1000000+200)) || (strandNo1==1 && pos1>(1000000+200) && pos1<(1000000+400))))
                    
                    gapFileNo=checkPos(contigNo1,pos1, strandNo1);
                    if(gapFileNo>=0)
                    {
                    //    gapFileNo=0;
                        writeSam(mixedReads1[i], gapFiles[gapFileNo]);
                        writeSam(mixedReads2[0], gapFiles[gapFileNo]);
                        
                    }
                }
                
                for(int i=0;i<mixedReads1.size();i++)
                    delete mixedReads1[i];
                for(int i=0;i<mixedReads2.size();i++)
                    delete mixedReads2[i];
                
                
                mixedReads1.clear();
                mixedReads2.clear();
                
                return;
                
            }
            //Case3: If one of them read 2 is not unmapped and other one(read 1) is unmapped, write in singleton file
            //Because read2 is unmapped, find which gap it falls in
            else if(((read1->flag) & 4) != 0 && ((read2->flag)& 4) == 0)
            {
            //cout<<"Case -3 ####################\n";
                for(int i=0;i<mixedReads1.size();i++)
                {
                    writeSam(mixedReads1[i], single);
                }
                for(int i=0;i<mixedReads2.size();i++)
                {
                    writeSam(mixedReads2[i], single);
                }
                
                for(int i=0;i<mixedReads2.size();i++)
                {
                    contigNo2=mixedReads2[i]->contigNo;
                    pos2=mixedReads2[i]->pos;
                    strandNo2=(mixedReads2[i]->flag&16)>>4;
                    
                //    if(contigNo2==0 && ((strandNo2==0 && pos2>(1000000-200) && pos2<(1000000+200)) || (strandNo2==1 && pos2>(1000000+200) && pos2<(1000000+400))))
                    gapFileNo=checkPos(contigNo2,pos2, strandNo2);
                    if(gapFileNo>=0)
                    {
                    //    gapFileNo=0;
                        writeSam(mixedReads2[i], gapFiles[gapFileNo]);
                        writeSam(mixedReads1[0], gapFiles[gapFileNo]);
                        
                    }
                    
                
                }
                
                
                for(int i=0;i<mixedReads1.size();i++)
                    delete mixedReads1[i];
                for(int i=0;i<mixedReads2.size();i++)
                    delete mixedReads2[i];
                
                mixedReads1.clear();
                mixedReads2.clear();
                return;
            }
            //case4: Those that are both mapped but from different references?
            else if(strcmp(read1->rname,read2->rname)!=0)
			{
			    //cout<<"Case -4 ####################\n";
                read1->ih=ih1;
                read2->ih=ih2;
                readCounts[read1->contigNo][read2->contigNo]++;
                writeSam(read1,linking);
                writeSam(read2, linking);
                //	printSam(read1);
                //	printSam(read2);
                //	getchar();
                
			}
		}
	}
    
	for(int i=0;i<mixedReads1.size();i++)
		delete mixedReads1[i];
	for(int i=0;i<mixedReads2.size();i++)
		delete mixedReads2[i];
	
    
	mixedReads1.clear();
	mixedReads2.clear();
    
}

SAM *getSAM(char *line)
{
	SAM *sam=new SAM;
	char *temp;
    
	temp=strtok(line,"\t\n ");
	strcpy(sam->qname,temp);
    
	temp=strtok(NULL,"\t\n ");
	sam->flag=atoi(temp);
    
	temp=strtok(NULL,"\t\n ");
	strcpy(sam->rname,temp);
    
	temp=strtok(NULL,"\t\n ");
	sam->pos=atoi(temp);
    
	temp=strtok(NULL,"\t\n ");
	sam->mapq=atoi(temp);
    
	temp=strtok(NULL,"\t\n ");
	strcpy(sam->cigar,temp);
	
	temp=strtok(NULL,"\t\n ");
	strcpy(sam->rnext,temp);
    
	temp=strtok(NULL,"\t\n ");
	sam->pnext=atoi(temp);
    
	temp=strtok(NULL,"\t\n ");
	sam->tlen=atoi(temp);
    
	temp=strtok(NULL,"\t\n ");
	strcpy(sam->seq,temp);
	
	temp=strtok(NULL,"\t\n ");
	strcpy(sam->qual,temp);
    
	sam->nm=-1;
    
	while((temp=strtok(NULL,"\t\n "))!=NULL)
	{
		if(temp[0]=='M' && temp[1]=='D')
		{
			strcpy(sam->md,temp);
		}
		if(temp[0]=='N' && temp[1]=='M')
		{
			sam->nm=atoi(&temp[5]);
		}
	}
    
	sam->contigNo=getContigNo(sam->rname);
    
	return sam;
}

void printHelp()
{
    
	cout<<"Invalid parameters"<<endl;
	exit(1);
    
}

int main(int argc, char *argv[])
{
	int MAX_FILE_READ = 0;
	noContigs=0;
	long int contigLength=0;
	unsigned long read;
	long int bufferLength=1024;
    long int tempContigLength=0;
    int nStart=0;//How many gaps are there?
    long int nStartPos=0;//What are those gaps start positions?
    int nCount=0,afterCount=0;//ncount is number of NNNN in a particular gap
	
    char *line= new char[MAX_REC_LEN];
    //char *templine= new char[MAX_REC_LEN];
	const char * mapFileName="result.sam";
	char *	contigFileName=argv[1];
    char *contig;
	char *newcontig;
	char *contigName;
	
	
	contig=new char[bufferLength];
	contig[0]='\0';
	
	
	MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);
	
	contigFile=fopen(contigFileName, "r");
    
    if (contigFile == NULL)
	{
		printf("Can't open contig file\n");
		exit(1);
	}
    
	
    //======================================Start reading the contig file=======================================
    
    //It contains 19 scaffolds in total with gaps
    int cnt=0;
	while(fgets(line, MAX_FILE_READ, contigFile)!=NULL)
	{
	    cnt++;
		if(line[0]==';')//There is no ';' in contig file, why this?
		{
			continue;
		}
		else if(line[0]=='>')//fasta format name of the contig
		{
			contigName = new char[strlen(line)];
			strcpy(contigName,line+1);
			contigName[strlen(contigName)-1]='\0';
			contigNames.push_back(strtok(contigName," \t\n"));//In case of multiple contigs, 
			                                                    //we are pushing each name into a vector of char array
			
			if(contigLength>0)//This block is accessed after completing the full reading of each contig, not at the start
			{
				noContigs++;
				contig=(char *)realloc(contig, (contigLength+1)*sizeof(char));
				contigs.push_back(contig);
				contigLengths.push_back(contigLength);
				contigLength=0;
				bufferLength=1024;
				contig=new char[bufferLength];
				contig[0]='\0';
			}
		}
		else
		{
			read=strlen(line);
			tempContigLength=contigLength;
			if(read<MAX_FILE_READ-1)
			{
				contigLength+=(read-1);
			}
			else
			{
				contigLength+=MAX_FILE_READ-1;
				read++;
				
			}
			if(contigLength>bufferLength)
			{
				bufferLength=max(bufferLength*2,contigLength+1);
				newcontig=new char[bufferLength];
				strcpy(newcontig,contig);
				line[read-1]='\0';
				strcpy(newcontig+tempContigLength, line);
				delete []contig;
				contig=newcontig;
			}
			else
			{
				line[read-1]='\0';
				strcpy(contig+tempContigLength, line);
			}
            
		}
		//if(cnt == 3)break;
		
	}
    
	noContigs++;
	contig=(char *)realloc(contig, (contigLength+1)*sizeof(char));
	contigs.push_back(contig);
	contigLengths.push_back(contigLength);
	
    
	//cout<<"Total no. of contigs = "<<noContigs<<endl;
    fclose(contigFile);
    
    //=================================End of reading the contig file==========================================//
    
    
    //=================================Gather Info about gaps in the contig====================================//
    
    gapInfoFile=fopen("gapInfo.txt","w");
    
    for(int i=0;i<noContigs;i++)//For each contig
    {
        for(long int j=0;j<contigLengths[i];j++)//For the entire contig length
        {
            if(contigs[i][j]=='N' || contigs[i][j]=='n')//If there is a gap at that position
            {
                if(nStart==0)
                {
                    nStart=1;
                    nCount=1;
                    nStartPos=j;
                }
                else
                {
                    nCount++;
                }
            }
            else if(nStart==1)//if there was ever any gap
            {
                if(nCount>=1)//Why >=10? if the gaplen was >=10, ony then consider the gap
                {
                    Gap * g=new Gap();
                    g->contigNo=i;
                    g->gapStart=nStartPos;
                    g->gapLength=nCount;
                    gaps.push_back(g);
                    
                    fprintf(gapInfoFile,"%d\t%ld\t%d\n",g->contigNo,g->gapStart,g->gapLength);
                }
                nStart=0;
            }
            
        }
        
    }
    
    //===============================================Done with gaps=================================//
    
    initHashTable();//Fill a hashtable of contigs for some future purpose[Insignificant]? hashtable has only name and number of the contigs
    
	contigReadCounts=new double[noContigs];//1-D array of doubles
    
	for(int i=0;i<noContigs;i++)
	{
		contigReadCounts[i]=0;
	}

	statFile=fopen("stat.txt","w");
    scaffoldFile=fopen("scaffold.txt","w");
    
    
	readCounts=new int*[noContigs];//2D array of ints, row->contigs, col->contigs? 19*19
	
	for (int i = 0; i < noContigs; i++)
		readCounts[i] = new int[noContigs];
    
    
    
	for(int i=0;i<noContigs;i++)
	{
		for(int j=0;j<noContigs;j++)
		{
			readCounts[i][j]=0;
		}
	}
    
	//=======Start working with the SAM file output by Bowtie2====================
	
	//=================================Preprocessing works========================
	
	mapFile=fopen(mapFileName, "r");
    
	if (mapFile == NULL)
	{
		printf("Can't open map file\n");
		exit(1);
	}
    
	outFile=fopen("myout.sam","w");
    linkingFile=fopen("linking.sam","w");
    singletonFile=fopen("singletons.sam","w");
   
    
    gapFiles=new FILE*[gaps.size()];
    
    char *gapFileName=new char[200];
    char *temp=new char[200];
    
    for(int i=0;i<gaps.size();i++)
	{
        strcpy(gapFileName,"Gaps/gaps_");
        myitoa(i,temp,10);
        strcat(gapFileName,temp);
        strcat(gapFileName,".sam");
        gapFiles[i]=fopen(gapFileName,"w");

    }
    
	unFile=fopen("unmapped.txt","w");
    
    //	char *temp,nhstring[200];
    
    
	int it=0; 
    int end=0;
    
	char preqname1[100];
	char preqname2[100];
	
	strcpy(preqname1,"*");
	strcpy(preqname2,"*");
    
    char qname[200];
    int flag_segment;
    
    //===============================Starting the read of SAM file result.sam line by line================
        
	while(fgets(line, MAX_FILE_READ, mapFile)!=NULL)
	{
		if(line[0]=='@')//Skip headers
			continue;
        
		it++;
        

		
		SAM *read1=getSAM(line);

        //
        //The following condition means that
        //this read is "NOT" part of a pair that aligned in a pair end  fashion
        //So take extra care for them to find a mate that maps??!!
        while((read1->flag & 2) == 0)
        {

            strcpy(qname,read1->qname);
            flag_segment=read1->flag & 192;//Only interested in MSBth and (MSB-1)th position of the flag according to the AND op.
            //192 = 1100 0000
            //MSB = Indicates if this read is mate 2 in the pair
            //MSB-1 = Indicates if this read is mate 1 in the pair
            
            //Paired end reads are arranged one by one in SAM file
            //their qnames are same always
            while(strcmp(qname,read1->qname)==0 && (read1->flag & 192) == flag_segment)
            //Continue reading the entire file and search if
            {
                mixedReads1.push_back(read1);//
                if(fgets(line, MAX_FILE_READ, mapFile)!=NULL)
                {
                    //Enters here only once
                    //printf("Here %d?\n",it);
                    read1=getSAM(line);
                    
                }
            }
            strcpy(qname,read1->qname);
            flag_segment=read1->flag & 192;
            
            while(strcmp(qname,read1->qname)==0 && (read1->flag & 192) == flag_segment)
            //Continue reading the entire file and search if the same read appears elsewhere?
            {
                mixedReads2.push_back(read1);
                if(fgets(line, MAX_FILE_READ, mapFile)!=NULL)
                {
                    //Enters here only once                    
                    //printf("THere? %d\n",it);
                    read1=getSAM(line);
                    
                }
                else
                {                    
                    //Never in this block
                    end=1;
                    break;
                }
            }
            //They are called mixed because they are not part of the same pair that mapped properly
            printMixedVectors(linkingFile,singletonFile,unFile);
            
            if(end==1)
            {
                //Never in this block
                break;
            }
            
        }
        
        if(end==1)
        {
            //Never in this block
            break;
        }


		fgets(line, MAX_FILE_READ, mapFile);
        
		SAM *read2=getSAM(line);
        
        
		if(strcmp(read1->qname,preqname1)!=0 || strcmp(read2->qname,preqname2)!=0)
		{
			strcpy(preqname1,read1->qname);
			strcpy(preqname2,read2->qname);
			printVectors(outFile, unFile);
			reads1.push_back(read1);
			reads2.push_back(read2);
		}
		else
		{
			reads1.push_back(read1);
			reads2.push_back(read2);
			//cout<<"ever?? here??\n";//Never Here
		}
		//if(it%10000==0)
			//cout<<"Iteration:- "<<it<<endl;
	}
    //======================Ending the read of SAM file result.sam line by line====================================
   
    //printvector runs on reads1 and reads2 vector, but on first call they are empty, see at line 1139
    //after the call, we pushback??why??
    //But anyways, is that why there is this extra one at the last??
    printVectors(outFile, unFile);
    
    
    for(int i=0;i<noContigs;i++)
        fprintf(scaffoldFile,";%s",contigNames[i]);
    fprintf(scaffoldFile,";\n");
    
    
	for(int i=0;i<noContigs;i++)
	{
        fprintf(scaffoldFile,"%s;",contigNames[i]);
		for(int j=0;j<noContigs;j++)
		{
            if(i==j)
                fprintf(scaffoldFile,"%d;",0);
			else if(readCounts[i][j]>10)
                fprintf(scaffoldFile,"%d;",readCounts[i][j]);
            else
                fprintf(scaffoldFile,"%d;",0);
            //			if(i!=j && readCounts[i][j]>0)
            //				cout<<contigNames[i]<<" "<<contigNames[j]<<" "<<readCounts[i][j]<<endl;
		}
		fprintf(scaffoldFile,"\n");
	}
    
	fclose(mapFile);
	fclose(outFile);
	fclose(unFile);
	fclose(scaffoldFile);
    fclose(gapInfoFile);
	
    //    cout<<totalCount<<"\t"<<unCount<<"\t"<<maxReadLength<<"\t"<<MAX_FRAGMENT_SIZE<<endl;
    
    fprintf(statFile,"%ld %ld %ld %ld",totalCount, unCount, maxReadLength, MAX_FRAGMENT_SIZE);//doesn't match with figbird.cpp
    fclose(statFile);
    
	FILE *countFile=fopen("counts.txt","w");
    
	for(int i=0;i<noContigs;i++)
	{
		fprintf(countFile,"%d %ld %lf\n",i, contigLengths[i], contigReadCounts[i]);
        //cout<<contigNames[i]<<"\t"<<contigLengths[i]<<"\t"<<contigReadCounts[i]<<endl;
	}
	fclose(countFile);    
    //cout<<"Done with main.cpp successfully"<<endl;
	return 0;
}

