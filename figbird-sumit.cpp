//
//  main.cpp
//  gapFiller
//
//  Created by Atif Rahman on 2/7/16.
//  Copyright (c) 2016 Atif Rahman. All rights reserved.
//


#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <cmath>
#include <cfloat>
#include<cstdlib>
using namespace std;

#define _USE_MATH_DEFINES
#include <pthread.h>
#include <math.h>
#include "cgal.h"
#include <time.h>

#define MAX_REC_LEN 1024

int noContigs=0;
int noReads=0;
long int contigLength=0;
int maxReadLength=0;

long int totalContigLength=0;
vector<char*> contigs;
vector<char*> contigNames;
vector<long int> contigLengths;
double *contigReadCounts;
double *nc;



FILE *contigFile;
FILE *mapFile;
FILE *summaryFile;
FILE *outFile;
FILE *countFile;

FILE *filledContigFile;

char *contigFileName;
char const *mapFileName;

long int *insertCounts;
int maxInsertSize=0;
int MAX_INSERT_SIZE;
double insertSizeMean;
double insertSizeVar;
double insertSizeSD;
int insertSizeMode;
int insertCutoffMax=0;
int insertCutoffMin=0;
int insertThresholdMax=0;
int insertThresholdMin=0;
int tolerance;
int insertCountMax;

int priorSize=0;
int f;
long int *insertCountsMapped;
long int *insertCountsUnmapped;

long int errorTypes[5][5];
long int baseCounts[5];
long int *errorPos;
long int *inPos;
long int *inLengths;
long int *delPos;
long int *delLengths;
long int *readLengths;

long int *effectiveLengths;

double errorTypeProbs[5][5];
double baseErrorRates[5];
double *errorPosDist;
double *inPosDist;
double *inLengthDist;
double *delPosDist;
double *delLengthDist;
double *insertLengthDist;
double *insertLengthDistSmoothed;

int windowSize=12;

double *noErrorProbs;
int tmpCount=0;
int toAlign;

long int erroredReads=0;
long int uniqueMappedReads=0;
long int discardedReads=0;
long int totalCount,unCount;

char tempCigar[500], tempMD[500];
char noErrorCigar[500], noErrorMD[500];

int maxDistance=200;

int charCodes[256];

pthread_mutex_t overlaps_mutex = PTHREAD_MUTEX_INITIALIZER;

int isJump=0;
int isConservative=0;
int isGaussian=0;
double gaussianMean;
double gaussianSD;
double inputMean;

int MIN_GAP=-500;
int MAX_GAP=5000;
int MIN_READ=0;
int MIN_READ_JOIN=1;
int NO_INTERVALS=0;


double *insertProbsSums;

long int gapProbs[1000];
int gapProbCutOff;


FILE *gapFile;
FILE *gapInfoFile;
FILE *gapOutFile;
//=======================================================================

struct MAP
{
	double errorProb;
	int insertSize;
    long int pos1;
    long int pos2;
    long int contigNo1;
    long int contigNo2;
    int readLength1;
    int readLength2;
    bool isSameStrand;
};

struct InsertTable
{
	int insertSize;
    int count;
};


void initInsertCounts(int max)
{
	maxInsertSize=max;
	insertCounts=new long int[maxInsertSize];
	for(int i=0;i<maxInsertSize;i++)
	{
		insertCounts[i]=1;
	}
}

void updateInsertCounts(int index)
{
    
	if(index<=0)
		return;
	if(index<maxInsertSize)
	{
		insertCounts[index]++;
	}
	else
	{
		
		if(index>MAX_INSERT_SIZE)
		{
			discardedReads++;
			return;
		}
		int tempInsertSize=max(maxInsertSize*2,index);
		long int *tempCounts=new long int[maxInsertSize];
		for(int i=0;i<maxInsertSize;i++)
		{
			tempCounts[i]=insertCounts[i];
		}
		insertCounts=new long int[tempInsertSize];
		for(int i=0;i<maxInsertSize;i++)
		{
			insertCounts[i]=tempCounts[i];
		}
		for(int i=maxInsertSize;i<tempInsertSize;i++)
		{
			insertCounts[i]=1;
		}
        
		insertCounts[index]++;
		maxInsertSize=tempInsertSize;
		delete []tempCounts;
		
	}
    
}

void initErrorTypes(int readLength)
{
	for(int i=0;i<5;i++)
		for(int j=0;j<5;j++)
			errorTypes[i][j]=1;
    
	for(int i=0;i<5;i++)
		baseCounts[i]=1;
    
	errorPos=new long int[readLength];
	inPos=new long int[readLength];
	inLengths=new long int[readLength];
	delPos=new long int[readLength];
	delLengths=new long int[readLength];
	readLengths=new long int[readLength];
	
	for(int i=0;i<readLength;i++)
	{
		errorPos[i]=1;
		inPos[i]=1;
		inLengths[i]=1;
		delPos[i]=1;
		delLengths[i]=1;
		readLengths[i]=0;
	}
}


int getLength(char *read)
{
	int i=0;
	while(read[i])
	{
		if(read[i]=='A')
			baseCounts[0]++;
		else if(read[i]=='C')
			baseCounts[1]++;
		else if(read[i]=='G')
			baseCounts[2]++;
		else if(read[i]=='T')
			baseCounts[3]++;
		else
			baseCounts[4]++;
        
		i++;
	}
	
	return i;
}

long int getContigNo(char *contigName)
{
    
	return atol(contigName);
    
    for(long int i=0;i<contigNames.size();i++)
    {
		if(strcmp(contigNames[i],contigName)==0)
			return i;
	}
	return -1;
    
}

void processErrorTypes(char *cigar, char *md, char *read, int strandNo)
{
    
	int readLength=getLength(read);
	readLengths[readLength-1]++;
    
	if(strcmp(md,noErrorCigar)!=0)
		erroredReads++;
	else
		return;
    
    
	unsigned long mdLength=strlen(md)-5;
	unsigned long tempLength=0;
    
	char *temp;
	int index=0,totalLength=0;
	
    
	int curIndex=0;
	int *inserts=new int[readLength];
	
	for(int i=0;i<readLength;i++)
	{
		inserts[i]=0;
	}
    
	unsigned long cigarLength=strlen(cigar);
    //	char *tempCigar=new char[cigarLength];
	char cigarChar;
    
	strcpy(tempCigar,cigar);
    
    
	temp=strtok(tempCigar,"IDMS^\t\n ");
    
	while(temp!=NULL)
	{
        
		tempLength=atoi(temp);
		totalLength+=strlen(temp);
		cigarChar=cigar[totalLength];
        
		if(cigarChar=='M')
		{
			index+=tempLength;
			curIndex+=tempLength;
		}
		else if(cigarChar=='I' || cigarChar=='S')
		{
			if(strandNo==0)
			{
				inPos[index]++;
				inLengths[tempLength-1]++;
                
			}
			else
			{
                
				inPos[readLength-index-1]++;
				inLengths[tempLength-1]++;
			}
            
			inserts[curIndex]=tempLength;
            
			index+=tempLength;
		}
		else if(cigarChar=='D' )
		{
			if(strandNo==0)
			{
				delPos[index]++;
				delLengths[tempLength-1]++;
                
			}
			else
			{
                
				delPos[readLength-index-1]++;
				delLengths[tempLength-1]++;
			}
		}
		totalLength++;
		temp=strtok(NULL,"IDMS^\t\n ");
	}
    
	
	strcpy(tempMD,md);
    
	strtok(tempMD,":");
	strtok(NULL,":");
    
	
	index=0,totalLength=0,tempLength=0;
    
	int f, t;
    
	while((temp=strtok(NULL,"ACGTN^\t\n "))!=NULL)
	{
		tempLength=strlen(temp);
        
        
		totalLength+=tempLength;
		
        
		if(totalLength<mdLength)
		{
			char from=md[5+totalLength];
			
			
			if(from=='^')
			{
				totalLength++;
				index+=atoi(temp);
				for(int i=totalLength;i<mdLength;i++)
				{
					from=md[5+totalLength];
					if(from=='A' || from=='C' || from=='G' || from=='T'|| from=='N')
						totalLength++;
					else
						break;
                    
				}
			}
			else if(from=='A' || from=='C' || from=='G' || from=='T'|| from=='N')
			{
				totalLength++;
				index+=atoi(temp)+1;
				
				
				
				curIndex=0;
				for(int i=0;i<index;i++)
				{
					curIndex+=inserts[i];
				}
				char to=read[index-1+curIndex];
                
				if(strandNo==0)
					errorPos[index-1+curIndex]++;
				else
					errorPos[readLength-index-curIndex]++;
                
                
				switch(from)
				{
					case 'A':
						f=0;
						break;
					case 'C':
						f=1;
						break;
					case 'G':
						f=2;
						break;
					case 'T':
						f=3;
						break;
					default:
						f=4;
				}
                
				switch(to)
				{
					case 'A':
						t=0;
						break;
					case 'C':
						t=1;
						break;
					case 'G':
						t=2;
						break;
					case 'T':
						t=3;
						break;
					default:
						t=4;
				}
                
				if(f==t)
				{
                    
				}
				else
                    
					errorTypes[f][t]++;
                
			}
			else
				break;
		}
		
	}
	delete []inserts;
    
}

double dnorm(double x,double mean, double variance)
{
	double val=1/sqrt(M_PI*2*variance);
	val*=exp(-((x-mean)*(x-mean))/(2*variance));
	return val;
}


void computeProbabilites()
{
	
	int errorCount=0;
    
	for(int i=0;i<5;i++)
	{
		errorCount=0;
		for(int j=0;j<5;j++)
		{
			errorCount+=errorTypes[i][j];
		}
		for(int j=0;j<5;j++)
		{
			errorTypeProbs[i][j]=(double)errorTypes[i][j]/errorCount;
		}
		
		baseErrorRates[i]=errorCount/(double)baseCounts[i];
	}
    
	double sum=0;
	for(int i=0;i<4;i++)
		sum+=baseErrorRates[i];
    
	for(int i=0;i<4;i++)
	{
		baseErrorRates[i]=4*baseErrorRates[i]/sum;
	}
    
	baseErrorRates[4]=1;
    
	for(int i=maxReadLength-1;i>0;i--)
	{
		readLengths[i-1]=readLengths[i]+readLengths[i-1];
	}
    
	errorPosDist=new double[maxReadLength];
    
	for(int i=0;i<maxReadLength;i++)
	{
		errorPosDist[i]=(double)errorPos[i]/readLengths[i];
	}
    
	inPosDist=new double[maxReadLength];
    
	for(int i=0;i<maxReadLength;i++)
	{
		inPosDist[i]=(double)inPos[i]/readLengths[i];
	}
    
	inLengthDist=new double[maxReadLength];
    
	int inCount=0;
    
	for(int i=0;i<maxReadLength;i++)
	{
		inCount+=inLengths[i];
	}
	
	for(int i=0;i<maxReadLength;i++)
	{
		inLengthDist[i]=(double)inLengths[i]/inCount;
	}
    
	delPosDist=new double[maxReadLength];
    
	for(int i=0;i<maxReadLength;i++)
	{
		delPosDist[i]=(double)delPos[i]/readLengths[i];
	}
    
	delLengthDist=new double[maxReadLength];
    
	int delCount=0;
    
	for(int i=0;i<maxReadLength;i++)
	{
		delCount+=delLengths[i];
	}
	
	for(int i=0;i<maxReadLength;i++)
	{
		delLengthDist[i]=(double)delLengths[i]/delCount;
	}
	
    
	insertLengthDist=new double[maxInsertSize];
    
	long int insCount=discardedReads;
    
	sum=0;
    
	for(int i=0;i<maxInsertSize;i++)
	{
		insCount+=(insertCounts[i]-1);
		sum+=i*(insertCounts[i]-1);
        
	}
	insertSizeMean=sum/insCount;

	sum=0;
    
	for(int i=0;i<maxInsertSize;i++)
	{
		insertLengthDist[i]=(double)insertCounts[i]/insCount;
        
		sum+=(insertCounts[i]-1)*(insertSizeMean-i)*(insertSizeMean-i);
	}
    
	insertSizeVar=sum/insCount;
    
    insertSizeSD=sqrt(insertSizeVar);
    
	noErrorProbs=new double[maxReadLength];
    
	double noErrorProb=1.0;
    
    
	for(int i=0;i<maxReadLength;i++)
	{
		noErrorProb*=(1-errorPosDist[i]-inPosDist[i]-delPosDist[i]);
		noErrorProbs[i]=noErrorProb;
  	}
    
	effectiveLengths=new long int[maxInsertSize];
    
	for(int i=0;i<maxInsertSize;i++)
	{
		effectiveLengths[i]=-1;
	}
    
	long int totalContigLength=0;
	for(int i=0;i<contigLengths.size();i++)
	{
		totalContigLength+=contigLengths[i];
	}
	effectiveLengths[0]=totalContigLength;
	
    
    insertCountsMapped=new long int[maxInsertSize];
    insertCountsUnmapped=new long int[maxInsertSize];
    
    
    for(long int i=0;i<maxInsertSize;i++)
    {
        insertCountsMapped[i]=0;
        insertCountsUnmapped[i]=0;
    }
    
    insertLengthDistSmoothed=new double[maxInsertSize];
    double windowSum=0;
    
    for(int i=0;i<windowSize;i++)
    {
        insertLengthDistSmoothed[i]=insertLengthDist[i];
        
    }
    
    for(int i=0;i<2*windowSize+1;i++)
    {
        windowSum+=insertLengthDist[i];
        
    }
    insertLengthDistSmoothed[windowSize]=windowSum/(2*windowSize+1);
    
    for(int i=windowSize+1;i<maxInsertSize-windowSize;i++)
    {
        windowSum-=insertLengthDist[i-windowSize-1];
        windowSum+=insertLengthDist[i+windowSize];
        insertLengthDistSmoothed[i]=windowSum/(2*windowSize+1);
    }
    for(int i=maxInsertSize-windowSize;i<maxInsertSize;i++)
    {
        insertLengthDistSmoothed[i]=insertLengthDist[i];
    }
    
	for(int i=0;i<maxInsertSize;i++)
	{
		insertLengthDistSmoothed[i]=insertLengthDistSmoothed[i]-1/(double)(insCount)+(1/(double)maxInsertSize)/(double)(insCount+1);
        
	}
    
    int count=0;
    
    for(int i=insertSizeMean;i<maxInsertSize;i++)
    {
        if(insertCounts[i]<=1)
        {
            count++;
            if(count==10)
            {
                insertCutoffMax=i;
                break;
            }
        }
        else
        {
            count=0;
        }
        
    }
    
    count=0;
    
    for(int i=insertSizeMean;i>=0;i--)
    {
        if(insertCounts[i]<=1)
        {
            count++;
            if(count==10)
            {
                insertCutoffMin=i;
                break;
            }
        }
        else
        {
            count=0;
        }
        
    }
    
    insertCountMax=0;
    for(int i=0;i<maxInsertSize;i++)
    {
        if(insertCounts[i]>insertCountMax)
        {
            insertCountMax=insertCounts[i];
            insertSizeMode=i;
            
        }
    }
    
    
    count=0;
    
    for(int i=insertSizeMean;i<maxInsertSize;i++)
    {
        if(insertCounts[i]<=max(insertCountMax/1000,2))
        {
            count++;
            if(count==2)
            {
                insertThresholdMax=i;
                break;
            }
        }
        else
        {
            count=0;
        }
        
    }
    
    count=0;
    
    for(int i=insertSizeMean;i>=0;i--)
    {
        if(insertCounts[i]<=max(insertCountMax/1000,2))
        {
            count++;
            if(count==2)
            {
                insertThresholdMin=i;
                break;
            }
        }
        else
        {
            count=0;
        }
        
    }

    //mean used instead of mode
    
    double insertSum=0;
    double insertCount=0;
	
    for(int i=insertCutoffMin;i<insertCutoffMax;i++)
    {
        insertCount+=insertCounts[i]-1;
        insertSum+=(insertCounts[i]-1)*i;
    }
    insertSizeMode=insertSum/insertCount;
    
	double leftSD, rightSD;
    
	insertSum=0;
    insertCount=0;
    
	for(int i=insertSizeMean+1;i<maxInsertSize;i++)
	{
		insertSum=insertSum+(insertCounts[i]-1)*(i-insertSizeMean)*(i-insertSizeMean);
		insertCount+=(insertCounts[i]-1);
	}
	rightSD=sqrt(insertSum/insertCount);

	insertSum=0;
    insertCount=0;
	for(int i=max((int)(insertSizeMean-10*rightSD),0);i<insertSizeMean;i++)
	{
		insertSum=insertSum+(insertCounts[i]-1)*(insertSizeMean-i)*(insertSizeMean-i);
		insertCount+=(insertCounts[i]-1);
	}
	leftSD=sqrt(insertSum/insertCount);
    
	if((rightSD>1000 || leftSD>1000) && isGaussian==0)
	{
		cerr<<"Switching to conservative mode"<<endl;
		isConservative=1;
	}

	insertThresholdMin=max((int)(insertSizeMean-2.5*leftSD),1);
	insertThresholdMax=min((int)(insertSizeMean+2.5*rightSD),maxInsertSize);
    
	if(isGaussian==1)
	{
		insertThresholdMin=max((int)(gaussianMean-2.5*gaussianSD),1);
		insertThresholdMax=min((int)(gaussianMean+2.5*gaussianSD),maxInsertSize);
		insertSizeMode=gaussianMean;
        /*
         for(int i=0;i<insertThresholdMin;i++)
         {
         insertLengthDistSmoothed[i]=(1/(double)maxInsertSize)/(double)(insCount+1);
         }
         for(int i=insertThresholdMin;i<insertThresholdMax;i++)
         {
         insertLengthDistSmoothed[i]=dnorm(i,gaussianMean,gaussianSD*gaussianSD)+(1/(double)maxInsertSize)/(double)(insCount+1);
         }
         for(int i=insertThresholdMax;i<maxInsertSize;i++)
         {
         insertLengthDistSmoothed[i]=(1/(double)maxInsertSize)/(double)(insCount+1);
         }
         */
		for(int i=0;i<maxInsertSize;i++)
		{
			insertLengthDistSmoothed[i]=dnorm(i,gaussianMean,gaussianSD*gaussianSD)+(1/(double)maxInsertSize)/(double)(insCount+1);
		}
        
	}

    insertCutoffMax=insertThresholdMax;
    insertCutoffMin=insertThresholdMin;
    
//	MAX_GAP=insertCutoffMax;
    
	
}

void processMapping(char *line)
{
	
	char * temp;
	char *qname, *rname, *mapq;
	int	pos,flag;
	char * cigar, * readString; // * md, *nhstring;
	long int contigNo;
    
    
	char md[500];
	char nhstring[500];
    
	int nh;
    
	int strandNo=0;
    
    
	qname=strtok(line,"\t");
	
    
	temp=strtok(NULL,"\t");
	flag=atoi(temp);
	
    
	strandNo=(flag&16)>>4;
    
    rname=strtok(NULL,"\t");
	
    
	temp=strtok(NULL,"\t");
	pos=atoi(temp);
	
	
	cigar=strtok(NULL,"\t");
	
	
	temp=strtok(NULL,"\t");
	
	readString=strtok(NULL,"\t");
    
	int insertSize=atoi(temp);
    
	while((temp=strtok(NULL,"\t\n"))!=NULL)
	{
		if(temp[0]=='M' && temp[1]=='D')
		{
			strcpy(md,temp);
		}
		else if(temp[0]=='I' && temp[1]=='H')
		{
			strcpy(nhstring,(temp+5));
			nh=atoi(nhstring) ;
		}
        
	}
    
	
	if(nh==1 && md[5]!='^')
	{
		contigNo=getContigNo(rname);
		if(isGaussian==1)
		{
			updateInsertCounts(insertSize);
		}
		else
		{
			if(contigLengths[contigNo]>inputMean)
				updateInsertCounts(insertSize);
		}
		processErrorTypes(cigar,md,readString,strandNo);
		uniqueMappedReads++;
	}
	
    
}

long int getEffectiveLength(int insertSize)
{
	if(insertSize<0)
		return effectiveLengths[0];
    
	if(insertSize>=maxInsertSize)
	{
		long int effectiveLength=0;
		for(int i=0;i<contigLengths.size();i++)
		{
			if(contigLengths[i]>=insertSize)
				effectiveLength+=(contigLengths[i]-insertSize+1);
		}
		return effectiveLength;
        
	}
	if(effectiveLengths[insertSize]==-1)
	{
		long int effectiveLength=0;
		for(int i=0;i<contigLengths.size();i++)
		{
			if(contigLengths[i]>=insertSize)
				effectiveLength+=(contigLengths[i]-insertSize+1);
		}
		effectiveLengths[insertSize]=effectiveLength;
	}
	return effectiveLengths[insertSize];
}

long double computeErrorProb(char *cigar, char *md, char *read, int strandNo)
{
    
    unsigned long readLength=strlen(read);
	
    
	long double errorProb=noErrorProbs[readLength-1];
    
    
    
	if(md[5]=='^')
		return errorProb;
	
	
	char tempMD[1000], tempCigar[1000];
    
	unsigned long mdLength=strlen(md)-5;
	unsigned long tempLength=0;
	
	char *temp;
	int index=0,totalLength=0;
	
    
	int curIndex=0;
	int *inserts=new int[readLength];
	
	for(int i=0;i<readLength;i++)
	{
		inserts[i]=0;
	}
    
    //	int cigarLength=strlen(cigar);
	char cigarChar;
    
	strcpy(tempCigar,cigar);
    
	temp=strtok(tempCigar,"IDM^\t\n ");
    
	while(temp!=NULL)
	{
        
		tempLength=atoi(temp);
		totalLength+=strlen(temp);
		cigarChar=cigar[totalLength];
        
		if(cigarChar=='M')
		{
			index+=tempLength;
			curIndex+=tempLength;
		}
		else if(cigarChar=='I')
		{
			unsigned long i;
			if(strandNo==0)
			{
				//look up insert probs
				i=index;
				
			}
			else
			{
				i=readLength-index-1;
			}
            
			errorProb=errorProb*inPosDist[i]*inLengthDist[tempLength-1]/(1-errorPosDist[i]-inPosDist[i]-delPosDist[i]);
            
			inserts[curIndex]=tempLength;
            
			index+=tempLength;
		}
		else if(cigarChar=='D')
		{
			unsigned long i;
			if(strandNo==0)
			{
				i=index;
                //	look up delete probs
			}
			else
			{
				i=readLength-index-1;
			}
            
			errorProb=errorProb*delPosDist[i]*delLengthDist[tempLength-1]/(1-errorPosDist[i]-inPosDist[i]-delPosDist[i]);
		}
		totalLength++;
		temp=strtok(NULL,"IDM^\t\n ");
	}
    
    
	strcpy(tempMD,md);
    
	strtok(tempMD,":");
	strtok(NULL,":");
    
	index=0,totalLength=0,tempLength=0;
    
	int f, t;
    
	while((temp=strtok(NULL,"ACGTN^\t\n "))!=NULL)
	{
		tempLength=strlen(temp);
        
		totalLength+=tempLength;
		
		if(totalLength<mdLength)
		{
			char from=md[5+totalLength];
            
			if(from=='^')
			{
				totalLength++;
				index+=atoi(temp);
				for(int i=totalLength;i<mdLength;i++)
				{
					from=md[5+totalLength];
					if(from=='A' || from=='C' || from=='G' || from=='T'|| from=='N')
						totalLength++;
					else
						break;
				}
			}
			else if(from=='A' || from=='C' || from=='G' || from=='T'|| from=='N')
			{
				totalLength++;
				index+=atoi(temp)+1;
				
                
				curIndex=0;
				for(int i=0;i<index;i++)
				{
					curIndex+=inserts[i];
				}
				char to=read[index-1+curIndex];
                
				int i;
				if(strandNo==0)
					i=index-1+curIndex;
				else
					i=readLength-index-curIndex;
                
                
				errorProb=errorProb*errorPosDist[i]/(1-errorPosDist[i]-inPosDist[i]-delPosDist[i]);
                
                
				switch(from)
				{
					case 'A':
						f=0;
						break;
					case 'C':
						f=1;
						break;
					case 'G':
						f=2;
						break;
					case 'T':
						f=3;
						break;
					default:
						f=4;
				}
                
				switch(to)
				{
					case 'A':
						t=0;
						break;
					case 'C':
						t=1;
						break;
					case 'G':
						t=2;
						break;
					case 'T':
						t=3;
						break;
					default:
						t=4;
				}
                
				if(f==t)
				{
					
                    
				}
				else
				{
					//errorTypeProb
					errorProb*=baseErrorRates[f]*errorTypeProbs[f][t];
				}
			}
			else
				break;
            
		}
		
	}
	
	delete []inserts;
	return errorProb;
}


double computeLikelihood(char const *file)
{
    mapFile=fopen(file, "r");
	char *line1= new char[MAX_REC_LEN];
	char *line2= new char[MAX_REC_LEN];
    
	char *qname1,*qname2,preqname1[500],preqname2[500];
    
	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line1[0]);
    
	long double sum=0.0;
	long double logsum=0.0;
    
	char * temp;
	char *rname1, *rname2;
	int	pos1,pos2,flag,strandNo1, strandNo2, insertSize1, insertSize2;
	char *cigar1, *cigar2, *readString1, *readString2, md1[1000], md2[1000];
    
    
	long double insertSizeProb;
    
	long double errorProb1, errorProb2;
	
	long double gapProb;
    
    int tempInsertSize=0;
    long double tempProb=0;
    
	preqname1[0]=preqname2[0]='*';
	preqname1[1]=preqname2[1]=0;
    
    
	int it=0;
    
	while(fgets(line1, MAX_FILE_READ, mapFile)!=NULL)
	{
		if(line1[0]=='@')
			continue;
        //????
		if(fgets(line2, MAX_FILE_READ, mapFile)==NULL)
			break;
        
		qname1=strtok(line1,"\t");
		
        temp=strtok(NULL,"\t");
        flag=atoi(temp);
		
        
		strandNo1=(flag&16)>>4;
        
        temp=strtok(NULL,"\t");
		
		temp=strtok(NULL,"\t");
		pos1=atoi(temp);
        
		
		cigar1=strtok(NULL,"\t");
        
		
		temp=strtok(NULL,"\t");
		insertSize1=atoi(temp);
        
        readString1=strtok(NULL,"\t");
        
        
		while((temp=strtok(NULL,"\t\n"))!=NULL)
		{
			if(temp[0]=='M' && temp[1]=='D')
			{
				strcpy(md1,temp);
			}
		}
        
        //        		cout<<insertSize1<<" "<<cigar1<<" "<<md1<<endl;
		
        //second of the pair
        
		qname2=strtok(line2,"\t");
		temp=strtok(NULL,"\t");
		flag=atoi(temp);
		
        
		strandNo2=(flag&16)>>4;
        
        temp=strtok(NULL,"\t");
		
		temp=strtok(NULL,"\t");
		pos2=atoi(temp);
        
		cigar2=strtok(NULL,"\t");
        
		temp=strtok(NULL,"\t");
		
		insertSize2=atoi(temp);
        
		readString2=strtok(NULL,"\t");
        
        
		while((temp=strtok(NULL,"\t\n"))!=NULL)
		{
			if(temp[0]=='M' && temp[1]=='D')
			{
				strcpy(md2,temp);
			}
		}
        
        //cout<<insertSize2<<" "<<cigar2<<" "<<md2<<endl;

		int insertSize=max(insertSize1, insertSize2);
		
        
		insertSizeProb=0;
        
		if(insertSize>=0 && insertSize<maxInsertSize)
		{
			insertSizeProb=insertLengthDist[insertSize];
		}
        
		if(insertSizeProb==0)
		{
			insertSizeProb=1/(double)uniqueMappedReads;
		}
		
        
		errorProb1=computeErrorProb(cigar1,md1,readString1,strandNo1);
		errorProb2=computeErrorProb(cigar2,md2,readString2,strandNo2);

		long int totalEffectiveLength=getEffectiveLength(insertSize);
        
        
		long double prob=(1/(long double)(totalEffectiveLength))*insertSizeProb*errorProb1*errorProb2;

        //        cout<<errorProb1<<" "<<errorProb2<<" "<<insertSizeProb<<" "<<prob<<endl;
        if(strcmp(qname1,preqname1)==0 && strcmp(qname2,preqname2)==0)
		{
            if(tempProb<prob)
            {
                tempProb=prob;
                tempInsertSize=insertSize;
                tempInsertSize=tempInsertSize<0?0:tempInsertSize;
                gapProb=insertSizeProb*errorProb2;
            }
            
			sum+=prob;
            
		}
		else if(strcmp("*",preqname1)!=0 && strcmp("*",preqname2)!=0)
		{
			if(sum<1e-320 || isnan(sum))
			{
				sum=1e-320;
			}
			logsum+=log(sum);
            
            
            int gapIndex=-log(gapProb);
            
            if(gapIndex<1000 && gapIndex>=0)
            {
            	gapProbs[gapIndex]++;
            }
            else
            {
            	gapProbs[999]++;
            }
            
            
            if(tempInsertSize>=maxInsertSize)
                insertCountsMapped[maxInsertSize-1]++;
            else
                insertCountsMapped[tempInsertSize]++;
			
			sum=prob;
            
            tempProb=prob;
            tempInsertSize=insertSize;
            tempInsertSize=tempInsertSize<0?0:tempInsertSize;
            
            gapProb=insertSizeProb*errorProb2;
		}
		else
		{
			sum=prob;
			
            tempProb=prob;
            tempInsertSize=insertSize;
            tempInsertSize=tempInsertSize<0?0:tempInsertSize;
            
			gapProb=insertSizeProb*errorProb2;
		}
        
		strcpy(preqname1,qname1);
		strcpy(preqname2,qname2);
		it++;
        
		
		if(isinf( logsum ))
		{
			cout<<it<<endl;
			exit(1);
		}
        
        
	}
	if(sum!=0)
    {
        if(sum<1e-320 || isnan(sum))
        {
            sum=1e-320;
        }
		logsum+=log(sum);
    }
    
	fclose(mapFile);
	
	return logsum;
}

void printHelp()
{
    
	cout<<"cgal v0.9.5-beta"<<endl;
	cout<<"----------------"<<endl;
	cout<<endl;
	cout<<"cgal - computes likelihood"<<endl;
	cout<<"Usage:"<<endl;
	cout<<"cgal [options] <contigfile>"<<endl;
	cout<<endl;
	cout<<"Required arguments:"<<endl;
	cout<<"<contigfile>\t Assembly file in FASTA format"<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
	cout<<"-h [--help]\t\t Prints this message"<<endl;
	cout<<endl;
	cout<<"Output: "<<endl;
	cout<<"(In file out.txt) <numberContigs> <totalLikelihood> <mappedLikelihood> <unmappedLikelihood> <noReads> <noReadsUnmapped>"<<endl;
	cout<<"<numberContigs>\t\t Number of contigs"<<endl;
	cout<<"<totalLikelihood>\t Total log likelihood value"<<endl;
	cout<<"<mappedLikelihood>\t Likelihood value of reads mapped by the mapping tool"<<endl;
	cout<<"<unmappedLikelihood>\t Likelihood value corresponding to reads not mapped by alignment tool"<<endl;
	cout<<"<noReads>\t\t Total number of paired-end reads"<<endl;
	cout<<"<noReadsUnmapped>\t Number of reads not mapped by the alignment tool"<<endl;
	cout<<endl;
	exit(1);
    
}

void reverse(char *reverse, char *read)
{
    
	char ch='A';
	int readLength=strlen(read);
	for(int i=1; i<=readLength;i++)
	{
		ch=read[i-1];
		if(ch=='A' || ch=='a')
			reverse[readLength-i]='T';
		else if(ch=='C' || ch=='c')
			reverse[readLength-i]='G';
		else if(ch=='G' || ch=='g')
			reverse[readLength-i]='C';
		else if(ch=='T' || ch=='t')
			reverse[readLength-i]='A';
		else
			reverse[readLength-i]='N';
	}
	reverse[readLength]='\0';
    
}


int getDistance(char *s, char * t, int sStart, int sEnd, int tStart, int tEnd, int ** dis)
{
	int m=sEnd-sStart+1;
    int n=tEnd-tStart+1;
	
    
	int val1;
	int val2;
	int val3;
    
	
    int i,j;
    
	for(i=0;i<=m;i++)
        dis[i][0]=i;
    
	for(j=0;j<=n;j++)
        dis[0][j]=j;
    
	for(i=1;i<=m;i++)
	{
		for(j=1;j<=n;j++)
		{
			val2=dis[i-1][j]+1;
			val3=dis[i][j-1]+1;
			val2=val2<val3?val2:val3;
			
            val1=dis[i-1][j-1];
            if(s[i-1+sStart]!=t[j-1]+tStart)
			{
				val1++;
                
			}
			val3=val1<val2?val1:val2;
			
            dis[i][j]=val3;
			
		}
		
	}
    
    
    
    return 	dis[m][n];
    
}

int getMismatch(char *s, char * t, int sStart, int sEnd, int tStart, int tEnd, double error)
{
	int m=sEnd-sStart+1;
    int n=tEnd-tStart+1;
    
	int mismatches=0;
    
	for(int i=0;i<m,i<n;i++)
	{
		if(s[sStart+i]!=t[tStart+i])
		{
			mismatches++;
			if(mismatches>error)
				return mismatches;
		}
	}
	return mismatches;
    
}


int getOverlap(char *s, char *t, double error)
{
    int sLength=strlen(s);
    int tLength=strlen(t);
    int distance=0;
    
    int **dis;
    
    dis=new int*[sLength+1];
    
    for(int i=0;i<=sLength;i++)
    {
        dis[i]=new int[tLength+1];
    }
    
    
    for(int i=sLength>tLength?sLength-tLength:0;i<sLength;i++)
    {
        distance=getDistance(s,t , i, sLength-1, 0, sLength-i-1,dis);
        
        if(distance<=error*(sLength-i))
        {
            for(int j=0;j<=sLength;j++)
            {
                delete [] dis[j];
            }
            
            delete [] dis;
            
            return (sLength-i);
        }
    }
    
    for(int i=0;i<=sLength;i++)
    {
        delete [] dis[i];
    }
    
    delete [] dis;
    
    return 0;
}

class GapFiller
{
    double **countsGap;
    double **probsGap;
    double **errorProbsGap;

	char *bestString;
	char *secondBestString;

    int contigNo;
    long int gapStart;
    int gapLength;

    int originalGap;

    long int startPos, endPos;
	char gapFileName[500];
    char *concensus;



    int maxGap;


public:

    void allocate(int maxGap)
    {
        this->maxGap=maxGap;
        
        int maxSize=maxGap+2*maxDistance;
        concensus=new char[maxGap+1];


        countsGap=new double *[maxSize];
        probsGap=new double *[maxSize];
        errorProbsGap=new double *[maxSize];


        bestString=new char[maxGap+1];
        secondBestString=new char[maxGap+1];

        for(long int i=0;i<maxSize;i++)
        {
            countsGap[i]=new double[5];
            probsGap[i]=new double[5];
            errorProbsGap[i]=new double[5];

        }
    }

	void computeProbsGap()
	{
	    double total=0;
	    double Ncount=0;
	    for(long int i=0;i<endPos-startPos;i++)
	    {
	        total=0;
	        for(int j=0;j<5;j++)
	        {
	            total=total+countsGap[i][j];
	        }
	        Ncount=countsGap[i][4];

	        for(int j=0;j<4;j++)
	        {
	            probsGap[i][j]=((countsGap[i][j]+Ncount/4)/total);
	        }
	        probsGap[i][4]=0;
	    }
	}
	
	void computeErrorProbsGap()
	{
	    double sum=0;
	    
	    for(long int i=0;i<endPos-startPos;i++)
	    {
	        for(int j=0;j<5;j++)
	        {
	            sum=0;
	            for(int k=0;k<4;k++)
	            {
	                if(j==k)
	                    continue;
	                
	                sum+=probsGap[i][k]*(errorTypeProbs[k][j]);
	            }
	            errorProbsGap[i][j]=sum;
	        }
	    }
	}
	
	void initGapFiller(int contigNo,long int gapStart, int originalGap, char * gapFileName)
	{
	    this->gapStart=gapStart;
	    this->contigNo=contigNo;
	    
	    this->originalGap=originalGap;
	    
	    strcpy(this->gapFileName,gapFileName);
	}
	
	void initialize(int gapLength)//called using gapestimate
	{
		this->gapLength=gapLength;
	    startPos=(this->gapStart-maxDistance)<0?0:(this->gapStart-maxDistance);
	    endPos=this->gapStart+this->gapLength+maxDistance;

	    //cout<<"endpos = "<<endPos<<"start = "<<startPos<<endl;
	    for(long int i=0;i<endPos-startPos;i++)
	    {
	        for(int j=0;j<5;j++)
	        {
	            countsGap[i][j]=0;
	            probsGap[i][j]=0;
	            errorProbsGap[i][j]=0;
	        }
	    }
        /*
	    for(long int i=startPos;i<gapStart;i++)
	    {
	        cout<<contigs[contigNo][i];
	        
	    }
	    cout<<endl;
	    for(long int i=gapStart;i<gapStart+gapLength;i++)
	    {
	        cout<<contigs[contigNo][i];
	        
	    }
	    cout<<endl;
	    for(long int i=gapStart+gapLength;i<endPos;i++)
	    {
	        cout<<contigs[contigNo][i];
	        
	    }
	    cout<<endl;
	    */
	    char ch;
	    for(long int i=0;i<maxDistance;i++)
	    {
	        ch=contigs[contigNo][i+startPos];
	        countsGap[i][charCodes[ch]]++;
	    }
	    for(long int i=maxDistance;i<maxDistance+gapLength;i++)
	    {
	        countsGap[i][4]++;
	    }
	    for(long int i=maxDistance+gapLength;i<maxDistance+gapLength+maxDistance;i++)
	    {
	        ch=contigs[contigNo][i+startPos-gapLength+originalGap];
	        countsGap[i][charCodes[ch]]++;
	    }

	    computeProbsGap();
	    computeErrorProbsGap();

	}

	double placeReads(int ge)
	{
	    //int f=0;//old=0;new =1
	    FILE * scaffoldMap=fopen(gapFileName, "r");
	    
		char *qname1,*qname2;
	    
	    char *line1= new char[MAX_REC_LEN];
		char *line2= new char[MAX_REC_LEN];
	
		int MAX_FILE_READ=MAX_REC_LEN/sizeof(line1[0]);
	    
	    
		char * temp;
		char *rname1, *rname2;
		int	pos1,pos2,flag1, flag2, strandNo1, strandNo2, insertSize1, insertSize2;
		char *cigar1, *cigar2, *readString1, *readString2, md1[1000], md2[1000];
	    int ih1, ih2;
	    int readLength1, readLength2;
	    int charCode,index;
		
	    char tempReadString[100];
	    
	    double maxProb=0,tempProb,maxLikelihood=0;
	    
	    int insertSize,mleInsertSize,tempInsertSize;
	    
	    int maxMatch, tempMatch, maxPos;
	    
	    int isReverse=0, readIndex;
	    
	    for(long int i=maxDistance;i<endPos-startPos-maxDistance;i++)
	    {
	        for(int j=0;j<5;j++)
	        {
	            countsGap[i][j]=0;
	        }
	    }

	    int count_of_gap=1;
	    int read_no = 14;
	    while(fgets(line1, MAX_FILE_READ, scaffoldMap)!=NULL)
		{
			if(line1[0]=='@')
				continue;
	        //????
			if(fgets(line2, MAX_FILE_READ, scaffoldMap)==NULL)
				break;

			qname1=strtok(line1,"\t");
			
	        temp=strtok(NULL,"\t");

	        flag1=atoi(temp);

			strandNo1=(flag1&16)>>4;

	        rname1=strtok(NULL,"\t");

			temp=strtok(NULL,"\t");

			pos1=atoi(temp);

			cigar1=strtok(NULL,"\t");

	        temp=strtok(NULL,"\t");

			insertSize1=atoi(temp);

			readString1=strtok(NULL,"\t");
	        
	/*
			while((temp=strtok(NULL,"\t\n"))!=NULL)
			{
				if(temp[0]=='M' && temp[1]=='D')
				{
					strcpy(md1,temp);
				}
	            if(temp[0]=='I' && temp[1]=='H')
				{
					ih1=atoi(&temp[5]);
				}
			}
	*/        
	        
			
	        //second of the pair
	        
			qname2=strtok(line2,"\t");
			temp=strtok(NULL,"\t");
			flag2=atoi(temp);
			
	        
			strandNo2=(flag2&16)>>4;
	        
	        rname2=strtok(NULL,"\t");
			
			temp=strtok(NULL,"\t");
			pos2=atoi(temp);

			cigar2=strtok(NULL,"\t");

			temp=strtok(NULL,"\t");
			
			insertSize2=atoi(temp);
	        
			readString2=strtok(NULL,"\t");
	        
	  /*
	        
			while((temp=strtok(NULL,"\t\n"))!=NULL)
			{
				if(temp[0]=='M' && temp[1]=='D')
				{
					strcpy(md2,temp);
				}
	            if(temp[0]=='I' && temp[1]=='H')
				{
					ih2=atoi(&temp[5]);
				}
			}
	        
	    */    
	        readLength1=strlen(readString1);
	        readLength2=strlen(readString2);
	
	        if(strandNo1==0)
	        {
	            reverse(tempReadString, readString2);
	            strcpy(readString2,tempReadString);
	            isReverse=1;
	        }
	        else
	        {
	        	isReverse=0;
			}

	        if(f==1)maxProb=-DBL_MAX;
	        else maxProb=0;

	        maxMatch=0;
	        maxPos=0;
	        double aa=0;

	       //cout<<"Pos1 = "<<pos1<<"\tPos2 = "<<pos2<<endl;//same
	        if(pos1<gapStart)
	        {
	            insertSize=gapStart-pos1+readLength2;
	            //cout<<"In if placereads, insertsize = "<<insertSize<<endl;
	            for(long int i=gapStart-readLength2+1;i<gapStart+gapLength;i++)
	            {
	                tempInsertSize=insertSize+i-gapStart;
	                if(tempInsertSize<(insertThresholdMin) || tempInsertSize > (insertThresholdMax))
                    {

                        continue;
                    }

                    if(f==1)tempProb=log(insertLengthDistSmoothed[tempInsertSize]);
	                else tempProb=(insertLengthDistSmoothed[tempInsertSize]);

	                tempMatch=0;

	                for(int j=0;j<readLength2;j++)
	                {
	                    charCode=charCodes[readString2[j]];//taken from read2
	                    index=i+j-startPos;//this index refers to the ref genome with gap
	                    
	                    if(isReverse==1)
	                	{
	                		readIndex=readLength2-1-j;//this index refers to the read
						}
						else
						{
							readIndex=j;
						}
	                    if(charCode<4)
	                    {
	                        aa =(probsGap[index][charCode]*(1-errorPosDist[readIndex])+errorPosDist[readIndex]*errorProbsGap[index][charCode]);
	                        if(f==1)tempProb += log(aa);
	                        else tempProb*= aa;
	                    }
	                    else
	                    {
	                        if(f==1)tempProb +=log(errorPosDist[readIndex]*errorProbsGap[index][4]);
	                        else tempProb *=(errorPosDist[readIndex]*errorProbsGap[index][4]);
	                    }
	                }

	                if(tempProb>maxProb)
	                {
	                    maxMatch=tempMatch;//Not used here
	                    maxProb=tempProb;
	                    mleInsertSize=insertSize+i-gapStart;
	                    maxPos=i-startPos;
	                }

	                for(int j=0;j<readLength2;j++)
	                {
	                    if(i-startPos+j>=maxDistance && i-startPos+j<maxDistance+gapLength)
	                    {
                            if(f==1)countsGap[i-startPos+j][charCodes[readString2[j]]]+=pow(10,tempProb);//antilog
                            else countsGap[i-startPos+j][charCodes[readString2[j]]]+= tempProb;
	                    }
	                }
	            }//for end
	        }
	        else
	        {

	            insertSize=pos1-(gapStart+gapLength)+readLength1;
	            //cout<<"In else placereads, insertsize = "<<insertSize<<endl;
	            for(long int i=gapStart-readLength2+1;i<gapStart+gapLength;i++)
	            {
	                tempInsertSize=insertSize+gapStart+gapLength-i;
	                if(tempInsertSize<insertThresholdMin || tempInsertSize > insertThresholdMax)
	                    {

	                        continue;
	                    }
	                
	                if(f==1)tempProb=log(insertLengthDistSmoothed[tempInsertSize]);
	                else tempProb=(insertLengthDistSmoothed[tempInsertSize]);

	                tempMatch=0;


	                for(int j=0;j<readLength2;j++)
	                {
	                    charCode=charCodes[readString2[j]];
	                    index=i+j-startPos;
	                    
	                    
	                    if(isReverse==1)
	                	{
	                		readIndex=readLength2-1-j;
						}
						else
						{
							readIndex=j;
						}

	                    if(charCode<4)
	                    {
	                        aa=(probsGap[index][charCode]*(1-errorPosDist[readIndex])+errorPosDist[readIndex]*errorProbsGap[index][charCode]);
	                        if(f==1)tempProb += log(aa);
	                        else tempProb*= aa;
	                    }
	                    else
	                    {
	                        if(f==1)tempProb +=log(errorPosDist[readIndex]*errorProbsGap[index][4]);
                            else tempProb *=(errorPosDist[readIndex]*errorProbsGap[index][4]);
	                    }
	                }

	                if(tempProb>maxProb)
	                {
	                    maxMatch=tempMatch;
	                    maxProb=tempProb;
	                    mleInsertSize=insertSize+gapStart+gapLength-i;
	                    maxPos=i-startPos;
	                }

	                for(int j=0;j<readLength2;j++)
	                {
	                    if(i-startPos+j>=maxDistance && i-startPos+j<maxDistance+gapLength)
	                    {
	                        if(f==1)countsGap[i-startPos+j][charCodes[readString2[j]]]+=pow(10,tempProb);
	                        else countsGap[i-startPos+j][charCodes[readString2[j]]]+= tempProb;
	                    }
	                }
	            }//for end
	        }
	        //maxLikelihood += log(maxProb);
            if(f==1)maxLikelihood += (maxProb);
            else
            {
                if(maxProb>0)maxLikelihood += log(maxProb);
            }
	        count_of_gap++;
	    }//while end
	    
	    delete [] line1;
	    delete [] line2;
	    
	    fclose(scaffoldMap);
	    //cout<<"Iteration "<<ge<<" , maxlikelihood = "<<maxLikelihood<<endl;
	    return maxLikelihood;
	}

	int checkDuplicate(double *arr)
	{
	    int flag=0;
	    int size = sizeof(arr) /  sizeof(arr[0]);
	    int count[10000]={0};
        for(int i = 0; i < size; i++)
         {
             if(count[(int)arr[i]] == 1)
             {
                flag=1;
                break;
             }
             else count[(int)arr[i]]++;
         }
        return flag;
	}

	void computeSequence(int check)//based on the 2Darray countsGap
	{
	    double max=0;
	    int maxIndex=0;

        for(long int i=maxDistance;i<endPos-startPos-maxDistance;i++)
	    {
	        max=0;
	        maxIndex=0;
	        float aaa = 1;
	        int flag=-1;
	        for(int j=0;j<=4;j++)
	        {
	            if(countsGap[i][j]>max)
                {
                    max=countsGap[i][j];
                    maxIndex=j;
                    if(max == 0)flag=1;//all elements in the array is 0
                    else flag=0;
                }
             }

	        ///Check =0 when, computesequence is called from  func to compute best string, we don't
	                            //want to put N for that call
	        ///Check =1 when, computesequence is called from finalize() func to compute final string, we
            	                            //want to put N for that call

	        int flag2=0;

	        if(check == 1)flag2 = checkDuplicate(countsGap[i]);//0 means no duplicate
	        if((flag == 0 && flag2 == 0) || !check)
            {
                if(maxIndex==0)
                {
                    concensus[i-maxDistance]='A';
                }
                else if(maxIndex==1)
                {
                    concensus[i-maxDistance]='C';

                }
                else if(maxIndex==2)
                {
                    concensus[i-maxDistance]='G';

                }
                else if(maxIndex==3)
                {
                    concensus[i-maxDistance]='T';

                }
                else
                {
                    concensus[i-maxDistance]='N';
                }
            }
            else
            {
                //cout<<"Countsgap row index = "<<i<<" , flag = "<<flag<<" , flag2 = "<<flag2<<endl;
                concensus[i-maxDistance]='N';
            }
	    }
	    concensus[endPos-startPos-2*maxDistance]='\0';
	}


	int getConcensus(char *s)
	{
	    int len;
	    len = this-> gapLength;


		for(int i=0;i<len;i++)
		{
			s[i]=concensus[i];
		}

		s[len]='\0';
		return len;
	}

	int getDiff(char * target, int length)
	{
	
	    int diff=0;
	    for(int i=0;i<length;i++)
	    {
	        if(toupper(target[i])!=concensus[i])
	        {
	            diff++;
	        }
	        
	    }
	    
	//    cout<<concensus<<endl;
	    
	    return diff;
	    
	//    cout<<diff<<endl;
	    
	    
	}

	void finalize(int gapLength)
    {

        FILE * scaffoldMap=fopen(gapFileName, "r");

        char *qname1,*qname2;

        char *line1= new char[MAX_REC_LEN];
        char *line2= new char[MAX_REC_LEN];

        int MAX_FILE_READ=MAX_REC_LEN/sizeof(line1[0]);


        char * temp;
        char *rname1, *rname2;
        int	pos1,pos2,flag1, flag2, strandNo1, strandNo2, insertSize1, insertSize2;
        char *cigar1, *cigar2, *readString1, *readString2, md1[1000], md2[1000];
        int ih1, ih2;
        int readLength1, readLength2;
        int fromCharCode, toCharCode, index, k;

        char tempReadString[100];

        double maxProb=0,tempProb,maxLikelihood=0;

        int insertSize,mleInsertSize,tempInsertSize;

        int maxMatch, tempMatch, maxPos;

        int maxSize=gapLength+2*maxDistance;

        int totalCount=0, discardedCount=0;

        char *gapString=new char[maxSize+1];

        for(long int i=0;i<maxDistance;i++)
        {
            gapString[i]=contigs[contigNo][i+startPos];

        }
        for(long int i=maxDistance;i<maxDistance+gapLength;i++)
        {
            gapString[i]=bestString[i-maxDistance];

        }
        for(long int i=maxDistance+gapLength;i<maxDistance+gapLength+maxDistance;i++)
        {
            gapString[i]=contigs[contigNo][i+startPos-gapLength+originalGap];

        }
        gapString[maxSize]='\0';
        //cout<<"GapString = \n"<< gapString<<endl;
        //cout<<"end pos = "<<endPos<<" , startpos = "<<startPos<<endl;
        for(long int i=maxDistance;i<endPos-startPos-maxDistance;i++)
        {
            for(int j=0;j<5;j++)
            {
                countsGap[i][j]=0;
            }
        }


        startPos=this->gapStart-maxDistance;
        endPos=this->gapStart+this->gapLength+maxDistance;

        this->gapLength=gapLength;

        while(fgets(line1, MAX_FILE_READ, scaffoldMap)!=NULL)
        {
            //cout<<"At the start of while"<<endl;
            if(line1[0]=='@')
                continue;
            //????
            if(fgets(line2, MAX_FILE_READ, scaffoldMap)==NULL)
                break;

            qname1=strtok(line1,"\t");

            temp=strtok(NULL,"\t");
            flag1=atoi(temp);


            strandNo1=(flag1&16)>>4;

            rname1=strtok(NULL,"\t");


            temp=strtok(NULL,"\t");
            pos1=atoi(temp);


            cigar1=strtok(NULL,"\t");


            temp=strtok(NULL,"\t");



            insertSize1=atoi(temp);


            readString1=strtok(NULL,"\t");




            //second of the pair

            qname2=strtok(line2,"\t");
            temp=strtok(NULL,"\t");
            flag2=atoi(temp);


            strandNo2=(flag2&16)>>4;

            rname2=strtok(NULL,"\t");

            temp=strtok(NULL,"\t");
            pos2=atoi(temp);



            cigar2=strtok(NULL,"\t");


            temp=strtok(NULL,"\t");

            insertSize2=atoi(temp);

            readString2=strtok(NULL,"\t");


            readLength1=strlen(readString1);
            readLength2=strlen(readString2);

            if(strandNo1==0)
            {
                reverse(tempReadString, readString2);
                strcpy(readString2,tempReadString);
                strandNo2=1;
            }
            else
            {
                strandNo2=0;
            }

            maxProb=0;
            maxMatch=0;
            maxPos=0;

            totalCount++;

            if(pos1<gapStart)
            {

                insertSize=gapStart-pos1+readLength2;
                for(long int i=gapStart-readLength2+1;i<gapStart+gapLength;i++)
                {

                    tempInsertSize=insertSize+i-gapStart;
                    if(tempInsertSize<insertThresholdMin || tempInsertSize > insertThresholdMax)
                        continue;


                    tempProb=insertLengthDistSmoothed[tempInsertSize];
                    tempMatch=0;

                    for(int j=0;j<readLength2;j++)
                    {

                        toCharCode=charCodes[readString2[j]];
                        index=i+j-startPos;
                        fromCharCode=charCodes[gapString[index]];

                        if(strandNo2==0)
                        {
                            k=j;
                        }
                        else
                        {
                            k=readLength2-j-1;
                        }

                        if(fromCharCode==toCharCode)
                        {
                            tempProb*=(1-errorPosDist[k]-inPosDist[k]-delPosDist[k]);
                        }
                        else
                        {
                            if(fromCharCode >=0 && fromCharCode<=4)tempProb*=errorPosDist[k]*errorTypeProbs[fromCharCode][toCharCode];
                        }

                    }
                    if(tempProb>maxProb)
               //     if(tempMatch>maxMatch)
                    {
                        maxMatch=tempMatch;
                        maxProb=tempProb;
                        mleInsertSize=insertSize+i-gapStart;
                        maxPos=i-startPos;

                    }


                }
               //cout<<"gap_cutoff = "<<-log(maxProb)<<endl;
                if(-log(maxProb)<=gapProbCutOff)
                {
                    for(int j=0;j<readLength2;j++)
                    {
                        if(maxPos+j>=maxDistance && maxPos+j<maxDistance+gapLength)
                        {
                            countsGap[maxPos+j][charCodes[readString2[j]]]+=1;
                        }
                    }
                }
                else
                {
                    discardedCount++;
                }
     //           cout<<(pos1-startPos)<<"\t"<<mleInsertSize<<"\t"<<maxProb<<"\t"<<maxMatch<<"\t"<<readString2<<endl;

            }
            else
            {
                insertSize=pos1-(gapStart+gapLength)+readLength1;

                for(long int i=gapStart-readLength2+1;i<gapStart+gapLength;i++)
                {
                    tempInsertSize=insertSize+gapStart+gapLength-i;
                    if(tempInsertSize<insertThresholdMin || tempInsertSize > insertThresholdMax)
                        continue;

                    //cout<<"After continue"<<endl;
                    tempProb=insertLengthDistSmoothed[tempInsertSize];
                    tempMatch=0;

                    for(int j=0;j<readLength2;j++)
                    {

                        toCharCode=charCodes[readString2[j]];
                        index=i+j-startPos;
                        fromCharCode=charCodes[gapString[index]];

                        if(strandNo2==0)
                        {
                            k=j;
                        }
                        else
                        {
                            k=readLength2-j-1;
                        }

                        //cout<<"Value of index = "<<index<<" , From = "<<fromCharCode<<", to = "<<toCharCode<<endl;
                        if(fromCharCode==toCharCode)
                        {
                            tempProb*=(1-errorPosDist[k]-inPosDist[k]-delPosDist[k]);
                            //cout<<"tempProb = "<<tempProb<<endl;
                        }
                        else
                        {
                            if(fromCharCode >=0 && fromCharCode<=4)tempProb*=errorPosDist[k]*errorTypeProbs[fromCharCode][toCharCode];
                            //cout<<"tempProb = "<<tempProb<<endl;
                        }
                    }
                    if(tempProb>maxProb)
                    {
                        maxMatch=tempMatch;
                        maxProb=tempProb;
                        mleInsertSize=insertSize+gapStart+gapLength-i;
                        maxPos=i-startPos;
                    }
                }
       //         cout<<(endPos-pos1+readLength1)<<"\t"<<mleInsertSize<<"\t"<<maxProb<<"\t"<<maxMatch<<"\t"<<readString2<<endl;
                //cout<<"Out offor loop"<<endl;
                //cout<<"else gap_cutoff = "<<-log(maxProb)<<endl;
                if(-log(maxProb)<=gapProbCutOff)
                {
                    for(int j=0;j<readLength2;j++)
                    {
                        if(maxPos+j>=maxDistance && maxPos+j<maxDistance+gapLength)
                        {
                            countsGap[maxPos+j][charCodes[readString2[j]]]+=1;
                        }
                    }
                }
                else
                {

                    discardedCount++;
                }
            }
            //cout<<"At the end of while"<<endl;
        }

        cout<<"Total_Reads = "<<totalCount<<"\t"<<"Discard = "<<discardedCount<<endl;

        delete [] line1;
        delete [] line2;

        computeSequence(1);

        delete [] gapString;

        fclose(scaffoldMap);

    }

	void printArr(int a,int b)
	{
	    /*printf("Printing probsgap array\n\n");
	    for(int i=a;i<b;i++)
	    {
	        cout<<"i = "<<i<<"==>";
	        for(int j=0;j<4;j++)
	            cout<<probsGap[i][j]<<" ";
	        cout<<endl;


	    }

*/
        printf("Printing countsgap array\n\n");
        for(int i=a;i<b;i++)
        {
            cout<<"i = "<<i<<"==>";
            for(int j=0;j<4;j++)
                cout<<countsGap[i][j]<<" ";
            cout<<endl;
        }
	    /*printf("Printing errorprobsgap array\n\n");
	    for(int i=a;i<b;i++)
        {
            cout<<"i = "<<i<<"==>";
            for(int j=0;j<4;j++)
                cout<<errorProbsGap[i][j]<<"    ";
            cout<<endl;
        }

        printf("Printing errorposdist array\n\n");
        for(int i=0;i<maxReadLength;i++)
        {
            cout<< errorPosDist[i]<<endl;
        }
        cout<<endl;

        */
	}

	void fillGap()
	{
	    float gp_frac,gp_frac2;
	    if(originalGap>300) gp_frac = .95;
	    else
	    {
	        if(originalGap<100)gp_frac=0.1;
	        else gp_frac = 0.85;
	     }

	    if(originalGap>200) gp_frac2 = 1.05;
        else
        {
            gp_frac2 = 1.5;
        }

		int gapMin=originalGap*gp_frac<10?1:originalGap*gp_frac;
		int gapMax=originalGap*gp_frac2;
		int gapEstimate=gapMin;
        int maxGapEstimate = gapMin;
        int secondMaxGapEstimate = gapMin;

        double maxLikelihood=-DBL_MAX;
        double secondMaxLikelihood=-DBL_MAX;

        double likelihood=0;

        bestString[0]='\0';
		secondBestString[0]='\0';

        int num_itr = 5;

        for(int j=0;j<(gapMax-gapMin+1);j++)
        {
            initialize(gapEstimate);

            for(int i=0;i<num_itr;i++)
            {
                likelihood = placeReads(i);

                computeProbsGap();
                computeErrorProbsGap();
            }
            if(likelihood>=maxLikelihood)
            {
                secondMaxLikelihood=maxLikelihood;
                secondMaxGapEstimate=maxGapEstimate;
                strcpy(secondBestString,bestString);
                maxLikelihood=likelihood;
                maxGapEstimate=gapEstimate;


                computeSequence(0);
                strcpy(bestString,concensus);
            }
            else if(likelihood>=secondMaxLikelihood)
            {
                secondMaxLikelihood=likelihood;
                secondMaxGapEstimate=gapEstimate;
                computeSequence(0);
                strcpy(secondBestString,concensus);
            }

            //cout<<"GapEstimate = "<<gapEstimate<<"\t"<<"Likelihood = "<<likelihood<<endl;
            gapEstimate++;
        }

		cout<<"----##############MaxGapEstimate = "<<maxGapEstimate<<"\t"<<"MaxLikelihood = "<<maxLikelihood<<endl;
		//cout<<"Best String = "<<bestString<<endl;
		//printArr(maxDistance,maxDistance+maxGapEstimate);

		finalize(maxGapEstimate);

		//printArr(maxDistance,maxDistance+maxGapEstimate);
	}
    
    void freeGapFiller()
        {
    		int maxSize=maxGap+2*maxDistance;

            for(long int i=0;i<maxSize;i++)
            {
                delete [] countsGap[i];
                delete [] probsGap[i];
                delete [] errorProbsGap[i];
            }

            delete [] countsGap;
            delete [] probsGap;
            delete [] errorProbsGap;

            delete [] concensus;

            delete [] bestString;
            delete [] secondBestString;
           }
};


int main(int argc, char *argv[])
{
	/*	input contig file name, read file name
     contig file - fasta format
     read file - fastq format
     */
    srand (SEED);
    //    srand (rand());

    //==========================================================================
    //Declaration of all variables

    char *line= new char[MAX_REC_LEN];
   	char *line1= new char[MAX_REC_LEN];
   	char *line2= new char[MAX_REC_LEN];

   	int read;
   	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);
   	long int bufferLength=1024;

   	char *contig=new char[bufferLength];
   	contig[0]='\0';
   	char *newcontig;
   	char *contigName;
   	contigLength=0;

    long int contigNo;
    double readCount;
    char *temp;

   	long int tempContigLength=0;

   	long int totalInserts=0;
    long int insertSum=0;

    int nonzeroCount=0;

	//======================Declarations complete====================================================

    //======================Start reading the ref genome with gaps file==============================

    //Open gapped contig file for reading
	contigFileName=argv[1];
	f = atoi(argv[2]);
	contigFile=fopen(contigFileName, "r");
	outFile=fopen("out.txt", "w");
	if (contigFile == NULL)
	{
		printf("Can't open contig file\n");
		exit(1);
	}
	while(fgets(line, MAX_FILE_READ, contigFile)!=NULL)
	{
		if(line[0]==';')
		{
			continue;
		}
		else if(line[0]=='>')
		{
			contigName=new char[strlen(line)];
			strcpy(contigName,line+1);
			contigName[strlen(contigName)-1]='\0';
			contigNames.push_back(strtok(contigName," \t\n"));

			if(contigLength>0)
			{
				noContigs++;
				contigs.push_back(contig);
				contigLengths.push_back(contigLength);
				totalContigLength+=contigLength;
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
				strcat(newcontig, line);
				delete []contig;
				contig=newcontig;
			}
			else
			{
				line[read-1]='\0';
				strcpy(contig+tempContigLength, line);
			}
		}
	}

	noContigs++;
	contigs.push_back(contig);
	contigLengths.push_back(contigLength);
	totalContigLength+=contigLength;
    
    fclose(contigFile);
    
    //=====================Finished with reading the contig file===========================
    
	for(long int i=0;i<noContigs;i++)
	{
		for(long int j=0;j<contigLengths[i];j++)
		{
			contigs[i][j]=toupper(contigs[i][j]);
		}
	}

	/*
     use bfast or some other tool to map reads and save mapping
     */
	
	mapFileName="myout.sam";
    
	mapFile=fopen("myout.sam", "r");
    
	summaryFile=fopen("stat.txt","r");
    
	if (mapFile == NULL)
	{
		printf("Can't open map file\n");
		exit(1);
	}
    
	int count=0;

	//doesn't match with main.cpp
	//fscanf(summaryFile,"%ld %ld %d %d %d",&totalCount, &unCount, &toAlign, &maxReadLength, &MAX_INSERT_SIZE);
	fscanf(summaryFile,"%ld %ld %d %d",&totalCount, &unCount, &maxReadLength, &MAX_INSERT_SIZE);

    //cout<<totalCount<<" "<<unCount<<" "<<" "<<maxReadLength<<" "<<MAX_INSERT_SIZE<<endl;
    MAX_INSERT_SIZE=MAX_INSERT_SIZE>20000?MAX_INSERT_SIZE:20000;

    insertCutoffMin=MAX_INSERT_SIZE;//20k

    //There are two variables,MAX_INSERT_SIZE and maxInsertSize,
    //at this fn, initInsertCounts, 2nd one is initialized by 1st one .i.e. 20k
	initInsertCounts(MAX_INSERT_SIZE);
    initErrorTypes(maxReadLength);
    
    
	myitoa(maxReadLength, noErrorCigar, 10);
	strcpy(noErrorMD,"MD:Z:");
	strcat(noErrorMD,noErrorCigar);
	strcat(noErrorCigar,"M");

    //So, noErrorMD = MD:Z:100                   noErrorCigar = 100M

    //============================Read myout.sam line by line and store those reads=========
	int noMatches=0;
    
	while(fgets(line1, MAX_FILE_READ, mapFile)!=NULL)
	{
		
		if(line1[0]=='@')
			continue;
        
		processMapping(line1);
		count+=1;
	}
    //cout<<"Myout.sam file has "<<count<<" lines"<<endl;
	fclose(mapFile);
	fclose(summaryFile);
    
    //===========================Done with myout.sam========================================

	computeProbabilites();
    double val1 = computeLikelihood(mapFileName);
   // cout<<"likelihood //val 1  = "<<val1<<endl;
    
    FILE * mappedInsertFile=fopen("mappedInserts.txt","w");
    
    for(long int i=0;i<maxInsertSize;i++)
    {
        if(insertCountsMapped[i]>0)
        {
            //    cout<<i<<","<<insertCountsAll[i]<<endl;
            fprintf(mappedInsertFile, "%ld,%ld\n",i,insertCountsMapped[i]);
            totalInserts+=insertCountsMapped[i];
            insertSum+=insertCountsMapped[i]*i;
            nonzeroCount++;
        }
    }

    double insertMean=insertSum/(double)totalInserts;
    //cout<<"Total Inserts = "<<totalInserts<<endl;
    
    fclose(mappedInsertFile);
    fclose(outFile);
    //	cout<<"after val 2"<<endl;
    //cout<<"insertCutoffMin = "<<insertCutoffMin<<" ...insertCutoffMax = "<<insertCutoffMax<<endl;
    
    cout<<"insertThresholdMin = "<<insertThresholdMin<<"...insertThresholdMax = "<<insertThresholdMax<<endl;
    
    //cout<<"insertSizeMode = "<<insertSizeMode<<endl;

    /*for(int i=0;i<maxInsertSize;i++)
    {
        cout<<i<<" "<<insertLengthDistSmoothed[i]<<endl;
    }*/

    for(int i=0;i<256;i++)
    {
        if(i=='A')
        {
            charCodes[i]=0;
        }
        else if(i=='C')
        {
            charCodes[i]=1;
        }
        else if(i=='G')
        {
            charCodes[i]=2;
        }
        else if(i=='T')
        {
            charCodes[i]=3;
        }
        else
        {
            charCodes[i]=4;
        }

    }

	long int gapProbSum=0;
	for(int i=0;i<1000;i++)
    {
    	gapProbSum+=gapProbs[i];
    }

	long int gapProbCount=0;
	
    for(int i=0;i<1000;i++)
    {
    	gapProbCount+=gapProbs[i];
    	if(gapProbCount>=0.99*gapProbSum)
    	{
    		gapProbCutOff=i;
    		break;
    	}
    }

    //cout<<"gap_prob_cutoff in main = "<<gapProbCutOff<<endl;
    time_t now, end;
    time(&now);

    //====================================Start filling the gaps from here===================
    int gapContigNo;
    long int gapStart;
    int gapLength;
    int gapStringLength;
    int diff;

    char gapFileName[500];
    char *gapString=new char[MAX_GAP];


    temp=new char[100];
    char *temp2=new char[100];

    int gapNo=0;

    gapInfoFile=fopen("gapInfo.txt","r");
    gapOutFile=fopen("gapout.txt","w");

    while(fgets(line1, MAX_FILE_READ, gapInfoFile)!=NULL)
    {
    	strcpy(gapFileName,"Gaps/gaps_");
    	//itoa(gapNo,temp,10);
    	sprintf(temp2,"%d",gapNo);

    	strcat(gapFileName,temp2);

    	strcat(gapFileName,".sam");

        
        temp=strtok(line1,"\t");
        gapContigNo=atoi(temp);

        temp=strtok(NULL,"\t");
        gapStart=atol(temp);

        temp=strtok(NULL,"\t\n");
        gapLength=atoi(temp);

        GapFiller gf;
        gf.allocate(gapLength*2);

		gf.initGapFiller(gapContigNo,gapStart,gapLength, gapFileName);

        gf.fillGap();

		gapStringLength = gf.getConcensus(gapString);

		fprintf(gapOutFile,"%d\t%d\t%ld\t%d\t%d\t%s\n",gapNo,gapContigNo,gapStart,gapLength,gapStringLength,gapString);

		gf.freeGapFiller();

		cout<<"Gap "<<gapNo<<" is done."<<endl;

		gapNo++;
    }
   
    fclose(gapInfoFile);
    fclose(gapOutFile);
    time(&end);
    
    double seconds = difftime(end,now);

    cout<<"Total time to fill gaps in seconds = "<<seconds<<endl;
    
    //==========================================Gap filling is complete=================================
    
    int nStart=0;
    long int nStartPos=0;
    int nCount=0,afterCount=0;

    long int offset=0, index=0;
    
    gapOutFile=fopen("gapout.txt","r");
    
    filledContigFile=fopen("filledContigs.fa","w");
    
    cout<<"# of Contig = "<<noContigs<<endl;

    for(int i=0;i<noContigs;i++)
    {
        //cout<<"Contig len = "<<contigLengths[i]<<endl;
    	contig=new char[contigLengths[i]+1];
    	//cout<<"Allocation Done"<<endl;
    	offset=0;
    	index=0;

    	fprintf(filledContigFile,">%s\n",contigNames[i]);

        for(long int j=0;j<contigLengths[i];j++)
        {
            if(contigs[i][j]=='N' || contigs[i][j]=='n')
            {
                //cout<<"found N"<<endl;
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
            else 
            {
            	if(nStart==1)
                {
                	if(nCount>=1)
                	{
                   		fscanf(gapOutFile,"%d\t%d\t%ld\t%d\t%d\t%s\n",&gapNo,&gapContigNo,&gapStart,&gapLength,&gapStringLength,gapString);
                    	//cout<<gapNo<<" "<<gapContigNo<<" "<<gapStart<<" "<<gapLength<<endl;
                    	contig[index]='\0'; 
						fprintf(filledContigFile,"%s\n",contig);
                    	fprintf(filledContigFile,"%s\n",gapString);
                    
                    	offset=offset+(gapStringLength-gapLength);
                    	contig=(char *)realloc(contig, (contigLengths[i]+1+offset-index)*sizeof(char));
                    	index=0;
                	}
                	else
                	{
                		contig[index]='\0';
                		fprintf(filledContigFile,"%s\n",contig);
                		index=0;
                		for(int k=0;k<nCount;k++)
                		{
                			contig[k]='N';
                		}
                		contig[nCount]='\0';
                		fprintf(filledContigFile,"%s\n",contig);
                	}
                	nStart=0;
                }
                contig[index]=contigs[i][j];
                index++;
            }
        }
        contig[index]='\0';
    	fprintf(filledContigFile,"%s\n",contig);
        //cout<<"Done with Contig"<<i<< "successfully"<<endl;

    }

	fclose(gapOutFile);
	fclose(filledContigFile);
    cout<<"Figbird.cpp ends successfully"<<endl;
    
	return 0;
}
