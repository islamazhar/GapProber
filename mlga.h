#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*
void itoa(int value, char *s, int base)
{
	sprintf(s,"%d",value);
}
*/
char * toupper(char *contig)
{
	long int i=0;
	char c;
	while (contig[i])
  	{
    		c=contig[i];
    		contig[i]=toupper(c);
    		i++;
  	}
	return contig;
}


