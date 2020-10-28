#include <stdio.h>
#include <stdlib.h>

#define D 20 //maximum depth

int main(int argc,char* argv[])
{
int d,size[D],layers,hidden,instart,outstart,i,j;
char *netfile;
FILE *fp;

if(argc==1)
	{
	printf("expected list of two or more integers (layer sizes) followed by netfile name\n");
	
	return 1;
	}
	
for(d=0;d<argc-2;++d)
	size[d]=atoi(argv[d+1]);

layers=argc-3;

netfile=argv[argc-1];

hidden=0;
for(d=1;d<layers;++d)
	hidden+=size[d];

fp=fopen(netfile,"w");	
fprintf(fp,"%d %d %d\n\n",size[0],hidden,size[layers]);

instart=0;
outstart=size[0];
for(d=1;d<layers+1;++d)
	{
	for(j=0;j<size[d];++j)
		{
		fprintf(fp,"%d %d\n",outstart+j,size[d-1]);
		
		for(i=0;i<size[d-1];++i)
			fprintf(fp,"%d ",instart+i);
			
		fprintf(fp,"\n");
		}
		
	instart+=size[d-1];
	outstart+=size[d];
	}
	
fclose(fp);

return 0;
}
