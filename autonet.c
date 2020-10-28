#include <stdio.h>
#include <stdlib.h>

#define D 20 //maximum depth

int main(int argc,char* argv[])
{
int l,e,d,size[D],sgnsize,abssize,layers,datanodes,encodenodes,codenodes,decodenodes,nodes,in0,out0,i,j;
char *netfile;
FILE *fp;


if(argc==1)
	{
	printf("expected list of negative integers (encoder layer sizes), list of positive integers (decoder layer sizes), netfile name\n");
	
	return 1;
	}
	
encodenodes=0;
decodenodes=0;
e=0;
d=0;

for(l=0;l<argc-2;++l)
	{
	sgnsize=atoi(argv[l+1]);
	abssize=abs(sgnsize);
	
	size[l]=abssize;
	
	if(sgnsize<0)
		{
		if(e++==0)
			datanodes=abssize;
		else
			encodenodes+=abssize;
		}
	else
		{
		if(d++==0)
			codenodes=abssize;
		else
			decodenodes+=abssize;
		}
	}
	
layers=e+d;
nodes=datanodes+encodenodes+codenodes+decodenodes;

netfile=argv[argc-1];

fp=fopen(netfile,"w");	
fprintf(fp,"%d %d\n%d %d\n\n",datanodes,encodenodes,codenodes,decodenodes);

in0=nodes-size[layers-1];
out0=0;

for(j=0;j<size[0];++j)
	{
	fprintf(fp,"%d %d\n",out0+j,size[layers-1]);
		
	for(i=0;i<size[layers-1];++i)
		fprintf(fp,"%d ",in0+i);
			
	fprintf(fp,"\n");
	}

in0=0;
out0=size[0];

for(l=1;l<layers;++l)
	{
	for(j=0;j<size[l];++j)
		{
		fprintf(fp,"%d %d\n",out0+j,size[l-1]);
		
		for(i=0;i<size[l-1];++i)
			fprintf(fp,"%d ",in0+i);
			
		fprintf(fp,"\n");
		}
		
	in0+=size[l-1];
	out0+=size[l];
	}
	
fclose(fp);

return 0;
}