#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int datanodes,encodenodes,codenodes,decodenodes,nodes,edges,maxindeg;

int *indeg,**inedge,**innode,*outdeg,**outedge;
double *x,*y,*w,*b;

double wnorm,margin;
int zig;

FILE *dataptr;
fpos_t datastart;


int getnet(char *netfile)
{
int i,d,in;
FILE *fp;

fp=fopen(netfile,"r");
if(!fp)
	{
	printf("netfile not found\n");
	return 0;
	}
	
fscanf(fp,"%d%d%d%d",&datanodes,&encodenodes,&codenodes,&decodenodes);
	
nodes=datanodes+encodenodes+codenodes+decodenodes;

indeg=malloc(nodes*sizeof(int));
outdeg=malloc(nodes*sizeof(int));
innode=malloc(nodes*sizeof(int*));
inedge=malloc(nodes*sizeof(int*));
outedge=malloc(nodes*sizeof(int*));

for(i=0;i<nodes;++i)
	outdeg[i]=0;

maxindeg=0;	
edges=0;
for(i=0;i<nodes;++i)
	{
	fscanf(fp,"%*d%d",&indeg[i]);
	
	if(indeg[i]>maxindeg)
		maxindeg=indeg[i];
	
	innode[i]=malloc(indeg[i]*sizeof(int));
	inedge[i]=malloc(indeg[i]*sizeof(int));
	
	for(d=0;d<indeg[i];++d)
		{
		inedge[i][d]=edges++;
		
		fscanf(fp,"%d",&innode[i][d]);
		
		++outdeg[innode[i][d]];
		}
	}
	
fclose(fp);

for(i=0;i<nodes;++i)
	{
	outedge[i]=malloc(outdeg[i]*sizeof(int));
	outdeg[i]=0;
	}

for(i=0;i<nodes;++i)
for(d=0;d<indeg[i];++d)
	{
	in=innode[i][d];
	outedge[in][outdeg[in]++]=inedge[i][d];
	}

return 1;
}


int opendata(char *datafile)
{
int datanodes1;

dataptr=fopen(datafile,"r");
if(!dataptr)
	{
	printf("datafile not found\n");
	return 0;
	}
		
fscanf(dataptr,"%d%*d",&datanodes1);
	
if(codenodes!=datanodes1)
	{
	printf("datafile is not compatible with netfile\n");
	return 0;
	}

fgetpos(dataptr,&datastart);
	
return 1;
}


void closedata()
{
fclose(dataptr);
}


void setup()
{
x=malloc(nodes*sizeof(double));
y=malloc(nodes*sizeof(double));

w=malloc(edges*sizeof(double));
b=malloc(nodes*sizeof(double));
}


int readdata(int *class)
{
int i;

for(i=datanodes+encodenodes;i<datanodes+encodenodes+codenodes;++i)
	{
	if(fscanf(dataptr,"%lf",&x[i])==EOF)
		return 0;
	}
		
fscanf(dataptr,"%d",class);

return 1;
}


double act(double z)
{
if(zig && fabs(z)<margin)
	return (z+margin)/(2.*margin);
else
	return z<=0. ? 0. : 1.;
}


void decode()
{
int i,j,d,start,stop;

start=datanodes+encodenodes+codenodes;
stop=start+decodenodes+datanodes;

for(i=start;i<stop;++i)
	{
	j=i%nodes;
	
	y[j]=0.;
	for(d=0;d<indeg[j];++d)
		y[j]+=x[innode[j][d]]*w[inedge[j][d]];

	y[j]/=wnorm;
	
	x[j]=act(y[j]-b[j]);
	}
}


int readparam(char *paramfile)
{
FILE *fp;
int e,i,d;

fp=fopen(paramfile,"r");

fscanf(fp,"%lf%lf%d",&wnorm,&margin,&zig);

e=0;
for(i=0;i<nodes;++i)
	{
	fscanf(fp,"%lf",&b[i]);
	
	for(d=0;d<indeg[i];++d)
		fscanf(fp,"%lf",&w[e++]);
	}

if(e==edges)
	return 1;
else
	{
	fprintf(stderr,"netfile and paramfile are incompatible\n");
	return 0;
	}
}


void printparam()
{
FILE *fp;
int i,d;

fp=fopen("paramcheck","w");

fprintf(fp,"%lf %lf %d\n\n",wnorm,margin,zig);

edges=0;
for(i=0;i<nodes;++i)
	{
	fprintf(fp,"%lf\n",b[i]);
	
	for(d=0;d<indeg[i];++d)
		fprintf(fp,"%lf ",w[edges++]);
		
	fprintf(fp,"\n\n");
	}

fclose(fp);
}


int main(int argc,char* argv[])
{
char *id,*datafile,*netfile,*paramfile,decodefile[50];
int i,class;
FILE *fp;

if(argc==5)
	{
	datafile=argv[1];
	netfile=argv[2];
	paramfile=argv[3];
	
	id=argv[4];
	}
else
	{
	fprintf(stderr,"expected 4 arguments: datafile, netfile, paramfile, id\n");
	return 1;
	}

if(!getnet(netfile))
	return 0;
	
if(!opendata(datafile))
	return 0;
	
setup();

if(!readparam(paramfile))
	return 0;

printparam();

sprintf(decodefile,"%s.decode",id);

fp=fopen(decodefile,"w");

fsetpos(dataptr,&datastart);

while(readdata(&class))
	{
	decode();
	
	for(i=0;i<datanodes;++i)
		fprintf(fp,"%11.8lf",x[i]);
		
	fprintf(fp," %d\n",class);
	}
		
fclose(fp);

closedata();

return 1;
}
		
		
		
		
		
	
	
	
	
