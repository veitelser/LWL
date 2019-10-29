#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int datanodes,encodenodes,codenodes,decodenodes,nodes,edges,maxindeg;

int *indeg,**inedge,**innode,*outdeg,**outedge;
double *x,*y,*w,*b;

double **codeval;
int codevalsamples;

double wnorm,margin;
int zig;

FILE *dataptr;
fpos_t datastart;


double urand()
{
return ((double)rand())/RAND_MAX;
}


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
	
if(datanodes!=datanodes1)
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
int i;

x=malloc(nodes*sizeof(double));
y=malloc(nodes*sizeof(double));

w=malloc(edges*sizeof(double));
b=malloc(nodes*sizeof(double));
	
codeval=malloc(codenodes*sizeof(double*));
for(i=0;i<codenodes;++i)
	codeval[i]=malloc(codevalsamples*sizeof(double));
}


void readdata()
{
int i;

for(i=0;i<datanodes;++i)
	if(fscanf(dataptr,"%lf",&x[i])==EOF)
		{
		fsetpos(dataptr,&datastart);
		i=0;
		}
		
fscanf(dataptr,"%*d");
}


double act(double z)
{
if(zig && fabs(z)<margin)
	return (z+margin)/(2.*margin);
else
	return z<=0. ? 0. : 1.;
}


void encode()
{
int i,d;

for(i=datanodes;i<datanodes+encodenodes+codenodes;++i)
	{
	y[i]=0.;
	for(d=0;d<indeg[i];++d)
		y[i]+=x[innode[i][d]]*w[inedge[i][d]];

	y[i]/=wnorm;
	
	x[i]=act(y[i]-b[i]);
	}
}


void makecode()
{
int k,i;

for(k=0;k<codevalsamples;++k)
	{
	readdata();
	encode();
	
	for(i=0;i<codenodes;++i)
		codeval[i][k]=x[datanodes+encodenodes+i];
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


int main(int argc,char* argv[])
{
char *id,*datafile,*netfile,*paramfile,codefile[50];
int s,codesamples,i;
double truefrac;
FILE *fp;

if(argc==8)
	{
	datafile=argv[1];
	netfile=argv[2];
	paramfile=argv[3];
	
	codevalsamples=atoi(argv[4]);
	
	codesamples=atoi(argv[5]);
	truefrac=atof(argv[6]);
	
	id=argv[7];
	}
else
	{
	fprintf(stderr,"expected 7 arguments: datafile, netfile, paramfile, codevalsamples, codesamples, truefrac, id\n");
	return 1;
	}

if(!getnet(netfile))
	return 0;
	
if(!opendata(datafile))
	return 0;
	
setup();

if(!readparam(paramfile))
	return 0;

sprintf(codefile,"%s.encode",id);

srand(time(NULL));

makecode();

fp=fopen(codefile,"w");

fprintf(fp,"%d 2\n\n",codenodes);

for(s=0;s<codesamples;++s)
	if(urand()<truefrac)
		{
		readdata();
		encode();
	
		for(i=datanodes+encodenodes;i<datanodes+encodenodes+codenodes;++i)
			fprintf(fp,"%11.8lf",x[i]);
			
		fprintf(fp," 1\n");
		}
	else
		{
		for(i=0;i<codenodes;++i)
			fprintf(fp,"%11.8lf",codeval[i][rand()%codevalsamples]);
		
		fprintf(fp," 0\n");
		}
		
fclose(fp);

closedata();

return 1;
}
		
		
		
		
		
	
	
	
	
