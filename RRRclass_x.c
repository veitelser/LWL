#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NEWTONITER 50
#define NEWTONEPS 1e-12

int datanodes,hiddennodes,classnodes,nodes,edges,maxindeg;

int *indeg,**inedge,**innode,*outdeg,**outedge;
double *xref,*xnew,**x,**xp1,**xp2;
double yref,ynew,**y,**yp1,**yp2;
double *wref,*wnew,**w,**wp1,*wp2;
double **b,*bp1,**bp2;

int *class;
double *classdist,*baddist,leastbad;
int ee;

int *perm,*rank,trainsize,testsize,batchsize,batchnum,netiter=0;

double beta,wnorm,ynorm,margin;

FILE *trainptr,*testptr;
fpos_t trainstart,teststart,*trainitem,*testitem;


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
	
fscanf(fp,"%d%d%d",&datanodes,&hiddennodes,&classnodes);

nodes=datanodes+hiddennodes+classnodes;

indeg=malloc(nodes*sizeof(int));
outdeg=malloc(nodes*sizeof(int));
innode=malloc(nodes*sizeof(int*));
inedge=malloc(nodes*sizeof(int*));
outedge=malloc(nodes*sizeof(int*));

for(i=0;i<nodes-classnodes;++i)
	outdeg[i]=0;

maxindeg=0;	
edges=0;
for(i=datanodes;i<nodes;++i)
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

for(i=0;i<nodes-classnodes;++i)
	{
	outedge[i]=malloc(outdeg[i]*sizeof(int));
	outdeg[i]=0;
	}

for(i=datanodes;i<nodes;++i)
for(d=0;d<indeg[i];++d)
	{
	in=innode[i][d];
	outedge[in][outdeg[in]++]=inedge[i][d];
	}

return 1;
}


int opendata(char *trainfile, char *testfile)
{
int datanodes1,datanodes2,classnodes1,classnodes2;

trainptr=fopen(trainfile,"r");
if(!trainptr)
	{
	printf("trainfile not found\n");
	return 0;
	}
	
testptr=fopen(testfile,"r");
if(!testptr)
	{
	printf("testfile not found\n");
	return 0;
	}
		
fscanf(trainptr,"%d",&datanodes1);
fscanf(testptr,"%d",&datanodes2);

fscanf(trainptr,"%d",&classnodes1);
fscanf(testptr,"%d",&classnodes2);
	
if(datanodes!=datanodes1 || datanodes!=datanodes2)
	{
	printf("trainfile or testfile is not compatible with netfile\n");
	return 0;
	}
	
if(classnodes1!=classnodes2)
	{
	printf("trainfile and testfile differ in class number\n");
	return 0;
	}

fgetpos(trainptr,&trainstart);
fgetpos(testptr,&teststart);
	
return 1;
}


void closedata()
{
fclose(trainptr);
fclose(testptr);
}


void setup()
{
int k;

xref=malloc(maxindeg*sizeof(double));
xnew=malloc(maxindeg*sizeof(double));
wref=malloc(maxindeg*sizeof(double));
wnew=malloc(maxindeg*sizeof(double));

x=malloc(batchsize*sizeof(double*));
xp1=malloc(batchsize*sizeof(double*));
xp2=malloc(batchsize*sizeof(double*));

y=malloc(batchsize*sizeof(double*));
yp1=malloc(batchsize*sizeof(double*));
yp2=malloc(batchsize*sizeof(double*));

w=malloc(batchsize*sizeof(double*));
wp1=malloc(batchsize*sizeof(double*));
wp2=malloc(edges*sizeof(double));

b=malloc(batchsize*sizeof(double*));
bp1=malloc(nodes*sizeof(double));
bp2=malloc(batchsize*sizeof(double*));

for(k=0;k<batchsize;++k)
	{
	x[k]=malloc(edges*sizeof(double));
	xp1[k]=malloc(edges*sizeof(double));
	xp2[k]=malloc(nodes*sizeof(double));
	
	y[k]=malloc(nodes*sizeof(double));
	yp1[k]=malloc(nodes*sizeof(double));
	yp2[k]=malloc(nodes*sizeof(double));
	
	w[k]=malloc(edges*sizeof(double));
	wp1[k]=malloc(edges*sizeof(double));
	
	b[k]=malloc(nodes*sizeof(double));
	bp2[k]=malloc(nodes*sizeof(double));
	}
	
class=malloc(batchsize*sizeof(int));

classdist=malloc(batchsize*sizeof(double));
baddist=malloc(ee*sizeof(double));
}


int readdata(int k,FILE *dataptr)
{
int i;

for(i=0;i<datanodes;++i)
	if(fscanf(dataptr,"%lf",&xp2[k][i])==EOF)
		return 0;
		
fscanf(dataptr,"%d",&class[k]);

return 1;
}


void initdata()
{
int itemcount;

trainsize=0;
while(readdata(0,trainptr))
	++trainsize;
	
trainitem=malloc(trainsize*sizeof(fpos_t));

perm=malloc(trainsize*sizeof(int));
rank=malloc(trainsize*sizeof(int));

fsetpos(trainptr,&trainstart);

itemcount=0;
do
	{
	fgetpos(trainptr,&trainitem[itemcount++]);
	}
while(readdata(0,trainptr));


testsize=0;
while(readdata(0,testptr))
	++testsize;
	
testitem=malloc(testsize*sizeof(fpos_t));

fsetpos(testptr,&teststart);

testsize=0;
do
	{
	fgetpos(testptr,&testitem[testsize++]);
	}
while(readdata(0,testptr));
}


double act(double z)
{
return z>0. ? z : 0.;
}


void evaldata(int k)
{
int i,d;

for(i=datanodes;i<nodes;++i)
	{
	yp2[k][i]=0.;
	for(d=0;d<indeg[i];++d)
		yp2[k][i]+=xp2[k][innode[i][d]]*wp2[inedge[i][d]];

	yp2[k][i]/=wnorm;
	
	xp2[k][i]=yp2[k][i]-bp1[i];
	
	if(i<datanodes+hiddennodes)
		xp2[k][i]=act(xp2[k][i]);
	}	
}


double urand()
{
return ((double)rand())/RAND_MAX;
}


void init(int item)
{
int k,e,i,d;

for(k=0;k<batchsize;++k)
	{
	fsetpos(trainptr,&trainitem[perm[item+k]]);
	
	readdata(k,trainptr);
		
	for(e=0;e<edges;++e)
		w[k][e]=wp2[e];
		
	for(i=datanodes;i<nodes;++i)
		b[k][i]=bp1[i];
		
	evaldata(k);
	
	for(i=datanodes;i<nodes;++i)
		{
		for(d=0;d<indeg[i];++d)
			x[k][inedge[i][d]]=xp2[k][innode[i][d]];
			
		y[k][i]=yp2[k][i];
		}
	}
}


double func(double g,double pdot,double qdot,double t)
{
double den;

den=1.-t*t;
return wnorm*wnorm*t/g+(pdot*(1.+t*t)+qdot*t)/(den*den)-wnorm*yref;
}


double dfunc(double g,double pdot,double qdot,double t)
{
double den;

den=1.-t*t;
return wnorm*wnorm/g+(2.*pdot*t+qdot)/(den*den)+4.*t*(pdot*(1.+t*t)+qdot*t)/(den*den*den);
}


void bilinproj(int n)
{
int j,iter;
double pdot,qdot,g,t,ta,tb,tc,f,df,den;

if(n<datanodes+hiddennodes)
	g=outdeg[n];
else
	g=ynorm;

pdot=0.;
qdot=0.;
for(j=0;j<indeg[n];++j)
	{
	pdot+=wref[j]*xref[j];
	qdot+=wref[j]*wref[j]+xref[j]*xref[j];
	}
	
ta=0.;
f=func(g,pdot,qdot,ta);
tb=f>0. ? -1. : 1.;

for(iter=1;iter<=NEWTONITER;++iter)
	{
	df=dfunc(g,pdot,qdot,ta);
	
	t=ta-f/df;
	
	tc=.5*(ta+tb);
	
	if(tb>ta)
		{
		if(t>tc)
			t=tc;
		
		f=func(g,pdot,qdot,t);
		if(f>0.)
			tb=ta;
		}
	else
		{
		if(t<tc)
			t=tc;
		
		f=func(g,pdot,qdot,t);
		if(f<0.)
			tb=ta;
		}
	
	if(fabs(t-ta)<NEWTONEPS)
		break;
		
	ta=t;
	}

ynew=yref-wnorm*t/g;

den=1.-t*t;

for(j=0;j<indeg[n];++j)
	{	
	wnew[j]=(wref[j]+t*xref[j])/den;
	xnew[j]=(xref[j]+t*wref[j])/den;
	}
}


void proj1()
{
int k,i,d,in;

for(k=0;k<batchsize;++k)
for(i=datanodes;i<nodes;++i)
	{
	yref=2.*yp2[k][i]-y[k][i];
	
	for(d=0;d<indeg[i];++d)
		{
		in=inedge[i][d];
		wref[d]=2.*wp2[in]-w[k][in];
		xref[d]=2.*xp2[k][innode[i][d]]-x[k][in];
		}
		
	bilinproj(i);
	
	yp1[k][i]=ynew;
	
	for(d=0;d<indeg[i];++d)
		{
		in=inedge[i][d];
		wp1[k][in]=wnew[d];
		xp1[k][in]=xnew[d];
		}
	}
	
for(i=datanodes;i<nodes;++i)
	{
	bp1[i]=0.;
	for(k=0;k<batchsize;++k)
		bp1[i]+=2.*bp2[k][i]-b[k][i];
		
	bp1[i]/=batchsize;
	}
}


void actproj(double *xa,double *ya,double *ba)
{
double s,x1,y1,b1,d1,x2,y2,b2,d2;

x1=0.;
d1=(*xa)*(*xa);

if(*ya<*ba)
	{
	y1=*ya;
	b1=*ba;
	}
else
	{
	s=(*ya-*ba)/2.;
	y1=*ya-s;
	b1=y1;
	d1+=2.*s*s;
	}

s=(*xa-*ya+*ba)/3.;

x2=*xa-s;

if(x2<0.)
	{
	x2=0.;
	
	s=(*ya-*ba)/2.;
	y2=*ya-s;
	b2=y2;
	d2=(*xa)*(*xa)+2.*s*s;
	}
else
	{
	y2=*ya+s;
	b2=*ba-s;
	d2=3.*s*s;
	}
	
if(d1<d2)
	{
	*xa=x1;
	*ya=y1;
	*ba=b1;
	}
	else
	{
	*xa=x2;
	*ya=y2;
	*ba=b2;
	}
}


void classproj()
{
int k,i,j,e,leastbadpos;
double s;

for(k=0;k<batchsize;++k)
	{
	classdist[k]=0.;
	
	for(j=0;j<classnodes;++j)
		{
		i=nodes-classnodes+j;
		
		yp2[k][i]=y[k][i];
		bp2[k][i]=b[k][i];
				
		if(class[k]==j)
			{
			s=(y[k][i]-b[k][i]-margin)/2.;
			if(s<0.)
				{
				yp2[k][i]-=s;
				bp2[k][i]+=s;
	
				classdist[k]+=2.*s*s;
				}
			}
		else
			{
			s=(y[k][i]-b[k][i])/2.;
			if(s>0.)
				{
				yp2[k][i]-=s;
				bp2[k][i]+=s;
	
				classdist[k]+=2.*s*s;
				}
			}
		}
	}

if(ee==0)
	{
	leastbad=1e10;
	return;
	}
	
baddist[0]=classdist[0];
leastbad=baddist[0];
leastbadpos=0;

for(e=1;e<ee;++e)
	{
	baddist[e]=classdist[e];
	if(baddist[e]<leastbad)
		{
		leastbad=baddist[e];
		leastbadpos=e;
		}
	}
	
for(k=ee;k<batchsize;++k)
	if(classdist[k]>leastbad)
		{
		baddist[leastbadpos]=classdist[k];
	
		leastbad=baddist[0];
		leastbadpos=0;
	
		for(e=1;e<ee;++e)
			if(baddist[e]<leastbad)
				{
				leastbad=baddist[e];
				leastbadpos=e;
				}
		}
		
for(k=0;k<batchsize;++k)
	if(classdist[k]>=leastbad)
		{
		for(j=0;j<classnodes;++j)
			{
			i=nodes-classnodes+j;
			
			yp2[k][i]=y[k][i];
			bp2[k][i]=b[k][i];
			}
		}
}



void proj2()
{
int i,k,d,in;
double norm;

for(k=0;k<batchsize;++k)
for(i=datanodes;i<nodes-classnodes;++i)
	{
	xp2[k][i]=0.;
	for(d=0;d<outdeg[i];++d)
		xp2[k][i]+=x[k][outedge[i][d]];
			
	xp2[k][i]/=outdeg[i];
	yp2[k][i]=y[k][i];
	bp2[k][i]=b[k][i];
		
	actproj(&xp2[k][i],&yp2[k][i],&bp2[k][i]);
	}
	
classproj();
	
for(i=datanodes;i<nodes;++i)
	{
	norm=0.;
	for(d=0;d<indeg[i];++d)
		{
		in=inedge[i][d];
		wp2[in]=0.;
		for(k=0;k<batchsize;++k)
			wp2[in]+=w[k][in];
			
		wp2[in]/=batchsize;
		
		norm+=wp2[in]*wp2[in];
		}
		
	norm=sqrt(norm)/wnorm;
	
	for(d=0;d<indeg[i];++d)
		{
		in=inedge[i][d];
		wp2[in]/=norm;
		}
	}
}


double RRR()
{
int k,i,d,e;
double g,err,diff;

err=0.;

proj2();
proj1();

for(k=0;k<batchsize;++k)
	{
	for(i=datanodes;i<nodes;++i)
		{
		for(d=0;d<indeg[i];++d)
			{
			diff=xp1[k][inedge[i][d]]-xp2[k][innode[i][d]];
			x[k][inedge[i][d]]+=beta*diff;
			
			err+=diff*diff;
			}
			
		if(i<datanodes+hiddennodes)
			g=outdeg[i];
		else
			g=ynorm;
			
		diff=yp1[k][i]-yp2[k][i];
		y[k][i]+=beta*diff;
	
		err+=g*diff*diff;
		
		diff=bp1[i]-bp2[k][i];
		b[k][i]+=beta*diff;
	
		err+=g*diff*diff;
		}
	
	for(e=0;e<edges;++e)
		{
		diff=wp1[k][e]-wp2[e];
		w[k][e]+=beta*diff;
		
		err+=diff*diff;
		}
	}

return sqrt(err/batchsize);
}


void randparam()
{
int i,in,d;
double ave,norm;

for(i=datanodes;i<nodes;++i)
	{
	ave=0.;
	for(d=0;d<indeg[i];++d)
		{
		in=inedge[i][d];
		wp2[in]=2.*urand()-1.;
		ave+=wp2[in];
		}
	
	ave/=indeg[i];
		
	norm=0.;
	for(d=0;d<indeg[i];++d)
		{
		in=inedge[i][d];
		wp2[in]-=ave;
		norm+=wp2[in]*wp2[in];
		}
		
	norm=sqrt(norm)/wnorm;
	
	for(d=0;d<indeg[i];++d)
		{
		in=inedge[i][d];
		wp2[in]/=norm;
		}
		
	bp1[i]=0.;
	}
}


int classerr(int k)
{
int maxclass,c;
double xmax;

evaldata(k);

xmax=xp2[k][nodes-classnodes];
maxclass=0;
	
for(c=1;c<classnodes;++c)
	if(xp2[k][nodes-classnodes+c]>xmax)
		{
		xmax=xp2[k][nodes-classnodes+c];
		maxclass=c;
		}
	
if(maxclass!=class[k])
	return 1;
else
	return 0;
}
	
	
void printbad(char *badfile)
{
int k,i;
FILE *fp;

fp=fopen(badfile,"a");

for(k=0;k<batchsize;++k)
	if(classdist[k]>=leastbad)
		{
		for(i=0;i<datanodes;++i)
			fprintf(fp,"%lf ",xp2[k][i]);
			
		fprintf(fp,"%d\n",class[k]);
		}
		
fclose(fp);
}
	

void train(int itermax,double *aveiter,double tol,double *errstart,double *errstop,double *trainerr,char *badfile)
{
int batchcount,iter,k,samples;
double err;

*errstart=0.;
*errstop=0.;
*aveiter=0.;
*trainerr=0.;

samples=0;
for(batchcount=0;batchcount<batchnum;++batchcount)
	{
	init(batchsize*batchcount);
	
	err=RRR();
	iter=1;
	
	*errstart+=err;
	
	while(iter<itermax && err>tol)
		{
		err=RRR();
		++iter;
		}
			
	*errstop+=err;
	
	*aveiter+=iter;
	
	netiter+=iter;
	
	for(k=0;k<batchsize;++k)
		if(classdist[k]<leastbad)
			{
			*trainerr+=classerr(k);
			++samples;
			}
		
	if(ee)
		printbad(badfile);
	}
	
*errstart/=batchnum;
*errstop/=batchnum;
*aveiter/=batchnum;
*trainerr/=samples;
}


double testerr(FILE *dataptr,fpos_t *datastart,int datasize)
{
int k,err;

err=0;

fsetpos(dataptr,datastart);

for(k=0;k<datasize;++k)
	{
	readdata(0,dataptr);
	err+=classerr(0);
	}

return ((double)err)/datasize;
}


void printparam(char *paramfile)
{
FILE *fp;
int i,d;

fp=fopen(paramfile,"w");

edges=0;
for(i=datanodes;i<nodes;++i)
	{
	fprintf(fp,"%lf\n",bp1[i]);
	
	for(d=0;d<indeg[i];++d)
		fprintf(fp,"%lf ",wp2[edges++]);
		
	fprintf(fp,"\n\n");
	}
	
fclose(fp);
}


int compare(const void *a,const void *b)
{
return rank[*(int*)a] - rank[*(int*)b];
}


void shuffle()
{
int k;

for(k=0;k<trainsize;++k)
	{
	rank[k]=rand();
	perm[k]=k;
	}
	
qsort(perm,trainsize,sizeof(int),compare);
}


double work()
{
return 1e-9*(double)edges*(double)batchsize*(double)netiter;
}


int main(int argc,char* argv[])
{
char *trainfile,*testfile,*netfile,*id,logfile[50],paramfile[50],badfile[50];
int c,e,epochs,itermax;
double tol,aveiter,errstart,errstop,trainerr;
FILE *fp;

if(argc==14)
	{
	trainfile=argv[1];
	testfile=argv[2];
	netfile=argv[3];
	
	batchsize=atoi(argv[4]);
	ee=atoi(argv[5]);
	epochs=atoi(argv[6]);
	
	itermax=atoi(argv[7]);
	tol=atof(argv[8]);
	
	beta=atof(argv[9]);
	wnorm=atof(argv[10]);
	ynorm=atof(argv[11]);
	margin=atof(argv[12]);
	
	id=argv[13];
	}
else
	{
	fprintf(stderr,"expected 13 arguments: trainfile, testfile, netfile, batchsize, ee, epochs, itermax, tol, beta, wnorm, ynorm, margin, id\n");
	return 1;
	}
	
if(!getnet(netfile))
	return 1;
	
if(!opendata(trainfile,testfile))
	return 1;

setup();

sprintf(logfile,"%s.log",id);
sprintf(paramfile,"%s.dat",id);
sprintf(badfile,"%s.bad",id);

srand(time(NULL));
randparam();

initdata();

batchnum=trainsize/batchsize;

fp=fopen(logfile,"w");

for(c=0;c<=12;++c)
	fprintf(fp,"%s ",argv[c]);
fprintf(fp,"\n\n");

fprintf(fp,"epoch      GWMs   aveiter  errstart   errstop  batcherr  trainerr   testerr\n\n");
fclose(fp);

for(e=1;e<=epochs;++e)
	{
	if(ee)
		{
		fp=fopen(badfile,"w");
		fclose(fp);
		}
		
	shuffle();
	train(itermax,&aveiter,tol,&errstart,&errstop,&trainerr,badfile);
	
	fp=fopen(logfile,"a");
	fprintf(fp,"%5d%10.2f%10.2e%10.2e%10.2e%10.6f%10.6f%10.6f\n",e,work(),aveiter,errstart,errstop,trainerr,testerr(trainptr,&trainstart,trainsize),testerr(testptr,&teststart,testsize));
	fclose(fp);
	
	if(errstart<tol)
		break;
	}
	
//printparam(paramfile);

closedata();

return 0;
}