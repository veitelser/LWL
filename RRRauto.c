#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NEWTONITER 50
#define NEWTONEPS 1e-12

int datanodes,encodenodes,codenodes,decodenodes,nodes,edges,maxindeg;

int *indeg,**inedge,**innode,*outdeg,**outedge;
double *xref,*xnew,**x,**xp1,**xp2;
double yref,ynew,**y,**yp1,**yp2;
double *wref,*wnew,**w,**wp1,*wp2;
double **b,*bp1,**bp2;

int *perm,*rank,datasize,databatch,codebatch,batchsize,batchnum,netiter=0;
int decodesamples;
double **codeval;

double beta,wnorm,margin;
int zig,exhaust;

FILE *dataptr;
fpos_t datastart,*dataitem;

FILE *errptr;


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
int k,i;

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
	
codeval=malloc(codenodes*sizeof(double*));
for(i=0;i<codenodes;++i)
	codeval[i]=malloc(databatch*sizeof(double));
}


int readdata(int k)
{
int i;

for(i=0;i<datanodes;++i)
	if(fscanf(dataptr,"%lf",&xp2[k][i])==EOF)
		return 0;
		
fscanf(dataptr,"%*d");

return 1;
}


void initdata()
{
int itemcount;

datasize=0;
while(readdata(0))
	++datasize;

dataitem=malloc(datasize*sizeof(fpos_t));

perm=malloc(datasize*sizeof(int));
rank=malloc(datasize*sizeof(int));

fsetpos(dataptr,&datastart);

itemcount=0;
do
	{
	fgetpos(dataptr,&dataitem[itemcount++]);
	}
while(readdata(0));
}


void exhaustcode(int k,int c)
{
int i;

for(i=0;i<codenodes;++i)
	{
	xp2[k][datanodes+encodenodes+i]=c%2;
	c/=2;
	}
}


void randcode(int k)
{
int i;

for(i=0;i<codenodes;++i)
	xp2[k][datanodes+encodenodes+i]=codeval[i][rand()%databatch];
}


double act(double z)
{
if(zig && fabs(z)<margin)
	return (z+margin)/(2.*margin);
else
	return z<=0. ? 0. : 1.;
}


void eval(int k)
{
int i,j,d,start,size;

if(k<databatch)
	{
	start=datanodes;
	size=datanodes;
	}
else
	{
	start=(datanodes+encodenodes+codenodes)%nodes;
	size=codenodes;
	}

j=start;
for(i=0;i<nodes;++i)
	{
	yp2[k][j]=0.;
	for(d=0;d<indeg[j];++d)
		yp2[k][j]+=xp2[k][innode[j][d]]*wp2[inedge[j][d]];

	yp2[k][j]/=wnorm;
	
	if(i<nodes-size)
		xp2[k][j]=act(yp2[k][j]-bp1[j]);
	
	j=(j+1)%nodes;
	}
}


void init(int item)
{
int k,e,i,d;

for(k=0;k<batchsize;++k)
	{
	if(k<databatch)
		{
		fsetpos(dataptr,&dataitem[perm[item+k]]);
		
		readdata(k);
			
		eval(k);
	
		for(i=0;i<codenodes;++i)
			codeval[i][k]=xp2[k][datanodes+encodenodes+i];
		}
	else
		{
		if(exhaust)
			exhaustcode(k,k-databatch);
		else
			randcode(k);
			
		eval(k);
		}
		
	for(e=0;e<edges;++e)
		w[k][e]=wp2[e];
		
	for(i=0;i<nodes;++i)
		b[k][i]=bp1[i];
	
	for(i=0;i<nodes;++i)
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

g=outdeg[n];

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
for(i=0;i<nodes;++i)
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
	
for(i=0;i<nodes;++i)
	{
	bp1[i]=0.;
	for(k=0;k<batchsize;++k)
		bp1[i]+=2.*bp2[k][i]-b[k][i];
		
	bp1[i]/=batchsize;
	}
}


void actproj(double *xa,double *ya,double *ba)
{
double s,y0,b0,d0,y1,b1,d1,xz,yz,bz,dz;

d0=(*xa)*(*xa);
s=(*ya-*ba+margin)/2.;

if(s<0.)
	{
	y0=*ya;
	b0=*ba;
	}
else
	{
	y0=*ya-s;
	b0=*ba+s;
	d0+=2.*s*s;
	}

d1=(1.-*xa)*(1.-*xa);
s=(*ya-*ba-margin)/2.;

if(s>0.)
	{
	y1=*ya;
	b1=*ba;
	}
else
	{
	y1=*ya-s;
	b1=*ba+s;
	d1+=2.*s*s;
	}

if(zig)
	{
	s=(*xa-(*ya-*ba)/(2.*margin)-.5)/(1.+.5/(margin*margin));
	
	xz=*xa-s;
	yz=*ya+s/(2.*margin);
	bz=*ba-s/(2.*margin);
	
	if(xz<0. || xz>1.)
		dz=1e10;
	else
		dz=(1.+.5/(margin*margin))*s*s;
	}
	
if(d0<d1)
	{
	if(zig && dz<d0)
		{
		*xa=xz;
		*ya=yz;
		*ba=bz;
		}
	else
		{
		*xa=0.;
		*ya=y0;
		*ba=b0;
		}
	}
else
	{
	if(zig && dz<d1)
		{
		*xa=xz;
		*ya=yz;
		*ba=bz;
		}
	else
		{
		*xa=1.;
		*ya=y1;
		*ba=b1;
		}
	}
}


void autoproj(double xauto,double *yout,double *bout)
{
double s;

if(xauto==0.)
	{
	s=(*yout-*bout+margin)/2.;
	if(s>0.)
		{
		*yout-=s;
		*bout+=s;
		}
	}
else if(xauto==1.)
	{
	s=(*yout-*bout-margin)/2.;
	if(s<0.)
		{
		*yout-=s;
		*bout+=s;
		}
	}
else
	{
	s=(*yout-*bout+margin*(1.-2.*xauto))/2.;
	*yout-=s;
	*bout+=s;
	}
}


void autonodeproj(int k,int i)
{
yp2[k][i]=y[k][i];
bp2[k][i]=b[k][i];
		
autoproj(xp2[k][i],&yp2[k][i],&bp2[k][i]);
}


void actnodeproj(int k,int i)
{
int d;

xp2[k][i]=0.;
for(d=0;d<outdeg[i];++d)
	xp2[k][i]+=x[k][outedge[i][d]];
			
xp2[k][i]/=outdeg[i];
yp2[k][i]=y[k][i];
bp2[k][i]=b[k][i];
		
actproj(&xp2[k][i],&yp2[k][i],&bp2[k][i]);
}


void proj2()
{
int k,i,d,in;
double norm;

for(k=0;k<databatch;++k)
for(i=0;i<nodes;++i)
	{
	if(i<datanodes)
		autonodeproj(k,i);
	else
		actnodeproj(k,i);
	}
	
for(k=databatch;k<databatch+codebatch;++k)
for(i=0;i<nodes;++i)
	{
	if(i>=datanodes+encodenodes && i<datanodes+encodenodes+codenodes)
		autonodeproj(k,i);
	else
		actnodeproj(k,i);
	}
	
for(i=0;i<nodes;++i)
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
	for(i=0;i<nodes;++i)
		{
		for(d=0;d<indeg[i];++d)
			{
			diff=xp1[k][inedge[i][d]]-xp2[k][innode[i][d]];
			x[k][inedge[i][d]]+=beta*diff;
			
			err+=diff*diff;
			}
			
		g=outdeg[i];
		
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

for(i=0;i<nodes;++i)
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


void train(int itermax,double *aveiter,double tol,double *errstart,double *errstop)
{
int batchcount,iter;
double err;

*aveiter=0.;
*errstart=0.;
*errstop=0.;

for(batchcount=0;batchcount<batchnum;++batchcount)
	{
	init(databatch*batchcount);
	
	err=RRR();
	iter=1;
	
	*errstart+=err;
	
	while(iter<itermax && err>tol)
		{
		err=RRR();
		
		fprintf(errptr,"%f\n",err);
		
		++iter;
		}
	
	*errstop+=err;
	
	*aveiter+=iter;
	
	netiter+=iter;
	}
	
*aveiter/=batchnum;
*errstart/=batchnum;
*errstop/=batchnum;
}


double autoerr(int k,int start,int size)
{
int i;
double err,diff;

err=0.;
for(i=start;i<start+size;++i)
	{
	diff=act(yp2[k][i]-bp1[i])-xp2[k][i];
	err+=diff*diff;
	}

return sqrt(err/size);
}


double dataerr(char *decodefile)
{
int k,i;
double err;
FILE *fp;

fp=fopen(decodefile,"a");

err=0.;

fsetpos(dataptr,&datastart);

for(k=0;k<datasize;++k)
	{
	readdata(0);
	eval(0);
	
	if(k<decodesamples)
		{
		for(i=0;i<datanodes;++i)
			fprintf(fp,"%lf ",act(yp2[0][i]-bp1[i]));
	
		fprintf(fp,"1\n");
		}
		
	err+=autoerr(0,0,datanodes);
	}

fclose(fp);

return err/datasize;
}


double codeerr(char *decodefile)
{
int c,i,samples;
double err;
FILE *fp;

if(codebatch==0)
	return 0.;
	
if(exhaust)
	samples=codebatch;
else
	samples=datasize;
	
fp=fopen(decodefile,"a");

err=0.;

for(c=0;c<samples;++c)
	{
	if(exhaust)
		exhaustcode(databatch,c);
	else
		randcode(databatch);
			
	eval(databatch);
	
	if(c<decodesamples)
		{
		for(i=0;i<datanodes;++i)
			fprintf(fp,"%lf ",act(yp2[databatch][i]-bp1[i]));
	
		fprintf(fp,"0\n");
		}
	
	err+=autoerr(databatch,datanodes+encodenodes,codenodes);
	}
	
fclose(fp);

return err/samples;
}


void printparam(char *paramfile)
{
FILE *fp;
int i,d;

fp=fopen(paramfile,"w");

fprintf(fp,"%lf %lf %d\n\n",wnorm,margin,zig);

edges=0;
for(i=0;i<nodes;++i)
	{
	fprintf(fp,"%lf\n",bp1[i]);
	
	for(d=0;d<indeg[i];++d)
		fprintf(fp,"%lf ",wp2[edges++]);
		
	fprintf(fp,"\n\n");
	}

fclose(fp);
}


double work()
{
return 1e-9*(double)edges*(double)batchsize*(double)netiter;
}


int compare(const void *a,const void *b)
{
return rank[*(int*)a] - rank[*(int*)b];
}


void shuffle()
{
int k;

for(k=0;k<datasize;++k)
	{
	rank[k]=rand();
	perm[k]=k;
	}
	
qsort(perm,datasize,sizeof(int),compare);
}


int main(int argc,char* argv[])
{
char *datafile,*netfile;
char *id,logfile[50],paramfile[50],decodefile[50];

char errfile[50];

int epochs,itermax;
int c,e;
double tol,aveiter,errstart,errstop,derr,cerr,bestdataerr,bestcodeerr,besterr2;
FILE *fp;


if(argc==14)
	{
	datafile=argv[1];
	netfile=argv[2];
	
	databatch=atoi(argv[3]);
	codebatch=atoi(argv[4]);
	epochs=atoi(argv[5]);
	
	itermax=atoi(argv[6]);
	tol=atof(argv[7]);
	
	beta=atof(argv[8]);
	wnorm=atof(argv[9]);
	margin=atof(argv[10]);
	zig=atoi(argv[11]);
	
	decodesamples=atoi(argv[12]);
	id=argv[13];
	}
else
	{
	fprintf(stderr,"expected 13 arguments: datafile, netfile, databatch, codebatch, epochs, itermax, tol, beta, wnorm, margin, zig, decodesamples, id\n");
	return 1;
	}
	
batchsize=databatch+codebatch;
	
if(!getnet(netfile))
	return 0;
	
if(zig==0 && codenodes<=10 && codebatch==1<<codenodes)
	exhaust=1;
else
	exhaust=0;

if(!opendata(datafile))
	return 0;
	
setup();

sprintf(logfile,"%s.log",id);
sprintf(paramfile,"%s.dat",id);
sprintf(decodefile,"%s.decode",id);

sprintf(errfile,"%s.err",id);

errptr=fopen(errfile,"w");

srand(time(NULL));
randparam();

initdata();

batchnum=datasize/databatch;

fp=fopen(logfile,"w");

for(c=0;c<argc-1;++c)
	fprintf(fp,"%s ",argv[c]);
fprintf(fp,"\n\n");

fprintf(fp,"epoch      GWMs   aveiter  errstart   errstop   dataerr   codeerr  dataerr*  codeerr*\n\n");
fclose(fp);

besterr2=1e10;

for(e=1;e<=epochs;++e)
	{
	shuffle();
	
	train(itermax,&aveiter,tol,&errstart,&errstop);
	
	fp=fopen(decodefile,"w");
	fclose(fp);
	
	derr=dataerr(decodefile);
	
	cerr=codeerr(decodefile);
	
	if(derr*derr+cerr*cerr<besterr2)
		{
		besterr2=derr*derr+cerr*cerr;
		bestdataerr=derr;
		bestcodeerr=cerr;
		}
	
	fp=fopen(logfile,"a");
	fprintf(fp,"%5d%10.2f%10.2e%10.2e%10.2e%10.6f%10.6f%10.6f%10.6f\n",e,work(),aveiter,errstart,errstop,derr,cerr,bestdataerr,bestcodeerr);
	fclose(fp);
	
	printparam(paramfile);
	
	if(besterr2==0.)
		break;
	}

closedata();

fclose(errptr);

return 0;
}
