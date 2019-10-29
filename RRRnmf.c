#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NEWTONITER 50
#define NEWTONEPS 1e-12
#define ERRFILE 1

int datanodes,codenodes,edges;
int batchsize,batchcount;

double **x,***z,***zp1,**zp2,***w,***wp1,**wp2;
double **graminv,*gramdiag,**wp2inv;
double *zref,*wref;

double beta,wnorm;

int netiter=0;

FILE *dataptr;
fpos_t datastart;

FILE *errptr;


int opendata(char *datafile)
{
dataptr=fopen(datafile,"r");
if(!dataptr)
	{
	printf("datafile not found\n");
	return 0;
	}
		
fscanf(dataptr,"%d",&datanodes);

fgetpos(dataptr,&datastart);

return 1;
}


void closedata()
{
fclose(dataptr);
}


void setup()
{
int i,j,k;

zref=malloc(codenodes*sizeof(double));
wref=malloc(codenodes*sizeof(double));

x=malloc(batchsize*sizeof(double*));

z=malloc(batchsize*sizeof(double**));
zp1=malloc(batchsize*sizeof(double**));
zp2=malloc(batchsize*sizeof(double*));

w=malloc(batchsize*sizeof(double**));
wp1=malloc(batchsize*sizeof(double**));
wp2=malloc(datanodes*sizeof(double*));

wp2inv=malloc(codenodes*sizeof(double*));
graminv=malloc(codenodes*sizeof(double*));
gramdiag=malloc(codenodes*sizeof(double));

for(k=0;k<batchsize;++k)
	{
	x[k]=malloc(datanodes*sizeof(double));
	
	z[k]=malloc(datanodes*sizeof(double*));
	zp1[k]=malloc(datanodes*sizeof(double*));
	zp2[k]=malloc(codenodes*sizeof(double));
	
	w[k]=malloc(datanodes*sizeof(double*));
	wp1[k]=malloc(datanodes*sizeof(double*));
	
	for(i=0;i<datanodes;++i)
		{
		z[k][i]=malloc(codenodes*sizeof(double));
		zp1[k][i]=malloc(codenodes*sizeof(double));
		
		w[k][i]=malloc(codenodes*sizeof(double));
		wp1[k][i]=malloc(codenodes*sizeof(double));
		}
	}
	
for(i=0;i<datanodes;++i)	
	wp2[i]=malloc(codenodes*sizeof(double));
	
for(j=0;j<codenodes;++j)
	{
	wp2inv[j]=malloc(datanodes*sizeof(double));
	graminv[j]=malloc(codenodes*sizeof(double));
	}
}
	
	
int readdata(int k)
{
int i;

for(i=0;i<datanodes;++i)
	if(fscanf(dataptr,"%lf",&x[k][i])==EOF)
		{
		fsetpos(dataptr,&datastart);
		return 0;
		}

return 1;
}


int choleskyinverse(int m,double **a,double *d)
{
int i,j,k;
double s;

for(i=0;i<m;++i)
	{
	if(a[i][i]==0.)
		return 0;
		
	d[i]=a[i][i];
	}

for(i=0;i<m;++i)
	{
	a[i][i]=d[i];
	for(j=0;j<i;++j)
		a[i][i]-=a[j][i]*a[j][i];
		
	a[i][i]=sqrt(a[i][i]);
	
	for(k=i+1;k<m;++k)
		{
		a[i][k]=a[k][i];
		for(j=0;j<i;++j)
			a[i][k]-=a[j][i]*a[j][k];
			
		if(a[i][i]==0.)
			return 0;
		
		a[i][k]/=a[i][i];
		}
	}
	
for(i=0;i<m;++i)
	{
	if(a[i][i]==0.)
		return 0;
		
	a[i][i]=1./a[i][i];
	}
	
for(k=m-1;k>0;--k)
for(i=k-1;i>=0;--i)
	{
	a[i][k]*=-a[k][k];
	for(j=i+1;j<k;++j)
		a[i][k]-=a[i][j]*a[j][k];
		
	a[i][k]*=a[i][i];
	}
	
for(j=0;j<m;++j)
for(i=j;i<m;++i)
	{
	s=.0;
	for(k=i;k<m;++k)
		s+=a[j][k]*a[i][k];
		
	a[j][i]=s;
	}
	
for(i=0;i<m;++i)
for(j=i+1;j<m;++j)
	a[j][i]=a[i][j];
	
return 1;
}

	
void makeencoder()
{
int i,j1,j2;

for(j1=0;j1<codenodes;++j1)
for(j2=0;j2<codenodes;++j2)
	{
	graminv[j1][j2]=0.;
	for(i=0;i<datanodes;++i)
		graminv[j1][j2]+=wp2[i][j1]*wp2[i][j2];
	}
	
if(choleskyinverse(codenodes,graminv,gramdiag))
	{
	for(i=0;i<datanodes;++i)
	for(j1=0;j1<codenodes;++j1)
		{
		wp2inv[j1][i]=0.;
		for(j2=0;j2<codenodes;++j2)
			wp2inv[j1][i]+=wp2[i][j2]*graminv[j1][j2];
		}
	}
}

	
void encode(int k)
{
int i,j;

for(j=0;j<codenodes;++j)
	{
	zp2[k][j]=0.;
	
	for(i=0;i<datanodes;++i)
		zp2[k][j]+=wp2inv[j][i]*x[k][i];

	if(zp2[k][j]<0.)
		zp2[k][j]=0.;
	}
}
		

double urand()
{
return ((double)rand())/RAND_MAX;
}


int init()
{
int k,i,j;

for(k=0;k<batchsize;++k)
	if(!readdata(k))
		return 0;
	
makeencoder();

for(k=0;k<batchsize;++k)
	{
	encode(k);
	
	for(i=0;i<datanodes;++i)
	for(j=0;j<codenodes;++j)
		{
		w[k][i][j]=wp2[i][j];
		z[k][i][j]=zp2[k][j];
		}
	}
	
return 1;
}


double func(double c,double pdot,double qdot,double t)
{
double den;

den=1.-t*t;
	
return (pdot*(1.+t*t)+qdot*t)/(den*den)-c;
}


double dfunc(double pdot,double qdot,double t)
{
double den;

den=1.-t*t;
	
return (2.*pdot*t+qdot)/(den*den)+4.*t*(pdot*(1.+t*t)+qdot*t)/(den*den*den);
}


void bilinproj(double c,double *wold,double *zold,double *wnew,double *znew)
{
int j,iter;
double pdot,qdot,t,ta,tb,f,df,den;

pdot=0.;
qdot=0.;
for(j=0;j<codenodes;++j)
	{
	pdot+=wold[j]*zold[j];
	qdot+=wold[j]*wold[j]+zold[j]*zold[j];
	}

ta=0.;
f=func(c,pdot,qdot,ta);
tb=f>0. ? -1. : 1.;

for(iter=1;iter<=NEWTONITER;++iter)
	{
	df=dfunc(pdot,qdot,ta);
		
	t=ta-f/df;
	
	if(tb>ta)
		{
		if(t>tb)
			t=.5*(ta+tb);
		
		f=func(c,pdot,qdot,t);
		if(f>0.)
			tb=ta;
		}
	else
		{
		if(t<tb)
			t=.5*(ta+tb);
		
		f=func(c,pdot,qdot,t);
		if(f<0.)
			tb=ta;
		}
	
	if(fabs(t-ta)<NEWTONEPS)
		break;
		
	ta=t;
	}

den=1.-t*t;

for(j=0;j<codenodes;++j)
	{	
	wnew[j]=(wold[j]+t*zold[j])/den;
	znew[j]=(zold[j]+t*wold[j])/den;
	}
}


void proj1()
{
int i,j,k;

for(k=0;k<batchsize;++k)
for(i=0;i<datanodes;++i)
	{
	for(j=0;j<codenodes;++j)
		{
		wref[j]=2.*wp2[i][j]-w[k][i][j];
		zref[j]=2.*zp2[k][j]-z[k][i][j];
		}
		
	bilinproj(x[k][i],wref,zref,wp1[k][i],zp1[k][i]);
	}
}


void proj2()
{
int i,j,k,imax;
double norm,wmax;

for(j=0;j<codenodes;++j)
	{
	wmax=-1e10;
	for(i=0;i<datanodes;++i)
		{
		wp2[i][j]=0.;
		for(k=0;k<batchsize;++k)
			wp2[i][j]+=w[k][i][j];
		
		wp2[i][j]/=batchsize;
	
		if(wp2[i][j]>wmax)
			{
			wmax=wp2[i][j];
			imax=i;
			}
				
		if(wp2[i][j]<0.)
			wp2[i][j]=0.;
		}
		
	if(wmax<=0.)
		wp2[imax][j]=1.;
	}
	
for(j=0;j<codenodes;++j)
	{
	norm=0.;
	for(i=0;i<datanodes;++i)
		norm+=wp2[i][j]*wp2[i][j];
		
	norm=sqrt(norm)/wnorm;
	for(i=0;i<datanodes;++i)
		wp2[i][j]/=norm;
	}
	
for(k=0;k<batchsize;++k)
for(j=0;j<codenodes;++j)
	{
	zp2[k][j]=0.;
	for(i=0;i<datanodes;++i)
		zp2[k][j]+=z[k][i][j];
		
	zp2[k][j]/=datanodes;
	
	if(zp2[k][j]<0.)
		zp2[k][j]=0.;
	}
}


double RRR()
{
int i,j,k;
double err,diff;

err=0.;

proj2();
proj1();

for(k=0;k<batchsize;++k)
for(i=0;i<datanodes;++i)
for(j=0;j<codenodes;++j)
	{
	diff=wp1[k][i][j]-wp2[i][j];
	w[k][i][j]+=beta*diff;
			
	err+=diff*diff;
	
	diff=zp1[k][i][j]-zp2[k][j];
	z[k][i][j]+=beta*diff;
			
	err+=diff*diff;
	}
	
return sqrt(err/(batchsize*edges));
}


void randparam()
{
int i,j;
double norm;

for(j=0;j<codenodes;++j)
	{
	norm=0.;
	for(i=0;i<datanodes;++i)
		{
		wp2[i][j]=urand();
		norm+=wp2[i][j]*wp2[i][j];
		}
		
	norm=wnorm/sqrt(norm);
	
	for(i=0;i<datanodes;++i)
		wp2[i][j]*=norm;
	}
}
	

double reconerr()
{
int samples,k,i,j;
double err,diff;

err=0.;

if(batchcount>1)
	{
	makeencoder();
	
	samples=0;
	
	while(readdata(0))
		{
		++samples;
		encode(0);
		
		for(i=0;i<datanodes;++i)
			{
			diff=x[0][i];
			for(j=0;j<codenodes;++j)
				diff-=wp2[i][j]*zp2[0][j];
		
			err+=diff*diff;
			}
		}
	}
else
	{
	samples=batchsize;
	
	for(k=0;k<batchsize;++k)
	for(i=0;i<datanodes;++i)
		{
		diff=x[k][i];
		for(j=0;j<codenodes;++j)
			diff-=wp2[i][j]*zp2[k][j];
		
		err+=diff*diff;
		}
	}

return sqrt(err/(samples*datanodes));
}


void train(int itermax,double *aveiter,double errstop,double *errstart)
{
int iter;
double err;

batchcount=0;

*aveiter=0.;
*errstart=0.;

while(init())
	{
	++batchcount;
	
	err=RRR();
	iter=1;
	
	if(ERRFILE)
		fprintf(errptr,"%e\n",err);
	
	*errstart+=err;
	
	while(err>errstop && iter<itermax)
		{
		err=RRR();
		++iter;
		
		if(ERRFILE)
			fprintf(errptr,"%e\n",err);
		}
			
	netiter+=iter;
	*aveiter+=iter;
	}
	
*aveiter/=batchcount;
*errstart/=batchcount;
}


void printparam(char *paramfile)
{
FILE *fp;
int i,j,k;

fp=fopen(paramfile,"w");

for(i=0;i<datanodes;++i)
	{
	for(j=0;j<codenodes;++j)
		fprintf(fp,"%lf ",wp2[i][j]);	
	fprintf(fp,"\n");
	}
	
fprintf(fp,"\n");

for(k=0;k<batchsize;++k)
	{
	for(j=0;j<codenodes;++j)
		fprintf(fp,"%lf ",zp2[k][j]);	
	fprintf(fp,"\n");
	}
	
fclose(fp);
}


double work()
{
return 1e-9*(double)edges*(double)batchsize*(double)netiter;
}


int main(int argc,char* argv[])
{
char *datafile,*id,paramfile[50],logfile[50];
int c,e,epochs,itermax;
double errstart,errstop,aveiter;
FILE *fp;
clock_t start;
double elapsed;

if(argc==10)
	{
	datafile=argv[1];
	codenodes=atoi(argv[2]);
	batchsize=atoi(argv[3]);
	epochs=atoi(argv[4]);
	itermax=atoi(argv[5]);
	beta=atof(argv[6]);
	wnorm=atof(argv[7]);
	errstop=atof(argv[8]);
	id=argv[9];
	}
else
	{
	fprintf(stderr,"expected 9 arguments: datafile, codenodes, batchsize, epochs, itermax, beta, wnorm, errstop, id\n");
	return 1;
	}

if(!opendata(datafile))
	return 0;
	
edges=codenodes*datanodes;

setup();

sprintf(logfile,"%s.log",id);
sprintf(paramfile,"%s.dat",id);

srand(time(NULL));
randparam();

fp=fopen(logfile,"w");

for(c=0;c<=8;++c)
	fprintf(fp,"%s ",argv[c]);
fprintf(fp,"\n\n");

fprintf(fp,"epoch      GOPs  aveiter  errstart  reconerr\n\n");
fclose(fp);

if(ERRFILE)
	errptr=fopen("rrrerr","w");

start=clock();

for(e=1;e<=epochs;++e)
	{
	train(itermax,&aveiter,errstop,&errstart);
	
	fp=fopen(logfile,"a");
	fprintf(fp,"%5d%10.2f%9.2f%10.2e%10.6f\n",e,work(),aveiter,errstart,reconerr());
	fclose(fp);
	
	printparam(paramfile);
	}

elapsed=((double)(clock()-start))/CLOCKS_PER_SEC;

fp=fopen(logfile,"a");
fprintf(fp,"\nsecs: %12.3f\n",elapsed);
fprintf(fp,"%d iterations\n",netiter);
fclose(fp);

closedata();

if(ERRFILE)
	fclose(errptr);
	
return 0;
}