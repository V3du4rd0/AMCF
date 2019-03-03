def init_ker(Np,wi,he):
	global TanNorCurv, distanceV,Evolution,BLin,MMxV,UPT,pixval,Au,Bu
	TanNorCurv="""
#include <math.h>

__global__ void TNK(double *M,double *T,  double *K, double *N, double mu)
{
 int Np=%(Np)s;
 
int row=blockIdx.y*blockDim.y+threadIdx.y;
int  col=blockIdx.x*blockDim.x+threadIdx.x;

double d1=0;
double t1=0;

double d2=0;
double t2=0;

double d_1=0;
double t_1=0;

double d_2=0;
double t_2=0;

int p1=0;
int p2=0;
int p_1=0;
int p_2=0;


if (row<Np && col<2 ){
	
	p1=fmodf(row+1,Np);
	p2=fmodf(row+2,Np);
	if (row==0){
		p_1=Np-1;
		p_2= Np-2;}
	else if(row==1){
		p_1=0;
		p_2=Np-1;}
	else{
		p_1=row-1;
		p_2=row-2;}
	
	
	d1=sqrt(  pow(M[ p1*2]-M[row*2] ,2) + pow(M[p1*2+1]-M[row*2+1] ,2));
	d2=sqrt(  pow(M[p2*2]-M[row*2] ,2) + pow(M[p2*2+1]-M[row*2+1] ,2));
	d_1=sqrt(  pow(M[p_1*2]-M[row*2] ,2) + pow(M[p_1*2+1]-M[row*2+1] ,2));
	d_2=sqrt(  pow(M[p_2*2]-M[row*2] ,2) + pow(M[p_2*2+1]-M[row*2+1] ,2));
	t1=(M[p1*2+col]-M[2*row+col])/d1;
	t2=(M[p2*2+col]-M[2*row+col])/d2;
	t_1=-(M[p_1*2+col]-M[2*row+col])/d_1;
	t_2=-(M[p_2*2+col]-M[2*row+col])/d_2;
	
	
T[row*2+col]=(-t2+4*t1+4.*t_1-t_2)/6.0;
K[row*2+col]=-(2.*mu/(d1+d_1))*(t1-t_1)-(1.-mu)*(2./(d2+d_2))*(t2-t_2);

	t1=(M[p1*2+1-col]-M[2*row+1-col])/d1;
	t2=(M[p2*2+1-col]-M[2*row+1-col])/d2;
	t_1=-(M[p_1*2+1-col]-M[2*row+1-col])/d_1;
	t_2=-(M[p_2*2+1-col]-M[2*row+1-col])/d_2;
	
	
N[row*2+col]=pow(-1.0,col+1)*(-t2+4*t1+4.*t_1-t_2)/6.0;

}}

"""%{'Np':Np}



	distanceV="""
#include <math.h>

__global__ void Dist(double *M, double d[%(Np)s])
{
 int Np=%(Np)s;
 
int row=blockIdx.y*blockDim.y+threadIdx.y;
int  col=blockIdx.x*blockDim.x+threadIdx.x;

if(row< Np && col<2){
	if(row<Np-1){
		d[row]=sqrt(pow(M[row*2]-M[(row+1)*2],2)+pow(M[row*2+1]-M[(row+1)*2+1],2));
	}
	else{
		d[row]=sqrt(pow(M[row*2]-M[0],2)+pow(M[row*2+1]-M[1],2));
	}
}
}
"""%{'Np':Np}

	BLin="""
#include <math.h>

__global__ void BLinear(double Dis[%(Np)s], double h, double C_lengpu, double b_gpu[%(Np)s])
{
int i =threadIdx.x;

b_gpu[i]=-h*(Dis[i]-C_lengpu*h)/(h*h*h);
if (i==0){
b_gpu[i]=0;
}
}
"""%{'Np':Np}
	MMxV="""
__global__ void MultiMxV(double M[ %(Np)s ][ %(Np)s ], double V[ %(Np)s ], double R[ %(Np)s ])
{
int i =threadIdx.x;
int k=0;
double Res=0;
for(k=0;k< %(Np)s ; k++){
	Res+=M[i][k]*V[k];
}
R[i]=Res;}
"""%{'Np':Np}

	UPT="""
 __global__ void Up(double *T, double A[%(Np)s])
 {
 int Np=%(Np)s;
 
 
int row=blockIdx.y*blockDim.y+threadIdx.y;
int  col=blockIdx.x*blockDim.x+threadIdx.x;

if (row<Np && col<2){
	T[row*2+col]*=A[row];
}
 
 }
 """%{'Np':Np}

	pixval="""
#include <math.h>
__global__ void pixv(double M[%(Np)s][2],double *Fpix, double Img[%(w)s][%(h)s])
{
 int Np=%(Np)s;
 int wi=%(w)s;
 int he=%(h)s;


int i=blockIdx.x*blockDim.x+threadIdx.x;

if(i<Np){
	Fpix[i]=0.0;
	
	double RSx=(100*M[i][0]+wi*0.5);
	double RSy=(-100*M[i][1]+he*0.5);
	int r=RSx; 
	int s=RSy; 
	
	
	if( r<wi && s<he && r>=0 && s>=0){
		Fpix[i]=Img[r][s];}
	else{
		Fpix[i]=255;}
	
}}	
"""%{'Np':Np,'w':wi,'h':he}

	Evolution="""
#include <math.h>

__global__ void Evolution(double *M, double *T, double *v,double dt)
{
 int Np=%(Np)s;
 
int row=blockIdx.y*blockDim.y+threadIdx.y;
int  col=blockIdx.x*blockDim.x+threadIdx.x;

if(row<Np && col<2){
	
		M[row*2+col]=M[row*2+col]+dt*(T[row*2+col]+v[row*2+col]);
	
	
}
}
"""%{'Np':Np}
	
	Bu="""
#include <math.h>
__host__ __device__ double dot(double p[2],double q[2])
{return(p[0]*q[0]+p[1]*q[1]);}
__host__ __device__ double norm(double p[2])
{return(sqrt(dot(p,p)));}
__host__ __device__ double Phi(double p[2],double q[2]) 
{
double PI=3.14159265358979;
double v[2];
v[0]=p[0]-q[0];
v[1]=p[1]-q[1];
return(  (-1./(2*PI))*log(norm(v)) );}
__host__ __device__ double func(double x[2],double x1[2],double x2[2],double k1,double k2,double xi)
{
double v[2];
v[0]=0.5*((1-xi)*x1[0]+(1+xi)*x2[0]);
v[1]=0.5*((1-xi)*x1[1]+(1+xi)*x2[1]);
double y[2];
y[0]=x2[0]-x1[0];
y[1]=x2[1]-x1[1];
return(  0.125*Phi(x,v)*((1-xi)*k1+(1+xi)*k2)*norm(y)  );}
extern "C"
__global__ void Bu(double M[%(Np)s][2],double Kur[%(Np)s][2],double N[%(Np)s][2],double q,double qp[2], double Pix[%(Np)s],double *F)
{
int i=blockIdx.x*blockDim.x+threadIdx.x;
int Np=%(Np)s;
int ii=0;
if(i<2*Np){
F[i]=0;

	if(i<Np){
		F[i]=q*Phi(M[i],qp);
	}
	else{
	ii=i-Np;
	for(int k=0; k<Np;k++){
		int kpos1=k+1;
		if(k==(Np-1)){
			kpos1=0;}

		double k1=dot(Kur[k],N[k])*Pix[k]/255.;
		double k2=dot(Kur[kpos1],N[kpos1])*Pix[k]/255.;
		F[i]+=(5./9.)*func(M[ii],M[k],M[kpos1],k1,k2,-sqrt(3./5.));
		F[i]+=(8./9.)*func(M[ii],M[k],M[kpos1],k1,k2,0.);
		F[i]+=(5./9.)*func(M[ii],M[k],M[kpos1],k1,k2,sqrt(3./5.));
	}}
	F[i]+=q*Phi(M[ii],qp)*Pix[ii]/255.;
}}
"""%{'Np':Np}

	Au="""
#include <math.h>
__host__ __device__ double dot(double p[2],double q[2])
{return(p[0]*q[0]+p[1]*q[1]);}
__host__ __device__ double norm(double p[2])
{return(sqrt(dot(p,p)));}
__host__ __device__ double Phi(double p[2],double q[2]) 
{
double PI=3.14159265358979;
double v[2];
v[0]=p[0]-q[0];
v[1]=p[1]-q[1];
return((-1./(2*PI) )*log(norm(v)) );}
__host__ __device__ double DPhi(double p[2],double q[2],double ny[2]) 
{
double PI=3.14159265358979;
double v[2];
v[0]=p[0]-q[0];
v[1]=p[1]-q[1];
return((-1./(2*PI) )*dot(v,ny)/pow(norm(v),2) );}
__host__ __device__ double func(double x[2],double x1[2],double x2[2], double n1[2],double n2[2], double xi,double c)
{
double y[2];
y[0]=0.5*((1-xi)*x1[0]+(1+xi)*x2[0]);
y[1]=0.5*((1-xi)*x1[1]+(1+xi)*x2[1]);
double ny[2];
ny[0]=0.5*((1-xi)*n1[0]+(1+xi)*n2[0]);
ny[1]=0.5*((1-xi)*n1[1]+(1+xi)*n2[1]);
double v[2];
v[0]=x2[0]-x1[0];
v[1]=x2[1]-x1[1];
return(0.25*DPhi(x,y,ny)*norm(v)*c);}
__host__ __device__ double func2(double x[2],double x1[2], double x2[2],double xi,double c)
{
double v[2];
v[0]=x2[0]-x1[0];
v[1]=x2[1]-x1[1];
double y[2];
y[0]=0.5*((1-xi)*x1[0]+(1+xi)*x2[0]);
y[1]=0.5*((1-xi)*x1[1]+(1+xi)*x2[1]);
return( 0.25*Phi(x,y)*norm(v)*c );}
extern "C"
__global__ void Au(double M[%(Np)s][2],double N[%(Np)s][2],double A[%(Np2)s][%(Np2)s],double q, double qp[2],double Pix[%(Np)s])
{
int j=blockIdx.x*blockDim.x+threadIdx.x;
int Np=%(Np)s;
int jj=0;
if(j<2*Np){
	if(j<Np){
	
	int jpos1=j+1;
	int jpos_1=j-1;
	if(j==0){
		jpos_1=Np-1;}
	if(j==(Np-1)){
		jpos1=0;}
	
	A[j][j+Np]=0.5;
	A[j][j+Np]+=(5./9.)*func(M[j],M[jpos_1],M[j],N[jpos_1],N[j],-sqrt(3./5.),1-sqrt(3./5.));
	A[j][j+Np]+=(8./9.)*func(M[j],M[jpos_1],M[j],N[jpos_1],N[j],0.,1);
	A[j][j+Np]+=(5./9.)*func(M[j],M[jpos_1],M[j],N[jpos_1],N[j],sqrt(3./5.),1+sqrt(3./5.));
	
	A[j][j+Np]+=(5./9.)*func(M[j],M[j],M[jpos1],N[j],N[jpos1],-sqrt(3./5.),1+sqrt(3./5.));
	A[j][j+Np]+=(8./9.)*func(M[j],M[j],M[jpos1],N[j],N[jpos1],0.,1);
	A[j][j+Np]+=(5./9.)*func(M[j],M[j],M[jpos1],N[j],N[jpos1],sqrt(3./5.),1-sqrt(3./5.));
	
	A[j][jpos1+Np]+=(5./9.)*func(M[j],M[j],M[jpos1],N[j],N[jpos1],-sqrt(3./5.),1-sqrt(3./5.));
	A[j][jpos1+Np]+=(8./9.)*func(M[j],M[j],M[jpos1],N[j],N[jpos1],0.,1);
	A[j][jpos1+Np]+=(5./9.)*func(M[j],M[j],M[jpos1],N[j],N[jpos1],sqrt(3./5.),1+sqrt(3./5.));
	
	A[j][jpos_1+Np]+=(5./9.)*func(M[j],M[jpos_1],M[j],N[jpos_1],N[j],-sqrt(3./5.),1+sqrt(3./5.));
	A[j][jpos_1+Np]+=(8./9.)*func(M[j],M[jpos_1],M[j],N[jpos_1],N[j],0.,1);
	A[j][jpos_1+Np]+=(5./9.)*func(M[j],M[jpos_1],M[j],N[jpos_1],N[j],sqrt(3./5.),1-sqrt(3./5.));
	
	for(int k=0;k<Np;k++){
		if(k!=j && k!=jpos_1){
		int kpos1=k+1;
		if(k==(Np-1)){
				kpos1=0;}
		
		A[j][k+Np]+=(5./9.)*func(M[j],M[k],M[kpos1],N[k],N[kpos1],-sqrt(3./5.),1+sqrt(3./5.));
		A[j][k+Np]+=(8./9.)*func(M[j],M[k],M[kpos1],N[k],N[kpos1],0.,1);
		A[j][k+Np]+=(5./9.)*func(M[j],M[k],M[kpos1],N[k],N[kpos1],-sqrt(3./5.),1-sqrt(3./5.));
		
		A[j][kpos1+Np]+=(5./9.)*func(M[j],M[k],M[kpos1],N[k],N[kpos1],-sqrt(3./5.),1-sqrt(3./5.));
		A[j][kpos1+Np]+=(8./9.)*func(M[j],M[k],M[kpos1],N[k],N[kpos1],0.,1);
		A[j][kpos1+Np]+=(5./9.)*func(M[j],M[k],M[kpos1],N[k],N[kpos1],-sqrt(3./5.),1+sqrt(3./5.));
	}}
	
	A[j][j]=(5./9.)*func2(M[j],M[j],M[jpos1],-sqrt(3./5.),1+sqrt(3./5.));	
	A[j][j]+=(8./9.)*func2(M[j],M[j],M[jpos1],0.,1);	
	A[j][j]+=(5./9.)*func2(M[j],M[j],M[jpos1],sqrt(3./5.),1-sqrt(3./5.));

	A[j][j]+=(5./9.)*func2(M[j],M[jpos_1],M[j],-sqrt(3./5.),1-sqrt(3./5.));	
	A[j][j]+=(8./9.)*func2(M[j],M[jpos_1],M[j],0.,1);	
	A[j][j]+=(5./9.)*func2(M[j],M[jpos_1],M[j],sqrt(3./5.),1+sqrt(3./5.));
	
	A[j][jpos_1]=(5./9.)*func2(M[j],M[jpos_1],M[j],-sqrt(3./5.),1+sqrt(3./5.));	
	A[j][jpos_1]+=(8./9.)*func2(M[j],M[jpos_1],M[j],0.,1);	
	A[j][jpos_1]+=(5./9.)*func2(M[j],M[jpos_1],M[j],sqrt(3./5.),1-sqrt(3./5.));
	
	A[j][jpos1]=(5./9.)*func2(M[j],M[j],M[jpos1],-sqrt(3./5.),1-sqrt(3./5.));	
	A[j][jpos1]+=(8./9.)*func2(M[j],M[j],M[jpos1],0.,1);	
	A[j][jpos1]+=(5./9.)*func2(M[j],M[j],M[jpos1],sqrt(3./5.),1+sqrt(3./5.));
	
	for(int k=0;k<Np; k++){
		if(k!=j && k!=(jpos_1)){
			int kpos1=k+1;
			
			if(k==(Np-1)){
				kpos1=0;}
				
			A[j][k]+=(5./9.)*func2(M[j],M[k],M[kpos1],-sqrt(3./5.),1+sqrt(3./5.));
			A[j][k]+=(8./9.)*func2(M[j],M[k],M[kpos1],0,1);
			A[j][k]+=(5./9.)*func2(M[j],M[k],M[kpos1],sqrt(3./5.),1-sqrt(3./5.));
			
			A[j][kpos1]+=(5./9.)*func2(M[j],M[k],M[kpos1],-sqrt(3./5.),1-sqrt(3./5.));
			A[j][kpos1]+=(8./9.)*func2(M[j],M[k],M[kpos1],0,1);
			A[j][kpos1]+=(5./9.)*func2(M[j],M[k],M[kpos1],sqrt(3./5.),1+sqrt(3./5.));
			
	}}
	
	
	}
	else{
	jj=j-Np;
	
	int jpos1=jj+1;
	int jpos_1=jj-1;
	if(jj==0){
		jpos_1=Np-1;}
	if(jj==(Np-1)){
		jpos1=0;}
	A[j][j]=(5./9.)*func(M[jj],M[jj],M[jpos1],N[jj],N[jpos1],-sqrt(3./5.),1+sqrt(3./5.));
	A[j][j]+=(8./9.)*func(M[jj],M[jj],M[jpos1],N[jj],N[jpos1],0.,1);
	A[j][j]+=(5./9.)*func(M[jj],M[jj],M[jpos1],N[jj],N[jpos1],sqrt(3./5.),1-sqrt(3./5.));
	
	A[j][j]+=(5./9.)*func(M[jj],M[jpos_1],M[jj],N[jpos_1],N[jj],-sqrt(3./5.),1-sqrt(3./5.));
	A[j][j]+=(8./9.)*func(M[jj],M[jpos_1],M[jj],N[jpos_1],N[jj],0.,1);
	A[j][j]+=(5./9.)*func(M[jj],M[jpos_1],M[jj],N[jpos_1],N[jj],sqrt(3./5.),1+sqrt(3./5.));
	
	A[j][j]+=0.5;
	
	A[j][jpos1+Np]=(5./9.)*func(M[jj],M[jj],M[jpos1],N[jj],N[jpos1],-sqrt(3./5.),1-sqrt(3./5.));
	A[j][jpos1+Np]+=(8./9.)*func(M[jj],M[jj],M[jpos1],N[jj],N[jpos1],0.,1);
	A[j][jpos1+Np]+=(5./9.)*func(M[jj],M[jj],M[jpos1],N[jj],N[jpos1],sqrt(3./5.),1+sqrt(3./5.));
	
	A[j][jpos_1+Np]=(5./9.)*func(M[jj],M[jpos_1],M[jj],N[jpos_1],N[jj],-sqrt(5./9.),1+sqrt(5./9.));
	A[j][jpos_1+Np]+=(8./9.)*func(M[jj],M[jpos_1],M[jj],N[jpos_1],N[jj],0.,1.);
	A[j][jpos_1+Np]+=(5./9.)*func(M[jj],M[jpos_1],M[jj],N[jpos_1],N[jj],sqrt(5./9.),1-sqrt(5./9.));
	
	for(int k=0;k<Np;k++){
		if (k!=jj && k!=(jpos_1)){
			int kpos1=k+1;
			
			if(k==(Np-1)){
				kpos1=0;}
				
			A[j][k+Np]+=(5./9.)*func(M[jj],M[k],M[kpos1],N[k],N[kpos1],-sqrt(3./5.),1+sqrt(3./5.));
			A[j][k+Np]+=(8./9.)*func(M[jj],M[k],M[kpos1],N[k],N[kpos1],0.,1);
			A[j][k+Np]+=(5./9.)*func(M[jj],M[k],M[kpos1],N[k],N[kpos1],sqrt(3./5.),1-sqrt(3./5.));
			
			A[j][kpos1+Np]+=(5./9.)*func(M[jj],M[k],M[kpos1],N[k],N[kpos1],-sqrt(3./5.),1-sqrt(3./5.));
			A[j][kpos1+Np]+=(8./9.)*func(M[jj],M[k],M[kpos1],N[k],N[kpos1],0.,1);
			A[j][kpos1+Np]+=(5./9.)*func(M[jj],M[k],M[kpos1],N[k],N[kpos1],sqrt(3./5.),1+sqrt(3./5.));
	}}	
	}
	
	

}}
"""%{'Np':Np,'Np2':2*Np}

	
