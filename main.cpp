#include <stdio.h>
#include <stdlib.h>
#include <math.h>



double Re=100;

double a_x = 0, b_x = 1;
double a_y = 0, b_y = b_x;
const int n = 10;
const int m = n;
double hx = (b_x-a_x)/(n);
double hy = (b_y-a_y)/(m);

const double tau=0.01;
double T=hx;
double CI=0;



void Solve1(double *u_temp,double *A,double *C,double *B,double *F,int N)
{
    double alfa[N+1];
    double beta[N+1];
    alfa[1]=-B[0]/C[0];
    beta[1]=F[0]/C[0];
    for(int i = 1; i<=N; i++)
    {
        alfa[i+1]= -B[i]/(A[i]*alfa[i]+C[i]);
        beta[i+1]=(F[i]-A[i]*beta[i])/(A[i]*alfa[i]+C[i]);
    }
    u_temp[N]=(F[N]-A[N]*beta[N])/(C[N]+A[N]*alfa[N]);
    for(int i = N-1;i>-1; i--)
    {
        u_temp[i] = alfa[i+1]*u_temp[i+1]+beta[i+1];
    }
}
double F_t()
{
    double p1 = 1;
    double p2 = 0;

    return (p2-p1)/(b_x-a_x);
}
void FindSpeed(double *u,int N, double h,double t)
{
    double *u_temp = new double[N+1];

    double *A=new double[N+1];
    double *B=new double[N+1];
    double *C=new double[N+1];
    double *F=new double[N+1];


    for(int i = 0; i<=N; i++)
        {
            if(i==0)
            {
                A[i]=0;
                B[i]=1;
                C[i]=0;
                F[i]=0;
            }
            else if(i==N)
            {
                A[i]=0;
                B[i]=1;
                C[i]=0;
                F[i]=0;
            }
            else
            {
                A[i] =  - 1.0 / (h*h) * (1.0/Re);
                B[i] =  1.0/tau + (2.0)/(h*h)* (1.0/Re);
                C[i] =  - 1.0 / (h*h)* (1.0/Re);
                F[i] =  CI + u[i]/tau;
            }
        }
        Solve1(u_temp,A,B,C,F,N);

        for(int i=0;i<=N;i++)
            u[i]=u_temp[i];


    delete []A;
    delete []B;
    delete []C;
    delete []F;
    delete []u_temp;
}
double Derivative_x(double **u,int i,int j,double h)
{
    return (u[i+1][j]-u[i-1][j])/(2*h);
}
double Derivative_y(double **u,int i,int j,double h)
{
    return (u[i][j+1]-u[i][j-1])/(2*h);
}
double Laplas_y(double **u,int i,int j,double h)
{
    return ((u[i][j+1]-2*u[i][j]+u[i][j-1])/(h*h));
}
double Laplas_x(double **u,int i,int j,double h)
{
    return ((u[i+1][j]-2*u[i][j]+u[i-1][j])/(h*h));
}
void null(double **u)
{
    for(int i=0;i<=n;i++)
        for(int j=0;j<=m;j++)
            u[i][j]=0;
}
double Norm(double **a,int lN,int lM)
{
    double max=fabs(a[0][0]);

    for(int i=0;i<=lN;i++)
    {
        for(int j=0;j<=lM;j++)
        {
            if(fabs(a[i][j])>max)
            {
                max=fabs(a[i][j]);
            }
        }
    }
    return max;
}
double U_touch(double t,double y)
{
    return - (4.0/(b_y*b_y)) * sin(t)*y*y + (4.0/b_y) * sin(t)*y + CI*sin(t);
}
double RightPart(double x,double y,double t)
{
    return - (4.0/(b_y*b_y)) * cos(t)*y*y + (4.0/b_y) * cos(t)*y + CI*cos(t) + (  (1.0/Re) * 2*(4.0/(b_y*b_y))*sin(t) );
}
double Aprox(double t_1,double t_n,double t_1_2,int i,int j,double *x,double *y)
{
    return ( U_touch(t_1,y[j]) - U_touch(t_n,y[j]) ) / tau
    - (1.0/Re) * ( ( U_touch(t_1,y[j+1]) - 2*U_touch(t_1,y[j]) + U_touch(t_1,y[j-1]) )/(hy*hy) + ( U_touch(t_n,y[j+1]) - 2*U_touch(t_n,y[j]) + U_touch(t_n,y[j-1]) )/(hy*hy)  )/2.0
    - RightPart(0,y[j],t_1_2);
}
double Aprox_ch(double **ut_1,double **ut_n,double **ut_1_2,double t_1_2,int i,int j,double *x,double *y)
{
    return ( ut_1[i][j] - ut_n[i][j] ) / tau
    - ( (1.0/Re) * ( ( ut_1_2[i+1][j] - 2*ut_1_2[i][j] + ut_1_2[i-1][j] ) / (hx*hx) ) )
    - ( (1.0/Re) * ( ( ut_1[i][j+1] - 2*ut_1[i][j] + ut_1[i][j-1]) /(hy*hy) + ( ut_n[i][j+1] - 2*ut_n[i][j] + ut_n[i][j-1] )/(hy*hy)  )/2.0 )
    - RightPart(0,y[j],t_1_2);
}
double Aprox_ch_1_2(double **ut_n,double **ut_1_2,double t_1_2,int i,int j,double *x,double *y)
{
    return ( ut_1_2[i][j] - ut_n[i][j] ) / (tau/2.0)
    - ( (1.0/Re) * ( ( ut_1_2[i+1][j] - 2*ut_1_2[i][j] + ut_1_2[i-1][j] ) / (hx*hx) ) )
    - ( (1.0/Re) * ( ( ut_n[i][j+1] - 2*ut_n[i][j] + ut_n[i][j-1] )/(hy*hy)  ) )
    - RightPart(0,y[j],t_1_2);
}

void SolveTransport(double *S,double **ne,double **u,double **u_n,int N,int M,double *x,double *y,double t)
{
    double *u_temp_x=new double[N+1];
    double *u_temp_y=new double[M+1];

    double **u_1_2=new double*[N+1];
    for(int i=0;i<=N;i++)
        u_1_2[i]=new double[M+1];

     for(int i=0;i<=N;i++)
        for(int j=0;j<=M;j++)
        {
            u_1_2[i][j]=u_n[i][j];
        }


        double *A=new double[M+1];
        double *B=new double[M+1];
        double *C=new double[M+1];
        double *F=new double[M+1];

        for(int i=1;i<N;i++)
        {
            for(int j = 0; j<=M; j++)
            {
                if(j == 0)
                {
                    A[j]=0;
                    B[j]=1;
                    C[j]=0;
                    F[j]=S[j];
                }
                else if(j == M)
                {
                    A[j]=0;
                    B[j]=1;
                    C[j]=0;
                    F[j]=S[j];
                }
                else
                {
                    A[j] =  - 1.0 / (hy*hy)*(1.0/Re);
                    B[j] =  2.0/tau+(2.0)/(hy*hy)*(1.0/Re);
                    C[j] =  - 1.0 / (hy*hy)*(1.0/Re);
                    F[j] = RightPart(x[i],y[j],t) + (1.0/Re)*Laplas_x(u_n,i,j,hx) + 2.0*u_n[i][j]/tau;
                }
            }
            Solve1(u_temp_y,A,B,C,F,M);
            double *check=new double[M+1];
            for(int j=1;j<M;j++)
            {
                check[j]= A[j]*u_temp_y[j-1]+B[j]*u_temp_y[j]+C[j]*u_temp_y[j+1]-F[j];
            }
            double max=0;
            for(int j=1;j<M;j++)
            {
                if(fabs(check[j])>max)
                    max=fabs(check[j]);
            }
            printf("CHECK = %lf\n",max);
            delete []check;

            for(int j=0;j<=M;j++)
                u_1_2[i][j]=u_temp_y[j];
        }

        for(int i=0;i<=N;i++)
        {
            for(int j=0;j<=M;j++)
            {
                u[i][j]=u_1_2[i][j];
                ne[i][j]=u_1_2[i][j];

            }
        }


        delete []A;
        delete []B;
        delete []C;
        delete []F;

        A=new double[N+1];
        B=new double[N+1];
        C=new double[N+1];
        F=new double[N+1];

        for(int j=1;j<M;j++)
        {
            for(int i = 0; i<=N; i++)
            {
                if(i==0)
                {
                    A[i]=0;
                    B[i]=1;
                    C[i]=0;
                    F[i]=S[j];
                }
                else if(i==N)
                {
                    A[i]=0;
                    B[i]=1;
                    C[i]=0;
                    F[i]=S[j];
                }
                else
                {
                    A[i] =  - 1.0 / (hx*hx)*(1.0/Re);
                    B[i] =  2.0/tau+(2.0)/(hx*hx)*(1.0/Re);
                    C[i] =  - 1.0 / (hx*hx)*(1.0/Re);
                    F[i] = RightPart(x[i],y[j],t)+ (1.0/Re)*Laplas_y(u_1_2,i,j,hy)+2.0*u_1_2[i][j]/tau;
                }
            }
        Solve1(u_temp_x,A,B,C,F,N);
        for(int i=0;i<=N;i++)
            u[i][j]=u_temp_x[i];
        }

        for(int i=1;i<N;i++)
        {
            u[i][0]=S[0];
            u[i][M]=S[M];
        }


    for(int i=0;i<=N;i++)
    {
        delete []u_1_2[i];
    }
    delete []u_1_2;
    delete []A;
    delete []B;
    delete []C;
    delete []F;
    delete []u_temp_x;
    delete []u_temp_y;
    u_1_2=NULL;
    A=NULL;
    B=NULL;
    C=NULL;
    F=NULL;
    u_temp_x=NULL;
    u_temp_y=NULL;
}



int main()
{
    printf("\n\n%lf %lf \n\n",hx,hy);

    double *x=new double[n+1];
    double *y=new double[m+1];


    double **u=new double*[n+1];
    for(int i=0;i<=n;i++)
        u[i]=new double[m+1];

    double **u_1_2=new double*[n+1];
    for(int i=0;i<=n;i++)
        u_1_2[i]=new double[m+1];

    double **_u=new double*[n+1];
    for(int i=0;i<=n;i++)
        _u[i]=new double[m+1];


    double **v=new double*[n+1];
    for(int i=0;i<=n;i++)
        v[i]=new double[m+1];



    double **u_n=new double*[n+1];
    for(int i=0;i<=n;i++)
        u_n[i]=new double[m+1];




    double **w=new double*[n+1];
    for(int i=0;i<=n;i++)
        w[i]=new double[m+1];


    x[0]=0;
    for(int i=1;i<=n;i++)
        x[i]=x[i-1]+hx;


    y[0]=0;
    for(int i=1;i<=m;i++)
        y[i]=y[i-1]+hy;




    for(int i=0;i<=n;i++)
        for(int j=0;j<=m;j++)
        {
            u[i][j]=0;
        }




    for(int i=0;i<=n;i++)
        for(int j=0;j<=m;j++)
        {
            u_1_2[i][j]=u[i][j];
            u_n[i][j]=u[i][j];
            _u[i][j]=u[i][j];
            v[i][j]=u[i][j];
        }

    /*
    int M=1000;
    double h = (b_y-a_y)/(M);
    double *u_one = new double[M+1];
    */
    double *u_one_g = new double[m+1];
    /*for(int i=0;i<=M;i++){
        u_one[i]=0;
    }*/

    for(int t=0;t<=T;t++)
    {
        for(int i=1;i<n;i++)
            for(int j=1;j<m;j++)
            {
                _u[i][j]=Aprox((t+1)*tau,(t)*tau,(t+1.0/2.0)*tau,i,j,x,y);
            }
        printf("Aprox = %lf\n",Norm(_u,n,m));

        /*
        FindSpeed(u_one,M,h,(t+1)*tau);


        u_one_g[0]=u_one[0];
        u_one_g[m]=u_one[M];
        u_one_g[1]=u_one[M/2];

        for(int i=0;i<=m;i++)
        {
            u_one_g[i] = u_one_g[i] + tau * F_t();
        }
        //*/

        //printf("\n");
        for(int i=0;i<=m;i++)
        {
            u_one_g[i] = U_touch((t+1)*tau,y[i]);
            //printf("%lf\n",u_one_g[i]);
        }
        //printf("\n");
        SolveTransport(u_one_g,u_1_2,u,u_n,n,m,x,y,(t+1.0/2.0)*tau);
        /*for(int i=0;i<=n;i++)
            for(int j=0;j<=m;j++)
            {
                u_1_2[i][j]=u[i][j]/2.0;
            }
        */


        for(int j=1;j<m;j++){
            for(int i=1;i<n;i++){
                w[i][j]=u[i][j]-u_one_g[j];
            }
        }

        printf("Norm(w) = %.14lf\n",Norm(w,n,m));
        /*
        printf("\n");
        for(int i=0;i<=m;i++)
        {
            printf("%lf\n",u[2][i]);
        }
        printf("\n");
        //*/
        for(int i=1;i<n;i++)
            for(int j=1;j<m;j++)
            {
                v[i][j]=Aprox_ch_1_2(u_n,u_1_2,(t+1.0/2.0)*tau,i,j,x,y);
            }
        printf("||aprox_ch_1_2|| = %.14lf\n",Norm(v,n,m));
        /*
        printf("\n U PPP\n");
        for(int j=0;j<=m;j++){
            for(int i=0;i<=n;i++){
                printf("%lf ",u[i][j]);
            }
            printf("\n");
        }

        printf("\n U Touch\n");
        for(int j=0;j<=m;j++){
            for(int i=0;i<=n;i++){
                printf("%lf ",u_one_g[j]);
            }
            printf("\n");

        }

        //*/


    }

    return 0;
}
