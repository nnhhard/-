#include <stdio.h>
#include <stdlib.h>
#include <math.h>



double Re=100;

double a_x = 0, b_x = 1;
double a_y = 0, b_y = b_x;
const int n = 1000;
const int m = n;
double hx = (b_x-a_x)/(n);
double hy = (b_y-a_y)/(m);

const double tau=hx;
double T=tau;
double CI=10;

int method = 1;


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
double Norm(double *a,int N)
{
    double max=fabs(a[0]);

    for(int i=0;i<=N;i++)
    {
        if(fabs(a[i])>max)
        {
            max=fabs(a[i]);
        }
    }
    return max;
}



double U_touch(double t,double y)
{
    return - (4.0/(b_y*b_y)) * sin(t)*y*y + (4.0/b_y) * sin(t)*y +
           CI*sin(t);
}
double RightPart(double x,double y,double t)
{
    if(method==0)
    {
        return - (4.0/(b_y*b_y)) * cos(t)*y*y + (4.0/b_y) * cos(t)*y +
            CI*cos(t) + (  (1.0/Re) * 2*(4.0/(b_y*b_y))*sin(t) );
    }
    else if(method==1)
    {
        return CI*sin(t);
    }
}
double Aprox(double t_1,double t_n,double t_1_2,int i,int j,double *x,double *y)
{
    return ( U_touch(t_1,y[j]) - U_touch(t_n,y[j]) ) / tau
    - (1.0/Re) * ( ( U_touch(t_1,y[j+1]) - 2*U_touch(t_1,y[j]) + U_touch(t_1,y[j-1]) )/(hy*hy) + ( U_touch(t_n,y[j+1]) - 2*U_touch(t_n,y[j]) + U_touch(t_n,y[j-1]) )/(hy*hy)  )/2.0
    - RightPart(0,y[j],t_1_2);
}
double Aprox_ne(double t_1,double t_n,double t_1_2,int i,int j,double *x,double *y)
{
    return ( U_touch(t_1,y[j]) - U_touch(t_n,y[j]) ) / tau
    - (1.0/Re) * ( ( U_touch(t_1,y[j+1]) - 2*U_touch(t_1,y[j]) + U_touch(t_1,y[j-1]) )/(hy*hy) )
    - RightPart(0,y[j],t_1);
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
void SolveTransport(double *S_1,double *S_1_2,double **u,double **u_n,int N,int M,double *x,double *y,double t,double **us_1_2)
{
    /*
    for(int j=0;j<=M;j++)
    {
        printf("%e ",S[j]);
    }
    printf("\n");
    */

    double *u_temp_x=new double[N+1];
    double *u_temp_y=new double[M+1];

    double **u_1_2=new double*[N+1];
    for(int i=0;i<=N;i++)
        u_1_2[i]=new double[M+1];

    double** check_aprox=new double*[N+1];
    for(int i=0;i<=N;i++)
    {
        check_aprox[i]=new double[M+1];
    }


    for(int i=0;i<=N;i++)
        for(int j=0;j<=M;j++)
        {
            u_1_2[i][j]=0;
        }

    for(int j=1;j<M;j++)
    {
        u_1_2[0][j]=S_1_2[j];
        u_1_2[N][j]=S_1_2[j];
    }


    double *check=new double[M+1];

    double *A=new double[M+1];
    double *B=new double[M+1];
    double *C=new double[M+1];
    double *F=new double[M+1];

    for(int j=0;j<=M;j++)
    {
        A[j]=0;
        B[j]=0;
        C[j]=0;
        F[j]=0;
        check[j]=0;
    }

    for(int i=1;i<N;i++)
    {
        for(int j = 0; j<=M; j++)
        {
            if(j == 0)
            {
                A[j]=0;
                B[j]=1;
                C[j]=0;
                F[j]=S_1_2[j];
            }
            else if(j == M)
            {
                A[j]=0;
                B[j]=1;
                C[j]=0;
                F[j]=S_1_2[j];
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

        for(int j=0;j<=M;j++)
        {
            if(j==0)
            {
                check[j]= B[j]*u_temp_y[j]+C[j]*u_temp_y[j+1]-F[j];
            }
            else if(j==M)
            {
                 check[j]= A[j]*u_temp_y[j-1]+B[j]*u_temp_y[j]-F[j];
            }
            else
                check[j]= A[j]*u_temp_y[j-1]+B[j]*u_temp_y[j]+C[j]*u_temp_y[j+1]-F[j];
        }
        //printf("CHECK_X = %.20lf\n",Norm(check,M));

        for(int j=0;j<=M;j++)
            u_1_2[i][j]=u_temp_y[j];

        }

        for(int i=1;i<N;i++)
        {
            for(int j=1;j<M;j++)
            {
                check_aprox[i][j]= ( (u_1_2[i][j]-u_n[i][j])/(tau/2.0) )
                - ( (1.0/Re) * ( Laplas_x(u_n,i,j,hx) + Laplas_y(u_1_2,i,j,hy) ) )
                - RightPart(x[i],y[j],t);
            }
        }

       // printf("Norm(check_aprox_1_2) = %.20lf\n",Norm(check_aprox,N,M));

        for(int i=0;i<=N;i++)
        {
            for(int j=0;j<=M;j++)
            {
                us_1_2[i][j]=u_1_2[i][j];
                u[i][j]=u_1_2[i][j];
            }
        }

        delete []A;
        delete []B;
        delete []C;
        delete []F;
        delete []check;


        A=new double[N+1];
        B=new double[N+1];
        C=new double[N+1];
        F=new double[N+1];
        check=new double[N+1];

        for(int j=1;j<M;j++)
        {
            for(int i = 0; i<=N; i++)
            {
                if(i==0)
                {
                    A[i]=0;
                    B[i]=1;
                    C[i]=0;
                    F[i]=S_1[j];
                }
                else if(i==N)
                {
                    A[i]=0;
                    B[i]=1;
                    C[i]=0;
                    F[i]=S_1[j];
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
        {
            if(i==0)
            {
                check[i]= B[i]*u_temp_x[i]+C[i]*u_temp_x[i+1]-F[i];
            }
            else if(i==N)
            {
                 check[i]= A[i]*u_temp_x[i-1]+B[i]*u_temp_x[i]-F[i];
            }
            else
                check[i]= A[i]*u_temp_x[i-1]+B[i]*u_temp_x[i]+C[i]*u_temp_x[i+1]-F[i];
        }
        //printf("CHECK_Y = %.20lf\n",Norm(check,M));
        for(int i=0;i<=N;i++)
            u[i][j]=u_temp_x[i];
        }

        for(int i=1;i<N;i++)
        {
            u[i][0]=S_1[0];
            u[i][M]=S_1[M];
        }

        for(int i=1;i<N;i++)
        {
            for(int j=1;j<M;j++)
            {
                check_aprox[i][j]= ( (u[i][j]-u_1_2[i][j])/(tau/2.0) )
                - ( (1.0/Re) * ( Laplas_x(u,i,j,hx) + Laplas_y(u_1_2,i,j,hy) ) )
                - RightPart(x[i],y[j],t);
            }
        }

      //  printf("Norm(check_aprox_1) = %.20lf\n",Norm(check_aprox,N,M));

        for(int i=1;i<N;i++)
        {
            for(int j=1;j<M;j++)
            {
                check_aprox[i][j]= ( (u[i][j]-u_n[i][j])/(tau) )
                - ( (1.0/Re) * ( Laplas_y(u_1_2,i,j,hy) ) )
                - ( (1.0/Re) * ( Laplas_x(u,i,j,hx) + Laplas_x(u_n,i,j,hx) ) / 2.0 )
                - RightPart(x[i],y[j],t);
            }
        }
      //  printf("Norm(check_aprox) = %.20lf\n",Norm(check_aprox,N,M));





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
void FindSpeed(double *u,double *u_n,int N, double h,double lTau,double t)
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
                F[i]=CI*sin(t);
            }
            else if(i==N)
            {
                A[i]=0;
                B[i]=1;
                C[i]=0;
                F[i]=CI*sin(t);
            }
            else
            {
                A[i] =  - 1.0 / (h*h) * (1.0/Re);
                B[i] =  1.0/lTau + (2.0)/(h*h)* (1.0/Re);
                C[i] =  - 1.0 / (h*h)* (1.0/Re);
                F[i] =  RightPart(0,0,t) + u_n[i]/lTau;
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

double *u_one;
double *u_one_;
double *u_one_n;
double *u_one_n_1_2;

void u_pr(double *u_one_g_1,double *u_one_g_1_2,double *y,int t,int lM)
{
    if(method==0)
    {
        for(int i=0;i<=lM;i++)
        {
            u_one_g_1[i] = U_touch((t+1)*tau,y[i]);
            u_one_g_1_2[i] = U_touch((t+1.0/2.0)*tau,y[i]);
        }
    }
    else if(method==1)
    {
        int K = 1;
        int M = lM * K;
        double h = (b_y-a_y)/(M);
        double Tau=tau/2.0;

        if(t==0)
        {
            u_one=new double[M+1];
            u_one_=new double[M+1];
            u_one_n_1_2=new double[M+1];
            u_one_n=new double[M+1];

            for(int i=0;i<=M;i++)
            {
                u_one[i]=0;
                u_one_n_1_2[i]=0;
                u_one_n[i]=0;

            }
        }

        double *check=new double[M+1];

        for(int i=0;i<=M;i++)
        {
            check[i]=0;
        }

        FindSpeed(u_one_n_1_2,u_one_n,M,h,Tau,(t+1.0/2.0)*tau);

        for(int i=1;i<M;i++)
        {
            check[i]= ( (u_one_n_1_2[i]-u_one_n[i])/Tau ) - 1.0/Re*( (u_one_n_1_2[i+1]-2*u_one_n_1_2[i]+u_one_n_1_2[i-1])/(h*h) ) - RightPart(0,0,(t+1.0/2.0)*tau);
        }
        printf("Norm(check 1d) = %e\n",Norm(check,M));

        for(int i=0;i<=lM;i++)
        {
            u_one_g_1_2[i]=u_one_n_1_2[K*i];
        }


        FindSpeed(u_one,u_one_n_1_2,M,h,Tau,(t+1)*tau);
        for(int i=1;i<M;i++)
        {
           check[i]= ( (u_one[i]-u_one_n_1_2[i])/Tau ) - 1.0/Re*( (u_one[i+1]-2*u_one[i]+u_one[i-1])/(h*h) ) - RightPart(0,0,(t+1.0)*tau);
        }
        printf("Norm(check 1d) = %e\n",Norm(check,M));


        FindSpeed(u_one_,u_one_n,M,h,tau,(t+1)*tau);
        for(int i=1;i<M;i++)
        {
           check[i]= ( (u_one_[i]-u_one_n[i])/tau ) - 1.0/Re*( (u_one_[i+1]-2*u_one_[i]+u_one_[i-1])/(h*h) ) - RightPart(0,0,(t+1.0)*tau);
        }
        printf("Norm(check 1d) = %e\n",Norm(check,M));


        ///Регулятор
        for(int i=0;i<=lM;i++)
        {
            u_one_g_1[i]=u_one[K*i];
            ///u_one_g_1[i]=u_one_[K*i];
        }

        for(int i=0;i<=M;i++)
        {
            u_one_n[i]=u_one[i];
            ///u_one_n[i]=u_one_[i];
        }
        ///Регулятор

        //printf("\n");
        for(int i=0;i<=lM;i++)
        {
            //printf("%lf %lf\n",u_one_g_1[i],u_one_g_1_2[i]);
        }
    }
}

int main()
{
    printf("\n\n%lf %lf \n\n",hx,hy);

    double *x=new double[n+1];
    double *y=new double[m+1];


    double **u=new double*[n+1];
    for(int i=0;i<=n;i++)
        u[i]=new double[m+1];

    double **r=new double*[n+1];
    for(int i=0;i<=n;i++)
        r[i]=new double[m+1];
    double **r_ne=new double*[n+1];
    for(int i=0;i<=n;i++)
        r_ne[i]=new double[m+1];


    double **_u=new double*[n+1];
    for(int i=0;i<=n;i++)
        _u[i]=new double[m+1];

    double **u_touch=new double*[n+1];
    for(int i=0;i<=n;i++)
        u_touch[i]=new double[m+1];

    double **us_1_2=new double*[n+1];
    for(int i=0;i<=n;i++)
        us_1_2[i]=new double[m+1];

    double **u_1_2=new double*[n+1];
    for(int i=0;i<=n;i++)
        u_1_2[i]=new double[m+1];

    double **check=new double*[n+1];
    for(int i=0;i<=n;i++)
        check[i]=new double[m+1];


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
            us_1_2[i][j]=0;
            r[i][j]=0;
            r_ne[i][j]=0;
            u_n[i][j]=0;
            _u[i][j]=0;
            u_1_2[i][j]=0;
        }


    double *u_one_g_1 = new double[m+1];
    double *u_one_g_1_2 = new double[m+1];

    for(int i=0;i<=m;i++){
        u_one_g_1[i]=0;
        u_one_g_1_2[i]=0;
    }



    for(int t=0;t<=T;t++)
    {
        /*
        for(int i=1;i<n;i++)
            for(int j=1;j<m;j++)
            {
                _u[i][j]=Aprox((t+1)*tau,(t)*tau,(t+1.0/2.0)*tau,i,j,x,y);
            }
        //printf("Aprox = %.20lf\n",Norm(_u,n,m));

         for(int i=1;i<n;i++)
            for(int j=1;j<m;j++)
            {
                _u[i][j]=Aprox_ne((t+1)*tau,(t)*tau,(t+1.0/2.0)*tau,i,j,x,y);
            }
        //printf("Aprox_ne = %.20lf\n",Norm(_u,n,m));
        */
        /*
        for(int i=0;i<=m;i++)
        {
            u_one_g_1[i] = U_touch((t+1)*tau,y[i]);
            u_one_g_1_2[i] = U_touch((t+1.0/2.0)*tau,y[i]);
        }
        */

        u_pr(u_one_g_1,u_one_g_1_2,y,t,m);

        SolveTransport(u_one_g_1,u_one_g_1_2,u,u_n,n,m,x,y,(t+1.0/2.0)*tau,us_1_2);

        for(int i=0;i<=n;i++)
        {
          for(int j=0;j<=m;j++)
          {
            r[i][j]=0;
            r_ne[i][j]=0;
          }
        }

        //точное решение
        for(int i=0;i<=n;i++)
            for(int j=0;j<=m;j++)
            {
                u_touch[i][j]=u_one_g_1[j];
            }



        //численное решение на шаге 1/2
        for(int i=0;i<=n;i++)
        {
            for(int j=0;j<=m;j++)
            {
                u_1_2[i][j]=u_one_g_1_2[j];
            }
        }
        for(int i=1;i<n;i++)
        {
            for(int j=1;j<m;j++)
            {
                r[i][j] = u_1_2[i][j] - us_1_2[i][j];
            }
        }
        printf("Norm(w_1_2) = %e\n",Norm(r,n,m));





        ///невязка схемы расщепления на точном решении
        for(int i=1;i<n;i++)
        {
            for(int j=1;j<m;j++)
            {
                r[i][j]= ( (u_touch[i][j]-u_n[i][j])/(tau) )
                - ( (1.0/Re) * ( Laplas_y(u_1_2,i,j,hy) ) )
                - ( (1.0/Re) * ( Laplas_x(u_touch,i,j,hx) + Laplas_x(u_n,i,j,hx) ) / 2.0 )
                - RightPart(x[i],y[j],(t+1.0/2.0)*tau);
            }
        }
        printf("Norm(r_split_ex) = %e\n",Norm(r,n,m));

        ///невязка схемы расщепления на численном решении
        for(int i=1;i<n;i++)
        {
            for(int j=1;j<m;j++)
            {
                r[i][j]= ( (u[i][j]-u_n[i][j])/(tau) )
                - ( (1.0/Re) * ( Laplas_y(us_1_2,i,j,hy) ) )
                - ( (1.0/Re) * ( Laplas_x(u,i,j,hx) + Laplas_x(u_n,i,j,hx) ) / 2.0 )
                - RightPart(x[i],y[j],(t+1.0/2.0)*tau);
            }
        }
        printf("Norm(r_split_num) = %e\n",Norm(r,n,m));


        ///невязка неявной схемы на точном решении
        for(int i=1;i<n;i++)
        {
            for(int j=1;j<m;j++)
            {
                r_ne[i][j]= ( (u_touch[i][j]-u_n[i][j])/(tau) )
                - ( (1.0/Re) * ( Laplas_y(u_touch,i,j,hy) ) )
                - ( (1.0/Re) * ( Laplas_x(u_touch,i,j,hx) ) )
                - RightPart(0,y[j],(t+1)*tau);
            }
        }
        printf("Norm(r_im_ex) = %e\n",Norm(r_ne,n,m));



        //невязка неявной схемы на численном решении
        for(int i=1;i<n;i++)
        {
            for(int j=1;j<m;j++)
            {
                r_ne[i][j]= ( (u[i][j]-u_n[i][j])/(tau) )
                - ( (1.0/Re) * ( Laplas_y(u,i,j,hy) ) )
                - ( (1.0/Re) * ( Laplas_x(u,i,j,hx) ) )
                - RightPart(0,y[j],(t+1)*tau);
            }
        }
        printf("Norm(r_im_num) = %e\n",Norm(r_ne,n,m));


        for(int j=0;j<=m;j++){
            for(int i=0;i<=n;i++){
             if ( (i==0) || (i==n) || (j==0) || (j==m) )
             {
                w[i][j]=0;
             }
            else {
                w[i][j]=u[i][j]-u_touch[i][j]; }
            }
        }
        printf("Norm(w) = %e\n",Norm(w,n,m));


        //Оператор от погрешности
        for(int i=1;i<n;i++)
        {
            for(int j=1;j<m;j++)
            {
                r_ne[i][j]=  (w[i][j]/(tau) )
                - ( (1.0/Re) * ( Laplas_y(w,i,j,hy) ) )
                - ( (1.0/Re) * ( Laplas_x(w,i,j,hx) ) );
            }
        }
        printf("Norm(Aw) = %e\n",Norm(r_ne,n,m));


        //Погрешность отдельных членов

        for(int j=1;j<m;j++)
        {
            for(int i=1;i<n;i++)
            {
                r_ne[i][j]= //( (u_touch[i][j]-u_n[i][j])/(tau) )
                //- ( (1.0/Re) * ( Laplas_y(u_touch,i,j,hy) ) )
                - ( (1.0/Re) * ( Laplas_x(u,i,j,hx) ) );
                //- ( (1.0/Re) * ( Laplas_x(us_1_2,i,j,hx) ) );
                //- RightPart(0,y[j],(t+1)*tau);
            }
        }
        printf("Norm(!) = %e\n",Norm(r_ne,n,m));

        return 0;





    }

    return 0;
}
