#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define M 3
#define N 3
#define pi 3.1415926
typedef struct complex{
    double real;
    double image;
}Complex;



void init(Complex *fx);
void DFT(Complex *fu,Complex fx[]);
void IDFT(Complex fu[],Complex *fx);


int main()
{
    Complex fx[M][N];
    Complex fu[M][N];

    init(&fx[0][0]);

    DFT(&fu[0][0],fx);

    IDFT(fu,&fx[0][0]);

    return 0;
}

void init(Complex *fx)
{
    printf("-----------------------------\n");
    printf("----------产生随机数---------\n");
    printf("-----------------------------\n");
    int i,j;
    srand((int)time(0));
    for(i=0;i<M;i++){
        for(j=0;j<N;j++){
            fx[i*M+j].real=(double)(rand()%10);
            fx[i*M+j].image=(double)(rand()%10);
            printf("fx[%2d][%2d]:%8.2lf+%8.2lfj\t",i,j,fx[i*M+j].real,fx[i*M+j].image);
        }
        printf("\n");
    }
}


void DFT(Complex *fu,Complex fx[])
{
    int u,v,x,y;
    double real,image;
    printf("-----------------------------\n");
    printf("---------傅里叶正变换--------\n");
    printf("-----------------------------\n");

    for(u=0;u<M;u++){
        for(v=0;v<N;v++){
            real=0;image=0;
            for(x=0;x<M;x++){
                for(y=0;y<N;y++){
                    real+=(fx[x*M+y].real*cos(-2*pi*((double)(u*x)/M+(double)(v*y)/N))-fx[x*M+y].image*sin(-2*pi*((double)(u*x)/M+(double)(v*y)/N)));
                    image+=(fx[x*M+y].real*sin(-2*pi*((double)(u*x)/M+(double)(v*y)/N))+fx[x*M+y].image*cos(-2*pi*((double)(u*x)/M+(double)(v*y)/N)));
                }
            }
            fu[u*M+v].real=real;
            fu[u*M+v].image=image;
            printf("fu[%2d][%2d]:%8.2lf+%8.2lfj\t",u,v,fu[u*M+v].real,fu[u*M+v].image);
        }
        printf("\n");
    }
}

void IDFT(Complex fu[],Complex *fx)
{
    int u,v,x,y;
    double real,image;
    printf("-----------------------------\n");
    printf("---------傅里叶逆变换--------\n");
    printf("-----------------------------\n");

    for(x=0;x<M;x++){
        for(y=0;y<N;y++){
            real=0;image=0;
            for(u=0;u<M;u++){
                for(v=0;v<N;v++){
                    real+=(fu[u*M+v].real*cos(2*pi*((double)(u*x)/M+(double)(v*y)/N))-fu[u*M+v].image*sin(2*pi*((double)(u*x)/M+(double)(v*y)/N)));
                    image+=(fu[u*M+v].real*sin(2*pi*((double)(u*x)/M+(double)(v*y)/N))+fu[u*M+v].image*cos(2*pi*((double)(u*x)/M+(double)(v*y)/N)));
                }
            }
            fx[x*M+y].real=real/(M*N);
            fx[x*M+y].image=image/(M*N);
            printf("fx[%2d][%2d]:%8.2lf+%8.2lfj\t",x,y,fx[x*M+y].real,fx[x*M+y].image);
        }
         printf("\n");
    }
}


