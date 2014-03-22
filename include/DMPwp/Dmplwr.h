#ifndef _Dmplwr
#define _Dmplwr

#include <string.h>
#include <cstdlib>
#include <cstdio>
#include <time.h>
#include <iostream>

#include "armadillo"
#define REALMIN 2.2251e-200
#define REALMAX 1.7977e200
#define PI 3.14

using namespace arma;
using namespace std;


struct Modeldmp
{
	int nbData,nbWeight;
	double dt,alpha_g,beta_g,alpha_z,x1_0,target_f,target_d,A,z, x1,x2,target,x2_d,x1_d,zd,x1_dd,move_time,tau;
	mat c,h,w,psi,t;
	mat demo,demo_d,demo_dd;
	mat rPos, rVel, rAcc;
};



class Dmplwr
{
public:
	Dmplwr();
	Dmplwr(int t,int s,double deltat);//模型参数初始化
  int Load(string path);//加载原始运动数据
	mat Derivative(mat data,double dt);//计算导数子函数
	int Learnlwr();//Dmp模型参数学习
	int Reprolwr(double dtar,double dvel,string path);//目标改变时，新数据生成，并输出在txt文件中

	Modeldmp model;

	 private:
};

#endif
