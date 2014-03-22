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
	Dmplwr(int t,int s,double deltat);//ģ�Ͳ�����ʼ��
  int Load(string path);//����ԭʼ�˶�����
	mat Derivative(mat data,double dt);//���㵼���Ӻ���
	int Learnlwr();//Dmpģ�Ͳ���ѧϰ
	int Reprolwr(double dtar,double dvel,string path);//Ŀ��ı�ʱ�����������ɣ��������txt�ļ���

	Modeldmp model;

	 private:
};

#endif
