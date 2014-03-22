#include "Dmplwr.h"

//--------------------------------------------------------------
Dmplwr::Dmplwr(){
}

//参数初始化--------------------------------------------------
Dmplwr::Dmplwr(int n, int s,double deltat){
  model.nbData = n; //Number of training sample of the model
  model.nbWeight =n*s; 
  model.dt = deltat;
  model.alpha_g = 15;
  model.beta_g = model.alpha_g/4;
  model.alpha_z = 5;
  model.move_time = (n-1)*model.dt;
  model.tau = 1/model.move_time;
  model.t = linspace(0,model.move_time,model.nbData);
  //计算核函数中心值及接收域宽度
  model.h = zeros(model.nbWeight,1);
  model.c = zeros(model.nbWeight,1);
  vec temp = linspace(0,1,model.nbWeight);

  model.c = exp(-model.alpha_z*temp);    
 model.h = Derivative(model.c,1)*0.65;
 model.h = 1/(model.h%model.h); 
   //cout<<"temp data"<<model.h<<endl;
  cout<<"Parameter initialization of the model is finished"<<endl;
};
//---------------------------------------------------------------
//加载训练样本
int Dmplwr::Load(string path){
  mat DataM;
  string UpdatedPath;

  UpdatedPath = path+"demo.txt";
  DataM.load(UpdatedPath,raw_ascii);
  //DataM = DataM;
  model.demo = DataM.col(0);
 
 //----------------------------------------------------------------  
  model.demo_d =Derivative(model.demo,model.dt) ;
 
  model.demo_dd =Derivative(model.demo_d,model.dt) ; 
  model.demo = trans(model.demo);
  model.demo_d = trans(model.demo_d);
  model.demo_dd = trans(model.demo_dd);
  
  //cout<<"temp data"<<model.demo_dd<<endl;
  cout<<"Loading of the demonstartions is finished"<<endl;
  return 1;
}

//求导子函数--------------------------------------------------------------
mat Dmplwr::Derivative(mat DataH, double dt)
{
	mat rData = zeros(DataH.n_rows,DataH.n_cols);
  rData.rows(0,DataH.n_rows-2) = DataH.rows(1, DataH.n_rows-1);
  rData -= DataH;
  rData.row(rData.n_rows-1) = rData.row(rData.n_rows-2);
  rData = rData/dt;
  return rData;
}
//采用LWR学习参数-------------------------------------------------------------------------------
int Dmplwr::Learnlwr(){
	model.x1_0 = model.demo(0,0);
	model.target_f = model.demo(0,model.nbData-1);
	model.target_d = model.demo_d(0,model.nbData-1);
	model.A = model.demo(0,model.nbData-1)-model.demo(0,0);

	//cout<<":1:"<<model.x1_0<<endl<<":2:"<<model.target_f<<endl<<":3:"<<model.target_d<<endl<<";4;"<<model.A<<endl;

	mat target = (model.target_f-model.move_time*model.target_d)*ones(model.nbData,1);
	target += model.t*model.target_d;
	target = trans(target);
	

	mat Z = zeros(model.nbData,1);
	Z(0,0) = 1;
	for(int i=1;i<model.nbData;i++)
	{
		Z(i,0) =Z(i-1,0)-model.alpha_z*Z(i-1,0)*model.tau*model.dt;
	}
     Z = trans(Z);
	 
	mat PSI =  -model.h*ones(1,model.nbData)%((ones(model.nbWeight,1)*Z - model.c*ones(1,model.nbData))%(ones(model.nbWeight,1)*Z - model.c*ones(1,model.nbData)));
	PSI = exp(PSI);
	
		mat fdemo = (model.demo_dd/(model.tau*model.tau)-model.alpha_g*(1-Z)%(model.beta_g*(target-model.demo)+(model.target_d-model.demo_d)/model.tau))/model.A;
		
		vec wDnom = sum(PSI%(ones(model.nbWeight,1)*(Z%Z)),1);
		vec wNom = sum(PSI%(ones(model.nbWeight,1)*(Z%fdemo)),1);
		model.w = wNom/(wDnom+REALMIN);
		//cout<<"model.w"<<model.w<<endl;
		return 1;
}
//目标改变，生成新位置数据，可自由设置改变的幅度-----------------------------------------------------------------------
int Dmplwr::Reprolwr(double dtar,double dvel,string path){
	model.z = 1;
	model.x1 = 0;
	model.x2 = model.demo_d(0)*model.move_time;

	model.target_d +=dvel;
	model.target =model.target_f-model.move_time*model.target_d+dtar;

	//cout<<model.target_d<<endl<<model.target<<endl;

	model.rPos = zeros(1,model.nbData); 
	model.rVel = zeros(1,model.nbData);
	model.rAcc = zeros(1,model.nbData);
	model.rPos.col(0) = model.demo.col(0);
	model.rVel.col(0) = model.demo_d.col(0);
	model.rAcc.col(0) = model.demo_dd.col(0);

	//cout<<model.rPos.col(0)<<endl<<model.rVel.col(0)<<endl<<model.rAcc.col(0)<<endl;
	for(int i = 0;i<model.nbData;i++)
	{
		model.target += (model.target_d*model.dt);
		model.psi = exp(-model.h%(model.z-model.c)%(model.z-model.c));
		mat ftemp = sum(model.psi%model.w*model.z,0);
		ftemp = ftemp/sum(model.psi+REALMIN,0);
		double f = ftemp(0,0);
		model.x2_d = (model.alpha_g*(1-model.z)*(model.beta_g*(model.target-model.x1)+((model.target_d/model.tau)-model.x2))+model.A*f)*model.tau;
		model.x1_d = model.x2*model.tau;
		model.zd = -model.alpha_z*model.z*model.tau;
		model.x1_dd = model.x2_d*model.tau;

		model.x2 +=model.dt*model.x2_d;
		model.x1 +=model.dt*model.x1_d;
		model.z +=model.dt*model.zd;

		model.rPos(0,i) = model.x1;
		model.rVel(0,i) = model.x1_d;
		model.rAcc(0,i) = model.x1_dd;
	}
	cout<<"model.rPos:"<<model.rPos<<endl;
	model.rPos = trans(model.rPos);
	model.rVel = trans(model.rVel);
	model.rAcc = trans(model.rAcc);

// Saving the reproduction data
    string UpdatedPath;
    UpdatedPath = path + "data_rPos" + ".txt";
	model.rPos.save(UpdatedPath, raw_ascii);
                     
	UpdatedPath = path + "data_rVel" + ".txt";
	model.rVel.save(UpdatedPath, raw_ascii);

	UpdatedPath = path + "data_rAcc" + ".txt";
	model.rAcc.save(UpdatedPath, raw_ascii);

	cout<<"Reproduction to new target is finished"<<endl;

	return 1;
}