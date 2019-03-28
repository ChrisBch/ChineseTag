// kill.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#define pi 3.1415926
#define MIL pi/3000


void Xiangguanxishu(double alfa, double beta, double gama, double n, double h,double* ph,double* ps,double* pw);
void Dandaozhuyuan(double H, double Dq, double V0,double* palfa, double* pomg,double* pvq, double* ptf, double* psita,
					  double* pfai, double* pd_aF,double* pd_tF,double* pd_ZF);
void Dandaopiancha(double H, double Dq, double V0, double* pd_av0, double* pd_aK,
				   double* pd_tv0,double* pd_tK,double* pd_tD);

void Nihe(double H,double Dq,double V0,double* pGy, double* pGw,double* pGv, double* pGt);
void Zuobiaozhuyuan(double dj,double lq,double dq,double Hq,double* pdq,double* pDq,double* pqq, double* pepsq);
double Touying(double SB,double SS,double SF,double s1,double s2,double s3, 
			   double sita,double landa,double qq, double epsq, double vq, double vm);
void Fangxiangyuxian(double vq,double vm,double sita,double qq,double landa,double Hq,double Dq,
		 double* pvxd, double* ps1,double* ps2,double* ps3,double* peVx1,double* peVx2);
void B_xg(double sig_N, double sig_Z,double sig_V0p,double d_av0, double d_tv0, double Dq, double dq,double vm,
		  double eVx1, double eVx2, double* psig11, double* psig12, double* psig21,double* psig22);
void R_xg(double sigh_f, double sigh_b,double sigs_f, double sigs_b,double sigw_f, double sigw_b,
		  double Dq, double dq,double* psig11, double* psig12, double* psig21,double* psig22);
void R_xgh(double sig_f,double sig_b, double Dq, double dq, double* psig11,double* psig12, double* psig21, double* psig22);
void R_xgs(double sig_f,double sig_b, double Dq, double dq, double* psig11,double* psig12, double* psig21, double* psig22);
void R_xgw(double sig_f,double sig_b, double Dq, double dq, double* psig11,double* psig12, double* psig21, double* psig22);


void G_xg(double sig_V0, double sig_K,double sig_Fx,double sig_Fz,double Dq, double dq,double vm, 
		  double d_av0, double d_aK, double d_aF, double d_tv0, double d_tK, double d_tF,double d_ZF,
		  double eVx1, double eVx2, double* psig11, double* psig22, double* psig12,double* psig21);
void Xitong_wucha(double Dq, double dq, double af, double abq, double* pa1, double* pa2);

void F_xiefangcha(double b_xg[], double r_xgh[], double r_xgs[], double r_xgw[], double h, double s, double w,
				  double* psig11,double* psig12,double* psig21,double* psig22);
void C_xiefangcha(double g_xg[], double r_xgh[], double r_xgs[], double r_xgw[], double h, double s, double w,
				  double* psig11,double* psig12,double* psig21,double* psig22);
void Jifenbianhuan(double f[],double c[],double* psita1, double* psita2, double* plc1, double* plc2, double* plf1, double* plf2);
void Jifenxian(double l, double lf1, double lf2, double sita1, double sita2, double a1, double a2,
			   double yc1, double yc2, double* pc, double* pd, double* pe, double* pf);

double lapfunction(double x);
double jiewei(double l1,double l2,/*两个积分限参数*/double mu, double sigma/*正态分布的均值和均方差*/);

double Normal(double t);

double Beijihanshu(double m, double p, double n, double l/*命中面积的边长*/, double w, 
				   double lc1,double lc2,double lf1,double lf2, double sita1,double sita2,
				   double a1,double a2, double yc1, double yc2);
double Hanshu_S(double m, double p, double n, double l, double w,double lc1,double lc2,
				double lf1,double lf2,double sita1,double sita2,double a1, double a2, double yc2, int N1);
double Xinpusen_jifen(double m, double p, double n, double l,double w, double lc1,
					  double lc2,double lf1,double lf2,double sita1,double sita2,double a1,double a2,int N1,int N2);


double Killprobability(double H,  double Dq, double Pq,double landa/*倾斜角*/, double Vm, double SB, double SS, double SF,double m, double p, double n, double w, double V0, double h, double af, double abq, 
					 double sig_N, double sig_Z, double sigh_f,double sigh_b, double sigs_f, double sigs_b, double sigw_f,double sigw_b,double hk,double sd,double wd, 
		             double sig_V0p, double sig_V0, double sig_K, double sig_Fx, double sig_Fz);


int main(int argc, char* argv[])
{
	double H=209.1, Dq=2000, Pq=500, landa=0, Vm=250;
	//书中参数数值
	//double SB=48.72, SS=24.05, SF=3.96;
	double SB=51.6, SS=23.1, SF=5.9; 
	double m=1, p=2, n=18, w=1.5, V0=1175, h=9, af=-0.18*MIL, abq=-1.6*MIL;
	double sig_N=1.5*MIL,  sig_Z=1.6*MIL,  sigh_f=2.69*MIL, sigh_b=0.88*MIL; 
	double sigs_f=0*MIL,  sigs_b=0*MIL,  sigw_f=0*MIL, sigw_b=0*MIL;
	double hk=0.25, sd=5, wd=5;
	double sig_V0p=0.2,  sig_V0=0.2,  sig_K=1.3,  sig_Fx=2.2,  sig_Fz=2.2;

	double K=Killprobability(H,Dq,Pq,landa,Vm,SB,SS,SF,m,p,n,w,V0,h,af,abq,
		sig_N,sig_Z,sigh_f,sigh_b,sigs_f,sigs_b,sigw_f,sigw_b,hk,sd,wd,sig_V0p,
		sig_V0,sig_K,sig_Fx,sig_Fz);

	printf("K=%10.7f\n",K); 

	return 0;
}

//计算弱相关误差的分解系数
void Xiangguanxishu(double alfa, double beta, double gama, double n, double h,double* ph,double* ps,double* pw)
{
	//alfa-火控计算机输出误差相关衰减系数
    //beta-火炮随动系统输出误差相关衰减系数
    //gama-稳定系统输出误差相关衰减系数
    // n-每管一次点射的弹数
    // h-单管射速
	//ph,ps,pw分别表示指向火控、随动和稳定的弱相关分解系数

	double tao=(n-1)/h;

	*ph=2/(alfa*tao)-2*(1-exp(-alfa*tao))/pow((alfa*tao),2);
	*ps=2/(beta*tao)-2*(1-exp(-beta*tao))/pow((beta*tao),2);
	*pw=2/(gama*tao)-2*(1-exp(-gama*tao))/pow((gama*tao),2);

}

//提前点弹道诸元计算公式
void Dandaozhuyuan(double H,double Dq,double V0,double* palfa, double* pomg,double* pvq, double* ptf, double* psita,
					  double* pfai, double* pd_aF,double* pd_tF,double* pd_ZF)
{
	//H-目标提前点高度
//Dq-目标提前点斜距离
//V0-弹丸炮口初速

	double vq;//-弹丸存速
	double sita;//-提前点的弹道倾斜角
	double tf;//-弹丸飞行时间
	double omg;//-着角
	double alfa;//-高角
	double fai;//-射角

	double g=9.8;//重力加速度
	double C43=1.893;
	double CHD=C43*Dq*pow((1-2.1905*1e-5*H/2),(4.4));
    double dq=sqrt(Dq*Dq-H*H);//目标提前点水平距离

	//计算拟合系数
    //epsq-目标提前点高低角
    double epsq=asin(H/Dq);
	double Gy,Gw,Gv,Gt;
    Nihe(H,Dq,V0,&Gy,&Gw,&Gv,&Gt);

	alfa=cos(epsq)*g*Dq/(2*V0*V0)*Gy;
    omg=cos(epsq)*g*Dq/(2*V0*V0)*Gw;
    fai=epsq+alfa;
    sita=epsq-omg;
    vq=V0*cos(fai)/cos(sita)*Gv;
    tf=cos(epsq)/cos(fai)*Dq/V0*Gt;

	//临时变量
    double t1=0;double t2=0;double t3=0;
//计算d_aF
    double d_aF;
	t1 = (0.2595097e-2)-((0.1444118e-5)+(0.1348652e-7)*fai-(0.1775943e-9)*pow(fai,2)-(0.1877511e-11)*dq*fai)*dq*fai;
	t2 = ((0.5699964e-9)+(0.6742182e-13)*dq*fai-(0.4077875e-15)*pow(dq,2))*pow(dq,2);
	d_aF = pow((960*C43)/(1.508*V0),1.2)*(t1-t2);

//计算d_tF
    double d_tF;
	t1 = (0.1655924e-3)-((0.2358256e-7)+(0.1893686e-12)*dq*fai)*dq*fai+((0.1388580e-10)*fai+(0.1150672e-12)*dq-(0.1789927e-14)*dq*fai+(0.9054865e-17)*pow(dq,2))*pow(dq,2);
	t2 = ((0.3999206e-6)-(0.6522549e-8)*fai+(0.5080177e-11)*dq*fai)*pow(fai,2);
	d_tF = pow((960*C43)/(1.508*V0),2.55)*(t1+t2);

//计算d_ZF
    double d_ZF;
	t1 = (0.8966708e-1)-(0.1648553e-3)*dq+((0.8872561e-7)+(0.1133668e-10)*dq+(0.9842440e-11)*H)*pow(H,2);
	t2 = ((0.1667704e-6)-(0.1934115e-11)*H+(0.7819438e-11)*dq)*pow(dq,2);
	d_ZF = 	t1 + t2;
/*
//计算d_aF
	double d_aF;
t1=(0.1444118e-5+0.1348652e-7*fai-0.1775943e-9*fai*fai-0.1877511e-11*Dq*fai)*Dq*fai;
t2=(0.5699964e-9+0.6742182e-13*Dq*fai-0.4077875e-15*Dq*Dq)*Dq*Dq;
d_aF=pow((960*C43/(1.508*V0)),1.2)*(0.2595097e-2-t1-t2);

//计算d_tF
double d_tF;
t1=(0.2358256e-7+0.1893686e-12*Dq*fai)*Dq*fai;
t2=(0.138858e-10*fai+0.1150672e-12*Dq-0.1789927e-14*Dq*fai+0.9054865e-17*Dq*Dq)*Dq*Dq;
t3=(0.3999206e-6-0.6522549e-8*fai+0.5080177e-11*Dq*fai)*fai*fai;
d_tF=pow((960*C43/(1.508*V0)),2.55)*(0.1655924e-3-t1+t2+t3);

//计算d_ZF
double d_ZF;
t1=0.1648553e-3*dq;
t2=(0.8872561e-7+0.1133668e-10*dq+0.984244e-11*H)*H*H;
t3=(0.1667704e-6-0.1934115e-11*H+0.7819438e-11*dq)*dq*dq;
d_ZF=pow((960*C43/(1.508*V0)),1.32)*(0.8966708e-1-t1+t2+t3);
*/
//最后赋值
    *palfa=alfa;
	*pomg=omg;
	*pvq=vq;
	*ptf=tf;
	*psita=sita;
	*pfai=fai;
	*pd_aF=d_aF;
	*pd_tF=d_tF;
	*pd_ZF=d_ZF;

	return ;
}

//计算弹道诸元所需的拟合函数计算公式 
void Nihe(double H,double Dq,double V0,double* pGy, double* pGw,double* pGv, double* pGt)
{
	double *p=0;

	double dq=sqrt(Dq*Dq-H*H);
	double C43=1.893;
	double CHD=C43*Dq*pow((1-2.1905*1e-5*H/2),(4.4));
	//定义一些临时变量
double t1=0;
double t2=0;
double t3=0;

t1=(0.4406926e-6-0.283558e-13*CHD*V0)*CHD*V0;
t2=(0.3241693e-8-0.3199712e-10*V0+0.1166426e-11*CHD-0.2123943e-14*CHD*V0+0.588093e-16*CHD*CHD)*CHD*CHD;
t3=(0.144344e-5+0.6348184e-9*CHD-0.1934912e-8*V0-0.2794775e-12*CHD*V0+0.7250672e-12*V0*V0)*V0*V0;
double fy=exp(0.2347343+t1-t2-t3);

//计算Gy
double Gy;
if((H<=800)||((H>800)&&(H<=1500)&&(Dq<=2000)))
    Gy=fy;

if ((H>800)&&(Dq>2000))
    Gy=fy*pow((Dq/2000),0.075);

if ((H>1500)&&(Dq<=2000))
    Gy=fy*pow((H/1500),0.115);

if ((H>1500)&&(Dq>2000))
    Gy=fy*pow((Dq/2000),0.075)*pow((H/1500),0.115);

if ((H>4000)&&(Dq>6500))
    Gy=fy*pow((Dq/2000),0.075)*pow((H/1500),0.115)*pow((Dq/6500),0.45);

//计算中间系数Cr
double Cr;
if (dq<=2500)
    Cr=CHD;

if (dq>2500)
    Cr=CHD*pow((dq/2500),0.085);

if ((dq>3000)&&(dq<=5000)&&(H>=1500)&&(Dq>=4000))
    Cr=CHD*pow((dq/2500),0.085)*pow((Dq/4000),0.045)*pow((H/1500),0.14);

if ((dq>5000)&&(H>=1500))
    Cr=CHD*pow((dq/2500),0.085)*pow((Dq/5000),0.085)*pow((H/1500),0.18);

//计算Gw
double Gw;
t1=0.2152844e-3*Cr;
t2=(0.1660815e-6-0.4358138e-9*V0+0.2247045e-12*V0*V0-0.4583237e-13*Cr*V0)*Cr*V0;
t3=(0.7358717e-8-0.5559608e-10*V0+0.2756331e-11*Cr-0.4376855e-14*Cr*V0+0.1413296e-15*Cr*Cr)*Cr*Cr;
Gw=exp(-0.703454e-3+t1+t2-t3);

//计算Gv
double Gv;
t1=(0.1677301e-3+0.7276473e-7*V0-0.981196e-8*CHD+0.2390793e-10*CHD*V0-0.114599e-11*CHD*CHD)*CHD;
t2=(0.2952592e-6+0.1309007e-9*CHD-0.2306377e-9*V0)*V0*V0;
Gv=exp(-0.6546777e-1-t1+t2);

//计算中间系数Ct
double Ct;
if (H>3500)
    Ct=CHD*pow((Dq/3500),0.018);
else
    Ct=CHD;

//计算Gt
double Gt;
t1=(0.3163075e-6-0.1850771e-13*Ct*V0)*Ct*V0;
t2=(0.3603319e-8-0.2294901e-10*V0+0.9587002e-12*Ct-0.1503861e-14*Ct*V0+0.436919e-16*Ct*Ct)*Ct*Ct;
t3=(0.9955841e-6+0.4341479e-9*Ct-0.1309399e-8*V0-0.1810468e-12*Ct*V0+0.4795298e-12*V0*V0)*V0*V0;
Gt=exp(0.1658128+t1-t2-t3);
     //最后赋值

    *pGy=Gy;
	*pGw=Gw;
	*pGv=Gv;
	*pGt=Gt;
    return ;
}

//提前点坐标诸元计算公式
void Zuobiaozhuyuan(double dj,double lq,double dq,double Hq,double* pdq,double* pDq,double* pqq, double* pepsq)
{
	//double dq=sqrt(dj*dj+lq*lq);
	double Dq=sqrt(dq*dq+Hq*Hq);
	double qq=atan(dj/lq);
	double epsq=asin(Hq/Dq);

	//最后赋值
    *pdq=dq;
	*pDq=Dq;
	*pqq=qq;
	*pepsq=epsq;
    return ;
}

//拉普拉斯函数 按资料中定义 
double lapfunction(double x)
{
	double y=0;
	double a[6]={0.0705230784,0.0422820123,0.0092705272,0.0001520143,0.0002765672,0.0000430638};

	double sum=0;
	
	for (int i=0;i<=5;i++)
	{
		sum=sum+a[i]*pow(x,i+1);
	}
	y=1-pow(1+sum,-16);
	return y;
}

//计算给定积分限，带有两个参数的正态分布截尾积分，计算单发条件毁伤概率的积分所用 资料中定义
double jiewei(double l1,double l2,/*两个积分限参数*/double mu, double sigma/*正态分布的均值和均方差*/)
{
	double y=0;
	double f1=lapfunction((mu+l1)/(sqrt(2)*sigma));
	double f2=lapfunction((mu+l2)/(sqrt(2)*sigma));
	y=1.0/2*(f1-f2);
	return y;
}

//书中定义的计算正态分布公式
double Normal(double t)
{
	double a[4]={0.196854,0.115194,0.000344,0.019527};//拟合系数
	
	double y;
	double sum=0;
	
	for (int i=0;i<4;i++)
		{
			sum=sum+a[i]*pow(fabs(t),i+1);
		}
	if (t>=0)
	{
		y=1-0.5*pow(1+sum,-4);
	}
	if (t<0)
	{
		y=0.5*pow(1+sum,-4);
	}

	return y;
}



//计算目标投影面积
double Touying(double SB,double SS,double SF,double s1,double s2,double s3, 
			   double sita,double landa,double qq, double epsq, double vq, double vm)
{
	double ST=0;
	double t1,t2,t3;
	/*资料中定义t1=Sx*fabs(s1*sin(landa)*cos(qq)+s2*cos(landa)-s3*sin(landa)*sin(qq));
t2=Sy*fabs(s1*sin(qq)*s3*cos(qq));
t3=Sz*fabs(-s1*cos(landa)*cos(qq)+s2*sin(landa)+s3*cos(landa)*sin(qq));*/

	//书中定义
	double Wpm=fabs(vq*cos(epsq-sita)+vm*(cos(landa)*cos(qq)*cos(epsq)-sin(landa)*sin(epsq)));

	t1=SB*vq*fabs(cos(sita)*sin(landa)*cos(qq)+sin(sita)*cos(landa))/Wpm;
	t2=SS*vq*fabs(cos(sita)*sin(qq))/Wpm;
	t3=SF*fabs(vq*(sin(sita)*sin(landa)-cos(sita)*cos(landa)*cos(qq))-vm)/Wpm;

    ST=t1+t2+t3;

    return ST;
}

//计算初速、气密和斜距离改变时的弹道诸元
void Dandaopiancha(double H, double Dq, double V0, double* pd_av0, double* pd_aK,double* pd_tv0,double* pd_tK,double* pd_tD)
{
	double d_av0,  d_aK, d_tv0, d_tK, d_tD;

	double epsq=asin(H/Dq);
	double g=9.8;//重力加速度

	double Gy,Gw,Gv,Gt;
	Nihe(H,Dq,V0,&Gy,&Gw,&Gv,&Gt);

	double 	alfa_0=cos(epsq)*g*Dq/(2*V0*V0)*Gy;//初始值
	double fai=epsq+alfa_0;
    double tf_0=cos(epsq)/cos(fai)*Dq/V0*Gt;//初始值

	//计算d_av0 d_tv0
	double V0_1=V0*(1+0.01);//改变后的初速
	Nihe(H,Dq,V0_1,&Gy,&Gw,&Gv,&Gt);
	double alfa_1=cos(epsq)*g*Dq/(2*V0_1*V0_1)*Gy;//改变后的值
	fai=epsq+alfa_1;
    double tf_1=cos(epsq)/cos(fai)*Dq/V0_1*Gt;//改变后的值
    d_av0=fabs(alfa_0-alfa_1);
	d_tv0=fabs(tf_0-tf_1);

	//计算d_aK d_tK
	double A=2.1905e-5/2;double B=4.4;
	double d_H=fabs(0.01/(-A*B*pow(1-A*H,B-1)));//计算使得空气密度h改变1%时高度的变化量
	double H_1=H+d_H;
	Nihe(H_1,Dq,V0,&Gy,&Gw,&Gv,&Gt);
	   alfa_1=cos(epsq)*g*Dq/(2*V0*V0)*Gy;//改变后的值
	fai=epsq+alfa_1;
	   tf_1=cos(epsq)/cos(fai)*Dq/V0*Gt;//改变后的值
	d_aK=fabs(alfa_0-alfa_1);
	d_tK=fabs(tf_0-tf_1);

	//计算d_tD
	double Dq_1=Dq+100;//改变后的初速
	Nihe(H,Dq_1,V0,&Gy,&Gw,&Gv,&Gt);
	   alfa_1=cos(epsq)*g*Dq_1/(2*V0*V0)*Gy;//改变后的值
	fai=epsq+alfa_1;
       tf_1=cos(epsq)/cos(fai)*Dq_1/V0*Gt;//改变后的值
	d_tD=fabs(tf_0-tf_1);

	//最后赋值
	*pd_av0=d_av0;
	*pd_tv0=d_tv0;
	*pd_aK=d_aK;
	*pd_tK=d_tK;
	*pd_tD=d_tD;
}

//计算提前点的相对速度和各个方向余弦
void Fangxiangyuxian(double vq,double vm,double sita,double qq,double landa,double Hq,double Dq,
		 double* pvxd, double* ps1,double* ps2,double* ps3,double* peVx1,double* peVx2)
{
	double vxd,s1,s2,s3,eVx1,eVx2;

	double epsq=asin(Hq/Dq);//提前点的高低角
	vxd=sqrt(vq*vq+vm*vm+2*vq*vm*(cos(sita)*cos(qq)*cos(landa)-sin(sita)*sin(landa)));

	//按照书中定义
	s1=(vq*cos(sita)+vm*cos(landa)*cos(qq))/vxd;
    s2=-vm*cos(landa)*sin(qq)/vxd;
    s3=(vq*sin(sita)-vm*sin(landa))/vxd;

	eVx1=(s1*sin(landa)+s3*cos(landa)*cos(qq))/(s1*cos(epsq)+s3*sin(epsq));
	eVx2=cos(landa)*sin(qq)-(s2*(sin(landa)*sin(epsq)-cos(landa)*cos(qq)*cos(epsq)))/(s1*cos(epsq)+s3*sin(epsq));

	*pvxd=vxd;
	*ps1=s1;
	*ps2=s2;
	*ps3=s3;
	*peVx1=eVx1;
	*peVx2=eVx2;
}

//计算不相关误差协方差矩阵系数
void B_xg(double sig_N, double sig_Z,double sig_V0p,double d_av0, double d_tv0, double Dq, double dq,
		  double vm, double eVx1, double eVx2, double* psig11, double* psig12, double* psig21,double* psig22)
//sig_N, sig_Z-射弹高低角均方差和方位角均方差
//Dq,dq-目标提前点斜距离和水平斜距离
{
	double sigc11,sigc12,sigc21,sigc22;
	double sigcp11,sigcp12,sigcp21,sigcp22;
	
	sigc11=Dq*Dq*sig_N*sig_N;
	sigc22=dq*dq*sig_Z*sig_Z;
	sigc12=0;sigc21=0;

	sigcp11=sig_V0p*sig_V0p*(Dq*d_av0+vm*d_tv0*eVx1)*(Dq*d_av0+vm*d_tv0*eVx1);
	sigcp22=(sig_V0p*vm*d_tv0*eVx2)*(sig_V0p*vm*d_tv0*eVx2);
	sigcp12=sig_V0p*sig_V0p*(Dq*d_av0+vm*d_tv0*eVx1)*(vm*d_tv0*eVx2);
	sigcp21=sigcp12;

	*psig11=sigc11+sigcp11;
	*psig22=sigc22+sigcp22;
	*psig12=sigc12+sigcp12;
	*psig21=sigc21+sigcp21;
}

//计算弱相关误差协方差矩阵系数
void R_xg(double sigh_f, double sigh_b,double sigs_f, double sigs_b,double sigw_f, double sigw_b,
		  double Dq, double dq,double* psig11, double* psig12, double* psig21,double* psig22)
//Dq,dq-目标提前点斜距离和水平斜距离
//sigh_f、sigs_f、sigw_f分别表示火控（指挥仪）、随动和稳定系统输出射角均方差
//sigh_b、sigs_b、sigw_b分别表示火控（指挥仪）、随动和稳定系统输出方位角均方差
{
	*psig12=0;
	*psig21=0;
	*psig11=Dq*Dq*(sigh_f*sigh_f+sigs_f*sigs_f+sigw_f*sigw_f);
	*psig22=dq*dq*(sigh_b*sigh_b+sigs_b*sigs_b+sigw_b*sigw_b);
}

//计算弱相关误差协方差矩阵系数中火控部分
void R_xgh(double sig_f,double sig_b, double Dq, double dq, double* psig11,double* psig12, double* psig21, double* psig22)
{
	*psig12=0;
	*psig21=0;
	*psig11=Dq*Dq*(sig_f*sig_f);
	*psig22=dq*dq*(sig_b*sig_b);
}

//计算弱相关误差协方差矩阵系数中随动部分
void R_xgs(double sig_f,double sig_b, double Dq, double dq, double* psig11,double* psig12, double* psig21, double* psig22)
{
	*psig12=0;
	*psig21=0;
	*psig11=Dq*Dq*(sig_f*sig_f);
	*psig22=dq*dq*(sig_b*sig_b);
}

//计算弱相关误差协方差矩阵系数中稳定部分
void R_xgw(double sig_f,double sig_b, double Dq, double dq, double* psig11,double* psig12, double* psig21, double* psig22)
{
	*psig12=0;
	*psig21=0;
	*psig11=Dq*Dq*(sig_f*sig_f);
	*psig22=dq*dq*(sig_b*sig_b);
}

//计算强相关误差协方差矩阵系数
void G_xg(double sig_V0, double sig_K,double sig_Fx,double sig_Fz,double Dq, double dq,double vm, 
		  double d_av0, double d_aK, double d_aF, double d_tv0, double d_tK, double d_tF, double d_ZF,
		  double eVx1, double eVx2, double* psig11, double* psig12, double* psig21,double* psig22)
//Dq,dq-目标提前点斜距离和水平斜距离
//
//
{
	double sig11,sig22,sig21,sig12;

	double t1=0,t2=0,t3=0,t4=0;//临时变量

	t1=sig_V0*sig_V0*(Dq*d_av0+vm*d_tv0*eVx1)*(Dq*d_av0+vm*d_tv0*eVx1);
	t2=sig_K*sig_K*(Dq*d_aK+vm*d_tK*eVx1)*(Dq*d_aK+vm*d_tK*eVx1);
	t3=sig_Fx*sig_Fx*(Dq*d_aF+vm*d_tF*eVx1)*(Dq*d_aF+vm*d_tF*eVx1);
	sig11=t1+t2+t3;

	t1=sig_V0*sig_V0*(vm*d_tv0*eVx2)*(vm*d_tv0*eVx2);
	t2=sig_K*sig_K*(vm*d_tK*eVx2)*(vm*d_tK*eVx2);
	t3=sig_Fx*sig_Fx*(vm*d_tF*eVx2)*(vm*d_tF*eVx2);
	t4=sig_Fz*sig_Fz*d_ZF*d_ZF;
	sig22=t1+t2+t3+t4;

	t1=sig_V0*sig_V0*(Dq*d_av0+vm*d_tv0*eVx1)*vm*d_tv0*eVx2;
	t2=sig_K*sig_K*(Dq*d_aK+vm*d_tK*eVx1)*vm*d_tK*eVx2;
	t3=sig_Fx*sig_Fx*(Dq*d_aF+vm*d_tF*eVx1)*vm*d_tF*eVx2;
	sig12=t1+t2+t3;
	sig21=sig12;

	*psig12=sig12;
	*psig21=sig21;
	*psig11=sig11;
	*psig22=sig22;
}

//计算高炮系统误差
void Xitong_wucha(double Dq, double dq, double af, double abq, double* pa1, double* pa2)
{
	double a1,a2;
	a1=Dq*af;
	a2=dq*abq;

	*pa1=a1;
	*pa2=a2;
}

//计算非重复误差协方差系数
void F_xiefangcha(double b_xg[], double r_xgh[], double r_xgs[], double r_xgw[], double h, double s, double w,
				  double* psig11,double* psig12,double* psig21,double* psig22)
//b_xg[4]-不相关误差协方差系数，其中下标0-3分别表示协方差矩阵的11、12、21、22四个元素，其他参数也同样如此
//r_xgh[4]-火控输出协方差系数
//h-火控衰减系数
{
	double sig11, sig12, sig21, sig22;
	
	sig11=b_xg[0]+(1-h)*r_xgh[0]+(1-s)*r_xgs[0]+(1-w)*r_xgw[0];
	sig12=b_xg[1]+(1-h)*r_xgh[1]+(1-s)*r_xgs[1]+(1-w)*r_xgw[1];
	sig21=b_xg[2]+(1-h)*r_xgh[2]+(1-s)*r_xgs[2]+(1-w)*r_xgw[2];
	sig22=b_xg[3]+(1-h)*r_xgh[3]+(1-s)*r_xgs[3]+(1-w)*r_xgw[3];

	*psig11=sig11;
	*psig12=sig12;
	*psig21=sig21;
	*psig22=sig22;
}

//计算重复误差协方差系数
void C_xiefangcha(double g_xg[], double r_xgh[], double r_xgs[], double r_xgw[], double h, double s, double w,
				  double* psig11,double* psig12,double* psig21,double* psig22)
{
	double sig11, sig12, sig21, sig22;
	
	sig11=g_xg[0]+h*r_xgh[0]+s*r_xgs[0]+w*r_xgw[0];
	sig12=g_xg[1]+h*r_xgh[1]+s*r_xgs[1]+w*r_xgw[1];
	sig21=g_xg[2]+h*r_xgh[2]+s*r_xgs[2]+w*r_xgw[2];
	sig22=g_xg[3]+h*r_xgh[3]+s*r_xgs[3]+w*r_xgw[3];

	*psig11=sig11;
	*psig12=sig12;
	*psig21=sig21;
	*psig22=sig22;
}

//计算积分变换参数
void Jifenbianhuan(double f[],double c[],double* psita1, double* psita2, double* plc1, double* plc2, double* plf1, double* plf2)
//f[4]-表示非重复误差协方差矩阵元素，即书中的矩阵I，其中下标0-3分别表示协方差矩阵的11、12、21、22四个元素，
//c[4]-表示重复误差协方差矩阵元素，即书中的矩阵II，其中下标0-3分别表示协方差矩阵的11、12、21、22四个元素，
{
	double sita1,sita2,lc1,lc2,lf1,lf2;

	sita1=1.0/2*atan((2*f[1])/(f[0]-f[3]));
	lf1=f[0]*cos(sita1)*cos(sita1)+f[1]*sin(2*sita1)+f[3]*sin(sita1)*sin(sita1);
	lf2=f[0]*sin(sita1)*sin(sita1)-f[1]*sin(2*sita1)+f[3]*cos(sita1)*cos(sita1);

	sita2=1.0/2*atan((2*c[1])/(c[0]-c[3]));
	lc1=c[0]*cos(sita2)*cos(sita2)+c[1]*sin(2*sita2)+c[3]*sin(sita2)*sin(sita2);
	lc2=c[0]*sin(sita2)*sin(sita2)-c[1]*sin(2*sita2)+c[3]*cos(sita2)*cos(sita2);

	*psita1=sita1;
	*psita2=sita2;
	*plc1=lc1;
	*plc2=lc2;
	*plf1=lf1;
	*plf2=lf2;
}

//计算积分变换后的积分限和其他参数
void Jifenxian(double l, double lf1, double lf2, double sita1, double sita2, double a1, double a2,
			   double yc1, double yc2, double* pc, double* pd, double* pe, double* pf)
{
	double c,d,e,f,t,s;

	t=-yc1*cos(sita1-sita2)-yc2*sin(sita1-sita2)-a1*cos(sita1)-a2*sin(sita1);
	s=yc1*sin(sita1-sita2)-yc2*cos(sita1-sita2)+a1*sin(sita1)-a2*cos(sita1);

	c=(t-l)/sqrt(lf1);
	d=(t+l)/sqrt(lf1);
	e=(s-l)/sqrt(lf2);
	f=(s+l)/sqrt(lf2);

	*pc=c;
	*pd=d;
	*pe=e;
	*pf=f;
}

//毁歼概率中积分变换后可用yc1和yc2作为变量表示的二元函数，即被积函数
double Beijihanshu(double m, double p, double n, double l/*命中面积的边长的一半*/, double w, 
				   double lc1,double lc2,double lf1,double lf2, double sita1,double sita2,
				   double a1,double a2, double yc1, double yc2)
// m p n 分别为火炮门数、单炮身管数和点射长度
//w 为单发毁伤概率，即1.34
//c d e f 为标准正态分布积分限，由数jifenxian给出
//lc1 lc2 为二维正态分布参数，由函数jifenbianhuan给出
{
	double y=0;

	double c,d,e,f;
	Jifenxian(l,lf1,lf2,sita1,sita2,a1,a2,yc1,yc2,&c,&d,&e,&f);
	//用函数Normal计算两个标准正态分布函数的数值
	double Z1,Z2,F;
	Z1=Normal(d)-Normal(c);
	Z2=Normal(f)-Normal(e);
	
	F=Z1*Z2/w;//即书中的hv

	double E=pow(1-F,m*p*n);//这里暂时只考虑各炮基线已经修正的情况

	//y=E*exp(-1.0/2*(yc1*yc1/lc1+yc2*yc2/lc2))/(2*pi*sqrt(lc1*lc2));
	y=E*exp(-0.5*(yc1*yc1/lc1+yc2*yc2/lc2))/(2*pi*sqrt(lc1*lc2));
	return y;
}

//计算书中的函数S(yc2)
double Hanshu_S(double m, double p, double n, double l, double w,double lc1,double lc2,
				double lf1,double lf2,double sita1,double sita2,double a1, double a2, double yc2, int N1)
{
	double y=0;

	//积分限
	double b=4*sqrt(lc2);
//步长
	double h1=2*b/N1;
//临时变量
	double t1,t2,t3,t4;
	
	t1=Beijihanshu(m, p, n, l, w,lc1,lc2,lf1,lf2, sita1,sita2,a1,a2, -b, yc2);
	
	t2=Beijihanshu(m, p, n, l, w,lc1,lc2,lf1,lf2, sita1,sita2,a1,a2, b, yc2);

	t3=0;t4=0;
	for (int ix=1;ix<=N1-1;ix++)
	{
		t3=t3+Beijihanshu(m, p, n, l, w,lc1,lc2,lf1,lf2, sita1,sita2,a1,a2, -b+ix*h1, yc2);
		t4=t4+Beijihanshu(m, p, n, l, w,lc1,lc2,lf1,lf2, sita1,sita2,a1,a2, -b-h1/2+ix*h1, yc2);
	}
	t3=t3*2;
	t4=t4+Beijihanshu(m, p, n, l, w,lc1,lc2,lf1,lf2, sita1,sita2,a1,a2, -b-h1/2+N1*h1, yc2);
	t4=t4*4;

	y=h1/6*(t1+t2+t3+t4);

	return y;
}

//计算辛普森数值积分，利用书中给出的求和公式表示，即毁伤概率中1减去的部分
double Xinpusen_jifen(double m, double p, double n, double l,double w,
					  double lc1,double lc2,double lf1,double lf2,
					  double sita1,double sita2,double a1, double a2,  int N1, int N2)
//a,b-积分限 由于二重积分的四个积分限为对称，即a和-a，b和-b 
//N1 N2-两个方向的等分数
{
	double y=0;
	
	//积分限
	double a=4*sqrt(lc1);

	//步长
	double h2=2*a/N2;

//临时变量
	double t1,t2,t3,t4;	
	
	t1=Hanshu_S(m,p,n,l,w,lc1,lc2,lf1,lf2,sita1,sita2,a1,a2,-a,N1);

	t2=Hanshu_S(m,p,n,l,w,lc1,lc2,lf1,lf2,sita1,sita2,a1,a2,a,N1);

	t3=0;t4=0;
	for (int ix=1;ix<=N2-1;ix++)
	{
		t3=t3+Hanshu_S(m,p,n,l,w,lc1,lc2,lf1,lf2,sita1,sita2,a1,a2,-a+ix*h2,N1);
		t4=t4+Hanshu_S(m,p,n,l,w,lc1,lc2,lf1,lf2,sita1,sita2,a1,a2,-a-h2/2+ix*h2,N1);
	}
	t3=t3*2;
	t4=t4+Hanshu_S(m,p,n,l,w,lc1,lc2,lf1,lf2,sita1,sita2,a1,a2,-a-h2/2+N2*h2,N1);
	t4=t4*4;
	
	y=h2/6*(t1+t2+t3+t4);
	return y;
}

//计算毁伤概率
double Killprobability(double H,  double Dq, double Pq,double landa/*倾斜角*/, double Vm, double SB, double SS, double SF,double m, double p, double n, double w, double V0, double h, double af, double abq, 
					 double sig_N, double sig_Z, double sigh_f,double sigh_b, double sigs_f, double sigs_b, double sigw_f,double sigw_b,double hk,double sd,double wd, 
		             double sig_V0p, double sig_V0, double sig_K, double sig_Fx, double sig_Fz)
{
	double K;//毁歼概率

	double g=9.8;//重力加速度

    double dq=sqrt(Dq*Dq-H*H);//水平斜距离
	double qq=asin(Pq/dq);//航路角 弧度
	double epsq=asin(H/Dq);//高低角 弧度

	double vq=725.4;//-弹丸存速（斜距离3500，高低角300密位，高度1081.6米时射表的数据472.4）
	double sita=12.9*MIL;//-提前点的弹道倾斜角（斜距离3500，高低角300密位，高度1081.6米时射表的数据39.4）
	double tf=2.167;//-弹丸飞行时间（斜距离3500，高低角300密位，高度1081.6米时射表的数据4.722）
	double omg=0;//-着角
	double alfa=9.4*MIL;//-高角（斜距离3500，高低角300密位，高度1081.6米时射表的数据21.4）
	double fai=0;//-射角
    double d_aF=-0.02*MIL;//（斜距离3500，高低角300密位，高度1081.6米时射表的数据0.12）
	double d_tF=0.0009;//（斜距离3500，高低角300密位，高度1081.6米时射表的数据0.0046）
	double d_ZF=0.23*MIL*dq;//（斜距离3500，高低角300密位，高度1081.6米时射表的数据0.48）
    double d_av0=0.2*MIL; //（斜距离3500，高低角300密位，高度1081.6米时射表的数据0.6）
	double d_aK=0.0*MIL;//（斜距离3500，高低角300密位，高度1081.6米时射表的数据0.2）
	double d_tv0=0.025;//（斜距离3500，高低角300密位，高度1081.6米时射表的数据0.065）
	double d_tK=0.006;//（斜距离3500，高低角300密位，高度1081.6米时射表的数据0.026）
	double d_tD=0;

	//dandaozhuyuan(H,Dq,V0,&alfa, &omg,&vq,&tf,&sita,&fai,&d_aF,&d_tF,&d_ZF);
	//dandaopiancha(H,Dq,V0,&d_av0,&d_aK,&d_tv0,&d_tK,&d_tD);

	double vxd,  s1, s2, s3, eVx1, eVx2;
	Fangxiangyuxian(vq,Vm,sita,qq,landa,H,Dq,&vxd,&s1,&s2,&s3,&eVx1,&eVx2);

	double h_xg,s_xg,w_xg;//弱相关误差的相关系数
	Xiangguanxishu(hk,sd,wd,n,h,&h_xg,&s_xg,&w_xg);

	//计算目标投影面积
	double ST=Touying(SB,SS,SF,s1,s2,s3, sita,landa,qq, epsq, vq, Vm);
	double l=sqrt(ST)/2;
	
	//计算各种随机误差
	double b_xs[4],h_xs[4],s_xs[4],w_xs[4],g_xs[4];
	B_xg(sig_N,sig_Z,sig_V0p,d_av0, d_tv0, Dq, dq,Vm, eVx1, eVx2, &b_xs[0],&b_xs[1],&b_xs[2],&b_xs[3]);
	R_xgh(sigh_f,sigh_b,Dq,dq,&h_xs[0],&h_xs[1],&h_xs[2],&h_xs[3]);
	R_xgs(sigs_f,sigs_b,Dq,dq,&s_xs[0],&s_xs[1],&s_xs[2],&s_xs[3]);
	R_xgw(sigw_f,sigw_b,Dq,dq,&w_xs[0],&w_xs[1],&w_xs[2],&w_xs[3]);
	G_xg(sig_V0,sig_K,sig_Fx,sig_Fz,Dq,dq,Vm,d_av0,d_aK,d_aF,d_tv0,d_tK,d_tF,d_ZF,eVx1,eVx2,&g_xs[0],&g_xs[1],
		&g_xs[2],&g_xs[3]);
	
	double f_sig11,f_sig12,f_sig21,f_sig22,f_XieFC[4];//非重复协方差矩阵
	double c_sig11,c_sig12,c_sig21,c_sig22,c_XieFC[4];//重复协方差矩阵
	F_xiefangcha(b_xs,h_xs,s_xs,w_xs,h_xg,s_xg,w_xg,&f_sig11,&f_sig12,&f_sig21,&f_sig22);
	f_XieFC[0]=f_sig11,f_XieFC[1]=f_sig12,f_XieFC[2]=f_sig21,f_XieFC[3]=f_sig22;

	C_xiefangcha(g_xs,h_xs,s_xs,w_xs,h_xg,s_xg,w_xg,&c_sig11,&c_sig12,&c_sig21,&c_sig22); 
	c_XieFC[0]=c_sig11,c_XieFC[1]=c_sig12,c_XieFC[2]=c_sig21,c_XieFC[3]=c_sig22;

	//计算积分变换后的参数
	double sita1,  sita2,  lc1,  lc2,  lf1,  lf2;
	Jifenbianhuan(f_XieFC,c_XieFC,&sita1,&sita2,&lc1,&lc2,&lf1,&lf2);

	//计算系统误差
	double a1,a2;
	Xitong_wucha(Dq,dq,af,abq,&a1,&a2);

	int N1=50,N2=50;
	K=1-Xinpusen_jifen(m,p,n,l,w,lc1,lc2,lf1,lf2,sita1,sita2,a1,a2,N1,N2); 

	return K;
}

//优化速度将辛普森积分中的结果制成一张表
void Init_table(double m, double p, double n, double l, double w, 
				   double lc1,double lc2,double lf1,double lf2, double sita1,double sita2,
				   double a1,double a2, double yc1, double yc2, int N1,int N2)
{
}
