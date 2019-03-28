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
double jiewei(double l1,double l2,/*���������޲���*/double mu, double sigma/*��̬�ֲ��ľ�ֵ�;�����*/);

double Normal(double t);

double Beijihanshu(double m, double p, double n, double l/*��������ı߳�*/, double w, 
				   double lc1,double lc2,double lf1,double lf2, double sita1,double sita2,
				   double a1,double a2, double yc1, double yc2);
double Hanshu_S(double m, double p, double n, double l, double w,double lc1,double lc2,
				double lf1,double lf2,double sita1,double sita2,double a1, double a2, double yc2, int N1);
double Xinpusen_jifen(double m, double p, double n, double l,double w, double lc1,
					  double lc2,double lf1,double lf2,double sita1,double sita2,double a1,double a2,int N1,int N2);


double Killprobability(double H,  double Dq, double Pq,double landa/*��б��*/, double Vm, double SB, double SS, double SF,double m, double p, double n, double w, double V0, double h, double af, double abq, 
					 double sig_N, double sig_Z, double sigh_f,double sigh_b, double sigs_f, double sigs_b, double sigw_f,double sigw_b,double hk,double sd,double wd, 
		             double sig_V0p, double sig_V0, double sig_K, double sig_Fx, double sig_Fz);


int main(int argc, char* argv[])
{
	double H=209.1, Dq=2000, Pq=500, landa=0, Vm=250;
	//���в�����ֵ
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

//������������ķֽ�ϵ��
void Xiangguanxishu(double alfa, double beta, double gama, double n, double h,double* ph,double* ps,double* pw)
{
	//alfa-��ؼ�������������˥��ϵ��
    //beta-�����涯ϵͳ���������˥��ϵ��
    //gama-�ȶ�ϵͳ���������˥��ϵ��
    // n-ÿ��һ�ε���ĵ���
    // h-��������
	//ph,ps,pw�ֱ��ʾָ���ء��涯���ȶ�������طֽ�ϵ��

	double tao=(n-1)/h;

	*ph=2/(alfa*tao)-2*(1-exp(-alfa*tao))/pow((alfa*tao),2);
	*ps=2/(beta*tao)-2*(1-exp(-beta*tao))/pow((beta*tao),2);
	*pw=2/(gama*tao)-2*(1-exp(-gama*tao))/pow((gama*tao),2);

}

//��ǰ�㵯����Ԫ���㹫ʽ
void Dandaozhuyuan(double H,double Dq,double V0,double* palfa, double* pomg,double* pvq, double* ptf, double* psita,
					  double* pfai, double* pd_aF,double* pd_tF,double* pd_ZF)
{
	//H-Ŀ����ǰ��߶�
//Dq-Ŀ����ǰ��б����
//V0-�����ڿڳ���

	double vq;//-�������
	double sita;//-��ǰ��ĵ�����б��
	double tf;//-�������ʱ��
	double omg;//-�Ž�
	double alfa;//-�߽�
	double fai;//-���

	double g=9.8;//�������ٶ�
	double C43=1.893;
	double CHD=C43*Dq*pow((1-2.1905*1e-5*H/2),(4.4));
    double dq=sqrt(Dq*Dq-H*H);//Ŀ����ǰ��ˮƽ����

	//�������ϵ��
    //epsq-Ŀ����ǰ��ߵͽ�
    double epsq=asin(H/Dq);
	double Gy,Gw,Gv,Gt;
    Nihe(H,Dq,V0,&Gy,&Gw,&Gv,&Gt);

	alfa=cos(epsq)*g*Dq/(2*V0*V0)*Gy;
    omg=cos(epsq)*g*Dq/(2*V0*V0)*Gw;
    fai=epsq+alfa;
    sita=epsq-omg;
    vq=V0*cos(fai)/cos(sita)*Gv;
    tf=cos(epsq)/cos(fai)*Dq/V0*Gt;

	//��ʱ����
    double t1=0;double t2=0;double t3=0;
//����d_aF
    double d_aF;
	t1 = (0.2595097e-2)-((0.1444118e-5)+(0.1348652e-7)*fai-(0.1775943e-9)*pow(fai,2)-(0.1877511e-11)*dq*fai)*dq*fai;
	t2 = ((0.5699964e-9)+(0.6742182e-13)*dq*fai-(0.4077875e-15)*pow(dq,2))*pow(dq,2);
	d_aF = pow((960*C43)/(1.508*V0),1.2)*(t1-t2);

//����d_tF
    double d_tF;
	t1 = (0.1655924e-3)-((0.2358256e-7)+(0.1893686e-12)*dq*fai)*dq*fai+((0.1388580e-10)*fai+(0.1150672e-12)*dq-(0.1789927e-14)*dq*fai+(0.9054865e-17)*pow(dq,2))*pow(dq,2);
	t2 = ((0.3999206e-6)-(0.6522549e-8)*fai+(0.5080177e-11)*dq*fai)*pow(fai,2);
	d_tF = pow((960*C43)/(1.508*V0),2.55)*(t1+t2);

//����d_ZF
    double d_ZF;
	t1 = (0.8966708e-1)-(0.1648553e-3)*dq+((0.8872561e-7)+(0.1133668e-10)*dq+(0.9842440e-11)*H)*pow(H,2);
	t2 = ((0.1667704e-6)-(0.1934115e-11)*H+(0.7819438e-11)*dq)*pow(dq,2);
	d_ZF = 	t1 + t2;
/*
//����d_aF
	double d_aF;
t1=(0.1444118e-5+0.1348652e-7*fai-0.1775943e-9*fai*fai-0.1877511e-11*Dq*fai)*Dq*fai;
t2=(0.5699964e-9+0.6742182e-13*Dq*fai-0.4077875e-15*Dq*Dq)*Dq*Dq;
d_aF=pow((960*C43/(1.508*V0)),1.2)*(0.2595097e-2-t1-t2);

//����d_tF
double d_tF;
t1=(0.2358256e-7+0.1893686e-12*Dq*fai)*Dq*fai;
t2=(0.138858e-10*fai+0.1150672e-12*Dq-0.1789927e-14*Dq*fai+0.9054865e-17*Dq*Dq)*Dq*Dq;
t3=(0.3999206e-6-0.6522549e-8*fai+0.5080177e-11*Dq*fai)*fai*fai;
d_tF=pow((960*C43/(1.508*V0)),2.55)*(0.1655924e-3-t1+t2+t3);

//����d_ZF
double d_ZF;
t1=0.1648553e-3*dq;
t2=(0.8872561e-7+0.1133668e-10*dq+0.984244e-11*H)*H*H;
t3=(0.1667704e-6-0.1934115e-11*H+0.7819438e-11*dq)*dq*dq;
d_ZF=pow((960*C43/(1.508*V0)),1.32)*(0.8966708e-1-t1+t2+t3);
*/
//���ֵ
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

//���㵯����Ԫ�������Ϻ������㹫ʽ 
void Nihe(double H,double Dq,double V0,double* pGy, double* pGw,double* pGv, double* pGt)
{
	double *p=0;

	double dq=sqrt(Dq*Dq-H*H);
	double C43=1.893;
	double CHD=C43*Dq*pow((1-2.1905*1e-5*H/2),(4.4));
	//����һЩ��ʱ����
double t1=0;
double t2=0;
double t3=0;

t1=(0.4406926e-6-0.283558e-13*CHD*V0)*CHD*V0;
t2=(0.3241693e-8-0.3199712e-10*V0+0.1166426e-11*CHD-0.2123943e-14*CHD*V0+0.588093e-16*CHD*CHD)*CHD*CHD;
t3=(0.144344e-5+0.6348184e-9*CHD-0.1934912e-8*V0-0.2794775e-12*CHD*V0+0.7250672e-12*V0*V0)*V0*V0;
double fy=exp(0.2347343+t1-t2-t3);

//����Gy
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

//�����м�ϵ��Cr
double Cr;
if (dq<=2500)
    Cr=CHD;

if (dq>2500)
    Cr=CHD*pow((dq/2500),0.085);

if ((dq>3000)&&(dq<=5000)&&(H>=1500)&&(Dq>=4000))
    Cr=CHD*pow((dq/2500),0.085)*pow((Dq/4000),0.045)*pow((H/1500),0.14);

if ((dq>5000)&&(H>=1500))
    Cr=CHD*pow((dq/2500),0.085)*pow((Dq/5000),0.085)*pow((H/1500),0.18);

//����Gw
double Gw;
t1=0.2152844e-3*Cr;
t2=(0.1660815e-6-0.4358138e-9*V0+0.2247045e-12*V0*V0-0.4583237e-13*Cr*V0)*Cr*V0;
t3=(0.7358717e-8-0.5559608e-10*V0+0.2756331e-11*Cr-0.4376855e-14*Cr*V0+0.1413296e-15*Cr*Cr)*Cr*Cr;
Gw=exp(-0.703454e-3+t1+t2-t3);

//����Gv
double Gv;
t1=(0.1677301e-3+0.7276473e-7*V0-0.981196e-8*CHD+0.2390793e-10*CHD*V0-0.114599e-11*CHD*CHD)*CHD;
t2=(0.2952592e-6+0.1309007e-9*CHD-0.2306377e-9*V0)*V0*V0;
Gv=exp(-0.6546777e-1-t1+t2);

//�����м�ϵ��Ct
double Ct;
if (H>3500)
    Ct=CHD*pow((Dq/3500),0.018);
else
    Ct=CHD;

//����Gt
double Gt;
t1=(0.3163075e-6-0.1850771e-13*Ct*V0)*Ct*V0;
t2=(0.3603319e-8-0.2294901e-10*V0+0.9587002e-12*Ct-0.1503861e-14*Ct*V0+0.436919e-16*Ct*Ct)*Ct*Ct;
t3=(0.9955841e-6+0.4341479e-9*Ct-0.1309399e-8*V0-0.1810468e-12*Ct*V0+0.4795298e-12*V0*V0)*V0*V0;
Gt=exp(0.1658128+t1-t2-t3);
     //���ֵ

    *pGy=Gy;
	*pGw=Gw;
	*pGv=Gv;
	*pGt=Gt;
    return ;
}

//��ǰ��������Ԫ���㹫ʽ
void Zuobiaozhuyuan(double dj,double lq,double dq,double Hq,double* pdq,double* pDq,double* pqq, double* pepsq)
{
	//double dq=sqrt(dj*dj+lq*lq);
	double Dq=sqrt(dq*dq+Hq*Hq);
	double qq=atan(dj/lq);
	double epsq=asin(Hq/Dq);

	//���ֵ
    *pdq=dq;
	*pDq=Dq;
	*pqq=qq;
	*pepsq=epsq;
    return ;
}

//������˹���� �������ж��� 
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

//������������ޣ�����������������̬�ֲ���β���֣����㵥���������˸��ʵĻ������� �����ж���
double jiewei(double l1,double l2,/*���������޲���*/double mu, double sigma/*��̬�ֲ��ľ�ֵ�;�����*/)
{
	double y=0;
	double f1=lapfunction((mu+l1)/(sqrt(2)*sigma));
	double f2=lapfunction((mu+l2)/(sqrt(2)*sigma));
	y=1.0/2*(f1-f2);
	return y;
}

//���ж���ļ�����̬�ֲ���ʽ
double Normal(double t)
{
	double a[4]={0.196854,0.115194,0.000344,0.019527};//���ϵ��
	
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



//����Ŀ��ͶӰ���
double Touying(double SB,double SS,double SF,double s1,double s2,double s3, 
			   double sita,double landa,double qq, double epsq, double vq, double vm)
{
	double ST=0;
	double t1,t2,t3;
	/*�����ж���t1=Sx*fabs(s1*sin(landa)*cos(qq)+s2*cos(landa)-s3*sin(landa)*sin(qq));
t2=Sy*fabs(s1*sin(qq)*s3*cos(qq));
t3=Sz*fabs(-s1*cos(landa)*cos(qq)+s2*sin(landa)+s3*cos(landa)*sin(qq));*/

	//���ж���
	double Wpm=fabs(vq*cos(epsq-sita)+vm*(cos(landa)*cos(qq)*cos(epsq)-sin(landa)*sin(epsq)));

	t1=SB*vq*fabs(cos(sita)*sin(landa)*cos(qq)+sin(sita)*cos(landa))/Wpm;
	t2=SS*vq*fabs(cos(sita)*sin(qq))/Wpm;
	t3=SF*fabs(vq*(sin(sita)*sin(landa)-cos(sita)*cos(landa)*cos(qq))-vm)/Wpm;

    ST=t1+t2+t3;

    return ST;
}

//������١����ܺ�б����ı�ʱ�ĵ�����Ԫ
void Dandaopiancha(double H, double Dq, double V0, double* pd_av0, double* pd_aK,double* pd_tv0,double* pd_tK,double* pd_tD)
{
	double d_av0,  d_aK, d_tv0, d_tK, d_tD;

	double epsq=asin(H/Dq);
	double g=9.8;//�������ٶ�

	double Gy,Gw,Gv,Gt;
	Nihe(H,Dq,V0,&Gy,&Gw,&Gv,&Gt);

	double 	alfa_0=cos(epsq)*g*Dq/(2*V0*V0)*Gy;//��ʼֵ
	double fai=epsq+alfa_0;
    double tf_0=cos(epsq)/cos(fai)*Dq/V0*Gt;//��ʼֵ

	//����d_av0 d_tv0
	double V0_1=V0*(1+0.01);//�ı��ĳ���
	Nihe(H,Dq,V0_1,&Gy,&Gw,&Gv,&Gt);
	double alfa_1=cos(epsq)*g*Dq/(2*V0_1*V0_1)*Gy;//�ı���ֵ
	fai=epsq+alfa_1;
    double tf_1=cos(epsq)/cos(fai)*Dq/V0_1*Gt;//�ı���ֵ
    d_av0=fabs(alfa_0-alfa_1);
	d_tv0=fabs(tf_0-tf_1);

	//����d_aK d_tK
	double A=2.1905e-5/2;double B=4.4;
	double d_H=fabs(0.01/(-A*B*pow(1-A*H,B-1)));//����ʹ�ÿ����ܶ�h�ı�1%ʱ�߶ȵı仯��
	double H_1=H+d_H;
	Nihe(H_1,Dq,V0,&Gy,&Gw,&Gv,&Gt);
	   alfa_1=cos(epsq)*g*Dq/(2*V0*V0)*Gy;//�ı���ֵ
	fai=epsq+alfa_1;
	   tf_1=cos(epsq)/cos(fai)*Dq/V0*Gt;//�ı���ֵ
	d_aK=fabs(alfa_0-alfa_1);
	d_tK=fabs(tf_0-tf_1);

	//����d_tD
	double Dq_1=Dq+100;//�ı��ĳ���
	Nihe(H,Dq_1,V0,&Gy,&Gw,&Gv,&Gt);
	   alfa_1=cos(epsq)*g*Dq_1/(2*V0*V0)*Gy;//�ı���ֵ
	fai=epsq+alfa_1;
       tf_1=cos(epsq)/cos(fai)*Dq_1/V0*Gt;//�ı���ֵ
	d_tD=fabs(tf_0-tf_1);

	//���ֵ
	*pd_av0=d_av0;
	*pd_tv0=d_tv0;
	*pd_aK=d_aK;
	*pd_tK=d_tK;
	*pd_tD=d_tD;
}

//������ǰ�������ٶȺ͸�����������
void Fangxiangyuxian(double vq,double vm,double sita,double qq,double landa,double Hq,double Dq,
		 double* pvxd, double* ps1,double* ps2,double* ps3,double* peVx1,double* peVx2)
{
	double vxd,s1,s2,s3,eVx1,eVx2;

	double epsq=asin(Hq/Dq);//��ǰ��ĸߵͽ�
	vxd=sqrt(vq*vq+vm*vm+2*vq*vm*(cos(sita)*cos(qq)*cos(landa)-sin(sita)*sin(landa)));

	//�������ж���
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

//���㲻������Э�������ϵ��
void B_xg(double sig_N, double sig_Z,double sig_V0p,double d_av0, double d_tv0, double Dq, double dq,
		  double vm, double eVx1, double eVx2, double* psig11, double* psig12, double* psig21,double* psig22)
//sig_N, sig_Z-�䵯�ߵͽǾ�����ͷ�λ�Ǿ�����
//Dq,dq-Ŀ����ǰ��б�����ˮƽб����
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

//������������Э�������ϵ��
void R_xg(double sigh_f, double sigh_b,double sigs_f, double sigs_b,double sigw_f, double sigw_b,
		  double Dq, double dq,double* psig11, double* psig12, double* psig21,double* psig22)
//Dq,dq-Ŀ����ǰ��б�����ˮƽб����
//sigh_f��sigs_f��sigw_f�ֱ��ʾ��أ�ָ���ǣ����涯���ȶ�ϵͳ�����Ǿ�����
//sigh_b��sigs_b��sigw_b�ֱ��ʾ��أ�ָ���ǣ����涯���ȶ�ϵͳ�����λ�Ǿ�����
{
	*psig12=0;
	*psig21=0;
	*psig11=Dq*Dq*(sigh_f*sigh_f+sigs_f*sigs_f+sigw_f*sigw_f);
	*psig22=dq*dq*(sigh_b*sigh_b+sigs_b*sigs_b+sigw_b*sigw_b);
}

//������������Э�������ϵ���л�ز���
void R_xgh(double sig_f,double sig_b, double Dq, double dq, double* psig11,double* psig12, double* psig21, double* psig22)
{
	*psig12=0;
	*psig21=0;
	*psig11=Dq*Dq*(sig_f*sig_f);
	*psig22=dq*dq*(sig_b*sig_b);
}

//������������Э�������ϵ�����涯����
void R_xgs(double sig_f,double sig_b, double Dq, double dq, double* psig11,double* psig12, double* psig21, double* psig22)
{
	*psig12=0;
	*psig21=0;
	*psig11=Dq*Dq*(sig_f*sig_f);
	*psig22=dq*dq*(sig_b*sig_b);
}

//������������Э�������ϵ�����ȶ�����
void R_xgw(double sig_f,double sig_b, double Dq, double dq, double* psig11,double* psig12, double* psig21, double* psig22)
{
	*psig12=0;
	*psig21=0;
	*psig11=Dq*Dq*(sig_f*sig_f);
	*psig22=dq*dq*(sig_b*sig_b);
}

//����ǿ������Э�������ϵ��
void G_xg(double sig_V0, double sig_K,double sig_Fx,double sig_Fz,double Dq, double dq,double vm, 
		  double d_av0, double d_aK, double d_aF, double d_tv0, double d_tK, double d_tF, double d_ZF,
		  double eVx1, double eVx2, double* psig11, double* psig12, double* psig21,double* psig22)
//Dq,dq-Ŀ����ǰ��б�����ˮƽб����
//
//
{
	double sig11,sig22,sig21,sig12;

	double t1=0,t2=0,t3=0,t4=0;//��ʱ����

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

//�������ϵͳ���
void Xitong_wucha(double Dq, double dq, double af, double abq, double* pa1, double* pa2)
{
	double a1,a2;
	a1=Dq*af;
	a2=dq*abq;

	*pa1=a1;
	*pa2=a2;
}

//������ظ����Э����ϵ��
void F_xiefangcha(double b_xg[], double r_xgh[], double r_xgs[], double r_xgw[], double h, double s, double w,
				  double* psig11,double* psig12,double* psig21,double* psig22)
//b_xg[4]-��������Э����ϵ���������±�0-3�ֱ��ʾЭ��������11��12��21��22�ĸ�Ԫ�أ���������Ҳͬ�����
//r_xgh[4]-������Э����ϵ��
//h-���˥��ϵ��
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

//�����ظ����Э����ϵ��
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

//������ֱ任����
void Jifenbianhuan(double f[],double c[],double* psita1, double* psita2, double* plc1, double* plc2, double* plf1, double* plf2)
//f[4]-��ʾ���ظ����Э�������Ԫ�أ������еľ���I�������±�0-3�ֱ��ʾЭ��������11��12��21��22�ĸ�Ԫ�أ�
//c[4]-��ʾ�ظ����Э�������Ԫ�أ������еľ���II�������±�0-3�ֱ��ʾЭ��������11��12��21��22�ĸ�Ԫ�أ�
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

//������ֱ任��Ļ����޺���������
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

//�ټ߸����л��ֱ任�����yc1��yc2��Ϊ������ʾ�Ķ�Ԫ����������������
double Beijihanshu(double m, double p, double n, double l/*��������ı߳���һ��*/, double w, 
				   double lc1,double lc2,double lf1,double lf2, double sita1,double sita2,
				   double a1,double a2, double yc1, double yc2)
// m p n �ֱ�Ϊ��������������������͵��䳤��
//w Ϊ�������˸��ʣ���1.34
//c d e f Ϊ��׼��̬�ֲ������ޣ�����jifenxian����
//lc1 lc2 Ϊ��ά��̬�ֲ��������ɺ���jifenbianhuan����
{
	double y=0;

	double c,d,e,f;
	Jifenxian(l,lf1,lf2,sita1,sita2,a1,a2,yc1,yc2,&c,&d,&e,&f);
	//�ú���Normal����������׼��̬�ֲ���������ֵ
	double Z1,Z2,F;
	Z1=Normal(d)-Normal(c);
	Z2=Normal(f)-Normal(e);
	
	F=Z1*Z2/w;//�����е�hv

	double E=pow(1-F,m*p*n);//������ʱֻ���Ǹ��ڻ����Ѿ����������

	//y=E*exp(-1.0/2*(yc1*yc1/lc1+yc2*yc2/lc2))/(2*pi*sqrt(lc1*lc2));
	y=E*exp(-0.5*(yc1*yc1/lc1+yc2*yc2/lc2))/(2*pi*sqrt(lc1*lc2));
	return y;
}

//�������еĺ���S(yc2)
double Hanshu_S(double m, double p, double n, double l, double w,double lc1,double lc2,
				double lf1,double lf2,double sita1,double sita2,double a1, double a2, double yc2, int N1)
{
	double y=0;

	//������
	double b=4*sqrt(lc2);
//����
	double h1=2*b/N1;
//��ʱ����
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

//��������ɭ��ֵ���֣��������и�������͹�ʽ��ʾ�������˸�����1��ȥ�Ĳ���
double Xinpusen_jifen(double m, double p, double n, double l,double w,
					  double lc1,double lc2,double lf1,double lf2,
					  double sita1,double sita2,double a1, double a2,  int N1, int N2)
//a,b-������ ���ڶ��ػ��ֵ��ĸ�������Ϊ�Գƣ���a��-a��b��-b 
//N1 N2-��������ĵȷ���
{
	double y=0;
	
	//������
	double a=4*sqrt(lc1);

	//����
	double h2=2*a/N2;

//��ʱ����
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

//������˸���
double Killprobability(double H,  double Dq, double Pq,double landa/*��б��*/, double Vm, double SB, double SS, double SF,double m, double p, double n, double w, double V0, double h, double af, double abq, 
					 double sig_N, double sig_Z, double sigh_f,double sigh_b, double sigs_f, double sigs_b, double sigw_f,double sigw_b,double hk,double sd,double wd, 
		             double sig_V0p, double sig_V0, double sig_K, double sig_Fx, double sig_Fz)
{
	double K;//�ټ߸���

	double g=9.8;//�������ٶ�

    double dq=sqrt(Dq*Dq-H*H);//ˮƽб����
	double qq=asin(Pq/dq);//��·�� ����
	double epsq=asin(H/Dq);//�ߵͽ� ����

	double vq=725.4;//-������٣�б����3500���ߵͽ�300��λ���߶�1081.6��ʱ��������472.4��
	double sita=12.9*MIL;//-��ǰ��ĵ�����б�ǣ�б����3500���ߵͽ�300��λ���߶�1081.6��ʱ��������39.4��
	double tf=2.167;//-�������ʱ�䣨б����3500���ߵͽ�300��λ���߶�1081.6��ʱ��������4.722��
	double omg=0;//-�Ž�
	double alfa=9.4*MIL;//-�߽ǣ�б����3500���ߵͽ�300��λ���߶�1081.6��ʱ��������21.4��
	double fai=0;//-���
    double d_aF=-0.02*MIL;//��б����3500���ߵͽ�300��λ���߶�1081.6��ʱ��������0.12��
	double d_tF=0.0009;//��б����3500���ߵͽ�300��λ���߶�1081.6��ʱ��������0.0046��
	double d_ZF=0.23*MIL*dq;//��б����3500���ߵͽ�300��λ���߶�1081.6��ʱ��������0.48��
    double d_av0=0.2*MIL; //��б����3500���ߵͽ�300��λ���߶�1081.6��ʱ��������0.6��
	double d_aK=0.0*MIL;//��б����3500���ߵͽ�300��λ���߶�1081.6��ʱ��������0.2��
	double d_tv0=0.025;//��б����3500���ߵͽ�300��λ���߶�1081.6��ʱ��������0.065��
	double d_tK=0.006;//��б����3500���ߵͽ�300��λ���߶�1081.6��ʱ��������0.026��
	double d_tD=0;

	//dandaozhuyuan(H,Dq,V0,&alfa, &omg,&vq,&tf,&sita,&fai,&d_aF,&d_tF,&d_ZF);
	//dandaopiancha(H,Dq,V0,&d_av0,&d_aK,&d_tv0,&d_tK,&d_tD);

	double vxd,  s1, s2, s3, eVx1, eVx2;
	Fangxiangyuxian(vq,Vm,sita,qq,landa,H,Dq,&vxd,&s1,&s2,&s3,&eVx1,&eVx2);

	double h_xg,s_xg,w_xg;//������������ϵ��
	Xiangguanxishu(hk,sd,wd,n,h,&h_xg,&s_xg,&w_xg);

	//����Ŀ��ͶӰ���
	double ST=Touying(SB,SS,SF,s1,s2,s3, sita,landa,qq, epsq, vq, Vm);
	double l=sqrt(ST)/2;
	
	//�������������
	double b_xs[4],h_xs[4],s_xs[4],w_xs[4],g_xs[4];
	B_xg(sig_N,sig_Z,sig_V0p,d_av0, d_tv0, Dq, dq,Vm, eVx1, eVx2, &b_xs[0],&b_xs[1],&b_xs[2],&b_xs[3]);
	R_xgh(sigh_f,sigh_b,Dq,dq,&h_xs[0],&h_xs[1],&h_xs[2],&h_xs[3]);
	R_xgs(sigs_f,sigs_b,Dq,dq,&s_xs[0],&s_xs[1],&s_xs[2],&s_xs[3]);
	R_xgw(sigw_f,sigw_b,Dq,dq,&w_xs[0],&w_xs[1],&w_xs[2],&w_xs[3]);
	G_xg(sig_V0,sig_K,sig_Fx,sig_Fz,Dq,dq,Vm,d_av0,d_aK,d_aF,d_tv0,d_tK,d_tF,d_ZF,eVx1,eVx2,&g_xs[0],&g_xs[1],
		&g_xs[2],&g_xs[3]);
	
	double f_sig11,f_sig12,f_sig21,f_sig22,f_XieFC[4];//���ظ�Э�������
	double c_sig11,c_sig12,c_sig21,c_sig22,c_XieFC[4];//�ظ�Э�������
	F_xiefangcha(b_xs,h_xs,s_xs,w_xs,h_xg,s_xg,w_xg,&f_sig11,&f_sig12,&f_sig21,&f_sig22);
	f_XieFC[0]=f_sig11,f_XieFC[1]=f_sig12,f_XieFC[2]=f_sig21,f_XieFC[3]=f_sig22;

	C_xiefangcha(g_xs,h_xs,s_xs,w_xs,h_xg,s_xg,w_xg,&c_sig11,&c_sig12,&c_sig21,&c_sig22); 
	c_XieFC[0]=c_sig11,c_XieFC[1]=c_sig12,c_XieFC[2]=c_sig21,c_XieFC[3]=c_sig22;

	//������ֱ任��Ĳ���
	double sita1,  sita2,  lc1,  lc2,  lf1,  lf2;
	Jifenbianhuan(f_XieFC,c_XieFC,&sita1,&sita2,&lc1,&lc2,&lf1,&lf2);

	//����ϵͳ���
	double a1,a2;
	Xitong_wucha(Dq,dq,af,abq,&a1,&a2);

	int N1=50,N2=50;
	K=1-Xinpusen_jifen(m,p,n,l,w,lc1,lc2,lf1,lf2,sita1,sita2,a1,a2,N1,N2); 

	return K;
}

//�Ż��ٶȽ�����ɭ�����еĽ���Ƴ�һ�ű�
void Init_table(double m, double p, double n, double l, double w, 
				   double lc1,double lc2,double lf1,double lf2, double sita1,double sita2,
				   double a1,double a2, double yc1, double yc2, int N1,int N2)
{
}
