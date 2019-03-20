#include <stdio.h> 
#include <math.h>

double CHD(double d, double i_43, double G, double H, double D_q, double C_43)
/*
* ����CHD�������CHD������double����ֵ 
* parameters 
* d: ����ֱ��  i_43: 1943���������ɵ���ϵ��   G: ��������   
* H: ��ǰ�㴦��Ŀ��߶�  D_q: ��ǰ�㴦��Ŀ��б���� C_43: ����ϵ��
* Note: �����֪C_43��ֵ�����뺯���У���d,i_43,G����ֵ��Ϊ��Ч 
*/
{
	double h,result;
	h = pow(1-(2.1905e-5)*H*0.5, 4.4);
	if (C_43 == 0.0)
	    C_43 = (i_43*pow(d,2)*1000)/G;
	result = C_43*h*D_q;
	return result;
}


double func_yita(double d, double i_43, double G, double H, double D_q, double v0, double C_43)
/*
* ����func_yitaΪ����g_yita����Ϻ�����һ���֣�����ֵ�ᱻ����Ϻ����õ� 
* parameters 
* d: ����ֱ��  i_43: 1943���������ɵ���ϵ��   G: ��������   
* H: ��ǰ�㴦��Ŀ��߶�  D_q: ��ǰ�㴦��Ŀ��б���� v0: ������ٶ� C_43: ����ϵ��
* Note: �����֪C_43��ֵ�����뺯���У���d,i_43,G����ֵ��Ϊ��Ч 
*/ 
{
	double C_HD, fir, sec, result;
	C_HD = CHD(d, i_43, G, H, D_q, C_43);
	fir = 0.2347343+((0.4406926e-6)-(0.2835580e-13)*C_HD*v0)*C_HD*v0-((0.3241693e-8)-(0.3199712e-10)*v0+(0.1166426e-11)*C_HD\
	-(0.2123943e-14)*C_HD*v0+(0.5880930e-16)*pow(C_HD,2))*pow(C_HD,2);
    sec = ((0.1443440e-5)+(0.6348184e-9)*C_HD-(0.1934912e-8)*v0-(0.2794775e-12)*C_HD*v0+(0.7250672e-12)*pow(v0,2))*pow(v0,2);
    result = exp(fir-sec);
    return result;
}


double G_yita(double d, double i_43, double G, double v0, double H, double D_q, double C_43=0.0)
/*
* ����G_yitaΪg_yita����Ϻ���������ֵΪdouble 
* parameters 
* d: ����ֱ��  i_43: 1943���������ɵ���ϵ��   G: ��������   
* H: ��ǰ�㴦��Ŀ��߶�  D_q: ��ǰ�㴦��Ŀ��б���� v0: ������ٶ� C_43: ����ϵ��
* Note: �����֪C_43��ֵ�����뺯���У���d,i_43,G����ֵ��Ϊ��Ч 
*/
{
	double fir, sec, result;
	result = 0.0;
	if((H>4000)&&(D_q>6500))
	{
		result = func_yita(d, i_43, G, H, D_q,v0, C_43)*pow(D_q/2000.0, 0.075)*pow(H/1500.0, 0.115)*pow(D_q/6500.0, 0.45);
	}
	else if((H>1500)&&(D_q>2000))
	{
		result = func_yita(d, i_43, G, H, D_q,v0, C_43)*pow(D_q/2000.0, 0.075)*pow(H/1500.0, 0.115);
	}
	else if((H>1500)&&(D_q<=2000))
	{
		result = func_yita(d, i_43, G, H, D_q,v0, C_43)*pow(H/1500.0, 0.115);
	}
	else if((H>800)&&(D_q>2000))
	{
		result = func_yita(d, i_43, G, H, D_q,v0, C_43)*pow(D_q/2000.0, 0.075);
	}
	else if((H<=800)||((H>800)&&(H<=1500)&&(D_q<=2000)))
	{
		result = func_yita(d, i_43, G, H, D_q,v0, C_43);
	}
	return result;	
}


double C_r(double d_q, double H, double D_q, double d, double i_43, double G, double C_43=0.0)
/*
* ����C_rΪ����g_omega����Ϻ�����һ���֣�����ֵ�ᱻ����Ϻ����õ� 
* parameters 
* d: ����ֱ��  i_43: 1943���������ɵ���ϵ��   G: ��������   d_q: Ŀ����ǰ���ˮƽ����
* H: ��ǰ�㴦��Ŀ��߶�  D_q: ��ǰ�㴦��Ŀ��б���� C_43: ����ϵ��
* Note: �����֪C_43��ֵ�����뺯���У���d,i_43,G����ֵ��Ϊ��Ч 
*/
{
	double C_HD, result;
	C_HD = CHD(d, i_43, G, H, D_q, C_43);
	result = 0.0;
	if ((d_q>5000)&&(H>=1500))
	{
		result = C_HD*pow(d_q/2500.0, 0.085)*pow(D_q/5000.0, 0.085)*pow(H/1500.0, 0.18);
	}
	else if((d_q>3000)&&(d_q<=5000)&&(H>=1500)&&(D_q>=4000))
	{
		result = C_HD*pow(d_q/2500.0, 0.085)*pow(D_q/5000.0, 0.045)*pow(H/1500.0, 0.14);
	}
	else if(d_q>2500)
	{
		result = C_HD*pow(d_q/2500.0, 0.085);
	}
	else if(d_q<=2500)
	{
		result = C_HD;
	}
	return result;
}


double G_omega(double d, double i_43, double G, double H,double v0, double d_q, double D_q, double C_43=0.0)
/*
* ����G_omegaΪg_omega����Ϻ���������ֵΪdouble 
* parameters 
* d: ����ֱ��  i_43: 1943���������ɵ���ϵ��   G: ��������  d_q: Ŀ����ǰ���ˮƽ���� 
* H: ��ǰ�㴦��Ŀ��߶�  D_q: ��ǰ�㴦��Ŀ��б���� v0: ������ٶ� C_43: ����ϵ��
* Note: �����֪C_43��ֵ�����뺯���У���d,i_43,G����ֵ��Ϊ��Ч 
*/
{
	double fir, sec, result, Cr;
	Cr = C_r(d_q, H, D_q, d, i_43, G);
	fir = (-0.7034540e-3)+(0.2152844e-3)*Cr+((0.1660815e-6)-(0.4358138e-9)*v0+(0.2247045e-12)*pow(v0,2)-(0.4583237e-13)*Cr*v0)*Cr*v0;
	sec = ((0.7358717e-8)-(0.5559608e-10)*v0+(0.2756331e-11)*Cr-(0.4376855e-14)*Cr*v0+(0.1413296e-15)*pow(Cr,2))*pow(Cr,2);
	result = exp(fir-sec);
	return result;
}


double G_v(double d, double i_43, double G, double H, double D_q,double v0, double C_43=0.0)
/*
* ����G_vΪg_v����Ϻ���������ֵΪdouble 
* parameters 
* d: ����ֱ��  i_43: 1943���������ɵ���ϵ��   G: �������� 
* H: ��ǰ�㴦��Ŀ��߶�  D_q: ��ǰ�㴦��Ŀ��б���� v0: ������ٶ� C_43: ����ϵ��
* Note: �����֪C_43��ֵ�����뺯���У���d,i_43,G����ֵ��Ϊ��Ч 
*/
{
	double C_HD, fir, sec, result;
	C_HD = CHD(d, i_43, G, H, D_q, C_43);
	fir = (-0.6546777e-1)-((0.1677301e-3)+(0.7276473e-7)*v0-(0.9811960e-8)*C_HD+(0.2390793e-10)*C_HD*v0-(0.1145990e-11)*pow(C_HD,2))*C_HD;
	sec = ((0.2952592e-6)+(0.1309007e-9)*C_HD-(0.2306377e-9)*v0)*pow(v0,2);
	result = exp(fir+sec);
	return result;
}


double C_t(double D_q, double d, double i_43, double G, double H, double C_43=0.0)
/*
* ����C_tΪ����g_t����Ϻ�����һ���֣�����ֵ�ᱻ����Ϻ����õ� 
* parameters 
* d: ����ֱ��  i_43: 1943���������ɵ���ϵ��   G: ��������   
* H: ��ǰ�㴦��Ŀ��߶�  D_q: ��ǰ�㴦��Ŀ��б���� C_43: ����ϵ��
* Note: �����֪C_43��ֵ�����뺯���У���d,i_43,G����ֵ��Ϊ��Ч 
*/
{
	double C_HD, result;
	C_HD = CHD(d, i_43, G, H, D_q, C_43);
	result = C_HD;
	if (H>3500)
	{
		result = result*pow(D_q/3500.0, 0.018);
	}
	return result;
}


double G_t(double D_q, double d, double i_43, double G, double H, double v0, double C_43=0.0)
/*
* ����G_tΪg_t����Ϻ���������ֵΪdouble 
* parameters 
* d: ����ֱ��  i_43: 1943���������ɵ���ϵ��   G: �������� 
* H: ��ǰ�㴦��Ŀ��߶�  D_q: ��ǰ�㴦��Ŀ��б���� v0: ������ٶ� C_43: ����ϵ��
* Note: �����֪C_43��ֵ�����뺯���У���d,i_43,G����ֵ��Ϊ��Ч 
*/
{
	double fir, sec, result, Ct;
	Ct = C_t(D_q, d, i_43, G, H);
	fir = 0.1658128+((0.3163075e-6)-(0.1850771e-13)*Ct*v0)*Ct*v0-((0.3603319e-8)-(0.2294901e-10)*v0+(0.9587002e-12)*Ct-(0.1503861e-14)*Ct*v0\
	+(0.4369190e-16)*pow(Ct,2))*pow(Ct,2);
	sec = ((0.9955841e-6)+(0.4341479e-9)*Ct-(0.1309399e-8)*v0-(0.1810468e-12)*Ct*v0+(0.4795298e-12)*pow(v0,2))*pow(v0,2);
	result = exp(fir-sec);
	return result;
}


double compute_alpha(double g, double e_q, double d, double i_43, double G, double v0, double H, double D_q, double C_43=0.0)
/*
* ����compute_alphaΪ����߽Ǻ���������ֵΪdouble 
* parameters 
* d: ����ֱ��  i_43: 1943���������ɵ���ϵ��   G: �������� 
* H: ��ǰ�㴦��Ŀ��߶�  D_q: ��ǰ�㴦��Ŀ��б���� v0: ������ٶ� C_43: ����ϵ��
* e_q: ��ǰ�㴦��Ŀ��ߵͽ� g: �������ٶ� 
* Note: �����֪C_43��ֵ�����뺯���У���d,i_43,G����ֵ��Ϊ��Ч 
*/
{
	double result;
	result = (cos(e_q)*g*D_q*G_yita(d, i_43, G, v0, H, D_q))/(2*pow(v0,2));
	return result;
}


double compute_omega(double g, double e_q, double d, double i_43, double G, double v0,double H, double d_q, double D_q, double C_43=0.0)
/*
* ����compute_omegaΪ�����ŽǺ���������ֵΪdouble 
* parameters 
* d: ����ֱ��  i_43: 1943���������ɵ���ϵ��   G: �������� 
* H: ��ǰ�㴦��Ŀ��߶�  D_q: ��ǰ�㴦��Ŀ��б���� v0: ������ٶ� C_43: ����ϵ��
* e_q: ��ǰ�㴦��Ŀ��ߵͽ� g: �������ٶ� d_q: Ŀ����ǰ���ˮƽ����
* Note: �����֪C_43��ֵ�����뺯���У���d,i_43,G����ֵ��Ϊ��Ч 
*/
{
	double result;
	result = (cos(e_q)*g*D_q*G_omega(d, i_43, G, v0, H, d_q, D_q))/(2*pow(v0,2));
	return result;
}


double compute_V_p(double g, double e_q, double d, double i_43, double G, double v0, double H, double D_q, double d_q, double C_43=0.0)
/*
* ����compute_V_pΪ������ǰ�㴦�ĵ�����ٵĺ���������ֵΪdouble 
* parameters 
* d: ����ֱ��  i_43: 1943���������ɵ���ϵ��   G: �������� 
* H: ��ǰ�㴦��Ŀ��߶�  D_q: ��ǰ�㴦��Ŀ��б���� v0: ������ٶ� C_43: ����ϵ��
* e_q: ��ǰ�㴦��Ŀ��ߵͽ� g: �������ٶ� d_q: Ŀ����ǰ���ˮƽ����
* Note: �����֪C_43��ֵ�����뺯���У���d,i_43,G����ֵ��Ϊ��Ч 
*/
{
	double theta, fai, result;
	theta = e_q - compute_omega(g, e_q, d, i_43, G, v0, H, d_q, D_q);
	fai = e_q + compute_alpha(g, e_q, d, i_43, G, v0, H, D_q);
	result = (v0*cos(fai)*G_v(d, i_43, G, H, D_q, v0))/cos(theta);
	return result;
}


double compute_t_f(double g, double e_q, double d, double i_43, double G, double v0, double H, double D_q, double C_43=0.0)
/*
* ����compute_t_fΪ���㵯�����ʱ��ĺ���������ֵΪdouble 
* parameters 
* d: ����ֱ��  i_43: 1943���������ɵ���ϵ��   G: �������� 
* H: ��ǰ�㴦��Ŀ��߶�  D_q: ��ǰ�㴦��Ŀ��б���� v0: ������ٶ� C_43: ����ϵ��
* e_q: ��ǰ�㴦��Ŀ��ߵͽ� g: �������ٶ� d_q: Ŀ����ǰ���ˮƽ����
* Note: �����֪C_43��ֵ�����뺯���У���d,i_43,G����ֵ��Ϊ��Ч 
*/
{
	double fai, result;
	fai = e_q + compute_alpha(g, e_q, d, i_43, G, v0, H, D_q);
	result = (cos(e_q)*D_q*G_t(D_q, d, i_43, G, H, v0))/(cos(fai)*v0);
	return result;
}


double delta_alpha_f(double g, double e_q, double d, double i_43, double G, double v0, double H, double D_q,double C_43 = 0.0)
/*
* ����delta_alpha_fΪ�����ݷ�1m/s����ĸ߽����ĺ���������ֵΪdouble 
* parameters 
* d: ����ֱ��  i_43: 1943���������ɵ���ϵ��   G: �������� 
* H: ��ǰ�㴦��Ŀ��߶�  D_q: ��ǰ�㴦��Ŀ��б���� v0: ������ٶ� C_43: ����ϵ��
* e_q: ��ǰ�㴦��Ŀ��ߵͽ� g: �������ٶ�
* Note: �����֪C_43��ֵ�����뺯���У���d,i_43,G����ֵ��Ϊ��Ч 
*/
{
	double fai, result, fir,sec;
	fai = e_q + compute_alpha(g, e_q, d, i_43, G, v0, H, D_q);
	fir = (0.2595097e-2)-((0.1444118e-5)+(0.1348652e-7)*fai-(0.1775943e-9)*pow(fai,2)-(0.1877511e-11)*D_q*fai)*D_q*fai;
	sec = ((0.5699964e-9)+(0.6742182e-13)*D_q*fai-(0.4077875e-15)*pow(D_q,2))*pow(D_q,2);
	result = pow((960*C_43)/(1.508*v0),1.2)*(fir-sec);
	return result;
}


double delta_t_F(double g, double e_q, double d, double i_43, double G, double v0, double H, double D_q,double C_43=0.0)
/*
* ����delta_t_FΪ�����ݷ�1m/s����ĵ������ʱ�����ĺ���������ֵΪdouble 
* parameters 
* d: ����ֱ��  i_43: 1943���������ɵ���ϵ��   G: �������� 
* H: ��ǰ�㴦��Ŀ��߶�  D_q: ��ǰ�㴦��Ŀ��б���� v0: ������ٶ� C_43: ����ϵ��
* e_q: ��ǰ�㴦��Ŀ��ߵͽ� g: �������ٶ�
* Note: �����֪C_43��ֵ�����뺯���У���d,i_43,G����ֵ��Ϊ��Ч 
*/
{
	double fai, result, fir,sec;
	fai = e_q + compute_alpha(g, e_q, d, i_43, G, v0, H, D_q);
	fir = (0.1655924e-3)-((0.2358256e-7)+(0.1893686e-12)*D_q*fai)*D_q*fai+((0.1388580e-10)*fai+(0.1150672e-12)*D_q-(0.1789927e-14)*D_q*fai+(0.9054865e-17)*pow(D_q,2))*pow(D_q,2);
	sec = ((0.3999206e-6)-(0.6522549e-8)*fai+(0.5080177e-11)*D_q*fai)*pow(fai,2);
	result = pow((960*C_43)/(1.508*v0),2.55)*(fir+sec);
	return result;
}


double delta_Z_F(double H, double d_q, double v0)
/*
* ����delta_Z_FΪ������1m/s����ķ������ĺ���������ֵΪdouble 
* parameters 
* H: ��ǰ�㴦��Ŀ��߶�  d_q: Ŀ����ǰ���ˮƽ���� v0: ������ٶ�
*/
{
	double fir, sec, result;
	fir = (0.8966708e-1)-(0.1648553e-3)*d_q+((0.8872561e-7)+(0.1133668e-10)*d_q+(0.9842440e-11)*H)*pow(H,2);
	sec = ((0.1667704e-6)-(0.1934115e-11)*H+(0.7819438e-11)*d_q)*pow(d_q,2);
	return result;
}
