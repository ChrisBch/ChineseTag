
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.141592654

//gaussrand�������ݸ����ľ�ֵ�ͷ�����ط�����̬�ֲ���һ������
//parameters: 
//mean: float��������ֵ��stdiv��float���������� 
//  
float gaussrand(float mean,float stdiv)
{
    static float U, V;
    float Z;
    float X;
    U = rand() / (RAND_MAX + 1.0);
    V = rand() / (RAND_MAX + 1.0);
    Z = sqrt(-2.0 * log(U))* sin(2.0 * PI * V);
    X = mean+Z*stdiv;
    return X;
}  
 
 //change this function to change g(x) 
//�����ϲ�ȷ���ȵĹ��ϸ��ʵĺ���������ֵΪ0�����ڹ��ϣ���1��û�й��ϣ� 
//parameters��
//sample��float���� m��yi��float������Ϊ����һ�������Ĳ��� 
//
int omega_f_mix(float X,float Y)
{
	int g, x1=1, x2=3, y1=1, y2=3; //��� 
	if ((X>=x1)&&(X<=x2)&&(Y>=y1)&&(Y<=y2))
	    g = 1;
	else
	    g = 0;
	return g;
 }

 float A_neg(float c, int Mx, int My, float mean, float stdiv, float m, float yi)
/*
 * Mx, My�ֱ�Ϊ������֪�����ĸ��ʲ�ȵĲ������Լ�������֪��ȷ��������assigned nominal probability measure�Ĳ����� 
 * mean,stdiv�ֱ�Ϊ��֪��̬�ֲ��ı����Ĳ��� 
 * m��yiΪnominal probability measure�Ĳ��� 
 */
 {
 
 	float samples_x[Mx];
 	float samples_y[My];
 	double sum_samples_inside = 0.0;
 	double sum_samples_outside = 0.0;
 	for(int i=0; i<Mx; i++)
 	    samples_x[i]=gaussrand(mean,stdiv);
 	for(int i=0; i<My; i++)
 	    samples_y[i]=gaussrand(m,yi);
 	//for(int i=0; i<My; i++)
 	//    samples_y[i]=
 	for(int i=0; i<My; i++)
 	{
 		for(int j=0; j<Mx; j++)
 		{
 			sum_samples_inside = sum_samples_inside + (-1)*c*omega_f_mix(samples_x[j],samples_y[i]);
		 }
		 sum_samples_inside = exp(sum_samples_inside/Mx);
		 sum_samples_outside = sum_samples_outside + sum_samples_inside;
		 sum_samples_inside = 0.0;
	 }
	 sum_samples_outside = (-1)*log(sum_samples_outside/My)/c;
	 return sum_samples_outside;
 }
 
 float A_pos(float c, int Mx, int My, float mean, float stdiv, float m, float yi)
 /*
 * Mx, My�ֱ�Ϊ������֪�����ĸ��ʲ�ȵĲ������Լ�������֪��ȷ��������assigned nominal probability measure�Ĳ����� 
 * mean,stdiv�ֱ�Ϊ��֪��̬�ֲ��ı����Ĳ��� 
 * m��yiΪnominal probability measure�Ĳ��� 
 */
 {
 	float samples_x[Mx];
 	float samples_y[My];
 	double sum_samples_inside = 0.0;
 	double sum_samples_outside = 0.0;
 	for(int i=0; i<Mx; i++)
 	    samples_x[i]=gaussrand(0.0,1.0);
 	for(int i=0; i<My; i++)
 	    samples_y[i]=gaussrand(m,yi);
 	for(int i=0; i<My; i++)
 	{
 		for(int j=0; j<Mx; j++)
 		{
 			sum_samples_inside = sum_samples_inside + c*omega_f_mix(samples_x[j],samples_y[i]);
		 }
		 sum_samples_inside = exp(sum_samples_inside/Mx);
		 sum_samples_outside = sum_samples_outside + sum_samples_inside;
		 sum_samples_inside = 0.0;
	 }
	 sum_samples_outside = log(sum_samples_outside/My)/c;
	 return sum_samples_outside;
 }
 
 /*
 double theta_pos(float c, int M, float mean, float stdiv)
 {
	 double samples[M];
	 for(int i=0; i<M; i++)
	 	samples[i]=gaussrand(mean,stdiv);
 	 double sum_sample = 0.0;
 	 for(int i=0;i<M;i++)
 	 	sum_sample = exp(c*omega_f_epi(samples[i])) + sum_sample;
 	 double result = (log(sum_sample/M))/c;
	 return result;
 }

 double theta_neg(float c, int M, float mean, float stdiv)
 {
	 double samples[M];
	 for(int i=0; i<M; i++)
	 	samples[i]=gaussrand(mean,stdiv);
 	 double sum_sample = 0.0;
 	 for(int i=0;i<M;i++)
 	 	sum_sample = exp((-1)*c*omega_f_epi(samples[i])) + sum_sample;
 	 double result = ((-1)*log(sum_sample/M))/c;
	 return result;
 }
*/

 float GaussianDistribution(float mean, float stdiv, float x)
 {
 	return exp(pow(x-mean,2)/(2*pow(stdiv,2)))/(sqrt(2*PI)*stdiv);
 }

 float relative_entropy(float mean_1, float mean_2, float stdiv_1, float stdiv_2)
 {
 	float IncreaseRate = 0.001;
 	float count = 0.0;
 	float low_bound,up_bound,Gaussian_1,Gaussian_2;
 	if ((mean_1-2.58*stdiv_1)>=(mean_2-2.58*stdiv_2))
 	    low_bound = mean_2-2.58*stdiv_2;
 	else
 	    low_bound = mean_1-2.58*stdiv_1;
 	if ((mean_1+2.58*stdiv_1)>=(mean_2+2.58*stdiv_2))
 	    up_bound = mean_1+2.58*stdiv_1;
 	else
 	    up_bound = mean_2+2.58*stdiv_2;
 	for(float i=low_bound;i<up_bound;i=i+IncreaseRate) //change here
 	{
 		Gaussian_1 = GaussianDistribution(mean_1,stdiv_1,i);
 		Gaussian_2 = GaussianDistribution(mean_2,stdiv_2,i);
 		count = Gaussian_1*log(Gaussian_1/Gaussian_2)*IncreaseRate + count;
	 }
	return count;
 }

float sup_relative_entropy(float mean_range_left, float mean_range_right, float stdiv_range_left, float stdiv_range_right, float mean, float stdiv)
{
 	float RateMean = 0.001;
 	float RateStdiv = 0.001;
 	float TempSup = relative_entropy(mean_range_left,mean,stdiv_range_right,stdiv);
 	if (mean_range_left==mean_range_right)
 	{
 	    for(float j=stdiv_range_left; j<stdiv_range_right; j=j+RateStdiv)
 		{
 		  if (relative_entropy(mean,mean,j,stdiv)<TempSup)
 		 	 TempSup =  relative_entropy(mean,mean,j,stdiv);
 		}
 	}
 	else
 	{
 	    for(float i=mean_range_left;i<mean_range_right;i=i+RateMean)
 		{
 	    	if (stdiv_range_left==stdiv_range_right)
 	    	{
 	    	    if (relative_entropy(i,mean,stdiv,stdiv)<TempSup)
 	    	    {
 	    	    	TempSup =  relative_entropy(i,mean,stdiv,stdiv);
 	    	    }
 	    	}
 	    	else
 	    	{
 	    	   for(float j=stdiv_range_left; j<stdiv_range_right; j=j+RateStdiv)
 	    	   {
 	    		 if (relative_entropy(i,mean,j,stdiv)<TempSup)
 	    		 	TempSup =  relative_entropy(i,mean,j,stdiv);
 	    	   }
 	    	}
 		}
 	}
   	return TempSup;
}
/*
//mode1 is upper bound, and mode2 is lower bound
float EpiUncer(float c, int M, int mode,float mean_range_left, float mean_range_right, float stdiv_range_left, float stdiv_range_right, float mean, float stdiv)
{
	float R_star;
	float theta;
	float result;
	if (mode==1)
	{
		R_star = sup_relative_entropy(mean_range_left, mean_range_right, stdiv_range_left, stdiv_range_right, mean, stdiv);
		theta = theta_pos(c,M,mean,stdiv);
		result = theta + (R_star/c);
	}
	else
	{
		R_star = sup_relative_entropy(mean_range_left, mean_range_right, stdiv_range_left, stdiv_range_right, mean, stdiv);
		theta = theta_neg(c,M,mean,stdiv);
		result = theta - (R_star/c);
	}
	return result;
}
*/

float MixEpiAle(float c, int Mx, int My, int mode,float mean_range_left, float mean_range_right, float stdiv_range_left, float stdiv_range_right, float mean_known, float stdiv_known, float mean_nomi, float stdiv_nomi)
{
/*
* cΪ������Խ����Խ�����ڱ�׼��mode=1ʱ������ȷ�磬����ʱ������ȷ�� 
* mean_left,mean_right,stdiv_left,stdiv_rightΪ��֪��ȷ�������Ĳ����仯��Χ���˴�Ĭ�ϸñ���Ϊ��̬�ֲ�������ͨ���趨����������ͬ��ֵ�ķ�ʽ��ֻ���� 
* һ�������仯������� mean_known, stdiv_known����֪�ĸ��ʷֲ��Ĳ��� 
* mean_nomi, stdiv_nomi��nominal probability measure�Ĳ��� 
*/
	float R_star;
	float A;
	float result;
	if (mode==1)
	{
		R_star = sup_relative_entropy(mean_range_left, mean_range_right, stdiv_range_left, stdiv_range_right, mean_nomi, stdiv_nomi);
		A = A_pos(c, Mx, My, mean_known, stdiv_known, mean_nomi, stdiv_nomi);
		result = A + R_star;
	}
	else
	{
	    R_star = sup_relative_entropy(mean_range_left, mean_range_right, stdiv_range_left, stdiv_range_right, mean_nomi, stdiv_nomi);
		A = A_neg(c, Mx, My, mean_known, stdiv_known, mean_nomi, stdiv_nomi);
		result = A - R_star;	
	}
	return result;
}

int main()
{
	for(int i=0; i<=10; i++)
	    printf("%f \n",MixEpiAle(50,100,100,1,0.0,1.0,0.0,0.0,2.0,1.0,2.0,2.0));
 	return 0;
}



