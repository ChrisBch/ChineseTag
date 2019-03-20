#include <iostream>
#include <math.h>
#include <vector>
#include <numeric>

using namespace std;

/*����get_A_para�õ����̵Ĳ���A
* parameters: 
*  v0: �����ͱ�����Ϊ������ٶ� 
*/ 
float get_A_para(float v0)
{
	float cs = 341.2; //�����׼�ٶȣ����Ը��� 
	float A = 11.260*(cs/v0) + 10.564*atan(2.3330*(cs/v0)-3.8846) + 7.280*log(pow(cs/v0,2)-3.3302*(cs/v0)+2.9562);
	return A;
}

/*����p����ⷽʽ p=4.737 * 10**(-4) * cH(y)��Ϊ����ϵ��������Ϊ0.003���ҽӽ����£� 
* ����AccelerDifferEquation���ݵ�����ʽ����v�������ظ�����vֵ 
* parameters:   
* a0: p�ĳ�ʼ����ֵ     c:�����׼�ٶ�    t: ����ʱ�̣��ɽ��ܸ�����    v0: ������ٶ�
* mode:��ѡ�����ֵ���������ʽ��trueΪ�̶�����������falseΪ�̶�������
* time: ������������ʼ�趨Ϊ100    error���������ʼ�趨Ϊ0.01 
*/ 
float AccelerDifferEquation(float a0, float c, float p, float t, float v0, bool mode, int time=100, float error=0.01)
{
	float A = get_A_para(v0);
	float q = ((-10.564)*atan(2.3330*a0-3.8846)-7.2870*log(pow(a0,2)-3.3302*a0+2.9562)+c*p*t+A)/11.260;
	if (mode == true)
	{
		for(int i=0;i<time;i++)
		{
			q = ((-10.564)*atan(2.3330*q-3.8846)-7.2870*log(pow(q,2)-3.3302*q+2.9562)+c*p*t+A)/11.260;
		}
	}
	else
	{
		float new_q;
		while(1)
	    {
		    new_q = ((-10.564)*atan(2.3330*q-3.8846)-7.2870*log(pow(q,2)-3.3302*q+2.9562)+c*p*t+A)/11.260;
		    if((abs(new_q-q)/q)<=error)
		        break;
		    q = new_q;
	    }
	}
	return 341.2/q; //������׼�ٶ� 
}

/*����DispHorizon����ˮƽλ�ƣ����ظ����ͱ��� 
* parameters:
* delta_t: ʱ�䲽��   t: ����ʱ�̣�����Ϊdelta_t���������� 
* mode, a0, c, p, v0, time, error�μ�AccelerDifferEquation���� 
*/
float DispHorizon(float delta_t, float t, bool mode, float a0, float c, float p, float v0, int time=100, float error=0.01)
{
	vector<float> mini_disp;
	for(int j=0; j<=t/delta_t; j++)
		mini_disp.push_back(AccelerDifferEquation(a0,c,p,(j+1)*delta_t,v0,mode,time,error)*delta_t);
	vector<float>(mini_disp).swap(mini_disp);
	float x = ((t/delta_t)/((t/delta_t)+1))*accumulate(mini_disp.begin(),mini_disp.end(),0.0);
	return x;
 } 
 
int main()
{
	float x = DispHorizon(0.001,0.233,true,0.5,341.2,0.003,890,50);
	cout<<x<<endl;
	return 0;
}

