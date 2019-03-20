#include <iostream>
#include <math.h>
#include <vector>
#include <numeric>

using namespace std;

/*函数get_A_para得到方程的参数A
* parameters: 
*  v0: 浮点型变量，为弹丸初速度 
*/ 
float get_A_para(float v0)
{
	float cs = 341.2; //地面标准速度，可以更改 
	float A = 11.260*(cs/v0) + 10.564*atan(2.3330*(cs/v0)-3.8846) + 7.280*log(pow(cs/v0,2)-3.3302*(cs/v0)+2.9562);
	return A;
}

/*参数p的求解方式 p=4.737 * 10**(-4) * cH(y)，为弹道系数（测试为0.003左右接近文章） 
* 函数AccelerDifferEquation根据迭代公式计算v，并返回浮点型v值 
* parameters:   
* a0: p的初始迭代值     c:地面标准速度    t: 任意时刻，可接受浮点型    v0: 弹丸初速度
* mode:可选择两种迭代结束方式，true为固定迭代次数，false为固定相对误差
* time: 迭代次数，初始设定为100    error：相对误差，初始设定为0.01 
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
	return 341.2/q; //空气标准速度 
}

/*函数DispHorizon计算水平位移，返回浮点型变量 
* parameters:
* delta_t: 时间步长   t: 任意时刻（必须为delta_t的整数倍） 
* mode, a0, c, p, v0, time, error参见AccelerDifferEquation函数 
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

