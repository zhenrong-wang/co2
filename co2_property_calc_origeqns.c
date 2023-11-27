/*
 *This code is distributed under the license: MIT License
 * Originally written by Zhenrong WANG
 * mailto: zhenrongwang@live.com
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "co2_property_calc.h"

const double coeff0[8][2]={
	8.37304456,	0,
	-3.70454304,0,
	2.5000000,0,
	1.99427042,3.15163,
	0.62105248,6.11190,
	0.41195293,6.77708,
	1.04028922,11.32384,
	0.08327678,27.08792,
};

const double coeffr1[7][3]={
	0.38856823203161,1,0,
	0.2938547594274e1,1,0.75,
	-0.55867188534934e1,1,1,
	-0.76753199592477,1,2,
	0.31729005580416,2,0.75,
	0.54803315897767,2,2,
	0.12279411220335,3,0.75,
};

const double coeffr2[27][4]={
	0.21658961543220e1,1,1.5,1,
	0.15841735109724e1,2,1.5,1,
	-0.23132705405503,4,2.5,1,
	0.58116916431436e-1,5,0,1,
	-0.55369137205382,5,1.5,1,
	0.48946615909422,5,2,1,
	-0.24275739843501e-1,6,0,1,
	0.62494790501678e-1,6,1,1,
	-0.12175860225246,6,2,1,
	-0.37055685270086,1,3,2,
	-0.16775879700426e-1,1,6,2,
	-0.11960736637987,4,3,2,
	-0.45619362508778e-1,4,6,2,
	0.35612789270346e-1,4,8,2,
	-0.74422727132052e-2,7,6,2,
	-0.17395704902432e-2,8,0,2,
	-0.21810121289527e-1,2,7,3,
	0.24332166559236e-1,3,12,3,
	-0.37440133423463e-1,3,16,3,
	0.14338715756878,5,22,4,
	-0.13491969083286,5,24,4,
	-0.23151225053480e-1,6,16,4,
	0.12363125492901e-1,7,24,4,
	0.21058321972940e-2,8,8,4,
	-0.33958519026368e-3,10,2,4,
	0.55993651771592e-2,4,28,5,
	-0.30335118055646e-3,8,14,6,
};

const double coeffr3[5][7]={
	-0.21365488688320e3,2,1,25,325,1.16,1,
	0.26641569149272e5,2,0,25,300,1.19,1,
	-0.24027212204557e5,2,1,25,300,1.19,1,
	-0.28341603423999e3,3,3,15,275,1.25,1,
	0.21247284400179e3,3,3,20,275,1.22,1,
};

const double coeffr4[3][8]={
	-0.66642276540751,3.5,0.875,0.300,0.700,0.3,10,275,
	0.72608632349897,3.5,0.925,0.3,0.7,0.3,10,275,
	0.55068668612842e-1,3,0.875,0.3,0.7,1,12.5,275,
};

double phi0(double del, double tau)
{
	double res,ai,thi;
	int i;
	res=log(del)+coeff0[0][0]+coeff0[1][0]*tau+coeff0[2][0]*log(tau);
	for(i=3;i<8;i++)
	{
		ai=coeff0[i][0];
		thi=coeff0[i][1];
		res+=ai*log(1-exp(-thi*tau));
	}
	return res;
}

double phi0d(double del)
{
	return 1/del;
}

double phi0dd(double del)
{
	return -1/(del*del);
}

double phi0t(double tau)
{
	double res, ai, thi;
	int i;
	res=coeff0[1][0]+coeff0[2][0]/tau;
	for(i=3;i<8;i++)
	{
		ai=coeff0[i][0];
		thi=coeff0[i][1];
		res+=ai*thi*(1/(1-exp(-thi*tau))-1);
	}
	return res;
}

double phi0tt(double tau)
{
	double res,ai,thi;
	int i;
	res=-coeff0[2][0]/(tau*tau);
	for(i=3;i<8;i++)
	{
		ai=coeff0[i][0];
		thi=coeff0[i][1];
		res-=ai*thi*thi*exp(-thi*tau)/((1-exp(-thi*tau))*(1-exp(-thi*tau)));
	}
	return res;
}

double deltd(double del, double delta, double aai, double theta, double bti, double bbi, double ai) 
{
	return (del-1)*(aai*theta*2*pow((del-1)*(del-1),0.5/bti-1)/bti+2*bbi*ai*pow((del-1)*(del-1),ai-1));
}

double deltdd(double dd, double del, double delta, double aai, double theta, double bti, double bbi, double ai)
{
//	double dd=deltd(del,delta,aai,theta,bti,bbi,ai);
//	printf("\n\n\n%lf\n",dd);
	return dd/(del-1)+(del-1)*(del-1)*(4*bbi*ai*(ai-1)*pow((del-1)*(del-1),ai-2)+2*aai*aai*pow((del-1)*(del-1),1/bti-2)/(bti*bti)+aai*theta*(4/bti)*(0.5/bti-1)*pow((del-1)*(del-1),0.5/bti-2));
}

double deltabid(double bi, double delta, double dd)
{
	return bi*pow(delta,bi-1)*dd;
}

double deltabidd(double bi, double delta, double dd, double ddd)
{
	return bi*(pow(delta,bi-1)*ddd+(bi-1)*pow(delta,bi-2)*pow(dd,2));
}

double deltabit(double bi, double theta, double delta)
{
	return -2*theta*bi*pow(delta,bi-1);
}

double deltabitt(double bi, double theta, double delta)
{
	return 2*bi*pow(delta,bi-1)+4*theta*theta*bi*(bi-1)*pow(delta,bi-2);
}

double deltabidt(double aai, double bi, double bti, double delta, double del, double theta, double dd)
{
	return -aai*bi*2*pow(delta,bi-1)*(del-1)*pow((del-1)*(del-1),0.5/bti-1)/bti-2*theta*bi*(bi-1)*pow(delta,bi-2)*dd;
}

double psid(double del, double psi, double cci)
{
	return -2*cci*(del-1)*psi;
}

double psidd(double del, double psi, double cci)
{
	return 2*cci*psi*(2*cci*(del-1)*(del-1)-1);
}

double psit(double tau, double psi, double ddi)
{
	return -2*ddi*(tau-1)*psi;
}

double psitt(double tau, double psi, double ddi)
{
	return 2*ddi*psi*(2*ddi*(tau-1)*(tau-1)-1);
}

double psidt(double tau, double del, double psi, double cci, double ddi)
{
	return 4*cci*ddi*(del-1)*(tau-1)*psi;
}

double phir(double del, double tau)
{
	int i;
	double ni,di,ti,ci,afi,bti,gmi,epi,ai,bi,aai,bbi,cci,ddi,delta,theta,psi;
	double res=0;
	for(i=0;i<7;i++)
	{
		ni=coeffr1[i][0];
		di=coeffr1[i][1];
		ti=coeffr1[i][2];
		res+=ni*pow(del,di)*pow(tau,ti);
	}
//	printf("........%lf\n",res);
	
	for(i=0;i<27;i++)
	{
		ni=coeffr2[i][0];
		di=coeffr2[i][1];
		ti=coeffr2[i][2];
		ci=coeffr2[i][3];
		res+=ni*pow(del,di)*pow(tau,ti)*exp(-pow(del,ci));
	}
//		printf("........%lf\n",res);
	for(i=0;i<5;i++)
	{
		ni=coeffr3[i][0];
		di=coeffr3[i][1];
		ti=coeffr3[i][2];
		afi=coeffr3[i][3];
		bti=coeffr3[i][4];
		gmi=coeffr3[i][5];
		epi=coeffr3[i][6];
		res+=ni*pow(del,di)*pow(tau,ti)*exp(-afi*(del-epi)*(del-epi)-bti*(tau-gmi)*(tau-gmi));
	}
//		printf("........%lf\n",res);
	for(i=0;i<3;i++)
	{
		ni=coeffr4[i][0];
		ai=coeffr4[i][1];
		bi=coeffr4[i][2];
		bti=coeffr4[i][3];
		aai=coeffr4[i][4];
		bbi=coeffr4[i][5];
		cci=coeffr4[i][6];
		ddi=coeffr4[i][7];
		
//		printf("************%lf,%lf,%lf,%lf,%lf,%Lf,%lf,%lf\n",ni,ai,bi,bti,aai,bbi,cci,ddi);
		delta=pow(1-tau+aai*pow((del-1)*(del-1),0.5/bti),2)+bbi*pow((del-1)*(del-1),ai);
		psi=exp(-cci*(del-1)*(del-1)-ddi*(tau-1)*(tau-1));
//		printf("*******%lf\n",delta);
		res+=ni*pow(delta,bi)*del*psi;
	}
//		printf("........%lf\n",res);
	return res;
}

double phird(double del, double tau)
{
	int i;
	double ni,di,ti,ci,afi,bti,gmi,epi,ai,bi,aai,bbi,cci,ddi,delta,psi,theta,dd,ddd;
	double res=0;
	for(i=0;i<7;i++)
	{
		ni=coeffr1[i][0];
		di=coeffr1[i][1];
		ti=coeffr1[i][2];
		res+=ni*di*pow(del,di-1)*pow(tau,ti);
	}
	
	for(i=0;i<27;i++)
	{
		ni=coeffr2[i][0];
		di=coeffr2[i][1];
		ti=coeffr2[i][2];
		ci=coeffr2[i][3];
		res+=ni*pow(tau,ti)*exp(-pow(del,ci))*pow(del,di-1)*(di-ci*pow(del,ci));
	}
	
	for(i=0;i<5;i++)
	{
		ni=coeffr3[i][0];
		di=coeffr3[i][1];
		ti=coeffr3[i][2];
		afi=coeffr3[i][3];
		bti=coeffr3[i][4];
		gmi=coeffr3[i][5];
		epi=coeffr3[i][6];
		res+=ni*pow(del,di)*pow(tau,ti)*exp(-afi*(del-epi)*(del-epi)-bti*(tau-gmi)*(tau-gmi))*(di/del-2*afi*(del-epi));
	}
//	printf("..............%lf\n",res);
	for(i=0;i<3;i++)
	{
		ni=coeffr4[i][0];
		ai=coeffr4[i][1];
		bi=coeffr4[i][2];
		bti=coeffr4[i][3];
		aai=coeffr4[i][4];
		bbi=coeffr4[i][5];
		cci=coeffr4[i][6];
		ddi=coeffr4[i][7];
		
		theta=(1-tau)+aai*pow((del-1)*(del-1),0.5/bti);
		delta=theta*theta+bbi*pow((del-1)*(del-1),ai);
		psi=exp(-cci*(del-1)*(del-1)-ddi*(tau-1)*(tau-1));
		dd=deltd(del,delta,aai,theta,bti,bbi,ai);
		ddd=deltdd(dd,del,delta,aai,theta,bti,bbi,ai);
//		printf("\n%lf\t%lf\t%lf\n",theta,delta,psi);
//		printf("\n%lf,%lf,%lf\n",psid(del,psi,cci),deltabid(bi,delta,dd))
		res+=ni*(pow(delta,bi)*(psi+del*psid(del,psi,cci))+deltabid(bi,delta,dd)*del*psi);
	}
//		printf("..............%lf\n\n",res);
	return res;
}

double phirdd(double del, double tau)
{
	int i;
	double ni,di,ti,ci,afi,bti,gmi,epi,ai,bi,aai,bbi,cci,ddi,delta,psi,theta;
	double dd,ddd;
	double res=0;
	for(i=0;i<7;i++)
	{
		ni=coeffr1[i][0];
		di=coeffr1[i][1];
		ti=coeffr1[i][2];
		res+=ni*di*(di-1)*pow(del,di-2)*pow(tau,ti);
	}
	
	for(i=0;i<27;i++)
	{
		ni=coeffr2[i][0];
		di=coeffr2[i][1];
		ti=coeffr2[i][2];
		ci=coeffr2[i][3];
		res+=ni*exp(-pow(del,ci))*(pow(del,di-2)*pow(tau,ti)*((di-ci*pow(del,ci))*(di-1-ci*pow(del,ci))-ci*ci*pow(del,ci)));
	}
	
	for(i=0;i<5;i++)
	{
		ni=coeffr3[i][0];
		di=coeffr3[i][1];
		ti=coeffr3[i][2];
		afi=coeffr3[i][3];
		bti=coeffr3[i][4];
		gmi=coeffr3[i][5];
		epi=coeffr3[i][6];
		res+=ni*pow(tau,ti)*exp(-afi*(del-epi)*(del-epi)-bti*(tau-gmi)*(tau-gmi))*(-2*afi*pow(del,di)+4*afi*afi*pow(del,di)*(del-epi)*(del-epi)-4*di*afi*pow(del,di-1)*(del-epi)+di*(di-1)*pow(del,di-2));
	}
	
	for(i=0;i<3;i++)
	{
		ni=coeffr4[i][0];
		ai=coeffr4[i][1];
		bi=coeffr4[i][2];
		bti=coeffr4[i][3];
		aai=coeffr4[i][4];
		bbi=coeffr4[i][5];
		cci=coeffr4[i][6];
		ddi=coeffr4[i][7];
		
		theta=(1-tau)+aai*pow((del-1)*(del-1),0.5/bti);
		delta=theta*theta+bbi*pow((del-1)*(del-1),ai);
		psi=exp(-cci*(del-1)*(del-1)-ddi*(tau-1)*(tau-1));
		dd=deltd(del,delta,aai,theta,bti,bbi,ai);
		ddd=deltdd(dd,del,delta,aai,theta,bti,bbi,ai);
//		printf("%lf,,,,,,%lf,,,,,\n",dd,ddd);
		res+=ni*(pow(delta,bi)*(2*psid(del,psi,cci)+del*psidd(del,psi,cci))+2*deltabid(bi,delta,dd)*(psi+del*psid(del,psi,cci))+deltabidd(bi,delta,dd,ddd)*del*psi);
	}
//	printf("..............%lf\n",res);
	return res;
}

double phirt(double del, double tau)
{
	int i;
	double ni,di,ti,ci,afi,bti,gmi,epi,ai,bi,aai,bbi,cci,ddi,delta,psi,theta;
	double dd,ddd;
	double res=0;
	for(i=0;i<7;i++)
	{
		ni=coeffr1[i][0];
		di=coeffr1[i][1];
		ti=coeffr1[i][2];
		res+=ni*ti*pow(tau,ti-1)*pow(del,di);
	}
	
	for(i=0;i<27;i++)
	{
		ni=coeffr2[i][0];
		di=coeffr2[i][1];
		ti=coeffr2[i][2];
		ci=coeffr2[i][3];
		res+=ni*ti*pow(tau,ti-1)*pow(del,di)*exp(-pow(del,ci));
	}
	
	for(i=0;i<5;i++)
	{
		ni=coeffr3[i][0];
		di=coeffr3[i][1];
		ti=coeffr3[i][2];
		afi=coeffr3[i][3];
		bti=coeffr3[i][4];
		gmi=coeffr3[i][5];
		epi=coeffr3[i][6];
		res+=ni*pow(del,di)*pow(tau,ti)*exp(-afi*(del-epi)*(del-epi)-bti*(tau-gmi)*(tau-gmi))*(ti/tau-2*bti*(tau-gmi));
	}
	
	for(i=0;i<3;i++)
	{
		ni=coeffr4[i][0];
		ai=coeffr4[i][1];
		bi=coeffr4[i][2];
		bti=coeffr4[i][3];
		aai=coeffr4[i][4];
		bbi=coeffr4[i][5];
		cci=coeffr4[i][6];
		ddi=coeffr4[i][7];
		
		theta=(1-tau)+aai*pow((del-1)*(del-1),0.5/bti);
		delta=theta*theta+bbi*pow((del-1)*(del-1),ai);
		psi=exp(-cci*(del-1)*(del-1)-ddi*(tau-1)*(tau-1));
		res+=ni*del*(deltabit(bi,theta,delta)*psi+pow(delta,bi)*psit(tau,psi,ddi));
	}
	return res;
}

double phirtt(double del, double tau)
{
	int i;
	double ni,di,ti,ci,afi,bti,gmi,epi,ai,bi,aai,bbi,cci,ddi,delta,psi,theta;
	double dd,ddd;
	double res=0;
	for(i=0;i<7;i++)
	{
		ni=coeffr1[i][0];
		di=coeffr1[i][1];
		ti=coeffr1[i][2];
		res+=ni*ti*(ti-1)*pow(tau,ti-2)*pow(del,di);
	}
	
	for(i=0;i<27;i++)
	{
		ni=coeffr2[i][0];
		di=coeffr2[i][1];
		ti=coeffr2[i][2];
		ci=coeffr2[i][3];
		res+=ni*ti*(ti-1)*pow(tau,ti-2)*pow(del,di)*exp(-pow(del,ci));
	}
	
	for(i=0;i<5;i++)
	{
		ni=coeffr3[i][0];
		di=coeffr3[i][1];
		ti=coeffr3[i][2];
		afi=coeffr3[i][3];
		bti=coeffr3[i][4];
		gmi=coeffr3[i][5];
		epi=coeffr3[i][6];
		res+=ni*pow(del,di)*pow(tau,ti)*exp(-afi*(del-epi)*(del-epi)-bti*(tau-gmi)*(tau-gmi))*(pow(ti/tau-2*bti*(tau-gmi),2)-ti/(tau*tau)-2*bti);
	}
	
	for(i=0;i<3;i++)
	{
		ni=coeffr4[i][0];
		ai=coeffr4[i][1];
		bi=coeffr4[i][2];
		bti=coeffr4[i][3];
		aai=coeffr4[i][4];
		bbi=coeffr4[i][5];
		cci=coeffr4[i][6];
		ddi=coeffr4[i][7];
		
		theta=(1-tau)+aai*pow((del-1)*(del-1),0.5/bti);
		delta=theta*theta+bbi*pow((del-1)*(del-1),ai);
		psi=exp(-cci*(del-1)*(del-1)-ddi*(tau-1)*(tau-1));
		res+=ni*del*(deltabitt(bi,theta,delta)*psi+2*deltabit(bi,theta,delta)*psit(tau,psi,ddi)+pow(delta,bi)*psitt(tau,psi,ddi));
	}
	return res;
}

double phirdt(double del, double tau)
{
	int i;
	double ni,di,ti,ci,afi,bti,gmi,epi,ai,bi,aai,bbi,cci,ddi,delta,psi,theta;
	double dd,ddd;
	double res=0;
	for(i=0;i<7;i++)
	{
		ni=coeffr1[i][0];
		di=coeffr1[i][1];
		ti=coeffr1[i][2];
		res+=ni*ti*di*pow(tau,ti-1)*pow(del,di-1);
	}
	
	for(i=0;i<27;i++)
	{
		ni=coeffr2[i][0];
		di=coeffr2[i][1];
		ti=coeffr2[i][2];
		ci=coeffr2[i][3];
		res+=ni*exp(-pow(del,ci))*pow(del,di-1)*ti*pow(tau,ti-1)*(di-ci*pow(del,ci));
	}
	
	for(i=0;i<5;i++)
	{
		ni=coeffr3[i][0];
		di=coeffr3[i][1];
		ti=coeffr3[i][2];
		afi=coeffr3[i][3];
		bti=coeffr3[i][4];
		gmi=coeffr3[i][5];
		epi=coeffr3[i][6];
		res+=ni*pow(del,di)*pow(tau,ti)*exp(-afi*(del-epi)*(del-epi)-bti*(tau-gmi)*(tau-gmi))*(di/del-2*afi*(del-epi))*(ti/tau-2*bti*(tau-gmi));
	}
	
	for(i=0;i<3;i++)
	{
		ni=coeffr4[i][0];
		ai=coeffr4[i][1];
		bi=coeffr4[i][2];
		bti=coeffr4[i][3];
		aai=coeffr4[i][4];
		bbi=coeffr4[i][5];
		cci=coeffr4[i][6];
		ddi=coeffr4[i][7];
		
		theta=(1-tau)+aai*pow((del-1)*(del-1),0.5/bti);
		delta=theta*theta+bbi*pow((del-1)*(del-1),ai);
		psi=exp(-cci*(del-1)*(del-1)-ddi*(tau-1)*(tau-1));
		dd=deltd(del,delta,aai,theta,bti,bbi,ai);
		ddd=deltdd(dd,del,delta,aai,theta,bti,bbi,ai);
		
		res+=ni*(pow(delta,bi)*(psit(tau,psi,ddi)+del*psidt(tau,del,psi,cci,ddi))+del*deltabid(bi,delta,dd)*psit(tau,psi,ddi)+deltabit(bi,theta,delta)*(psi+del*psid(del,psi,cci))+deltabidt(aai,bi,bti,delta,del,theta,dd)*del*psi);
	}
	return res;
}

double p_tr(double temp, double dens)
{
	double tau,del;
	tau=T_C/temp;
	del=dens/r_c;
	return dens*GAS_CONST_CO2*temp*(1+del*phird(del,tau));
}

double s_tr(double temp, double dens)
{
	double tau, del;
	tau=T_C/temp;
	del=dens/r_c;
	return GAS_CONST_CO2*(tau*(phi0t(tau)+phirt(del,tau))-phi0(del,tau)-phir(del,tau)); 
}

double u_tr(double temp, double dens)
{
	double tau,del;
	tau=T_C/temp;
	del=dens/r_c;
	return GAS_CONST_CO2*temp*tau*(phi0t(tau)+phirt(del,tau));
} 

double cv_tr(double temp, double dens)
{
	double tau,del;
	tau=T_C/temp;
	del=dens/r_c;
	return -GAS_CONST_CO2*tau*tau*(phi0tt(tau)+phirtt(del,tau));
}

double h_tr(double temp, double dens)
{
	double tau,del;
	tau=T_C/temp;
	del=dens/r_c;
	return GAS_CONST_CO2*temp*(1+tau*(phi0t(tau)+phirt(del,tau))+del*phird(del,tau));	
}

double cp_tr(double temp, double dens)
{
	double tau,del;
	tau=T_C/temp;
	del=dens/r_c;
	return GAS_CONST_CO2*(-tau*tau*(phi0tt(tau)+phirtt(del,tau))+pow(1+del*phird(del,tau)-del*tau*phirdt(del,tau),2)/(1+2*del*phird(del,tau)+del*del*phirdd(del,tau)));
}

double cs_tr(double temp, double dens)
{
	double tau,del;
	tau=T_C/temp;
	del=dens/r_c;
	return GAS_CONST_CO2*(-tau*tau*(phi0tt(tau)+phirtt(del,tau))+(1+del*phird(del,tau)-del*tau*phirdt(del,tau))*(1+del*phird(del,tau)-del*tau*phirdt(del,tau)-r_c*1.00000/(GAS_CONST_CO2*del))/(1+2*del*phird(del,tau)+del*del*phirdd(del,tau)));

}

double w_tr(double temp, double dens)
{
	double tau,del;
	tau=T_C/temp;
	del=dens/r_c;
	return sqrt(GAS_CONST_CO2*temp*(1+2*del*phird(del,tau)+del*del*phirdd(del,tau)-pow((1+del*phird(del,tau)-del*tau*phirdt(del,tau))/tau,2)/(phi0tt(tau)+phirtt(del,tau))));		
}

double miu_tr(double temp, double dens)
{
	double tau,del;
	tau=T_C/temp;
	del=dens/r_c;
	return -(del*phird(del,tau)+del*del*phirdd(del,tau)+del*tau*phirdt(del,tau))/(pow(1+del*phird(del,tau)-del*tau*phirdt(del,tau),2)-tau*tau*(phi0tt(tau)+phirtt(del,tau))*(1+2*del*phird(del,tau)+del*del*phirdd(del,tau)))/(GAS_CONST_CO2*dens);
}

double f_tr(double temp, double dens)
{
	double tau,del;
	tau=T_C/temp;
	del=dens/r_c;
	return exp(phir(del,tau)+del*phird(del,tau)-log(1+del*phird(del,tau)));
}