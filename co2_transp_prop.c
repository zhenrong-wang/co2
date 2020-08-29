#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "co2_property_calc.h"
#include "co2_property_calc.c"

const double tab1[5]={
	0.235156,-0.491266,5.211155e-2,5.347906e-2,-1.537102e-2,
};

const double dxx[5]={
	0.4071119e-2,0.7198037e-4,0.2411697e-16,0.2971072e-22,-0.1627888e-22,
};

const double lk[4]={
	1.51874307e-2,2.80674040e-2,2.28564190e-2,-7.41624210e-3
};

const double Bc[6][3]={
	1,1.00128e-2,4.30829e-3,
	2,5.60488e-2,-3.58563e-2,
	3,-8.11620e-2,6.71480e-2,
	4,6.24337e-2,-5.22855e-2,
	5,-2.06336e-2,1.74571e-2,
	6,2.53248e-3,-1.96414e-3,
};

double eta0(double temp)
{
	double sum=0;
	int i;
	for(i=0;i<5;i++)
	{
		sum+=tab1[i]*pow(log(temp/251.196),i);
	}
	return 1.00697*sqrt(temp)/exp(sum);
}

double deleta(double temp, double dens)
{
	return dxx[0]*dens+dxx[1]*dens*dens+dxx[2]*pow(dens,6)/pow(temp/251.196,3)+dxx[3]*pow(dens,8)+dxx[4]*pow(dens,8)/(temp/251.196);
}

double visc_calc(double temp, double dens)
{
	return eta0(temp)+deleta(temp,dens);
}

double lmd0(double temp)
{	
	double sum=0;
	int i;
	for(i=0;i<4;i++)
	{
		sum+=lk[i]/pow(temp/T_C,i);
	}
	return sqrt(temp/T_C)/sum;
}

double deltlmd(double temp, double dens)
{
	double sum=0;
	int i;
	for(i=0;i<6;i++)
	{
		sum+=(Bc[i][1]+Bc[i][2]*(temp/T_C))*pow(dens/r_c,Bc[i][0]);
	}
	return 1e3*sum;
}

double deltlmdc_emp(double temp, double dens)
{
	double dtc,drc;
	dtc=temp/T_C-1;
	drc=dens/r_c-1;
	return (-17.47-44.88*dtc)/(0.8563-0.4503*drc-7.197*dtc-exp(8.865*dtc+4.16*drc*drc+2.302*dtc*drc-drc*drc*drc));
}

double drdp(double temp, double dens)
{
	double del=dens/r_c;
	double tau=T_C/temp;
	return 1/(GAS_CONST_CO2*temp+2*dens*phird(del,tau)/r_c+dens*dens*phirdd(del,tau)/(r_c*r_c));
}

double deltlmdc(co2_prop* prop, double temp, double dens)
{
	double rd=1.02;
	double v=0.63;
	double gm=1.239;
	double gm_upper=0.052;
	double ep0=1.5e-10;
	double t_ref=456.19;
	double qd=0.25e10;
	double kb=1.380648813e-23;
	double ep,vis;
	double omg,omg0,pi;
	double cp,cv;
	cp=prop->spe_h_p;
	cv=prop->spe_h_v;
	ep=ep0*pow(P_C*dens/(gm_upper*r_c*r_c),v/gm)*pow(drdp(temp,dens)-drdp(t_ref,dens)*(t_ref/temp),v/gm);
	pi=3.14159265358;
	omg=((cp-cv)*atan(qd*ep)/cp+cv*qd*ep/cp)*2/pi;
	omg0=(1-exp(-1/(1/(qd*ep)+pow(qd*ep*r_c/dens,2)/3)))*2/pi;
	vis=visc_calc(temp,dens);
	printf("...........%lf,%lf,%lf,%.15lf,%lf,%lf\n",vis,cp,cv,ep,omg,omg0);
	return 1e9*dens*cp*rd*kb*temp*(omg-omg0)/(6*pi*vis*ep);
}

void calc_trans_prop(co2_prop* prop)
{
	double lm1,lm2,lm3;
	lm1=lmd0(prop->temp);
	lm2=deltlmd(prop->temp,prop->dens);
	lm3=deltlmdc_emp(prop->temp,prop->dens);
	prop->viscous=visc_calc(prop->temp,prop->dens);
//	printf("%lf,%lf,%lf,%lf\n\n",lm1,lm2,lm3,lm1+lm2+lm3);
	prop->thcond=lmd0(prop->temp)+deltlmd(prop->temp,prop->dens)+deltlmdc_emp(prop->temp,prop->dens);
}

int main()
{
	
	co2_prop pp,pp2;
//	co2_prop_calc_pt(&pp,&pp2,7.5e6,304);
//	calc_trans_prop(&pp);
//	printf("%lf,%lf\n",pp.viscous,pp.thcond);
	
//	co2_prop_calc_pt(&pp,&pp2,7e6,304);
//	calc_trans_prop(&pp);
//	printf("%lf,%lf\n",pp.viscous,pp.thcond);
	
//	co2_prop_calc_tr(&pp,310,400);
//	calc_trans_prop(&pp);
//	printf("%lf,%lf\n",pp.viscous,pp.thcond);
//	printf("....%lf,%lf,%lf,%lf\n",pp.pres,pp.temp,pp.spe_h_v,pp.spe_h_p);
//	printf("\n###%lf,%lf\n",lmd0(310)+deltlmd(310,400),deltlmdc(&pp,310,400)); 
	
//	co2_prop_calc_tr(&pp,250,1058);
//	printf("....%lf,%lf,%lf,%lf\n",pp.pres,pp.temp,pp.spe_h_v,pp.spe_h_p);
//	printf("\n###%lf,%lf\n",lmd0(250)+deltlmd(250,1058),deltlmdc(&pp,310,400)); 
	
	
//	co2_prop_calc_pt(&pp,&pp2,0.1e6,240);
//	printf("%lf\n",visc_calc(pp.temp,pp.dens));
/*	co2_prop pp,pp2;
	co2_prop_calc_pt(&pp,&pp2,0.1e6,220);
	printf("%lf,%lf\n",pp.temp,pp.dens);
	printf("%lf\n",visc_calc(pp.dens,pp.temp));
		co2_prop_calc_pt(&pp,&pp2,0.1e6,300);
	printf("%lf,%lf\n",pp.temp,pp.dens);
	printf("%lf\n",visc_calc(pp.dens,pp.temp));
		co2_prop_calc_pt(&pp,&pp2,0.1e6,800);
	printf("%lf,%lf\n",pp.temp,pp.dens);
	printf("%lf\n",visc_calc(pp.dens,pp.temp));
		co2_prop_calc_pt(&pp,&pp2,7e6,304);
	printf("%lf,%lf\n",pp.temp,pp.dens);
	printf("%lf\n",visc_calc(pp.dens,pp.temp));
		co2_prop_calc_pt(&pp,&pp2,15e6,220);
	printf("%lf,%lf\n",pp.temp,pp.dens);
	printf("%lf\n",visc_calc(pp.dens,pp.temp));
			co2_prop_calc_pt(&pp,&pp2,50e6,300);
	printf("%lf,%lf\n",pp.temp,pp.dens);
	printf("%lf\n",visc_calc(pp.dens,pp.temp));
				co2_prop_calc_pt(&pp,&pp2,75e6,800);
	printf("%lf,%lf\n",pp.temp,pp.dens);
	printf("%lf\n",visc_calc(pp.dens,pp.temp));*/
//	printf("%lf\n",visc_calc(304,254.3205));
} 