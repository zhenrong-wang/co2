/*
 *This code is distributed under the license: MIT License
 * Originally written by Zhenrong WANG
 * mailto: zhenrongwang@live.com
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "co2_property_calc.h"
#include "co2_property_calc_origeqns.h"
#define ERR 1e-6

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

int is_vld_prop(co2_prop *prop)
{
	if(isnan(prop->pres)||isnan(prop->temp)||isnan(prop->spe_vol)||isnan(prop->dens)||isnan(prop->spe_entr)||isnan(prop->spe_ener)||isnan(prop->spe_h_v)||isnan(prop->spe_enth)||isnan(prop->spe_h_p)||isnan(prop->speed_sound)||isnan(prop->Joule_Thompson))
		return -1;
	else if(prop->pres<0||prop->temp<0||prop->spe_vol<0||prop->dens<0||prop->spe_h_v<0||prop->spe_h_p<0||prop->speed_sound<0)
		return -1;
	else if(prop->speed_sound>2500||prop->spe_enth<-500e3)
		return -1;
	return 0;
}

void co2_prop_calc_tr(co2_prop *prop, double temp, double dens)
{	
	prop->pres=p_tr(temp,dens);
	prop->temp=temp;
	prop->spe_vol=1/dens;
	prop->dens=dens;
	prop->spe_entr=s_tr(temp,dens);
	prop->spe_ener=u_tr(temp,dens);
	prop->spe_h_v=cv_tr(temp,dens);
	prop->spe_enth=h_tr(temp,dens);
	prop->spe_h_p=cp_tr(temp,dens);
	prop->sat_liq_hc=cs_tr(temp,dens);
	prop->speed_sound=w_tr(temp,dens);
	prop->Joule_Thompson=miu_tr(temp,dens);
	prop->fugacity=f_tr(temp,dens);
	prop->vf=-1;
	prop->viscous=-1;
	prop->thcond=-1;
}

void swap_line(double left_matrix[3][3],double right_array[3],int i,int j)
{
    double temp;
    int k;
    temp=right_array[i];
    right_array[i]=right_array[j];
    right_array[j]=temp;
    for(k=0; k<3; k++)
    {
        temp=left_matrix[i][k];
        left_matrix[i][k]=left_matrix[j][k];
        left_matrix[j][k]=temp;
    }
}

int Gaussian_elimination(double *solution, double left_matrix[3][3],double right_hand_array[3])
{
    int i,j,k,kk;
    int max_index;
    double max_temp;
    double coeff_temp;
    double sigma_temp;
    int la,lb;

    for(i=0; i<2; i++)
    {
        max_temp=fabs(left_matrix[i][i]);
        max_index=i;
        for(j=i+1; j<3; j++)
        {
            if(fabs(left_matrix[j][i])>max_temp)
            {
                max_temp=fabs(left_matrix[j][i]);
                max_index=j;
            }
        }

        if(fabs(max_temp-0)<ERR)
        {
            return -1;
        }
        if(max_index!=i)
        {
            swap_line(left_matrix,right_hand_array,i,max_index);
        }

        for(j=i+1; j<3; j++)
        {
            coeff_temp=left_matrix[j][i]/left_matrix[i][i];
            right_hand_array[j]=right_hand_array[j]-coeff_temp*right_hand_array[i];
            for(k=i; k<3; k++)
            {
                left_matrix[j][k]=left_matrix[j][k]-coeff_temp*left_matrix[i][k];
            }
        }
    }
    
    if(fabs(left_matrix[2][2])<ERR)
    {
    	return -1;
    }
    solution[2]=right_hand_array[2]/left_matrix[2][2];
    solution[1]=(right_hand_array[1]-left_matrix[1][2]*solution[2])/left_matrix[1][1];
    solution[0]=(right_hand_array[0]-left_matrix[0][2]*solution[2]-left_matrix[0][1]*solution[1])/left_matrix[0][0];
//    printf("........%lf,%lf,%lf\n",solution[0],solution[1],solution[2]);
    return 0;
}

double find_max(double* solu)
{
	double max=fabs(solu[0]);
	double max2;
	if(fabs(solu[1])>fabs(solu[2]))
		max2=fabs(solu[1]);
	else
		max2=fabs(solu[2]);
	if(max>max2)
		return max;
	else
		return max2;
}

void calcfgh_dfgh_psat(double dfgh[3][3], double fgh[3], double tau, double d1, double d2, double ps)
{
//	int ii;
	fgh[0]=-(d1+d1*d1*phird(d1,tau)-ps*tau/(r_c*GAS_CONST_CO2*T_C));
	dfgh[0][0]=d1*d1*phirdt(d1,tau)-ps/(r_c*GAS_CONST_CO2*T_C);
	dfgh[0][1]=1+2*d1*phird(d1,tau)+d1*d1*phirdd(d1,tau);
	dfgh[0][2]=0;
	
	fgh[1]=-(d2+d2*d2*phird(d2,tau)-ps*tau/(r_c*GAS_CONST_CO2*T_C));
	dfgh[1][0]=d2*d2*phirdt(d2,tau)-ps/(r_c*GAS_CONST_CO2*T_C);
	dfgh[1][1]=0;
	dfgh[1][2]=1+2*d2*phird(d2,tau)+d2*d2*phirdd(d2,tau);
	
	fgh[2]=-(phir(d1,tau)-phir(d2,tau)-(ps*tau/(r_c*GAS_CONST_CO2*T_C))*(1/d2-1/d1)+log(d1/d2));
	dfgh[2][0]=phirt(d1,tau)-phirt(d2,tau)-ps*(1/d2-1/d1)/(r_c*GAS_CONST_CO2*T_C);
	dfgh[2][1]=phird(d1,tau)-ps*tau/(r_c*GAS_CONST_CO2*T_C*d1*d1)+1/d1;
	dfgh[2][2]=-phird(d2,tau)+ps*tau/(r_c*GAS_CONST_CO2*T_C*d2*d2)-1/d2;
//	printf("%lf,%lf,%lf\n",fgh[0],fgh[1],fgh[2]);
//					for(ii=0;ii<3;ii++)
//				{
//					printf("%lf,%lf,%lf\n",dfgh[ii][0],dfgh[ii][1],dfgh[ii][2]);
//				}
}

int calc_tsat_p(double *tsat, double *densl, double *densg, double psat) 
{
	double solu_iter[3]={0,0,0};
	double left_matrix[3][3]={
		0,0,0,
		0,0,0,
		0,0,0,
	};
	int flag=0;
	int iter=0;
	int solu_flag=-1;
	int ii,jj,kk;
	double right_array[3]={0,0,0};
	double tau_ini[10]={T_C/217,T_C/227,T_C/237,T_C/247,T_C/257,T_C/267,T_C/277,T_C/287,T_C/297,T_C/304};
	double d1_ini[3]={1000/r_c,800/r_c,600/r_c};
	double d2_ini[9]={20/r_c,60/r_c,100/r_c,150/r_c,200/r_c,250/r_c,300/r_c,350/r_c,400/r_c};
	double tau_nxt,tau_prev,d1_prev,d1_nxt,d2_prev,d2_nxt;
	int i,j,k;
	for(i=0;i<10;i++)
	{
		for(j=0;j<3;j++)
		{
			for(k=0;k<9;k++)
			{
//				printf("........%lf\t%lf\t%lf\n\n",tau_ini[i],d1_ini[j],d2_ini[k]);
				calcfgh_dfgh_psat(left_matrix,right_array,tau_ini[i],d1_ini[j],d2_ini[k],psat);
//				printf("%lf,%lf,%lf\n",right_array[0],right_array[1],right_array[2]);
//				for(ii=0;ii<3;ii++)
//				{
//					printf("%lf,%lf,%lf\n",left_matrix[ii][0],left_matrix[ii][1],left_matrix[ii][2]);
//				}
				flag=Gaussian_elimination(solu_iter,left_matrix,right_array);
//				printf("..........%lf,%lf,%lf\n",solu_iter[0],solu_iter[1],solu_iter[2]);
				if(flag==-1)
					continue;
				tau_nxt=tau_ini[i]+solu_iter[0];
				d1_nxt=d1_ini[j]+solu_iter[1];
				d2_nxt=d2_ini[k]+solu_iter[2];
//				printf("%lf,%lf,%lf\n",tau_nxt,d1_nxt,d2_nxt);
//				printf("%lf,%lf,%lf\n",tau_ini[i],d1_ini[j],d2_ini[k]);
				iter=0;
				do
				{
					tau_prev=tau_nxt;
					d1_prev=d1_nxt;
					d2_prev=d2_nxt;
					calcfgh_dfgh_psat(left_matrix,right_array,tau_prev,d1_prev,d2_prev,psat);
//									printf("%lf,%lf,%lf\n",right_array[0],right_array[1],right_array[2]);
					flag=Gaussian_elimination(solu_iter,left_matrix,right_array);
					if(flag==-1)
						break;
					tau_nxt=tau_prev+solu_iter[0];
					d1_nxt=d1_prev+solu_iter[1];
					d2_nxt=d2_prev+solu_iter[2];
					iter++;
				}while(iter<100&&find_max(solu_iter)>ERR);
				if(iter==100||flag==-1)
				{
//					printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%d,%d\n",iter,flag);
					continue;
				}
				else if(!isnan(tau_nxt)&&!isnan(d1_nxt)&&!isnan(d2_nxt))
				{
		 			*tsat=T_C/tau_nxt;
			 		*densl=d1_nxt*r_c;
			 		*densg=d2_nxt*r_c;
			 		return 0;
//			 		printf("{{{{{{{{{{%d,%lf,%lf,%lf\n\n",iter,tau_nxt,d1_nxt,d2_nxt);	
				}
				else
					continue;				
			}
		}
	}
	return -1; 
//	printf("######%d,%d,%d\n",i,j,k);
//	printf("######%lf,%lf,%lf\n",T_C/solu[0],solu[1]*r_c,solu[2]*r_c);
}


void calcfgh_dfgh_tsat(double dfgh[3][3], double fgh[3], double ps, double d1, double d2, double ts)
{
	int ii;
	double tau=T_C/ts;
	double coeff=r_c*GAS_CONST_CO2*ts;
	fgh[0]=-coeff*(d1+d1*d1*phird(d1,tau))+ps;
	dfgh[0][0]=-1;
	dfgh[0][1]=coeff*(1+2*d1*phird(d1,tau)+d1*d1*phirdd(d1,tau));
	dfgh[0][2]=0;
	
	fgh[1]=-coeff*(d2+d2*d2*phird(d2,tau))+ps;
	dfgh[1][0]=-1;
	dfgh[1][1]=0;
	dfgh[1][2]=coeff*(1+2*d2*phird(d2,tau)+d2*d2*phirdd(d2,tau));
	
	fgh[2]=-coeff*(phir(d1,tau)-phir(d2,tau))+ps*(1/d2-1/d1)-coeff*log(d1/d2);
	dfgh[2][0]=1/d1-1/d2;
	dfgh[2][1]=coeff*(phird(d1,tau)+1/d1)-ps/(d1*d1);
	dfgh[2][2]=coeff*(-phird(d2,tau)-1/d2)+ps/(d2*d2);
//	printf("%lf,%lf,%lf\n",fgh[0],fgh[1],fgh[2]);
//					for(ii=0;ii<3;ii++)
//				{
//					printf("%.12lf,%lf,%lf\n",dfgh[ii][0],dfgh[ii][1],dfgh[ii][2]);
//				}
}

int calc_psat_t(double *psat, double *densl, double *densg, double tsat) 
{
	double solu_iter[3]={0,0,0};
	double left_matrix[3][3]={
		0,0,0,
		0,0,0,
		0,0,0,
	};
	int flag=0;
	int iter=0;
	int solu_flag=-1;
	int ii,jj,kk;
	double right_array[3]={0,0,0};
	double ps_ini[10]={0.6e6,0.8e6,1.2e6,2.2e6,3.2e6,4.2e6,5.2e6,6.2e6,7.2e6,7.3e6};
	double d1_ini[3]={1000/r_c,800/r_c,600/r_c};
	double d2_ini[9]={20/r_c,60/r_c,100/r_c,150/r_c,200/r_c,250/r_c,300/r_c,350/r_c,400/r_c};
	double ps_nxt,ps_prev,d1_prev,d1_nxt,d2_prev,d2_nxt;
	int i,j,k;
	for(i=0;i<1;i++)
	{
		for(j=0;j<1;j++)
		{
			for(k=0;k<1;k++)
			{
				calcfgh_dfgh_tsat(left_matrix,right_array,ps_ini[i],d1_ini[j],d2_ini[k],tsat);
				flag=Gaussian_elimination(solu_iter,left_matrix,right_array);
//								printf("%lf,%lf,%lf\n",right_array[0],right_array[1],right_array[2]);
//				for(ii=0;ii<3;ii++)
//				{
//					printf("%lf,%lf,%lf\n",left_matrix[ii][0],left_matrix[ii][1],left_matrix[ii][2]);
//				}
				if(flag==-1)
					continue;
				ps_nxt=ps_ini[i]+solu_iter[0];
				d1_nxt=d1_ini[j]+solu_iter[1];
				d2_nxt=d2_ini[k]+solu_iter[2];
				iter=0;
				do
				{
					ps_prev=ps_nxt;
					d1_prev=d1_nxt;
					d2_prev=d2_nxt;
					calcfgh_dfgh_tsat(left_matrix,right_array,ps_prev,d1_prev,d2_prev,tsat);
					flag=Gaussian_elimination(solu_iter,left_matrix,right_array);
					if(flag==-1)
						break;
					ps_nxt=ps_prev+solu_iter[0];
					d1_nxt=d1_prev+solu_iter[1];
					d2_nxt=d2_prev+solu_iter[2];
					iter++;
				}while(iter<100&&find_max(solu_iter)>ERR);
				if(iter==100||flag==-1)
					continue;
				else if(!isnan(ps_nxt)&&!isnan(d1_nxt)&&!isnan(d2_nxt))
				{
		 			*psat=ps_nxt;
			 		*densl=d1_nxt*r_c;
			 		*densg=d2_nxt*r_c;
			 		return 0;	
				}
				else
					continue;				
			}
		}
	}
	return -1;
}

void calcfd_dfd_pt(double *f, double *df, double del, double pres, double tau)
{
	*f=del+del*del*phird(del,tau)-pres*tau/(r_c*GAS_CONST_CO2*T_C);
	*df=1+2*del*phird(del,tau)+del*del*phirdd(del,tau);
}

int dens_pt(double *dens, double pres, double temp)
{
	int i,iter;
	double tau=T_C/temp;
	double f,df; 
	co2_prop pp;
	double del_ini[12]={
		0.1/r_c,1300/r_c,10/r_c,1200/r_c,200/r_c,1000/r_c,300/r_c,1400/r_c,400/r_c,800/r_c,467/r_c,600/r_c,
	};
	double del_ini2;
	double del_prev,del_nxt;
	for(i=0;i<12;i++)
	{
//		if(pres>0.7e6)
//			del_ini2=del_ini[i+1];/
//		else
			del_ini2=del_ini[i];
//			printf(".......%lf,%lf\n",pres,del_ini2);
		iter=0;
		calcfd_dfd_pt(&f,&df,del_ini2,pres,tau);
		del_nxt=del_ini2-f/df;
		do
		{
			del_prev=del_nxt;
			calcfd_dfd_pt(&f,&df,del_prev,pres,tau);
			del_nxt=del_prev-f/df;
			iter++;
		}while(iter<100&&fabs(del_nxt-del_prev)>ERR);
		if(iter==100)
			continue;
		else
		{
			if(!isnan(del_nxt))
			{
				co2_prop_calc_tr(&pp,temp,del_nxt*r_c);
				if(is_vld_prop(&pp)==-1||fabs(pres-pp.pres)/pres>1e-3)
					continue;
				else
				{
					*dens=del_nxt*r_c;
					return 0;
				}
			}
			else
				continue;
		}
	}
	return -1;
}

void calcfd_dfd_pr(double *f, double *df, double tau, double pres, double del)
{
	*f=del+del*del*phird(del,tau)-pres*tau/(r_c*GAS_CONST_CO2*T_C);
	*df=del*del*phirdt(del,tau)-pres/(r_c*GAS_CONST_CO2*T_C);
}

int temp_pr(double *temp, double pres, double dens)
{
	int i,iter;
	double del=dens/r_c;
	double f,df; 
	co2_prop pp;
	double tau_ini[10]={T_C/217,T_C/304,T_C/317,T_C/417,T_C/517,T_C/617,T_C/717,T_C/817,T_C/917,T_C/1017};
	double tau_prev,tau_nxt;
	for(i=0;i<10;i++)
	{
//		if(pres>0.7e6)
//			del_ini2=del_ini[i+1];/
//		else
//			del_ini2=del_ini[i];
//			printf(".......%lf,%lf\n",pres,del_ini2);
		iter=0;
		calcfd_dfd_pr(&f,&df,tau_ini[i],pres,del);
		tau_nxt=tau_ini[i]-f/df;
		do
		{
			tau_prev=tau_nxt;
			calcfd_dfd_pr(&f,&df,tau_prev,pres,del);
			tau_nxt=tau_prev-f/df;
			iter++;
//			printf("\n%lf.......\n",f/df);
		}while(iter<100&&fabs(tau_nxt-tau_prev)>ERR);
		if(iter==100)
			continue;
		else
		{
			if(!isnan(tau_nxt))
			{
				co2_prop_calc_tr(&pp,T_C/tau_nxt,dens);
				if(is_vld_prop(&pp)==0&&fabs(pres-pp.pres)/pres<1e-4)
				{
//					printf("\n.........\n");
					*temp=T_C/tau_nxt;
					return 0;
				}
				else
					continue;
			}
			else
				continue;
		}
	}
//	printf("\n%d....",i);
	return -1;
}

void calcfd_dfd_th(double *f, double *df, double del, double tau, double h)
{
	*f=1+tau*(phi0t(tau)+phirt(del,tau))+del*phird(del,tau)-h*tau/(GAS_CONST_CO2*T_C);
	*df=tau*phirdt(del,tau)+del*phirdd(del,tau)+phird(del,tau);
}
int dens_th(double* dens, double temp, double h)
{
	int i,iter;
	double tau=T_C/temp;
	double f,df; 
	co2_prop pp;
	double del_ini[12]={
		0.1/r_c,1300/r_c,10/r_c,1200/r_c,200/r_c,1000/r_c,300/r_c,1400/r_c,400/r_c,800/r_c,467/r_c,600/r_c,
	};
	double del_prev,del_nxt;
	for(i=0;i<12;i++)
	{
//		if(pres>0.7e6)
//			del_ini2=del_ini[i+1];/
//		else
//			del_ini2=del_ini[i];
//			printf(".......%lf,%lf\n",pres,del_ini2);
		iter=0;
		calcfd_dfd_th(&f,&df,del_ini[i],tau,h);
		del_nxt=del_ini[i]-f/df;
		do
		{
			del_prev=del_nxt;
			calcfd_dfd_th(&f,&df,del_prev,tau,h);
			del_nxt=del_prev-f/df;
			iter++;
//			printf("\n%lf.......\n",f/df);
		}while(iter<100&&fabs(del_nxt-del_prev)>ERR);
		if(iter==100)
			continue;
		else
		{
			if(!isnan(del_nxt))
			{
				co2_prop_calc_tr(&pp,temp,del_nxt*r_c);
				if(is_vld_prop(&pp)==0&&fabs((h-pp.spe_enth)/h)<1e-4)
				{
//					printf("\n.........\n");
					*dens=del_nxt*r_c;
					return 0;
				}
				else
					continue;
			}
			else
				continue;
		}
	}
//	printf("\n%d....",i);
	return -1;
}

void calcfd_dfd_ts(double *f, double *df, double del, double tau, double s)
{
	*f=tau*(phi0t(tau)+phirt(del,tau))-phi0(del,tau)-phir(del,tau)-s/GAS_CONST_CO2;
	*df=tau*phirdt(del,tau)-phi0d(del)-phird(del,tau);
//	printf("\n%lf,%lf\n",*f,*df);
}
int dens_ts(double* dens, double temp, double s)
{
	int i,iter;
	double tau=T_C/temp;
	double f,df; 
	co2_prop pp;
	double del_ini[12]={
		0.1/r_c,1300/r_c,10/r_c,1200/r_c,200/r_c,1000/r_c,300/r_c,1400/r_c,400/r_c,800/r_c,467/r_c,600/r_c,
	};
	double del_prev,del_nxt;
	for(i=0;i<12;i++)
	{
//		if(pres>0.7e6)
//			del_ini2=del_ini[i+1];/
//		else
//			del_ini2=del_ini[i];
//			printf(".......%lf,%lf\n",pres,del_ini2);
		iter=0;
		calcfd_dfd_ts(&f,&df,del_ini[i],tau,s);
		del_nxt=del_ini[i]-f/df;
		do
		{
			del_prev=del_nxt;
			calcfd_dfd_ts(&f,&df,del_prev,tau,s);
			del_nxt=del_prev-f/df;
			iter++;
//			printf("\n%lf.......\n",f/df);
		}while(iter<100&&fabs(del_nxt-del_prev)>ERR);
		if(iter==100)
			continue;
		else
		{
			if(!isnan(del_nxt))
			{
				co2_prop_calc_tr(&pp,temp,del_nxt*r_c);
				if(is_vld_prop(&pp)==0&&fabs((s-pp.spe_entr)/s)<1e-4)
				{
//					printf("\n.........\n");
					*dens=del_nxt*r_c;
					return 0;
				}
				else
					continue;
			}
			else
				continue;
		}
	}
//	printf("\n%d....",i);
	return -1;
}


void calcfg_dfg_hs(double *f, double *g, double *dfdt, double *dfdd, double *dgdt, double *dgdd, double tau, double del, double h, double s)
{
	*f=1+tau*(phi0t(tau)+phirt(del,tau))+del*phird(del,tau)-h*tau/(GAS_CONST_CO2*T_C);
	*dfdt=tau*(phi0tt(tau)+phirtt(del,tau))+phi0t(tau)+phirt(del,tau)+del*phirdt(del,tau)-h/(GAS_CONST_CO2*T_C);
	*dfdd=tau*phirdt(del,tau)+del*phirdd(del,tau)+phird(del,tau);
	*g=tau*(phi0t(tau)+phirt(del,tau))-phi0(del,tau)-phir(del,tau)-s/GAS_CONST_CO2;
	*dgdt=tau*(phi0tt(tau)+phirtt(del,tau));
	*dgdd=tau*phirdt(del,tau)-phi0d(del)-phird(del,tau);
}

int td_hs(double *temp, double *dens, double h, double s)
{
	double tau_prev,tau_nxt;
	double del_prev,del_nxt,t2,d2;
	int i=0;
	double f,g,dfdt,dgdt,dfdd,dgdd;
	int t_index,d_index;
	int flag=0;
	double dt,dd,maxd;
	co2_prop prop;
	double tau_ini[10]={T_C/217,T_C/304,T_C/317,T_C/417,T_C/517,T_C/617,T_C/717,T_C/817,T_C/917,T_C/1017};
	double del_ini[12]={
		1/r_c,10/r_c,100/r_c,200/r_c,300/r_c,400/r_c,467/r_c,600/r_c,800/r_c,1000/r_c,1200/r_c,1400/r_c
	};
	for(t_index=0;t_index<10;t_index++)
	{
		for(d_index=0;d_index<12;d_index++)
		{
//			printf("\t%lf,%lf\n",pi_ini[p_index],tau_ini[t_index]);
			calcfg_dfg_hs(&f,&g,&dfdt,&dfdd,&dgdt,&dgdd,tau_ini[t_index],del_ini[d_index],h,s);
//			printf("%\n!!!!!!%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\n\n",f,g,dfdt,dfdd,dgdt,dgdd);
			if(fabs(dfdd*dgdt-dfdt*dgdd)<ERR)
			{
				continue;
			}
			dd=(g*dfdt-f*dgdt)/(dfdd*dgdt-dgdd*dfdt);
			dt=(g*dfdd-f*dgdd)/(dfdt*dgdd-dgdt*dfdd);
			tau_nxt=tau_ini[t_index]+dt;
			del_nxt=del_ini[d_index]+dd;
			flag=0;
			i=0;
//	printf("\n%.12lf\t%.12lf\n",dt,dd);
			do
			{
				tau_prev=tau_nxt;
				del_prev=del_nxt;
				calcfg_dfg_hs(&f,&g,&dfdt,&dfdd,&dgdt,&dgdd,tau_prev,del_prev,h,s);
				if(fabs(dfdd*dgdt-dfdt*dgdd)<ERR)
				{
					flag=-1;
					break;
				}
				dd=(g*dfdt-f*dgdt)/(dfdd*dgdt-dgdd*dfdt);
				dt=(g*dfdd-f*dgdd)/(dfdt*dgdd-dgdt*dfdd);
				tau_nxt=tau_prev+dt;
				del_nxt=del_prev+dd;
				if(fabs(dt)>fabs(dd))
					maxd=fabs(dt);
				else
					maxd=fabs(dd);
				i++;
//				printf("........%lf,%lf\n",dt,dd);
			}while(maxd>ERR&&i<100);
			if(i==100||flag==-1)
			{
				continue;
			}	
			else
			{
				if(!isnan(tau_nxt)&&!isnan(del_nxt))
				{
					t2=T_C/tau_nxt;
					d2=r_c*del_nxt;
					co2_prop_calc_tr(&prop,t2,d2);
					if(is_vld_prop(&prop)==0&&fabs((h-prop.spe_enth)/h)<1e-4||fabs((s-prop.spe_entr)/h)<1e-4)
					{
						*temp=t2;
						*dens=d2;
						return 0;
					}
					else
						continue;
				}
			}		
		}
	}
//	printf("::::::::%d,%d",t_index,d_index);
	return -1;
}

int hs_find_tsat(double *t, double *x, double spe_enth, double spe_entr)
{
	double tt1,tt2,ps1,ps2,xs1,xs2,xh1,xh2,ttm,psm,xhm,xsm;
	double dl1,dg1,dl2,dg2,dml,dmg;
	co2_prop pl,pg,pml,pmg;
	int flag1,flag2,flag3;
	int i=0;
	
	tt1=T_t;
	tt2=T_C;
	if(calc_psat_t(&ps1,&dl1,&dg1,tt1)==-1)
	{
		*t=-1;
		*x=-1;
		return -1;
	}
	co2_prop_calc_tr(&pl,tt1,dl1);
	co2_prop_calc_tr(&pg,tt1,dg1);
	xs1=(spe_entr-pl.spe_entr)/(pg.spe_entr-pl.spe_entr);
	xh1=(spe_enth-pl.spe_enth)/(pg.spe_enth-pg.spe_enth);
	
	if(calc_psat_t(&ps2,&dl2,&dg2,tt2)==-1)
	{
		*t=-1;
		*x=-1;
		return -1;
	}
	co2_prop_calc_tr(&pl,tt2,dl2);
	co2_prop_calc_tr(&pg,tt2,dg2);
	xs2=(spe_entr-pl.spe_entr)/(pg.spe_entr-pl.spe_entr);
	xh2=(spe_enth-pl.spe_enth)/(pg.spe_enth-pg.spe_enth);
//	printf("\n%lf,%lf\n",xh2,xs2);
	if((xh1-xs1)*(xh2-xs2)>0)
	{
		*t=-1;
		*x=-1;
		return -1;
	}
	do
	{
		ttm=0.5*(tt1+tt2);
//		printf("\n............%lf\n",ttm);
		if(calc_psat_t(&psm,&dml,&dmg,ttm)==-1)
		{
			*t=-1;
			*x=-1;
			return -1;
		}
		co2_prop_calc_tr(&pl,ttm,dml);
		co2_prop_calc_tr(&pg,ttm,dmg);
		xsm=(spe_entr-pl.spe_entr)/(pg.spe_entr-pl.spe_entr);
		xhm=(spe_enth-pl.spe_enth)/(pg.spe_enth-pg.spe_enth);
//			printf("\n...........%lf,%lf\n",xhm,xsm);
		i++;
		if(fabs(xhm-xsm)<1e-6)
		{
			*t=ttm;
			*x=0.5*(xhm+xsm);
			return 0;
		}
		else if((xhm-xsm)*(xh1-xs1)>0)
		{
			tt1=ttm;
		}
		else if((xhm-xsm)*(xh2-xs2)>0)
		{
			tt2=ttm;
		}
	}while(i<200);
	if(i==200)
	{
		*t=-1;
		*x=-1;
		return -1;
	}
}


void calcfg_dfg_ph(double *f, double *g, double *dfdt, double *dfdd, double *dgdt, double *dgdd, double tau, double del, double p, double h)
{
	*f=del+del*del*phird(del,tau)-p*tau/(r_c*GAS_CONST_CO2*T_C); 
	*dfdt=del*del*phirdt(del,tau)-p/(r_c*GAS_CONST_CO2*T_C);
	*dfdd=1+del*del*phirdd(del,tau)+2*del*phird(del,tau);
	
	*g=1+tau*(phi0t(tau)+phirt(del,tau))+del*phird(del,tau)-h*tau/(GAS_CONST_CO2*T_C);
	*dgdt=tau*(phi0tt(tau)+phirtt(del,tau))+phi0t(tau)+phirt(del,tau)+del*phirdt(del,tau)-h/(GAS_CONST_CO2*T_C);
	*dgdd=tau*phirdt(del,tau)+del*phirdd(del,tau)+phird(del,tau);
}

int td_ph(double *temp, double *dens, double p, double h)
{
	double tau_prev,tau_nxt;
	double del_prev,del_nxt,t2,d2;
	int i=0;
	double f,g,dfdt,dgdt,dfdd,dgdd;
	int t_index,d_index;
	int flag=0;
	double dt,dd,maxd,di;
	co2_prop prop;
	double tau_ini[10]={T_C/217,T_C/304,T_C/317,T_C/417,T_C/517,T_C/617,T_C/717,T_C/817,T_C/917,T_C/1017};
	double del_ini[12]={
		0.1/r_c,1300/r_c,10/r_c,1200/r_c,200/r_c,1000/r_c,300/r_c,1400/r_c,400/r_c,800/r_c,467/r_c,600/r_c,
	};
	for(t_index=0;t_index<10;t_index++)
	{
		for(d_index=0;d_index<12;d_index++)
		{
//			if(p>0.7e6)
//				di=del_ini[i+1];
//			else
				di=del_ini[i];
//			printf("\t%lf,%lf\n",pi_ini[p_index],tau_ini[t_index]);
			calcfg_dfg_ph(&f,&g,&dfdt,&dfdd,&dgdt,&dgdd,tau_ini[t_index],di,p,h);
//			printf("%\n!!!!!!%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\n\n",f,g,dfdt,dfdd,dgdt,dgdd);
			if(fabs(dfdd*dgdt-dfdt*dgdd)<ERR)
			{
				continue;
			}
			dd=(g*dfdt-f*dgdt)/(dfdd*dgdt-dgdd*dfdt);
			dt=(g*dfdd-f*dgdd)/(dfdt*dgdd-dgdt*dfdd);
			tau_nxt=tau_ini[t_index]+dt;
			del_nxt=del_ini[d_index]+dd;
			flag=0;
			i=0;
//	printf("\n%.12lf\t%.12lf\n",dt,dd);
			do
			{
				tau_prev=tau_nxt;
				del_prev=del_nxt;
				calcfg_dfg_ph(&f,&g,&dfdt,&dfdd,&dgdt,&dgdd,tau_prev,del_prev,p,h);
				if(fabs(dfdd*dgdt-dfdt*dgdd)<ERR)
				{
					flag=-1;
					break;
				}
				dd=(g*dfdt-f*dgdt)/(dfdd*dgdt-dgdd*dfdt);
				dt=(g*dfdd-f*dgdd)/(dfdt*dgdd-dgdt*dfdd);
				tau_nxt=tau_prev+dt;
				del_nxt=del_prev+dd;
				if(fabs(dt)>fabs(dd))
					maxd=fabs(dt);
				else
					maxd=fabs(dd);
				i++;
//				printf("........%lf,%lf\n",dt,dd);
			}while(maxd>ERR&&i<100);
			if(i==100||flag==-1)
			{
				continue;
			}	
			else
			{
				if(!isnan(tau_nxt)&&!isnan(del_nxt))
				{
					t2=T_C/tau_nxt;
					d2=r_c*del_nxt;
					co2_prop_calc_tr(&prop,t2,d2);
//					printf("\n~~~~~~~~~~~~%d\n",is_vld_prop(&prop));
					if(is_vld_prop(&prop)==0&&fabs(p-prop.pres)/p<1e-4&&fabs((h-prop.spe_enth)/h)<1e-4)
					{
						*temp=t2;
						*dens=d2;
						return 0;
					}
					else
						continue;
				}
			}		
		}
	}
//	printf("::::::::%d,%d",t_index,d_index);
	return -1;
}

void calcfg_dfg_ps(double *f, double *g, double *dfdt, double *dfdd, double *dgdt, double *dgdd, double tau, double del, double p, double s)
{
	*f=del+del*del*phird(del,tau)-p*tau/(r_c*GAS_CONST_CO2*T_C); 
	*dfdt=del*del*phirdt(del,tau)-p/(r_c*GAS_CONST_CO2*T_C);
	*dfdd=1+del*del*phirdd(del,tau)+2*del*phird(del,tau);
	
	*g=tau*(phi0t(tau)+phirt(del,tau))-phi0(del,tau)-phir(del,tau)-s/GAS_CONST_CO2;
	*dgdt=tau*(phi0tt(tau)+phirtt(del,tau));
	*dgdd=tau*phirdt(del,tau)-phi0d(del)-phird(del,tau);
}

int td_ps(double *temp, double *dens, double p, double s)
{
	double tau_prev,tau_nxt;
	double del_prev,del_nxt,t2,d2;
	int i=0;
	double f,g,dfdt,dgdt,dfdd,dgdd;
	int t_index,d_index;
	int flag=0;
	double dt,dd,maxd,di;
	co2_prop prop;
	double tau_ini[10]={T_C/217,T_C/304,T_C/317,T_C/417,T_C/517,T_C/617,T_C/717,T_C/817,T_C/917,T_C/1017};
	double del_ini[12]={
		0.1/r_c,1300/r_c,10/r_c,1200/r_c,200/r_c,1000/r_c,300/r_c,1400/r_c,400/r_c,800/r_c,467/r_c,600/r_c,
	};
	for(t_index=0;t_index<10;t_index++)
	{
		for(d_index=0;d_index<12;d_index++)
		{
//			if(p>0.7e6)
//				di=del_ini[i+1];
//			else
				di=del_ini[i];
//			printf("\t%lf,%lf\n",pi_ini[p_index],tau_ini[t_index]);
			calcfg_dfg_ps(&f,&g,&dfdt,&dfdd,&dgdt,&dgdd,tau_ini[t_index],di,p,s);
//			printf("%\n!!!!!!%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\n\n",f,g,dfdt,dfdd,dgdt,dgdd);
			if(fabs(dfdd*dgdt-dfdt*dgdd)<ERR)
			{
				continue;
			}
			dd=(g*dfdt-f*dgdt)/(dfdd*dgdt-dgdd*dfdt);
			dt=(g*dfdd-f*dgdd)/(dfdt*dgdd-dgdt*dfdd);
			tau_nxt=tau_ini[t_index]+dt;
			del_nxt=del_ini[d_index]+dd;
			flag=0;
			i=0;
//	printf("\n%.12lf\t%.12lf\n",dt,dd);
			do
			{
				tau_prev=tau_nxt;
				del_prev=del_nxt;
				calcfg_dfg_ps(&f,&g,&dfdt,&dfdd,&dgdt,&dgdd,tau_prev,del_prev,p,s);
				if(fabs(dfdd*dgdt-dfdt*dgdd)<ERR)
				{
					flag=-1;
					break;
				}
				dd=(g*dfdt-f*dgdt)/(dfdd*dgdt-dgdd*dfdt);
				dt=(g*dfdd-f*dgdd)/(dfdt*dgdd-dgdt*dfdd);
				tau_nxt=tau_prev+dt;
				del_nxt=del_prev+dd;
				if(fabs(dt)>fabs(dd))
					maxd=fabs(dt);
				else
					maxd=fabs(dd);
				i++;
//				printf("........%lf,%lf\n",dt,dd);
			}while(maxd>ERR&&i<100);
			if(i==100||flag==-1)
			{
				continue;
			}	
			else
			{
				if(!isnan(tau_nxt)&&!isnan(del_nxt))
				{
					t2=T_C/tau_nxt;
					d2=r_c*del_nxt;
					co2_prop_calc_tr(&prop,t2,d2);
					if(is_vld_prop(&prop)==0&&fabs(p-prop.pres)/p<1e-4&&fabs((s-prop.spe_entr)/s)<1e-4)
					{
						*temp=t2;
						*dens=d2;
						return 0;
					}
					else
						continue;
				}
			}		
		}
	}
//	printf("::::::::%d,%d",t_index,d_index);
	return -1;
}

void print_prop_co2(co2_prop *prop)
{
	printf("%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n",prop->pres,prop->temp,prop->dens,prop->spe_ener,prop->spe_enth,prop->spe_entr,prop->spe_h_v,prop->spe_h_p,prop->speed_sound);
	printf("\n");
}

int co2_prop_calc_pt(co2_prop *prop, co2_prop *prop_bkp, double pres, double temp)
{
	double dens,tsat,densl,densg,psat,densl2,densg2;
	int flag,flag2;
	
	prop->viscous=-1;
	prop->thcond=-1;
	prop_bkp->viscous=-1;
	prop_bkp->thcond=-1;
	
	if(pres<0&&temp<0)
	{
		return -1;
	}
	
	else if(pres<0&&temp>0)
	{
		flag=calc_psat_t(&psat,&densl,&densg,temp);
		if(flag==-1)
			return -1;
		co2_prop_calc_tr(prop,temp,densl);
		co2_prop_calc_tr(prop_bkp,temp,densg);
		prop->vf=0;
		prop_bkp->vf=1;
		return 1;
	}
	
	else if(pres>0&&temp<0)
	{
		flag=calc_tsat_p(&tsat,&densl,&densg,pres);
		if(flag==-1)
			return -1;
		co2_prop_calc_tr(prop,tsat,densl);
		co2_prop_calc_tr(prop_bkp,tsat,densg);
		prop->pres=pres;
		prop_bkp->pres=pres;
		prop->vf=0;
		prop_bkp->vf=1;
		return 1;
	}
	
	else
	{
		flag=calc_tsat_p(&tsat,&densl,&densg,pres);
//		flag2=calc_psat_t(&psat,&densl2,&densg2,temp);
//		printf("\n.....%lf,%lf,%lf,%d\n",tsat,temp,psat,flag);
		if(flag==0&&fabs(tsat-temp)/tsat<1e-3)
		{
			co2_prop_calc_tr(prop,temp,densl);
			co2_prop_calc_tr(prop_bkp,temp,densg);
			prop->pres=pres;
			prop->temp=temp;
			prop_bkp->pres=pres;
			prop_bkp->temp=temp;
			prop->vf=0;
			prop_bkp->vf=1;
//		print_prop_co2(prop);
//		print_prop_co2(prop_bkp);
			return 1;
		}
	
		else
		{
			flag=dens_pt(&dens,pres,temp);
			if(flag==-1)
				return -1;
//			printf("\n....%lf\n",dens);
			co2_prop_calc_tr(prop,temp,dens);
//			printf("........%.12lf\n",prop->pres);
			return 0;
		}
	}
}

int co2_prop_calc_ph(co2_prop *prop, double pres, double spe_enth)
{
	double temp,dens,hl,hg,dl,dg,tsat,x;
	int flag,flag2;
	co2_prop pl,pg;
	
	prop->viscous=-1;
	prop->thcond=-1;
	
	if(pres<P_C)
	{
		flag=calc_tsat_p(&tsat,&dl,&dg,pres);
		if(flag==0)
		{
			co2_prop_calc_tr(&pl,tsat,dl);
			co2_prop_calc_tr(&pg,tsat,dg);
			if(spe_enth<pg.spe_enth&&spe_enth>pl.spe_enth)
			{
				x=(spe_enth-pl.spe_enth)/(pg.spe_enth-pl.spe_enth);
				prop->pres=pres;
				prop->temp=tsat;
				prop->spe_vol=x*pg.spe_vol+(1-x)*pl.spe_vol;
				prop->dens=1/prop->spe_vol;
				prop->spe_entr=x*pg.spe_entr+(1-x)*pl.spe_entr;
				prop->spe_ener=x*pg.spe_ener+(1-x)*pl.spe_ener;
				prop->spe_h_v=x*pg.spe_h_v+(1-x)*pl.spe_h_v;
				prop->spe_enth=spe_enth;
				prop->spe_h_p=x*pg.spe_h_p+(1-x)*pl.spe_h_p;
				prop->sat_liq_hc=x*pg.sat_liq_hc+(1-x)*pl.sat_liq_hc;
				prop->speed_sound=x*pg.speed_sound+(1-x)*pl.speed_sound;
				prop->Joule_Thompson=x*pg.Joule_Thompson+(1-x)*pl.Joule_Thompson;
				prop->fugacity=x*pg.fugacity+(1-x)*pl.fugacity;
				prop->vf=x;
				return 0;
			}
		}
	}
	
	flag=td_ph(&temp,&dens,pres,spe_enth);
	if(flag==-1)
		return -1;
	co2_prop_calc_tr(prop,temp,dens);
	return 0;
}

int co2_prop_calc_ps(co2_prop *prop, double pres, double spe_entr)
{
	double temp,dens,hl,hg,dl,dg,tsat,x;
	int flag,flag2;
	co2_prop pl,pg;
	
	prop->viscous=-1;
	prop->thcond=-1;
		
	if(pres<P_C)
	{
		flag=calc_tsat_p(&tsat,&dl,&dg,pres);
		if(flag==0)
		{
			co2_prop_calc_tr(&pl,tsat,dl);
			co2_prop_calc_tr(&pg,tsat,dg);
			if(spe_entr<pg.spe_entr&&spe_entr>pl.spe_entr)
			{
				x=(spe_entr-pl.spe_entr)/(pg.spe_entr-pl.spe_entr);
				prop->pres=pres;
				prop->temp=tsat;
				prop->spe_vol=x*pg.spe_vol+(1-x)*pl.spe_vol;
				prop->dens=1/prop->spe_vol;
				prop->spe_entr=spe_entr;
				prop->spe_ener=x*pg.spe_ener+(1-x)*pl.spe_ener;
				prop->spe_h_v=x*pg.spe_h_v+(1-x)*pl.spe_h_v;
				prop->spe_enth=x*pg.spe_enth+(1-x)*pl.spe_enth;
				prop->spe_h_p=x*pg.spe_h_p+(1-x)*pl.spe_h_p;
				prop->sat_liq_hc=x*pg.sat_liq_hc+(1-x)*pl.sat_liq_hc;
				prop->speed_sound=x*pg.speed_sound+(1-x)*pl.speed_sound;
				prop->Joule_Thompson=x*pg.Joule_Thompson+(1-x)*pl.Joule_Thompson;
				prop->fugacity=x*pg.fugacity+(1-x)*pl.fugacity;
				prop->vf=x;
				return 0;
			}
		}
	}
	
	flag=td_ps(&temp,&dens,pres,spe_entr);
	if(flag==-1)
		return -1;
	co2_prop_calc_tr(prop,temp,dens);
	return 0;
}

 int co2_prop_calc_hs(co2_prop *prop, double spe_enth, double spe_entr)
{
	double temp,dens,hl,hg,dl,dg,tsat,x,psat;
	int flag,flag2;
	co2_prop pl,pg;

	prop->viscous=-1;
	prop->thcond=-1;
	
	flag=td_hs(&temp,&dens,spe_enth,spe_entr);
	if(flag==0)
	{
		co2_prop_calc_tr(prop,temp,dens);
		return 0;
	}
	flag=hs_find_tsat(&tsat,&x,spe_enth,spe_entr);
	if(flag==-1)
		return -1;
	if(calc_psat_t(&psat,&dl,&dg,tsat)==-1)
		return -1;
	if(co2_prop_calc_pt(&pl,&pg,psat,tsat)==-1)
		return -1;
	prop->pres=psat;
	prop->temp=tsat;
	prop->spe_vol=x*pg.spe_vol+(1-x)*pl.spe_vol;
	prop->dens=1/prop->spe_vol;
	prop->spe_entr=spe_entr;
	prop->spe_ener=x*pg.spe_ener+(1-x)*pl.spe_ener;
	prop->spe_h_v=x*pg.spe_h_v+(1-x)*pl.spe_h_v;
	prop->spe_enth=spe_enth;
	prop->spe_h_p=x*pg.spe_h_p+(1-x)*pl.spe_h_p;
	prop->sat_liq_hc=x*pg.sat_liq_hc+(1-x)*pl.sat_liq_hc;
	prop->speed_sound=x*pg.speed_sound+(1-x)*pl.speed_sound;
	prop->Joule_Thompson=x*pg.Joule_Thompson+(1-x)*pl.Joule_Thompson;
	prop->fugacity=x*pg.fugacity+(1-x)*pl.fugacity;
	prop->vf=x;
	return 0;
}

int co2_prop_calc_pr(co2_prop *prop, double pres, double dens)
{
	double densl,densg,tsat,x,temp;
	int flag;
	co2_prop pl,pg;

	prop->viscous=-1;
	prop->thcond=-1;
		
	if(pres<P_C)
	{
		flag=calc_tsat_p(&tsat,&densl,&densg,pres);
		if(flag==0)
		{
			if(dens>densg&&dens<densl)
			{
				co2_prop_calc_tr(&pl,tsat,densl);
				co2_prop_calc_tr(&pg,tsat,densg);
				x=densg*(densl-dens)/(dens*(densl-densg));
				prop->pres=pres;
				prop->temp=tsat;
				prop->spe_vol=1/dens;
				prop->dens=dens;
				prop->spe_entr=x*pg.spe_entr+(1-x)*pl.spe_entr;
				prop->spe_ener=x*pg.spe_ener+(1-x)*pl.spe_ener;
				prop->spe_h_v=x*pg.spe_h_v+(1-x)*pl.spe_h_v;
				prop->spe_enth=x*pg.spe_enth+(1-x)*pl.spe_enth;
				prop->spe_h_p=x*pg.spe_h_p+(1-x)*pl.spe_h_p;
				prop->sat_liq_hc=x*pg.sat_liq_hc+(1-x)*pl.sat_liq_hc;
				prop->speed_sound=x*pg.speed_sound+(1-x)*pl.speed_sound;
				prop->Joule_Thompson=x*pg.Joule_Thompson+(1-x)*pl.Joule_Thompson;
				prop->fugacity=x*pg.fugacity+(1-x)*pl.fugacity;
				prop->vf=x;
				return 0;				
			}
		}
	}

	flag=temp_pr(&temp,pres,dens);
	if(flag==-1)
		return -1;
	co2_prop_calc_tr(prop,temp,dens);
	return 0;
}

 int co2_prop_calc_th(co2_prop *prop, double temp, double h)
{
	double densl,densg,psat,x,dens;
	int flag;
	co2_prop pl,pg;

	prop->viscous=-1;
	prop->thcond=-1;
		
	if(temp<T_C)
	{
		flag=calc_psat_t(&psat,&densl,&densg,temp);
		if(flag==0)
		{
			co2_prop_calc_tr(&pl,temp,densl);
			co2_prop_calc_tr(&pg,temp,densg);
			if(h>pl.spe_enth&&h<pg.spe_enth)
			{
				x=(h-pl.spe_enth)/(pg.spe_enth-pl.spe_enth);
				prop->pres=psat;
				prop->temp=temp;
				prop->spe_vol=x*pg.spe_vol+(1-x)*pl.spe_vol;
				prop->dens=1/prop->spe_vol;
				prop->spe_entr=x*pg.spe_entr+(1-x)*pl.spe_entr;
				prop->spe_ener=x*pg.spe_ener+(1-x)*pl.spe_ener;
				prop->spe_h_v=x*pg.spe_h_v+(1-x)*pl.spe_h_v;
				prop->spe_enth=h;
				prop->spe_h_p=x*pg.spe_h_p+(1-x)*pl.spe_h_p;
				prop->sat_liq_hc=x*pg.sat_liq_hc+(1-x)*pl.sat_liq_hc;
				prop->speed_sound=x*pg.speed_sound+(1-x)*pl.speed_sound;
				prop->Joule_Thompson=x*pg.Joule_Thompson+(1-x)*pl.Joule_Thompson;
				prop->fugacity=x*pg.fugacity+(1-x)*pl.fugacity;
				prop->vf=x;
				return 0;				
			}
		}
	}
	flag=dens_th(&dens,temp,h);
//	printf("\n%d.....",flag);
	if(flag==-1)
		return -1;
	co2_prop_calc_tr(prop,temp,dens);
	return 0;
}

int co2_prop_calc_ts(co2_prop *prop, double temp, double s)
{
	double densl,densg,psat,x,dens;
	int flag;
	co2_prop pl,pg;
//	printf("\nHERE\n");

	prop->viscous=-1;
	prop->thcond=-1;
		
	if(temp<T_C)
	{
		flag=calc_psat_t(&psat,&densl,&densg,temp);
//		printf("%lf,%lf,%lf\n",densl,densg,psat);
		if(flag==0)
		{
			co2_prop_calc_tr(&pl,temp,densl);
			co2_prop_calc_tr(&pg,temp,densg);
			if(s>pl.spe_entr&&s<pg.spe_entr)
			{
				x=(s-pl.spe_entr)/(pg.spe_entr-pl.spe_entr);
				prop->pres=psat;
				prop->temp=temp;
				prop->spe_vol=x*pg.spe_vol+(1-x)*pl.spe_vol;
				prop->dens=1/prop->spe_vol;
				prop->spe_entr=s;
				prop->spe_ener=x*pg.spe_ener+(1-x)*pl.spe_ener;
				prop->spe_h_v=x*pg.spe_h_v+(1-x)*pl.spe_h_v;
				prop->spe_enth=x*pg.spe_enth+(1-x)*pl.spe_enth;
				prop->spe_h_p=x*pg.spe_h_p+(1-x)*pl.spe_h_p;
				prop->sat_liq_hc=x*pg.sat_liq_hc+(1-x)*pl.sat_liq_hc;
				prop->speed_sound=x*pg.speed_sound+(1-x)*pl.speed_sound;
				prop->Joule_Thompson=x*pg.Joule_Thompson+(1-x)*pl.Joule_Thompson;
				prop->fugacity=x*pg.fugacity+(1-x)*pl.fugacity;
				prop->vf=x;
				return 0;				
			}
		}
	}
//	printf("!!!!!!!!!!!!!!!!!!!");
	flag=dens_ts(&dens,temp,s);
	if(flag==-1)
		return -1;
	co2_prop_calc_tr(prop,temp,dens);
	return 0;
}

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


/*int main()
{
	co2_prop prop,prop2;
//	co2_prop_calc_tr(&prop,T_C+0.00000001,r_c);
//	print_prop_co2(&prop);
//	co2_prop_calc_tr(&prop,T_C+0.00000001,446.456975);
//	print_prop_co2(&prop);
	co2_prop_calc_pt(&prop,&prop2,0.1e6,300);
	print_prop_co2(&prop);
	co2_prop_calc_pt(&prop,&prop2,0.101325e6,298.15);
	print_prop_co2(&prop);
}

/*int main()
{
	co2_prop pp1;
	co2_prop_calc_tr(&pp1,325,773.4575);
	print_prop_co2(&pp1);	
		co2_prop_calc_tr(&pp1,325,773.46);
	print_prop_co2(&pp1);
		co2_prop_calc_tr(&pp1,700,149.27);
	print_prop_co2(&pp1);
		co2_prop_calc_tr(&pp1,370,143.82);
	print_prop_co2(&pp1);
			co2_prop_calc_tr(&pp1,245,1089.93);
	print_prop_co2(&pp1);
			co2_prop_calc_tr(&pp1,325,212.9);
	print_prop_co2(&pp1);

	
	co2_prop_calc_tr(&pp1,490,468.31);
	print_prop_co2(&pp1);
		co2_prop_calc_tr(&pp1,300,865.82);
	print_prop_co2(&pp1);
			co2_prop_calc_tr(&pp1,300,801.62);
	print_prop_co2(&pp1);
				co2_prop_calc_tr(&pp1,300,753.17);
	print_prop_co2(&pp1);
					co2_prop_calc_tr(&pp1,300,733.9);
	print_prop_co2(&pp1);
					co2_prop_calc_tr(&pp1,300,706.06);
	print_prop_co2(&pp1);
						co2_prop_calc_tr(&pp1,300,91.965);
	print_prop_co2(&pp1);
}

/*int main()
{
		co2_prop prop,prop2;
	double pres,temp,densl,densg,dens;
	int flag;
	dens_pt(&dens,30e6,280);
	printf("%lf\n",dens);
	calc_psat_t(&pres,&densl,&densg,224);
	printf("\n\n%lf\n",pres);
	
		calc_tsat_p(&temp,&densl,&densg,0.5e6);
	printf("\n\n%lf\n",temp);
//	co2_prop_calc_tr(&prop,280,720.1141);
//	print_prop_co2(&prop);	
//		co2_prop_calc_tr(&prop,280,1031.04);
//	print_prop_co2(&prop);	
}

/*int main()
{
	co2_prop prop,prop2;
	double pres,temp,densl,densg,dens;
	int flag;
	
	flag=td_hs(&temp,&dens,22.645e3,206.55);
	co2_prop_calc_tr(&prop,temp,dens);
	print_prop_co2(&prop);
	
	flag=td_hs(&temp,&dens,881.58e3,575.27);
	co2_prop_calc_tr(&prop,temp,dens);
	print_prop_co2(&prop);
	
		flag=td_hs(&temp,&dens,-379.22e3,-2038);
	co2_prop_calc_tr(&prop,temp,dens);
	print_prop_co2(&prop);

		flag=td_hs(&temp,&dens,-88944,-237.57);
	co2_prop_calc_tr(&prop,temp,dens);
	print_prop_co2(&prop);
	
			flag=td_hs(&temp,&dens,275.86e3,-683.89);
	co2_prop_calc_tr(&prop,temp,dens);
	print_prop_co2(&prop);

	
/	printf("$$$$$%lf\n",log(2.718281829));
//	co2_prop_calc_tr(&prop,250,1067.96);
//	print_prop_co2(&prop);
	flag=calc_tsat_p(&temp,&densl,&densg,0.51796e6);
	co2_prop_calc_tr(&prop,temp,densl);
	co2_prop_calc_tr(&prop2,temp,densg);
//	printf(";;;;;;;;;;;%d\n",flag);
	print_prop_co2(&prop);
	print_prop_co2(&prop2);	
	
	flag=calc_psat_t(&pres,&densl,&densg,290);
	co2_prop_calc_tr(&prop,290,densl);
	co2_prop_calc_tr(&prop2,290,densg);
//	printf(";;;;;;;;;;;%d\n",flag);
	print_prop_co2(&prop);
	print_prop_co2(&prop2);	
	for(temp=190;temp<1110;temp=temp+10)
	{
		flag=dens_pt(&dens,0.05e6,temp);
		if(flag==-1)
		{
			printf("\n! ERRRRR!\n");
		}	
		co2_prop_calc_tr(&prop,temp,dens);
		print_prop_co2(&prop);
	}
}*/
