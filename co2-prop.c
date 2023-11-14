#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "co2_property_calc.c"

int welcome(void)
{
    time_t rtime;
    struct tm* timeinfo=NULL;
    time(&rtime);
    timeinfo=localtime(&rtime);
	printf("/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */\n");
	printf("/*                    CO2 Properties Calculation.                     */\n");
	printf("/*                   VERSION 1.0(beta) LICENSE: MIT                   */\n");
	printf("/* Author: Zhenrong WANG. zhenrongwang@live.com; K495458966(wechat).  */\n");
	printf("/* Main function: calculating properties by given two values:         */\n");
	printf("/*       1-pt, 2-ph, 3-ps, 4-hs, 5-tr, 6-pr, 7-th, 8-ts.              */\n"); 
	printf("/*       Input File: _input.dat; Output File: _properties.dat.        */\n");
	printf("/* Range: 216.15K<T<1100K     && P<=800MPa;                           */\n");
	printf("/* NOTE: For detailed information, please read the help doc.          */\n");
	printf("/*       Any bugs found, please contact Zhenrong WANG.                */\n");
	printf("/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */\n");
	printf("\n# CURRENT DATE AND TIME: %s\n",asctime(timeinfo));
	return 0;
}

int main(int argc, char** argv)
{
	co2_prop prop,prop_bkup;
	FILE* fin;
	FILE* fout;
	int type,flag,flag2,i=0;
	char sep1,sep2,enter;
	double v1,v2;

	welcome();
	fin=fopen("_input.dat","r");
	if(fin==NULL)
	{
		printf("\n! FATAL ERROR: input file not found.\n! Program abort.\n! Please press any key to exit.");
		printf("\n@ Any problems found, please contact the author.\n@ Zhenrong Wang, zhenrongwang@live.com, K495458966(wechat).\n@ All rights reserved.\n");		
		fflush(stdin);
		getchar();
		return -1;
	}
	fout=fopen("_properties.dat","w");
	if(fout==NULL)
	{
		printf("\n! FATAL ERROR: cannot creat output file.\n! Program abort.\n! Please press any key to exit.");
		printf("\n@ Any problems found, please contact the author.\n@ Zhenrong Wang, zhenrongwang@live.com, K495458966(wechat).\n@ All rights reserved.\n");
		fclose(fin);
		fflush(stdin);
		getchar();
		return -1;
	}
	
	printf("\n# Please press ENTER to start calculation:");
	fflush(stdin);
	getchar(); 
	
	fprintf(fout,"\tPRES\t\tTEMP\t\tDENS\t\tu\t\th\t\ts\t\tCv\t\tCp\t\tVsound\t\tV_FRAC\t\tVISC\t\tTHCOND\n");
	fprintf(fout,"\tMPa\t\tK\t\tkg/m^3\t\tkJ/kg\t\tkJ/kg\t\tkJ/(kg*K)\tkJ/(kg*K)\tkJ/(kg*K)\tm/s\t\t-\t\t10^-6Pa.s\tmW/(m.K)\n");
	while(!feof(fin))
	{
		i++;
		fscanf(fin,"%d%c%lf%c%lf%c",&type,&sep1,&v1,&sep2,&v2,&enter);
		printf("\n%d\t%.8lf\t%.8lf",type,v1,v2);
	
		if(type==1)
			flag=co2_prop_calc_pt(&prop,&prop_bkup,v1,v2);
		else if(type==2)
			flag=co2_prop_calc_ph(&prop,v1,v2);
		else if(type==3)
			flag=co2_prop_calc_ps(&prop,v1,v2);
		else if(type==4)
			flag=co2_prop_calc_hs(&prop,v1,v2);
		else if(type==5)
		{
			flag=0;
			co2_prop_calc_tr(&prop,v1,v2);
		}
		else if(type==6)
			flag=co2_prop_calc_pr(&prop,v1,v2);
		else if(type==7)
			flag=co2_prop_calc_th(&prop,v1,v2);
		else if(type==8)
			flag=co2_prop_calc_ts(&prop,v1,v2);
		else
		{
			continue;
		}
		if(flag==-1)
		{
			printf("\n! WARNING: Calculation error at line # %d of the input file.\n", i);
			fprintf(fout,"%d\t-1\n",i);
			type=-1;
			continue;
		}
		else if(flag==1)
		{
			calc_trans_prop(&prop);
			calc_trans_prop(&prop_bkup);
			fprintf(fout,"%d\t%8.4lf*\t%8.4lf\t%10.4lf\t%8.4lf\t%8.4lf\t%8.4lf\t%8.4lf\t%10.5lf\t%10.5lf\t%8.4lf\t%8.6lf\t%8.6lf\n",i,prop.pres/1e6,prop.temp,prop.dens,prop.spe_ener/1000,prop.spe_enth/1000,prop.spe_entr/1000,prop.spe_h_v/1000,prop.spe_h_p/1000,prop.speed_sound,prop.vf,prop.viscous,prop.thcond);
			fprintf(fout,"\t%8.4lf*\t%8.4lf\t%10.4lf\t%8.4lf\t%8.4lf\t%8.4lf\t%8.4lf\t%10.5lf\t%10.5lf\t%8.4lf\t%8.6lf\t%8.6lf\n",prop_bkup.pres/1e6,prop_bkup.temp,prop_bkup.dens,prop_bkup.spe_ener/1000,prop_bkup.spe_enth/1000,prop_bkup.spe_entr/1000,prop_bkup.spe_h_v/1000,prop_bkup.spe_h_p/1000,prop_bkup.speed_sound,prop_bkup.vf,prop_bkup.viscous,prop_bkup.thcond);
			type=-1;
		}
		else
		{
			calc_trans_prop(&prop);
			fprintf(fout,"%d\t%8.4lf\t%8.4lf\t%10.4lf\t%8.4lf\t%8.4lf\t%8.4lf\t%8.4lf\t%10.5lf\t%10.5lf\t%8.4lf\t%8.6lf\t%8.6lf\n",i,prop.pres/1e6,prop.temp,prop.dens,prop.spe_ener/1000,prop.spe_enth/1000,prop.spe_entr/1000,prop.spe_h_v/1000,prop.spe_h_p/1000,prop.speed_sound,prop.vf,prop.viscous,prop.thcond);
			type=-1;
		}	
	}	
	
	printf("\n\n# Calculation finished. Please check _properties.dat for results.\n# Press any key to exit.\n\n@ Any problems found, please contact the author.\n@ Zhenrong Wang, zhenrongwang@live.com, K495458966(wechat).\n@ All rights reserved.\n");
	fclose(fin);
	fclose(fout);
	fflush(stdin);
	getchar();
	return 0;
	
}

