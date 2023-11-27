/*
 * This code is distributed under the license: MIT License
 * Originally written by Zhenrong WANG
 * mailto: zhenrongwang@live.com
 */

// units: pressure-Pa, temp-K, spe_vol-m3/kg, dens-kg/m3, spe_ener-J/kg, spe_entr-J/kgK, spe_enth-J/kg, spe_h-J/kgK

#ifndef CO2_PROPERTY_CALC_H_
#define CO2_PROPERTY_CALC_H_

#define MAX_PRES 800e6
#define MAX_TEMP 1100
#define MIN_PRES 0
#define MIN_TEMP 216

#define T_C 304.1282 
#define P_C 7.377324711212e6
#define r_c 467.6 

#define T_t 216.592
#define P_t 0.51795e6

#define T_0 298.15
#define p_0 0.101325e6
#define h_0 0
#define s_0 0

#define GAS_CONST_CO2 188.9241

typedef struct{
	double pres;
	double temp;
	double spe_vol;
	double dens;
	double spe_entr;
	double spe_ener;
	double spe_h_v;
	double spe_enth;
	double spe_h_p;
	double sat_liq_hc;
	double speed_sound;
	double Joule_Thompson;
	double fugacity;
	double vf;
	double viscous;
	double thcond; 
//	double vc2nd;
//	double vc3rd;
}co2_prop;

int is_vld_prop(co2_prop *prop);
void co2_prop_calc_tr(co2_prop *prop, double temp, double dens);
void swap_line(double left_matrix[3][3],double right_array[3],int i,int j);
int Gaussian_elimination(double *solution, double left_matrix[3][3],double right_hand_array[3]);
double find_max(double* solu);
void calcfgh_dfgh_psat(double dfgh[3][3], double fgh[3], double tau, double d1, double d2, double ps);
int calc_tsat_p(double *tsat, double *densl, double *densg, double psat);
void calcfgh_dfgh_tsat(double dfgh[3][3], double fgh[3], double ps, double d1, double d2, double ts);
int calc_psat_t(double *psat, double *densl, double *densg, double tsat);
void calcfd_dfd_pt(double *f, double *df, double del, double pres, double tau);
int dens_pt(double *dens, double pres, double temp);
void calcfd_dfd_pr(double *f, double *df, double tau, double pres, double del);
int temp_pr(double *temp, double pres, double dens);
void calcfd_dfd_th(double *f, double *df, double del, double tau, double h);
int dens_th(double* dens, double temp, double h);
void calcfd_dfd_ts(double *f, double *df, double del, double tau, double s);
int dens_ts(double* dens, double temp, double s);
void calcfg_dfg_hs(double *f, double *g, double *dfdt, double *dfdd, double *dgdt, double *dgdd, double tau, double del, double h, double s);
int td_hs(double *temp, double *dens, double h, double s);
int hs_find_tsat(double *t, double *x, double spe_enth, double spe_entr);
void calcfg_dfg_ph(double *f, double *g, double *dfdt, double *dfdd, double *dgdt, double *dgdd, double tau, double del, double p, double h);
int td_ph(double *temp, double *dens, double p, double h);
void calcfg_dfg_ps(double *f, double *g, double *dfdt, double *dfdd, double *dgdt, double *dgdd, double tau, double del, double p, double s);
int td_ps(double *temp, double *dens, double p, double s);
void print_prop_co2(co2_prop *prop);
int co2_prop_calc_pt(co2_prop *prop, co2_prop *prop_bkp, double pres, double temp);
int co2_prop_calc_ph(co2_prop *prop, double pres, double spe_enth);
int co2_prop_calc_ps(co2_prop *prop, double pres, double spe_entr);
int co2_prop_calc_hs(co2_prop *prop, double spe_enth, double spe_entr);
int co2_prop_calc_pr(co2_prop *prop, double pres, double dens);
int co2_prop_calc_th(co2_prop *prop, double temp, double h);
int co2_prop_calc_ts(co2_prop *prop, double temp, double s);
double eta0(double temp);
double deleta(double temp, double dens);
double visc_calc(double temp, double dens);
double lmd0(double temp);
double deltlmd(double temp, double dens);
double deltlmdc_emp(double temp, double dens);
double drdp(double temp, double dens);
double deltlmdc(co2_prop* prop, double temp, double dens);
void calc_trans_prop(co2_prop* prop);

#endif