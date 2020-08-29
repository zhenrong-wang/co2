//units: pressure-Pa, temp-K, spe_vol-m3/kg, dens-kg/m3, spe_ener-J/kg, spe_entr-J/kgK, spe_enth-J/kg, spe_h-J/kgK


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

#endif