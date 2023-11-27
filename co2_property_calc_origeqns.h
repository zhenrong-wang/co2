/*
 *This code is distributed under the license: MIT License
 * Originally written by Zhenrong WANG
 * mailto: zhenrongwang@live.com
 */

#ifndef CO2_PROPERTY_CALC_ORIGEQNS_H_
#define CO2_PROPERTY_CALC_ORIGEQNS_H_

double phi0(double del, double tau);
double phi0d(double del);
double phi0dd(double del);
double phi0t(double tau);
double phi0tt(double tau);
double deltd(double del, double delta, double aai, double theta, double bti, double bbi, double ai);
double deltdd(double dd, double del, double delta, double aai, double theta, double bti, double bbi, double ai);
double deltabid(double bi, double delta, double dd);
double deltabidd(double bi, double delta, double dd, double ddd);
double deltabit(double bi, double theta, double delta);
double deltabitt(double bi, double theta, double delta);
double deltabidt(double aai, double bi, double bti, double delta, double del, double theta, double dd);
double psid(double del, double psi, double cci);
double psidd(double del, double psi, double cci);
double psit(double tau, double psi, double ddi);
double psitt(double tau, double psi, double ddi);
double psidt(double tau, double del, double psi, double cci, double ddi);
double phir(double del, double tau);
double phird(double del, double tau);
double phirdd(double del, double tau);
double phirt(double del, double tau);
double phirtt(double del, double tau);
double phirdt(double del, double tau);
double p_tr(double temp, double dens);
double s_tr(double temp, double dens);
double u_tr(double temp, double dens);
double cv_tr(double temp, double dens);
double h_tr(double temp, double dens);
double cp_tr(double temp, double dens);
double cs_tr(double temp, double dens);
double w_tr(double temp, double dens);
double miu_tr(double temp, double dens);
double f_tr(double temp, double dens);

#endif