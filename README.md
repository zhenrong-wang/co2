# CO2-Carbon-Dioxide-Property-Calculation

# 1. Background

The thermophysical properties of CO2 are very important to researchers and engineers in various domains (energy, power, chemistry, materials ...).

Usually, researchers and engineers rely on the [NIST Refprop](https://www.nist.gov/srd/refprop) to calculate the properties. However, Refprop is commercial and not open-source. Therefore, I wrote this code to provide an alternative to the refprop.

# 2. Brief Intro

Program Name: High-resolution CO2 (Carbon Dioxide) Properties Calculation

Purpose: Calculate the **very-high-resolution** thermophysical properties at given point(s)

License: MIT

Technical Reference:

- R. Span, W. Wagner, A New Equation of State for Carbon Dioxide Covering the Fluid Region from the Triple –Point Temperature to 1100K at Pressures up to 800MPa, J. Phys. Chem. Ref. Data, Vol 25, No.6, 1996
- Marcia L Huber, Marc J Assel, R. A. Perkins, Reference Correlation of the Thermal Conductivity of Carbon Dioxide from the Triple Point to 1100K and up to 200MPa. Joursssnal of Physical and Chemical References Date. Feb. 2016
- William A. Wakeham, The Viscosity of Carbon Dioxide, Journal of Physical and Chemical Reference Data, Jan, 1998.

# 3. How-To: Build, Run, and Use

## 3.1 Build

### 3.1.1 Prerequisites

You need a C compiler to build. 

- For Microsoft Windows users, [mingw](https://sourceforge.net/projects/mingw/) is a good choice
- For GNU/Linux Distro or other *nix users, the [GNU Compiler Collections](https://gcc.gnu.org/), known as gcc, is a perfect one
- For macOS users, [clang](https://clang.llvm.org/) is easy to install and use (brew is not needed to install clang on macOS).

### 3.2 Build Guide

1. Use git to clone this code: `git clone https://github.com/zhenrong-wang/co2.git`
2. Build command example: `gcc co2-prop.c -o my-co2.exe -lm`

## 3.2 Run

- Copy the executable to a working directory. and open a terminal/command prompt
- Create an input file named '**_input.dat**' with some input data (See Section 3.3 below)
- Run the executable in the terminal or command prompt by the command `./my-co2.exe` (for *nix users) or `.\my-co2.exe` (for Microsoft Windows users)

## 3.3 Use

For the **_input.dat** file, you need to input the data in lines. A single line refers to a point. Strict format:

`POINT_TYPE(int),PARAM1(float/double),PARAM2 (float/double)` 

**Please separate the params by a comma ','

POINT_TYPE: 

- 1 Pressure and Temperature
- 2 Pressure and Specific Enthalpy
- 3 Pressure and Specific Entropy
- 4 Specific Enthalpy and Specific Entropy
- 5 Temperature and Density
- 6 Pressure and Density
- 7 Temperature and Specific Enthalpy
- 8 Temperature and Specific Entropy

Example: `1,100000,300`

Example: `1,1e5,3e2`

Units:

- Pressure: Pa
- Temperature: K
- Density: kg/m3
- Specific Enthalpy: J/kg
- Specific Entropy: J/(kg.K)

## 3.4 Output

The program writes out a file named '**_properties.dat**' in the working directory. Here are some extra properties:

- Vsound: Velocity of sound
- V_FRAC: The mass ratio of gas when the point is in dual-phase
- VISC: Viscosity
- THCOND: Thermal Conductivity ratio

# 4 Bugs and Communications

This code was written years ago, I didn't continuously develop it. Therefore, bugs may occur. 

If you are interested in this project, please submit issues to this repo. I'd be glad to communicate on any issues.

Or, you can also email me: zhenrongwang@live.com
