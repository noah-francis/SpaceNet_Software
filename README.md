## SpaceNet Matlab Application
---
Overall description of reposity goes here.

#### **TLE_DATA_V_1_1.mlapp**
Version 1.1 of the TLE Data app

#### **TLE_DATA_V_1_2.mlapp**
......

#### **TLE_DATA_V_1_3.mlapp**
......

#### **gibbs.m**
MATLAB function that runs Gibb's Method to produce a velocity vector from three consecutive position vectors given in the ECEF frame.
  
**_Inputs_**:

'r_vec' -- array of three time-consecutive position vectors $[km]$

'r_mag' -- array of corresponding vector magnitides $[km]$

'mu'    -- Standard Gravitational Parameter of the body orbited $[km^3 / s^2]$

**_Outputs_**:

'v_vec'        -- velocity vector at the second time $[km/s]$

'coplanar_val' -- dot product of unit vectors that should be coplanar for Gibb's Method to work properly $[-]$ (closer to zero is better)

#### **read_STK**
Read and parse an STK ECEF satellite position file (assumes data directory)


**_Inputs_**:

'filename' -- string of the file name

**_Outputs_**:

'data' -- struct containing the STK file info