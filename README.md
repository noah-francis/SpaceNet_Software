# SpaceNet Software
CU Boulder Aerospace Engineering Senior Projects Team SpaceNet software development.

---
## Making .mlapp Changes
1) Make a copy of the latest version of the app
2) Make changes on copied version of the app on separate branch
3) Update version number. First number in version corresponds to a graphical update, second corresponds to a software update.
4) When finished, stage and commit new .mlapp version
5) Then push branch to github with 'git push origin <your_branch_name>'
6) Note: editing a past version or even running will cause merging errors, so make sure to make copies.

## .mlapp Version History
#### **TLE_DATA_V_1_1.mlapp**
New format includes TLE comparison tab

#### **TLE_DATA_V_1_2.mlapp**
Gibb's method and orbital element calculation functions implemented

#### **TLE_DATA_V_1_3.mlapp**
Implemented a TLE text file header parser (needs work)

#### **TLE_DATA_V_1_4.mlapp**
Implemented first functional sensitivity analysis 

#### **TLE_DATA_V_1_5.mlapp**
Has true velocity vector compared to perturbed velocity vector plotted in Sensitivity Analysis tab