%nprocshared=16
%mem=32GB
# wb97xd 6-31g* opt=modredundant freq=noraman

truncated_active_site_radius_5_conf_1_geom_opt

-1 1
 H  38.93600000  34.34600000  -3.27400000
 C  39.83300000  33.42000000  -2.56400000
 H  39.90900000  32.07700000  -3.26400000
 C  39.40500000  33.19500000  -1.10300000
 C  37.89900000  32.91600000  -0.86300000
 C  37.13100000  34.19600000  -0.53400000
 O  36.45400000  34.29500000   0.49200000
 N  37.26500000  35.19300000  -1.39000000
 H  30.97300000  31.86900000   8.80000000
 C  30.63300000  32.69900000   7.65300000
 C  31.78400000  32.96900000   6.70200000
 O  32.90100000  33.32000000   7.14400000
 N  31.52900000  32.78500000   5.40000000
 C  32.38500000  33.30900000   4.34500000
 C  31.47400000  33.85900000   3.23600000
 O  30.69000000  33.11800000   2.64000000
 C  33.38400000  32.24300000   3.77800000
 C  34.09300000  32.77000000   2.53400000
 C  34.41400000  31.90200000   4.79800000
 N  31.56700000  35.15300000   2.96900000
 C  30.70400000  35.81700000   1.97700000
 C  30.32700000  37.21300000   2.44900000
 O  30.94000000  37.73300000   3.37300000
 N  29.29000000  37.80400000   1.87200000
 C  28.99700000  39.20800000   2.09900000
 C  28.11000000  39.38000000   3.30900000
 O  27.25200000  38.54400000   3.58900000
 H  28.30900000  39.83100000   0.89500000
 N  28.31000000  40.48900000   4.01200000
 C  27.58100000  40.81100000   5.23600000
 H  26.05400000  40.64300000   5.13800000
 H  27.91300000  42.24900000   5.59500000
 H  35.32000000  34.60400000  11.25600000
 C  35.39900000  35.27200000   9.93000000
 C  34.24200000  36.26900000   9.73400000
 O  34.17200000  37.29400000  10.42900000
 C  35.57200000  34.26500000   8.74400000
 S  36.10200000  35.00000000   7.02600000
 N  33.31800000  35.97900000   8.82300000
 C  32.16000000  36.87300000   8.63800000
 H  31.30900000  36.91700000   9.91600000
 C  31.32100000  36.47000000   7.41700000
 C  30.05300000  37.31700000   7.16400000
 C  30.39100000  38.81500000   6.91400000
 C  29.09700000  36.74700000   6.04100000
 H  32.23100000  39.88600000  12.12000000
 C  31.68100000  41.10900000  11.59100000
 H  30.18200000  41.22200000  11.89400000
 C  31.94400000  41.16300000  10.09000000
 C  33.42000000  41.19900000   9.73500000
 C  33.64600000  41.35000000   8.26600000
 O  32.68900000  41.47800000   7.49800000
 N  34.92200000  41.34500000   7.84700000
 H  30.25800000  46.91800000   3.79400000
 C  30.00400000  45.85800000   2.81100000
 H  28.51900000  45.65600000   2.51800000
 C  30.61300000  44.52300000   3.25200000
 C  30.54000000  43.43100000   2.17900000
 C  31.23500000  42.13100000   2.57000000
 O  31.01400000  41.62800000   3.68000000
 O  31.97100000  41.57400000   1.72700000
 H  41.27900000  35.53000000   5.54700000
 C  40.79900000  35.93700000   4.20200000
 C  39.65600000  36.92300000   4.36800000
 O  38.53200000  36.51300000   4.64800000
 C  40.26000000  34.80200000   3.29700000
 C  39.63200000  35.43600000   2.01200000
 C  41.35800000  33.81300000   2.88800000
 N  39.93800000  38.21200000   4.19700000
 C  38.91400000  39.25000000   4.36900000
 C  39.34800000  40.59600000   3.80900000
 O  40.54400000  40.96700000   3.80400000
 C  38.53800000  39.41200000   5.82800000
 C  39.63000000  40.00300000   6.67000000
 N  39.52200000  41.23800000   7.27400000
 C  40.84900000  39.52200000   7.01300000
 C  40.63100000  41.49600000   7.94300000
 N  41.45200000  40.46800000   7.80400000
 N  38.34000000  41.30800000   3.33600000
 C  38.47000000  42.62100000   2.78200000
 C  38.10300000  43.68700000   3.81800000
 O  38.51600000  44.83400000   3.68600000
 C  37.57400000  42.69300000   1.55700000
 O  38.07100000  41.77500000   0.58900000
 C  37.51300000  44.09300000   0.92800000
 N  37.36000000  43.30400000   4.86200000
 C  36.79100000  44.26400000   5.80800000
 H  37.29300000  44.04100000   7.22600000
 C  35.23500000  44.24500000   5.74800000
 C  34.67600000  44.77900000   4.43800000
 C  34.65300000  43.98600000   3.28600000
 C  34.18600000  46.09000000   4.35400000
 C  34.14700000  44.49400000   2.05200000
 C  33.70600000  46.61100000   3.14700000
 C  33.68100000  45.81400000   2.01100000
 O  33.21400000  46.35800000   0.84500000
 N  33.00000000  39.52000000   3.38700000
 C  33.85900000  39.90900000   4.53000000
 C  35.24900000  40.39200000   4.09200000
 O  35.73300000  40.08500000   2.98500000
 C  33.94800000  38.78200000   5.57600000
 C  34.90800000  37.59700000   5.30400000
 C  34.71000000  36.89400000   3.96600000
 O  33.58100000  36.60500000   3.53000000
 N  35.83200000  36.60600000   3.29500000
 O  35.90200000  41.12000000   4.87200000
 O  31.28300000  42.91200000   6.06600000
 O  33.61300000  38.76000000   0.74200000
 O  36.16700000  38.93200000   0.59400000

X 2 F 
X 10 F 
X 14 F 
X 21 F 
X 25 F 
X 30 F 
X 34 F 
X 40 F 
X 47 F 
X 55 F 
X 63 F 
X 70 F 
X 80 F 
X 87 F 
X 98 F 


