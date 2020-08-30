#---------------Transient optimization of a gas pipeline system----------------
using JuMP
using Ipopt
using XLSX
xf = XLSX.readxlsx("tra.xlsx")
Profile_P = xf["P and M profiles"]["B2:B96"]
Profile_M = xf["P and M profiles"]["C2:C96"]
D4 = xf["Demand variations"]["B1:G49"] 
w1= 1; #w1 is the weight of two objectives, which can be adjusted
# the unit of flow rate is kg/s, the unit used in paper is 104m3/d
e1 = 0.35;
e2 = 0.3;
D1=10;D2=10;D3=20;
d=1.016-17.4*2*10^-3;
Ke=0.017*10^-3;
T=293;
ZMui=16.3776;
a=356.6; #wave speed
ratespeed=9030;  # rated rotational speed of compressors
minflow=2.9113;  # surge flow rate of compressors
maxflow=6.3333;  # stone flow rate of compressors
QVSH=[-243.41 684.18 8677.2];  #coefficient of compressors
QVSN=[-2.38 18.96 47.72]; #coefficient of compressors
A=0.25*3.14*d*d; #pipe segment area
dt=0.5*60*60 #set the time step as 0.5 hour
nt=49 # the time period is 48h
nt1=13;
ns=6;#number of scenarios
dx=10*1000 #set the space step as 10 km
A01=200.326827108680;B01=0.0444;C01=2212409.40302635;D01=37369729.1764816;E01=1308725230.88976;ZMui=16.3776;
aa1=7.7925;bb1=0.0053;cc1=302203.851412505;dd1=241.3268;rr1=0.0054;arr1=6.90617229223050*10^-5;Rm = 8.3143;
gp = Model(Ipopt.Optimizer)
hed_st=[10000,3700,3900,4000,3000];
#----------------------defining variable-------------------------
@variable(gp, lam>=0, start=10) #defining the friction factor
Pl = 4*10^3*ones(95)
Pl[17]=6.5*10^3
Pl[34]=6.5*10^3
Pl[51]=6.5*10^3
Pl[68]=6.5*10^3
#first stage
@variable(gp, 6.5*10^3<=P1[1:5, 1:nt1]<=10*10^3, start=10000) #defining pressure variations
@variable(gp, M1[1:nt1]>=0, start=360) #defining mass flow rate variations
@variable(gp, 30<=UGS1[1:nt1]<=60, start=40) #defining gas withdrawn variations kg/s
#second stage
@variable(gp, Pl[i]<=P[i in 1:95, 1:nt, 1:ns]<=10*10^3, start=Profile_P[i]) #defining pressure variations
@variable(gp, M[i in 1:95, 1:nt, 1:ns]>=0, start=Profile_M[i]) #defining mass flow rate variations
@variable(gp, 30<=UGS[1:nt, 1:ns]<=60, start=40) #defining gas withdrawn variations kg/s
@variable(gp, 0<=Z[1:4, 1:nt, 1:ns]<=1, start=0.9)
@variable(gp, den[1:4, 1:nt, 1:ns]>=0, start=4)
@variable(gp, kv[1:4, 1:nt, 1:ns]>=0, start=1.4)
@variable(gp, hed[i=1:5, 1:nt, 1:ns]>=0, start=hed_st[i])
@variable(gp, 6000<=rspeed1[1:nt, 1:ns]<=10500,start=8250)
@variable(gp, 5160<=rspeed2[1:nt, 1:ns]<=9030,start=6454)
@variable(gp, 5160<=rspeed3[1:nt, 1:ns]<=9030,start=6454)
@variable(gp, 5160<=rspeed4[1:nt, 1:ns]<=9030,start=6403)
@variable(gp, 5160<=rspeed5[1:nt, 1:ns]<=9030,start=5953)
@variable(gp, 0<=xiaolv[1:5, 1:nt, 1:ns]<=100, start=85)
@variable(gp, w)
#----------------------------Objective ---------------------- -----------------
@NLobjective(gp, Min, w)
#---------------------------Robust constraint---------------------------------
@NLconstraint(gp, [k = 1:ns], w - (w1*((sum(4*100*hed[1,j,k]*9.8*(M[1,j,k]/4)/1000/xiaolv[1,j,k] for j in 1:nt-1)
    + sum(2*100*hed[2,j,k]*9.8*(M[18,j,k]/2)/1000/xiaolv[2,j,k] for j in 1:nt-1)
    + sum(2*100*hed[3,j,k]*9.8*(M[35,j,k]/2)/1000/xiaolv[3,j,k] for j in 1:nt-1)
    + sum(2*100*hed[4,j,k]*9.8*(M[52,j,k]/2)/1000/xiaolv[4,j,k] for j in 1:nt-1)
    + sum(2*100*hed[5,j,k]*9.8*(M[69,j,k]/2)/1000/xiaolv[5,j,k] for j in 1:nt-1))*0.5*e1 +
    e2*sum(UGS[j,k]*30*60/(0.0417*16.3776) for j in 1:nt-1))-(1-w1)*(sum((P[i,nt,k]+P[i+1,nt,k])/2 for i in 1:16)*(10^7)*A/a/a +
    sum((P[i,nt,k]+P[i+1,nt,k])/2 for i in 18:33)*(10^7)*A/a/a + sum((P[i,nt,k]+P[i+1,nt,k])/2 for i in 35:50)*(10^7)*A/a/a +
    sum((P[i,nt,k]+P[i+1,nt,k])/2 for i in 52:67)*(10^7)*A/a/a + sum((P[i,nt,k]+P[i+1,nt,k])/2 for i in 69:94)*(10^7)*A/a/a)/(0.0417*16.3776)/10)/1000>=0)

#-------------------------- Initial values------------------------------------
@constraint(gp, [i = 1:95, k = 1:ns], M[i, 1, k] - Profile_M[i] == 0)
@constraint(gp, [i = 1:95, k = 1:ns], P[i, 1, k] - Profile_P[i] == 0)
@constraint(gp, [k = 1:ns], UGS[1, k] - 40 == 0)
#-------------------------- Working efficiency of compressors------------------------------------

@NLconstraint(gp, [j=1:nt, k = 1:ns], -13.91*(M[1,j,k]/(2.2599*ZMui)/4)^2*(8000/rspeed1[j,k])^2+
68.03*(M[1,j,k]/(2.2599*ZMui)/4)*(8000/rspeed1[j,k]) + 2.08 - xiaolv[1,j,k]==0)

@NLconstraint(gp, [j=1:nt, k = 1:ns], QVSN[1]*(M[18,j,k]/(ZMui*den[1,j,k])/2)^2*(ratespeed/rspeed2[j,k])^2+
QVSN[2]*(M[18,j,k]/(ZMui*den[1,j,k])/2)*(ratespeed/rspeed2[j,k]) + QVSN[3]- xiaolv[2,j,k]==0)

@NLconstraint(gp, [j=1:nt, k = 1:ns], QVSN[1]*(M[35,j,k]/(ZMui*den[2,j,k])/2)^2*(ratespeed/rspeed3[j,k])^2+
QVSN[2]*(M[35,j,k]/(ZMui*den[2,j,k])/2)*(ratespeed/rspeed3[j,k]) + QVSN[3]- xiaolv[3,j,k]==0)

@NLconstraint(gp, [j=1:nt, k = 1:ns], QVSN[1]*(M[52,j,k]/(ZMui*den[3,j,k])/2)^2*(ratespeed/rspeed4[j,k])^2+
QVSN[2]*(M[52,j,k]/(ZMui*den[3,j,k])/2)*(ratespeed/rspeed4[j,k]) + QVSN[3]- xiaolv[4,j,k]==0)

@NLconstraint(gp, [j=1:nt, k = 1:ns], QVSN[1]*(M[69,j,k]/(ZMui*den[4,j,k])/2)^2*(ratespeed/rspeed5[j,k])^2+
QVSN[2]*(M[69,j,k]/(ZMui*den[4,j,k])/2)*(ratespeed/rspeed5[j,k]) + QVSN[3]- xiaolv[5,j,k]==0)

#----------------------pyhsical parameter at inlet of the compressor-------------------------
@NLconstraint(gp, [j=1:nt, k = 1:ns], kv[1,j,k]-1.4275==0)
@NLconstraint(gp, [j=1:nt, k = 1:ns], kv[2,j,k]-1.4275==0)
@NLconstraint(gp, [j=1:nt, k = 1:ns], kv[3,j,k]-1.4275==0)
@NLconstraint(gp, [j=1:nt, k = 1:ns], kv[4,j,k]-1.4336==0)
#CS2 in
@NLconstraint(gp, [j=1:nt, k=1:ns], den[1,j,k]*Rm*293+(B01*Rm*293-A01-C01/(293^2)+D01/(293^3)-E01/(293^4))*(den[1,j,k]^2)+(bb1*Rm*(293)-aa1-dd1/(293))*(den[1,j,k]^3)+
arr1*(aa1+dd1/293)*(den[1,j,k]^6) +(cc1*(den[1,j,k]^3)/(293^2))*(1+rr1*(den[1,j,k]^2))*exp(-rr1*(den[1,j,k]^2))-P[17,j,k] == 0)

@NLconstraint(gp, [j=1:nt, k=1:ns], den[1,j,k]*Rm*(293)*Z[1,j,k]==P[17,j,k])

#CS3 in
@NLconstraint(gp, [j=1:nt, k=1:ns], den[2,j,k]*Rm*293+(B01*Rm*293-A01-C01/(293^2)+D01/(293^3)-E01/(293^4))*(den[2,j,k]^2)+(bb1*Rm*(293)-aa1-dd1/(293))*(den[2,j,k]^3)+
arr1*(aa1+dd1/293)*(den[2,j,k]^6) +(cc1*(den[2,j,k]^3)/(293^2))*(1+rr1*(den[2,j,k]^2))*exp(-rr1*(den[2,j,k]^2))-P[34,j,k] == 0)

@NLconstraint(gp, [j=1:nt, k=1:ns], den[2,j,k]*Rm*(293)*Z[2,j,k]==P[34,j,k])

#CS4 in
@NLconstraint(gp, [j=1:nt, k=1:ns], den[3,j,k]*Rm*293+(B01*Rm*293-A01-C01/(293^2)+D01/(293^3)-E01/(293^4))*(den[3,j,k]^2)+(bb1*Rm*(293)-aa1-dd1/(293))*(den[3,j,k]^3)+
arr1*(aa1+dd1/293)*(den[3,j,k]^6) +(cc1*(den[3,j,k]^3)/(293^2))*(1+rr1*(den[3,j,k]^2))*exp(-rr1*(den[3,j,k]^2))-P[51,j,k] == 0)

@NLconstraint(gp, [j=1:nt, k=1:ns], den[3,j,k]*Rm*(293)*Z[3,j,k]==P[51,j,k])

#CS5 in
@NLconstraint(gp, [j=1:nt, k=1:ns], den[4,j,k]*Rm*293+(B01*Rm*293-A01-C01/(293^2)+D01/(293^3)-E01/(293^4))*(den[4,j,k]^2)+(bb1*Rm*(293)-aa1-dd1/(293))*(den[4,j,k]^3)+
arr1*(aa1+dd1/293)*(den[4,j,k]^6) +(cc1*(den[4,j,k]^3)/(293^2))*(1+rr1*(den[4,j,k]^2))*exp(-rr1*(den[4,j,k]^2))-P[68,j,k] == 0)

@NLconstraint(gp, [j=1:nt, k=1:ns], den[4,j,k]*Rm*(293)*Z[4,j,k]==P[68,j,k])
#----------------------Nodes-------------------------
@NLconstraint(gp, lam+2*log10(Ke/3.7/d)==0) #equation of the friction factor, actual lam=1/(lam^2)

#CS1 Pin=5000KPa,Tin=293K,Z=0.9082
@NLconstraint(gp,[j=1:nt, k=1:ns],(M[1,j,k]/(2.2599*ZMui)/4)*8000-1.79615*rspeed1[j,k]>=0)

@NLconstraint(gp,[j=1:nt, k=1:ns],3.6336*rspeed1[j,k]-(M[1,j,k]/(2.2599*ZMui)/4)*8000>=0)

@NLconstraint(gp,[j=1:nt, k=1:ns],-1423.7*(M[1,j,k]/(2.2599*ZMui)/4)^2+
4713*(M[1,j,k]/(2.2599*ZMui)/4)*rspeed1[j,k]/8000 + 6635.9*(rspeed1[j,k]/8000)^2-hed[1,j,k]==0)

@NLconstraint(gp,[j=1:nt, k=1:ns],P[1,j,k]-5*10^3*((1+(1.346-1)*9.8*ZMui*hed[1,j,k]/1.346/0.9082/Rm/293/1000)^(1.346/(1.346-1)))==0)

#pipe segment from the CS1 to the CS2
@NLconstraint(gp, [i = 1:16, j = 1:(nt-1), k=1:ns], (P[i,j+1,k] + P[i+1,j+1,k] - P[i,j,k] - P[i+1,j,k])*1000/2/dt +
(a^2/A)*(M[i+1,j,k] + M[i+1,j+1,k] - M[i,j,k] - M[i,j+1,k])/2/dx == 0)

@NLconstraint(gp, [i = 1:16, j = 1:(nt-1), k=1:ns], (P[i+1,j,k] + P[i+1,j+1,k] - P[i,j,k] - P[i,j+1,k])*1000/2/dx +
(1/A)*(M[i,j+1,k] + M[i+1,j+1,k] - M[i,j,k] - M[i+1,j,k])/2/dt +
(1/lam^2)*a*a*((M[i,j,k] + M[i,j+1,k] + M[i+1,j,k] + M[i+1,j+1,k])/4)*abs((M[i,j,k] + M[i,j+1,k] + M[i+1,j,k] + M[i+1,j+1,k])/4)/
2/d/A/A/((P[i,j,k] + P[i,j+1,k] + P[i+1,j,k] + P[i+1,j+1,k])*1000/4) == 0)

#CS2
@constraint(gp, [j = 2:nt, k=1:ns], M[17,j,k] - M[18,j,k] == 0)

@NLconstraint(gp,[j=1:nt, k=1:ns], (M[18,j,k]/(ZMui*den[1,j,k])/2)*ratespeed-minflow*rspeed2[j,k]>=0)

@NLconstraint(gp,[j=1:nt, k=1:ns], maxflow*rspeed2[j,k]-(M[18,j,k]/(ZMui*den[1,j,k])/2)*ratespeed>=0)

@NLconstraint(gp,[j=1:nt, k=1:ns], QVSH[1]*(M[18,j,k]/(ZMui*den[1,j,k])/2)^2+
QVSH[2]*(M[18,j,k]/(ZMui*den[1,j,k])/2)*rspeed2[j,k]/ratespeed + QVSH[3]*(rspeed2[j,k]/ratespeed)^2-hed[2,j,k]==0)

@NLconstraint(gp,[j=1:nt, k=1:ns], P[18,j,k]-P[17,j,k]*((1+(kv[1,j,k]-1)*9.8*ZMui*hed[2,j,k]/kv[1,j,k]/Z[1,j,k]/Rm/293/1000)^(kv[1,j,k]/(kv[1,j,k]-1)))==0)

#Nodes from the CS2 to the CS3
@NLconstraint(gp, [i = 18:33, j = 1:(nt-1), k=1:ns], (P[i,j+1,k] + P[i+1,j+1,k] - P[i,j,k] - P[i+1,j,k])*1000/2/dt +
(a^2/A)*(M[i+1,j,k] + M[i+1,j+1,k] - M[i,j,k] - M[i,j+1,k])/2/dx == 0)

@NLconstraint(gp, [i = 18:33, j = 1:(nt-1), k=1:ns], (P[i+1,j,k] + P[i+1,j+1,k] - P[i,j,k] - P[i,j+1,k])*1000/2/dx +
(1/A)*(M[i,j+1,k] + M[i+1,j+1,k] - M[i,j,k] - M[i+1,j,k])/2/dt +
(1/lam^2)*a*a*((M[i,j,k] + M[i,j+1,k] + M[i+1,j,k] + M[i+1,j+1,k])/4)*abs((M[i,j,k] + M[i,j+1,k] + M[i+1,j,k] + M[i+1,j+1,k])/4)/
2/d/A/A/((P[i,j,k] + P[i,j+1,k] + P[i+1,j,k] + P[i+1,j+1,k])*1000/4) == 0)

#CS3
@constraint(gp, [j = 2:nt, k=1:ns], M[34,j,k] - M[35,j,k] == 0)

@NLconstraint(gp,[j=1:nt, k=1:ns], (M[35,j,k]/(ZMui*den[2,j,k])/2)*ratespeed-minflow*rspeed3[j,k]>=0)

@NLconstraint(gp,[j=1:nt, k=1:ns], maxflow*rspeed3[j,k]-(M[35,j,k]/(ZMui*den[2,j,k])/2)*ratespeed>=0)

@NLconstraint(gp,[j=1:nt, k=1:ns], QVSH[1]*(M[35,j,k]/(ZMui*den[2,j,k])/2)^2+
QVSH[2]*(M[35,j,k]/(ZMui*den[2,j,k])/2)*rspeed3[j,k]/ratespeed + QVSH[3]*(rspeed3[j,k]/ratespeed)^2-hed[3,j,k]==0)

@NLconstraint(gp,[j=1:nt, k=1:ns], P[35,j,k]-P[34,j,k]*((1+(kv[2,j,k]-1)*9.8*ZMui*hed[3,j,k]/kv[2,j,k]/Z[2,j,k]/Rm/293/1000)^(kv[2,j,k]/(kv[2,j,k]-1)))==0)

#Nodes from the CS3 to the CS4
@NLconstraint(gp, [i = 35:50, j = 1:(nt-1), k=1:ns], (P[i,j+1,k] + P[i+1,j+1,k] - P[i,j,k] - P[i+1,j,k])*1000/2/dt +
(a^2/A)*(M[i+1,j,k] + M[i+1,j+1,k] - M[i,j,k] - M[i,j+1,k])/2/dx == 0)

@NLconstraint(gp, [i = 35:50, j = 1:(nt-1), k=1:ns], (P[i+1,j,k] + P[i+1,j+1,k] - P[i,j,k] - P[i,j+1,k])*1000/2/dx +
(1/A)*(M[i,j+1,k] + M[i+1,j+1,k] - M[i,j,k] - M[i+1,j,k])/2/dt +
(1/lam^2)*a*a*((M[i,j,k] + M[i,j+1,k] + M[i+1,j,k] + M[i+1,j+1,k])/4)*abs((M[i,j,k] + M[i,j+1,k] + M[i+1,j,k] + M[i+1,j+1,k])/4)/
2/d/A/A/((P[i,j,k] + P[i,j+1,k] + P[i+1,j,k] + P[i+1,j+1,k])*1000/4) == 0)

#CS4
@constraint(gp, [j = 2:nt, k=1:ns], M[51,j,k] - M[52,j,k] == D1)

@NLconstraint(gp,[j=1:nt, k=1:ns], (M[52,j,k]/(ZMui*den[3,j,k])/2)*ratespeed-minflow*rspeed4[j,k]>=0)

@NLconstraint(gp,[j=1:nt, k=1:ns], maxflow*rspeed4[j,k]-(M[52,j,k]/(ZMui*den[3,j,k])/2)*ratespeed>=0)

@NLconstraint(gp,[j=1:nt, k=1:ns], QVSH[1]*(M[52,j,k]/(ZMui*den[3,j,k])/2)^2+
QVSH[2]*(M[52,j,k]/(ZMui*den[3,j,k])/2)*rspeed4[j,k]/ratespeed + QVSH[3]*(rspeed4[j,k]/ratespeed)^2-hed[4,j,k]==0)

@NLconstraint(gp,[j=1:nt, k=1:ns], P[52,j,k]-P[51,j,k]*((1+(kv[3,j,k]-1)*9.8*ZMui*hed[4,j,k]/kv[3,j,k]/Z[3,j,k]/Rm/293/1000)^(kv[3,j,k]/(kv[3,j,k]-1)))==0)

#Nodes from the CS4 to the CS5
@NLconstraint(gp, [i = 52:67, j = 1:(nt-1), k=1:ns], (P[i,j+1,k] + P[i+1,j+1,k] - P[i,j,k] - P[i+1,j,k])*1000/2/dt +
(a^2/A)*(M[i+1,j,k] + M[i+1,j+1,k] - M[i,j,k] - M[i,j+1,k])/2/dx == 0)

@NLconstraint(gp, [i = 52:67, j = 1:(nt-1), k=1:ns], (P[i+1,j,k] + P[i+1,j+1,k] - P[i,j,k] - P[i,j+1,k])*1000/2/dx +
(1/A)*(M[i,j+1,k] + M[i+1,j+1,k] - M[i,j,k] - M[i+1,j,k])/2/dt +
(1/lam^2)*a*a*((M[i,j,k] + M[i,j+1,k] + M[i+1,j,k] + M[i+1,j+1,k])/4)*abs((M[i,j,k] + M[i,j+1,k] + M[i+1,j,k] + M[i+1,j+1,k])/4)/
2/d/A/A/((P[i,j,k] + P[i,j+1,k] + P[i+1,j,k] + P[i+1,j+1,k])*1000/4) == 0)

#CS5
@constraint(gp, [j = 2:nt, k=1:ns], M[68,j,k] - M[69,j,k] == D2)

@NLconstraint(gp,[j=1:nt, k=1:ns], (M[69,j,k]/(ZMui*den[4,j,k])/2)*ratespeed-minflow*rspeed5[j,k]>=0)

@NLconstraint(gp,[j=1:nt, k=1:ns], maxflow*rspeed5[j,k]-(M[69,j,k]/(ZMui*den[4,j,k])/2)*ratespeed>=0)

@NLconstraint(gp,[j=1:nt, k=1:ns], QVSH[1]*(M[69,j,k]/(ZMui*den[4,j,k])/2)^2+
QVSH[2]*(M[69,j,k]/(ZMui*den[4,j,k])/2)*rspeed5[j,k]/ratespeed + QVSH[3]*(rspeed5[j,k]/ratespeed)^2-hed[5,j,k]==0)

@NLconstraint(gp,[j=1:nt, k=1:ns], P[69,j,k]-P[68,j,k]*((1+(kv[4,j,k]-1)*9.8*ZMui*hed[5,j,k]/kv[4,j,k]/Z[4,j,k]/Rm/293/1000)^(kv[4,j,k]/(kv[4,j,k]-1)))==0)

#Nodes from the CS5 to the Customer 3
@NLconstraint(gp, [i = 69:78, j = 1:(nt-1), k=1:ns], (P[i,j+1,k] + P[i+1,j+1,k] - P[i,j,k] - P[i+1,j,k])*1000/2/dt +
(a^2/A)*(M[i+1,j,k] + M[i+1,j+1,k] - M[i,j,k] - M[i,j+1,k])/2/dx == 0)

@NLconstraint(gp, [i = 69:78, j = 1:(nt-1), k=1:ns], (P[i+1,j,k] + P[i+1,j+1,k] - P[i,j,k] - P[i,j+1,k])*1000/2/dx +
(1/A)*(M[i,j+1,k] + M[i+1,j+1,k] - M[i,j,k] - M[i+1,j,k])/2/dt +
(1/lam^2)*a*a*((M[i,j,k] + M[i,j+1,k] + M[i+1,j,k] + M[i+1,j+1,k])/4)*abs((M[i,j,k] + M[i,j+1,k] + M[i+1,j,k] + M[i+1,j+1,k])/4)/
2/d/A/A/((P[i,j,k] + P[i,j+1,k] + P[i+1,j,k] + P[i+1,j+1,k])*1000/4) == 0)

#Customer 3
@NLconstraint(gp, [j = 1:(nt-1), k=1:ns], (P[79,j+1,k] + P[80,j+1,k] - P[79,j,k] - P[80,j,k])*1000/2/dt +
(a^2/A)*(M[80,j,k] + M[80,j+1,k] - (M[79,j,k]-D3) - (M[79,j+1,k]-D3))/2/dx == 0)

@NLconstraint(gp, [j = 1:(nt-1), k=1:ns], (P[80,j,k] + P[80,j+1,k] - P[79,j,k] - P[79,j+1,k])*1000/2/dx +
(1/A)*((M[79,j+1,k]-D3) + M[80,j+1,k] - (M[79,j,k]-D3) - M[80,j,k])/2/dt +
(1/lam^2)*a*a*(((M[79,j,k]-D3) + (M[79,j+1,k]-D3) + M[80,j,k] + M[80,j+1,k])/4)*abs(((M[79,j,k]-D3) + (M[79,j+1,k]-D3) + M[80,j,k] + M[80,j+1,k])/4)/
2/d/A/A/((P[79,j,k] + P[79,j+1,k] + P[80,j,k] + P[80,j+1,k])*1000/4) == 0)

#Nodes from the Customer 3 to the underground gas storage
@NLconstraint(gp, [i = 80:86, j = 1:(nt-1), k=1:ns], (P[i,j+1,k] + P[i+1,j+1,k] - P[i,j,k] - P[i+1,j,k])*1000/2/dt +
(a^2/A)*(M[i+1,j,k] + M[i+1,j+1,k] - M[i,j,k] - M[i,j+1,k])/2/dx == 0)

@NLconstraint(gp, [i = 80:86, j = 1:(nt-1), k=1:ns], (P[i+1,j,k] + P[i+1,j+1,k] - P[i,j,k] - P[i,j+1,k])*1000/2/dx +
(1/A)*(M[i,j+1,k] + M[i+1,j+1,k] - M[i,j,k] - M[i+1,j,k])/2/dt +
(1/lam^2)*a*a*((M[i,j,k] + M[i,j+1,k] + M[i+1,j,k] + M[i+1,j+1,k])/4)*abs((M[i,j,k] + M[i,j+1,k] + M[i+1,j,k] + M[i+1,j+1,k])/4)/
2/d/A/A/((P[i,j,k] + P[i,j+1,k] + P[i+1,j,k] + P[i+1,j+1,k])*1000/4) == 0)

#Underground Gas storage
@NLconstraint(gp, [j = 1:(nt-1), k=1:ns], (P[87,j+1,k] + P[88,j+1,k] - P[87,j,k] - P[88,j,k])*1000/2/dt +
(a^2/A)*(M[88,j,k] + M[88,j+1,k] - (M[87,j,k]+ UGS[j,k]) - (M[87,j+1,k]+ UGS[j+1,k]))/2/dx == 0)

@NLconstraint(gp, [j = 1:(nt-1), k=1:ns], (P[88,j,k] + P[88,j+1,k] - P[87,j,k] - P[87,j+1,k])*1000/2/dx +
(1/A)*((M[87,j+1,k]+ UGS[j+1,k]) + M[88,j+1,k] - (M[87,j,k]+ UGS[j,k]) - M[88,j,k])/2/dt +
(1/lam^2)*a*a*((M[87,j,k]+ UGS[j,k] + M[87,j+1,k]+ UGS[j+1,k] + M[88,j,k] + M[88,j+1,k])/4)*abs((M[87,j,k]+ UGS[j,k] + M[87,j+1,k]+ UGS[j+1,k] + M[88,j,k] + M[88,j+1,k])/4)/
2/d/A/A/((P[87,j,k] + P[87,j+1,k] + P[88,j,k] + P[88,j+1,k])*1000/4) == 0)

#Nodes from the the underground gas storage to customer 4
@NLconstraint(gp, [i = 88:94, j = 1:(nt-1), k=1:ns], (P[i,j+1,k] + P[i+1,j+1,k] - P[i,j,k] - P[i+1,j,k])*1000/2/dt +
(a^2/A)*(M[i+1,j,k] + M[i+1,j+1,k] - M[i,j,k] - M[i,j+1,k])/2/dx == 0)

@NLconstraint(gp, [i = 88:94, j = 1:(nt-1), k=1:ns], (P[i+1,j,k] + P[i+1,j+1,k] - P[i,j,k] - P[i,j+1,k])*1000/2/dx +
(1/A)*(M[i,j+1,k] + M[i+1,j+1,k] - M[i,j,k] - M[i+1,j,k])/2/dt +
(1/lam^2)*a*a*((M[i,j,k] + M[i,j+1,k] + M[i+1,j,k] + M[i+1,j+1,k])/4)*abs((M[i,j,k] + M[i,j+1,k] + M[i+1,j,k] + M[i+1,j+1,k])/4)/
2/d/A/A/((P[i,j,k] + P[i,j+1,k] + P[i+1,j,k] + P[i+1,j+1,k])*1000/4) == 0)

#Boundary condition
@NLconstraint(gp, [j = 2:nt, k = 1:ns], M[95,j,k] - D4[j,k] == 0)

#First stage constraints
@constraint(gp, [j = 2:nt1, k = 1:ns], P1[1, j] - P[1,j,k] == 0)
@constraint(gp, [j = 2:nt1, k = 1:ns], P1[2, j] - P[18,j,k] == 0)
@constraint(gp, [j = 2:nt1, k = 1:ns], P1[3, j] - P[35,j,k] == 0)
@constraint(gp, [j = 2:nt1, k = 1:ns], P1[4, j] - P[52,j,k] == 0)
@constraint(gp, [j = 2:nt1, k = 1:ns], P1[5, j] - P[69,j,k] == 0)

@constraint(gp, [j = 2:nt1, k = 1:ns], M1[j] - M[1,j,k] == 0)

@constraint(gp, [j = 2:nt1, k = 1:ns], UGS1[j] - UGS[j,k] == 0)

@NLconstraint(gp, [j = 2:nt-1, k = 1:ns], UGS[j,k] - UGS[j-1,k] <= 5)

@NLconstraint(gp, [j = 2:nt-1, k = 1:ns], UGS[j,k] - UGS[j-1,k] >= -5)

optimize!(gp)

println("gas withdrawn cost:")
println(sum(e2*sum(value(UGS[j,k])*30*60/(0.0417*16.3776) for j in 1:nt-1) for k in 1:ns)/ns)

println("line pack:")
println(sum(sum(value((P[i,nt,k]+P[i+1,nt,k])/2) for i in 1:16)*(10^7)*A/a/a +
sum(value((P[i,nt,k]+P[i+1,nt,k])/2) for i in 18:33)*(10^7)*A/a/a + sum(value((P[i,nt,k]+P[i+1,nt,k])/2) for i in 35:50)*(10^7)*A/a/a +
sum(value((P[i,nt,k]+P[i+1,nt,k])/2) for i in 52:67)*(10^7)*A/a/a + sum(value((P[i,nt,k]+P[i+1,nt,k])/2) for i in 69:94)*(10^7)*A/a/a for k in 1:ns)/ns/(0.0417*16.3776))

println("electricity cost:")
println(sum((sum(4*100*value(hed[1,j,k])*9.8*(value(M[1,j,k])/4)/1000/value(xiaolv[1,j,k]) for j in 1:nt-1)
+ sum(2*100*value(hed[2,j,k])*9.8*(value(M[18,j,k])/2)/1000/value(xiaolv[2,j,k]) for j in 1:nt-1)
+ sum(2*100*value(hed[3,j,k])*9.8*(value(M[35,j,k])/2)/1000/value(xiaolv[3,j,k]) for j in 1:nt-1)
+ sum(2*100*value(hed[4,j,k])*9.8*(value(M[52,j,k])/2)/1000/value(xiaolv[4,j,k]) for j in 1:nt-1)
+ sum(2*100*value(hed[5,j,k])*9.8*(value(M[69,j,k])/2)/1000/value(xiaolv[5,j,k]) for j in 1:nt-1))*0.5*e1  for k in 1:ns)/ns)
