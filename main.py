import numpy as np
from numpy import sqrt as sqrt
from numpy import linspace
from numpy import array 
from Classes.EngineStates import *
from Equations.Equations import *
import matplotlib.pyplot as plt

# Initial Values for Static Values
h = 11  # Km
g = 9.81  # m/s/s
M0 = 0.85  #
BPR = 10  #
Yc = 1.4  #
Yh = 1.333  #
R = 287
FPR = 1.5  #
CPR = 36  #
S = 285.0  # m^2
CombustorPressureLoss = 0.96  #
IntakeEfficiency = 0.98  #
MechEfficiency = 0.99  #
PolytropicFan = 0.89  #
PolytropicCompressor = 0.9  #
PolytropicTurbine = 0.9  #
CombustorEfficiency = 0.99  #
NozzleEfficiency = 0.99  #
HFuel = 431000.0  # Kj/kg
Cpc = 1005
Cph = 1148

PAmbient = 22696.8  # Pa
TAmbient = 216.775245  # K
DensityAmbient = 0.364748106  # Kg/ M^3

# Converting Imperial to Metric
WTO = 370000.0  # lbf
WTO = WTO * 4.448
# WTO = 1645760 N
Weight = WTO * 0.8  # 296000 lbf, 1316608 N

# Initiaize all engine States all will be after the component mentioned
Ambient = EngineStates()
Intake = EngineStates()
Fan = EngineStates()
FanNozzle = EngineStates()
Compressor = EngineStates()
Combustor = EngineStates()
Turbine = EngineStates()
TurbineNozzle = EngineStates()

EngineState = {
    "Ambient": Ambient, "Intake": Intake, "Fan": Fan, "FanNozzle": FanNozzle, "Compressor": Compressor,
    "Combustor": Combustor,
    "Turbine": Turbine, "TurbineNozzle": TurbineNozzle
}
EngineState["Ambient"].Stage = "Ambient"
EngineState["Intake"].Stage = "Intake"
EngineState["Fan"].Stage = "Fan"
EngineState["FanNozzle"].Stage = "FanNozzle"
EngineState["Compressor"].Stage = "Compressor"
EngineState["Combustor"].Stage = "Combustor"
EngineState["Turbine"].Stage = "Turbine"
EngineState["TurbineNozzle"].Stage = "TurbineNozzle"

q = DynamicPressure(Yc, PAmbient, M0)
Cl = LiftCoefficient(Weight, S, q)
CD = DragCoefficient(Cl)
D = Drag(CD, S, q)
A0 = SpeedofSound(R, Yc, TAmbient)
V0 = A0 * M0
PFC = PolytropicTurbineCoefficent(PolytropicFan, Yc)
PCC = PolytropicTurbineCoefficent(PolytropicCompressor, Yc)

T02 = inlet.Temp(V0, Cpc, TAmbient)
P02 = inlet.Pressure(Yc, V0, TAmbient, PAmbient, Cpc, IntakeEfficiency)

T02_1 = fan.temperature(PolytropicFan, T02, FPR, Yc)
P02_1 = fan.pressure(P02, FPR)

T03 = compressor.temperature(PolytropicCompressor, T02_1, CPR, Yc)
P03 = compressor.pressure(P02_1, CPR)

T04 = 1560
P04 = P03 * CombustorPressureLoss

T05 = turbine.temperature(BPR, Cpc, Cph, MechEfficiency, T02, T02, T02_1, T03, T04)
# T05 = turbine.temperatureWork(Cph,Cpc,T02,T03,T04)
P05 = turbine.pressure(T04, T05, P04, Yh, PolytropicTurbine)

P9 = PAmbient
T9 = nozzle.temperature(NozzleEfficiency, T05, PAmbient, P05, Yh)
M9 = nozzle.Mach(T05, T9, Yh)
V9 = M9 * np.sqrt(Yh * R * T9)
print(P05 / P9)
print("Mach is:", M9)
print("Velocity is:", V9)


P19 = PAmbient
T19 = nozzle.temperature(NozzleEfficiency, T02_1, PAmbient, P02_1, Yc)
M19 = nozzle.Mach(T02_1, T19, Yc)
V19 = M19 * np.sqrt(Yc * R * T19)

massTotal, mdot_c, mdot_f = combustor.massflow(BPR, V0, V9, V19, D)
f0, MassFlowFuel = combustor.FAR(Cph, Cpc, T03, T04, HFuel, CombustorEfficiency, mdot_c)



F_m00 = F_m0(D, mdot_f, mdot_c)
TSFC0 = TSFC(D,massTotal,MassFlowFuel)
fa_ratio0, mdot_h0 = fa_ratio(T04, T03, CombustorEfficiency, Cpc, HFuel, mdot_f, mdot_c)
##
nt0, np0, no0 = perf_eff(MassFlowFuel, mdot_c, mdot_f, V9, V19, V0, HFuel)
#print(mdot_c,mdot_f,mdot_h0,MassFlowFuel)
print(F_m00,nt0,np0,no0)

###################
#Initial values for Graphs
BPR_values = linspace(5, 20, num=100,dtype=None) #20 values between 5-20 BPR
FanPR_values = linspace(1.2, 2.0, num=100,dtype=float) #20 values between 1.2-2.0 Fan pressure ratio
CompPR_values = linspace(20, 40, num=100,dtype=float) #20 values between 20-40 compressor pressure ratio

#Data Lists for Graphs
#Data when BPR changes 
BPR_fs_data = []
BPR_TSFC_data = []
BPR_f_data = []
BPR_eta_t_data = []
BPR_eta_p_data = []
BPR_eta_o_data = []

#Data when Fan Pressure Ratio Changes
Fan_fs_data = []
Fan_TSFC_data = []
Fan_f_data = []
Fan_eta_t_data = []
Fan_eta_p_data = []
Fan_eta_o_data = []

#data lists for compressor pressure ratios
Comp_fs_data = []
Comp_TSFC_data = []
Comp_f_data = []
Comp_eta_t_data = []
Comp_eta_p_data = []
Comp_eta_o_data = []

#Class that runs Calcs
def solve_engine(BPR, fan_PR, comp_PR):

    #obtaining the values given from the for loops
    BPR = BPR #bypass pressure ratio
    FPR = fan_PR #pressure ratio of fan
    CPR = comp_PR #compressor ratio
    T02 = inlet.Temp(V0, Cpc, TAmbient)
    P02 = inlet.Pressure(Yc, V0, TAmbient, PAmbient, Cpc, IntakeEfficiency)

    T02_1 = fan.temperature(PolytropicFan, T02, FPR, Yc)
    P02_1 = fan.pressure(P02, FPR)

    T03 = compressor.temperature(PolytropicCompressor, T02_1, CPR, Yc)
    P03 = compressor.pressure(P02_1, CPR)

    T04 = 1560
    P04 = P03 * CombustorPressureLoss

    T05 = turbine.temperature(BPR, Cpc, Cph, MechEfficiency, T02, T02, T02_1, T03, T04)
    # T05 = turbine.temperatureWork(Cph,Cpc,T02,T03,T04)
    P05 = turbine.pressure(T04, T05, P04, Yh, PolytropicTurbine)

    P9 = PAmbient
    T9 = nozzle.temperature(NozzleEfficiency, T05, PAmbient, P05, Yh)
    M9 = nozzle.Mach(T05, T9, Yh)
    V9 = M9 * sqrt(Yh * R * T9)

    #mdot_f = 0
    #mdot_c = 0
    P19=PAmbient
    T19 = nozzle.temperature(NozzleEfficiency, T02_1, PAmbient, P02_1, Yh)
    M19 = nozzle.Mach(T02_1, T19, Yh)
    V19 = M19 * sqrt(Yh * R * T19)

    massTotal, mdot_c, mdot_f = combustor.massflow(BPR, V0, V9, V19, D)
    f1, MassFlowFuel = combustor.FAR(Cph, Cpc, T03, T04, HFuel, CombustorEfficiency, mdot_c)
    ##print(D,mdot_f,mdot_c)
    F_m01 = F_m0(D,mdot_f,mdot_c)
    TSFC1 = TSFC(D,massTotal,MassFlowFuel)
    fa_ratio1,mdot_h = fa_ratio(T04,T03,CombustorEfficiency,Cpc,HFuel,mdot_f,mdot_c)
    n_t1,n_p1,n_o1 = perf_eff(MassFlowFuel,mdot_c,mdot_f,V9,V19,V0,HFuel)
    #f=0

    x = array([F_m01, TSFC1, f1, n_t1, n_p1, n_o1]) #array with holding the calculated parameters

    return x


#Loop that Calculates values when BPR changes
for i in BPR_values:
    x = solve_engine(i, 1.5, 36)

    BPR_fs_data.append(x[0])
    BPR_TSFC_data.append(x[1])
    BPR_f_data.append(x[2])
    BPR_eta_t_data.append(x[3])
    BPR_eta_p_data.append(x[4])
    BPR_eta_o_data.append(x[5])

#Loop to calculate values when fan Pressure ratio changes

for i in FanPR_values:
    x = solve_engine(10,i,36)

    Fan_fs_data.append(x[0])
    Fan_TSFC_data.append(x[1])
    Fan_f_data.append(x[2])
    Fan_eta_t_data.append(x[3])
    Fan_eta_p_data.append(x[4])
    Fan_eta_o_data.append(x[5])

#Loop to calculate values when Compressor pressure ratio changes

for i in CompPR_values:
    x =  solve_engine(10,1.5,i)

    Comp_fs_data.append(x[0])
    Comp_TSFC_data.append(x[1])
    Comp_f_data.append(x[2])
    Comp_eta_t_data.append(x[3])
    Comp_eta_p_data.append(x[4])
    Comp_eta_o_data.append(x[5])


#Graph for BPR

fig1, ax1 = plt.subplots(nrows = 2, ncols=3)

ax1[0,0].plot(BPR_values, BPR_fs_data)
ax1[0,1].plot(BPR_values, BPR_f_data)
ax1[0,2].plot(BPR_values,BPR_TSFC_data)

ax1[1,0].plot(BPR_values,BPR_eta_t_data)
ax1[1,1].plot(BPR_values,BPR_eta_p_data)
ax1[1,2].plot(BPR_values,BPR_eta_o_data)



#Graph for Fan Pressure Ratio
fig2, ax2 = plt.subplots(nrows = 2, ncols=3)

ax2[0,0].plot(FanPR_values, Fan_fs_data)
ax2[0,1].plot(FanPR_values, Fan_f_data)
ax2[0,2].plot(FanPR_values,Fan_TSFC_data)

ax2[1,0].plot(FanPR_values,Fan_eta_t_data)
ax2[1,1].plot(FanPR_values,Fan_eta_p_data)
ax2[1,2].plot(FanPR_values,Fan_eta_o_data)



#Graph for Compressor Pressure Ratio
fig3,ax3 = plt.subplots(nrows = 2, ncols = 3)

ax3[0,0].plot(CompPR_values, Comp_fs_data)
ax3[0,1].plot(CompPR_values, Comp_f_data)
ax3[0,2].plot(CompPR_values,Comp_TSFC_data)

ax3[1,0].plot(CompPR_values,Comp_eta_t_data)
ax3[1,1].plot(CompPR_values,Comp_eta_p_data)
ax3[1,2].plot(CompPR_values,Comp_eta_o_data)



plt.show()