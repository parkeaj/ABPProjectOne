import numpy as np
from Classes.EngineStates import *
from Equations.Equations import *

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

mdot_f = 0
mdot_c = 0
P19 = PAmbient
T19 = nozzle.temperature(NozzleEfficiency, T02_1, PAmbient, P02_1, Yh)
M19 = nozzle.Mach(T02_1, T19, Yh)
V19 = M19 * np.sqrt(Yh * R * T19)







F_m0 = F_m0(D, mdot_f, mdot_c)
TSFC = TSFC(D, mdot_f)
fa_ratio, mdot_h = fa_ratio(T04, T03, CombustorEfficiency, Cpc, HFuel, mdot_f, mdot_c)
nt, np, no = perf_eff(mdot_h, mdot_c, mdot_f, V9, V19, V0, HFuel)
