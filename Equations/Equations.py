import sympy as s
import numpy as np


def DragCoefficient(Cl):
    Cd = 0.056 * Cl ** 2 - 0.004 * Cl + 0.0140
    return Cd


def Drag(DragCoefficient, WingArea, DynamicPressure):
    Drag = DragCoefficient * DynamicPressure * WingArea
    return Drag


def LiftCoefficient(lift, dynamicpressure, wingarea):
    Cl = lift / (dynamicpressure * wingarea)
    return Cl


def DynamicPressure(y, pressure, mach):
    DynamicPressure = 0.5 * mach * mach * pressure * y
    return DynamicPressure


def SpeedofSound(R, y, temperature):
    a = np.sqrt(R * temperature * y)
    return a


def PolytropicCompressCoefficent(PolytropicE, y):
    Coefficent = PolytropicE * (y - 1) / y
    return Coefficent


def PolytropicTurbineCoefficent(PolytropicE, y):
    Coefficent = (y - 1) / (y * PolytropicE)
    return Coefficent


def F_m0(drag, mdot_f, mdot_c):
    F_m0 = drag / (mdot_c + mdot_f)
    return F_m0


def TSFC(drag, mdot_f):
    TSFC = mdot_f / drag
    return TSFC


def fa_ratio(t04, t03, nb, cp, Q, mdot_f, mdot_c):
    f = ((t04 / t03) - 1) / ((nb * Q / cp * t03) - (t04 / t03))
    mdot_h = f * (mdot_f + mdot_c)
    return f, mdot_h


def perf_eff(mdot_h, mdot_c, mdot_f, c9, c19, Ca, Qf):
    mdot_a = mdot_f + mdot_c
    nt = (0.5 * (mdot_c * c9 ** 2 + mdot_f * c19 ** 2 - mdot_a * Ca ** 2)) / (mdot_h * Qf)
    np = (Ca * (mdot_f * (c19 - Ca) + mdot_c * (c9 - Ca))) / (
            0.5 * (mdot_c * c9 ** 2 + mdot_f * c19 ** 2 - mdot_a * Ca ** 2))
    no = nt * np
    return nt, np, no


"""
def PolytropicCompressor(PR, Coefficent, P01, T01):
    P02 = PR * P01
    T02 = T01 * (PR ** Coefficent)
    DeltaT = T02 - T01
    return P02, T02, DeltaT
def PolytropicTurbine(Cpc,Cph,nm,DeltaT, Coefficent, P01, T01):
    P02 = PR * P01
    T02 = T01 * (PR ** Coefficent)
    DeltaT = T02 - T01
    return P02, T02, DeltaT

"""


class inlet:
    @staticmethod
    def Temp(velocity, Cp, T0):
        T02 = T0 + 0.98 * (velocity ** 2 / (2 * Cp))
        return T02

    @staticmethod
    def Pressure(y, velocity, T0, P0, Cp, IntakeEfficiency):
        temp1 = IntakeEfficiency * velocity ** 2 / (2 * Cp * T0)
        temp2 = (1 + temp1) ** (y / (y - 1))
        P02 = temp2 * P0
        return P02


class fan:
    @staticmethod
    def pressure(P02, fanpressureratio):
        after = P02 * fanpressureratio
        return after

    @staticmethod
    def temperature(polytropicfan, T02, fanpressureratio, y):
        temp = (fanpressureratio) ** ((y - 1) / (y * polytropicfan))
        after = T02 * temp
        return after

    @staticmethod
    def Work(temperatureafter, temperaturebefore, cp):
        workdone = cp * (temperatureafter - temperaturebefore)
        return workdone


class compressor:
    @staticmethod
    def pressure(P02, compressorpressureratio):
        P03 = P02 * compressorpressureratio
        return P03

    @staticmethod
    def temperature(polytropiccompressor, T02, compressorpressureratio, y):
        temp = (compressorpressureratio) ** ((y - 1) / (y * polytropiccompressor))
        T03 = T02 * temp
        return T03

    @staticmethod
    def Work(temperatureafter, temperaturebefore, cp):
        workdone = cp * (temperatureafter - temperaturebefore)
        return workdone


class combustor:
    @staticmethod
    def pressure(PR, P03):
        P04 = PR * P03
        return P04

    @staticmethod
    def massflow(BPR, Ca, Cc, Ch, Drag):
        Temp1 = - (BPR + 1) * Drag
        Temp2 = BPR * (Ca - Cc) + Ca - Ch
        m = Temp1 / Temp2
        MassFlowCore = m / (BPR + 1)
        MassFlowNonCore = m * BPR / (BPR + 1)
        return m, MassFlowCore, MassFlowNonCore
    @staticmethod
    def FAR(Cph,Cpc,T04,T03,HFuel,Nb,MassFlowCore):
        Temp1 = Cph * T04 - Cpc * T03
        Temp2 = Nb * (HFuel - Cph * T04)

        f = Temp1/Temp2
        MassFlowFuel = MassFlowCore * f

        return f, MassFlowFuel

class turbine:

    @staticmethod
    def temperature(BPR, Cpc, Cph, nm, T01, T02, T02_1, T03, T04):
        # temp = (BPR * Cpc * (T02_1 - T02) -   Cpc * (T02_1 - T03) + Cph * nm * T04)
        temp = -1 * (BPR * (T02 - T02_1) + T02 - T03) * Cpc
        T05 = temp / (Cph * nm)
        T05 = T05 - T04
        T05 = T05 * -1
        return T05

    @staticmethod
    def temperatureWork(Cph, Cpc, T02, T03, T04):
        Temp1 = Cpc * (T03 - T02)
        Temp2 = -1 * Temp1 / Cph
        T05 = Temp2 + T04
        return T05

    @staticmethod
    def pressure(T04, T05, P04, y, ntinf):
        PTC = PolytropicTurbineCoefficent(ntinf, y)
        temp = T05 / T04
        P05 = temp ** (PTC) * P04
        return P05


class nozzle:
    @staticmethod
    def temperature(Ni, T05, P9, P06, y):
        T9 = (Ni * T05 * (1 - (P9 / P06) ** ((y - 1) / y)))
        return T9

    @staticmethod
    def Mach(T05, T9, y):
        temp1 = 2 * (T05 - T9)
        temp2 = T9 * (y - 1)
        M = np.sqrt(temp1 / temp2)
        return M
