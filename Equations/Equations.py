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
        T02 = T0 + velocity ** 2 / (2 * Cp)
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


class turbine:

    @staticmethod
    def temperature(BPR, Cpc, Cph, nm, T01, T02, T02_1, T03, T04):
        temp = (BPR * Cpc * (T02_1 - T02) - Cpc * (T02_1 - T03) - Cph * nm * T04 + Cpc * (T02_1 - T02)) * -1
        T05 = temp / Cph / nm
        return T05

    """def temperature(WorkCompressor, WorkFan, Cph, T04):
        T05 = (WorkCompressor + WorkFan) * - 1 / Cph + T04
        return T05"""

    @staticmethod
    def pressure(T04, T05, P04, y, ntinf):
        PTC = PolytropicTurbineCoefficent(ntinf, y)
        temp = T04 / T05
        P05 = temp ** (1 / PTC) * P04
        return P05


class nozzle:
    @staticmethod
    def temperature():
        T06 = 1
        return T06

    @staticmethod
    def pressure():
        P06 = 1
        return P06
