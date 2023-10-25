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
        after = P02 * compressorpressureratio
        return after

    @staticmethod
    def temperature(polytropiccompressor, T02, compressorpressureratio, y):
        temp = (compressorpressureratio) ** ((y - 1) / (y * polytropiccompressor))
        after = T02 * temp
        return after

    @staticmethod
    def Work(temperatureafter, temperaturebefore, cp):
        workdone = cp * (temperatureafter - temperaturebefore)
        return workdone
