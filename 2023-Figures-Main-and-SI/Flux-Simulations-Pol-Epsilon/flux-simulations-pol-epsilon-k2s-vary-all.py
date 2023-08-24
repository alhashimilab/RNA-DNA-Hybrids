# -*- coding: UTF-8 -*-
import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.optimize import curve_fit
from scipy import stats
from scipy.integrate import odeint


def dxdt(x, T, *args):
    # define differential equations system

    # A = GS
    # B = ES1
    # C = ES2
    # D = Mis

    #k1 = GS to ES1
    #k2 = GS to ES2
    #k3 = ES1 to ES2
    #k4 = ES1/ES2 to Mis
    
    A, B, C, D= x
    k1, k_1, k2, k_2, k3, k_3, k4, k_4 = args

    dAdt = -k1*A + k_1*B - k2*A + k_2*C
    dBdt = k1*A - k_1*B - k3*B - k4*B + k_3*C + k_4*D
    dCdt = k2*A - k_2*C + k3*B - k_3*C - k4*C + k_4*D
    dDdt = k4*B + k4*C - k_4*D - k_4*D

    return dAdt, dBdt, dCdt, dDdt

def sim(params=(5.74, 6127.67, 10.82, 3642, 2437, 769, 268, 100)):
    # This functions solves the kinetic differential equations at equilibrium
    # The parameters are taken from the input csv file, input params above are examples
    time = np.linspace(0, 10000000, 2)  # 'Solve' the diff eqs for equilibirum by picking a time point far in teh future, this is way faster in bulk then fully solving for d/dt=0'
    k1, k_1, k2, k_2, k3, k_3, k4, k_4 = params
    C0 = np.array([1, 0, 0, 0])
    Data = odeint(dxdt, C0, time, args = params)[-1]
    GS, ES1, ES2, Mis = Data

    # calculate flux through ES1 and through ES2
    Ftt = ((1/(GS*k1)) + (1/(ES1*k4)))**-1
    Fta = ((1/(GS*k2)) + (1/(ES2*k_3)) + (1/(ES1*k4)))**-1
    Faa = ((1/(GS*k2)) + (1/(ES2*k4)))**-1
    Fat = ((1/(GS*k1)) + (1/(ES1*k3)) + (1/(ES2*k4)))**-1
    FT = Ftt + Fta
    FA = Faa + Fat
    # fractional flux
    fT = FT / (FT + FA)
    fA = FA / (FT + FA)

    return [GS, ES1, ES2, Mis, FT, FA, fT, fA]


def read_table(csv_path):
    csv_table = []
    relevant_titles = ['kGStES1_pH', 'kES1tGS_pH', 'kGStES2_pH', 'kES2tGS_pH', 'kES1tES2_pH', 'kES2tES1_pH']

    with open(csv_path) as csv:
        lines = csv.readlines()

    titles = lines[0].split(',')
    find_titles = []
    for x in titles:
       find_titles.append(x in relevant_titles)

    for line in lines[1:]:
        csv_table += [np.array(line.strip().split(','))[find_titles].astype(float)]
    return csv_table

def create_vector(length, V_0, multip):
    k = [V_0]
    x = V_0
    for i in range(1, length):
        x = x*multip
        k.append(x)
    return np.asarray(k)

if __name__ == "__main__":
    parent_path = '/Users/orsula1/OrsDocs/OrsDocs_since_Jan2023/2023-RNA-DNA-Hybrid/2023-Manuscript-Figures/2023-For-Data-Deposition/2023-Figures-Main-and-SI/'

    params_file = 'params_for_python_flux/params_for_flux_calc_all.csv'
    path = os.path.join(parent_path,params_file)
    params_table = read_table(path)

    # create variable k2
    k2s = create_vector(115, 0.001, 1.235)

    out_folder = 'Flux-Simulations-Pol-Epsilon/k2s_all/'
    out_path = os.path.join(parent_path,out_folder)
    filenames = ['CGC_DNA_8p4','CGC_RNA_8p4','dGrU_pH_8p0','dTrG_pH_8p0','CGC_DNA_7p4','CGC_RNA_7p4','dGrU_pH_7p4','dTrG_pH_7p4','CGC_DNA_6p9','CGC_RNA_6p9','dGrU_pH_6p9','dTrG_pH_6p9']

    for params,filename in zip(params_table,filenames):
        out_file = out_path + filename + '.csv'
        with open(out_file, 'w') as out:
            l = 'GS,ES1,ES2,Mis,FT,FA,fT,fA\n'
            out.writelines(l)
            for k2 in k2s:
                comp_params = np.concatenate((params,np.array([k2,100])))
                sim_results = sim(params=tuple(comp_params))
                l = '{},{},{},{},{},{},{},{}\n'.format(*sim_results)
                out.writelines(l)