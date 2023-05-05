#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 17 13:54:28 2021
@author: sylvain
"""

import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv
import time

import os
os.chdir(os.path.dirname(os.path.realpath(__file__)))


if __name__ == "__main__":

    # data = np.loadtxt("log-benchmark-pi.txt", skiprows=1, delimiter="\t")

    data = read_csv("log-benchmark-pi.txt", delimiter="\t")

    for i in range(data.shape[0]):
        x = data["time"][i].split(':')
        data["time"][i] = float(x[1])




    for k in range(1,5):
        p =data.loc[data["nb-process"] == k]
        plt.scatter(p["nb-iter"], p["time"], label=f"nb-processes = {k}")

    plt.grid()

    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.xlabel("N nombre d'itérations")
    plt.ylabel("temps execution (s)")


    plt.legend()

    plt.show()
