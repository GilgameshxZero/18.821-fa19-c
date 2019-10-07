#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 12:12:13 2019

@author: tamiquedebrito
"""

import numpy as np

import matplotlib.pyplot as plt
"""
Example:
"""
graph_perts(Poly([6j - 1, 3, 1 + 5j, 4, -6, 1j]), 100, .5)


class Poly:
    # Coefficients an array with highest degree first
    def __init__(self, coeff):
        self.N = len(coeff) - 1
        self.coeff = coeff

    def evaluate(self, x):
        return sum(
            [self.coeff[i] * x**(self.N - i) for i in range(self.N + 1)])

    def __str__(self):
        y = str(self.coeff[0]) + "z^" + str(self.N)
        for i in range(1, len(self.coeff)):
            y += "+" + str(self.coeff[i]) + "z^" + str(self.N - i)
        return y


def perturb_coeff(coeff, epsilon, which="all"):
    # Perturbs a subset of the elements of the list "coeff" by a normal distribution with mean 0 and variance "epsilon".
    # The subset of elements which are perturbed are determined by "which", where which="all" means that all are perturbed,
    #     and otherwise, which is a list of booleans so that whether index i in "which" is True/False determines whether coeff[i]
    #     is perturbed.
    pert = coeff[:]
    if which == "all":
        for i in range(1, len(coeff)):
            p = np.random.normal(scale=epsilon, size=2)
            pert[i] = coeff[i] + p[0] + p[1] * 1j
        return pert
    else:
        assert type(which) == list
        for i in range(1, len(coeff)):
            pert[i] = coeff[i]
            if which[i]:
                p = np.random.normal(scale=epsilon, size=2)
                pert[i] = coeff[i] + p[0] + p[1] * 1j
        return pert


def complex_list_to_points(L):
    # Converts a list of complex numbers to two lists of corresponding real/imaginary parts.
    # This is for the purpose of being able to graph the complex numbers.
    return [z.real for z in L], [z.imag for z in L]


def pert_roots(P, epsilon, which="all"):
    # Gets the resulting roots from randomly perturbing a given set of coefficients according to "epsilon" and "which".
    return np.roots(perturb_coeff(P.coeff, epsilon, which))


def N_pert_roots(P, N, epsilon, which):
    # Returns a list of "N" randomly perturbed root-sets
    rootList = []
    for _ in range(N):
        rootList.append(pert_roots(P, epsilon, which))
    return rootList


def allpairs(rootList):
    col = []
    for l in rootList:
        col.append(l)
    return complex_list_to_points(col)


def graph_perts(P, N, epsilon, which="all"):
    # Graph the results of N perturbations
    x, y = allpairs(N_pert_roots(P, N, epsilon, which))
    plt.plot(x, y, 'ko')


def graph_perts_single_root(P, N, i, epsilon, which="all"):
    # Graph the results of N perturbations
    single_roots = [r[i] for r in N_pert_roots(P, N, epsilon, which)]
    x, y = complex_list_to_points(single_roots)
    plt.plot(x, y, 'ko')
    return x, y
