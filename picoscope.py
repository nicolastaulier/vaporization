#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import numpy as np

###############################################################################
def read(filename):
    """Reads a CSV file with a specified delimiter (;) and processes the data
    to extract time (in microseconds) and amplitude values (in mV).
    The first three rows are skipped, and the remaining rows are used to populate
    the time and amplitude arrays, which are converted in seconds and volts.

    :param filename: The name of the CSV file to read.
    :type filename: str
    :return: A tuple containing two numpy arrays: time values (t) and amplitude values (A).
    :rtype: tuple of numpy.ndarray
    """
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile,delimiter=';')
        data = [[x.strip() for x in row] for row in reader]
        M = 3
        N = len(data) - M
        t = np.zeros(N)
        A = np.zeros(N)
        for i in range(0,N):
            t[i] = float(data[i+M][0].replace(',','.'))*1e-6
            A[i] = float(data[i+M][1].replace(',','.'))*1e-3
        return (t,A)
###############################################################################

