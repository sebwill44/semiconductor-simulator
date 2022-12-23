# importing the modules
import numpy as np
import matplotlib.pyplot as plt

def sn(res):
    if res > 10 or abs(res)< 0.001:
        return np.format_float_scientific(res, precision=10, sign=True)
    else:
        return np.format_float_positional(res, precision=10, sign=True)