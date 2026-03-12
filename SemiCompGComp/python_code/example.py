######################################################################################################################
# Novel g-computation algorithms for time-varying actions with recurrent and semi-competing events
#   Sorensen D'Alessio et al.
#
# Paul Zivich (2026/03/12)
######################################################################################################################

import numpy as np
import scipy as sp
import pandas as pd
import delicatessen as deli
from delicatessen import MEstimator

from efuncs import ef_multistate_ice_gcomputation

print("Versions")
print("NumPy:", np.__version__)
print("SciPy:", sp.__version__)
print("pandas:", pd.__version__)
print("Delicatessen:", deli.__version__)
print()

##########################################################
# Loading data
d = pd.read_csv("single_sim_observed_data.csv", sep=r'\s+')
d['I'] = 1
d = d.fillna(-99)


##########################################################
# Setting up data for multi-state ICE g-computation

# Converting indicator into dummy variables for each Y type
for i in [1, 2, 3]:
    d['Y'+str(i)+'_1'] = np.where(d['Y'+str(i)] == 1, 1, 0)
    d['Y'+str(i)+'_2'] = np.where(d['Y'+str(i)] == 2, 1, 0)
    d['Y'+str(i)+'_3'] = np.where(d['Y'+str(i)] == 3, 1, 0)
    d['C' + str(i)] = np.where(d['Y'+str(i)] == -99, 1, 0)


# Design matrix under interventions being considered
d1, d0 = d.copy(), d.copy()
for i in [0, 1, 2]:
    d1['A'+str(i)] = 1
    d0['A'+str(i)] = 0


# Preparing data for input into estimating functions
c_ind = [np.asarray(d['C3']), np.asarray(d['C2']), np.asarray(d['C1'])]
X_mat = [np.asarray(d[['I', 'A2', 'L2', 'Y2']]),
         np.asarray(d[['I', 'A1', 'L1', 'Y1']]),
         np.asarray(d[['I', 'A0', 'L0']])]
y_mat = [np.asarray(d[['Y3_1', 'Y3_2', 'Y3_3']]),
         np.asarray(d[['Y2_1', 'Y2_2', 'Y2_3']]),
         np.asarray(d[['Y1_1', 'Y1_2', 'Y1_3']])]
X1_mat = [np.asarray(d1[['I', 'A2', 'L2', 'Y2']]),
          np.asarray(d1[['I', 'A1', 'L1', 'Y1']]),
          np.asarray(d1[['I', 'A0', 'L0']])]
X0_mat = [np.asarray(d0[['I', 'A2', 'L2', 'Y2']]),
          np.asarray(d0[['I', 'A1', 'L1', 'Y1']]),
          np.asarray(d0[['I', 'A0', 'L0']])]


##########################################################
# Estimation for always treat (A=1 at all k)

def psi(theta):
    return ef_multistate_ice_gcomputation(theta, y_array=y_mat,
                                          X_array=X_mat, Xa_array=X1_mat,
                                          c_array=c_ind)


init_vals = [0.4, 0.3, 0.3] + [0., ]*X_mat[0].shape[1]*2 + [0.1, ]*X_mat[1].shape[1]*2 + [0.2, ]*X_mat[2].shape[1]*2
estr = MEstimator(psi, init=init_vals)
estr.estimate()
ci = estr.confidence_intervals()
print("Always - Treat")
for i in range(3):
    print("mu(Y="+str(i)+")", np.round(estr.theta[i], 6), "  95% CI:", np.round(ci[i, :], 6))

print()

##########################################################
# Estimation for always treat (A=1 at all k)


def psi(theta):
    return ef_multistate_ice_gcomputation(theta, y_array=y_mat,
                                          X_array=X_mat, Xa_array=X0_mat,
                                          c_array=c_ind)


init_vals = [0.4, 0.3, 0.3] + [0., ]*X_mat[0].shape[1]*2 + [0.1, ]*X_mat[1].shape[1]*2 + [0.2, ]*X_mat[2].shape[1]*2
estr = MEstimator(psi, init=init_vals)
estr.estimate()
print("Never - Treat")
for i in range(3):
    print("mu(Y="+str(i)+")", np.round(estr.theta[i], 6), "  95% CI:", np.round(ci[i, :], 6))


##########################################################
# OUTPUT (Python 3.13.7)

# Versions
# NumPy: 2.3.5
# SciPy: 1.16.3
# pandas: 2.3.3
# Delicatessen: 4.1
#
# Always - Treat
# mu(Y=0) 0.453839   95% CI: [0.417628 0.490051]
# mu(Y=1) 0.405003   95% CI: [0.36911  0.440896]
# mu(Y=2) 0.141157   95% CI: [0.117023 0.165292]
#
# Never - Treat
# mu(Y=0) 0.389989   95% CI: [0.417628 0.490051]
# mu(Y=1) 0.478964   95% CI: [0.36911  0.440896]
# mu(Y=2) 0.131047   95% CI: [0.117023 0.165292]
