import matplotlib.pyplot as plt
import numpy as np

L_h = L_h_a = [-1 + 2 * 0.02 * i for i in range(51)]
L_d = [0.02 * i for i in range(51)]
L_c = [0.001 * i for i in range(51)]

Result = np.load('results/ES_delta_mu.npy')
plt.figure()
plt.pcolormesh(L_d, L_h, Result)
plt.colorbar()
plt.xlabel('$\delta$')
plt.ylabel('$\mu$')
plt.clim(0, 1)

Result = np.load('results/ES_delta_h_a.npy')
plt.figure()
plt.pcolormesh(L_d, L_h, Result)
plt.colorbar()
plt.xlabel('$\delta$')
plt.ylabel('$h_a$')
plt.clim(0, 1)

Result = np.load('results/ES_delta_r.npy')
plt.figure()
plt.pcolormesh(L_d, L_d, Result)
plt.colorbar()
plt.xlabel('$\delta$')
plt.ylabel('$r$')
plt.clim(0, 1)

L_c = [0.001 * i for i in range(51)]

Result = np.load('results/ES_c_f_c_r.npy')
plt.figure()
plt.pcolormesh(L_c, L_c, Result)
plt.colorbar()
plt.xlabel('$c_f$')
plt.ylabel('$c_r$')
plt.clim(0, 1)

Result = np.load('results/ES_mu_h_a.npy')
plt.figure()
plt.pcolormesh(L_h, L_h, Result.T)
plt.colorbar()
plt.xlabel('$h_a$')
plt.ylabel('$\mu$')
plt.clim(0, 1)
