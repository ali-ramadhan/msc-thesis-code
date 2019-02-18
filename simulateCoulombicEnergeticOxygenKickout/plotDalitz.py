import sys
import numpy as np
import matplotlib.pyplot as plt
 
momentumFilename = str(sys.argv[1])
title = str(sys.argv[2])

# E_min = float(sys.argv[2])
# E_max = float(sys.argv[3])

momenta = np.loadtxt(momentumFilename)
#momenta = momenta*1e-22

p_O = momenta[:,0:3]
p_C = momenta[:,3:6]
p_S = momenta[:,6:9]

amu = 1.66053886e-27
e   = 1.602176565e-19

m_O = 15.999*amu
m_C = 12.011*amu
m_S = 32.065*amu

E_O = (1/e) * (np.sum(p_O**2, 1) / (2*m_O))
E_C = (1/e) * (np.sum(p_C**2, 1) / (2*m_C))
E_S = (1/e) * (np.sum(p_S**2, 1) / (2*m_S))
E_total = E_O + E_C + E_S

epsilon_O = E_O / E_total
epsilon_C = E_C / E_total
epsilon_S = E_S / E_total

X = (epsilon_O - epsilon_S) / np.sqrt(3)
Y = epsilon_C - (1.0/3.0)

# Crude code for making a KER cut.
# mask = []
# for i in range(len(X)):
# 	if E_total[i] < E_min or E_total[i] > E_max:
# 		mask = mask + [i]

# X = np.delete(X, mask, None)
# Y = np.delete(Y, mask, None)

# Plot data
fig1 = plt.figure()
plt.plot(X, Y,'.r')
plt.xlabel(r'$\epsilon_{O^+} - \epsilon_{S^+}/\sqrt{3}$', fontsize=22)
plt.ylabel(r'$\epsilon_{C^+} - 1/3$', fontsize=22)
 
# Estimate the 2D histogram
nbins = 150
H, xedges, yedges = np.histogram2d(X, Y, bins=nbins)

weightsH = np.ones_like(X) / np.amax(H)
H, xedges, yedges = np.histogram2d(X, Y, bins=nbins, weights=weightsH)

# H needs to be rotated and flipped
H = np.rot90(H)
H = np.flipud(H)
 
# Mask zeros
Hmasked = np.ma.masked_where(H==0, H) # Mask pixels with a value of zero
 
# Plot 2D histogram using pcolor
fig2 = plt.figure()
plt.pcolormesh(xedges, yedges, Hmasked)

plt.gca().set_aspect('equal', adjustable='box')

plt.title(title, fontsize=36, y=1.02)
plt.xlabel(r'$(\epsilon_{O} - \epsilon_{S})/\sqrt{3}$', fontsize=28)
plt.ylabel(r'$\epsilon_{C} - 1/3$', fontsize=28)

cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts [arbitrary units]', fontsize=28)
cbar.ax.tick_params(labelsize=16)

plt.tick_params(axis='both', which='major', labelsize=16)
plt.tick_params(axis='both', which='minor', labelsize=16)

plt.gcf().subplots_adjust(bottom=0.15)

plt.xlim([-0.2,0.35])
plt.ylim([-0.35,0.5])

plt.scatter([0.1381], [-0.2213], s=250, marker='x', edgecolor='red', linewidth=3)

fig2.savefig("%s Dalitz.pdf"% str(sys.argv[1]))