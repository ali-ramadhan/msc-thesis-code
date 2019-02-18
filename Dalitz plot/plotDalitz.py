import numpy as np
import matplotlib.pyplot as plt
 
finX = open("momentaAr6+X_rotated.txt", "r")
finY = open("momentaAr6+Y_rotated.txt", "r")

X = []
Y = []

for line in finX.readlines():
	X.append(float(line.strip()))

for line in finY.readlines():
	Y.append(float(line.strip()))

# Plot data
fig1 = plt.figure()
plt.plot(X, Y,'.r')
plt.xlabel(r'$\frac{\epsilon_{O^+} - \epsilon_{S^+}}{\sqrt{3}}$', fontsize=22)
plt.ylabel(r'$\epsilon_{C^+} - \frac{1}{3}$', fontsize=22)
 
# Estimate the 2D histogram
nbins = 150
H, xedges, yedges = np.histogram2d(X, Y, bins=nbins)
 
# H needs to be rotated and flipped
H = np.rot90(H)
H = np.flipud(H)
 
# Mask zeros
Hmasked = np.ma.masked_where(H==0, H) # Mask pixels with a value of zero
 
# Plot 2D histogram using pcolor
fig2 = plt.figure()
plt.pcolormesh(xedges, yedges, Hmasked)
plt.xlabel(r'$\frac{\epsilon_{O^+} - \epsilon_{S^+}}{\sqrt{3}}$', fontsize=22)
plt.ylabel(r'$\epsilon_{C^+} - \frac{1}{3}$', fontsize=22)
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')

plt.savefig('momentaAr6+_rotated.pdf')

finX.close()
finY.close()