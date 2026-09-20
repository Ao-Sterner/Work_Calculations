import numpy as np
import matplotlib.pyplot as plt


#c partial molar volume vbar for methanol at 298.15 K
###Constants:
a0=0.0771323
a1=-0.0184755
a2=0.00120936
an0=223.7
akT=.41145e-20
vbar0=(25*25*25)/an0

print("\nWork Calculations from Hydrocarbon graphs, 9.5 Angstroms: Intrusion\n")

#pressure in bar = i in increments of 1 bar
pbulk = np.array([50.184917, 75.872596, 100.6008, 125.87796, 151.19697, 201.952, 252.85251, 303.88897, 612.6353, 925.17334])
nmol = np.array([1, 0.3, 44.1, 46, 47, 48.9, 49.64, 50.6, 55.1, 58.6])

an = an0 + (a0*pbulk) + (a1*pbulk**1.2) + (a2*pbulk**1.4)


#vbar will be expressed in angstrom*3 and number density in nm**-3
vbar = (25*25*25)/an
rho = 1/(vbar*1000)

vol = vbar*nmol

for i in range(len(vol) - 1):
    x = [vol[i], vol[i+1]]
    y = [pbulk[i], pbulk[i+1]]
    
    plt.fill_between(x, y, alpha=0.6)



plt.scatter(vol, pbulk, label='Original points', color='blue')
plt.xlabel('<N>*Vbar')
plt.ylabel('Pressure (bar)')


print("<N>*Vbar points: ",vol,"\n")


plt.scatter(vol,pbulk,label = 'My points')
plt.xlabel('<N>*Vbar'), plt.ylabel('Pressure (bar)')


W = np.trapezoid(pbulk, vol)
print(f"Trapezoidal integral from {min(vol):.3f} to {max(vol):.3f} = ", W, " Joules\n\n")


# If you want work over a region of volume
##EXAMPLE:
Work = np.trapezoid(pbulk[2:4], vol[2:4])  # slice arrays

print(f"Trapezoidal integral from {vol[2]:.3f} to {vol[3]:.3f} = {Work} Joules\n\n")


plt.show()