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

print("\nWork Calculations from Graphene graphs, 7 Angstroms: Intrusion\n")

#pressure in bar = i in increments of 1 bar
pbulk = np.array([1.0002392, 201.95198, 406.34201, 612.6352, 820.60525, 1030.111, 1241.0519, 1453.351, 1666.9466, 1881.7877, 2097.8309, 2315.0386, 
                  2533.3779, 2752.8198, 2973.3379, 3194.9086, 3417.5103])
nmol = np.array([34.45, 34.5, 35, 35.35, 35.45, 35.5, 36, 36, 36.5, 36.5, 36.9, 36.95, 36.99, 37.3, 37.2, 37.3, 37.9])

an = an0 + (a0*pbulk) + (a1*pbulk**1.2) + (a2*pbulk**1.4)


#vbar will be expressed in angstrom*3 and number density in nm**-3
vbar = (25*25*25)/an
rho = 1/(vbar*1000)

vol = vbar*nmol

for i in range(len(vol) - 1):
    x = [vol[i], vol[i+1]]
    y = [pbulk[i], pbulk[i+1]]
    
    plt.fill_between(x, y, alpha=0.6)



degree = 3 

coefficients = np.polyfit(vol, pbulk, degree)
vol_smooth = np.linspace(min(vol), max(vol), 100)
pbulk_fitted = np.polyval(coefficients, vol_smooth)

plt.scatter(vol, pbulk, label='Original points', color='blue')
plt.plot(vol_smooth, pbulk_fitted, color='blue', label=f'Degree {degree} fit', linewidth=2)
plt.xlabel('<N>*Vbar')
plt.ylabel('Pressure (bar)')



#print(f"Coefficients: {coefficients}")

integral_coefficients = np.polyint(coefficients)

print("<N>*Vbar points: ",vol,"\n")



# Work integral from min(vol) to max(vol)
integral_value = np.polyval(integral_coefficients, max(vol)) - \
                 np.polyval(integral_coefficients, min(vol))

print(f"Definite integral from {min(vol):.3f} to {max(vol):.3f} = {integral_value:.6f} Joules")


plt.scatter(vol,pbulk,label = 'My points')
plt.xlabel('<N>*Vbar'), plt.ylabel('Pressure (bar)')


W = np.trapezoid(pbulk, vol)
print(f"Trapezoidal integral from {min(vol):.3f} to {max(vol):.3f} = ", W, " Joules\n\n")

plt.show()