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

print("\nWork Calculations from Hydrocarbon graphs, 6 Angstroms: Intrusion\n")

#pressure in bar = i in increments of 1 bar
pbulk = np.array([509.26782, 716.42083, 925.17334, 1135.4076, 1347.0362, 1559.9907, 1666.9472, 1774.215, 1989.6622, 2097.8317, 2206.2922, 2424.0697, 
                  2642.9639, 2862.9472, 3083.9946, 3194.9102])
nmol = np.array([0.73325, 1.417, 3.70975, 12.95775, 21.46275, 26.2545, 27.6, 28.112, 29.9795, 30.9, 31.06, 31.93625, 32.6225, 32.918, 33.52775, 33.79425])



an = an0 + (a0*pbulk) + (a1*pbulk**1.2) + (a2*pbulk**1.4)


#vbar will be expressed in angstrom*3 and number density in nm**-3
vbar = (25*25*25)/an
rho = 1/(vbar*1000)

vol = vbar*nmol

for i in range(len(vol) - 1):
    x = [vol[i], vol[i+1]]
    y = [pbulk[i], pbulk[i+1]]
    
    plt.fill_between(x, y, alpha=0.6)



degree = 7 

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