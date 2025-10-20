import numpy as np
from scipy.integrate import dblquad
"""
cylinde de r = 2.5 mm sur 8 mm de hauteur
cone sur 10 mm de hauteur qui va de r = 2.5 mm à r = 1.15 mm
demi-sphere sur 0.9 mm (rayon ?)
épaisseur du tube 0.5 mm

la focale tape sur le tube où le diamètre extérieur est 5.41, donc un diamètre intérieur de 4.41
"""

# Paramètres (en mm)
a = 12.9*0.5  # demi-grand axe de l'ellipsoïde
b = 1.9*0.5   # demi-petit axe de l'ellipsoïde
r = 4.41*0.5       # rayon du cylindre

# Limite en x pour l'intersection
x_max = min(a, np.sqrt(a**2 * b**2 * (r**2 - b**2) / (a**2 - b**2)))

# Borne supérieure pour y, en tenant compte de l'intersection réelle
def y_bound(x):
    y_cylinder = np.sqrt(r**2 - x**2)
    y_ellipsoid = b * np.sqrt(1 - x**2/a**2)
    return np.minimum(y_cylinder, y_ellipsoid)

# Fonction à intégrer : hauteur locale de l'ellipsoïde selon z
def integrand(y, x):
    return 2 * b * np.sqrt(1 - x**2/a**2 - y**2/b**2)

# Bornes d'intégration corrigées
def y_lower(x):
    return -y_bound(x)

def y_upper(x):
    return y_bound(x)

# Intégration numérique
volume_mm3, _ = dblquad(integrand, -x_max, x_max, y_lower, y_upper)

# Conversion en microlitres (1 mm³ = 1 µL)
volume_µL = volume_mm3

print(f"Volume commun : {volume_mm3:.4f} mm³")
print(f"Volume commun : {volume_µL:.4f} µL")





