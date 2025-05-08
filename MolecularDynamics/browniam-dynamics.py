import numpy as np
import matplotlib.pyplot as plt

# Parámetros físicos y de simulación
k = 100.0              # Constante del resorte (kcal/mol/Å²)
r_eq = 0.96            # Distancia de equilibrio O-H (en Å)
dt = 0.001             # Paso de tiempo
gamma = 1.0            # Coeficiente de fricción
kT = 0.001             # Energía térmica (kcal/mol)
steps = 10000

# Posiciones iniciales (Oxígeno en el centro)
positions = np.array([
    [0.0, 0.0],         # Oxígeno
    [1.5, 0.0],         # Hidrógeno 1 (posición inicial alejada)
    [-1.0, 1.0]         # Hidrógeno 2
])

# Función para calcular la fuerza del resorte
def spring_force(r_vec, r0, k):
    dist = np.linalg.norm(r_vec)
    if dist == 0:
        return np.zeros_like(r_vec)
    direction = r_vec / dist
    return -k * (dist - r0) * direction

# Simulación
for step in range(steps):
    forces = np.zeros_like(positions)

    # Enlaces O-H1 y O-H2
    for i in [1, 2]:  # Hidrógenos
        r_vec = positions[i] - positions[0]
        f = spring_force(r_vec, r_eq, k)
        forces[i] += f
        forces[0] -= f

    # Movimiento Browniano (Euler-Maruyama)
    noise = np.random.normal(scale=np.sqrt(2 * gamma * kT * dt), size=positions.shape)
    positions += dt / gamma * forces + noise

# Visualización final
plt.figure(figsize=(5, 5))
plt.plot(positions[[0, 1], 0], positions[[0, 1], 1], 'b-', label='O-H1')
plt.plot(positions[[0, 2], 0], positions[[0, 2], 1], 'r-', label='O-H2')
plt.plot(positions[:, 0], positions[:, 1], 'ko', label='Átomos')
plt.legend()
plt.title("Molécula de agua relajada (2D)")
plt.xlabel("x (Å)")
plt.ylabel("y (Å)")
plt.axis('equal')
plt.grid(True)
plt.show()
