import numpy as np
import matplotlib.pyplot as plt

# Parámetros físicos y de simulación
k = 100.0              # Constante del resorte (kcal/mol/Å²)
r_eq = 0.96            # Longitud de enlace O-H (Å)
dt = 0.001             # Paso de tiempo
gamma = 1.0            # Coeficiente de fricción
kT = 0.001             # Energía térmica (kcal/mol)
steps = 10000

# Posiciones iniciales
positions = np.array([
    [0.0, 0.0],         # Oxígeno
    [1.5, 0.0],         # H1 desplazado
    [-1.0, 1.0]         # H2
])

trajectory = [positions.copy()]

# Función de fuerza por resorte
def spring_force(r_vec, r0, k):
    dist = np.linalg.norm(r_vec)
    if dist == 0:
        return np.zeros_like(r_vec)
    direction = r_vec / dist
    return -k * (dist - r0) * direction

# Simulación
for step in range(steps):
    forces = np.zeros_like(positions)

    # O-H1 y O-H2
    for i in [1, 2]:
        r_vec = positions[i] - positions[0]
        f = spring_force(r_vec, r_eq, k)
        forces[i] += f
        forces[0] -= f

    # Movimiento tipo Langevin overdamped
    noise = np.random.normal(scale=np.sqrt(2 * gamma * kT * dt), size=positions.shape)
    positions += dt / gamma * forces + noise

    trajectory.append(positions.copy())

# Convertir trayectoria a array
trajectory = np.array(trajectory)

# Graficar trayectoria en tiempo para cada átomo (coordenada x)
plt.figure(figsize=(10, 6))
labels = ['Oxígeno', 'Hidrógeno 1', 'Hidrógeno 2']
for i in range(3):
    plt.plot(trajectory[:, i, 0], label=f'{labels[i]} - x')
plt.xlabel('Tiempo (pasos)')
plt.ylabel('Posición x (Å)')
plt.title('Evolución de posición x con el tiempo')
plt.legend()
plt.grid()
plt.show()