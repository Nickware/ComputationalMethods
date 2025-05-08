import numpy as np
import matplotlib.pyplot as plt

# Parámetros físicos (unidades arbitrarias)
m_O = 16.0  # masa del oxígeno
m_C = 12.0  # masa del carbono
k = 500.0   # constante del resorte
r_eq = 1.16  # longitud de equilibrio C=O

# Vector de masas
masses = np.array([m_O, m_C, m_O])

# Condiciones iniciales
positions = np.array([-r_eq, 0.0, r_eq])  # posiciones en equilibrio
velocities = np.array([0.0, 0.0, 0.0])    # sin velocidad inicial

# Introduce una perturbación
positions[2] += 0.2  # desplaza O2

# Parámetros de integración
dt = 0.001
steps = 5000

# Función de fuerzas corregida
def compute_forces(pos):
    f = np.zeros_like(pos)
    
    # Resorte O1-C
    dx1 = pos[1] - pos[0]
    f12 = k * (dx1 - r_eq)
    f[0] += f12
    f[1] -= f12
    
    # Resorte C-O2
    dx2 = pos[2] - pos[1]
    f23 = k * (dx2 - r_eq)
    f[1] += f23
    f[2] -= f23
    
    return f

# Algoritmo Velocity Verlet
trajectory = [positions.copy()]
forces_old = compute_forces(positions)

for _ in range(steps):
    # Paso 1: Actualizar posiciones
    positions += velocities * dt + 0.5 * forces_old / masses * dt**2
    
    # Paso 2: Calcular nuevas fuerzas
    forces_new = compute_forces(positions)
    
    # Paso 3: Actualizar velocidades
    velocities += 0.5 * (forces_old + forces_new) / masses * dt
    
    # Guardar estado y preparar siguiente iteración
    trajectory.append(positions.copy())
    forces_old = forces_new.copy()

trajectory = np.array(trajectory)

# Graficar resultados
plt.figure(figsize=(10, 5))
plt.plot(trajectory[:, 0], label='Oxígeno 1 (O1)', color='red')
plt.plot(trajectory[:, 1], label='Carbono (C)', color='black')
plt.plot(trajectory[:, 2], label='Oxígeno 2 (O2)', color='blue')
plt.xlabel("Paso de tiempo")
plt.ylabel("Posición (1D)")
plt.title("Oscilaciones de CO₂ con Velocity Verlet")
plt.legend()
plt.grid(True)
plt.show()
