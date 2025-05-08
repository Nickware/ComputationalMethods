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

# Introduce una perturbación para ver movimiento oscilatorio
positions[2] += 0.2  # desplaza O2

# Parámetros de integración
dt = 0.001
steps = 5000

# Función para calcular fuerzas (CORRECCIÓN APLICADA)
def compute_forces(pos):
    f = np.zeros_like(pos)
    
    # Resorte entre O1 y C (CORRECCIÓN: eliminado signo negativo)
    dx1 = pos[1] - pos[0]
    f12 = k * (dx1 - r_eq)  # ← Aquí estaba el error principal
    f[0] += f12
    f[1] -= f12
    
    # Resorte entre C y O2 (CORRECCIÓN: eliminado signo negativo)
    dx2 = pos[2] - pos[1]
    f23 = k * (dx2 - r_eq)  # ← Aquí estaba el error principal
    f[1] += f23
    f[2] -= f23
    
    return f

# Inicialización de Leapfrog (medio paso atrás)
forces = compute_forces(positions)
velocities -= 0.5 * dt * forces / masses

# Guardar trayectoria
trajectory = [positions.copy()]

# Bucle principal de integración
for _ in range(steps):
    positions += dt * velocities
    forces = compute_forces(positions)
    velocities += dt * forces / masses
    trajectory.append(positions.copy())

trajectory = np.array(trajectory)

# Graficar resultados
plt.figure(figsize=(10, 5))
plt.plot(trajectory[:, 0], label='Oxígeno 1 (O1)', color='red')
plt.plot(trajectory[:, 1], label='Carbono (C)', color='black')
plt.plot(trajectory[:, 2], label='Oxígeno 2 (O2)', color='blue')
plt.xlabel("Paso de tiempo")
plt.ylabel("Posición (1D)")
plt.title("Oscilaciones de la molécula de CO₂ con Leapfrog (CORREGIDO)")
plt.legend()
plt.grid(True)
plt.show()
