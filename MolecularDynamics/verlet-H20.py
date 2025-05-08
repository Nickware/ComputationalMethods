import numpy as np
import matplotlib.pyplot as plt

# Parámetros físicos
k = 100.0          # Constante del resorte (kcal/mol/Å²)
r_eq = 0.96        # Longitud de enlace O-H (Å)
dt = 0.001         # Paso de tiempo
steps = 10000
gamma = 1.0        # Coeficiente de fricción (disipación)

# Posiciones iniciales (Oxígeno en el centro)
positions = np.array([
    [0.0, 0.0],       # Oxígeno
    [1.5, 0.0],       # Hidrógeno 1 (desplazado)
    [-1.0, 1.0]       # Hidrógeno 2
])

# Velocidades iniciales
velocities = np.zeros_like(positions)

# Masas (aproximadas)
masses = np.array([16.0, 1.0, 1.0]).reshape(-1, 1)

# Almacenar historia
pos_history = [positions.copy()]
vel_history = [velocities.copy()]

# Función de fuerza por resorte
def spring_force(r_vec, r0, k):
    dist = np.linalg.norm(r_vec)
    if dist == 0:
        return np.zeros_like(r_vec)
    direction = r_vec / dist
    return -k * (dist - r0) * direction

# Simulación con dinámica molecular y disipación (Langevin simplificado)
for step in range(steps):
    forces = np.zeros_like(positions)

    # Fuerzas del oxígeno a cada hidrógeno
    for i in [1, 2]:
        r_vec = positions[i] - positions[0]
        f = spring_force(r_vec, r_eq, k)
        forces[i] += f
        forces[0] -= f

    # Agregar fricción (fuerza viscosa proporcional a velocidad)
    forces -= gamma * velocities

    # Actualización por integración tipo Verlet con velocidad (Velocity-Verlet)
    accelerations = forces / masses
    velocities += accelerations * dt
    positions += velocities * dt

    pos_history.append(positions.copy())
    vel_history.append(velocities.copy())

# Convertir a arrays para análisis
pos_history = np.array(pos_history)  # shape: (steps, 3, 2)
vel_history = np.array(vel_history)

# Graficar evolución de posición x en el tiempo para cada átomo
t = np.arange(steps + 1) * dt
plt.figure(figsize=(10, 6))
for i, label in enumerate(['Oxígeno', 'Hidrógeno 1', 'Hidrógeno 2']):
    plt.plot(t, pos_history[:, i, 0], label=f'{label} - x')
plt.xlabel('Tiempo (ps)')
plt.ylabel('Posición x (Å)')
plt.title('Relajación de posiciones - coordenada x')
plt.legend()
plt.grid()
plt.show()

# Graficar estado de fase (posición vs velocidad en x)
plt.figure(figsize=(8, 6))
for i, label in enumerate(['Oxígeno', 'Hidrógeno 1', 'Hidrógeno 2']):
    plt.plot(pos_history[:, i, 0], vel_history[:, i, 0], label=label)
plt.xlabel('Posición x (Å)')
plt.ylabel('Velocidad x (Å/ps)')
plt.title('Estado de fase - coordenada x')
plt.legend()
plt.grid()
plt.show()