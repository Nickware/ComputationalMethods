import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from IPython.display import clear_output
import time

# ======================
# PARÁMETROS FÍSICOS
# ======================
m_C = 12.011  # masa del carbono [uma]
m_H = 1.008   # masa del hidrógeno [uma]
k_bond = 450.0  # constante de enlace [kcal/mol/Å²]
r_eq = 1.09    # distancia de equilibrio [Å]
dt = 0.002    # paso de tiempo reducido [fs]
steps = 3000   # pasos totales

# ======================
# CONFIGURACIÓN INICIAL
# ======================
np.random.seed(42)
positions = np.zeros((5, 3))  # [C, H1, H2, H3, H4]
positions[1:] = np.random.uniform(-0.3, 0.3, (4, 3))  # H cerca del centro
velocities = np.zeros_like(positions)
masses = np.array([m_C] + [m_H]*4)  # Definición de masas

# ======================
# FUNCIÓN DE FUERZAS
# ======================
def compute_forces(pos):
    forces = np.zeros_like(pos)
    C = pos[0]
    
    # Fuerzas de enlace C-H (armónico)
    for i in range(1, 5):
        r_vec = pos[i] - C
        r = max(np.linalg.norm(r_vec), 0.01)  # Evitar división por cero
        f_mag = -k_bond * (r - r_eq)
        forces[i] = f_mag * (r_vec/r)
        forces[0] -= f_mag * (r_vec/r)
    
    # Amortiguamiento viscoso (γ = 0.1)
    forces -= 0.1 * velocities
    
    return forces

# ======================
# SIMULACIÓN
# ======================
trajectory = []
energies = []
dist_history = []

plt.figure(figsize=(12, 5))
start_time = time.time()

for step in range(steps):
    # Velocity Verlet
    forces = compute_forces(positions)
    positions += velocities * dt + 0.5 * forces/masses[:, None] * dt**2
    new_forces = compute_forces(positions)
    velocities += 0.5 * (forces + new_forces)/masses[:, None] * dt
    
    # Energías y almacenamiento
    pe = 0.5 * k_bond * sum((np.linalg.norm(positions[i]-positions[0])-r_eq)**2 for i in range(1,5))
    ke = 0.5 * sum(m * np.sum(v**2) for m, v in zip(masses, velocities))
    
    if step % 10 == 0:
        trajectory.append(positions.copy())
        energies.append((pe, ke, pe+ke))
        dist_history.append([np.linalg.norm(positions[i]-positions[0]) for i in range(1,5)])
    
    # Visualización cada 50 pasos
    if step % 50 == 0:
        clear_output(wait=True)
        elapsed = time.time() - start_time
        print(f"Paso {step}/{steps} | Tiempo: {elapsed:.1f}s | PE: {pe:.2f} kcal/mol")
        
        # Gráfico 3D
        ax = plt.subplot(131, projection='3d')
        ax.scatter(*positions[0], c='black', s=200, label='C')
        ax.scatter(*positions[1:].T, c='red', s=100, label='H')
        for i in range(1,5):
            ax.plot(*np.array([positions[0], positions[i]]).T, 'b-', alpha=0.4)
        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-1.5, 1.5)
        ax.set_zlim(-1.5, 1.5)
        ax.set_title(f'Paso {step}')
        
        # Distancias C-H
        plt.subplot(132)
        plt.plot(np.array(dist_history))
        plt.axhline(r_eq, color='black', linestyle='--')
        plt.title('Distancias C-H')
        plt.ylabel('Å')
        plt.xlabel('Paso (x10)')
        
        # Energías
        plt.subplot(133)
        if energies:
            pe, ke, te = zip(*energies)
            plt.plot(pe, label='Potencial', c='red')
            plt.plot(ke, label='Cinética', c='blue')
            plt.plot(te, label='Total', c='black', linestyle='--')
            plt.title('Energías')
            plt.xlabel('Paso (x10)')
            plt.legend()
        
        plt.tight_layout()
        plt.pause(0.001)

# ======================
# ANÁLISIS FINAL
# ======================
final_pos = trajectory[-1]
distances = [np.linalg.norm(final_pos[i]-final_pos[0]) for i in range(1,5)]
angles = []
for i in range(1,5):
    for j in range(i+1,5):
        vec1 = final_pos[i] - final_pos[0]
        vec2 = final_pos[j] - final_pos[0]
        cos_theta = np.dot(vec1, vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))
        angles.append(np.degrees(np.arccos(np.clip(cos_theta, -1, 1))))

print("\n=== RESULTADOS FINALES ===")
print(f"Distancias C-H: {np.array(distances).round(3)} Å")
print(f"Ángulos H-C-H: {np.array(angles).round(1)}°")
print(f"\nDesviación media del tetraedro: {np.mean(np.abs(np.array(angles)-109.47)):.1f}°")
print(f"Distancia promedio final: {np.mean(distances):.3f} ± {np.std(distances):.3f} Å")