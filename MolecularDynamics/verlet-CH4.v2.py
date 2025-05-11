import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ==============================================
# PARÁMETROS FÍSICOS (UNIDADES REALISTAS)
# ==============================================
m_C = 12.011  # masa del carbono [uma]
m_H = 1.008   # masa del hidrógeno [uma]
k_bond = 450.0  # constante de enlace [kcal/mol/Å²]
r_eq = 1.09    # distancia de equilibrio [Å]
dt = 0.001    # paso de tiempo reducido [fs]
steps = 10000  # más pasos para mejor convergencia
damping = 0.3  # coeficiente de amortiguamiento

# ==============================================
# INICIALIZACIÓN DEL SISTEMA
# ==============================================
def init_positions():
    """Configuración inicial tetraédrica con pequeña perturbación"""
    positions = np.zeros((5, 3))
    # Posiciones tetraédricas ideales
    positions[1] = [ 1,  1,  1]  # H1
    positions[2] = [-1, -1,  1]  # H2
    positions[3] = [-1,  1, -1]  # H3
    positions[4] = [ 1, -1, -1]  # H4
    # Normalizar a distancia de equilibrio y añadir perturbación
    for i in range(1, 5):
        positions[i] = positions[i]/np.linalg.norm(positions[i]) * r_eq
        positions[i] += np.random.normal(0, 0.1, 3)  # Pequeña perturbación aleatoria
    return positions

# ==============================================
# CÁLCULO DE FUERZAS MEJORADO
# ==============================================
def compute_forces(pos, vel):
    forces = np.zeros_like(pos)
    C = pos[0]
    
    # Fuerzas de enlace C-H con protección numérica
    for i in range(1, 5):
        r_vec = pos[i] - C
        r = np.linalg.norm(r_vec)
        if r > 0.01:  # Evitar división por cero
            f_mag = -k_bond * (r - r_eq)
            forces[i] = f_mag * (r_vec/r)
            forces[0] -= f_mag * (r_vec/r)
    
    # Amortiguamiento más fuerte para mejor minimización
    forces -= damping * vel
    
    return forces

# ==============================================
# SIMULACIÓN CON VISUALIZACIÓN EN TIEMPO REAL
# ==============================================
def run_simulation_with_visualization():
    # Configuración inicial
    positions = init_positions()
    velocities = np.zeros_like(positions)
    masses = np.array([m_C] + [m_H]*4)
    
    # Preparar gráficos
    plt.ion()  # Modo interactivo
    fig = plt.figure(figsize=(18, 6))
    
    # Gráfico 3D
    ax1 = fig.add_subplot(131, projection='3d')
    ax1.set_title('Estructura Molecular')
    ax1.set_xlim(-2, 2)
    ax1.set_ylim(-2, 2)
    ax1.set_zlim(-2, 2)
    
    # Gráfico de distancias
    ax2 = fig.add_subplot(132)
    ax2.set_title('Distancias C-H')
    ax2.set_xlabel('Paso')
    ax2.set_ylabel('Distancia (Å)')
    ax2.axhline(r_eq, color='k', linestyle='--')
    
    # Gráfico de energía
    ax3 = fig.add_subplot(133)
    ax3.set_title('Energías del Sistema')
    ax3.set_xlabel('Paso')
    ax3.set_ylabel('Energía (kcal/mol)')
    
    # Datos para gráficos
    dist_history = [[] for _ in range(4)]
    pe_history = []
    ke_history = []
    te_history = []
    
    # Bucle principal de simulación
    for step in range(steps):
        # Velocity Verlet
        forces = compute_forces(positions, velocities)
        positions += velocities * dt + 0.5 * forces/masses[:, None] * dt**2
        new_forces = compute_forces(positions, velocities)
        velocities += 0.5 * (forces + new_forces)/masses[:, None] * dt
        
        # Calcular propiedades
        pe = 0.5 * k_bond * sum((np.linalg.norm(positions[i]-positions[0])-r_eq)**2 for i in range(1,5))
        ke = 0.5 * sum(m * np.sum(v**2) for m, v in zip(masses, velocities))
        
        # Actualizar datos históricos
        for i in range(4):
            dist_history[i].append(np.linalg.norm(positions[i+1]-positions[0]))
        pe_history.append(pe)
        ke_history.append(ke)
        te_history.append(pe + ke)
        
        # Actualizar gráficos cada 100 pasos
        if step % 100 == 0:
            # Limpiar gráficos
            ax1.cla()
            ax2.cla()
            ax3.cla()
            
            # Estructura 3D
            ax1.scatter(*positions[0], c='black', s=200, label='C')
            ax1.scatter(*positions[1:].T, c='red', s=100, label='H')
            for i in range(1,5):
                ax1.plot(*np.array([positions[0], positions[i]]).T, 'b-', alpha=0.5)
            ax1.set_title(f'Estructura Molecular (Paso {step})')
            ax1.set_xlim(-2, 2)
            ax1.set_ylim(-2, 2)
            ax1.set_zlim(-2, 2)
            ax1.legend()
            
            # Distancias
            for i in range(4):
                ax2.plot(dist_history[i], label=f'C-H{i+1}')
            ax2.axhline(r_eq, color='k', linestyle='--', label='Equilibrio')
            ax2.set_title('Evolución de Distancias C-H')
            ax2.legend()
            ax2.grid(True)
            
            # Energías
            ax3.plot(pe_history, 'r-', label='Potencial')
            ax3.plot(ke_history, 'b-', label='Cinética')
            ax3.plot(te_history, 'k--', label='Total')
            ax3.set_title('Energías del Sistema')
            ax3.legend()
            ax3.grid(True)
            
            plt.tight_layout()
            plt.draw()
            plt.pause(0.001)
    
    plt.ioff()
    return positions, dist_history, pe_history, ke_history, te_history

# ==============================================
# ANÁLISIS FINAL
# ==============================================
def analyze_final_state(final_pos):
    # Calcular distancias
    distances = [np.linalg.norm(final_pos[i]-final_pos[0]) for i in range(1,5)]
    
    # Calcular ángulos
    angles = []
    for i in range(1,5):
        for j in range(i+1,5):
            vec1 = final_pos[i] - final_pos[0]
            vec2 = final_pos[j] - final_pos[0]
            cos_theta = np.dot(vec1, vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))
            angles.append(np.degrees(np.arccos(np.clip(cos_theta, -1, 1))))
    
    return distances, angles

# ==============================================
# EJECUCIÓN PRINCIPAL
# ==============================================
if __name__ == "__main__":
    print("Iniciando simulación de minimización de energía...")
    
    # Ejecutar simulación con visualización
    final_pos, dist_hist, pe, ke, te = run_simulation_with_visualization()
    
    # Análisis final
    distances, angles = analyze_final_state(final_pos)
    
    print("\n=== RESULTADOS FINALES ===")
    print(f"Distancias C-H: {np.array(distances).round(4)} Å")
    print(f"Ángulos H-C-H: {np.array(angles).round(2)}°")
    print(f"Desviación media del tetraedro: {np.mean(np.abs(np.array(angles)-109.47)):.2f}°")
    print(f"Energía total final: {te[-1]:.4f} kcal/mol")
    
    # Mostrar gráficos finales
    plt.show()