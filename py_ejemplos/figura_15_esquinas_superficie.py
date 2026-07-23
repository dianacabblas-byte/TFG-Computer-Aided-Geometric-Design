import numpy as np
import matplotlib.pyplot as plt

def calcular_vector_nodos(n, grado, a, b):
    num_internos = n - grado
    h = (b - a) / (num_internos + 1)
    inicio = [a] * (grado + 1)
    medio = [a + i * h for i in range(1, num_internos + 1)]
    fin = [b] * (grado + 1)
    return inicio + medio + fin

def cox_de_boor(u, puntos_control, U, grado):
    n = len(puntos_control) - 1
    r = grado
    while r <= n:
        if u < U[r + 1] or r == n:
            break
        r += 1
    Q = np.array(puntos_control, dtype=float)
    for j in range(1, grado + 1):
        for i in range(r - grado, r - j + 1):
            dif = U[i + grado + 1] - U[i + j]
            if abs(dif) < 1e-9:
                Q[i,:] = 0.0
            else:
                factor_1 = (U[i + grado + 1] - u) / dif
                factor_2 = (u - U[i + j]) / dif
                Q[i,:] = factor_1 * Q[i,:] + factor_2 * Q[i+1,:]
    return Q[r - grado]

def evaluar_superficie(red_control, grado_u, grado_v, num_u, num_v):
    a, b = 0.0, 1.0
    c, d = 0.0, 1.0
    n = len(red_control) - 1
    m = len(red_control[0]) - 1
    
    U_u = calcular_vector_nodos(n, grado_u, a, b)
    V_v = calcular_vector_nodos(m, grado_v, c, d)
    
    tiempos_u = np.linspace(a, b, num_u + 1)
    tiempos_v = np.linspace(c, d, num_v + 1)
    
    X = np.zeros((num_u + 1, num_v + 1))
    Y = np.zeros((num_u + 1, num_v + 1))
    Z = np.zeros((num_u + 1, num_v + 1))
    
    curvas_v = np.zeros((num_v + 1, n + 1, 3))
    
    for j, val_v in enumerate(tiempos_v):
        for i in range(n + 1):
            curvas_v[j, i] = cox_de_boor(val_v, red_control[i], V_v, grado_v)
            
    for i, val_u in enumerate(tiempos_u):
        for j, val_v in enumerate(tiempos_v):
            puntos_intermedios = curvas_v[j]
            punto_final = cox_de_boor(val_u, puntos_intermedios, U_u, grado_u)
            X[i, j] = punto_final[0]
            Y[i, j] = punto_final[1]
            Z[i, j] = punto_final[2]
            
    return X, Y, Z

if __name__ == "__main__":
    try:
        # 1. Definición de una red de control de ejemplo (4x4)
        red_control = np.array([
            [[-2, -2, 0], [-2, -1, 4], [-2,  1, 4], [-2,  2, 0]],
            [[-1, -2, 3], [-1, -1, 6], [-1,  1, 6], [-1,  2, 3]],
            [[ 1, -2, 3], [ 1, -1, 6], [ 1,  1, 6], [ 1,  2, 3]],
            [[ 2, -2, 0], [ 2, -1, 4], [ 2,  1, 4], [ 2,  2, 0]]
        ], dtype=float)
        
        grado_u, grado_v = 3, 3
        res_u, res_v = 40, 40
        
        # 2. Cálculo de la superficie
        X_surf, Y_surf, Z_surf = evaluar_superficie(red_control, grado_u, grado_v, res_u, res_v)
        
        print("Generando gráfica de esquinas...")
        
        # 3. Dibujo
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Pintamos la superficie
        ax.plot_surface(X_surf, Y_surf, Z_surf, cmap='viridis', alpha=0.7, edgecolor='none')
        
        # Pintamos la malla de alambre de la red de control
        X_ctrl = red_control[:, :, 0]
        Y_ctrl = red_control[:, :, 1]
        Z_ctrl = red_control[:, :, 2]
        ax.plot_wireframe(X_ctrl, Y_ctrl, Z_ctrl, color='gray', alpha=0.5, linestyle='--')
        
        # Pintamos los puntos de control diferenciando las esquinas
        n = len(red_control) - 1
        m = len(red_control[0]) - 1
        
        for i in range(n + 1):
            for j in range(m + 1):
                esquina = (i == 0 and j == 0) or (i == 0 and j == m) or (i == n and j == 0) or (i == n and j == m)
                if esquina:
                    # Esquinas en rojo y gigantes
                    ax.scatter(red_control[i, j, 0], red_control[i, j, 1], red_control[i, j, 2], 
                               color='red', s=120, zorder=10)
                else:
                    # Resto de puntos en gris
                    ax.scatter(red_control[i, j, 0], red_control[i, j, 1], red_control[i, j, 2], 
                               color='gray', s=30, alpha=0.8)
        
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_box_aspect([1, 1, 0.6])
        ax.legend(fontsize=11)
        
        plt.tight_layout()
        plt.show()

    except Exception as e:
        print("\n[ERROR] Hubo un problema durante la ejecución:", e)