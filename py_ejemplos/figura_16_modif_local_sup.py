import numpy as np
import matplotlib.pyplot as plt

def calcular_vector_nodos(n, grado, a, b):
    num_internos = n - grado
    h = (b - a) / (num_internos + 1)
    
    # Repetimos el punto a (grado + 1) veces para que pase exactamente por él
    inicio = [a] * (grado + 1)
    
    # Puntos internos del intervalo
    medio = [a + i * h for i in range(1, num_internos + 1)]
    
    # Repetimos el punto b (grado + 1) veces para que pase exactamente por él
    fin = [b] * (grado + 1)
    
    return inicio + medio + fin # Devolvemos el vector ensamblado

def cox_de_boor(u, puntos_control, U, grado):
    n = len(puntos_control) - 1
    r = grado
    
    # Buscamos r tq el valor esté en el intervalo [u_r, u_(r+1)]
    while r <= n:
        if u < U[r + 1] or r == n:
            break
        r += 1
        
    # Matriz auxiliar para calcular la evaluación
    Q = np.array(puntos_control, dtype=float)

    # Aplicamos Cox-de Boor
    for j in range(1, grado + 1):
        for i in range(r - grado, r - j + 1):
            
            # Denominador de las formulas:
            dif = U[i + grado + 1] - U[i + j]
            if abs(dif) < 1e-9:  # si la diferencia es nula
                Q[i,:] = 0.0
            else:
                factor_1 = (U[i + grado + 1] - u) / dif
                factor_2 = (u - U[i + j]) / dif
                
                Q[i,:] = factor_1 * Q[i,:] + factor_2 * Q[i+1,:]

    return Q[r - grado]

def evaluar_superficie(red_control, grado_u, grado_v, num_u, num_v):
    """Aplica Cox-de Boor en las dos dimensiones (u y v)"""
    a, b = 0.0, 1.0
    c, d = 0.0, 1.0
    
    n = len(red_control) - 1      # Filas
    m = len(red_control[0]) - 1   # Columnas
    
    # Cálculos previos de nodos e intervalos
    U_u = calcular_vector_nodos(n, grado_u, a, b)
    V_v = calcular_vector_nodos(m, grado_v, c, d)
    
    tiempos_u = np.linspace(a, b, num_u + 1)
    tiempos_v = np.linspace(c, d, num_v + 1)
    
    # Pre-creamos las matrices para la superficie final
    X = np.zeros((num_u + 1, num_v + 1))
    Y = np.zeros((num_u + 1, num_v + 1))
    Z = np.zeros((num_u + 1, num_v + 1))
    
    # Matriz auxiliar
    curvas_v = np.zeros((num_v + 1, n + 1, 3))
    
    # PASO 1 (Dirección V): Evaluamos cada fila de la red para cada instante de V
    for j, val_v in enumerate(tiempos_v):
        for i in range(n + 1):
            curvas_v[j, i] = cox_de_boor(val_v, red_control[i], V_v, grado_v)
            
    # PASO 2 (Dirección U): Evaluamos las columnas precalculadas para cada instante de U
    for i, val_u in enumerate(tiempos_u):
        for j, val_v in enumerate(tiempos_v):
            
            # Cogemos la columna intermedia de la matriz auxiliar
            puntos_intermedios = curvas_v[j]
            
            # Aplicamos el paso final
            punto_final = cox_de_boor(val_u, puntos_intermedios, U_u, grado_u)
            
            # Asignamos el punto directamente a la celda de la matriz
            X[i, j] = punto_final[0]
            Y[i, j] = punto_final[1]
            Z[i, j] = punto_final[2]
        
    return X, Y, Z

if __name__ == "__main__":
    
    try:
        # 1. Definición de la red de control 
        red_orig = np.array([
            [[0, 0, 0], [0, 1, 0], [0, 2, 0], [0, 3, 0], [0, 4, 0], [0, 5, 0], [0, 6, 0]],
            [[1, 0, 0], [1, 1, 0], [1, 2, 0], [1, 3, 0], [1, 4, 0], [1, 5, 0], [1, 6, 0]],
            [[2, 0, 0], [2, 1, 0], [2, 2, 0], [2, 3, 0], [2, 4, 0], [2, 5, 0], [2, 6, 0]],
            [[3, 0, 0], [3, 1, 0], [3, 2, 0], [3, 3, 0], [3, 4, 0], [3, 5, 0], [3, 6, 0]],
            [[4, 0, 0], [4, 1, 0], [4, 2, 0], [4, 3, 0], [4, 4, 0], [4, 5, 0], [4, 6, 0]],
            [[5, 0, 0], [5, 1, 0], [5, 2, 0], [5, 3, 0], [5, 4, 0], [5, 5, 0], [5, 6, 0]],
            [[6, 0, 0], [6, 1, 0], [6, 2, 0], [6, 3, 0], [6, 4, 0], [6, 5, 0], [6, 6, 0]]
        ])
        
        red_mod = np.copy(red_orig)
        red_mod[3, 3, 2] = 8.0 
        
        n = len(red_orig) - 1
        m = len(red_orig[0]) - 1
        grado_u, grado_v = 3, 3
        
        # 2. Cálculo de la superficie
        res_u, res_v = 50, 50

        X_orig, Y_orig, Z_orig = evaluar_superficie(red_orig, grado_u, grado_v, res_u, res_v)
        X_mod, Y_mod, Z_mod = evaluar_superficie(red_mod, grado_u, grado_v, res_u, res_v)

        print("Generando gráfica...")
        
        # 3. Dibujo
        fig = plt.figure(figsize=(15, 6))

        # --- Subplot Izquierdo: Original ---
        ax1 = fig.add_subplot(121, projection='3d')
        ax1.plot_surface(X_orig, Y_orig, Z_orig, cmap='viridis', alpha=0.5, edgecolor='none')
        
        X_ctrl = red_orig[:, :, 0]
        Y_ctrl = red_orig[:, :, 1]
        Z_ctrl_orig = red_orig[:, :, 2]
        Z_ctrl_mod = red_mod[:, :, 2]
        
        ax1.plot_wireframe(X_ctrl, Y_ctrl, Z_ctrl_orig, color='black', alpha=0.5, linestyle='--')
        ax1.scatter(X_ctrl, Y_ctrl, Z_ctrl_orig, color='gray', s=20)
        
        ax1.set_title('Superficie Original', fontsize=12)
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_zlabel('z')
        ax1.set_zlim(-1, 8)
        ax1.set_box_aspect([1, 1, 0.6])

        # --- Subplot Derecho: Modificado ---
        ax2 = fig.add_subplot(122, projection='3d')
        ax2.plot_surface(X_mod, Y_mod, Z_mod, cmap='viridis', alpha=0.5, edgecolor='none')
        
        ax2.plot_wireframe(X_ctrl, Y_ctrl, Z_ctrl_mod, color='black', alpha=0.5, linestyle='--')
        
        colores = ['red' if (i==3 and j==3) else 'gray' for i in range(n+1) for j in range(m+1)]
        tamanos = [100 if (i==3 and j==3) else 20 for i in range(n+1) for j in range(m+1)]
        
        ax2.scatter(X_ctrl.flatten(), Y_ctrl.flatten(), Z_ctrl_mod.flatten(), color=colores, s=tamanos, zorder=10)

        ax2.set_title('Modificación de un punto de la malla', fontsize=12)
        ax2.set_xlabel('x')
        ax2.set_ylabel('y')
        ax2.set_zlabel('z')
        ax2.set_zlim(-1, 8)
        ax2.set_box_aspect([1, 1, 0.6])

        plt.tight_layout()
        plt.show()

    except Exception as e:
        print("\n[ERROR] Hubo un problema durante la ejecución:", e)