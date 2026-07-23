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
        red_control = np.array([
            [[ 0,  0, 0], [ 2,  0, 3], [ 2,  0, 6], [ 0,  0, 9]], 
            [[ 2,  4, 0], [ 4,  4, 3], [ 4,  4, 6], [ 2,  4, 9]], 
            [[ 4, -4, 0], [ 6, -4, 3], [ 6, -4, 6], [ 4, -4, 9]], 
            [[ 6,  4, 0], [ 8,  4, 3], [ 8,  4, 6], [ 6,  4, 9]], 
            [[ 8, -4, 0], [10, -4, 3], [10, -4, 6], [ 8, -4, 9]], 
            [[10,  0, 0], [12,  0, 3], [12,  0, 6], [10,  0, 9]]  
        ])
        
        grado_u = 3 
        grado_v = 2 
        
        
        # 2. Cálculo de la superfície
        res_u, res_v = 40, 40
        X_surf, Y_surf, Z_surf = evaluar_superficie(red_control, grado_u, grado_v, res_u, res_v)

        print("Generando gráfica...")
        
        
        # 3. Dibujo
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        
        # Superficie 
        ax.plot_surface(X_surf, Y_surf, Z_surf, cmap='viridis', alpha=0.7, edgecolor='none')


        X_ctrl = red_control[:, :, 0]
        Y_ctrl = red_control[:, :, 1]
        Z_ctrl = red_control[:, :, 2]
        
        # Red de control
        for j in range(len(red_control[0])):
            P_curva = red_control[:, j, :]
            if j == 0:
                ax.plot(P_curva[:, 0], P_curva[:, 1], P_curva[:, 2], color='gray', linestyle='--', marker='o', markersize=4, label='Red de Control')
            else:
                ax.plot(P_curva[:, 0], P_curva[:, 1], P_curva[:, 2], color='gray', linestyle='--', marker='o', markersize=4)
                
        for i in range(len(red_control)):
            P_linea = red_control[i, :, :]
            ax.plot(P_linea[:, 0], P_linea[:, 1], P_linea[:, 2], color='gray', linestyle='--')
            
        # Curvas generatrices
        U_u = calcular_vector_nodos(len(red_control) - 1, grado_u, 0.0, 1.0)
        vect_evaluacion = np.linspace(0.0, 1.0, 50)
        
        for j in range(len(red_control[0])):
            P_curva = red_control[:, j, :]
            curva_puntos = np.zeros((len(vect_evaluacion), 3))
            for k, u_val in enumerate(vect_evaluacion):
                curva_puntos[k] = cox_de_boor(u_val, P_curva, U_u, grado_u)
            
            if j == 0:
                ax.plot(curva_puntos[:, 0], curva_puntos[:, 1], curva_puntos[:, 2], color='crimson', linewidth=1.5, zorder=5, label='Curvas B-spline')
            else:
                ax.plot(curva_puntos[:, 0], curva_puntos[:, 1], curva_puntos[:, 2], color='crimson', linewidth=1.5, zorder=5)

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        
        # Proporción de la caja 3D
        ax.set_box_aspect([1, 1, 0.8])
        
        # Leyenda
        ax.legend(loc='upper right', fontsize=8)

        plt.tight_layout()
        plt.show()

    except Exception as e:
        print("\n[ERROR] Hubo un problema durante la ejecución:", e)