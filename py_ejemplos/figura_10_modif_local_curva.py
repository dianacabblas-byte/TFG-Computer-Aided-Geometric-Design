import numpy as np
import matplotlib.pyplot as plt

def calcular_vector_nodos(n, grado, a, b):
    ptos_internos = n - grado
    h = (b - a) / (ptos_internos + 1)
    
    inicio = [a] * (grado + 1)
    medio = [a + i * h for i in range(1, ptos_internos + 1)]
    fin = [b] * (grado + 1)
    
    return inicio + medio + fin

def cox_de_boor(u, puntos_control, U, grado):
    n = len(puntos_control) - 1
    r = grado
    
    # Buscamos r tal que u esté en el intervalo [U_r, U_(r+1)]
    while r <= n:
        if u < U[r + 1] or r == n: 
            break
        r += 1
        
    # Matriz auxiliar para calcular la evaluación
    Q = np.array(puntos_control, dtype=float)
    
    # Aplicamos Cox-de Boor
    for j in range(1, grado + 1):
        for i in range(r - grado, r - j + 1):
            
            # Denominador de las formulas
            dif = U[i + grado + 1] - U[i + j]
            if abs(dif) < 1e-9:
                Q[i, :] = 0.0
            else:
                factor_1 = (U[i + grado + 1] - u) / dif
                factor_2 = (u - U[i + j]) / dif
                
                Q[i, :] = factor_1 * Q[i, :] + factor_2 * Q[i+1, :]
                
    return Q[r - grado, :]


if __name__ == "__main__":
    
    try:
        # 1. Definición de los parámetros base
        grado = 3
        a, b = 0, 1
        numero_puntos = 300
        vect_evaluacion = np.linspace(a, b, numero_puntos)
        
        # Puntos de control 
        puntos_orig = np.array([
            [-2.5, -1.0],
            [-2.0, -1.2],
            [-1.5, -1.2],
            [-1.0, -0.9],
            [-0.5, -0.5],
            [ 0.0,  0.0],
            [ 0.5,  0.5],
            [ 1.0,  0.9],
            [ 1.5,  1.2],
            [ 2.0,  1.2],
            [ 2.5,  1.0]
        ])
        
        # Copiamos y modificamos P5
        puntos_mod = np.copy(puntos_orig)
        indice_modificado = 5
        puntos_mod[indice_modificado, 1] = 0.8  # Lo subimos hacia arriba
        
        n_puntos = len(puntos_orig) - 1
        U = calcular_vector_nodos(n_puntos, grado, a, b)
        
        curva_orig = np.zeros((numero_puntos, 2))
        for i, u_val in enumerate(vect_evaluacion):
            curva_orig[i] = cox_de_boor(u_val, puntos_orig, U, grado)
            
        curva_mod = np.zeros((numero_puntos, 2))
        for i, u_val in enumerate(vect_evaluacion):
            curva_mod[i] = cox_de_boor(u_val, puntos_mod, U, grado)
            
        print("Generando gráfica...")
        
        # 2. Dibujo
        fig, ax = plt.subplots(figsize=(10, 5))
        
        # Las curvas principales
        ax.plot(curva_orig[:, 0], curva_orig[:, 1], '-', color='tab:blue', linewidth=2.5, alpha=0.9, label='Curva B-spline original', zorder=3)
        ax.plot(curva_mod[:, 0], curva_mod[:, 1], '--', color='tab:red', linewidth=2.5, alpha=0.9, label='Curva tras modificar un punto de control', zorder=4)
        
        # Polígonos de control
        ax.plot(puntos_orig[:, 0], puntos_orig[:, 1], '--', color='silver', alpha=0.5, zorder=1)
        ax.plot(puntos_mod[:, 0], puntos_mod[:, 1], '--', color='lightcoral', alpha=0.3, zorder=1)
        
        # Puntos estáticos
        puntos_estaticos = np.delete(puntos_orig, indice_modificado, axis=0)
        ax.plot(puntos_estaticos[:, 0], puntos_estaticos[:, 1], 'o', color='silver', markersize=4, zorder=2)
        
        # Puntos modificados resaltados
        ax.plot(puntos_orig[indice_modificado, 0], puntos_orig[indice_modificado, 1], 'o', color='tab:blue', markeredgecolor='black', markersize=6, zorder=5)
        ax.plot(puntos_mod[indice_modificado, 0], puntos_mod[indice_modificado, 1], 'o', color='tab:red', markeredgecolor='black', markersize=6, zorder=5)
        
        ax.set_title('Efecto local de la modificación de un punto de control en B-splines', fontsize=12)
        ax.set_xlabel('x') 
        ax.set_ylabel('y')
        
        ax.legend(loc='upper left', fontsize=8)
        
        # Cuadrícula
        ax.grid(color='gray', linestyle='--', linewidth=0.2, alpha=0.5)
        
        plt.tight_layout()
        plt.show()
        
    except Exception as e:
        print("\n[ERROR] Hubo un problema durante la ejecución:", e)