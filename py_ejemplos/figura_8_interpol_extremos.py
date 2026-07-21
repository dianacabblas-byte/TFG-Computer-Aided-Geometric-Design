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
        # 1. Definición de parámetros generales
        grado = 3
        a, b = 0, 1
        numero_puntos = 150
        vect_evaluacion = np.linspace(a, b, numero_puntos)
        
        print("Calculando curva abierta...")
        
        # --- Curva 1: Abierta  ---
        puntos_abierta = np.array([
            [0, 0], 
            [2, 5], 
            [4, -5], 
            [6, 5], 
            [8, -5], 
            [10, 0]
        ])
        n_abierta = len(puntos_abierta) - 1
        U_abierta = calcular_vector_nodos(n_abierta, grado, a, b)
        
        curva_abierta = np.zeros((numero_puntos, 2))
        for i, u_val in enumerate(vect_evaluacion):
            curva_abierta[i] = cox_de_boor(u_val, puntos_abierta, U_abierta, grado)

        print("Calculando curva cerrada...")
        
        # --- Curva 2: Cerrada ---
        puntos_cerrada = np.array([
            [ 2, 6],  
            [ 5, 9], 
            [ 8, 6], 
            [10, 8],
            [12, 6],  
            [10, 4], 
            [ 8, 6],
            [ 5, 3], 
            [ 2, 6]
        ])
        n_cerrada = len(puntos_cerrada) - 1
        U_cerrada = calcular_vector_nodos(n_cerrada, grado, a, b)
        
        curva_cerrada = np.zeros((numero_puntos, 2))
        for i, u_val in enumerate(vect_evaluacion):
            curva_cerrada[i] = cox_de_boor(u_val, puntos_cerrada, U_cerrada, grado)

        print("Generando gráfica...")
        
        # 2. Dibujo
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Panel 1: Curva Abierta
        ax1.plot(puntos_abierta[:, 0], puntos_abierta[:, 1], 'o--', color='silver', label='Polígono de control')
        ax1.plot(curva_abierta[:, 0], curva_abierta[:, 1], '-', color='royalblue', linewidth=1.5, zorder=3, label='Curva')
        
        # Resaltar extremos con puntos rojos
        ax1.plot(puntos_abierta[0, 0], puntos_abierta[0, 1], 'o', color='crimson', markersize=6, zorder=5, label='Extremos interpolados')
        ax1.plot(puntos_abierta[-1, 0], puntos_abierta[-1, 1], 'o', color='crimson', markersize=6, zorder=5)
        
        ax1.set_title('Curva B-spline abierta', fontsize=12)
        ax1.grid(color='gray', linestyle='--', linewidth=0.2, alpha=0.5)
        
        # Forzamos los límites para crear espacio vacío abajo (-8)
        ax1.set_xlim(-3, 13)
        ax1.set_ylim(-8, 6)
        ax1.set_aspect('equal')
        ax1.legend(loc='lower left', fontsize=8)
        
        # Panel 2: Curva Cerrada 
        ax2.plot(puntos_cerrada[:, 0], puntos_cerrada[:, 1], 'o--', color='silver', label='Polígono de Control')
        ax2.plot(curva_cerrada[:, 0], curva_cerrada[:, 1], '-', color='mediumseagreen', linewidth=1.5, zorder=3, label='Curva')
        
        # Resaltar cierre con punto rojo
        ax2.plot(puntos_cerrada[0, 0], puntos_cerrada[0, 1], 'o', color='crimson', markersize=6, zorder=5, label='Punto de cierre ($P_0 = P_8$)')
        
        ax2.set_title('Curva B-spline cerrada', fontsize=12)
        ax2.grid(color='gray', linestyle='--', linewidth=0.2, alpha=0.5)
        
        # Aplicamos el ZOOM
        ax2.set_xlim(-6, 6)       
        ax2.set_ylim(-2.5, 8.0)  
        ax2.set_aspect('equal')
        ax2.legend(loc='lower right', fontsize=8)
        
        plt.tight_layout()
        plt.show()
        
    except Exception as e:
        print("\n[ERROR] Hubo un problema durante la ejecución:", e)