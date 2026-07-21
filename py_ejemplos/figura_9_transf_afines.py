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
        # 1. Definición de parámetros base
        grado = 3
        a, b = 0.0, 1.0
        numero_puntos = 150
        vect_evaluacion = np.linspace(a, b, numero_puntos)
     
        # --- Curva original ---
        puntos_originales = np.array([
            [0, 0], 
            [1, 3], 
            [3, 4], 
            [5, 2], 
            [6, 0], 
            [4,-2]
            ])
        n_orig = len(puntos_originales) - 1
        U_orig = calcular_vector_nodos(n_orig, grado, a, b)
        
        # Evaluamos la curva
        curva_orig = np.zeros((numero_puntos, 2))
        for i, u_val in enumerate(vect_evaluacion):
            curva_orig[i] = cox_de_boor(u_val, puntos_originales, U_orig, grado)
        
        # 2. Aplicamos la transformación afín: Rotación 45º y traslación [5, 3]
        theta = np.radians(45)
        matriz_rot = np.array([
            [np.cos(theta), -np.sin(theta)], 
            [np.sin(theta),  np.cos(theta)]
        ])
        vector_trasl = np.array([5, 3])
        
        # Transformamos los puntos de control
        puntos_trans = np.dot(puntos_originales, matriz_rot.T) + vector_trasl
        
        # Reconstruimos la curva usando los nuevos puntos
        curva_trans = np.zeros((numero_puntos, 2))
        for i, u_val in enumerate(vect_evaluacion):
            curva_trans[i] = cox_de_boor(u_val, puntos_trans, U_orig, grado)
            
        print("Generando gráfica...")
        
        # 3. Dibujo
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Panel 1: Gráfica Original
        ax1.plot(puntos_originales[:, 0], puntos_originales[:, 1], 'o--', color='silver', label='Polígono Original $P_i$', zorder=2)
        ax1.plot(curva_orig[:, 0], curva_orig[:, 1], '-', color='royalblue', linewidth=1.5, label='Curva Original $Q(u)$', zorder=3)
        ax1.set_title('Configuración Original', fontsize=12)
        ax1.grid(color='gray', linestyle='--', linewidth=0.2, alpha=0.5)
        
        # Panel 2: Gráfica Transformada
        ax2.plot(puntos_trans[:, 0], puntos_trans[:, 1], 's--', color='orange', label=r'Polígono Transformado $\Phi(P_i)$', zorder=2)
        ax2.plot(curva_trans[:, 0], curva_trans[:, 1], '-', color='crimson', linewidth=1.5, label=r'Curva Transformada $\Phi(Q(u))$', zorder=3)
        ax2.set_title(r'Traslación $v=(5,3)$ y Rotación (45º)', fontsize=12)
        ax2.grid(color='gray', linestyle='--', linewidth=0.2, alpha=0.5)
        
       
        # 4. Forzamos el mismo tamaño visual
        # Centro de la caja original y transformada
        cx_orig = (np.max(puntos_originales[:, 0]) + np.min(puntos_originales[:, 0])) / 2.0
        cy_orig = (np.max(puntos_originales[:, 1]) + np.min(puntos_originales[:, 1])) / 2.0
        cx_trans = (np.max(puntos_trans[:, 0]) + np.min(puntos_trans[:, 0])) / 2.0
        cy_trans = (np.max(puntos_trans[:, 1]) + np.min(puntos_trans[:, 1])) / 2.0
        
        span_x = max(np.ptp(puntos_originales[:, 0]), np.ptp(puntos_trans[:, 0])) / 2.0
        span_y = max(np.ptp(puntos_originales[:, 1]), np.ptp(puntos_trans[:, 1])) / 2.0
        span_max = max(span_x, span_y)
        margen = 0.8 
        
        # Unidades extra en blanco
        espacio_leyenda = 3.0 
        
        # PANEL 1 (Azul): Espacio en blanco a la izquierda (- espacio_leyenda)
        ax1.set_xlim(cx_orig - span_max - margen - espacio_leyenda, cx_orig + span_max + margen)
        ax1.set_ylim(cy_orig - span_max - margen, cy_orig + span_max + margen)
        ax1.set_aspect('equal') 
        
        # PANEL 2 (Roja): Espacio en blanco a la derecha (+ espacio_leyenda)
        ax2.set_xlim(cx_trans - span_max - margen, cx_trans + span_max + margen + espacio_leyenda)
        ax2.set_ylim(cy_trans - span_max - margen, cy_trans + span_max + margen)
        ax2.set_aspect('equal')
        
        # Leyendas
        ax1.legend(loc='lower left', fontsize=8)
        ax2.legend(loc='lower right', fontsize=8)
        
        plt.tight_layout()
        plt.show()
        
    except Exception as e:
        print("\n[ERROR] Hubo un problema durante la ejecución:", e)