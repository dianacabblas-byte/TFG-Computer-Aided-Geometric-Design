import matplotlib.pyplot as plt
import numpy as np

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
    
    # Definición de los parámetros de la curva
    n = 7           # 8 funciones de base (N_0 a N_7)
    grado = 3       # p = 3 (cúbico)
    a = 0.0
    b = 5.0
    
    numero_puntos = 500

    try:
        # 1. Truco matemático: Usar la matriz identidad para aislar las funciones base
        puntos_control = np.eye(n + 1) 
        
        # 2. Cálculos previos
        U = calcular_vector_nodos(n, grado, a, b)
        vect_evaluacion = np.linspace(a, b, numero_puntos + 1)
        
        # Matriz para guardar los resultados
        funciones_base = np.zeros((len(vect_evaluacion), n + 1))
        
        # 3. Evaluación de la curva
        for i, u_val in enumerate(vect_evaluacion):
            funciones_base[i] = cox_de_boor(u_val, puntos_control, U, grado)
            
        print("Generando gráfica...")
        
        # 4. Gráfica
        plt.figure(figsize=(10, 5))
        
        # Pintamos cada función base
        for i in range(n + 1):
            plt.plot(vect_evaluacion, funciones_base[:, i], label=f'$N_{{{i},{grado}}}$', linewidth=1.5, zorder=3)
            
            # Solo sombreamos el soporte (área bajo la curva) de N_3,3 en gris claro
            if i == 3:
                plt.fill_between(vect_evaluacion, funciones_base[:, i], color='lightgray', alpha=0.6, zorder=2)
                
        # Líneas verticales en los nodos internos
        nodos_internos = set(U) - {a, b} 
        for nodo in nodos_internos:
            plt.axvline(x=nodo, color='gray', linestyle='--', alpha=0.4, zorder=1)
            
        # Texto del soporte de N_3,3 (que va de 0 a 4, centro en x=2)
        plt.text(2.0, 0.92, 'Soporte de $N_{3,3}$', ha='center', va='center', fontsize=11, 
                 bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'), zorder=4)
                 
        plt.title(f'Funciones base B-spline de grado {grado}')
        plt.xlabel('$u$')
        plt.ylabel(f'$N_{{i,{grado}}}(u)$')
        
        # Leyenda en dos columnas
        plt.legend(loc='upper right', ncol=2, fontsize=10)
        
        # Cuadrícula
        plt.grid(color='gray', linestyle='--', linewidth=0.2, alpha=0.5)
        
        plt.xlim(a, b)
        plt.ylim(0, 1.05)
    
        plt.show()
        
    except Exception as e:
        print("\n[ERROR] Hubo un problema durante la ejecución:", e)