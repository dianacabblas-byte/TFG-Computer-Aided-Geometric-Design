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
        a, b = 0, 1
        numero_puntos = 500
        
        puntos = np.array([
            [0, 0],  
            [2, 8],  
            [4, 1],  
            [6, 8],  
            [8, 0]
        ])

        n_puntos = len(puntos) - 1
        U = calcular_vector_nodos(n_puntos, grado, a, b)
        
        vect_evaluacion = np.linspace(a, b, numero_puntos)
        curva = np.zeros((numero_puntos, 2))
        for i, u_val in enumerate(vect_evaluacion):
            curva[i] = cox_de_boor(u_val, puntos, U, grado)

        # 2. Cálculo exacto de intersecciones
        y_linea = 4.0
        
        intersecciones_poligono = []
        for i in range(len(puntos) - 1):
            x1, y1 = puntos[i]
            x2, y2 = puntos[i+1]
            # Detectar cruce por cambio de signo (o contacto exacto)
            if (y1 - y_linea) * (y2 - y_linea) < 0 or y1 == y_linea:
                x_int = x1 + (y_linea - y1) * (x2 - x1) / (y2 - y1)
                intersecciones_poligono.append(x_int)

        intersecciones_curva = []
        for i in range(len(curva) - 1):
            x1, y1 = curva[i]
            x2, y2 = curva[i+1]
            # Detectar cruce en la curva evaluada
            if (y1 - y_linea) * (y2 - y_linea) < 0:
                x_int = x1 + (y_linea - y1) * (x2 - x1) / (y2 - y1)
                intersecciones_curva.append(x_int)

        print("Generando gráfica...")
        
        # 3. Dibujo
        fig, ax = plt.subplots(figsize=(14, 6))

        # Dibujo de la recta secante
        ax.axhline(y=y_linea, color='gray', linestyle=':', linewidth=1.5, zorder=1)

        # Polígono y Curva
        ax.plot(puntos[:, 0], puntos[:, 1], '--', color='silver', linewidth=1.5, label='Polígono de control', zorder=2)
        ax.plot(puntos[:, 0], puntos[:, 1], 'o', color='silver', markersize=5, zorder=2)
        ax.plot(curva[:, 0], curva[:, 1], '-', color='blue', linewidth=1.5, alpha=0.9, label='Curva B-spline', zorder=3)

        # Puntos de intersección
        ax.plot(intersecciones_poligono, [y_linea]*len(intersecciones_poligono), 's', 
                color='silver', markeredgecolor='gray', zorder=4,
                label=f'Cortes en polígono ($N={len(intersecciones_poligono)}$)')

        ax.plot(intersecciones_curva, [y_linea]*len(intersecciones_curva), 'o', 
                color='royalblue', markeredgecolor='black', zorder=5, 
                label=f'Cortes en curva ($N={len(intersecciones_curva)}$)')

        ax.set_title('Propiedad de disminución de la variación', fontsize=12)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_aspect('equal')

        # Márgenes amplios
        ax.set_xlim(-3, 13) 
        ax.set_ylim(-1, 9)

        # Leyenda
        ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
        
        # Cuadrícula
        ax.grid(color='gray', linestyle='--', linewidth=0.2, alpha=0.5)

        plt.tight_layout()
        plt.show()

    except Exception as e:
        print("\n[ERROR] Hubo un problema durante la ejecución:", e)