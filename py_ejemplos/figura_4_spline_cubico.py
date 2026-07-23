import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

def funcion_runge(x):
    return 1.0 / (1.0 + 25.0 * x**2)

def interpolacion_spline(x_nodos, y_nodos, x_eval):
    spline = CubicSpline(x_nodos, y_nodos)
    return spline(x_eval)

if __name__ == "__main__":
    # Definición del intervalo [a,b] para el eje X
    a = -1.0
    b = 1.0
    
    # Número de nodos de interpolación
    numero_nodos = 12
    
    # Número de puntos para dibujar la curva
    numero_puntos = 500
    
    try:
        # 1. Definición de los nodos discretos
        x_nodos = np.linspace(a, b, numero_nodos)
        y_nodos = funcion_runge(x_nodos)
        
        # 2. Cálculos previos (vector de evaluación continuo)
        vect_evaluacion = np.linspace(a, b, numero_puntos + 1)
        
        # 3. Evaluación de la curva Spline
        y_spline = interpolacion_spline(x_nodos, y_nodos, vect_evaluacion)
        
        print("Generando gráfica...")
        
        # 4. Gráfica
        plt.figure(figsize=(8, 5))
        
        # Pintamos el Spline
        plt.plot(vect_evaluacion, y_spline, color='tab:blue', linewidth=1.5, label='Spline cúbico interpolante', zorder=2)
        
        # Pintamos los nodos discretos
        plt.plot(x_nodos, y_nodos, color='black', marker='o', linestyle='none', label='Datos interpolados', zorder=5)
        
        plt.xlabel('$x$') 
        plt.ylabel('$y$')
        
        # Leyenda arriba a la derecha, como en tu imagen
        plt.legend(loc='upper right', fontsize=8)
        
        # Añadimos una cuadrícula suave y profesional
        plt.grid(color='gray', linestyle='-', linewidth=0.2, alpha=0.5)
        
        plt.show()
        
    except Exception as e:
        print("\n[ERROR] Hubo un problema durante la ejecución:", e)