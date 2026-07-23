import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import BarycentricInterpolator

def funcion_runge(x):
    return 1.0 / (1.0 + 25.0 * x**2)

def interpolacion_polinomica(x_nodos, y_nodos, x_eval):
    """Calcula el polinomio interpolador global usando el método baricéntrico."""
    polinomio = BarycentricInterpolator(x_nodos, y_nodos)
    return polinomio(x_eval)

if __name__ == "__main__":
    # Definición del intervalo [a,b] para el eje X
    a = -1.0
    b = 1.0
    
    # Número de nodos de interpolación (16 nodos = polinomio de grado 15)
    numero_nodos = 16
    
    # Número de puntos sobre los que vamos a evaluar las curvas
    numero_puntos = 500
    
    try:
        # 1. Definición de los nodos
        x_nodos = np.linspace(a, b, numero_nodos)
        y_nodos = funcion_runge(x_nodos)
        
        # 2. Cálculos previos (vector de evaluación)
        vect_evaluacion = np.linspace(a, b, numero_puntos + 1)
        
        # 3. Evaluación de las curvas
        y_real = funcion_runge(vect_evaluacion)
        y_interp = interpolacion_polinomica(x_nodos, y_nodos, vect_evaluacion)
        
        print("Generando gráfica...")
        
        # 4. Gráfica
        plt.figure(figsize=(8, 5))
        
        # Pintamos la función original
        # Usamos r'' para que Matplotlib dibuje la fórmula matemática en formato LaTeX
        plt.plot(vect_evaluacion, y_real, color='black', linestyle='--', linewidth=1.5, label=r'$f(x) = \frac{1}{1 + 25x^2}$', zorder=2)
        
        # Pintamos la interpolación global
        plt.plot(vect_evaluacion, y_interp, color='red', linewidth=1.5, label=f'Interpolación polinómica (grado {numero_nodos - 1})', zorder=3)
        
        # Pintamos los nodos discretos
        plt.plot(x_nodos, y_nodos, color='tab:blue', marker='o', linestyle='none', label='Nodos equidistantes', zorder=4)
        
        plt.xlabel('x') 
        plt.ylabel('y')
        
        # Ajustamos un poco los límites del eje Y para que quede el mismo aire arriba y abajo que en tu foto
        plt.ylim(-0.9, 3.2)
        
        # Leyenda arriba en el centro
        plt.legend(loc='upper center', fontsize=8)
        
        # Añadimos una cuadrícula
        plt.grid(color='gray', linestyle='--', linewidth=0.3, alpha=0.7)
        
        plt.show()
        
    except Exception as e:
        print("\n[ERROR] Hubo un problema durante la ejecución:", e)