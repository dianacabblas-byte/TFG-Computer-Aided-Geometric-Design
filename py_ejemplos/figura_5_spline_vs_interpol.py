import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import BarycentricInterpolator, CubicSpline

def funcion_runge(x):
    return 1.0 / (1.0 + 25.0 * x**2)

def interpolacion_polinomica(x_nodos, y_nodos, x_eval):
    """Calcula el polinomio interpolador global usando el método baricéntrico."""
    polinomio = BarycentricInterpolator(x_nodos, y_nodos)
    return polinomio(x_eval)

def interpolacion_spline(x_nodos, y_nodos, x_eval):
    """Calcula la interpolación usando un Spline Cúbico."""
    spline = CubicSpline(x_nodos, y_nodos)
    return spline(x_eval)

if __name__ == "__main__":
    # Definición del intervalo [a,b] para el eje X
    a = -1.0
    b = 1.0
    
    # Número de nodos de interpolación
    numero_nodos = 16
    
    # Número de puntos para dibujar las curvas
    numero_puntos = 500
    
    try:
        # 1. Definición de los nodos discretos
        x_nodos = np.linspace(a, b, numero_nodos)
        y_nodos = funcion_runge(x_nodos)
        
        # 2. Cálculos previos (vector de evaluación continuo)
        vect_evaluacion = np.linspace(a, b, numero_puntos + 1)
        
        # 3. Evaluación de las dos curvas competidoras
        y_polinomio = interpolacion_polinomica(x_nodos, y_nodos, vect_evaluacion)
        y_spline = interpolacion_spline(x_nodos, y_nodos, vect_evaluacion)
        
        print("Generando gráfica...")
        
        # 4. Gráfica
        plt.figure(figsize=(8, 5))
        
        # Pintamos el Spline
        plt.plot(vect_evaluacion, y_spline, color='tab:blue', linewidth=2, alpha=0.9, label='Spline cúbico', zorder=2)
        
        # Pintamos el Polinomio
        plt.plot(vect_evaluacion, y_polinomio, color='red', linestyle='--', linewidth=1.5, alpha=0.7, label=f'Polinomio global (grado {numero_nodos - 1})', zorder=3)
        
        # Pintamos los datos
        plt.plot(x_nodos, y_nodos, color='black', marker='o', markersize=4, linestyle='none', label='Datos', zorder=5)
        
        plt.title('Spline cúbico frente a interpolación polinómica global')
        plt.xlabel('x') 
        plt.ylabel('y')
        
        # Ajustamos los límites de Y para que se vean las oscilaciones del polinomio
        plt.ylim(-0.9, 3.2)
        
        # Leyenda arriba en el centro
        plt.legend(loc='upper center', fontsize=8)
        
        # Añadimos la cuadrícula 
        plt.grid(color='gray', linestyle='--', linewidth=0.2, alpha=0.5)
        
        plt.show()
        
    except Exception as e:
        print("\n[ERROR] Hubo un problema durante la ejecución:", e)