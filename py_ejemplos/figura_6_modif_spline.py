import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

def interpolacion_spline(x_nodos, y_nodos, x_eval):
    """Calcula la interpolación usando un Spline Cúbico."""
    spline = CubicSpline(x_nodos, y_nodos)
    return spline(x_eval)

if __name__ == "__main__":
    # Definición del intervalo [a,b] para el eje X
    a = -2.5
    b = 2.5
    
    # Número de nodos de interpolación
    numero_nodos = 11
    
    # Número de puntos para dibujar las curvas
    numero_puntos = 500
    
    try:
        # 1. Definición de los nodos discretos originales
        x_nodos = np.linspace(a, b, numero_nodos)
        
        # Datos para generar la curva
        y_original = np.array([-1.0, -1.2, -1.2, -1.0, -0.55, 0.0, 0.55, 1.0, 1.2, 1.2, 1.0])
        
        # 2. Modificamos un único nodo (el central en x = 0.0)
        y_modificado = y_original.copy()
        indice_modificado = 5 # El nodo central
        y_modificado[indice_modificado] = 0.55 # Subimos el punto
        
        # 3. Cálculos previos (vector de evaluación continuo)
        vect_evaluacion = np.linspace(a, b, numero_puntos + 1)
        
        # 4. Evaluación de las dos curvas competidoras
        y_spline_orig = interpolacion_spline(x_nodos, y_original, vect_evaluacion)
        y_spline_mod = interpolacion_spline(x_nodos, y_modificado, vect_evaluacion)
        
        print("Generando gráfica...")
        
        # 5. Gráfica
        plt.figure(figsize=(8, 5))
        
        # Pintamos el Spline original (Línea azul sólida)
        plt.plot(vect_evaluacion, y_spline_orig, color='tab:blue', linewidth=2, alpha=0.9, label='Spline original', zorder=2)
        
        # Pintamos el Spline modificado (Línea roja discontinua)
        plt.plot(vect_evaluacion, y_spline_mod, color='tab:red', linestyle='--', linewidth=2, alpha=0.9, label='Spline tras modificar un nodo', zorder=3)
        
        # Pintamos los puntos base (nodos que no han cambiado)
        plt.plot(x_nodos, y_original, color='tab:red', marker='o', markersize=4, linestyle='none', zorder=4)
        
        # Resaltamos especialmente el punto original antes de moverlo
        plt.plot(x_nodos[indice_modificado], y_original[indice_modificado], color='tab:blue', marker='o', markeredgecolor='black', zorder=5)
        
        # Resaltamos especialmente el nuevo punto movido
        plt.plot(x_nodos[indice_modificado], y_modificado[indice_modificado], color='tab:red', marker='o', markeredgecolor='black', zorder=5)
        
        plt.title('Efecto global de la modificación de un dato en el spline cúbico')
        plt.xlabel('x') 
        plt.ylabel('y')
        
        # Leyenda arriba a la izquierda
        plt.legend(loc='upper left', fontsize=8)
        
        # Añadimos la cuadrícula 
        plt.grid(color='gray', linestyle='--', linewidth=0.2, alpha=0.5)
        
        plt.show()
        
    except Exception as e:
        print("\n[ERROR] Hubo un problema durante la ejecución:", e)