import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline


def interpolacion_exacta(x_datos, y_datos, x_eval):
    """Calcula una interpolación exacta (spline cúbico) que pasa por todos los puntos."""
    spline = CubicSpline(x_datos, y_datos)
    return spline(x_eval)

def aproximacion_minimos_cuadrados(x_datos, y_datos, x_eval, grado=3):
    """Calcula una aproximación polinómica usando mínimos cuadrados."""
    coeficientes = np.polyfit(x_datos, y_datos, grado)
    polinomio = np.poly1d(coeficientes)
    return polinomio(x_eval)

if __name__ == "__main__":
    # Definición del intervalo [a,b] para el eje X
    a = -3.0
    b = 3.0
    
    # Número de puntos sobre los que vamos a evaluar las curvas
    numero_puntos = 500
    
    try:
        # 1. Definición de los datos discretos originales 
        x_datos = np.linspace(a, b, 11)
       
        y_datos = np.array([-0.25, -0.5, -0.9, -1.1, -0.6, 0.2, 0.55, 0.75, 1.05, 0.85, 0.05])
        
        # 2. Cálculos previos (vector de evaluación)
        vect_evaluacion = np.linspace(a, b, numero_puntos + 1)
        
        # 3. Evaluación de las curvas
        y_interp = interpolacion_exacta(x_datos, y_datos, vect_evaluacion)
        y_aprox = aproximacion_minimos_cuadrados(x_datos, y_datos, vect_evaluacion, grado=3)
        
        print("Generando gráfica...")
        
        # 4. Gráfica
        plt.figure(figsize=(8, 5)) # Un poco más ancha que alta, igual que tu foto
        
        # Pintamos los datos discretos originales
        plt.plot(x_datos, y_datos, color='black', marker='o', linestyle='none', label='Datos discretos', zorder=5)
        
        # Pintamos la interpolación exacta
        plt.plot(vect_evaluacion, y_interp, color='tab:orange', linewidth=1.5, label='Interpolación exacta')
        
        # Pintamos la aproximación por mínimos cuadrado
        plt.plot(vect_evaluacion, y_aprox, color='tab:green', linewidth=1.5, label='Aproximación suave (mínimos cuadrados)')
        
        plt.title('Comparación conceptual: interpolación frente a aproximación')
        plt.xlabel('$x$') 
        plt.ylabel('$y$')
        
        # Leyenda arriba a la izquierda
        plt.legend(loc='upper left', fontsize=10)
        
        # Añadimos una cuadrícula
        plt.grid(color='gray', linestyle='-', linewidth=0.2, alpha=0.5)
        
        plt.show()
        
    except Exception as e:
        print("\n[ERROR] Hubo un problema durante la ejecución:", e)
