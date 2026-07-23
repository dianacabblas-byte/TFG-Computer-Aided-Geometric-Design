import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splprep, splev

if __name__ == "__main__":
    
    try:
        # 1. Definimos los puntos de datos (D_i) en R^2
        x_datos = np.array([0.0, 1.5, 2.0, 1.0, 0.5, 1.5, 3.0])
        y_datos = np.array([0.0, 1.0, 3.0, 3.5, 2.0, 0.5, 0.0])
        
        
        # 2. Interpolación
        # splprep nos devuelve dos cosas:
        # - interp: Contiene (Vector de nodos U, Puntos de control P_i, Grado p)
        # - t: Los instantes de tiempo t_i en cada punto de dato
        interp, t = splprep([x_datos, y_datos], s=0, k=3)

        # 3. Aproximacion
        aprox, _ = splprep([x_datos, y_datos], s=1.5, k=3)

        # Generamos una partición del parámetro para trazar la curva
        res = 100
        t_final = np.linspace(0, 1, res)
        
        # Evaluamos las curvas B-spline usando la partición final
        x_interp, y_interp = splev(t_final, interp)
        x_aprox, y_aprox = splev(t_final, aprox)

        print("Generando gráfica...")
        
        # 4. Dibujo
        fig = plt.figure(figsize=(9, 7))
        ax = fig.add_subplot(111)

        # Puntos de datos
        ax.plot(x_datos, y_datos, 'ko', zorder=5, label='Datos discretos ($D_i$)')

        # Curvas
        ax.plot(x_interp, y_interp, '-', color='tab:orange', linewidth=1.5, label='Interpolación exacta')
        ax.plot(x_aprox, y_aprox, '-', color='tab:green', linewidth=1.5, label='Aproximación suave (mínimos cuadrados)')

        # Configuración de texto y ejes
        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$')
        
        # Estilo de la cuadrícula y proporción de los ejes
        ax.grid(True, linestyle='--', alpha=0.6)
        ax.axis('equal') 

        plt.tight_layout()
        plt.show()

    except Exception as e:
        print("\n[ERROR] Hubo un problema durante la ejecución:", e)