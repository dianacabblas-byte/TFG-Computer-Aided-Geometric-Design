import matplotlib.pyplot as plt
import numpy as np

def chaikin_corte_abierto(puntos):
    """Realiza una iteración del algoritmo de Chaikin (Corte de Esquinas)."""
    nuevos_puntos = [puntos[0]] # Mantenemos fijo el primer punto
    
    for i in range(len(puntos) - 1):
        p0 = puntos[i]
        p1 = puntos[i+1]
        
        # Generamos los nuevos puntos a 1/4 y 3/4 del segmento
        q0 = 0.75 * p0 + 0.25 * p1
        q1 = 0.25 * p0 + 0.75 * p1
        
        nuevos_puntos.extend([q0, q1])
        
    nuevos_puntos.append(puntos[-1]) # Mantenemos fijo el último punto
    return np.array(nuevos_puntos)

if __name__ == "__main__":
    
    try:
        # 1. Polígono de control original
        puntos_orig = np.array([
            [0.0, 0.0],
            [2.0, 7.0],
            [5.0, 9.0],
            [7.5, 2.0],
            [10.0, 6.0]
        ])
        
        # 2. Primer corte
        puntos_corte_1 = chaikin_corte_abierto(puntos_orig)
        
        # 3. Curva final (Aplicamos el algoritmo en bucle 6 veces)
        curva_final = puntos_orig.copy()
        for _ in range(6):
            curva_final = chaikin_corte_abierto(curva_final)
            
        print("Generando gráfica...")
        
        # 4. Dibujo de la gráfica
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # 4.1 Áreas recortadas (Triángulos de las esquinas)
        # Extraemos los puntos para sombrearlos
        triangulo_1 = np.array([puntos_corte_1[2], puntos_orig[1], puntos_corte_1[3]])
        triangulo_2 = np.array([puntos_corte_1[4], puntos_orig[2], puntos_corte_1[5]])
        triangulo_3 = np.array([puntos_corte_1[6], puntos_orig[3], puntos_corte_1[7]])
        
        # Dibujamos el primer triángulo con la etiqueta para la leyenda
        ax.fill(triangulo_1[:, 0], triangulo_1[:, 1], color='lightgray', alpha=0.7, zorder=1, label='Área "recortada"')
        # Dibujamos los demás sin etiqueta
        ax.fill(triangulo_2[:, 0], triangulo_2[:, 1], color='lightgray', alpha=0.7, zorder=1)
        ax.fill(triangulo_3[:, 0], triangulo_3[:, 1], color='lightgray', alpha=0.7, zorder=1)
        
        # Polígono original
        ax.plot(puntos_orig[:, 0], puntos_orig[:, 1], '--o', color='gray', markersize=5, zorder=2, label='Polígono original')
        
        # Polígono tras el primer corte
        ax.plot(puntos_corte_1[:, 0], puntos_corte_1[:, 1], '-o', color='royalblue', markersize=5, linewidth=1.5, zorder=3, label='Polígono tras el primer corte')
        
        # Curva final
        ax.plot(curva_final[:, 0], curva_final[:, 1], '-', color='crimson', linewidth=2.0, zorder=4, label='Curva final')
        
        ax.set_xlabel('x', fontsize=12)
        ax.set_ylabel('y', fontsize=12)
        ax.set_aspect('equal')
        
        # Leyenda
        ax.legend(loc='lower right', fontsize=9)
        
        plt.tight_layout()
        plt.show()
        
    except Exception as e:
        print("\n[ERROR] Hubo un problema durante la ejecución:", e)