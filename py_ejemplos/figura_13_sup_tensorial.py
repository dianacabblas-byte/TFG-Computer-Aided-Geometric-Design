import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    
    try:
        
        # 1. Puntos de control
        # Matriz 2x2 para P_00, P_01, P_10, P_11
        P = np.zeros((2, 2, 3))
        P[0, 0] = [0, 0, 0]
        P[0, 1] = [0, 1, 2] 
        P[1, 0] = [1, 0, 2] 
        P[1, 1] = [1, 1, 0]

    
        # 2. Cálculo de la superfície
        res = 30
        s_vals = np.linspace(0, 1, res)
        t_vals = np.linspace(0, 1, res)

        X_sup = np.zeros((res, res))
        Y_sup = np.zeros((res, res))
        Z_sup = np.zeros((res, res))

        for i, s in enumerate(s_vals):
            for j, t in enumerate(t_vals):
                # Bases unidimensionales lineales
                f0, f1 = 1 - s, s
                g0, g1 = 1 - t, t
                
                # Producto tensorial directo
                punto = (f0 * g0) * P[0, 0] + (f0 * g1) * P[0, 1] + (f1 * g0) * P[1, 0] + (f1 * g1) * P[1, 1]
                
                X_sup[i, j] = punto[0]
                Y_sup[i, j] = punto[1]
                Z_sup[i, j] = punto[2]

        print("Generando gráfica...")
        
    
        # 3. Dibujo
        fig = plt.figure(figsize=(9, 7))
        ax = fig.add_subplot(111, projection='3d')

        # Superficie
        ax.plot_surface(X_sup, Y_sup, Z_sup, cmap='viridis', alpha=0.9, edgecolor='none')

        # Red de control 
        X_ctrl = P[:, :, 0]
        Y_ctrl = P[:, :, 1]
        Z_ctrl = P[:, :, 2]
        
        ax.plot_wireframe(X_ctrl, Y_ctrl, Z_ctrl, color='gray', linewidth=2, linestyle='--')

        # Puntos de control
        ax.scatter(X_ctrl, Y_ctrl, Z_ctrl, color='crimson', s=50, zorder=5, label='Puntos de Control ($P_{i,j}$)')

        ax.set_title('Superficie Tensorial Bilineal', fontsize=12)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        
        # Leyenda
        ax.legend(loc='upper right', fontsize=8)
        
        # Proporción de la caja 3D
        ax.set_box_aspect([1, 1, 0.8])

        plt.tight_layout()
        plt.show()

    except Exception as e:
        print("\n[ERROR] Hubo un problema durante la ejecución:", e)