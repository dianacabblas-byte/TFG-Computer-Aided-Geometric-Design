import matplotlib.pyplot as plt
import numpy as np

def calcular_vector_nodos(n, grado, a, b):
    num_internos = n - grado
    h = (b - a) / (num_internos + 1)
    
    # Repetimos el punto a (grado + 1) veces para que pase exactamente por él
    inicio = [a] * (grado + 1)
    
    # Puntos internos del intervalo
    medio = [a + i * h for i in range(1, num_internos + 1)]
    
    # Repetimos el punto b (grado + 1) veces para que pase exactamente por él
    fin = [b] * (grado + 1)
    
    return inicio + medio + fin # Devolvemos el vector ensamblado

def cox_de_boor(u, puntos_control, U, grado):
    n = len(puntos_control) - 1
    r = grado
    
    # Buscamos r tq el valor esté en el intervalo [u_r, u_(r+1)]
    while r <= n:
        if u < U[r + 1] or r == n:
            break
        r += 1
        
    # Matriz auxiliar para calcular la evaluación
    Q = np.array(puntos_control, dtype=float)

    # Aplicamos Cox-de Boor
    for j in range(1, grado + 1):
        for i in range(r - grado, r - j + 1):
            
            # Denominador de las formulas:
            dif = U[i + grado + 1] - U[i + j]
            if abs(dif) < 1e-9:  # si la diferencia es nula
                Q[i,:] = 0.0
            else:
                factor_1 = (U[i + grado + 1] - u) / dif
                factor_2 = (u - U[i + j]) / dif
                
                Q[i,:] = factor_1 * Q[i,:] + factor_2 * Q[i+1,:]

    return Q[r - grado]

def evaluar_superficie(red_control, grado_u, grado_v, num_u, num_v):
    """Aplica Cox-de Boor en las dos dimensiones (u y v)"""
    a, b = 0.0, 1.0
    c, d = 0.0, 1.0
    
    n = len(red_control) - 1      # Filas
    m = len(red_control[0]) - 1   # Columnas
    
    # Cálculos previos de nodos e intervalos
    U_u = calcular_vector_nodos(n, grado_u, a, b)
    V_v = calcular_vector_nodos(m, grado_v, c, d)
    
    tiempos_u = np.linspace(a, b, num_u + 1)
    tiempos_v = np.linspace(c, d, num_v + 1)
    
    # Pre-creamos las matrices para la superficie final
    X = np.zeros((num_u + 1, num_v + 1))
    Y = np.zeros((num_u + 1, num_v + 1))
    Z = np.zeros((num_u + 1, num_v + 1))
    
    # Matriz auxiliar
    curvas_v = np.zeros((num_v + 1, n + 1, 3))
    
    # PASO 1 (Dirección V): Evaluamos cada fila de la red para cada instante de V
    for j, val_v in enumerate(tiempos_v):
        for i in range(n + 1):
            curvas_v[j, i] = cox_de_boor(val_v, red_control[i], V_v, grado_v)
            
    # PASO 2 (Dirección U): Evaluamos las columnas precalculadas para cada instante de U
    for i, val_u in enumerate(tiempos_u):
        for j, val_v in enumerate(tiempos_v):
            
            # Cogemos la columna intermedia de la matriz auxiliar
            puntos_intermedios = curvas_v[j]
            
            # Aplicamos el paso final
            punto_final = cox_de_boor(val_u, puntos_intermedios, U_u, grado_u)
            
            # Asignamos el punto directamente a la celda de la matriz
            X[i, j] = punto_final[0]
            Y[i, j] = punto_final[1]
            Z[i, j] = punto_final[2]
        
    return X, Y, Z

def leer_red_archivo(archivo_path, dim):
    with open(archivo_path, 'r') as archivo:
        # Leemos la 'n' y la 'm', las dimensiones de la red 
        n = int(archivo.readline().strip()) - 1
        m = int(archivo.readline().strip()) - 1
        
        # Pre-creamos la matriz con el tamaño concreto
        # Dimensiones: (Filas) x (Columnas) x (Coordenadas)
        red_control = np.zeros((n + 1, m + 1, dim))
        
        # Leemos los puntos del archivo
        valores = list(map(float, archivo.read().split()))
        
        # Rellenamos la matriz 3D
        idx = 0 # Llevamos un contador de por qué punto de la lista larga vamos
        for i in range(n + 1):
            for j in range(m + 1):
                # Cortamos las coordenadas de este punto y las guardamos
                red_control[i, j, :] = valores[idx * dim : (idx + 1) * dim]
                idx += 1
            
    return red_control

if __name__ == "__main__":
    # Número de puntos sobre los que vamos a evaluar la superficie
    numero_puntos_u = 50
    numero_puntos_v = 50
    
    # Dimensión de los puntos
    dim_puntos = 3
    
    print("---CÁLCULO Y DIBUJO DE SUPERFICIES B-SPLINE 3D---")
    archivo_entrada = input("Introduzca el nombre del archivo con la red (ej. red.txt): ")
    
    try:
        # 1. Lectura de la red de control del archivo
        red_control = leer_red_archivo(archivo_entrada, dim_puntos)
        
        grado_u = int(input("Introduzca el grado en la dirección u (no puede superar n): "))
        grado_v = int(input("Introduzca el grado en la dirección v (no puede superar m): "))
        
        print("Calculando malla de la superficie")
        
        # 2. Evaluación de la superficie
        X, Y, Z = evaluar_superficie(red_control, grado_u, grado_v, numero_puntos_u, numero_puntos_v)
        
        print("¡Superficie calculada! Abriendo visor 3D...")
        
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Pintamos la superficie generada
        ax.plot_surface(X, Y, Z, cmap='viridis', alpha=0.8, edgecolor='none')
        
        # Extraemos las coordenadas de la red de control original
        X_red = red_control[:, :, 0]
        Y_red = red_control[:, :, 1]
        Z_red = red_control[:, :, 2]
        
        # Pintamos los puntos de control (rojos)
        ax.scatter(X_red, Y_red, Z_red, color='red', s=40, label='Puntos de Control')
        
        # Pintamos la línea de la red (malla de alambre gris)
        ax.plot_wireframe(X_red, Y_red, Z_red, color='gray', alpha=0.5, linestyle='--')

        ax.set_title('Superficie B-Spline')
        ax.set_xlabel('Eje X')
        ax.set_ylabel('Eje Y')
        ax.set_zlabel('Eje Z')
        ax.legend()
        
        # Para que los 3 ejes mantengan la proporción: 
        ax.set_box_aspect([1,1,1])
        
        # Mostramos la gráfica en pantalla
        plt.show()

    except FileNotFoundError:
        print(f"\n[ERROR] No se ha encontrado el archivo '{archivo_entrada}'.")
        print("Asegúrate de poner el nombre exacto o la ruta completa.")
    except Exception as e:
        print("\n[ERROR] Hubo un problema durante la ejecución:", e)