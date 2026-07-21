import matplotlib.pyplot as plt
import numpy as np

def calcular_vector_nodos(n, grado, a, b):
    ptos_internos = n - grado
    h = (b - a) / (ptos_internos + 1)
    
    # Repetimos el punto a (grado + 1) veces para que pase exactamente por él
    inicio = [a] * (grado + 1)
    
    # Puntos internos del intervalo
    medio = [a + i * h for i in range(1, ptos_internos + 1)]
    
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
    Q = np.array(puntos_control)

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

def leer_puntos_control(archivo_path, dim):
    with open(archivo_path, 'r') as archivo:
       # Leemos la 'n' (índice del número de puntos) de la primera línea
       n = int(archivo.readline().strip()) - 1
       
       # Pre-creamos la matriz con el tamaño concreto.
       # Dimensiones: (Cantidad de puntos) x (Coordenadas)
       puntos_control = np.zeros((n + 1, dim))
       
       # Leemos los puntos del archivo
       valores = list(map(float, archivo.read().split()))
       
       # Rellenamos la matriz 
       for i in range(n + 1):
           # Insertamos el trozo de números de este punto a su fila
           puntos_control[i,:] = valores[i * dim : (i + 1) * dim]
            
    return puntos_control

if __name__ == "__main__":
    # Definición del intervalo [a,b]
    a = 0.0 
    b = 1.0
    
    # Número de puntos sobre los que vamos a evaluar la curva
    numero_puntos = 500
    
    # Dimensión de los puntos (si queremos graficar en 3D, poner 3 y adaptar el archivo de puntos)
    dim_puntos = 2
    
    print("--CÁLCULO Y DIBUJO DE CURVAS B-SPLINE--")
    archivo_entrada = input("Introduce el nombre del archivo con los puntos (ej. puntos.txt): ")
    grado = int(input("Introduce el grado (tiene que ser igual o menor al núm. de puntos): "))
    
    try:
        # 1. Lectura de los puntos de control del archivo
        puntos_control = leer_puntos_control(archivo_entrada, dim_puntos)
        n = len(puntos_control) - 1
        
        # 2. Cálculos previos
        U = calcular_vector_nodos(n, grado, a, b)
        vect_evaluacion = np.linspace(a, b, numero_puntos + 1)
        
        # 3. Evaluación de la curva
        curva_final = np.zeros((len(vect_evaluacion), dim_puntos))
        
        for i, t in enumerate(vect_evaluacion):
            curva_final[i] = cox_de_boor(t, puntos_control, U, grado)
            
        print("¡Curva calculada con éxito! Generando gráfica...")
        
        # Extraemos las 'X' e 'Y' de los puntos de control originales
        x_control = puntos_control[:, 0]
        y_control = puntos_control[:, 1]
        
        # Extraemos la curva_final:
        x_curva = curva_final[:, 0]
        y_curva = curva_final[:, 1]
        
        plt.figure(figsize=(8, 6))
        
        # Pintamos los puntos de control (con línea de puntos gris)
        plt.plot(x_control, y_control, color='gray', linestyle='--', marker='o', label='Puntos de Control')
        
        # Pintamos la curva B-Spline suave (en azul)
        plt.plot(x_curva, y_curva, color='blue', linewidth=1.5, label=f'Curva B-Spline (Grado {grado})')
        
        plt.title('Curva B-Spline')
        plt.xlabel('Eje X')
        plt.ylabel('Eje Y')
        plt.legend(loc='lower right', fontsize=8)
        plt.grid(True)
        plt.axis('equal') # Ejes proporcionados para que no se deforme la curva
        
        # Mostramos la gráfica en pantalla
        plt.show()
        
    except FileNotFoundError:
        print(f"\n[ERROR] No se ha encontrado el archivo '{archivo_entrada}'.")
        print("Asegúrate de poner el nombre exacto o la ruta completa.")
    except Exception as e:
        print("\n[ERROR] Hubo un problema durante la ejecución:", e)