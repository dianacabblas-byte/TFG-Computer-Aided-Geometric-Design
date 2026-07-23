import urllib.request
import matplotlib.pyplot as plt
import numpy as np

def calcular_vector_nodos(n, grado, a, b):
    ptos_internos = n - grado
    h = (b - a) / (ptos_internos + 1)
    
    inicio = [a] * (grado + 1)
    medio = [a + i * h for i in range(1, ptos_internos + 1)]
    fin = [b] * (grado + 1)
    
    return inicio + medio + fin

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
            
            # Denominador de las formulas
            dif = U[i + grado + 1] - U[i + j]
            if abs(dif) < 1e-9:  
                Q[i, :] = 0.0
            else:
                factor_1 = (U[i + grado + 1] - u) / dif
                factor_2 = (u - U[i + j]) / dif
                
                Q[i, :] = factor_1 * Q[i, :] + factor_2 * Q[i+1, :]

    return Q[r - grado, :]

def evaluar_superficie(red_control, grado_u, grado_v, num_u, num_v):
    a, b = 0.0, 1.0
    c, d = 0.0, 1.0
    
    n = len(red_control) - 1
    m = len(red_control[0]) - 1
    
    # Cálculos previos
    U_u = calcular_vector_nodos(n, grado_u, a, b)
    U_v = calcular_vector_nodos(m, grado_v, c, d)
    
    vect_eval_u = np.linspace(a, b, num_u + 1)
    vect_eval_v = np.linspace(c, d, num_v + 1)
    
    # Pre-creamos las matrices
    X = np.zeros((num_u + 1, num_v + 1))
    Y = np.zeros((num_u + 1, num_v + 1))
    Z = np.zeros((num_u + 1, num_v + 1))
    
    curvas_v = np.zeros((num_v + 1, n + 1, 3))
    
    # PASO 1 (Dirección V)
    for j, v in enumerate(vect_eval_v):
        for i in range(n + 1):
            curvas_v[j, i, :] = cox_de_boor(v, red_control[i, :, :], U_v, grado_v)
            
    # PASO 2 (Dirección U)
    for i, u in enumerate(vect_eval_u):
        for j, v in enumerate(vect_eval_v):
            puntos_intermedios = curvas_v[j, :, :]
            punto_final = cox_de_boor(u, puntos_intermedios, U_u, grado_u)
            
            X[i, j] = punto_final[0]
            Y[i, j] = punto_final[1]
            Z[i, j] = punto_final[2]
        
    return X, Y, Z

if __name__ == "__main__":

    try:
        # 1. Descarga de los puntos originales desde GitHub
        url = 'https://raw.githubusercontent.com/DarrenTsung/bezier-surfaces/master/teapot.bez'
        response = urllib.request.urlopen(url)
        lineas = response.read().decode('utf-8').splitlines()
        
        # 2. Parseo del archivo 
        num_parches = int(lineas[0].strip())
        parches = []
        
        idx = 1
        for _ in range(num_parches):
            while idx < len(lineas) and lineas[idx].strip() == "":
                idx += 1
                
            red_control = np.zeros((4, 4, 3))
            for i in range(4):
                numeros = list(map(float, lineas[idx].split()))
                for j in range(4):
                    # Troceamos la lista usando la misma matemática que en las curvas
                    red_control[i, j, :] = numeros[j*3 : (j+1)*3]
                idx += 1
                
            parches.append(red_control)
        
        print(f"Calculando {num_parches} parches B-spline (esto tardará unos segundos)...")
        
        # 3. Gráfica
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Iteramos por cada red de control y la dibujamos
        for p, red in enumerate(parches):
            X, Y, Z = evaluar_superficie(red, 3, 3, 15, 15)
            ax.plot_surface(X, Y, Z, color='tab:olive', edgecolor='black', linewidth=0.1, alpha=0.9)
        
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        
        # Proporciones reales de la caja para que no se deforme
        rango_x = ax.get_xlim()[1] - ax.get_xlim()[0]
        rango_y = ax.get_ylim()[1] - ax.get_ylim()[0]
        rango_z = ax.get_zlim()[1] - ax.get_zlim()[0]
        ax.set_box_aspect([rango_x, rango_y, rango_z])

        plt.show()
        
    except Exception as e:
        print("\n[ERROR] Hubo un problema durante la ejecución:", e)