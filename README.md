# TFG: Aplicación de Splines Polinómicos en el Diseño Geométrico Asistido por Ordenador

Este repositorio contiene la implementación práctica en Python de los conceptos matemáticos desarrollados en el Trabajo de Fin de Grado. Abarca desde los fundamentos de la interpolación y aproximación clásica, hasta el modelado de curvas y superficies B-Spline.

## 📂 Estructura del Proyecto

El repositorio está dividido en dos partes principales:

1. **Scripts Generales (`curvabspline.py` y `superficiebspline.py`):** 
   Son las implementaciones matemáticas base de las B-splines. Estos algoritmos están preparados para leer los puntos de control desde archivos de texto externos y generar la curva o superficie resultante.
   
2. **Carpeta `py_ejemplos/`:**
   Contiene los scripts específicos (Figura 1 a Figura 17) que se han utilizado para generar todas las imágenes, gráficas y demostraciones visuales que ilustran la memoria del TFG (fenómeno de Runge, modificación local, interpolación de esquinas, curvas isoparamétricas, etc.).

## ⚙️ Formato de los Archivos de Entrada (.txt)

Para utilizar los scripts generales, es necesario proporcionar los puntos de control mediante archivos de texto. En el repositorio se incluyen `puntos.txt` y `red.txt` como ejemplos.

### Para Curvas B-Spline (`puntos.txt`)
El archivo debe contener en la primera línea el número total de puntos, seguido de una lista en columna con las coordenadas de cada uno.

**Formato:**
```text
[Número total de puntos]
x0 y0
x1 y1
x2 y2
...
```

### Para Superficies B-Spline (`red.txt`)
El archivo para superficies es muy intuitivo, ya que su estructura visual imita la propia malla (grid) de la red de control. Las dos primeras líneas definen las dimensiones (filas y columnas). A continuación, cada línea del archivo representa una fila entera de la red de control, agrupando las coordenadas `x y z` de cada punto de esa fila.

**Formato:**
```text
[Número de filas]
[Número de columnas]
x00 y00 z00   x01 y01 z01   x02 y02 z02 ...
x10 y10 z10   x11 y11 z11   x12 y12 z12 ...
x20 y20 z20   x21 y21 z21   x22 y22 z22 ...
...
```

## 🚀 Cómo ejecutar los ejemplos

Para visualizar cualquiera de las figuras creadas para la memoria del TFG, basta con ejecutar su script correspondiente usando Python. Por ejemplo:

```bash
python py_ejemplos/figura_15_esquinas_superficie.py
```
