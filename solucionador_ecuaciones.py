#DAVID ALEXANDER FONSECA PEREZ
#Programa para resolver ecuaciones no lineales utilizando los metodos de:
#1. Biseccion
#2. Newton-Raphson
#3. Secante

"""
Programa para resolver ecuaciones no lineales utilizando los metodos de:
1. Biseccion
2. Newton-Raphson
3. Secante

Este programa permite encontrar raices de ecuaciones no lineales y comparar
la eficiencia y precision de ambos metodos.

Author: [David Alexander Fonseca Perez]
Date: [2025-04-18]
Version: 1.0
"""

import math
import re
import sys
import argparse
import time
import tracemalloc
import matplotlib.pyplot as plt
from typing import Callable, List, Tuple, Dict, Optional, Any, Union

def analizar_ecuacion(ecuacion_str: str) -> Callable[[float], float]:
    """
    Convierte una cadena de texto que representa una ecuacion en una funcion evaluable.
    
    Args:
        ecuacion_str: Cadena de texto que representa la ecuacion a resolver.
    
    Returns:
        Una funcion que evalua la ecuacion para un valor dado de x.
    """
    # Reemplazar operaciones matematicas comunes con sus equivalentes en Python
    ecuacion_str = ecuacion_str.replace("^", "**")
    
    # Insertar operadores de multiplicación donde sea necesario (ej: 4x -> 4*x)
    ecuacion_str = re.sub(r'(\d)([a-zA-Z])', r'\1*\2', ecuacion_str)
    
    # Reemplazar funciones matematicas comunes
    ecuacion_str = re.sub(r'sin\(', 'math.sin(', ecuacion_str)
    ecuacion_str = re.sub(r'cos\(', 'math.cos(', ecuacion_str)
    ecuacion_str = re.sub(r'tan\(', 'math.tan(', ecuacion_str)
    ecuacion_str = re.sub(r'exp\(', 'math.exp(', ecuacion_str)
    ecuacion_str = re.sub(r'log\(', 'math.log(', ecuacion_str)
    ecuacion_str = re.sub(r'ln\(', 'math.log(', ecuacion_str)
    ecuacion_str = re.sub(r'sqrt\(', 'math.sqrt(', ecuacion_str)
    
    # Crear la funcion
    try:
        def f(x):
            try:
                return eval(ecuacion_str)
            except Exception as e:
                print(f"Error al evaluar la ecuacion: {e}")
                return float('nan')
        
        # Probar la funcion para verificar que sea valida
        f(0)
        return f
    except Exception as e:
        print(f"Error al analizar la ecuacion: {e}")
        print("Asegurate de que la ecuacion este correctamente escrita.")
        return None


def derivada(f: Callable[[float], float], x: float, h: float = 1e-6) -> float:
    """
    Calcula la derivada numerica de una funcion en un punto dado.
    
    Args:
        f: Funcion a derivar.
        x: Punto donde se evalua la derivada.
        h: Tamano del paso para la aproximacion numerica.
    
    Returns:
        El valor de la derivada en el punto x.
    """
    return (f(x + h) - f(x)) / h


def calcular_error(actual: float, anterior: float) -> Tuple[float, float]:
    """
    Calcula el error absoluto y relativo entre dos aproximaciones sucesivas.
    
    Args:
        actual: Valor actual.
        anterior: Valor anterior.
    
    Returns:
        Una tupla con el error absoluto y el error relativo.
    """
    error_abs = abs(actual - anterior)
    error_rel = error_abs / abs(actual) if actual != 0 else float('inf')
    return error_abs, error_rel


def metodo_biseccion(f: Callable[[float], float], a: float, b: float, 
                    tol: float = 1e-6, max_iter: int = 100) -> Tuple[float, int, List[float], List[float], List[float], float, float]:
    """
    Implementa el metodo de biseccion para encontrar una raiz de f en [a, b].
    
    Args:
        f: Funcion cuya raiz se busca.
        a: Extremo izquierdo del intervalo.
        b: Extremo derecho del intervalo.
        tol: Tolerancia para el criterio de parada.
        max_iter: Numero maximo de iteraciones.
    
    Returns:
        Una tupla con la raiz aproximada, el numero de iteraciones realizadas,
        una lista con los valores intermedios, una lista con los errores absolutos,
        una lista con los errores relativos, el tiempo de ejecucion y el uso de memoria en bytes.
    """
    # Verificar que f(a) y f(b) tengan signos opuestos
    if f(a) * f(b) >= 0:
        print(f"Error: f(a) = {f(a)} y f(b) = {f(b)} deben tener signos opuestos.")
        return None, 0, [], [], [], 0, 0
    
    # Inicializacion
    tracemalloc.start()
    tiempo_inicio = time.time()
    contador_iter = 0
    valores_c = [a]  # Incluir el valor inicial
    errores_abs = [float('inf')]
    errores_rel = [float('inf')]
    
    # Primera iteracion
    c = (a + b) / 2
    valores_c.append(c)
    error_abs, error_rel = calcular_error(c, a)
    errores_abs.append(error_abs)
    errores_rel.append(error_rel)
    
    while (b - a) / 2 > tol and contador_iter < max_iter:
        # Calcular el punto medio
        c_prev = c
        c = (a + b) / 2
        valores_c.append(c)
        
        # Calcular errores
        error_abs, error_rel = calcular_error(c, c_prev)
        errores_abs.append(error_abs)
        errores_rel.append(error_rel)
        
        # Evaluar la funcion en el punto medio
        fc = f(c)
        
        # Verificar si c es una raiz
        if abs(fc) < tol:
            tiempo_fin = time.time()
            memoria_actual, memoria_pico = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            return c, contador_iter + 1, valores_c, errores_abs, errores_rel, tiempo_fin - tiempo_inicio, memoria_pico
        
        # Actualizar el intervalo
        if f(a) * fc < 0:
            b = c
        else:
            a = c
        
        contador_iter += 1
    
    # Retornar el punto medio del intervalo final
    c = (a + b) / 2
    valores_c.append(c)
    error_abs, error_rel = calcular_error(c, valores_c[-2])
    errores_abs.append(error_abs)
    errores_rel.append(error_rel)
    
    tiempo_fin = time.time()
    memoria_actual, memoria_pico = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return c, contador_iter + 1, valores_c, errores_abs, errores_rel, tiempo_fin - tiempo_inicio, memoria_pico


def metodo_newton_raphson(f: Callable[[float], float], x0: float, 
                         tol: float = 1e-6, max_iter: int = 100) -> Tuple[float, int, List[float], List[float], List[float], float, float]:
    """
    Implementa el metodo de Newton-Raphson para encontrar una raiz de f.
    
    Args:
        f: Funcion cuya raiz se busca.
        x0: Aproximacion inicial.
        tol: Tolerancia para el criterio de parada.
        max_iter: Numero maximo de iteraciones.
    
    Returns:
        Una tupla con la raiz aproximada, el numero de iteraciones realizadas,
        una lista con los valores intermedios, una lista con los errores absolutos,
        una lista con los errores relativos, el tiempo de ejecucion y el uso de memoria en bytes.
    """
    # Inicializacion
    tracemalloc.start()
    tiempo_inicio = time.time()
    x = x0
    contador_iter = 0
    valores_x = [x0]
    errores_abs = [float('inf')]
    errores_rel = [float('inf')]
    
    while contador_iter < max_iter:
        # Calcular la derivada
        df = derivada(f, x)
        
        # Verificar que la derivada no sea cero
        if abs(df) < 1e-10:
            print("Error: La derivada es cercana a cero. El metodo puede no converger.")
            tiempo_fin = time.time()
            memoria_actual, memoria_pico = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            return x, contador_iter, valores_x, errores_abs, errores_rel, tiempo_fin - tiempo_inicio, memoria_pico
        
        # Calcular la nueva aproximacion
        x_nuevo = x - f(x) / df
        valores_x.append(x_nuevo)
        
        # Calcular errores
        error_abs, error_rel = calcular_error(x_nuevo, x)
        errores_abs.append(error_abs)
        errores_rel.append(error_rel)
        
        # Verificar el criterio de parada
        if abs(x_nuevo - x) < tol:
            tiempo_fin = time.time()
            memoria_actual, memoria_pico = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            return x_nuevo, contador_iter + 1, valores_x, errores_abs, errores_rel, tiempo_fin - tiempo_inicio, memoria_pico
        
        # Actualizar x
        x = x_nuevo
        contador_iter += 1
    
    print("Advertencia: Se alcanzo el numero maximo de iteraciones.")
    tiempo_fin = time.time()
    memoria_actual, memoria_pico = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return x, contador_iter, valores_x, errores_abs, errores_rel, tiempo_fin - tiempo_inicio, memoria_pico


def metodo_secante(f: Callable[[float], float], x0: float, x1: float, tol: float = 1e-6, max_iter: int = 100):
    """
    Implementa el metodo de la secante para encontrar una raiz de f.
    
    Args:
        f: Funcion cuya raiz se busca.
        x0: Primer valor inicial.
        x1: Segundo valor inicial.
        tol: Tolerancia para el criterio de parada.
        max_iter: Numero maximo de iteraciones.
    
    Returns:
        Una tupla con la raiz aproximada, el numero de iteraciones realizadas,
        una lista con los valores intermedios, una lista con los errores absolutos,
        una lista con los errores relativos, el tiempo de ejecucion y el uso de memoria en bytes.
    """
    valores = [x0, x1]
    errores_abs = [None]  # El primer error no se puede calcular
    errores_rel = [None]
    x_prev = x0
    x_curr = x1
    raiz = None
    start_time = time.time()
    tracemalloc.start()
    for i in range(max_iter):
        try:
            f_x_prev = f(x_prev)
            f_x_curr = f(x_curr)
            if f_x_curr - f_x_prev == 0:
                print("\n[Secante] División por cero en la iteración", i)
                break
            x_next = x_curr - f_x_curr * (x_curr - x_prev) / (f_x_curr - f_x_prev)
        except Exception as e:
            print(f"\n[Secante] Error en la iteración {i}: {e}")
            break
        valores.append(x_next)
        error_abs, error_rel = calcular_error(x_next, x_curr)
        errores_abs.append(error_abs)
        errores_rel.append(error_rel)
        if error_abs < tol:
            raiz = x_next
            break
        x_prev, x_curr = x_curr, x_next
    else:
        raiz = valores[-1]
    end_time = time.time()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    tiempo_ejecucion = end_time - start_time
    memoria_usada = peak
    return raiz, len(valores) - 2, valores[2:], errores_abs[1:], errores_rel[1:], tiempo_ejecucion, memoria_usada


def imprimir_tabla(nombre_metodo: str, valores: List[float], errores_abs: List[float], 
               errores_rel: List[float], f: Callable[[float], float]) -> None:
    """
    Imprime una tabla con los valores intermedios y la evaluacion de la funcion.
    
    Args:
        nombre_metodo: Nombre del metodo utilizado.
        valores: Lista de valores intermedios.
        errores_abs: Lista de errores absolutos.
        errores_rel: Lista de errores relativos.
        f: Funcion evaluada.
    """
    print(f"\nTabla de iteraciones para el metodo de {nombre_metodo}:")
    print("-" * 80)
    print(f"{'Iteracion':^10} | {'Valor de x':^15} | {'f(x)':^15} | {'Error Abs':^15} | {'Error Rel':^15}")
    print("-" * 80)
    
    for i, (x, err_abs, err_rel) in enumerate(zip(valores, errores_abs, errores_rel)):
        fx = f(x)
        print(f"{i:^10} | {x:^15.8f} | {fx:^15.8e} | {err_abs:^15.8e} | {err_rel:^15.8e}")


def comparar_metodos(resultados_biseccion, resultados_newton, resultados_secante, ecuacion_str: str) -> None:
    """
    Compara los resultados de los metodos de biseccion, Newton-Raphson y secante.
    
    Args:
        resultados_biseccion: Resultados del metodo de biseccion.
        resultados_newton: Resultados del metodo de Newton-Raphson.
        resultados_secante: Resultados del metodo de la secante.
        ecuacion_str: Ecuacion resuelta.
    """
    if resultados_biseccion[0] is None or resultados_newton[0] is None or resultados_secante[0] is None:
        print("\nNo se pueden comparar los metodos porque al menos uno fallo.")
        return
    
    print("\n" + "=" * 80)
    print("COMPARACION DE METODOS".center(80))
    print("=" * 80)
    
    print(f"\nEcuacion resuelta: {ecuacion_str}")
    
    # Crear tabla comparativa
    print("\nResumen comparativo:")
    print("-" * 80)
    print(f"{'Criterio':^20} | {'Biseccion':^25} | {'Newton-Raphson':^25} | {'Secante':^25}")
    print("-" * 80)
    
    # Raiz encontrada
    print(f"{'Raiz aproximada':^20} | {resultados_biseccion[0]:^25.10f} | {resultados_newton[0]:^25.10f} | {resultados_secante[0]:^25.10f}")
    
    # Diferencia entre raices
    diff_biseccion_newton = abs(resultados_biseccion[0] - resultados_newton[0])
    diff_biseccion_secante = abs(resultados_biseccion[0] - resultados_secante[0])
    diff_newton_secante = abs(resultados_newton[0] - resultados_secante[0])
    print(f"{'Diferencia raices':^20} | {diff_biseccion_newton:^51.10e} | {diff_biseccion_secante:^51.10e} | {diff_newton_secante:^51.10e}")
    
    # Numero de iteraciones
    print(f"{'Iteraciones':^20} | {resultados_biseccion[1]:^25d} | {resultados_newton[1]:^25d} | {resultados_secante[1]:^25d}")
    
    # Error final
    print(f"{'Error Abs final':^20} | {resultados_biseccion[3][-1]:^25.10e} | {resultados_newton[3][-1]:^25.10e} | {resultados_secante[3][-1]:^25.10e}")
    print(f"{'Error Rel final':^20} | {resultados_biseccion[4][-1]:^25.10e} | {resultados_newton[4][-1]:^25.10e} | {resultados_secante[4][-1]:^25.10e}")
    
    # Tiempo de ejecucion
    print(f"{'Tiempo (s)':^20} | {resultados_biseccion[5]:^25.6f} | {resultados_newton[5]:^25.6f} | {resultados_secante[5]:^25.6f}")
    
    # Uso de memoria
    print(f"{'Uso de memoria (bytes)':^20} | {resultados_biseccion[6]:^25d} | {resultados_newton[6]:^25d} | {resultados_secante[6]:^25d}")
    
    # Conclusiones
    print("\nConclusiones:")
    
    # Metodo mas rapido
    if resultados_biseccion[5] < resultados_newton[5] and resultados_biseccion[5] < resultados_secante[5]:
        print("- El metodo de Biseccion fue mas rapido en tiempo de ejecucion.")
    elif resultados_newton[5] < resultados_biseccion[5] and resultados_newton[5] < resultados_secante[5]:
        print("- El metodo de Newton-Raphson fue mas rapido en tiempo de ejecucion.")
    elif resultados_secante[5] < resultados_biseccion[5] and resultados_secante[5] < resultados_newton[5]:
        print("- El metodo de la Secante fue mas rapido en tiempo de ejecucion.")
    else:
        print("- Ambos metodos tuvieron tiempos de ejecucion similares.")
    
    # Metodo con menos iteraciones
    if resultados_biseccion[1] < resultados_newton[1] and resultados_biseccion[1] < resultados_secante[1]:
        print("- El metodo de Biseccion requirio menos iteraciones.")
    elif resultados_newton[1] < resultados_biseccion[1] and resultados_newton[1] < resultados_secante[1]:
        print("- El metodo de Newton-Raphson requirio menos iteraciones.")
    elif resultados_secante[1] < resultados_biseccion[1] and resultados_secante[1] < resultados_newton[1]:
        print("- El metodo de la Secante requirio menos iteraciones.")
    else:
        print("- Ambos metodos requirieron el mismo numero de iteraciones.")
    
    # Metodo con menor error
    if resultados_biseccion[3][-1] < resultados_newton[3][-1] and resultados_biseccion[3][-1] < resultados_secante[3][-1]:
        print("- El metodo de Biseccion tuvo menor error absoluto final.")
    elif resultados_newton[3][-1] < resultados_biseccion[3][-1] and resultados_newton[3][-1] < resultados_secante[3][-1]:
        print("- El metodo de Newton-Raphson tuvo menor error absoluto final.")
    elif resultados_secante[3][-1] < resultados_biseccion[3][-1] and resultados_secante[3][-1] < resultados_newton[3][-1]:
        print("- El metodo de la Secante tuvo menor error absoluto final.")
    else:
        print("- Ambos metodos tuvieron errores absolutos finales similares.")
    
    # Método con menor uso de memoria
    if resultados_biseccion[6] < resultados_newton[6] and resultados_biseccion[6] < resultados_secante[6]:
        print("- El metodo de Biseccion utilizó menos memoria.")
    elif resultados_newton[6] < resultados_biseccion[6] and resultados_newton[6] < resultados_secante[6]:
        print("- El metodo de Newton-Raphson utilizó menos memoria.")
    elif resultados_secante[6] < resultados_biseccion[6] and resultados_secante[6] < resultados_newton[6]:
        print("- El metodo de la Secante utilizó menos memoria.")
    else:
        print("- Ambos metodos utilizaron cantidades similares de memoria.")
    
    # Observaciones adicionales
    print("\nObservaciones adicionales:")
    if resultados_newton[1] < resultados_biseccion[1] and resultados_newton[5] < resultados_biseccion[5]:
        print("- El metodo de Newton-Raphson mostro una convergencia mas eficiente en este caso.")
    elif resultados_biseccion[1] < resultados_newton[1] and resultados_biseccion[5] < resultados_newton[5]:
        print("- El metodo de Biseccion mostro una convergencia mas eficiente en este caso.")
    elif resultados_secante[1] < resultados_biseccion[1] and resultados_secante[5] < resultados_biseccion[5]:
        print("- El metodo de la Secante mostro una convergencia mas eficiente en este caso.")
    else:
        print("- Cada metodo tiene sus ventajas en diferentes aspectos de la convergencia.")
    
    # Recomendacion
    print("\nRecomendacion:")
    if resultados_newton[1] < resultados_biseccion[1] and resultados_newton[3][-1] < resultados_biseccion[3][-1]:
        print("- Para esta ecuacion especifica, el metodo de Newton-Raphson es preferible.")
    elif resultados_biseccion[1] < resultados_newton[1] and resultados_biseccion[3][-1] < resultados_newton[3][-1]:
        print("- Para esta ecuacion especifica, el metodo de Biseccion es preferible.")
    elif resultados_secante[1] < resultados_biseccion[1] and resultados_secante[3][-1] < resultados_biseccion[3][-1]:
        print("- Para esta ecuacion especifica, el metodo de la Secante es preferible.")
    else:
        print("- Ambos metodos son adecuados para esta ecuacion, pero con diferentes caracteristicas de convergencia.")
    
    # Crear visualizaciones gráficas
    crear_graficas_comparativas(resultados_biseccion, resultados_newton, resultados_secante, ecuacion_str)


def crear_graficas_comparativas(resultados_biseccion, resultados_newton, resultados_secante, ecuacion_str: str) -> None:
    """
    Crea gráficas comparativas entre los métodos de bisección, Newton-Raphson y secante.
    
    Args:
        resultados_biseccion: Resultados del método de bisección.
        resultados_newton: Resultados del método de Newton-Raphson.
        resultados_secante: Resultados del método de la secante.
        ecuacion_str: Ecuación resuelta.
    """
    # Configurar el estilo de las gráficas
    plt.style.use('seaborn-v0_8-darkgrid')
    
    # Crear figura con subplots
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'Comparación de Métodos para la ecuación: {ecuacion_str}', fontsize=16)
    
    # 1. Gráfica de convergencia (valores vs iteraciones)
    axs[0, 0].plot(range(len(resultados_biseccion[2])), resultados_biseccion[2], 'b.-', label='Bisección')
    axs[0, 0].plot(range(len(resultados_newton[2])), resultados_newton[2], 'r.-', label='Newton-Raphson')
    axs[0, 0].plot(range(len(resultados_secante[2])), resultados_secante[2], 'g.-', label='Secante')
    axs[0, 0].set_title('Convergencia de los métodos')
    axs[0, 0].set_xlabel('Iteración')
    axs[0, 0].set_ylabel('Valor de x')
    axs[0, 0].legend()
    axs[0, 0].grid(True)
    
    # 2. Gráfica de error absoluto vs iteraciones
    axs[0, 1].semilogy(range(len(resultados_biseccion[3])), resultados_biseccion[3], 'b.-', label='Bisección')
    axs[0, 1].semilogy(range(len(resultados_newton[3])), resultados_newton[3], 'r.-', label='Newton-Raphson')
    axs[0, 1].semilogy(range(len(resultados_secante[3])), resultados_secante[3], 'g.-', label='Secante')
    axs[0, 1].set_title('Error absoluto vs Iteraciones')
    axs[0, 1].set_xlabel('Iteración')
    axs[0, 1].set_ylabel('Error absoluto (escala log)')
    axs[0, 1].legend()
    axs[0, 1].grid(True)
    
    # 3. Gráfica de error relativo vs iteraciones
    axs[1, 0].semilogy(range(len(resultados_biseccion[4])), resultados_biseccion[4], 'b.-', label='Bisección')
    axs[1, 0].semilogy(range(len(resultados_newton[4])), resultados_newton[4], 'r.-', label='Newton-Raphson')
    axs[1, 0].semilogy(range(len(resultados_secante[4])), resultados_secante[4], 'g.-', label='Secante')
    axs[1, 0].set_title('Error relativo vs Iteraciones')
    axs[1, 0].set_xlabel('Iteración')
    axs[1, 0].set_ylabel('Error relativo (escala log)')
    axs[1, 0].legend()
    axs[1, 0].grid(True)
    
    # 4. Gráfica de barras comparando tiempo y memoria
    criterios = ['Tiempo (s)', 'Memoria (KB)']
    biseccion_vals = [resultados_biseccion[5], resultados_biseccion[6]/1024]  # Convertir bytes a KB
    newton_vals = [resultados_newton[5], resultados_newton[6]/1024]  # Convertir bytes a KB
    secante_vals = [resultados_secante[5], resultados_secante[6]/1024]  # Convertir bytes a KB
    
    x = range(len(criterios))
    width = 0.35
    
    axs[1, 1].bar([i - width for i in x], biseccion_vals, width, label='Bisección', color='blue')
    axs[1, 1].bar(x, newton_vals, width, label='Newton-Raphson', color='red')
    axs[1, 1].bar([i + width for i in x], secante_vals, width, label='Secante', color='green')
    
    axs[1, 1].set_title('Comparación de recursos')
    axs[1, 1].set_xticks(x)
    axs[1, 1].set_xticklabels(criterios)
    axs[1, 1].legend()
    
    # Ajustar el diseño
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    
    # Crear un nombre de archivo seguro basado en la ecuación
    nombre_archivo = ecuacion_str.replace(' ', '_')
    # Reemplazar caracteres no seguros para nombres de archivo
    nombre_archivo = re.sub(r'[^\w\-_\.]', '', nombre_archivo)
    # Limitar la longitud del nombre
    if len(nombre_archivo) > 30:
        nombre_archivo = nombre_archivo[:30]
    # Añadir un sufijo si el nombre está vacío
    if not nombre_archivo:
        nombre_archivo = "ecuacion"
    
    nombre_archivo = f"comparacion_{nombre_archivo}.png"
    
    # Guardar la figura
    plt.savefig(nombre_archivo)
    
    print(f"\nSe ha generado una gráfica comparativa y guardado como '{nombre_archivo}'")
    
    # Mostrar la figura (opcional, dependiendo del entorno)
    try:
        plt.show()
    except Exception as e:
        print(f"Nota: No se pudo mostrar la gráfica interactivamente. Error: {e}")
        print(f"La gráfica ha sido guardada como '{nombre_archivo}'")
    finally:
        plt.close()


def entrada_segura(mensaje: str, valor_predeterminado=None):
    """
    Funcion segura para manejar la entrada del usuario, con proteccion contra EOF.
    
    Args:
        mensaje: Mensaje para mostrar al usuario.
        valor_predeterminado: Valor por defecto si hay un error o entrada vacia.
    
    Returns:
        La entrada del usuario o el valor por defecto.
    """
    try:
        valor = input(mensaje)
        if valor.strip() == "" and valor_predeterminado is not None:
            return valor_predeterminado
        return valor
    except EOFError:
        print("Error: Se detecto EOF al leer la entrada.")
        if valor_predeterminado is not None:
            print(f"Usando valor por defecto: {valor_predeterminado}")
            return valor_predeterminado
        else:
            print("No hay valor por defecto disponible. Saliendo...")
            sys.exit(1)
    except KeyboardInterrupt:
        print("\nOperacion cancelada por el usuario.")
        sys.exit(0)


def resolver_con_argumentos(args):
    """
    Resuelve una ecuacion utilizando argumentos de linea de comandos.
    
    Args:
        args: Argumentos de linea de comandos parseados.
    """
    # Parsear la ecuacion
    f = analizar_ecuacion(args.equation)
    if f is None:
        sys.exit(1)
    
    # Parametros comunes
    tol = args.tolerance
    max_iter = args.max_iterations
    
    # Variables para almacenar resultados
    resultados_biseccion = None
    resultados_newton = None
    resultados_secante = None
    
    # Ejecutar el metodo seleccionado
    if args.method in ['biseccion', 'both']:
        print("\n" + "-" * 60)
        print("METODO DE BISECCION".center(60))
        print("-" * 60)
        
        # Verificar que se proporcionaron los limites del intervalo
        if args.a is None or args.b is None:
            print("Error: Para el metodo de biseccion se requieren los limites del intervalo (a y b).")
            if args.method != 'both':
                sys.exit(1)
        else:
            a, b = args.a, args.b
            
            # Verificar que f(a) y f(b) tengan signos opuestos
            fa, fb = f(a), f(b)
            if fa * fb >= 0:
                print(f"Error: f({a}) = {fa} y f({b}) = {fb} deben tener signos opuestos.")
                if args.method == 'both':
                    print("\nContinuando con el metodo de Newton-Raphson...")
                else:
                    sys.exit(1)
            else:
                # Ejecutar el metodo de biseccion
                resultados_biseccion = metodo_biseccion(f, a, b, tol, max_iter)
                
                if resultados_biseccion[0] is not None:
                    print(f"\nRaiz encontrada: {resultados_biseccion[0]:.10f}")
                    print(f"Valor de f(raiz): {f(resultados_biseccion[0]):.10e}")
                    print(f"Iteraciones realizadas: {resultados_biseccion[1]}")
                    print(f"Error absoluto final: {resultados_biseccion[3][-1]:.10e}")
                    print(f"Error relativo final: {resultados_biseccion[4][-1]:.10e}")
                    print(f"Tiempo de ejecucion: {resultados_biseccion[5]:.6f} segundos")
                    print(f"Uso de memoria: {resultados_biseccion[6]} bytes")
                    
                    # Imprimir tabla de iteraciones
                    imprimir_tabla("Biseccion", resultados_biseccion[2], resultados_biseccion[3], 
                               resultados_biseccion[4], f)
    
    if args.method in ['newton', 'both']:
        print("\n" + "-" * 60)
        print("METODO DE NEWTON-RAPHSON".center(60))
        print("-" * 60)
        
        # Verificar que se proporciono la aproximacion inicial
        if args.x0 is None:
            print("Error: Para el metodo de Newton-Raphson se requiere una aproximacion inicial (x0).")
            sys.exit(1)
        
        # Ejecutar el metodo de Newton-Raphson
        resultados_newton = metodo_newton_raphson(f, args.x0, tol, max_iter)
        
        print(f"\nRaiz encontrada: {resultados_newton[0]:.10f}")
        print(f"Valor de f(raiz): {f(resultados_newton[0]):.10e}")
        print(f"Iteraciones realizadas: {resultados_newton[1]}")
        print(f"Error absoluto final: {resultados_newton[3][-1]:.10e}")
        print(f"Error relativo final: {resultados_newton[4][-1]:.10e}")
        print(f"Tiempo de ejecucion: {resultados_newton[5]:.6f} segundos")
        print(f"Uso de memoria: {resultados_newton[6]} bytes")
        
        # Imprimir tabla de iteraciones
        imprimir_tabla("Newton-Raphson", resultados_newton[2], resultados_newton[3], 
                   resultados_newton[4], f)
    
    if args.method in ['secante', 'both']:
        print("\n" + "-" * 60)
        print("METODO DE LA SECANTE".center(60))
        print("-" * 60)
        
        # Verificar que se proporcionaron los valores iniciales
        if args.x0 is None or args.x1 is None:
            print("Error: Para el metodo de la secante se requieren dos valores iniciales (x0 y x1).")
            sys.exit(1)
        
        # Ejecutar el metodo de la secante
        resultados_secante = metodo_secante(f, args.x0, args.x1, tol, max_iter)
        
        print(f"\nRaiz encontrada: {resultados_secante[0]:.10f}")
        print(f"Valor de f(raiz): {f(resultados_secante[0]):.10e}")
        print(f"Iteraciones realizadas: {resultados_secante[1]}")
        print(f"Error absoluto final: {resultados_secante[3][-1]:.10e}")
        print(f"Error relativo final: {resultados_secante[4][-1]:.10e}")
        print(f"Tiempo de ejecucion: {resultados_secante[5]:.6f} segundos")
        print(f"Uso de memoria: {resultados_secante[6]} bytes")
        
        # Imprimir tabla de iteraciones
        imprimir_tabla("Secante", resultados_secante[2], resultados_secante[3], 
                   resultados_secante[4], f)
    
    # Comparar metodos si se ejecutaron ambos
    if args.method == 'both' and resultados_biseccion is not None and resultados_newton is not None and resultados_secante is not None:
        comparar_metodos(resultados_biseccion, resultados_newton, resultados_secante, args.equation)


def modo_interactivo():
    """Funcion para el modo interactivo del programa."""
    print("\n" + "=" * 70)
    print(" SOLUCIONADOR DE ECUACIONES NO LINEALES ".center(70))
    print("=" * 70)
    
    print("\n DESCRIPCION:")
    print("  Este programa te ayuda a encontrar las raices (soluciones) de ecuaciones")
    print("  no lineales utilizando tres metodos numericos diferentes:")
    print("  1. Metodo de Biseccion: mas lento pero muy estable")
    print("  2. Metodo de Newton-Raphson: mas rapido pero requiere buena aproximacion inicial")
    print("  3. Metodo de la Secante: rapido pero requiere dos valores iniciales")
    
    print("\n INSTRUCCIONES PARA INGRESAR LA ECUACION:")
    print("  • Usa 'x' como variable (por ejemplo: x^2 - 4)")
    print("  • Operaciones disponibles: + (suma), - (resta), * (multiplicacion),")
    print("    / (division), ^ (potencia)")
    print("  • Funciones matematicas disponibles:")
    print("    - sin(x): seno de x")
    print("    - cos(x): coseno de x")
    print("    - tan(x): tangente de x")
    print("    - exp(x): exponencial (e^x)")
    print("    - log(x): logaritmo natural (base e)")
    print("    - sqrt(x): raiz cuadrada")
    
    print("\n EJEMPLOS DE ECUACIONES:")
    print("  • x^2 - 4 = 0       (soluciones: x = 2 y x = -2)")
    print("  • sin(x) = 0        (soluciones: x = 0, pi, 2pi, ...)")
    print("  • exp(x) - 2 = 0    (solucion: x = ln(2))")
    print("  • x^3 - 2*x - 5 = 0 (solucion numerica)")
    
    while True:
        try:
            # Solicitar la ecuacion
            print("\n" + "-" * 70)
            print(" PASO 1: INGRESA LA ECUACION".center(70))
            print("-" * 70)
            print("Escribe la ecuacion igualada a cero. Por ejemplo: x^2 - 4")
            print("(Para salir del programa escribe 'salir', 'exit' o 'q')")
            ecuacion_str = entrada_segura("\n Tu ecuacion: ")
            
            # Verificar si el usuario quiere salir
            if ecuacion_str.lower() in ['salir', 'exit', 'q', 'quit']:
                print("\n Gracias por usar el programa! ")
                break
            
            # Parsear la ecuacion
            f = analizar_ecuacion(ecuacion_str)
            if f is None:
                continue
            
            # Solicitar el metodo
            print("\n" + "-" * 70)
            print(" PASO 2: SELECCIONA EL METODO".center(70))
            print("-" * 70)
            print("1. Metodo de Biseccion (requiere intervalo [a,b] donde f(a) y f(b) tengan signos opuestos)")
            print("2. Metodo de Newton-Raphson (requiere un valor inicial aproximado)")
            print("3. Metodo de la Secante (requiere dos valores iniciales)")
            print("4. Ambos metodos (comparar resultados)")
            
            method = entrada_segura("\n Selecciona una opcion (1, 2, 3 o 4): ")
            
            if method not in ['1', '2', '3', '4']:
                print(" Opcion no valida. Por favor, selecciona 1, 2, 3 o 4.")
                continue
            
            # Parametros comunes
            print("\n" + "-" * 70)
            print(" PASO 3: CONFIGURA LOS PARAMETROS".center(70))
            print("-" * 70)
            
            print("\n TOLERANCIA:")
            print("  La tolerancia determina cuando el metodo se detiene.")
            print("  Un valor mas pequeno da resultados mas precisos pero puede tardar mas.")
            print("  Valor recomendado: 0.000001 (1e-6)")
            tol = float(entrada_segura("\n Ingresa la tolerancia (presiona ENTER para usar 1e-6): ", "1e-6"))
            
            print("\n ITERACIONES MAXIMAS:")
            print("  Numero maximo de pasos que realizara el algoritmo.")
            print("  Si el metodo no converge en este numero de pasos, se detendra.")
            print("  Valor recomendado: 100")
            max_iter = int(entrada_segura(" Ingresa el numero maximo de iteraciones (presiona ENTER para usar 100): ", "100"))
            
            # Variables para almacenar resultados
            resultados_biseccion = None
            resultados_newton = None
            resultados_secante = None
            
            # Ejecutar el metodo seleccionado
            if method in ['1', '4']:  # Biseccion
                print("\n" + "-" * 70)
                print(" METODO DE BISECCION".center(70))
                print("-" * 70)
                
                print("\n EXPLICACION:")
                print("  El metodo de biseccion requiere un intervalo [a,b] donde:")
                print("  • f(a) y f(b) deben tener signos opuestos")
                print("  • La raiz debe estar dentro de este intervalo")
                
                a = float(entrada_segura("\n Ingresa el extremo izquierdo del intervalo (a): "))
                b = float(entrada_segura(" Ingresa el extremo derecho del intervalo (b): "))
                
                # Verificar que f(a) y f(b) tengan signos opuestos
                fa, fb = f(a), f(b)
                if fa * fb >= 0:
                    print(f"\n ERROR: f({a}) = {fa} y f({b}) = {fb} deben tener signos opuestos.")
                    print("   Esto significa que la funcion debe cambiar de signo en el intervalo [a,b].")
                    print("   Intenta con otros valores donde uno sea positivo y otro negativo.")
                    if method == '4':
                        print("\n Continuando con el metodo de Newton-Raphson...")
                    else:
                        continue
                else:
                    print("\n Calculando solucion mediante biseccion...")
                    # Ejecutar el metodo de biseccion
                    resultados_biseccion = metodo_biseccion(f, a, b, tol, max_iter)
                    
                    if resultados_biseccion[0] is not None:
                        print("\n RESULTADOS DEL METODO DE BISECCION:")
                        print(f"  • Raiz encontrada: {resultados_biseccion[0]:.10f}")
                        print(f"  • Valor de f(raiz): {f(resultados_biseccion[0]):.10e}")
                        print(f"  • Iteraciones realizadas: {resultados_biseccion[1]}")
                        print(f"  • Error absoluto final: {resultados_biseccion[3][-1]:.10e}")
                        print(f"  • Error relativo final: {resultados_biseccion[4][-1]:.10e}")
                        print(f"  • Tiempo de ejecucion: {resultados_biseccion[5]:.6f} segundos")
                        print(f"  • Uso de memoria: {resultados_biseccion[6]} bytes")
                        
                        # Preguntar si desea ver la tabla de iteraciones
                        mostrar_tabla = entrada_segura("\n ¿Deseas ver la tabla de iteraciones? (s/n): ")
                        if mostrar_tabla.lower() in ['s', 'si', 'sí', 'y', 'yes']:
                            # Imprimir tabla de iteraciones
                            imprimir_tabla("Biseccion", resultados_biseccion[2], resultados_biseccion[3], 
                                       resultados_biseccion[4], f)
            
            if method in ['2', '4']:  # Newton-Raphson
                print("\n" + "-" * 70)
                print(" METODO DE NEWTON-RAPHSON".center(70))
                print("-" * 70)
                
                print("\n EXPLICACION:")
                print("  El metodo de Newton-Raphson requiere:")
                print("  • Un valor inicial x₀ cercano a la raiz")
                print("  • Converge muy rapido si el valor inicial es bueno")
                print("  • Puede diverger si el valor inicial esta lejos de la raiz")
                
                x0 = float(entrada_segura("\n Ingresa la aproximacion inicial (x₀): "))
                
                print("\n Calculando solucion mediante Newton-Raphson...")
                # Ejecutar el metodo de Newton-Raphson
                resultados_newton = metodo_newton_raphson(f, x0, tol, max_iter)
                
                print("\n RESULTADOS DEL METODO DE NEWTON-RAPHSON:")
                print(f"  • Raiz encontrada: {resultados_newton[0]:.10f}")
                print(f"  • Valor de f(raiz): {f(resultados_newton[0]):.10e}")
                print(f"  • Iteraciones realizadas: {resultados_newton[1]}")
                print(f"  • Error absoluto final: {resultados_newton[3][-1]:.10e}")
                print(f"  • Error relativo final: {resultados_newton[4][-1]:.10e}")
                print(f"  • Tiempo de ejecucion: {resultados_newton[5]:.6f} segundos")
                print(f"  • Uso de memoria: {resultados_newton[6]} bytes")
                
                # Preguntar si desea ver la tabla de iteraciones
                mostrar_tabla = entrada_segura("\n ¿Deseas ver la tabla de iteraciones? (s/n): ")
                if mostrar_tabla.lower() in ['s', 'si', 'sí', 'y', 'yes']:
                    # Imprimir tabla de iteraciones
                    imprimir_tabla("Newton-Raphson", resultados_newton[2], resultados_newton[3], 
                               resultados_newton[4], f)
            
            if method in ['3', '4']:  # Secante
                print("\n" + "-" * 70)
                print(" METODO DE LA SECANTE".center(70))
                print("-" * 70)
                
                print("\n EXPLICACION:")
                print("  El metodo de la secante requiere:")
                print("  • Dos valores iniciales x₀ y x₁")
                print("  • Converge rapido si los valores iniciales son buenos")
                print("  • Puede diverger si los valores iniciales estan lejos de la raiz")
                
                x0 = float(entrada_segura("\n Ingresa el primer valor inicial (x₀): "))
                x1 = float(entrada_segura("\n Ingresa el segundo valor inicial (x₁): "))
                
                print("\n Calculando solucion mediante la secante...")
                # Ejecutar el metodo de la secante
                resultados_secante = metodo_secante(f, x0, x1, tol, max_iter)
                
                print("\n RESULTADOS DEL METODO DE LA SECANTE:")
                print(f"  • Raiz encontrada: {resultados_secante[0]:.10f}")
                print(f"  • Valor de f(raiz): {f(resultados_secante[0]):.10e}")
                print(f"  • Iteraciones realizadas: {resultados_secante[1]}")
                print(f"  • Error absoluto final: {resultados_secante[3][-1]:.10e}")
                print(f"  • Error relativo final: {resultados_secante[4][-1]:.10e}")
                print(f"  • Tiempo de ejecucion: {resultados_secante[5]:.6f} segundos")
                print(f"  • Uso de memoria: {resultados_secante[6]} bytes")
                
                # Preguntar si desea ver la tabla de iteraciones
                mostrar_tabla = entrada_segura("\n ¿Deseas ver la tabla de iteraciones? (s/n): ")
                if mostrar_tabla.lower() in ['s', 'si', 'sí', 'y', 'yes']:
                    # Imprimir tabla de iteraciones
                    imprimir_tabla("Secante", resultados_secante[2], resultados_secante[3], 
                               resultados_secante[4], f)
            
            # Comparar metodos si se ejecutaron ambos
            if method == '4' and resultados_biseccion is not None and resultados_newton is not None and resultados_secante is not None:
                # Preguntar si desea ver la comparacion
                mostrar_comparacion = entrada_segura("\n ¿Deseas ver la comparacion entre los metodos? (s/n): ")
                if mostrar_comparacion.lower() in ['s', 'si', 'sí', 'y', 'yes']:
                    comparar_metodos(resultados_biseccion, resultados_newton, resultados_secante, ecuacion_str)
            
            # Preguntar si desea continuar
            continuar_opcion = entrada_segura("\n ¿Deseas resolver otra ecuacion? (s/n): ")
            if continuar_opcion.lower() not in ['s', 'si', 'sí', 'y', 'yes']:
                print("\n Gracias por usar el programa! ")
                break
                
        except Exception as e:
            print(f"\n ERROR: {e}")
            print("   Por favor, intenta de nuevo.")

def principal():
    """Funcion principal del programa."""
    # Configurar el parser de argumentos
    parser = argparse.ArgumentParser(description='Solucionador de ecuaciones no lineales')
    parser.add_argument('-e', '--equation', help='Ecuacion a resolver (f(x) = 0)')
    parser.add_argument('-m', '--method', choices=['biseccion', 'newton', 'secante', 'both'], 
                        help='Metodo a utilizar: biseccion, newton, secante o both')
    parser.add_argument('-a', type=float, help='Extremo izquierdo del intervalo para biseccion')
    parser.add_argument('-b', type=float, help='Extremo derecho del intervalo para biseccion')
    parser.add_argument('-x0', type=float, help='Aproximacion inicial para Newton-Raphson')
    parser.add_argument('-x1', type=float, help='Segundo valor inicial para la secante')
    parser.add_argument('-t', '--tolerance', type=float, default=1e-6, 
                        help='Tolerancia para el criterio de parada (default: 1e-6)')
    parser.add_argument('-i', '--max-iterations', type=int, default=100, 
                        help='Numero maximo de iteraciones (default: 100)')
    parser.add_argument('-f', '--file', help='Archivo con parametros de entrada')
    
    args = parser.parse_args()
    
    # Verificar si se proporcionaron argumentos para el modo no interactivo
    if args.equation and args.method:
        resolver_con_argumentos(args)
    elif args.file:
        # Leer parametros desde un archivo
        try:
            with open(args.file, 'r') as f:
                lines = f.readlines()
                
            # Parsear el archivo
            params = {}
            for line in lines:
                line = line.strip()
                if line and not line.startswith('#'):
                    key, value = line.split('=', 1)
                    params[key.strip()] = value.strip()
            
            # Configurar los argumentos
            if 'equation' in params:
                args.equation = params['equation']
            if 'method' in params:
                args.method = params['method']
            if 'a' in params:
                args.a = float(params['a'])
            if 'b' in params:
                args.b = float(params['b'])
            if 'x0' in params:
                args.x0 = float(params['x0'])
            if 'x1' in params:
                args.x1 = float(params['x1'])
            if 'tolerance' in params:
                args.tolerance = float(params['tolerance'])
            if 'max_iterations' in params:
                args.max_iterations = int(params['max_iterations'])
            
            # Verificar si se tienen los parametros minimos necesarios
            if args.equation and args.method:
                resolver_con_argumentos(args)
            else:
                print("Error: El archivo de parametros debe contener al menos 'equation' y 'method'.")
                sys.exit(1)
        except Exception as e:
            print(f"Error al leer el archivo de parametros: {e}")
            sys.exit(1)
    else:
        # Modo interactivo
        modo_interactivo()

if __name__ == "__main__":
    try:
        principal()
    except KeyboardInterrupt:
        print("\n\nPrograma interrumpido por el usuario.")
        sys.exit(0)
    except Exception as e:
        print(f"\nError inesperado: {e}")
        sys.exit(1)
