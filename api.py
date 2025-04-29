from flask import Flask, request, jsonify
from solucionador_ecuaciones import analizar_ecuacion, metodo_biseccion, metodo_newton_raphson, metodo_secante, comparar_metodos, crear_graficas_comparativas
from flask_cors import CORS
import traceback
import math

app = Flask(__name__)
CORS(app)

def limpiar_json(obj):
    if isinstance(obj, dict):
        return {k: limpiar_json(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [limpiar_json(v) for v in obj]
    elif isinstance(obj, float):
        if math.isinf(obj) or math.isnan(obj):
            return None
        return obj
    else:
        return obj

@app.route('/resolver', methods=['POST'])
def resolver():
    data = request.json
    try:
        ecuacion = data.get('ecuacion')
        metodo = data.get('metodo')
        tol = float(data.get('tol', 1e-6))
        max_iter = int(data.get('max_iter', 100))
        
        if not ecuacion or not metodo:
            return jsonify({'error': 'Faltan parámetros requeridos: ecuacion y metodo'}), 400

        f = analizar_ecuacion(ecuacion)

        # Ejecutar el método seleccionado
        if metodo == 'biseccion':
            a = float(data.get('a'))
            b = float(data.get('b'))
            resultado = metodo_biseccion(f, a, b, tol, max_iter)
            raiz, iteraciones, valores, errores_abs, errores_rel, tiempo, memoria = resultado
        elif metodo == 'newton-raphson':
            x0 = float(data.get('x0'))
            resultado = metodo_newton_raphson(f, x0, tol, max_iter)
            raiz, iteraciones, valores, errores_abs, errores_rel, tiempo, memoria = resultado
        elif metodo == 'secante':
            x0 = float(data.get('x0'))
            x1 = float(data.get('x1'))
            resultado = metodo_secante(f, x0, x1, tol, max_iter)
            raiz, iteraciones, valores, errores_abs, errores_rel, tiempo, memoria = resultado
        else:
            return jsonify({'error': 'Método no soportado. Usa biseccion, newton-raphson o secante.'}), 400

        # Determinar convergencia
        convergio = False
        mensaje_estado = ""
        if errores_abs:
            if errores_abs[-1] < tol:
                convergio = True
                mensaje_estado = "Convergió correctamente."
            else:
                mensaje_estado = "No convergió: máximo de iteraciones alcanzado o error no menor que tolerancia."
        else:
            mensaje_estado = "No se pudo calcular el error."

        # Resumen final
        respuesta = {
            'resumen': {
                'metodo': metodo,
                'ecuacion': ecuacion,
                'raiz': raiz,
                'f_raiz': f(raiz) if raiz is not None else None,
                'iteraciones': iteraciones,
                'error_absoluto_final': errores_abs[-1] if errores_abs else None,
                'error_relativo_final': errores_rel[-1] if errores_rel else None,
                'convergio': convergio,
                'mensaje_estado': mensaje_estado,
                'tiempo_ejecucion_seg': tiempo,
                'memoria_bytes': memoria,
            },
            'valores_intermedios': valores,
            'errores_absolutos': errores_abs,
            'errores_relativos': errores_rel,
            'detalle_iteraciones': [
                {
                    'iteracion': i + 1,
                    'valor': valores[i],
                    'f_x': f(valores[i]),
                    'error_absoluto': errores_abs[i] if i < len(errores_abs) else None,
                    'error_relativo': errores_rel[i] if i < len(errores_rel) else None,
                } for i in range(len(valores))
            ],
            'comparar_metodos_result': None,
            'graficas_generadas': [],
        }
        return jsonify(limpiar_json(respuesta))
    except Exception as e:
        return jsonify({'error': str(e), 'trace': traceback.format_exc()}), 500

@app.route('/')
def home():
    return 'API de solucionador de ecuaciones funcionando', 200

if __name__ == '__main__':
    app.run(debug=True, port=5000)
