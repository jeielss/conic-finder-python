import matplotlib.pyplot as plt
import numpy as np
import radiciacao as raiz
import eliminacao_gaussiana as gauss


def determinante_O2(A):  #funcao que calcula determinante de ordem 2
    return A[0][0] * A[1][1] - A[1][0] * A[0][1]


def erro(v):
    if ("Erro" in str(v)):
        return str(v)
    return False


def f(C_q, C_l, T_i, x):
    return C_q * x * x + 2 * C_l * x + T_i


def df(C_q, C_l, x):  # derivada da funcao parabola
    return 2 * C_q * x + 2 * C_l


def metodo_newton_conica(
    x0,
    C1,
    C2,
    Ti,
    erro,
    itmax,
    div_zero=False
):  #metodo de newton que recebe coeficientes da conica e repassa para funcao e a derivada
    cont_it = 0
    x_at = x0
    x_ant = x0
    er = erro + 1
    while (cont_it < itmax) & (er >= erro):
        x_ant = x_at
        funcao = f(C1, C2, Ti, x_at)
        derivada = df(C1, C2, x_at)
        if (abs(derivada) <= erro):
            if (div_zero): return x_at
            else: return "Erro: Divisao por zero"

        x_at = x_at - funcao / derivada
        #if(abs(x_at - x_ant)-dx_temp < erro): return x_ant
        er = abs((x_at - x_ant))

        cont_it += 1

    if (cont_it != itmax): return x_at
    else:
        return "Erro: Max Iteracoes"


#funcao da conica, recebe todos os coeficientes, o termo independente e o x
#retorna 2 possiveis valores para y(esses valores podem ser identicos)
def f_conica(A, B, C, D, E, F, x, div_zero=False, tol=10e-14):
    result = A * x * x + 2 * D * x + F
    b = (B * x + E)
    chute = 0
    if (abs(C) > tol): chute = -b / C  #vertice
    elif (abs(b) > tol): chute = -result / 2 * b
    return [
        metodo_newton_conica(chute + 1, C, b, result, tol, 100, div_zero),
        metodo_newton_conica(chute - 1, C, b, result, tol, 100, div_zero)
    ]


def rotacao_ponto(ponto, rad):
    x, y = ponto[0], ponto[1]
    novo_x = x * np.cos(rad) - y * np.sin(rad)
    novo_y = x * np.sin(rad) + y * np.cos(rad)
    return novo_x, novo_y


def conica_sem_B(A, B, C, D, E, F):
    if (A == C): ang = np.pi / 4
    else: ang = np.arctan(2 * B / (A - C)) / 2
    cos = np.cos(ang)
    sin = np.sin(ang)
    _A = A * cos * cos + 2 * B * cos * sin + C * sin * sin
    _B = 0
    _C = A * sin * sin - 2 * B * sin * cos + C * cos * cos
    _D = (D * cos + E * sin)
    _E = (-D * sin + E * cos)
    return ang, _A, _B, _C, _D, _E, F


def dist_entre_pontos(ponto1, ponto2):
    cateto1 = ponto1[0] - ponto2[0]
    cateto2 = ponto1[1] - ponto2[1]
    return raiz.raiz_n(cateto1 * cateto1 + cateto2 * cateto2, 2)


def plot_con(coeficientes,
             tipo_con,
             pontos=[],
             estimativas=[],
             p_graf=1000,
             tol=10e-14):

    A, B, C, D, E, F = coeficientes
    #Gerando o grafico da conica
    #O grafico e feito em diversas curvas, ja que alguns valores de x podem ter 2 valores correspondentes em y e vice-versa
    x_c, x_b, y_c, y_b, t = [], [], [], [], [
    ]  #criando as variaveis para armazenar os pontos

    c_q = (B * B - A * C)
    c_l = (2 * B * E - 2 * D * C) / 2
    t_i = (E * E - C * F)
    fator_x0 = 0
    if (abs(c_q) > tol): fator_x0 = -c_l / c_q
    elif (abs(c_q) > tol): fator_x0 = -t_i / (2 * c_l)
    encontro = [
        metodo_newton_conica(fator_x0 - 1, c_q, c_l, t_i, tol, 1000,
                             True), "Erro",
        metodo_newton_conica(fator_x0 + 1, c_q, c_l, t_i, tol, 1000, True),
        "Erro"
    ]

    if not erro(encontro[0]):
        encontro[1] = f_conica(A, B, C, D, E, F, encontro[0], True, tol)
        encontro[1] = encontro[1][0]
    if not erro(encontro[2]):
        encontro[3] = f_conica(A, B, C, D, E, F, encontro[2], True, tol)
        encontro[3] = encontro[3][0]
    fig, ax = plt.subplots()
    tam = pontos.shape[0]
    #print("Encontro:", encontro, "\nFator_x0:", fator_x0)
    #print("%.10fx² + %.10fxy + %.10fy² + %.10fx + %.10fy + %.10f = 0" % (A,2*B,C,2*D,2*E,F))
    angulo, _A, _B, _C, _D, _E, _F = conica_sem_B(A, B, C, D, E, F)

    #definindo dominio para pegar a curva
    if (tipo_con == "Parabola"):
        if (abs(_C) <= tol):
            vertice = (-_D / _A, ((-_D * _D) / _A + _F) / (-2 * _E))
            foco_orbitado = (vertice[0], vertice[1] - _E / (2 * _A))
            distancia_focal = abs(_E / (2 * _A))
        else:
            vertice = (((-_E * _E) / _C + _F) / (-2 * _D), -_E / _C)
            foco_orbitado = (vertice[0] - _D / (2 * _C), vertice[1])
            distancia_focal = abs(_D / (2 * _C))
        vertice = rotacao_ponto(vertice, angulo)
        foco_orbitado = rotacao_ponto(foco_orbitado, angulo)
        dom_x = [
            foco_orbitado[0] - 6 * distancia_focal,
            foco_orbitado[0] + 6 * distancia_focal
        ]
        dom_y = [
            foco_orbitado[1] - 6 * distancia_focal,
            foco_orbitado[1] + 6 * distancia_focal
        ]
    else:
        centro = (-_D / (_A), -_E / (_C))
        termo = (-F + _D * _D / _A + _E * _E / _C)
        fatores_a_b = (termo / _C, termo / _A)
        foco_orbitado = centro
        foco_2 = centro
        a = raiz.raiz_n(abs(fatores_a_b[0]), 2)
        b = raiz.raiz_n(abs(fatores_a_b[1]), 2)
        if (tipo_con == "Hiperbole"): c = raiz.raiz_n(a * a + b * b, 2)
        elif (tipo_con == "Elipse"):
            if (a > b):
                tmp = b
                b = a
                a = tmp
            c = raiz.raiz_n(-a * a + b * b, 2)

        lado = -1
        cont_vezes_perto_foco = 0
        if (fatores_a_b[0] > fatores_a_b[1]):
            for p in range(tam):
                x, y = rotacao_ponto((pontos[p, 0], pontos[p, 1]), -angulo)
                if (y > centro[1]): cont_vezes_perto_foco += 1
            if (cont_vezes_perto_foco > tam / 2): lado = 1
            foco_orbitado = (centro[0], centro[1] + c * lado)
            foco_2 = (centro[0], centro[1] - c * lado)

        else:
            for p in range(tam):
                x, y = rotacao_ponto((pontos[p, 0], pontos[p, 1]), -angulo)
                if (x > centro[0]): cont_vezes_perto_foco += 1
            if (cont_vezes_perto_foco > tam / 2): lado = 1
            foco_orbitado = (centro[0] + c * lado, centro[1])
            foco_2 = (centro[0] - c * lado, centro[1])

        foco_orbitado = rotacao_ponto(foco_orbitado, angulo)
        foco_2 = rotacao_ponto(foco_2, angulo)
        centro = rotacao_ponto(centro, angulo)
        if (tipo_con == "Hiperbole"):
            menor = pontos[0, 0]
            maior = pontos[0, 0]
            for p in range(tam):
                x = abs(pontos[p, 0])
                y = abs(pontos[p, 1])
                if (abs(x - foco_orbitado[0]) > maior):
                    maior = abs(x - foco_orbitado[0])
            if (maior < 1): maior += 1
            dom_x = [foco_orbitado[0] - maior, foco_orbitado[0] + maior]
            dom_y = [foco_orbitado[1] - maior, foco_orbitado[1] + maior]
        else:
            if (a > b): maior = a * 1.1
            else: maior = b * 1.1
            dom_x = [centro[0] - (maior), centro[0] + (maior)]
            dom_y = [centro[1] - (maior), centro[1] + (maior)]

    t = np.linspace(
        dom_x[0], dom_x[1], p_graf
    )  #Lista com "p_graf" pontos dentro do intervalo que foi informado pelo usuario

    dist_ponto = 1.1 * abs(t[0] - t[1])
    cont = 0
    for i in t:
        i_ant = t[cont - 1]
        if (not erro(encontro[1])):
            if (encontro[0] > i_ant) and (encontro[0] < i):
                t = np.insert(t, cont, encontro[0])
                cont += 1
        if (not erro(encontro[3])):
            if (encontro[2] > i_ant) and (encontro[2] < i):
                t = np.insert(t, cont, encontro[2])
            cont += 1
    i_ant = t[0]
    for i in t:  #loop para aplicar a funcao da conica em cada ponto de t

        if (tol > 10e-6): tol_graf = 10e-6
        else: tol_graf = tol
        valores_de_f = f_conica(
            A, B, C, D, E, F, i, 1, tol_graf
        )  #funcao da conica que retorna dois valores para o x informado
        ponto_y0 = (i, valores_de_f[0])
        ponto_y1 = (i, valores_de_f[1])
        if not erro(
                valores_de_f):  #os pontos validos sao adicionados as listas
            if (
                    abs(i - i_ant) >= dist_ponto
            ):  # caso tenha uma diferença entre pontos maior que 1, ele imprime uma parte curva
                ax.plot(x_c, y_c, color='b')
                ax.plot(x_b, y_b, color='b')
                x_c, x_b, y_c, y_b = [], [], [], []
            if (tipo_con != "Hiperbole"):
                x_c += [i]
                y_c += [valores_de_f[0]]
                x_b += [i]
                y_b += [valores_de_f[1]]
                i_ant = i
            else:
                if (dist_entre_pontos(foco_orbitado, ponto_y0) <=
                        dist_entre_pontos(foco_2, ponto_y0)):
                    x_c += [i]
                    y_c += [valores_de_f[0]]
                    i_ant = i
                if (dist_entre_pontos(foco_orbitado, ponto_y1) <=
                        dist_entre_pontos(foco_2, ponto_y1)):
                    x_b += [i]
                    y_b += [valores_de_f[1]]
                    i_ant = i

    ax.plot(x_c, y_c, color='b')
    ax.plot(x_b, y_b, color='b')
    cor = 'g'
    for p in range(tam):
        x = pontos[p, 0]
        y = pontos[p, 1]
        ax.plot(x, y, marker='o', color=cor, linewidth=0.5)
        if (cor == 'g'): cor = 'y'
        elif (cor == 'y'): cor = 'b'
        else: cor = 'g'
        if (p == 0): ax.annotate("    Inicio", xy=(x, y))
        elif (p == tam - 1): ax.annotate("    Fim", xy=(x, y))
    distancia_x_1dia = pontos[tam - 1, 0] - pontos[tam - 2, 0]

    for p in range(estimativas.shape[0]):
        x = estimativas[p, 0]
        y = estimativas[p, 1]
        ax.plot(x, y, marker='x', color='r')
        if (p == 0): ax.annotate("Previsoes \n\n", xy=(x, y))
    #ax.plot(centro[0],centro[1], marker = 'o')
    #ax.annotate(" Centro\n  (%.2f, %.2f)\n" %(centro), xy=(centro))
    #ax.plot(vertice[0], vertice[1], marker = 'o')
    #ax.plot(ponto0[0], ponto0[1], marker = 'o')
    ax.plot(foco_orbitado[0], foco_orbitado[1], marker='o')
    ax.plot(foco_2[0], foco_2[1], marker='o')
    #ax.set(xlim=(centro[0] -10, centro[0] + 10),ylim=(centro[1] - 10, centro[1] + 10))

    ax.set(xlim=(dom_x), ylim=(dom_y))
    plt.show()
    return estimativas
def orbita(A, mostrar=True, tol=10e-4):
    pontos = np.matrix(A)
    ordem = pontos.shape[0]
    matriz = np.zeros((ordem, 5))
    b = ordem * [1]
    for i in range(ordem):
        matriz[i, 0] = pontos[i, 0] * pontos[i, 0]
        matriz[i, 1] = pontos[i, 0] * pontos[i, 1] * 2
        matriz[i, 2] = pontos[i, 1] * pontos[i, 1]
        matriz[i, 3] = pontos[i, 0] * 2
        matriz[i, 4] = pontos[i, 1] * 2

    transposta = np.transpose(matriz)
    matriz = np.dot(transposta, matriz)
    b = np.dot(transposta, b)
    matriz_estendida = np.hstack((matriz, np.atleast_2d(b).T))
    res = gauss.escalonamento(matriz_estendida)
    res += [-1]
    A, B, C, D, E, F = res
    tipo_con = determinante_O2([[A, B], [B, C]])
    if (abs(tipo_con) < tol): tipo_con = "Parabola"
    elif (tipo_con < 0): tipo_con = "Hiperbole"
    else: tipo_con = "Elipse"

    estimativas = np.ndarray((3, 2), dtype=float)
    tam = pontos.shape[0]
    distancia_x_1dia = (pontos[tam - 1, 0] - pontos[tam - 6, 0]) / 5
    for n in range(1, 4):
        x = pontos[tam - 1, 0] + n * distancia_x_1dia
        y = f_conica(A, B, C, D, E, F, x)
        if (erro(y[0]) and erro(y[1])):
            x = pontos[tam - 1, 0] - n * distancia_x_1dia
            y = f_conica(A, B, C, D, E, F, x)
        if (not erro(y)):
            ponto_y0 = (x, y[0])
            ponto_y1 = (x, y[1])
            ponto_final = (pontos[tam - 1, 0], pontos[tam - 1, 1])
            sentido_y = pontos[tam - 1, 1] - pontos[tam - 2, 1]
            dist_y0 = ponto_y0[1] - ponto_final[1]
            dist_y1 = ponto_y1[1] - ponto_final[1]

            if (dist_y0 * dist_y1 > 0):
                if (dist_entre_pontos(ponto_y0, ponto_final) <
                        dist_entre_pontos(ponto_y1, ponto_final)):
                    y = y[0]
                else:
                    y = y[1]
            elif (dist_y0 * sentido_y > 0):
                y = y[0]
            else:
                y = y[1]

        elif (not erro(y[0])):
            y = y[0]
        else:
            y = y[1]

        estimativas[n - 1][0] = x
        estimativas[n - 1][1] = y
    if (mostrar):
        plot_con(res, tipo_con, pontos, estimativas, 10000, tol)
        ang, res[0], res[1], res[2], res[3], res[4], res[5] = conica_sem_B(
            A, B, C, D, E, F)
        #plot_con(res, tipo_con,pontos,1000,[-10,10],0.001*tol)
    if (tipo_con == "Parabola"): tipo_con = "parabolica"
    elif (tipo_con == "Elipse"): tipo_con = "eliptica"
    else: tipo_con = "hiperbolica"
    return [tipo_con, estimativas]


print("\nElipses")
EL = np.load("elipticas.npy", allow_pickle=True)
acertos = 0
erros = 0
for A in range(0, 1000):
    if orbita(EL[A], mostrar=0)[0] == "eliptica":
        acertos += 1
    else:
        erros += 1
print("Erros: " + str(erros))
print("Acertos: " + str(acertos))

print("\nHiperboles")
EL=np.load("hiperbolicas.npy",allow_pickle=True)
acertos=0
erros=0
for A in range(0,1000):
    if orbita(EL[A],mostrar=0)[0]=="hiperbolica":
        acertos+=1
    else:
        erros+=1
print("Erros: " + str(erros))
print("Acertos: " + str(acertos)) 

print("\nParabolas")
EL=np.load("parabolicas.npy",allow_pickle=True)
acertos=0
erros=0
for A in range(0,1000):
    if orbita(EL[A],mostrar=0)[0]=="parabolica":
        acertos+=1
    else:
        erros+=1
print("Erros: " + str(erros))
print("Acertos: " + str(acertos))
