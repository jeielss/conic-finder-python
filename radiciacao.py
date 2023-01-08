def f(c, n, x):
    return x**n - c

def df_exata(n,x):
    return n*(x**(n-1))
def df_prox(c,n,x_at,x_ant, tol=10e-10):
    dx = (x_at- x_ant)
    if(abs(dx) < tol): return "Erro"
    return (f(c,n,x_at) - f(c,n,x_ant))/dx
def metodo_bissecao_raiz(c,n, tol = 10e-10):
    #if(abs(c) < tol) return 0;
    if(c < 0) & (n % 2 != 0): return "Nao existe"
    inter = [-c, c]
    valor_fd, valor_fe = 0,0
    while(True):
        valor_f1 = f(c,n,inter[0])
        valor_f2 = f(c,n,inter[1])
        
        if(abs(valor_f1) < tol): return inter[0]
        elif(abs(valor_f2) < tol): return inter[1]
        dist = abs(inter[1] - inter[0])
        if(valor_f2 > 0):
            if(c > 0):
                inter[1] = inter[0] + dist/2
            else:
                inter[0] = inter[1] - dist/2
        else:
            if(c > 0):
                inter[0] = inter[1] - dist/2
            else:
                inter[1] = inter[0] + dist/2
        
    
       
def metodo_newton_raiz(x0, c, n, erro, itmax):
    cont_it = 0
    x_at = x0
    x_ant = x0 - 1
    er = abs((x_at-x_ant))
    while(cont_it < itmax) & (er >= erro):
        x_ant = x_at
        funcao = f(c,n,x_ant)
        derivada = df_exata(n,x_ant)
        if(abs(derivada) <= erro): return "Erro"
        x_at = x_ant - funcao/derivada
        er = abs((x_at-x_ant))
        cont_it+=1
    if(cont_it != itmax): return x_at
    else: return "Erro"

def metodo_secante_raiz(x0, c, n, erro, itmax):
    cont_it = 0
    x_at = x0
    x_ant = x0 - 1
    er = abs((x_at-x_ant))
    while(cont_it < itmax) & (er > erro):
       
        funcao = f(c,n,x_at)
        derivada = df_prox(c,n,x_at, x_ant)
        x_ant = x_at
        x_at = x_ant - funcao/derivada
        er = abs((x_at-x_ant))
        cont_it+=1
        print(x_ant)
    if(cont_it != itmax): return x_at
    else: return "Erro"

def raiz_n(x,n, tol=10e-10):
    return metodo_newton_raiz(x, x, n, tol, 1000)
def raiz_secante(x,n, tol=10e-10):
    return metodo_secante_raiz(x, x, n, tol, 1000)
