import numpy as np
import time
#inicio das funcoes para a eliminacao gaussiana
def subt_linhas(A, l1, l2, quo,tol=10e-14):
    #Realiza a subtracao: l1 - quo * l2 (elemento a elemento)
    for j in range(A.shape[1]):
        A[l1,j] -= quo * (A[l2,j])
        if(abs(A[l1,j]) < tol): A[l1,j] = 0 #Checa se o numero e maior que a tolerancia, senao for torna-o 0
def troca_linhas(A, l1, l2): #Funcao que troca de lugar l1 e l2
    cp_l1 = A[l1].copy()
    cp_l2 = A[l2].copy()
    
    A[l1] = cp_l2
    A[l2] = cp_l1
def checar_linhas_nulas(A,tol=10e-14):
    l_nul = 0 #variavel para quantidade de linhas nulas
    t = A.shape
    #Dois loops para checar todos os termos da matriz linha a linha
    for i in range(t[0]):
        ele_n_nulo = 0 #variavel para quantidade de elemento nao nulos
        for j in range(t[1] - 1): #Nao checa o termo independente
            if(abs(A[i,j]) > tol): ele_n_nulo = 1
        if(ele_n_nulo == 0): l_nul += 1
    return l_nul
def retroSubs_dp(A, chute=1, tol=10e-14):
    t = A.shape
    l_nul = checar_linhas_nulas(A) #fator de correcao aplicado para evitar erros com linhas nulas
    variaveis = (t[0] - l_nul)*[0] #lista que armazena os valores das variaveis
    
    variaveis_att= t[0]*[0] #lista que armazena quais variaveis da matriz foram atualizadas
    variaveis = t[0]*[0] #lista que armazena o valor das variaveis
    #dois loops para acessar os valores da matriz escalonada do final para o comeco, exceto a coluna do termo independente
    for i in reversed(range(t[0] - l_nul)):
        for j in reversed(range(t[1] - 1)):
            coeficiente = A[i,j] #obtendo o coeficiente da variavel
            if(abs(coeficiente) > tol): #Se o coeficiente for diferente de 0 
                if((i != j)): #e se nao for o coeficiente da propria variavel
                    if(variaveis_att[j] == 0): variaveis[j] = chute; #chuta o valor "1" para uma variavel que ainda nao foi atualizada, mas que outra variavel depende de seu valor
                    variaveis[i] -= variaveis[j]*coeficiente #como e um valor dependente de outra variavel realiza a multiplicao e subtrai o valor
                variaveis_att[j] = 1 #definindo que a variavel foi atualizada
        variaveis[i] += A[i, t[1] - 1] #somando o valor independente
                
    return variaveis
def retroSubs(A, tmp,tol=10e-14):
    t = A.shape
    retorno = [];
    l_nul = checar_linhas_nulas(A) #fator de correcao aplicado para evitar erros com linhas nulas
    variaveis = (t[0] - l_nul)*[0] #lista que armazena os valores das variaveis
    if(t[0] - l_nul != t[1] - 1): #Verificando se temos linhas nao nulas suficientes para obter o valor das variaveis
        #Realizando retrosubstituicao incompleta
        variaveis_dp = A.copy() 
        variaveis_att = t[1] * [0] #lista que armazena quais variaveis da matriz foram atualizadas
        #dois loops para acessar os valores da matriz escalonada do final para o comeco
        for i in reversed(range(t[0] - l_nul)): 
            for j in reversed(range(t[1])):
                variaveis_dp[i, j] = A[i,j]/A[i,i] #dividindo a linha pelo valor do item da linha que esta na diagonal principal
                if(variaveis_att[j] == 1): # Se a variavel da coluna ja foi atualizada, subtrai o valor que representa do termo independente e torna o seu coeficiente 0
                    for v in reversed(range(t[1])):
                        if(v == j): variaveis_dp[i, v] = 0 #tornando o coeficiente 0
                        elif(abs(variaveis_dp[i,j]) > tol): 
                            a = variaveis_dp[j, v] #obtendo o valor da variavel
                            b = variaveis_dp[i,j] #coeficiente da variavel na linha atual
                            variaveis_dp[i, v] -= a*b #subtrai o valor do termo independente
            variaveis_att[i] = 1 #guardando a informacao de que foi atualizada
        return variaveis_dp
    else:
        #dois loops para acessar os valores da matriz escalonada do final para o comeco
        for i in reversed(range(t[0] - l_nul)):
            for j in reversed(range(t[1])):
                if(j == (t[1] - 1)): variaveis[i] += A[i,j] #caso seja o valor independente soma no valor da variavel
                elif(j!=i): variaveis[i] -= variaveis[j] * A[i,j] #caso seja um valor dependente de outra variavel realiza a multiplicao e subtrai o valor
            
            
            variaveis[i] /= A[i,i] #divide pelo coficiente que multiplica a variavel para obter o valor dela
            if(abs(variaveis[i]) <= tol): variaveis[i] = 0 # caso a variavel esteja com valor abaixo da tolerancia torna-o 0
        for v in range(t[0] - l_nul):
            retorno += [variaveis[v]] 
        #print("Metodo Gaussiano - Tempo de execucao: " + str(time.time() - tmp) + "s") #Impressao do tempo de execucao
        return retorno
                
def escalonamento(A,tol=10e-14):
    te = time.time(); # Marca o inicio do tempo do processo de escalonamento
    A = np.matrix(A)
    t = A.shape #Guarda o tamanho da matriz
    pivo = 0
    l_pivo = -1 #variavel para linha do pivo atual
    for i in range(t[0]): #loop para passar por todas as linhas
        if(t[1] > i):pivo = A[i,i]
        cont = i + 1
        while (abs(pivo) <= tol) & (cont < t[0]): #loop para encotrar o pivo
            A[i,i] = 0; #definindo pra 0 ja que menor que tolerancia
            troca_linhas(A, i, cont)
            if(t[1] > i): pivo = A[i,i]
            cont+=1
        l_pivo = i
        
        #loop para subtrair a linha do pivo das outras linhas de baixo, de modo que os elementos abaixo do pivo sejam zerados
        if(t[1] > i):
            for l in range(i,t[0]): 
                if(l != i) & (abs(pivo) > tol): subt_linhas(A, l, i, (A[l, l_pivo]/ pivo)) 
    
    # Codigo pra imprimir a matriz escalonada
    #print("Matriz Escalonada:")
    #print(A, end = "\n")
    
    return retroSubs(A,te,tol) #Retornando o resultado da função de retrosubstituicao

#Fim da eliminacao gaussiana