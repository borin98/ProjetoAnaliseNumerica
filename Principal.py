import numpy as np
import random

from FatoracaoQR import FatoracaoQRGS   # biblioteca necessaria para rodar o algorítmo
from numpy.linalg import inv, LinAlgError

"""
    
    Funçoes que faz a verificacao do condicionamento

"""

def cond2 ( k1, A, normaA ) :

    linhasA, colunasA = A.shape

    w = np.zeros ( shape = ( linhasA, 1 ) )
    for i in range(0, linhasA):

        w [ i ] [ 0 ] = random.uniform (0, 1000 * linhasA)


    print("\nValor K(A)1 = {}\n".format(k1))

    v = np.dot(inv(A), w)

    for k in range ( 0, linhasA ) :

        soma = abs ( w [ k ] [ 0 ] )

    nomrAW = normaMatriz1 ( v )

    print("Valor A*A^(-1)*w/w = {}\n".format(normaA*nomrAW/soma))


def cond1 ( k1, A ) :

    linhasA, colunasA = A.shape

    for i in range ( 0, colunasA - 1 ) :

        print("\nValor K(A)1 = {}\n".format ( k1 ) )

        for k in range ( 0, colunasA ) :

            soma1 = abs ( A [ i ] [ k ] )
            soma2 = abs ( A [ i + 1 ] [ k ] )

        print("Valor A{}/A{} = {}\n".format( i, i + 1, soma1/soma2 ) )

"""

    Função que inicializa as colunas com valores aleatorios

"""

def numeroAle( n ) :

    U = np.zeros ( shape = ( n, 1 ) )

    for i in range ( 0, n ) :

        ale = random.uniform ( 0, 10*n )
        U [ i ] = ale

    return U

"""
    
    Função que cria a norma 1 de uma matriz

"""
def normaMatriz1 ( A ) :

    linhasA, colunasA = A.shape

    Max = np.zeros ( shape = ( linhasA, 1 ) )

    for i in range ( 0, linhasA ) :

        soma = 0

        for j in range ( 0, colunasA ) :

            soma += A [ i ] [ j ]

        Max [ i ] = soma

    return Max.max (  )

"""
    
    Função que cria a matriz identidade

"""
def matrizIdentidade ( n ) :

    id = np.zeros ( shape = ( n, n ) )

    for i in range ( 0, n ) :

        id [ i ] [ i ] = 1

    return id

"""

    Função que inicializa a matriz de input

"""
def criaA ( n,m ) :

    A = np.zeros ( shape = ( n, m ) )

    print("\n")

    for i in range ( 0, n ) :

        for j in range ( 0, m ) :

            A [ i ] [ j ] = ( float ) ( input ( "Digite o valor do elemento da linha {} coluna {} : ".format ( i, j ) ) )

        print("\n")

    return A

"""

    Função que cria a matriz 

"""
def criaB ( n ) :

    B = np.zeros(shape = ( n, 1 ) )

    print("\n")

    for i in range(0, n):

        B [ i ] [ 0 ] = ( float ) ( input("Digite o valor do elemento da linha {} : ".format ( i ) ) )

    print("\n")

    return B

"""
    
    Função que cria a matriz de hilbert

"""
def matrizHilbert ( n ) :

    H = np.zeros( shape = ( n, n ) )

    for i in range ( 0, n ) :

        for j in range ( 0, n ) :

            H [ i ] [ j ] = 1/ ( ( i + 1 ) + ( j + 1 ) - 1 )

    return H

"""

    Algorítmo que faz a retrosubstituição do sistema 

"""
def retroSubs ( A, B ) :

    n, _ = B.shape

    x = np.zeros ( shape = ( n, 1 ) )
    x [ n - 1 ] = B [ n - 1 ] / A [ n - 1 ] [ n - 1 ] # caso base

    for i in range ( n - 2, -1, -1 ) :

        soma = B [ i ]

        for j in range ( i + 1, n, 1 ) :

            soma -= A [ i ] [ j ] * x [ j ]

        x  [ i ] = soma / A [ i ] [ i ]

    return x

def main (warning=None) :

    op = ( int ) ( input ("Deseja criar uma matriz de hilbert ? ( 0 para nao e 1 para sim ) : ") )
    ale = ( int ) ( input ( "Deseja criar uma matriz A e B com valores aleatorios ? ( 0 para nao e 1 para sim ) : " ) )

    while ( ale != 0 and ale != 1 ) :

        print("\nValor inválido !!!\n")
        ale = (int)(input("Deseja criar uma matriz A e B com valores aleatorios ? ( 0 para nao e 1 para sim ) : "))

    while ( op != 1 and op != 0 ) :

        print("\nValor inválido !!!\n")
        op = (int)(input("Deseja criar uma matriz de hilbert ? ( 0 para nao e 1 para sim ) : "))

    n = ( int ) ( input ( "\nDigite um valor de n : " ) )

    if ( op == 0 ) :

        while ( n <= 0 ) :

            print ( "\nValor inválido !!!\n" )
            n = ( int ) ( input ( "Digite um valor de n : " ) )

        m = ( int ) ( input ( "Digite um valor valor de m : " ) )

        while ( m <= 0 ) :

            print ( "\nValor inválido !!!\n" )
            m = ( int ) ( input ( "Digite um valor de n : " ) )

        if ( ale == 0 ) :

            A = criaA ( n, m )

        else :

            A = np.zeros ( shape = ( n, m ) )

            for i in range ( 0, n ) :

                for j in range ( 0, m ) :

                    A [ i ] [ j ] = random.uniform ( 0, 1000*n )

    else :

        A = matrizHilbert ( n )

    print("\n---------- Vetor B -----------\n")

    if ( ale == 0 ) :

        B = criaB ( n )

    else :

        B = np.zeros ( shape = ( n, 1 ) )

        for i in range ( 0, n ) :

            B [ i ] [ 0 ] = random.uniform ( 0, 1000*n )

    FQRGS = FatoracaoQRGS ( A )

    Q, R = FQRGS.matrixQR (  )

    # realiza as operacoes somente se for possivel fazer a fatoracao QR por GS
    if ( Q.size != 0 ) :

        B = np.dot(Q.T, B)
        x = retroSubs(R, B)

        print("matriz Q : \n{}\n".format(Q))
        print("matriz R : \n{}\n".format(R))
        print("vetor x de resposta : \n{}\n".format(x))

        print("----------------------------------------------------------------\n")

        print("Fazendo Ax = B, com B = Id {} x {} e x = ( Q^T )( R^(-1) )\n")

        id = matrizIdentidade(n)
        print("Matriz identidade do problema : \n{}\n".format(id))

        x = np.dot( inv ( R ), Q.T )
        print("Produto x = ( R^(-1) )( Q^T ) ) : \n{}\n".format(x))

        prod = np.dot(A, x)
        print("Produto Ax : \n{}\n".format(prod))

        print("----------------------------------------------------------------\n")

        try :

            inv ( A )

        except LinAlgError :

            print ( "A matriz A não é invertível !!!!\n" )

            # no caso em que a matrix contem somente zeros
            if ( not np.any ( A ) ) :

                warning.warn ( "A matriz A é nula.\n Terminando o programa" )
                return 0

            warning.warn("Seu condicionamento será calculado fazendo A' = A^(T)*A :")

            A = np.dot ( A.T, A )

        normA = normaMatriz1(A)
        print("valor da norma de A : {}\n".format(normA))

        normInvA = normaMatriz1(inv(A))
        print("valor da norma de A^(-1) : {}\n".format(normInvA))

        print("----------------------------------------------------------------\n")

        cond = normA * normInvA
        print("valor de condicionamento do sistema : {}\n".format(cond))

        print("----------------------------------------------------------------\n")
        print("Comparando os valores de condicionamento\n")
        cond1 ( cond, A )

        print("----------------------------------------------------------------\n")
        print("Comparando os valores de condicionamento com valores de w aleatorio\n")
        cond2 ( cond, A, normA )


    else :

        print ( "Não é possível fazer a fatoração QR via QS .\nTerminando o programa" )
        return 0


if __name__ == '__main__':

    main (  )