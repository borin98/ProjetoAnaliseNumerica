import numpy as np
from numpy.linalg import matrix_rank

class FatoracaoQRGS ( object ) :

    """
        Construtor do objeto

     Param : A -> matrix de entrada do script principal

    """
    def __init__( self , A ) :

        self._A = A
        self._m, self._n = A.shape
        self._rank = matrix_rank ( A )
        self._q = np.zeros ( shape = ( self._m, self._n ) )
        self._r = np.zeros ( shape = ( self._n, self._n ) )

    """
        
        Funcao Interna para normalizar os vetores
    
    """
    def _normVetores ( self ) :

        for i in range ( 0, self._n ) :

            self._q [ :, i ] /= np.sqrt ( np.dot ( self._q [ :, i ], self._q [ :, i ] ) )

    """

        Função Interna para fazer o calculo da projeção no somatorio

    """
    def _projecao( self, u, v ) :

        return ( u * np.dot ( v, u ) /np.dot ( u, u ) )

    """
    
        Funcao Interna para fazer a fatoracao GS e normalizada
    
    """
    def _fatoracaoGS ( self ) :

        self._q [ :, 0 ] = self._A [ :, 0 ]    # caso base de GS

        for i in range ( 1, self._n ) :

            for j in range ( 0, i ) :

                self._q [ :, i ] = self._A [ :, i ] - self._projecao ( self._q [ :, j ], self._A [ :, i ] )

        self._normVetores (  )

    """
    
        Função interna que faz a fatoração
        da matriz R
    
    """

    def _fatoracaoR ( self ) :

        for i in range ( 0, self._n ) :

            for j in range ( i, self._n ) :

                self._r [ i, j ] = np.dot ( self._q [ :, i ], self._A [ :, j ] )

    """

       Função utilizada para retornar a matriz Q com os vetores normalizados
       da fatoração QR utilizando GS

    """
    def matrixQR ( self ) :

         _, colunasA = self._A.shape

         if ( colunasA != self._rank ) :

            print("Nao é possivel fatorar a matriz")

            return np.array ( [] ), np.array ( [] )

         self._fatoracaoGS (  )

         self._fatoracaoR (  )

         return self._q, self._r