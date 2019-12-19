/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ica;

import java.util.Random;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

/**
 * Clase que contiene todos los metodos necesarios para la implementacion del 
 * algoritmo FastICA
 *  
 *  Ademas contiene la matriz con los datos que son necesarios para realizar los 
 * calculos
 * 
 * @author jota
 */
public class FastICA 
{
    // Matriz W que es calculada 
    private double [][] W; 
    
    private double promedio[ ] ;
    //Corresponde al objeto de observaciones que tenemos 
    //
    //  Las observaciones son 
    //   Fila : Observacion del microfono X
    //   Columna : Observacion en el tiempo a 
    private double [][] X ; 
    private double [][] Xt ; 
    
    private double [][] Xcenter ; 
    
    private double [][] Xwhite ; 
    private double [][] XwhiteTras ; 
    
    public double [][] cov ;  
    
    // Matriz final con su traspuesta
    public double [][]S ;
    public double [][]Strans ;
    
    public double [][] Winv ;
    public double [][] WinvTras ;
    
    
    /**
     * Constructor del objeto. 
     */
    public FastICA(  double [][] X  )
    {
        this.X = X ;
        Xcenter = new double[ X.length][ X[0].length]  ;
        Xwhite = new double[ X.length][ X[0].length]  ;
        XwhiteTras = new double[ X[0].length][ X.length]  ;
        
        promedio = new double[ X[0].length] ;
        cov = new double[X[0].length][X[0].length];
        
        W = new double [ X[0].length] [ X[0].length];
        
        Winv = new double [ X[0].length] [ X[0].length];
        
        Xt = new double [X[0].length][X.length];
        for (int x=0; x < X.length; x++) 
        {
            for (int y=0; y < X[x].length; y++) 
            { 
                Xt[y][x] = X[x][y];
            }
        }
    }
    
    
    /**
     * Metodo que realiza el centrado de los datos.
     * 
     * Simplemente este metodo calcula el promedio por cada fila y posteriormente
     * lo va restando a cada fila de nuestra matriz.
     */
    public void centrado()
    {
       // Primero calculamos el promeido de cada fila
        for( int i = 0 ; i < this.X[i].length;i++)
        {
            for( int j = 0 ; j < X.length ; j++)
                promedio[i] += X[j][i] ;
            promedio[i] = promedio[i] / X.length;
            
        }   
        
        /** 
         * Ahora restamos el promedio de cada fila a nuestra matriz X
         */
        
        for( int i = 0 ; i < this.X[i].length;i++)
        {
            for( int j = 0 ; j < X.length ; j++)
                Xcenter[j][i]  =  X[j][i] - promedio[i]; 
            
        }
        
    }
    
    /**
     * Realiza el calculo de la covarianza.
     * 
     * El resultao de la matriz queda en el arreglo cov     
     */
    public void covarianza ()
    {
        double sum = 0 ;
        for( int i = 0 ; i < this.X[i].length;i++)
        {
            for( int j = 0 ; j < X.length ; j++)
            {
                sum += Math.pow( X[j][i] - promedio[i] , 2 );
            }
            cov[i][i] = sum / ( X.length - 1 ) ; 
            sum = 0 ;
        }  
        
        
        // Calculamos 1 con n
        for( int a = 0 ; a < X[0].length ; a ++  )
        {
            for( int i = (a+1)  ; i < this.X[i].length;i++)
            {
                for( int j = 0 ; j < X.length ; j++)
                {
                    sum +=  (X[j][a] - promedio[a] ) * (X[j][i] - promedio[i] );
                }
                cov[i][a] = sum / ( X.length - 1 ) ; 
                cov[a][i] = sum / ( X.length - 1 ) ; 
                sum = 0 ;
            }  
        }
        
    }   
    
    
    /**
     * Metodo que realiza el blanqueo de una matriz.
     * 
     * El proceso de blanqueo esta dado por la busqueda de un V calculado de la
     * siguiente manera
     * 
     *    V = E x D^(- 0.5) x Et
     * 
     * 
     */
    
    public void white()
    {
        int mues = X[0].length;
        covarianza();
       
       
        // COn la matriz de covarianza obtenemos los valores singulares  de la matriz
        RealMatrix rm = new Array2DRowRealMatrix( cov );
        /** 
         *  Realizamos la descomposicion de valores singulares
         * 
         * La descomposicion en valores singulares nos genera 3 matrices donde 
         * 
         *               A = U x S x Vt 
         * Y ademas
         * 
         *  Ut x U = I
         *  Vt x V = I 
         * 
         */ 
        
        SingularValueDecomposition svd = new SingularValueDecomposition( rm );
        RealMatrix Urm = svd.getU();
        RealMatrix Srm = svd.getS();
        RealMatrix Vrm = svd.getV();
        
        
        /** 
         * Supongamos que A es una matriz cuadrada. 
         * Un eigenvalues (valor propio) de A es un número r que cuando se 
         * resta de cada una de las entradas diagonales de A, se convierte 
         * A en una matriz singular. Restar un escalar r de cada entrada 
         * diagonal de A es lo mismo que restar r veces la matriz de identificación I de A.
         */
        
        double [][] eigenvalues = new double [ mues ][ mues ];
        for( int i =0; i<3 ; i++)
        {
            eigenvalues[i][i] = 1.0 / Math.sqrt(  Srm.getEntry(i, i) ) ;
        }
        
        
        // Obtenemos los datos de U y calculamos la traspuesta de U 
        double [][] Umt = Urm.getData();
        double [][] Umtt = new double [Umt.length][Umt[0].length ];
        for (int x=0; x < Umt.length; x++) 
        {
            for (int y=0; y < Umt[x].length; y++) 
            { 
                Umtt[y][x] = Umt[x][y];
            }
        }
        
        
        // Calculamos  E x D^(- 0.5) x Et    la matriz D el calculo interior esta en eigenvalues 
        double [][]  f1 = multiplyMatrices( eigenvalues , Umtt );
        double [][] whiteM = multiplyMatrices( Umt , f1  );

        
        
        // Por ultimo multiplicamos nuestra matriz de blanqueo con la matriz de datos X para 
        // obtener nuestra matriz blanqueada
        Xwhite = dotMatrices(whiteM, Xt);
        
        XwhiteTras = new double [ Xwhite[0].length][ Xwhite.length];
        for( int i = 0 ; i < Xwhite[0].length ; i++)
            for( int j = 0 ; j < Xwhite.length ; j++)
                XwhiteTras[i][j] = Xwhite[j][i];
    }
    
    
    /** 
     * Metodo principal que realiza el calculo del algoritmo FastICA.
     * 
     * 
     */
    public double [][] ICA( double  alpha  , double  lim )
    {
        // Almacenamos el nro de muestras que nos da la matriz X
        int nroMst = X.length;
        
        // Tenemos el nro de señales 
        int nroSig = W.length;
        
        // Lo primero es la generacion de la matriz W de forma aleatoria
        // Sabemos que W es una matriz cuadrada
        Random rd = new Random();
        for( int i = 0 ; i< W.length ; i ++)
            for( int j = 0 ; j< W.length ; j ++)
                W[i][j] = rd.nextDouble();
        
        
        // Recorremos por cada fila de la matriz W, que en realidad representa 
        // cuantas muestras tenemos ( es una matriz cuadrada de ponderacion)
        double [][]wTempNormal  = null ;
        double [][]wTempTrasp  = null ;
        for( int fila = 0  ; fila < nroSig; fila ++  )
        {
            
            double [] w = new double [ W.length ];
            double cuad = 0 ;
            
            // vector que almacena la nueva fila de W
            double[] wNew =  new double [ nroSig ] ;
            
            
            // Copiamos el vector FILA de la matriz W en el vector w
            // Ademas calculamos el valro de w segun 
            // w / sqrt( (w^2).sum())
            for( int q = 0 ; q < W.length ; q ++)
            {
                w[q] = W[fila][q];
                cuad +=w[q] * w[q] ; 
            }
            cuad = Math.sqrt( cuad ); 
            for( int q = 0 ; q < W.length ; q ++)
                w[q] = w[q]/cuad;
            
            
            // el contador de iteraciones para finalizar
            int iter = 0;
            //margen de error 
            double  err = 100 ; 
            
            // Ejecutamos el ciclo por 10mil veces o hasta que alcancemos un error bajo
            
            while(( iter < 10000) && ( err > lim) )
            {
                //producto punto entre la matriz w traspuesta y la señal blanqueada
                
                
                // Multiplicar ws por alpha y de ahi tangente hiperbolica con la transpuesta por cada elemento
                // de la matriz
                // por cada elemento del vector ws, calcula la tanh y la eleva al cuadrado, resta 1 y lo multiplica por alpha 
                double[] wg =  new double [ nroMst ] ;
                double[] wg_ =  new double [ nroMst ] ;
                
                double wg_Prom = 0;
                
                // Segun FastICA obtenemos 
                //  w <- E{ z * g( wt * z ) } - E{ g'(wt * z)} * w 
                //
                // Nuestro g es tanh 
                double [] ws = dotMatrizVector( w , Xwhite ); // Xwhite = z 
                for( int i = 0 ; i < nroMst ; i++ )
                {
                    wg[i] = Math.tanh( ws[i] * alpha ); 
                    
                    wg_[i] = (1 - Math.pow(Math.tanh(ws[i]) , 2 )) * alpha;
                    wg_Prom += wg_[i]; 
                }
                wg_Prom = wg_Prom / nroMst;
                // 
                
                // Se multiplica cada elemento del vector wg y la matriz de señales 
                // Esto es para calcular w / || w ||
                double [] temp = new double [nroSig];
                for( int i = 0 ; i < nroSig ; i ++ )
                {
                    for( int j = 0 ; j < nroMst ; j ++ )
                    {
                        temp[ i ] += wg[ j ] * Xwhite[ i ][ j ];
                    }
                    temp[ i ] = temp[ i ] / nroMst ; 
                    
                    wNew[ i ] = temp [ i ]  - wg_Prom * w[ i ]; 
                }
               
                
                // Aqui va la correlacion 
                // corr = np.dot(wNew, W[:c].T) 
                 // Es vector nuevo menos el producto punto de wNew y lo que llevamos acumulado de W y ese resultado otra vez por w 
                // Si es la primera fila no tenemos que calcular nada, porque el vector Winv esta vacio.
                double [][]corr;
                if( fila > 0 )
                {   
                    corr  =  multiplyMatrices( wNew , wTempTrasp ) ;
                    corr  =  multiplyMatrices(  corr  , wTempNormal ) ;
                    // print( corr );
                    
                    for ( int i  = 0 ; i < nroSig ; i ++)
                    {
                        wNew[i] = wNew[i] - corr[0][i];
                    }
                }
                
                
                // Tenemos el vector wNew que es la nueva W 
                cuad = 0 ;
                for( int q = 0 ; q < W.length ; q ++)
                {
                    cuad +=wNew[q] * wNew[q] ; 
                }
                cuad = Math.sqrt( cuad ); 
                for( int q = 0 ; q < W.length ; q ++)
                    wNew[q] = wNew[q]/cuad;
                
                
                // Falta calcular el limite si esta en cero
                // lim = np.abs(np.abs((wNew * w).sum()) - 1)
                double result = 0 ; 
                for( int i = 0 ; i < nroSig ; i ++ )
                {
                    result += wNew[i] * w[i];
                }
                
                err = Math.abs ( Math.abs( result ) - 1 );
                    
                
                
                // Copiamos wNew a w 
                for( int i = 0 ; i < nroSig ; i ++ )
                    w[i] = wNew[i];
                iter++;
                
                
            }//fin de while de iteraciones para encontrar w
            System.out.println( "Fila generada : " + fila + "    Nro de iteraciones : " + iter  );
            iter = 0 ; 
            //print( wNew );  
            /**
             * La matriz final W que tenemos es Winv que es la matriz W invertida
             * la que nos devuelve nuestra S
             */
            wTempTrasp = new double [ nroSig ][ fila +1 ]; 
            wTempNormal = new double [ fila +1 ][ nroSig  ]; 
                    
            for( int i = 0 ; i < nroSig ; i++)
            {
                Winv[ fila ][ i ] = wNew[ i ];
            }
            //Mantener una copia de wTemp en backup para la correlacion
            for( int i = 0 ; i < nroSig ; i++)
                for( int j = 0 ; j < (fila+1) ; j++)
                {
                    wTempTrasp[i][j] = Winv[j][i];
                    wTempNormal[j][i] = Winv[j][i];
                }
            
        }//Fin de for de recorrido de las distintas muestras
        
        
        //
        //Copiamos Winv a WinvTras que es la traspuesta 
        WinvTras = new double[ Winv[0].length ] [Winv.length];
        for(int i = 0 ; i < Winv[0].length; i ++  )
            for(int j = 0 ; j < Winv.length; j++  )
                WinvTras[i][j] = Winv[j][i];
        return Winv;
        
    }
    
    
    /**
     * Genera la matrix S con los elementos separados
     */
    public void generateUnMixed( )
    {
        
        double [][] un = multiplyMatrices(XwhiteTras, WinvTras );
        
        Strans  = new double [ un[ 0 ].length ][un.length];
        
        // Restamos el promedio 
        for(int i = 0 ; i < un.length ; i ++)
        {
            for( int j  = 0 ; j < un[0].length ; j++)
            {
                un[i][j] = un[i][j] - promedio[j];
                Strans[j][i] = un[i][j];
            }
        }
        
        S = un ;
        //print(un);
    }

    public double[][] getStrans() {
        return Strans;
    }

    public double[][] getW() {
        return W;
    }

    public void setW(double[][] W) {
        this.W = W;
    }

    public double[] getPromedio() {
        return promedio;
    }

    public void setPromedio(double[] promedio) {
        this.promedio = promedio;
    }

    public double[][] getX() {
        return X;
    }

    public void setX(double[][] X) {
        this.X = X;
    }

    public double[][] getXcenter() {
        return Xcenter;
    }

    public void setXcenter(double[][] Xcenter) {
        this.Xcenter = Xcenter;
    }

    public double[][] getXwhite() {
        return Xwhite;
    }

    public void setXwhite(double[][] Xwhite) {
        this.Xwhite = Xwhite;
    }

    public double[][] getXwhiteTras() {
        return XwhiteTras;
    }

    public void setXwhiteTras(double[][] XwhiteTras) {
        this.XwhiteTras = XwhiteTras;
    }

    public double[][] getCov() {
        return cov;
    }

    public void setCov(double[][] cov) {
        this.cov = cov;
    }

    public double[][] getS() {
        return S;
    }

    public void setS(double[][] S) {
        this.S = S;
    }

    public double[][] getWinv() {
        return Winv;
    }

    public void setWinv(double[][] Winv) {
        this.Winv = Winv;
    }

    public double[][] getWinvTras() {
        return WinvTras;
    }

    public void setWinvTras(double[][] WinvTras) {
        this.WinvTras = WinvTras;
    }
    
    
    
    
    /** 
     * Metodos para calculo de multiplicacion de matrices
     * @param firstMatrix
     * @param secondMatrix
     * @return 
     */
    
    double[][] multiplyMatrices(double[] firstMatrix, double[][] secondMatrix) 
    {
        double[][] result = new double[1][firstMatrix.length];
 
        for( int i  = 0 ; i < firstMatrix.length ; i ++ )
            result[0][i] = firstMatrix[i];
        return multiplyMatrices( result , secondMatrix) ;
        
    }
    
    double[][] multiplyMatrices(double[][] firstMatrix, double[][] secondMatrix) 
    {
        double[][] result = new double[firstMatrix.length][secondMatrix[0].length];
 
        for (int row = 0; row < result.length; row++) {
            for (int col = 0; col < result[row].length; col++) {
                result[row][col] = multiplyMatricesCell(firstMatrix, secondMatrix, row, col);
            }
        }
 
        return result;
    }
    double multiplyMatricesCell(double[][] firstMatrix, double[][] secondMatrix, int row, int col  ) {
        return multiplyMatricesCell( firstMatrix,  secondMatrix, row,  col, 0  );
    }
    
    double multiplyMatricesCell(double[][] firstMatrix, double[][] secondMatrix, int row, int col, int fore  ) {
        double cell = 0;
        if(fore == 0)
            fore = secondMatrix.length;
        for (int i = 0; i < fore; i++) {
            cell += firstMatrix[row][i] * secondMatrix[i][col];
        }
        return cell;
    }
 
    /** 
     * Calcula el producto punto de dos matrices
     * 
     * Asume que firstMatrix es n x n y que secondMatrix es m x n 
     * 
     * @param firstMatrix
     * @param secondMatrix
     * @return 
     */
    double[][] dotMatrices(double[][] firstMatrix, double[][] secondMatrix) 
    {
        double[][] result = new double[ secondMatrix.length ][ secondMatrix[0].length ];
        int contador = 0 ; 
        
        // recorremos las filas 
        for (int row = 0; row < secondMatrix.length; row++) 
        {
            for (int col = 0; col < secondMatrix[0].length; col++) 
            {
                // Ya sabemos que tenemos una fila y columna, ahora realizamos la multiplicacion
                //System.out.println( "Fila :" + row + " -  Col : " + col + " - cont :" + contador);
                result[row][col] = multiplyDot(firstMatrix, secondMatrix, row, col);
                
            }
        }
 
        return result;
    }
    
    double multiplyDot(double[][] firstMatrix, double[][] secondMatrix, int row, int col   ) 
    {
        double cell = 0;
        for (int i = 0; i < secondMatrix.length ; i++) {
            
            //int r = row%2;
            cell += firstMatrix[row][i] * secondMatrix[i][col];
        }
        return cell;
    }
    
    
    
    /** 
     * Producto pinto para un vector y una matriz
     * @param firstMatrix
     * @param secondMatrix
     * @return 
     */
    double[][] dotMatrices(double[] firstMatrix, double[][] secondMatrix) 
    {
        double[][] result = new double[ secondMatrix.length ][ secondMatrix[0].length ];
        int contador = 0 ; 
        
        // recorremos las filas 
        for (int row = 0; row < secondMatrix.length; row++) 
        {
            for (int col = 0; col < secondMatrix[0].length; col++) 
            {
                // Ya sabemos que tenemos una fila y columna, ahora realizamos la multiplicacion
                //System.out.println( "Fila :" + row + " -  Col : " + col + " - cont :" + contador);
                result[row][col] = multiplyDot(firstMatrix, secondMatrix, row, col);
                
            }
        }
 
        return result;
    }
    
    double multiplyDot(double[] firstMatrix, double[][] secondMatrix, int row, int col   ) 
    {
        double cell = 0;
        int vec = 0 ;

        for (int i = 0; i < secondMatrix.length ; i++) {
            
            if ( vec >= firstMatrix.length )
                vec = 0;
            cell += firstMatrix[vec] * secondMatrix[i][col];
            vec++;
        }
        return cell;
    }
    
    
    
    void print ( double ar[])
    {
        System.out.print( " [ " );
        for( int i = 0 ; i < ar.length ; i ++ )
            System.out.print( ar[i] + "  :   ");
        System.out.print( " ] \n" );
        
    }
    
    void print( double ar[][])
    {
        
        for( int i = 0 ; i < ar.length ; i ++ )
        {    
            System.out.print( " [ " );
            for( int j = 0 ; j < ar[i].length ; j ++ )
                System.out.print( ar[i][j] + "  :   ");
            System.out.print( " ] \n" );
        }
    }
    /** 
     * Metodo que calcula el producto punto entre un vector y una matriz
     */
    double[]  dotMatrizVector( double[] firstMatrix, double[][] secondMatrix ) 
    {
        double [] retorno = new double [ secondMatrix[0].length ];
        for( int i = 0 ; i < secondMatrix[0].length ; i++)
        {
            double elem = 0 ; 
            for( int j  = 0 ; j < firstMatrix.length ; j ++ )
            {
                elem += firstMatrix[j] * secondMatrix[j][i];
            }
            retorno[ i ] = elem;
        }
        return retorno;
    }
    
    /** 
     * Genera una matriz traspuesta
     * @param rec
     * @return 
     */
    public double [][] generateTraspuesta( double [][] rec) 
    {
        double [][]retorno   = new double [ rec[ 0 ].length ][rec.length];
        
        for(int i = 0 ; i < rec.length ; i ++)
        {
            for( int j  = 0 ; j < rec[0].length ; j++)
            {
                retorno[j][i] = rec[i][j];
            }
        }
        return retorno;
    }
}
