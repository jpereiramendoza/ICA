/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ica;

/**
 *
 * @author jota
 */

import java.awt.Color;
import java.text.DecimalFormat;
import java.util.Random;


/**
 * Clase principal que realiza el llamado a FastICA y genera los distintos graficos
 * 
 * @author jota
 */
public class PrincipalICA 
{
    private double [][] W;  
    private double [] y;        
    private double[][]I;            
    private double eta = 0.01;
        
    public PrincipalICA(int N){        
       
    }
    
    
    /**
     * Tiene el numero de muestras que tiene nuestras observaciones
     */
    public static int nroMuestras = 5000;       
    
    
    /** 
     * Main principal. Ejecuta el programa
     * @param args 
     */
    public static void main(String args[])
    {
        
        // Obtenemos el objeto X que tiene las senales mezcladas 
        double [][]X = Data.X ; 
        
        // Solo lo tenemos como referencia para imprimir al final y realizar la comparacion.
        double [][]Sinit = Data.S;
        
        
        // Creamos el objeto FastICA con la matriz X de objetos  (Por simpleza la matriz esta 
        // traspuesta de manera que sea mucho mas facil agregar las lineas 
        FastICA fi = new FastICA(X);
        
        /**
         * Primero realizamos el centrado de la matriz X
         */
        fi.centrado();
        
        /** 
         * Y despues "blanqueamos la matriz 
         */
        fi.white();
        
        /**
         * Realizamos el calculo de ICA con FastICA. Los parametros corresponden 
         * a alpha , y el error que vamos a aceptar para determinar si el proceso 
         * concluyo de forma exitosa
         */
        double [][] W = fi.ICA( 1 , 0.000001 );
        
        
        /**
         * Por ultimo generamos la matriz S. Todos los datos generados y utilizados
         * para la ejecucion del metodo, se encuentran dentro de la clase, y se pueden
         * revisar con los getters implementados
         */
        fi.generateUnMixed();
        
        double [][]S = fi.getS();
        
        

        new Graficos("Matriz S Inicial  ", fi.generateTraspuesta(Sinit) , new Color(148 , 139, 61 ) );   
        //new Graficos("Matriz X Mezclada ",fi.generateTraspuesta(  fi.getX()) , new Color(218 , 124, 48 )  );
        new Graficos("Resultado final ", fi.getStrans() , new Color(57 , 106, 177 )  );   
    }
}

