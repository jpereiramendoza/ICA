
package ica;

/**
 *
 * @author 
 */

import java.awt.*;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import javax.swing.*;
import org.jfree.chart.*;
import org.jfree.chart.plot.*;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.*;

public class Graficos 
{    
    private String name;           
    private JFrame frame;   
    
    private int count = 0 ; 
    
    private Color color ; 
    
    private Color c[] = { Color.DARK_GRAY ,Color.BLUE , Color.GREEN };
    public Graficos(String name, double y[][] , Color color  ) 
    {
        
        this.color = color;
        this.name = name;  
        
        frame = new JFrame( name );                                                        
        
        JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));        
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);        
        
        for(int row=0; row< y.length; row++)
        {           
            XYDataset d = createDataset( y ,  row);
            JFreeChart chart = createGrafico(d, row);
            ChartPanel chartPanel = new ChartPanel(chart);        
            chartPanel.setPreferredSize( new Dimension(800,100) );
            panel.add(chartPanel);                
        }  
        
        frame.getContentPane().add(new JScrollPane(panel));        
        frame.setSize(1070, 520);        
        frame.setVisible(true);

        frame.addComponentListener(new ComponentAdapter() {
            public void componentResized(ComponentEvent componentEvent) {
                System.out.println( frame.getWidth() + " : " + frame.getHeight());
            }
        });        
    }
            
    
    private XYDataset createDataset(double[][] y, int row) 
    {
        XYSeries series = new XYSeries("");
        for(int i=0; i<y[0].length;i++)
        {
            series.add(i,y[row][i]);   
        }        
        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series);                                    
        return dataset;        
    }
        
    private XYDataset createDataset(double[][][] y, int i, int j) {
        XYSeries series = new XYSeries("");
        int K = y[0][0].length;
        for(int k=0; k<K;k++){
            series.add(k,y[i][j][k]);   
        }        
        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series);                                    
        return dataset;      
    }
    
    private JFreeChart createGrafico(final XYDataset dataset, int row) {
               
        final JFreeChart chart = ChartFactory.createXYLineChart(
            name + String.valueOf( row + 1 ),                         //chart title
            "",                        // x axis label
            "",        // y axis label
            dataset,                    // data
            PlotOrientation.VERTICAL,
            false,                      // include legend
            false,                      // tooltips
            false                       // urls
        );  
        
        chart.getTitle().setFont(new Font("Arial", Font.PLAIN, 20));
        
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer( );
        renderer.setSeriesPaint( 0 , color );
        renderer.setSeriesStroke( 0 , new BasicStroke( 1.0f ) );
        renderer.setSeriesShapesVisible(0, false);
        
        chart.setBackgroundPaint(Color.WHITE);        
        XYPlot plot = chart.getXYPlot();  
        
        plot.setBackgroundPaint(Color.WHITE);             
        plot.setDomainGridlinePaint(Color.LIGHT_GRAY);
        plot.setRangeGridlinePaint(Color.LIGHT_GRAY);         
        plot.setRenderer( renderer );
        
        return chart;        
    }
}
