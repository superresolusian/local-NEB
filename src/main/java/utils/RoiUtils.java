package utils;

import ij.gui.PolygonRoi;
import ij.gui.Roi;

import java.awt.*;
import java.util.Random;

import static utils.ArrayUtils.unravel;

public class RoiUtils {

    private static Random random = new Random();

    public static Color[] getRandomColorPair(boolean opposite){
        float color1_r = random.nextFloat();
        float color1_g = random.nextFloat();
        float color1_b = random.nextFloat();

        Color color1 = new Color(color1_r, color1_g, color1_b);
        Color color2;
        if(opposite){
            color2 = new Color(1-color1_r, 1-color1_g, 1-color1_b);
        }
        else{
            float color2_r = random.nextFloat();
            float color2_g = random.nextFloat();
            float color2_b = random.nextFloat();

            color2 = new Color(color2_r, color2_g, color2_b);
        }

        Color[] colors = new Color[2];
        colors[0] = color1;
        colors[1] = color2;

        return colors;
    }

    public static PolygonRoi indicesToPolyline(int[] indices, int w, int xtranslate, int ytranslate){

        int npoints = indices.length;
        float[] xpoints = new float[npoints];
        float[] ypoints = new float[npoints];

        for(int i=0; i<npoints; i++){
            int[] xy = unravel(indices[i], w);
            xpoints[i] = (float) xy[0]+xtranslate;
            ypoints[i] = (float) xy[1]+ytranslate;
        }
        PolygonRoi polygonRoi = new PolygonRoi(xpoints, ypoints, Roi.POLYLINE);
        return polygonRoi;
    }

}
