package utils;

import java.awt.*;

public class ArrayUtils {

    public static double getLimitArray(double[] array, boolean isMin){
        int n = array.length;
        if(isMin){
            double min = Double.POSITIVE_INFINITY;
            for(int i=0; i<n; i++){
                min = Math.min(array[i], min);
            }
            return min;
        }
        double max = Double.NEGATIVE_INFINITY;
        for(int i=0; i<n; i++){
            max = Math.max(array[i], max);
        }
        return max;
    }

    public static int getIndexOfMax(double[] array){
        int nElements = array.length;
        int bestInd = -1;
        double bestMax = Double.NEGATIVE_INFINITY;

        for(int i=0; i<nElements; i++){
            if(array[i]>bestMax){
                bestMax = array[i];
                bestInd = i;
            }
        }
        return bestInd;
    }

    public static float[] intToFloat(int[] array){
        int n = array.length;
        float[] floatArray = new float[n];
        for(int i=0; i<n; i++) floatArray[i] = (float) array[i];
        return  floatArray;
    }

    public static int ravel(int x, int y, int w) {
        return y * w + x;
    }

    public static int[] unravel(int ind, int w) { return new int[]{ind % w, ind / w}; }


}
