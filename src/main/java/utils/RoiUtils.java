package utils;

import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.frame.RoiManager;
import ij.process.FloatProcessor;

import java.awt.*;
import java.util.ArrayList;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static utils.ArrayUtils.unravel;

public class RoiUtils {

    final protected static int ROI_NAME = 0;
    final protected static int FRAME = 1;
    final protected static int CENTRE_IDX = 2;

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

    public static ArrayList<Roi> getRoisOfType(RoiManager rm, int roiType){
        ArrayList<Roi> roiList = new ArrayList<>();

        int nRois = rm.getCount();

        for(int i=0; i<nRois; i++){
            Roi thisRoi = rm.getRoi(i);
            int thisRoiType = thisRoi.getType();
            if(thisRoiType==roiType) roiList.add(thisRoi);
        }

        return roiList;
    }

    public static int getNumericPropertyFromCircleName(String name, int property){

        Matcher matcher = Pattern.compile("\\d+").matcher(name);

        ArrayList<Integer> allMatches = new ArrayList<>();
        while (matcher.find()) {
            allMatches.add(Integer.parseInt(matcher.group()));
        }

        return allMatches.get(property);
    }

    public static class RoiWithInfo {
        public Roi roi;
        public String name;
        public int roiNumber;
        public int frame;
        public int centre;

        public RoiWithInfo(Roi roi){
            this.roi = roi;
            this.name = roi.getName();
            this.roiNumber = getNumericPropertyFromCircleName(name, ROI_NAME);
            this.frame = getNumericPropertyFromCircleName(name, FRAME);
            this.centre = getNumericPropertyFromCircleName(name, CENTRE_IDX);
        }

    }

    public static double getMeanSignalFromRoi(FloatProcessor fp, Roi roi){
        Point[] points = roi.getContainedPoints();
        int nPoints = points.length;
        double meanVal = 0;
        for(int i=0; i<nPoints; i++){
            Point p = points[i];
            meanVal += fp.getf(p.x, p.y)/nPoints;
        }
        return meanVal;
    }

}
