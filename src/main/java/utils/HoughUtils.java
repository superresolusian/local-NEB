package utils;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.measure.CurveFitter;
import ij.plugin.HyperStackConverter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.*;

import static java.lang.Math.*;
import static java.lang.Math.ceil;
import static utils.ArrayUtils.getIndexOfMax;
import static utils.ArrayUtils.getLimitArray;

public class HoughUtils {

    public static double[] getRadius(int maxRadius, int minRadius, double[] xVals, double[] yVals, double r2Threshold){

        // step 1: set up arrays
        double offset = getLimitArray(yVals, true);
        double amplitude = getLimitArray(yVals, false);
        double centre = (maxRadius+minRadius)/2;
        double sigma = (maxRadius+minRadius)/3;
        double[] params0 = new double[]{offset, amplitude, centre, sigma};

        Plot plot = new Plot("test", "radius", "hough", xVals, yVals);
        plot.show();

        // step 2: do fit
        CurveFitter cf = new CurveFitter(xVals, yVals);
        cf.setInitialParameters(params0);
        cf.doFit(CurveFitter.GAUSSIAN);

        double[] params = cf.getParams();
        double offset_ = params[0];
        double amplitude_ = params[1];
        double centre_ = params[2];
        double sigma_ = params[3];

        if(cf.getRSquared()<0.9) {
            if(cf.getRSquared()<r2Threshold){
                return new double[] {-1, cf.getRSquared()};
            }
            int ind = getIndexOfMax(yVals);
            return new double[] {xVals[ind], cf.getRSquared()};
        }

        return new double[] {centre_, cf.getRSquared()};
    }

    public static double getRadius(double[] xVals, double[] yVals, int localMaxRadius){
        int idxBest = -1;
        double radiusBest = -1;
        double thisMax = -1;
        for(int i=0; i<xVals.length; i++){
            double y = yVals[i];
            if(y>thisMax){
                thisMax = y;
                idxBest = i;
                radiusBest = xVals[i];
            }
        }

        double sum = 0;
        double wSum = 0;

        for(int i=max(idxBest-localMaxRadius, 0); i<=min(idxBest+localMaxRadius, yVals.length-1); i++){
            double val = yVals[i];
            double radius = xVals[i];
            sum += val;
            wSum += val*radius;
        }

        return wSum/sum;
    }

    public static int[] getOptimumRange(ImageStack ims, int nHoughs, int minRadius, int maxRadius, Polygon poly, double increment){

        int w = ims.getWidth();
        int h = ims.getHeight();

        double[] houghVals = new double[nHoughs];
        int[] xCs = poly.xpoints;
        int[] yCs = poly.ypoints;
        int nPoints = xCs.length;

        int averaging = 1;

        double[] radiusValues = new double[nHoughs];

        for(int i=0; i<nHoughs; i++){radiusValues[i] = minRadius+i*increment;} // as have decided two houghs per pixel

        for (int r = 1; r <= nHoughs; r++) {

            FloatProcessor fp = ims.getProcessor(r).convertToFloatProcessor();

            for(int n=0; n<nPoints; n++) {

                int xC = xCs[n];
                int yC = yCs[n];

                int counter = 0;
                double avVal = 0;

                for (int y = max(0, yC - averaging); y <= min(yC + averaging, h - 1); y++) {
                    for (int x = max(0, xC - averaging); x <= min(xC + averaging, w - 1); x++) {
                        avVal += fp.getf(x, y);
                        counter++;
                    }
                }

                avVal /= counter;
                houghVals[r - 1] += avVal;
            }

            houghVals[r-1]/=nPoints;

        }

        // fit Gaussian
        // step 1: set up arrays
        double offset = getLimitArray(houghVals, true);
        double amplitude = getLimitArray(houghVals, false);
        double centre = (maxRadius+minRadius)/2;
        double sigma = (maxRadius+minRadius)/3;
        double[] params0 = new double[]{offset, amplitude, centre, sigma};

        // step 2: do fit
        CurveFitter cf = new CurveFitter(radiusValues, houghVals);
        cf.setInitialParameters(params0);
        cf.doFit(CurveFitter.GAUSSIAN);

        double[] params = cf.getParams();
        double offset_ = params[0];
        double amplitude_ = params[1];
        double centre_ = params[2];
        double sigma_ = params[3];

        IJ.log("Fit parameters are offset="+offset_+", amplitude="+amplitude_+", centre="+centre_+", sigma="+sigma_);

        double idealRange = sigma_*3;
        double idealMin = max(1, centre_-(idealRange/2));
        double idealMax = min(w/2, centre_+(idealRange/2));

        return new int[]{(int) floor(idealMin), (int) ceil(idealMax)};

    }

    public static ImageStack getHoughAccumulators(ImageProcessor ip, int minRadius, int maxRadius, double increment, double threshold){
        int w = ip.getWidth();
        int h = ip.getHeight();

        float[] pixels = (float[]) ip.convertToFloatProcessor().getPixels();

        int wAcc = w + 2 * maxRadius;
        int hAcc = h + 2 * maxRadius;
        int nHoughs = (int) ((maxRadius - minRadius) * (1 / increment) + 1);

        ImageStack imsAccumulators = new ImageStack(wAcc, hAcc, nHoughs);

        for (int rad_ = 0; rad_ < nHoughs; rad_++) {
            int rad = (int) (minRadius + rad_ * increment);

            float[] accumulator = new float[wAcc * hAcc];

            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {

                    int p = x + y * w;

                    if (pixels[p] < threshold) continue;

                    for (int theta = 0; theta < 360; theta++) {

                        double a = x - rad * cos(toRadians(theta));
                        double b = y - rad * sin(toRadians(theta));

                        int a_ = (int) (a + maxRadius);
                        int b_ = (int) (b + maxRadius);

                        // weight by original pixel intensity
                        accumulator[a_ + b_ * wAcc] += pixels[p] / 360; // val?
                    }

                }

            }

            imsAccumulators.setProcessor(new FloatProcessor(wAcc, hAcc, accumulator), rad_ + 1);
        }

        imsAccumulators = imsAccumulators.crop(maxRadius - 1, maxRadius - 1, 0, w, h, nHoughs);

        return imsAccumulators;
    }

    public static Polygon peaksLocalMax(ImageProcessor ip, int minDistance, double threshold){

        int w = ip.getWidth();
        int h = ip.getHeight();

        ImageProcessor ipMax = ip.duplicate();
        Polygon maxima = new Polygon();

        //non maximum filter
        for(int y=minDistance; y<h-minDistance-1; y++){
            for(int x=minDistance; x<w-minDistance-1; x++){

                float maxVal = 0;

                for(int y_=y-minDistance; y_<=y+minDistance; y_++){
                    for(int x_=x-minDistance; x_<=x+minDistance; x_++){
                        maxVal = max(maxVal, ip.getf(x_,y_));
                    }
                }

                ipMax.setf(x, y, maxVal);
            }
        }

        for(int y=minDistance; y<h-minDistance-1; y++){
            for(int x=minDistance; x<w-minDistance-1; x++){
                float valMax = ipMax.getf(x,y);
                float valIm = ip.getf(x,y);
                if(valIm>threshold && valIm==valMax) maxima.addPoint(x, y);
            }
        }

        return maxima;

    }
}
