package utils;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.EllipseRoi;
import ij.gui.PolygonRoi;
import ij.plugin.frame.RoiManager;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

import java.awt.*;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import static java.lang.Math.abs;
import static java.lang.Math.toRadians;
import static utils.HoughUtils.*;
import static utils.ImageUtils.*;
import static utils.PixelPathUtils.*;
import static utils.RoiUtils.getRandomColorPair;

public class TrackerForSingleCrop {

    static final double rad30deg = toRadians(30);

    private ImageStack ims;
    private int w, h, nFrames;

    private double blurSize;
    private boolean invertLUT;
    private int maxFrameGap;
    private double exclusionFraction;

    private int minRadius, maxRadius;
    private double increment;
    private double thresholdModifier;

    private Color[] colors = getRandomColorPair(true);

    private boolean showHough, showSegmentation;
    public ImagePlus impHough, impSegmentations;
    public ImageStack imsHough, imsSegmentations;

    public RoiManager rm;
    int strokeWidth;

    public void setImageStack(ImageStack ims){
        this.ims = ims;
        this.w = ims.getWidth();
        this.h = ims.getHeight();
        this.nFrames = ims.getSize();
    }

    public void setTrackerOptions(double blurSize, boolean invertLUT, int maxFrameGap, double exclusionFraction) {
        this.blurSize = blurSize;
        this.invertLUT = invertLUT;
        this.maxFrameGap = maxFrameGap;
        this.exclusionFraction = exclusionFraction;
    }

    public void setHoughOptions(int minRadius, int maxRadius, double increment, double thresholdModifier){
        this.minRadius = minRadius;
        this.maxRadius = maxRadius;
        this.increment = increment;
        this.thresholdModifier = thresholdModifier;
    }

    public void setRoiOptions(RoiManager rm, int strokeWidth){
        this.rm = rm;
        this.strokeWidth = strokeWidth;
    }

    public void setDebugOptions(boolean showHough, boolean showSegmentation){
        this.showHough = showHough;
        this.showSegmentation = showSegmentation;
    }

    public void setupPreview(){
        this.imsHough = new ImageStack(w, h);
        this.impHough = new ImagePlus("Projected Hough transforms");
        this.imsSegmentations = new ImageStack(w, h);
        this.impSegmentations = new ImagePlus("Segmentations");
    }

    public void updatePreview(boolean isHough, ImageProcessor ip, int n){
        if (isHough) {
            imsHough.addSlice(ip);
            if (n == 1) {
                impHough.setStack(imsHough);
                impHough.show();
            } else {
                impHough.setSlice(n);
            }
        }
        else {

        }
    }


    public void doTracker(){

        setupPreview();

        //// determine division angle
        ImageStack imsNormalised = normaliseToMax(ims);
        FloatProcessor fpRoiAverageIntensity = averageProject(imsNormalised);
        ImageProcessor ipSkeletonAverage = getSkeletonImage(fpRoiAverageIntensity, blurSize, invertLUT);

        short[] pixelsSegmentedAverage = connectedComponents(ipSkeletonAverage);
        ShortProcessor spSegmentedAverage = new ShortProcessor(w, h);
        spSegmentedAverage.setPixels(pixelsSegmentedAverage);

        new ImagePlus("segmented average", spSegmentedAverage).show();

        int nSegmentsAverage = (int) getMaxOfImage(spSegmentedAverage);

        // get largest segment
        int bestArea = 0;
        int biggestSegment = 0;

        ArrayList<PixelPathUtils.regionProps> propsList = new ArrayList<>();

        for(int s=0; s<nSegmentsAverage; s++) {
            // get properties for each segment
            PixelPathUtils.regionProps thisProps = getProps(s + 1, spSegmentedAverage);

            // check size of segment
            int area = thisProps.area;
            if(area>bestArea){
                bestArea = area;
                biggestSegment = s;
            }

            propsList.add(thisProps);
        }

        PixelPathUtils.regionProps bestProps = propsList.get(biggestSegment);
        ArrayList<int[]> endpoints = bestProps.endpointsList;
        int nEndpoints = endpoints.size();

        double optimumAngle = 0;

        if(nEndpoints==2){
            int[] ep1 = endpoints.get(0);
            int[] ep2 = endpoints.get(1);
            optimumAngle = getAngleBetweenPoints(ep1, ep2);
            IJ.log("Angle between endpoints is "+Math.toDegrees(optimumAngle));
        }
        else{
            ArrayList<int[]> endpointsPruned = pruneEndpoints(bestProps);
            if(endpointsPruned.size()<2){
                IJ.log("Can't handle this skeleton");
                return; //TODO: continue for >1 Roi
            }
            int[] ep1 = endpointsPruned.get(0);
            int[] ep2 = endpointsPruned.get(1);
            optimumAngle = getAngleBetweenPoints(ep1, ep2);
            IJ.log("Angle between endpoints (pruned) is "+Math.toDegrees(optimumAngle));
        }

        boolean doContinueCounter = false;
        int continueCounter = 0;

        int[] bestCentre1 = new int[0], bestCentre2 = new int[0];
        int[] previousCentre1 = new int[0], previousCentre2 = new int[0];

        double bestRadius1 = 0, bestRadius2 = 0;
        double previousRadius1 = 0, previousRadius2 = 0;

        for (int n = 1; n <= nFrames; n++) {
            if (continueCounter > maxFrameGap) {
                IJ.log("max frame gap reached");
                break;
            }

            //// get circles in each frame
            ImageProcessor ipRoi1 = ims.getProcessor(n).duplicate();
            double threshold = getThresholdForIp(ipRoi1);
            IJ.log("image threshold is " + threshold);

            ImageStack imsAccumulators = getHoughAccumulators(ipRoi1, minRadius, maxRadius, increment, threshold);
            FloatProcessor fpMaxAccumulator = maxProject(imsAccumulators);
            double houghThreshold = getThresholdForIp(fpMaxAccumulator);
            IJ.log("hough threshold is " + threshold);

            int nHoughs = imsAccumulators.getSize();

            if(showHough) updatePreview(true, fpMaxAccumulator, n);

            Polygon maxima = peaksLocalMax(fpMaxAccumulator, minRadius, houghThreshold * thresholdModifier);
            int nMaxima = maxima.npoints;
            IJ.log("Frame " + n + ", nMaxima = " + nMaxima);

            if (nMaxima < 2) {
                if (doContinueCounter) continueCounter++;
                continue;
            }

            //get radii; crude
            ArrayList<int[]> maximaCoordinates = new ArrayList<int[]>();
            LinkedHashMap<int[], Double> radii = new LinkedHashMap<int[], Double>();

            double[] radiusVals = new double[nHoughs];
            for(int i=0; i<nHoughs; i++) radiusVals[i] = minRadius + increment*i;


            for (int i = 0; i < nMaxima; i++) {
                int thisX = maxima.xpoints[i];
                int thisY = maxima.ypoints[i];
                maximaCoordinates.add(new int[]{thisX, thisY});

                double[] houghVals = getHoughValsFromPeakSingleChannel(imsAccumulators, thisX, thisY, nHoughs, 1);
                double bestRadius = getRadius(radiusVals, houghVals, 2);
                IJ.log("centre at ("+thisX+","+thisY+") has radius "+bestRadius);

                radii.put(maximaCoordinates.get(i), bestRadius);
            }


            // see if good pair exists
            int[] bestPair = new int[2];
            boolean foundGoodPair = false;
            double bestAngleDifference = Double.MAX_VALUE;

            for (int i = 0; i < nMaxima - 1; i++) {
                int[] testCentre1 = maximaCoordinates.get(i);
                double r1 = radii.get(testCentre1);

                for (int j = i + 1; j < nMaxima; j++) {
                    int[] testCentre2 = maximaCoordinates.get(j);
                    double r2 = radii.get(testCentre2);

                    if (getDistance(testCentre1, testCentre2) < (r1+r2)) {
                        IJ.log("centres are too close!");
                        continue; // too close
                    }

                    // test if this pair of maxima have satisfactory angle
                    double angle = getAngleBetweenPoints(testCentre1, testCentre2);
                    if ((angle + rad30deg) < optimumAngle || (angle - rad30deg) > optimumAngle) continue;

                    foundGoodPair = true;
                    if (abs(angle - optimumAngle) < bestAngleDifference) {
                        bestAngleDifference = abs(angle - optimumAngle);
                        bestPair[0] = i;
                        bestPair[1] = j;
                    }
                }
            }

            if (!foundGoodPair) {
                if (doContinueCounter) continueCounter++;
                IJ.log("didn't find a good pair");
                continue;
            }

            if (doContinueCounter) {
                previousCentre1 = bestCentre1;
                previousCentre2 = bestCentre2;
                previousRadius1 = bestRadius1;
                previousRadius2 = bestRadius2;
            }

            bestCentre1 = maximaCoordinates.get(bestPair[0]);
            bestCentre2 = maximaCoordinates.get(bestPair[1]);

            if (!doContinueCounter) {
                // convention: nucleus 1 is left, nucleus 2 is right
                boolean centre1IsLeft = bestCentre1[0] < bestCentre2[0];
                if (!centre1IsLeft) {
                    // switch nuclei to keep direction constant
                    int[] tempBestCentre1 = bestCentre1;
                    int[] tempBestCentre2 = bestCentre2;
                    bestCentre1 = tempBestCentre2;
                    bestCentre2 = tempBestCentre1;
                }
            }

            if (doContinueCounter) {
                double dist11 = getDistance(previousCentre1, bestCentre1);
                double dist21 = getDistance(previousCentre2, bestCentre1);
                double dist12 = getDistance(previousCentre1, bestCentre2);
                double dist22 = getDistance(previousCentre2, bestCentre2);

                boolean condition1 = dist11 < previousRadius1*2 && dist22 < previousRadius2*2;
                boolean condition2 = dist12 < previousRadius1*2 && dist21 < previousRadius2*2;

                if (!condition1 && !condition2) {
                    IJ.log("Nuclei positions deviated too much");
                    continueCounter++;
                    continue;
                }

                if (condition2) {
                    // switch nuclei to keep direction constant
                    int[] tempBestCentre1 = bestCentre1;
                    int[] tempBestCentre2 = bestCentre2;
                    bestCentre1 = tempBestCentre2;
                    bestCentre2 = tempBestCentre1;
                }
            }

            doContinueCounter = true;
            continueCounter = 0;

            bestRadius1 = radii.get(bestCentre1);
            bestRadius2 = radii.get(bestCentre2);

            IJ.log(("best centre 1 is at (" + bestCentre1[0] + "," + bestCentre1[1] + "), radius = "+bestRadius1));
            IJ.log(("best centre 2 is at (" + bestCentre2[0] + "," + bestCentre2[1] + "), radius = "+bestRadius2));

            EllipseRoi ellipseRoi1 = new EllipseRoi(bestCentre1[0] - bestRadius1, bestCentre1[1] - bestRadius1, bestCentre1[0] + bestRadius1, bestCentre1[1] + bestRadius1, 1);
            ellipseRoi1.setPosition(n);
            ellipseRoi1.setStrokeColor(colors[0]);
            ellipseRoi1.setName("frame " + n + " centre 1");
            rm.addRoi(ellipseRoi1);

            EllipseRoi ellipseRoi2 = new EllipseRoi(bestCentre2[0] - bestRadius2, bestCentre2[1] - bestRadius2, bestCentre2[0] + bestRadius2, bestCentre2[1] + bestRadius2, 1);
            ellipseRoi2.setPosition(n);
            ellipseRoi2.setStrokeColor(colors[1]);
            ellipseRoi2.setName("frame " + n + " centre 2");
            rm.addRoi(ellipseRoi2);

            double radius1exc = bestRadius1 * exclusionFraction;
            double radius2exc = bestRadius2 * exclusionFraction;

            // get skeletons for this frame
            ImageProcessor ipSkeleton = getSkeletonImage(ipRoi1, blurSize, invertLUT);
            short[] pixelsSegmented = connectedComponents(ipSkeleton);
            ShortProcessor spSegmented = new ShortProcessor(w, h);
            spSegmented.setPixels(pixelsSegmented);

            if(showSegmentation) updatePreview(false, spSegmented, n);

            // find skeletons with endpoints inside circles
            int[] bridge = getBridge(spSegmented, bestCentre1, bestRadius1, bestCentre2, bestRadius2, 2, bestRadius1+bestRadius2);

            int[] bridgeNoCircle = pathExcludeCircle(bridge, bestCentre1, radius1exc, bestCentre2, radius2exc, w);

            PolygonRoi polygonRoi = indicesToPolyline(bridgeNoCircle, w);
            polygonRoi.setStrokeWidth(strokeWidth);
            polygonRoi.setPosition(n);
            polygonRoi.setStrokeColor(colors[0]);
            polygonRoi.setName("Frame " + n + " bridge");
            rm.addRoi(polygonRoi);

    }
}

}
