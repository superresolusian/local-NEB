package gui;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.*;
import ij.plugin.ChannelSplitter;
import ij.plugin.frame.RoiManager;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import utils.PixelPathUtils;

import java.awt.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Random;

import static java.lang.Math.*;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static utils.HoughUtils.*;
import static utils.ImageUtils.*;
import static utils.ImageUtils.getThresholdForIp;
import static utils.ImageUtils.maxProject;
import static utils.PixelPathUtils.*;
import static utils.PixelPathUtils.getAngleBetweenPoints;
import static utils.PixelPathUtils.getDistance;
import static utils.RoiUtils.getRandomColorPair;

public class TestTrackerCrop_ extends BaseGUI_ {


    int minRadius, maxRadius, strokeWidth;
    double minLength, maxGap, maxFrameGap, linkingDistance, thresholdModifier, exclusionFraction, blurSize;
    boolean showHough, showSegmentation, invertLUT, doFit;

    int nFrames;

    Random random = new Random();
    RoiManager rm = new RoiManager();
    ChannelSplitter CS = new ChannelSplitter();

    static final double rad30deg = toRadians(30);

    @Override
    public boolean beforeSetupDialog(String arg) {
        autoOpenImp = true;
        return true;
    }

    @Override
    public void setupDialog() {
        gd = new NonBlockingGenericDialog("Bridge tracking");
        gd.addMessage("~~~ Chough parameters ~~~");
        gd.addNumericField("Minimum radius (pixels)", getPrefs("minRadius", 5), 0);
        gd.addNumericField("Maximum radius (pixels)", getPrefs("maxRadius", 20), 0);
        gd.addNumericField("Threshold modifier (leave at 0.75 unless low SNR)", getPrefs("thresholdModifier", 0.75), 2);
        gd.addMessage("~~~ Miscellaneous ~~~");
        gd.addNumericField("Blur size", getPrefs("blurSize", 2), 0);
        gd.addNumericField("Line width for plotting (pixels)", getPrefs("strokeWidth", 3), 0);
        gd.addNumericField("Maximum empty frame gap", getPrefs("maxFrameGap", 3), 0);
        gd.addNumericField("Bridge exclusion fraction relative to radius", getPrefs("exclusionFraction", 1.0), 2);
        gd.addCheckbox("Do fit", getPrefs("doFit", false));
        gd.addCheckbox("Show Hough transforms", getPrefs("showHough", false));
        gd.addCheckbox("Show segmentation", getPrefs("showSegmentation", false));
        gd.addCheckbox("Invert LUT for skeletonization", getPrefs("invertLUT", true));
    }

    @Override
    public boolean loadSettings() {
        minRadius = (int) gd.getNextNumber();
        maxRadius = (int) gd.getNextNumber();
        thresholdModifier = gd.getNextNumber();

        blurSize = gd.getNextNumber();
        strokeWidth = (int) gd.getNextNumber();
        maxFrameGap = gd.getNextNumber();
        exclusionFraction = gd.getNextNumber();

        doFit = gd.getNextBoolean();
        showHough = gd.getNextBoolean();
        showSegmentation = gd.getNextBoolean();
        invertLUT = gd.getNextBoolean();


        setPrefs("minRadius", minRadius);
        setPrefs("maxRadius", maxRadius);
        setPrefs("thresholdModifier", thresholdModifier);

        setPrefs("blurSize", blurSize);
        setPrefs("strokeWidth", strokeWidth);
        setPrefs("maxFrameGap", maxFrameGap);
        setPrefs("exclusionFraction", exclusionFraction);

        setPrefs("doFit", doFit);
        setPrefs("showHough", showHough);
        setPrefs("showSegmentation", showSegmentation);
        setPrefs("invertLUT", invertLUT);

        return true;
    }

    @Override
    public void execute() throws InterruptedException, IOException {
        // image preparation
        nFrames = imp.getNFrames();
        nChannels = imp.getNChannels();

        ImageStack imsChannel1, imsChannel2 = new ImageStack();

        if(nChannels==2) {
            imsChannel1 = CS.getChannel(imp, 1);
            imsChannel2 = CS.getChannel(imp, 2);
        }
        else imsChannel1 = imp.getImageStack();

        RoiManager thisManager = new RoiManager().getInstance();

        if(thisManager==null){
            throw new IOException("Expecting ROI manager containing ROIs!");
        }

        rm = thisManager;

        int nRois = rm.getCount();

        if(nRois==0){
            throw new IOException("Expecting ROI manager containing ROIs!");
        }

        double increment = 1;

        // measure background
        Roi roiBG = rm.getRoi(0);
        rm.rename(0, "background");
        double[] backgroundCh1 = getBackground(imsChannel1, roiBG);
        double[] backgroundCh2 = new double[nFrames];
        if(nChannels==2) backgroundCh2 = getBackground(imsChannel2, roiBG);


        Color[] colors = getRandomColorPair(true);
        int w = imp.getWidth();
        int h = imp.getHeight();

        //// determine division angle
        ImageStack imsNormalised = normaliseToMax(imsChannel1);
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

        ImageStack imsProjectedHoughs = new ImageStack(w, h);
        ImagePlus impProjectedHoughs = new ImagePlus("Projected Houghs");

        ImageStack imsSegmentations = new ImageStack(w, h);
        ImagePlus impSegmentations = new ImagePlus("Segmentations");

        for (int n = 1; n <= nFrames; n++) {
            if (continueCounter > maxFrameGap) {
                IJ.log("max frame gap reached");
                break;
            }

            //// get circles in each frame
            ImageProcessor ipRoi1 = imsChannel1.getProcessor(n).duplicate();
            double threshold = getThresholdForIp(ipRoi1);
            IJ.log("image threshold is " + threshold);

            ImageStack imsAccumulators = getHoughAccumulators(ipRoi1, minRadius, maxRadius, increment, threshold);
            FloatProcessor fpMaxAccumulator = maxProject(imsAccumulators);
            double houghThreshold = getThresholdForIp(fpMaxAccumulator);
            IJ.log("hough threshold is " + threshold);

            int nHoughs = imsAccumulators.getSize();

            if (showHough) {
                imsProjectedHoughs.addSlice(fpMaxAccumulator);
                if (n == 1) {
                    impProjectedHoughs.setStack(imsProjectedHoughs);
                    impProjectedHoughs.show();
                } else {
                    impProjectedHoughs.setSlice(n);
                }
            }

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

            IJ.log(("best centre 1 is at (" + bestCentre1[0] + "," + bestCentre1[1] + ")"));
            IJ.log(("best centre 2 is at (" + bestCentre2[0] + "," + bestCentre2[1] + ")"));

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

            if (showSegmentation) {
                imsSegmentations.addSlice(spSegmented);
                if (n == 1) {
                    impSegmentations.setStack(imsSegmentations);
                    impSegmentations.show();
                } else {
                    impSegmentations.setSlice(n);
                }
            }

            // find skeletons with endpoints inside circles
            ArrayList<PixelPathUtils.regionProps> frameProps = new ArrayList<>();
            int nSegmentsFrame = (int) getMaxOfImage(spSegmented);

            for(int s=0; s<nSegmentsFrame; s++) frameProps.add(getProps(s+1, spSegmented));


            double angle = getAngleBetweenPoints(bestCentre1, bestCentre2);
            Line line = new Line(bestCentre1[0] + radius1exc * cos(angle), bestCentre1[1] - radius1exc * sin(angle),
                    bestCentre2[0] - radius2exc * cos(angle), bestCentre2[1] + radius2exc * sin(angle));

            line.setStrokeWidth(strokeWidth);
            Polygon poly = line.getPolygon();

            line = new Line(bestCentre1[0] + radius1exc * cos(angle), bestCentre1[1] - radius1exc * sin(angle),
                    bestCentre2[0] - radius2exc * cos(angle), bestCentre2[1] + radius2exc * sin(angle));

            line.setStrokeColor(colors[0]);
            line.setPosition(n);
            line.setStrokeWidth(strokeWidth);
            line.setName("Frame " + n + " connector");
            rm.addRoi(line);
        }

    }


    public static void main(String[] args) throws Exception {
        // set the plugins.dir property to make the plugin appear in the Plugins menu
        // see: https://stackoverflow.com/a/7060464/1207769
        Class<?> clazz = TestTrackerCrop_.class;
        java.net.URL url = clazz.getProtectionDomain().getCodeSource().getLocation();
        java.io.File file = new java.io.File(url.toURI());
        System.setProperty("plugins.dir", file.getAbsolutePath());

        // start ImageJ
        new ImageJ();

        // run the plugin
        IJ.runPlugIn(clazz.getName(), "");
    }


    @Override
    public boolean dialogItemChanged(GenericDialog genericDialog, AWTEvent awtEvent) {
        return true;
    }
}
