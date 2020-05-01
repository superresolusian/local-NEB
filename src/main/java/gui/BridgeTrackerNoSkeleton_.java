package gui;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.*;
import ij.measure.ResultsTable;
import ij.plugin.ChannelSplitter;
import ij.plugin.frame.RoiManager;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

import java.awt.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Random;

import static java.lang.Math.*;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static utils.HoughUtils.getHoughAccumulators;
import static utils.HoughUtils.peaksLocalMax;
import static utils.ImageUtils.*;
import static utils.PixelPathUtils.*;
import static utils.RoiUtils.getRandomColorPair;

public class BridgeTrackerNoSkeleton_ extends BaseGUI_ {

        /*
    Assumptions in data:
    - single channel
    - multiple timepoints

    - expect roimanager with rectangular rois
    - first roi should be for background estimation
     */

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

        IJ.log("Number of ROIs to examine is: " + nRois);

        double increment = 1;

        // determine background

        Roi roiBG = rm.getRoi(0);
        rm.rename(0, "background");
        Rectangle roiBoundsBG = roiBG.getBounds();

        int xTopLeftBG = roiBoundsBG.x;
        int yTopLeftBG = roiBoundsBG.y;
        int wBG = roiBoundsBG.width;
        int hBG = roiBoundsBG.height;
        int nPixelsBG = wBG*hBG;

        ImageStack imsBackgroundCh1 = imsChannel1.crop(xTopLeftBG, yTopLeftBG, 0, wBG, hBG, nFrames);
        ImageStack imsBackgroundCh2 = new ImageStack();
        if (nChannels == 2) imsBackgroundCh2 = imsChannel2.crop(xTopLeftBG, yTopLeftBG, 0, wBG, hBG, nFrames);

        double[] backgroundCh1 = new double[nFrames];
        double[] backgroundCh2 = new double[nFrames];

        for(int n=1; n<=nFrames; n++){
            ImageProcessor ip1 = imsBackgroundCh1.getProcessor(n);
            ImageProcessor ip2 = new FloatProcessor(wBG, hBG);
            if(nChannels==2) ip2 = imsBackgroundCh2.getProcessor(n);

            double avBG1 = 0, avBG2 = 0;
            for(int y=0; y<hBG; y++){
                for(int x=0; x<wBG; x++){
                    avBG1 += ip1.getf(x, y)/nPixelsBG;
                    if(nChannels==2) avBG2 += ip2.getf(x, y)/nPixelsBG;
                }
            }
            backgroundCh1[n-1] = avBG1;
            if(nChannels==2) backgroundCh2[n-1] = avBG2;
        }

        ResultsTable rt = new ResultsTable();

        // loop through other rois

        for(int r=1; r<nRois; r++) {
            Color[] colors = getRandomColorPair(true);

            IJ.log("ROI " + r);
            Roi roi = rm.getRoi(r);
            rm.rename(r, "ROI " + r);
            Rectangle roiBounds = roi.getBounds();

            int xTopLeft = roiBounds.x;
            int yTopLeft = roiBounds.y;
            int w = roiBounds.width;
            int h = roiBounds.height;

            // extract roi from each channel
            ImageStack imsRoiCh1 = imsChannel1.crop(xTopLeft, yTopLeft, 0, w, h, nFrames);
            ImageStack imsRoiCh2 = new ImageStack();
            if (nChannels == 2) imsRoiCh2 = imsChannel2.crop(xTopLeft, yTopLeft, 0, w, h, nFrames);

            new ImagePlus("ROI " + r, imsRoiCh1).show();

            //// determine division angle
            ImageStack imsNormalised = normaliseToMax(imsRoiCh1);
            FloatProcessor fpRoiAverageIntensity = averageProject(imsNormalised);
            ImageProcessor ipSkeletonAverage = getSkeletonImage(fpRoiAverageIntensity, blurSize, invertLUT);

            short[] pixelsSegmentedAverage = connectedComponents(ipSkeletonAverage);
            ShortProcessor spSegmentedAverage = new ShortProcessor(w, h);
            spSegmentedAverage.setPixels(pixelsSegmentedAverage);

            int nSegmentsAverage = (int) getMaxOfImage(spSegmentedAverage);

            // get largest segment
            int bestArea = 0;
            int biggestSegment = 0;

            ArrayList<regionProps> propsList = new ArrayList<>();

            for(int s=0; s<nSegmentsAverage; s++) {
                // get properties for each segment
                regionProps thisProps = getProps(s + 1, spSegmentedAverage);

                // check size of segment
                int area = thisProps.area;
                if(area>bestArea){
                    bestArea = area;
                    biggestSegment = s;
                }

                propsList.add(thisProps);
            }

            regionProps bestProps = propsList.get(biggestSegment);
            ArrayList<int[]> endpoints = bestProps.endpointsList;
            int nEndpoints = endpoints.size();

            double optimumAngle = 0;

            if(nEndpoints==2){
                int[] ep1 = endpoints.get(0);
                int[] ep2 = endpoints.get(1);
                optimumAngle = getAngleBetweenPoints(ep1, ep2);
            }
            else{
                // get list of pixels belonging to skeleton
                // TODO: replace with pruning?
                Rectangle boundingBox = bestProps.boundingBox;
                int y0 = boundingBox.y;
                int x0 = boundingBox.x;
                int x1 = x0+boundingBox.width;
                int y1 = y0+boundingBox.height;

                ArrayList<Integer> nodeList = new ArrayList<Integer>();
                for (int y = y0; y <= y1; y++) {
                    for (int x = x0; x <= x1; x++) {
                        int p = y * w + x;
                        if (pixelsSegmentedAverage[p] == biggestSegment + 1) {
                            nodeList.add(p);
                        }
                    }
                }
                LinkedHashMap<Integer, int[]> nodeNeighbours = findNodeNeighbours(spSegmentedAverage, nodeList, biggestSegment, bestArea);

                int[] beste1 = new int[2];
                int[] beste2 = new int[2];
                int[] beste1_ = new int[2];
                int[] beste2_ = new int[2];
                int longestPath = 0;
                int nextLongestPath = 0;

                for (int e1 = 0; e1 < nEndpoints - 1; e1++) {

                    int[] start = endpoints.get(e1); // get pixel index of endpoint e1
                    int p1 = start[0] + start[1] * w;

                    for (int e2 = e1 + 1; e2 < nEndpoints; e2++) {

                        int[] end = endpoints.get(e2); // get pixel index of endpoint e2
                        int p2 = end[0] + end[1] * w;

                        int[] path = findPath(p1, p2, nodeNeighbours);
                        if (path.length > longestPath) {
                            longestPath = path.length;
                            beste1 = start;
                            beste2 = end;
                        }
                        else if(path.length > nextLongestPath){
                            nextLongestPath = path.length;
                            beste1_ = start;
                            beste2_ = end;
                        }
                    }
                }

                double angle1 = getAngleBetweenPoints(beste1, beste2);
                double angle2 = getAngleBetweenPoints(beste1_, beste2_);

                optimumAngle = toRadians(0.5*(angle1+angle2));
            }

            boolean doContinueCounter = false;
            int continueCounter = 0;
            int[] bestCentre1 = new int[0], bestCentre2 = new int[0];
            int[] previousCentre1 = new int[0], previousCentre2 = new int[0];
            int bestRadius1 = 0, bestRadius2 = 0;
            int previousRadius1 = 0, previousRadius2 = 0;

            ImageStack imsProjectedHoughs = new ImageStack(w, h);
            ImagePlus impProjectedHoughs = new ImagePlus();

            for (int n = 1; n <= nFrames; n++) {
                if (continueCounter > maxFrameGap){
                    IJ.log("max frame gap reached");
                    break;
                }

                //// get circles in each frame
                ImageProcessor ipRoi1 = imsRoiCh1.getProcessor(n).duplicate();
                ImageProcessor ipRoi2 = new FloatProcessor(w, h);
                if (nChannels == 2) ipRoi2 = imsRoiCh2.getProcessor(n).duplicate();
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
                        impProjectedHoughs = new ImagePlus("Hough transforms for ROI " + r, imsProjectedHoughs);
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
                LinkedHashMap<int[], Integer> radii = new LinkedHashMap<int[], Integer>();

                for (int i = 0; i < nMaxima; i++) {
                    int thisX = maxima.xpoints[i];
                    int thisY = maxima.ypoints[i];
                    maximaCoordinates.add(new int[]{thisX, thisY});

                    double bestRadiusVal = 0;
                    int bestRadius = 0;
                    for (int j = 0; j < nHoughs; j++) {
                        ImageProcessor ip = imsAccumulators.getProcessor(j + 1);

                        double radVal = ip.getf(thisX, thisY);

                        if (radVal > bestRadiusVal) {
                            bestRadiusVal = radVal;
                            bestRadius = (int) (minRadius + j * increment);
                        }
                    }
                    radii.put(maximaCoordinates.get(i), bestRadius);
                }

                // see if good pair exists

                int[] bestPair = new int[2];
                boolean foundGoodPair = false;
                double bestAngleDifference = Double.MAX_VALUE;

                for (int i = 0; i < nMaxima - 1; i++) {
                    int[] testCentre1 = maximaCoordinates.get(i);
                    int r1 = radii.get(testCentre1);

                    for (int j = i + 1; j < nMaxima; j++) {
                        int[] testCentre2 = maximaCoordinates.get(j);
                        int r2 = radii.get(testCentre2);

                        if (getDistance(testCentre1, testCentre2) < r1 * 3) {
                            IJ.log("centres are too close!");
                            continue; // too close
                        }

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

                    boolean condition1 = dist11 < previousRadius1 && dist22 < previousRadius2;
                    boolean condition2 = dist12 < previousRadius1 && dist21 < previousRadius2;

                    if (!condition1 && !condition2) {
                        IJ.log("conditions not met");
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

                EllipseRoi ellipseRoi1 = new EllipseRoi(xTopLeft + bestCentre1[0] - bestRadius1, yTopLeft + bestCentre1[1] - bestRadius1, xTopLeft + bestCentre1[0] + bestRadius1, yTopLeft + bestCentre1[1] + bestRadius1, 1);
                ellipseRoi1.setPosition(n);
                ellipseRoi1.setStrokeColor(colors[0]);
                ellipseRoi1.setName("ROI " + r + " frame " + n + " centre 1");
                rm.addRoi(ellipseRoi1);

                EllipseRoi ellipseRoi2 = new EllipseRoi(xTopLeft + bestCentre2[0] - bestRadius2, yTopLeft + bestCentre2[1] - bestRadius2, xTopLeft + bestCentre2[0] + bestRadius2, yTopLeft + bestCentre2[1] + bestRadius2, 1);
                ellipseRoi2.setPosition(n);
                ellipseRoi2.setStrokeColor(colors[1]);
                ellipseRoi2.setName("ROI " + r + " frame " + n + " centre 2");
                rm.addRoi(ellipseRoi2);

                double radius1exc = bestRadius1 * exclusionFraction;
                double radius2exc = bestRadius2 * exclusionFraction;

                double averageBridgeIntensityCh1 = 0, averageBridgeIntensityCh2 = 0;
                double bridgeLength = 0;

                double angle = getAngleBetweenPoints(bestCentre1, bestCentre2);
                Line line = new Line(bestCentre1[0] + radius1exc * cos(angle), bestCentre1[1] - radius1exc * sin(angle),
                        bestCentre2[0] - radius2exc * cos(angle), bestCentre2[1] + radius2exc * sin(angle));

                line.setStrokeWidth(strokeWidth);
                Polygon poly = line.getPolygon();
                Rectangle polyBounds = poly.getBounds();

                int nPoints = 0;

                for (int y = polyBounds.y; y <= polyBounds.y + polyBounds.height; y++) {
                    for (int x = polyBounds.x; x <= polyBounds.x + polyBounds.width; x++) {
                        if (poly.contains(x, y)) {
                            averageBridgeIntensityCh1 += ipRoi1.getf(x, y);
                            if (nChannels == 2) averageBridgeIntensityCh2 += ipRoi2.getf(x, y);
                            nPoints++;
                        }
                    }
                }

                averageBridgeIntensityCh1 /= nPoints;
                averageBridgeIntensityCh2 /= nPoints;

                line = new Line(xTopLeft + bestCentre1[0] + radius1exc * cos(angle), yTopLeft + bestCentre1[1] - radius1exc * sin(angle),
                        xTopLeft + bestCentre2[0] - radius2exc * cos(angle), yTopLeft + bestCentre2[1] + radius2exc * sin(angle));

                line.setStrokeColor(colors[0]);
                line.setPosition(n);
                line.setStrokeWidth(strokeWidth);
                line.setName("ROI " + (r) + " frame " + n + " connector");
                rm.addRoi(line);

                bridgeLength = line.getLength();

                // populate results table
                double area1 = Math.PI * bestRadius1 * bestRadius1;
                double meanIntensity1Ch1 = 0, meanIntensity1Ch2 = 0;
                for (int y = (int) max(0, bestCentre1[1] - bestRadius1); y <= min(h - 1, bestCentre1[1] + bestRadius1); y++) {
                    for (int x = (int) max(0, bestCentre1[0] - bestRadius1); x <= min(w - 1, bestCentre1[0] + bestRadius1); x++) {
                        double dist = (x - bestCentre1[0]) * (x - bestCentre1[0]) + (y - bestCentre1[1]) * (y - bestCentre1[1]);
                        if (dist > (bestRadius1 * bestRadius1)) continue;
                        meanIntensity1Ch1 += ipRoi1.getf(x, y) / area1;
                        if (nChannels == 2) meanIntensity1Ch2 += ipRoi1.getf(x, y) / area1;
                    }
                }
                double area2 = Math.PI * bestRadius2 * bestRadius2;
                double meanIntensity2Ch1 = 0, meanIntensity2Ch2 = 0;
                for (int y = (int) max(0, bestCentre2[1] - bestRadius2); y <= min(h - 1, bestCentre2[1] + bestRadius2); y++) {
                    for (int x = (int) max(0, bestCentre2[0] - bestRadius2); x <= min(w - 1, bestCentre2[0] + bestRadius2); x++) {
                        double dist = (x - bestCentre2[0]) * (x - bestCentre2[0]) + (y - bestCentre2[1]) * (y - bestCentre2[1]);
                        if (dist > (bestRadius2 * bestRadius2)) continue;
                        meanIntensity2Ch1 += ipRoi2.getf(x, y) / area2;
                        if (nChannels == 2) meanIntensity2Ch2 += ipRoi2.getf(x, y) / area2;
                    }
                }

                double background1 = backgroundCh1[n-1];
                double background2 = 0;
                if(nChannels==2) background2 = backgroundCh2[n-1];

                double meanIntensityCh1 = (meanIntensity1Ch1 + meanIntensity2Ch1) / 2;
                double meanIntensityCh2 = (meanIntensity1Ch2 + meanIntensity2Ch2) / 2;
                double meanRadius = (bestRadius1 + bestRadius2) / 2;
                double distanceBetweenNuclei = sqrt((bestCentre1[0] - bestCentre2[0]) * (bestCentre1[0] - bestCentre2[0]) + (bestCentre1[1] - bestCentre2[1]) * (bestCentre1[1] - bestCentre2[1]));
                double separationIndex = distanceBetweenNuclei / meanRadius;

                double meanIntensityCh1BG = meanIntensityCh1-background1;
                double meanIntensityCh2BG = meanIntensityCh2-background2;
                double bridgeIntensityCh1BG = averageBridgeIntensityCh1-background1;
                double bridgeIntensityCh2BG = averageBridgeIntensityCh2-background2;


                rt.incrementCounter();
                rt.addValue("ROI", r);
                rt.addValue("Frame", n);
                rt.addValue("Radius 1", bestRadius1);
                rt.addValue("Radius 2", bestRadius2);
                rt.addValue("Distance between nucleus centres", distanceBetweenNuclei);
                rt.addValue("Mean intensity nuclei Ch1", meanIntensityCh1);
                rt.addValue("Mean intensity nuclei Ch1 - bg", meanIntensityCh1BG);
                if (nChannels == 2){
                    rt.addValue("Mean intensity nuclei Ch2", meanIntensityCh2);
                    rt.addValue("Mean intensity nuclei Ch2 - bg", meanIntensityCh2BG);
                }
                rt.addValue("Bridge length", bridgeLength);
                rt.addValue("Mean bridge intensity Ch1", averageBridgeIntensityCh1);
                rt.addValue("Mean bridge intensity Ch1 - BG", bridgeIntensityCh1BG);
                if (nChannels == 2){
                    rt.addValue("Mean bridge intensity Ch2", averageBridgeIntensityCh2);
                    rt.addValue("Mean bridge intensity Ch2 - BG", bridgeIntensityCh2BG);
                }

                rt.addValue("Separation index", separationIndex);
                rt.addValue("Normalised bridge intensity Ch1", bridgeIntensityCh1BG/meanIntensityCh1BG);
                if (nChannels == 2)
                    rt.addValue("Normalised bridge intensity Ch2", bridgeIntensityCh2BG/meanIntensityCh2BG);

            }

        }

        rt.show("Results");

    }

    public static void main(String[] args) throws Exception {
        // set the plugins.dir property to make the plugin appear in the Plugins menu
        // see: https://stackoverflow.com/a/7060464/1207769
        Class<?> clazz = BridgeTrackerNoSkeleton_.class;
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
