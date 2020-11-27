package gui;

import ij.*;
import ij.gui.*;
import ij.measure.Calibration;
import ij.measure.CurveFitter;
import ij.measure.ResultsTable;
import ij.plugin.ChannelSplitter;
import ij.plugin.Concatenator;
import ij.plugin.HyperStackConverter;
import ij.plugin.ZProjector;
import ij.plugin.filter.MaximumFinder;
import ij.plugin.frame.RoiManager;
import ij.process.AutoThresholder;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;

import java.awt.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.ConcurrentModificationException;

import static ij.process.ImageStatistics.getStatistics;
import static java.lang.Math.*;
import static utils.HoughUtils.getRadius;

public class StandardCircleFinder_ extends BaseGUI_ {

    Calibration calibration;
    ChannelSplitter CS = new ChannelSplitter();
    ZProjector zProjector = new ZProjector();
    MaximumFinder MF = new MaximumFinder();
    CurveFitter CF;
    Concatenator CC = new Concatenator();
    HyperStackConverter HSC = new HyperStackConverter();

    double pixelSize;
    String pixelUnit;

    RoiManager rm;

    public String unitChoice;
    public double maxRadius, minRadius;
    public int maxRadiusPx, minRadiusPx;
    public boolean overrideAutoThreshold;
    public int nHoughs;
    public int channel;

    int w, h, nFrames;

    double autoThreshold;
    double tolModifier;
    double threshold;

    boolean showThresholded, showAccumulators, showMaxProjected;

    boolean doManualCuration = false;


    ImagePlus imp_;


    @Override
    public boolean beforeSetupDialog(String arg) {
        autoOpenImp = true;
        return true;
    }

    @Override
    public void setupDialog() {
        calibration = imp.getCalibration();
        pixelSize = calibration.pixelWidth;
        pixelUnit = calibration.getUnit();

        channel = imp.getChannel();
        ImageProcessor ip;
        if(imp.getNFrames()>1 && imp.getNSlices()>1){
            IJ.error("This version does not yet know how to handle data with multiple T and Z, sorry! \n" +
                    "Please split out a frame series or slice series and then try again");
        }
        if(imp.getNChannels()>1){
            IJ.log("Warning - this version of CHough only works on single channels. CHough will run on current active channel (channel "+channel+")");
            ip = CS.getChannel(imp, channel).getProcessor(1).convertToFloatProcessor();
        }
        else{
            ip = imp.getProcessor();
        }

        autoThreshold = getThreshold(ip);

        gd = new NonBlockingGenericDialog("Find circles in single channel dataset");
        gd.addMessage("From image calibration, pixel size is "+String.format("%.3f", pixelSize)+" "+pixelUnit);

        gd.addRadioButtonGroup("Units for circle radius", new String[]{"Pixels", pixelUnit}, 1, 2, getPrefs("unitChoice", "Pixels"));
        gd.addNumericField("Minimum circle radius", getPrefs("minRadius", 10), 2);
        gd.addNumericField("Maximum circle radius ", getPrefs("maxRadius", 15), 2);

        gd.addMessage("By default, an automatic intensity threshold will be calculated per image unless overridden below");
        gd.addNumericField("Constant intensity threshold", autoThreshold, 0);
        gd.addCheckbox("Override auto threshold", getPrefs("overrideAutoThreshold", false));

        gd.addMessage("If you aren't getting enough circles detected, you can try decreasing the value below (e.g. 0.75)");
        gd.addNumericField("Tolerance modifier for peak detections", getPrefs("tolModifier", 1), 2);

        gd.addMessage("Debugging/curiosity options: check which intermediate stages you want to display");
        gd.addCheckbox("Thresholded images", getPrefs("showThresholded", false));
        gd.addCheckbox("Accumulator images", getPrefs("showAccumulators", false));
        gd.addCheckbox("Max-projected accumulators", getPrefs("showMaxProjected", true));

        gd.addMessage("Warning! Upon execution, this plugin will remove any currently stored Rois in the Roi Manager");
    }

    @Override
    public boolean loadSettings() {
        unitChoice = gd.getNextRadioButton();
        setPrefs("unitChoice", unitChoice);
        minRadius = gd.getNextNumber();
        setPrefs("minRadius", minRadius);
        maxRadius = gd.getNextNumber();
        setPrefs("maxRadius", maxRadius);

        if(unitChoice=="Pixels"){
            maxRadiusPx = (int) ceil(maxRadius);
            minRadiusPx = (int) floor(minRadius);
        }
        else{
            maxRadiusPx = (int) ceil(maxRadius/pixelSize);
            minRadiusPx = (int) floor(minRadius/pixelSize);
        }

        IJ.log("max radius (px) = "+maxRadiusPx+", min radius (px) = "+minRadiusPx);

        threshold = gd.getNextNumber();
        setPrefs("threshold", threshold);
        overrideAutoThreshold = gd.getNextBoolean();
        setPrefs("overrideAutoThreshold", overrideAutoThreshold);
        tolModifier = gd.getNextNumber();
        setPrefs("tolModifier", tolModifier);

        showThresholded = gd.getNextBoolean();
        setPrefs("showThresholded", showThresholded);
        showAccumulators = gd.getNextBoolean();
        setPrefs("showAccumulators", showAccumulators);
        showMaxProjected = gd.getNextBoolean();
        setPrefs("showMaxProjected", showMaxProjected);

        prefs.savePreferences();

        return true;
    }

    @Override
    public void execute() throws InterruptedException, IOException {
        if(imp_!=null){
            imp_.close();
        }
        imp_=null;

        RoiManager thisManager = new RoiManager().getInstance();

        if(thisManager!=null){
            rm = thisManager;
            rm.close();
        }

        rm = new RoiManager();

        // pull image info and stacks

        w = imp.getWidth();
        h = imp.getHeight();

        ImageStack ims;
        if(imp.getNChannels()>1) ims = CS.getChannel(imp, channel);
        else ims = imp.getImageStack();

        nFrames = ims.getSize();

        // hough set-up
        int nRadii = maxRadiusPx-minRadiusPx + 1;
        double increment = 1;
        if(nRadii<20) increment = 0.5;

        nHoughs = (int) ((maxRadiusPx-minRadiusPx)/increment) + 1;

        int wAcc = w+2*maxRadiusPx;
        int hAcc = h+2*maxRadiusPx;

        Rectangle rectangle = new Rectangle(maxRadiusPx - 1, maxRadiusPx - 1, w, h);

        // intermediate stages containers
        ImageStack imsThreshold = new ImageStack(w, h);
        ImagePlus impThreshold = null;
        ArrayList<ImageStack> houghsList = new ArrayList<>();
        ImagePlus impAccumulators = null;
        ImageStack imsMaxAccumulators = new ImageStack(w, h);
        ImagePlus impMaxAccumulators = null;
        ArrayList<double[]> circleParametersPx = new ArrayList<>();


        // do transforms

        for(int n=1; n<=nFrames; n++){
            float[] pixels =  (float[]) ims.getProcessor(n).convertToFloatProcessor().getPixels();
            IJ.showStatus("Frame "+n+": Calculating transform...");

            double threshold_;
            if(!overrideAutoThreshold) threshold_ = getThreshold(ims.getProcessor(n));
            else threshold_ = threshold;

            float[] pixelsThreshold = new float[w*h];
            ImageStack imsAccumulators = new ImageStack(w, h, nHoughs);

            // main hough loop
            for(int r_=0; r_<nHoughs; r_++){
                double r = minRadiusPx+(increment*r_);
                float[] accumulator = new float[wAcc*hAcc];
                IJ.showProgress(r_, nHoughs);

                for(int y=0; y<h; y++){
                    for(int x=0; x<w; x++){
                        int p = x + y*w;

                        if(pixels[p]>threshold_){
                            pixelsThreshold[p] = pixels[p];

                            for(int theta = 0; theta<360; theta++){
                                int a = (int) (x - r * Math.cos(toRadians(theta)));
                                int b = (int) (y - r * Math.sin(toRadians(theta)));

                                int a_ = a + maxRadiusPx;
                                int b_ = b + maxRadiusPx;

                                // weighted by original pixel intensity
                                accumulator[a_ + b_ * wAcc] += pixels[p]/360;
                            }
                        }

                    }
                }

                // update houghs image
                FloatProcessor fpAccumulator = new FloatProcessor(wAcc, hAcc, accumulator);
                fpAccumulator.setRoi(rectangle);
                // remove hough border and put into stack
                imsAccumulators.setProcessor(fpAccumulator.crop(), r_ +1);
                if(unitChoice=="Pixels") imsAccumulators.setSliceLabel("r=" + r + " pixels", r_ + 1);
                else imsAccumulators.setSliceLabel("r="+(r*pixelSize)+" "+pixelUnit, r_ + 1);
            }

            // update threshold image
            imsThreshold.addSlice(new FloatProcessor(w, h, pixelsThreshold));
            if(showThresholded){
                if(n==1){
                    impThreshold = new ImagePlus("Thresholded images", imsThreshold);
                    impThreshold.show();
                }
                else{
                    impThreshold.setStack(imsThreshold);
                    impThreshold.setSlice(n);
                }
            }

            // store accumulators image
            houghsList.add(imsAccumulators);

            if(showAccumulators){
                if(n==1){
                    impAccumulators = new ImagePlus("Accumulators per frame", imsAccumulators);
                    impAccumulators.show();
                }
                else{
                    impAccumulators = CC.concatenate(impAccumulators, new ImagePlus("", imsAccumulators.duplicate()),false);
                    impAccumulators = HSC.toHyperStack(impAccumulators, 1, nHoughs, n);
                    impAccumulators.setTitle("Accumulators per frame");
                    impAccumulators.show();
                }
            }

            // project accumulators image
            ImageProcessor ipMaxAccumulators = zProjector.run(new ImagePlus("", imsAccumulators), "max").getProcessor();
            imsMaxAccumulators.addSlice(ipMaxAccumulators);

            // update maxprojected image
            if(showMaxProjected){
                if(n==1){
                    impMaxAccumulators = new ImagePlus("Max-projected accumulators", imsMaxAccumulators);
                    impMaxAccumulators.show();
                }
                else{
                    impMaxAccumulators.setStack(imsMaxAccumulators);
                    impMaxAccumulators.setSlice(n);
                }
            }

            // find peaks in max projection
            Polygon maxima = MF.getMaxima(ipMaxAccumulators, getThreshold(ipMaxAccumulators)*tolModifier, true);
            int[] peaksX = maxima.xpoints;
            int[] peaksY = maxima.ypoints;
            int nPoints = maxima.npoints;
            IJ.log("Frame "+n+": "+nPoints+" candidate circles found");
//            PointRoi maxPoints = new PointRoi(peaksX, peaksY, nPoints);
//            maxPoints.setPosition(n);
//            maxPoints.setName("Circle centres frame "+n);
//            rm.addRoi(maxPoints);

            // find radii
            double[] radiusArray = new double[nHoughs];
            for(int h=0; h<nHoughs; h++){
                double radiusPx = minRadiusPx + increment*h;
                radiusArray[h] = radiusPx;
            }

            for(int i=0; i<nPoints; i++){
                int xC = peaksX[i];
                int yC = peaksY[i];

                // grab accumulator values for peak position - average in a 3x3 neighbourhood
                double[] houghVals = new double[nHoughs];

                for(int j=0; j<nHoughs; j++){
                    double val = 0;
                    double counter = 0;
                    FloatProcessor fp = imsAccumulators.getProcessor(j+1).convertToFloatProcessor();

                    for(int y = max(0, yC-1); y<=min(yC+1, h-1); y++){
                        for (int x = max(0, xC - 1); x <= min(xC + 1, w - 1); x++) {
                            val += fp.getf(x, y);
                            counter++;
                        }
                    }
                    val /= counter;
                    houghVals[j] = val;
                }

                // find local max of radius vs accumulator value
                double radiusMax = getRadius(radiusArray, houghVals, 2);

                // find max by Gaussian fit
                CF = new CurveFitter(radiusArray, houghVals);
                CF.doFit(CurveFitter.GAUSSIAN);
                double[] params = CF.getParams();
                double radiusGauss = params[2];
                double errorGauss = params[3];

                double fractionalDiffMaxGauss = abs(radiusGauss-radiusMax)/radiusGauss;

                // purge bad results
                if(radiusGauss<minRadiusPx || radiusGauss>maxRadiusPx) continue;
                if(errorGauss>0.5*(maxRadiusPx-minRadiusPx)) continue;
                if(fractionalDiffMaxGauss>0.2) continue;

                double[] circleInfoPx = new double[]{n, xC, yC, radiusGauss, errorGauss, radiusMax};
                circleParametersPx.add(circleInfoPx);
            }

        }

        // make circles with option of manual curation
        int nCirclesFound = circleParametersPx.size();

        NonBlockingGenericDialog gdConfirmCuration = new NonBlockingGenericDialog("Manual curation");
        gdConfirmCuration.addMessage("There were "+nCirclesFound+" circles found in your image. \n" +
                "Would you like to manually curate each of these circles?");
        gdConfirmCuration.addRadioButtonGroup("Do manual curation?", new String[]{"Yes", "No"},
                1, 2, "No");
        gdConfirmCuration.showDialog();

        String curation = gdConfirmCuration.getNextRadioButton();
        if(curation=="Yes") doManualCuration = true;

        ResultsTable rt = new ResultsTable();

        int rejectionCount = 0;

        for(int c=0; c<nCirclesFound; c++){
            double[] circleParamsPx = circleParametersPx.get(c);

            int cT = (int) circleParamsPx[0];
            float cX = (float) circleParamsPx[1];
            float cY = (float) circleParamsPx[2];
            float cR = (float) circleParamsPx[3];
            double cE = circleParamsPx[4];
            double cR_ = circleParamsPx[5];

            //EllipseRoi ellipseRoi = new EllipseRoi(cX - cR, cY - cR, cX + cR, cY + cR, 1);
            OvalRoi ellipseRoi = new OvalRoi(cX-cR, cY-cR, cR*2, cR*2);
            ellipseRoi.setStrokeColor(Color.cyan);
            ellipseRoi.setStrokeWidth(cE);
            ellipseRoi.setPosition(cT);
            imp.setSlice(cT);
            imp.setRoi(ellipseRoi);

            boolean reject = false;

            if(doManualCuration) {
                NonBlockingGenericDialog gdCheck = new NonBlockingGenericDialog("Check circle " + (c + 1) + "/" + nCirclesFound);
                gdCheck.addMessage("Circle at (" + cX + "," + cY + "):\n Radius (px) = " + cR);
                gdCheck.addCheckbox("Reject nucleus? ", false);
                gdCheck.showDialog();

                reject = gdCheck.getNextBoolean();
            }

            if(!reject) {

                rt.incrementCounter();
                rt.addValue("Stack position", cT);
                rt.addValue("x-position (px)", cX);
                rt.addValue("y-position (px)", cY);
                if(unitChoice!="Pixels") {
                    rt.addValue("x-position ("+pixelUnit+")", cX*pixelSize);
                    rt.addValue("y-position ("+pixelUnit+")", cY*pixelSize);
                }
                rt.addValue("radius from Gaussian fit (px)", cR);
                rt.addValue("uncertainty of Gaussian fit (σ, px)", cE);
                rt.addValue("radius from local maximum (px)", cR_);
                if(unitChoice!="Pixels") {
                    rt.addValue("radius from Gaussian fit ("+pixelUnit+")", cR*pixelSize);
                    rt.addValue("uncertainty of Gaussian fit (σ, "+pixelUnit+")", cE*pixelSize);
                    rt.addValue("radius from local maximum  ("+pixelUnit+")", cR_*pixelSize);
                }

                rm.addRoi(ellipseRoi);
            }
            else{
                rejectionCount++;
                ellipseRoi.setStrokeColor(Color.RED);
                rm.addRoi(ellipseRoi);
            }

        }
        rt.show("Hough results ("+rejectionCount+" nuclei manually rejected (Rois colored red))");
        IJ.log("Analysis completed!");
        if(nFrames>1){
            IJ.log("Tip: for multi-frame data, go to More >> Options... >> 'Associate Show All ROIs with slices' in the RoiManager");
        }

    }

    @Override
    public boolean dialogItemChanged(GenericDialog genericDialog, AWTEvent awtEvent) {
        return true;
    }


    public static void main(String[] args) throws Exception {
        // set the plugins.dir property to make the plugin appear in the Plugins menu
        // see: https://stackoverflow.com/a/7060464/1207769
        Class<?> clazz = StandardCircleFinder_.class;
        java.net.URL url = clazz.getProtectionDomain().getCodeSource().getLocation();
        java.io.File file = new java.io.File(url.toURI());
        System.setProperty("plugins.dir", file.getAbsolutePath());

        // start ImageJ
        new ImageJ();

        // run the plugin
        IJ.runPlugIn(clazz.getName(), "");
    }

    public static int getThreshold(ImageProcessor ip) {

        AutoThresholder autoThresholder = new AutoThresholder();
        ImageStatistics stats = getStatistics(ip);
        double min = stats.min;
        double max = stats.max;
        double interval = (max - min) / 256;

        int threshold = autoThresholder.getThreshold(AutoThresholder.Method.Otsu, ip.getHistogram(256));
        return (int) (threshold * interval + min);
    }
}
