package gui;

import ij.IJ;
import ij.ImageJ;
import ij.ImageStack;
import ij.gui.*;
import ij.plugin.ChannelSplitter;
import ij.plugin.frame.RoiManager;
import utils.TrackerForSingleCrop;

import java.awt.*;
import java.io.IOException;
import static java.lang.Math.*;

public class TestTrackerCrop_ extends BaseGUI_ {


    int minRadius, maxRadius, strokeWidth, maxFrameGap;
    double thresholdModifier, exclusionFraction, blurSize;
    boolean showHough, showSegmentation, invertLUT, doFit;

    int nFrames;

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
        maxFrameGap = (int) gd.getNextNumber();
        exclusionFraction = gd.getNextNumber();

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

        ImageStack imsChannel1 = new ImageStack();

        if(nChannels==2) imsChannel1 = CS.getChannel(imp, 1);
        else imsChannel1 = imp.getImageStack();

//        RoiManager thisManager = new RoiManager().getInstance();
//
//        if(thisManager==null){
//            throw new IOException("Expecting ROI manager containing ROIs!");
//        }
//
//        rm = thisManager;
//
//        int nRois = rm.getCount();
//
//        if(nRois==0){
//            throw new IOException("Expecting ROI manager containing ROIs!");
//        }
//
//        double increment = 1;
//
//        // measure background
//        Roi roiBG = rm.getRoi(0);
//        rm.rename(0, "background");
//        double[] backgroundCh1 = getBackground(imsChannel1, roiBG);
//        double[] backgroundCh2 = new double[nFrames];
//        if(nChannels==2) backgroundCh2 = getBackground(imsChannel2, roiBG);

        TrackerForSingleCrop tracker = new TrackerForSingleCrop();
        tracker.setImageStack(imsChannel1);
        tracker.setTrackerOptions(blurSize, invertLUT, maxFrameGap, exclusionFraction);
        tracker.setHoughOptions(minRadius, maxRadius, 1.0, thresholdModifier);
        tracker.setDebugOptions(showHough, showSegmentation);
        tracker.setRoiOptions(rm, strokeWidth);

        tracker.doTracker();

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
