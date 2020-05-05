package gui;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Roi;
import ij.io.DirectoryChooser;
import ij.measure.ResultsTable;
import ij.plugin.ChannelSplitter;
import ij.plugin.frame.RoiManager;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import utils.PixelPathUtils;
import utils.TrackerForSingleCrop;

import java.awt.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import static utils.HoughUtils.getHoughAccumulators;
import static utils.HoughUtils.peaksLocalMax;
import static utils.ImageUtils.*;
import static utils.PixelPathUtils.*;
import static utils.RoiUtils.getRandomColorPair;
import static utils.IOUtils.saveRoisAsZip;

public class BridgeTracker_ extends BaseGUI_ {

    /*
    Assumptions in data:
    - single channel
    - multiple timepoints

    - expect roimanager with rectangular rois
     */

    int minRadius, maxRadius, strokeWidth, nFrames, maxFrameGap;
    double minLength, maxGap, linkingDistance, thresholdModifier, exclusionFraction, blurSize;
    boolean showHough, showSegmentation, invertLUT, saveRoisToFile;

    String savePath;
    String imageTitle;

    DirectoryChooser DC;
    ChannelSplitter CS = new ChannelSplitter();
    RoiManager rm;

    @Override
    public boolean beforeSetupDialog(String arg) {
        autoOpenImp = true;
        return true;
    }

    @Override
    public void setupDialog() {
        gd = new NonBlockingGenericDialog("Bridge tracking");
        gd.addMessage("Please have a Roi manager open where the first Roi is a background region,\n" +
                "and other Rois fully contain division events");
        gd.addMessage("~~~ Circle parameters ~~~");
        gd.addNumericField("Minimum radius (pixels)", getPrefs("minRadius", 5), 0);
        gd.addNumericField("Maximum radius (pixels)", getPrefs("maxRadius", 20), 0);
        gd.addNumericField("Threshold modifier (leave at 0.75 unless low SNR)", getPrefs("thresholdModifier", 0.75), 2);
        gd.addMessage("~~~ Bridge parameters ~~~");
        gd.addNumericField("Blur size", getPrefs("blurSize", 2), 0);
        gd.addNumericField("Minimum line length (pixels)", getPrefs("minLength", 10), 0);
        gd.addNumericField("Maximum gap size (pixels)", getPrefs("maxGap", 10), 0);
        gd.addMessage("~~~ Miscellaneous ~~~");
        gd.addNumericField("Line width for plotting (pixels)", getPrefs("strokeWidth", 3), 0);
        gd.addNumericField("Maximum bridge to centre linking distance (pixels)", getPrefs("linkingDistance", 10), 0);
        gd.addNumericField("Maximum empty frame gap", getPrefs("maxFrameGap", 3), 0);
        gd.addNumericField("Bridge exclusion fraction relative to radius", getPrefs("exclusionFraction", 1.0), 2);
        gd.addCheckbox("Show Hough transforms", getPrefs("showHough", false));
        gd.addCheckbox("Show segmentation", getPrefs("showSegmentation", false));
        gd.addCheckbox("Invert LUT for skeletonization", getPrefs("invertLUT", true));
        gd.addCheckbox("Save ROIs to file", getPrefs("saveRoisToFile", true));
    }

    @Override
    public boolean loadSettings() {
        minRadius = (int) gd.getNextNumber();
        maxRadius = (int) gd.getNextNumber();
        thresholdModifier = gd.getNextNumber();

        blurSize = gd.getNextNumber();
        minLength = gd.getNextNumber();
        maxGap = gd.getNextNumber();

        strokeWidth = (int) gd.getNextNumber();
        linkingDistance = gd.getNextNumber();
        maxFrameGap = (int) gd.getNextNumber();
        exclusionFraction = gd.getNextNumber();

        showHough = gd.getNextBoolean();
        showSegmentation = gd.getNextBoolean();
        invertLUT = gd.getNextBoolean();
        saveRoisToFile = gd.getNextBoolean();


        setPrefs("minRadius", minRadius);
        setPrefs("maxRadius", maxRadius);
        setPrefs("thresholdModifier", thresholdModifier);

        setPrefs("blurSize", blurSize);
        setPrefs("minLength", minLength);
        setPrefs("maxGap", maxGap);

        setPrefs("strokeWidth", strokeWidth);
        setPrefs("linkingDistance", linkingDistance);
        setPrefs("maxFrameGap", maxFrameGap);

        setPrefs("exclusionFraction", exclusionFraction);

        setPrefs("showHough", showHough);
        setPrefs("showSegmentation", showSegmentation);
        setPrefs("showSegmentation", showSegmentation);
        setPrefs("invertLUT", invertLUT);
        setPrefs("saveRoisToFile", saveRoisToFile);

        return true;
    }

    @Override
    public void execute() throws InterruptedException, IOException {

        if(saveRoisToFile){
            DC = new DirectoryChooser("Select directory to save ROIs. zip file");
            savePath = DC.getDirectory();
            imageTitle = imp.getTitle();
        }

        // image preparation
        nFrames = imp.getNFrames();
        nChannels = imp.getNChannels();

        ImageStack imsChannel1 = new ImageStack();

        if(nChannels==2) imsChannel1 = CS.getChannel(imp, 1);
        else imsChannel1 = imp.getImageStack();

        RoiManager thisManager = new RoiManager().getInstance();

        if(thisManager==null){
            throw new IOException("Expecting ROI manager containing ROIs!");
        }

        rm = thisManager;

        int nRois = rm.getCount();

        // measure background
        rm.rename(0, "background");

        for(int r=1; r<nRois; r++) {

            IJ.log("ROI " + r);
            Roi roi = rm.getRoi(r);
            rm.rename(r, "ROI " + r);
            Rectangle roiBounds = roi.getBounds();

            int xTopLeft = roiBounds.x;
            int yTopLeft = roiBounds.y;
            int w = roiBounds.width;
            int h = roiBounds.height;

            // extract roi
            ImageStack imsRoiCh1 = imsChannel1.crop(xTopLeft, yTopLeft, 0, w, h, nFrames);

            new ImagePlus("ROI " + r, imsRoiCh1).show();

            TrackerForSingleCrop tracker = new TrackerForSingleCrop();
            tracker.setImageStack(imsRoiCh1);
            tracker.setTrackerOptions(blurSize, invertLUT, maxFrameGap, exclusionFraction);
            tracker.setHoughOptions(minRadius, maxRadius, 1.0, thresholdModifier);
            tracker.setDebugOptions(showHough, showSegmentation);
            tracker.setRoiOptions(rm, strokeWidth, r, xTopLeft, yTopLeft);

            tracker.doTracker();

            if(saveRoisToFile){
                saveRoisAsZip(rm, savePath, imageTitle);
            }
        }

    }

    public static void main(String[] args) throws Exception {
        // set the plugins.dir property to make the plugin appear in the Plugins menu
        // see: https://stackoverflow.com/a/7060464/1207769
        Class<?> clazz = BridgeTracker_.class;
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
