package gui;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import utils.PixelPathUtils;

import java.awt.*;
import java.io.IOException;
import java.util.ArrayList;

import static utils.ImageUtils.getMaxOfImage;
import static utils.ImageUtils.getSkeletonImage;
import static utils.PixelPathUtils.*;

public class TestSkeleton_ extends BaseGUI_ {

    private double blurSize = 2;
    private boolean invertLUT = true;

    @Override
    public boolean beforeSetupDialog(String arg) {
        autoOpenImp = true;
        return true;
    }

    @Override
    public void setupDialog() {
    }

    @Override
    public boolean loadSettings() {
        return true;
    }

    @Override
    public void execute() throws InterruptedException, IOException {

        ImageProcessor ip = imp.getProcessor();
        String title = imp.getTitle();

        ImageProcessor ipSkeleton = getSkeletonImage(ip.duplicate(), blurSize, invertLUT);
        short[] pixelsSegmented = connectedComponents(ipSkeleton);
        ShortProcessor spSegmented = new ShortProcessor(imp.getWidth(), imp.getHeight());
        spSegmented.setPixels(pixelsSegmented);

        new ImagePlus(title+": skeletons", ipSkeleton).show();
        new ImagePlus(title+": segmented", spSegmented).show();

        int nSegments = (int) getMaxOfImage(spSegmented);

        IJ.log("n skeletons is "+nSegments);

        int largestArea = 0;
        double divisionAngle = 0;
        int longestSegment = 0;

        for(int s=0; s<nSegments; s++){
            PixelPathUtils.regionProps thisProps = getProps(s+1, spSegmented);
            IJ.log("Skeleton "+(s+1)+" has area "+thisProps.area+", "+thisProps.endpointsList.size()+" endpoints and "+thisProps.branchpointsList.size()+" branchpoints");
            if(thisProps.endpointsList.size()>2) {
                ArrayList<int[]> prunedEndpoints = pruneEndpoints(thisProps);
                IJ.log("Skeleton " + (s + 1) + " after pruning: endpoints are (" + prunedEndpoints.get(0)[0] + ", " + prunedEndpoints.get(0)[1] + ") and (" + prunedEndpoints.get(1)[0] + ", " + prunedEndpoints.get(1)[1] + ")");
            }
            // use bounding box of longest skeleton to get angle
            ArrayList<int[]> endPoints = pruneEndpoints(thisProps);
            if(thisProps.area>largestArea){
                divisionAngle = getAngleBetweenPoints(endPoints.get(0), endPoints.get(1));
                largestArea = thisProps.area;
                longestSegment = s+1;
            }
        }

        IJ.log("Longest segment is "+longestSegment+" with division angle of "+Math.toDegrees(divisionAngle));


    }

    public static void main(String[] args) throws Exception {
        // set the plugins.dir property to make the plugin appear in the Plugins menu
        // see: https://stackoverflow.com/a/7060464/1207769
        Class<?> clazz = TestSkeleton_.class;
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
