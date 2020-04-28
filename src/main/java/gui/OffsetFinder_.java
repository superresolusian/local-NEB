package gui;

import ij.*;
import ij.gui.EllipseRoi;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.PointRoi;
import ij.measure.ResultsTable;
import ij.plugin.ChannelSplitter;
import ij.plugin.HyperStackConverter;
import ij.plugin.filter.MaximumFinder;
import ij.plugin.frame.RoiManager;
import ij.process.FloatProcessor;

import java.awt.*;
import java.io.IOException;
import java.util.ArrayList;

import static java.lang.Math.*;
import static utils.ArrayUtils.intToFloat;
import static utils.HoughUtils.getOptimumRange;
import static utils.HoughUtils.getRadius;
import static utils.ImageUtils.getFpStats;
import static utils.ImageUtils.groupedZProject;
import static utils.ImageUtils.maxProject;

public class OffsetFinder_ extends BaseGUI_  {

    ChannelSplitter CS = new ChannelSplitter();
    MaximumFinder MF = new MaximumFinder();
    RoiManager rm;

    public int maxRadius, minRadius;
    public int nChannels, nHoughs;

    int tolerance1, tolerance2;
    int w, h, nFrames;
    int averaging = 1;

    double autoThreshold1, autoThreshold2, autoTolerance1, autoTolerance2;
    double pixelSize, maxDisplacement, maxDisplacementPx, threshold1, threshold2;
    double[] thresholds, tolerances;

    float[] pixelsZOffset;

    ArrayList<float[]> xPositions, yPositions;

    boolean doManualCuration, doZCorrection;
    boolean[] isPeakGood;

    String[] imageTitles;
    Color[] colors = new Color[]{Color.magenta, Color.cyan, Color.yellow};

    ImagePlus impOffset;
    ImagePlus imp_;

    @Override
    public boolean beforeSetupDialog(String s) {
        autoOpenImp = true;

        String[] tempImageTitles = WindowManager.getImageTitles();
        int nImages = WindowManager.getImageCount();
        imageTitles = new String[nImages+1];
        imageTitles[0] = "No z-correction";
        for(int n=1; n<=nImages; n++) {
            imageTitles[n] = tempImageTitles[n-1];
        }

        return true;
    }

    @Override
    public void setupDialog() {
        if(imp.getNChannels()!=2){
            IJ.error("Expecting image with 2 colour channels");
        }
        FloatProcessor fp1 = CS.getChannel(imp, 1).getProcessor(1).convertToFloatProcessor();
        float[] fp1Stats = getFpStats(fp1);
        autoThreshold1 = fp1Stats[0] + fp1Stats[1];

        FloatProcessor fp2 = CS.getChannel(imp, 2).getProcessor(1).convertToFloatProcessor();
        float[] fp2Stats = getFpStats(fp2);
        autoThreshold2 = fp2Stats[0] + fp2Stats[1];

        autoTolerance1 = autoThreshold1;
        autoTolerance2 = autoThreshold2/2;


        gd = new NonBlockingGenericDialog("Radius difference between channels");
        gd.addNumericField("Pixel size (nm)", getPrefs("pixelSize", 100), 1);
        gd.addNumericField("Maximum displacement (nm)", getPrefs("maxDisplacement", 250), 0);
        gd.addNumericField("Minimum circle radius (pixels)", getPrefs("minRadius", 10), 0);
        gd.addNumericField("Maximum circle radius (pixels)", getPrefs("maxRadius", 15), 0);

        gd.addNumericField("Intensity threshold in channel 1", autoThreshold1, 0);
        gd.addNumericField("Intensity threshold in channel 2", autoThreshold2, 0);
        gd.addMessage("Please make sure you have set the above values before proceeding :-)");
        gd.addNumericField("Tolerance for peak detections in channel 1", autoTolerance1,0);
        gd.addNumericField("Tolerance for peak detections in channel 2", autoTolerance2,0);
        gd.addChoice("Image for z correction", imageTitles, imageTitles[0]);
        gd.addCheckbox("Manually curate results?", getPrefs("doManualCuration", true));
    }

    @Override
    public boolean loadSettings() {
        pixelSize = gd.getNextNumber();
        setPrefs("pixelSize", pixelSize);
        maxDisplacement = gd.getNextNumber();
        setPrefs("maxDisplacement", maxDisplacement);
        minRadius = (int) gd.getNextNumber();
        setPrefs("minRadius", minRadius);
        maxRadius = (int) gd.getNextNumber();
        setPrefs("maxRadius", maxRadius);
        nHoughs = (maxRadius-minRadius)*2 + 1;
        threshold1 = gd.getNextNumber();
        setPrefs("threshold1", threshold1);
        threshold2 = gd.getNextNumber();
        setPrefs("threshold2", threshold2);
        tolerance1 = (int) gd.getNextNumber();
        setPrefs("tolerance1", tolerance1);
        tolerance2 = (int) gd.getNextNumber();
        setPrefs("tolerance2", tolerance2);

        String path = gd.getNextChoice();
        if(path.equals(imageTitles[0])) doZCorrection = false;
        else{
            doZCorrection = true;
            impOffset = WindowManager.getImage(path);
        }

        doManualCuration = gd.getNextBoolean();
        setPrefs("doManualCuration", doManualCuration);

        prefs.savePreferences();

        thresholds = new double[2];
        thresholds[0] = threshold1;
        thresholds[1] = threshold2;


        tolerances = new double[2];
        tolerances[0] = tolerance1;
        tolerances[1] = tolerance2;

        maxDisplacementPx = maxDisplacement/pixelSize;

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

        // grab correction map if required
        if(doZCorrection){
            assert(impOffset.getWidth()==w && impOffset.getHeight()==h);
            FloatProcessor fp = impOffset.getProcessor().convertToFloatProcessor();
            pixelsZOffset = (float[]) fp.getPixels();
        }

        // pull image info and stacks

        w = imp.getWidth();
        h = imp.getHeight();
        nFrames = imp.getNFrames();
        if(nFrames>1){
            IJ.log("This plugin currently only works for single frame images, sorry! Will proceed using frame 1.");
            nFrames = 1;
        }
        nChannels = imp.getNChannels();

        ArrayList<ImageStack> imageStacks = new ArrayList<ImageStack>(nChannels);
        ImageStack ims;

        for(int ch=0; ch<nChannels; ch++){
            ims = CS.getChannel(imp, ch+1);
            imageStacks.add(ims);
        }

        // hough set-up
        int wAcc = w+2*maxRadius;
        int hAcc = h+2*maxRadius;
        double increment = 0.5;

        ImageStack imsRadius = new ImageStack(w, h, nHoughs*nChannels);
        Rectangle rectangle = new Rectangle(maxRadius - 1, maxRadius - 1, w, h);

        // do transforms
        ImageStack imsThreshold = new ImageStack(w, h, nChannels);

        int stackCounter = 1;

        for(int ch=1; ch<=nChannels; ch++) {
            ims = imageStacks.get(ch-1);

            float[] pixels = (float[]) ims.getProcessor(1).convertToFloatProcessor().getPixels();

            IJ.showStatus("Channel "+ch+": Calculating transform...");

            int n = 1;
            float[] pixelsThreshold = new float[w*h];

            long start = System.currentTimeMillis();

            // main hough loop
            for (int r_ = 0; r_ < nHoughs; r_++) {
                double r = minRadius+(increment*r_);
                float[] accumulator = new float[wAcc * hAcc];
                IJ.showProgress(n, nHoughs);

                for (int y = 0; y < h; y++) {
                    for (int x = 0; x < w; x++) {
                        int p = x + y * w;

                        if (pixels[p] > thresholds[ch-1]) {
                            pixelsThreshold[p] = pixels[p];
                            for (int theta = 0; theta < 360; theta++) {
                                int a = (int) (x - r * Math.cos(toRadians(theta)));
                                int b = (int) (y - r * Math.sin(toRadians(theta)));

                                int a_ = a + maxRadius;
                                int b_ = b + maxRadius;
                                // weighted by original pixel intensity
                                accumulator[a_ + b_ * wAcc] += pixels[p]/360;
                            }
                        }

                    }
                }

                // update threshold image
                imsThreshold.setProcessor(new FloatProcessor(w, h, pixelsThreshold), ch);

                FloatProcessor fpHough = new FloatProcessor(wAcc, hAcc, accumulator);
                fpHough.setRoi(rectangle);
                // remove hough border
                imsRadius.setProcessor(fpHough.crop(), stackCounter);
                imsRadius.setSliceLabel("r=" + r, stackCounter);
                stackCounter++;
                n++;
            }
            IJ.showProgress(1);

            long stop = System.currentTimeMillis();
            IJ.log("Analysis took for channel "+ch+" took "+(stop-start)+"ms");

        }


        ImagePlus impHough = new ImagePlus("Hough transforms", imsRadius);
        ImagePlus impHoughHS = HyperStackConverter.toHyperStack(impHough, nChannels, nHoughs, nFrames, "xyztc", "color");
        impHoughHS.show();

        ImagePlus impThreshold = new ImagePlus("Thresholded pixels", imsThreshold);
        ImagePlus impThresholdHS = HyperStackConverter.toHyperStack(impThreshold, nChannels, 1, nFrames, "xyztc", "composite");
        impThresholdHS.show();

        // get projected stack of first channel
        ImageStack ims_ = new ImageStack(w, h, nChannels);
        for(int ch=1; ch<=nChannels; ch++){
            ImageStack imsTemp = CS.getChannel(impHoughHS, ch);
            ims_.setProcessor(groupedZProject(imsTemp, nHoughs).getProcessor(1),ch);
        }
        imp_ = new ImagePlus("Projected", ims_);
        imp_.show();

        ArrayList<ImageStack> houghStacks = new ArrayList<ImageStack>(nChannels);
        for(int i=0; i<nChannels; i++){
            houghStacks.add(CS.getChannel(impHoughHS, i+1));
        }

        xPositions = new ArrayList<float[]>(nChannels);
        yPositions = new ArrayList<float[]>(nChannels);

        // loops through channels and find maxima
        // store in arraylists xpositions and ypositions
        for(int ch=1; ch<=nChannels; ch++){

            FloatProcessor fp = ims_.getProcessor(ch).convertToFloatProcessor();

            Polygon poly = MF.getMaxima(fp, tolerances[ch-1], true);
            //IJ.log("prominence for find maxima is "+tolerance*(thresholds[ch-1]/thresholds[0]));
            PointRoi roi = new PointRoi(poly);
            roi.setStrokeColor(colors[ch-1]);
            roi.setPointType(2);
            roi.setPosition(ch);

            rm.addRoi(roi);

            float[] xPos = intToFloat(poly.xpoints);
            float[] yPos = intToFloat(poly.ypoints);

            xPositions.add(ch-1,xPos);
            yPositions.add(ch-1,yPos);

            IJ.log("Peaks found! There are "+xPos.length+" peaks in channel "+ch);

        }
        int nPeaksChannel1 = xPositions.get(0).length;

        // set up for radius finding
        double[] radiusArray = new double[nHoughs];
        for(int n=0; n<nHoughs; n++){
            radiusArray[n] = minRadius + increment*n;
        }

        // holder for average radius per nucleus. [nChannels][nPeaks]
        ArrayList<double[][]> radiusList = new ArrayList<double[][]>();

        // find radii

        float zOffset = 0;

        double[][] radiusPerChannel = new double[nChannels][];

        // set up array for flagging poor peaks
        isPeakGood = new boolean[nPeaksChannel1];
        float[] xChannel1 = xPositions.get(0);
        float[] yChannel1 = yPositions.get(0);

        // test if peaks in channel 2 match well with peaks in channel 1
        for(int ch=1; ch<=nChannels; ch++) {
            float[] thisX = xPositions.get(ch - 1);
            float[] thisY = yPositions.get(ch - 1);
            double[] thisRadiusList = new double[xChannel1.length];

            for (int p = 0; p < xChannel1.length; p++) {
                isPeakGood[p] = true;

                // grab centres from channel 1
                float xC_ = xChannel1[p];
                float yC_ = yChannel1[p];

                if(xC_==0 && yC_==0){
                    //ij maximumfinder is a piece of crap so loads of detections at origin
                    //IJ.log("piece of crap");
                    isPeakGood[p] = false;
                    continue;
                }

                //IJ.log("x: "+xC_+", y: "+yC_);
                if(ch>1){
                    float bestDist = Float.MAX_VALUE;
                    int bestInd = 0;
                    //find closest neighbouring peak
                    for(int i=0; i<thisX.length; i++){
                        if(thisX[i]==0 && thisY[i]==0){
                            continue;
                        }
                        float dist = (float) sqrt((xC_-thisX[i])*(xC_-thisX[i]) + (yC_-thisY[i])*(yC_-thisY[i]));
                        if(dist<bestDist){
                            bestInd = i;
                            bestDist = dist;
                        }
                    }
                    if(bestDist>3){
                        isPeakGood[p] = false;
                        continue;
                    }
                    xC_ = thisX[bestInd];
                    yC_ = thisY[bestInd];

                }

                int xC = (int) xC_;
                int yC = (int) yC_;

                if(xC<minRadius || xC>w-minRadius-1 || yC <minRadius || yC>h-minRadius-1){
                    thisRadiusList[p] = -1;
                    isPeakGood[p] = false;
                    continue;
                }

                if(doZCorrection) {zOffset = pixelsZOffset[xC+yC*w];}

                double[] houghVals = new double[nHoughs];

                for (int r = 1; r <= nHoughs; r++) {

                    int ind = impHoughHS.getStackIndex(ch,r,1);
                    FloatProcessor fp = impHoughHS.getStack().getProcessor(ind).convertToFloatProcessor();

                    int counter = 0;
                    double avVal = 0;

                    for(int y = max(0, yC-averaging); y<=min(yC+averaging, h-1); y++){
                        for (int x = max(0, xC - averaging); x <= min(xC + averaging, w - 1); x++) {
                            avVal += fp.getf(x, y);
                            counter++;
                        }
                    }
                    avVal /= counter;

                    houghVals[r-1] = avVal;

                }

                // fit to z-profile to get radius
                //double[] thisRadiusAndR2 = getRadius(maxRadius, minRadius, radiusArray, houghVals, 0.8);
                //double thisRadius = thisRadiusAndR2[0];
                double thisRadius = getRadius(radiusArray, houghVals, 2);
                //IJ.log("peak "+p+" at ("+xC+", "+yC+") has radius "+thisRadius);
                thisRadiusList[p] = thisRadius;

                if(thisRadius==-1){
                    isPeakGood[p] = false;
                    continue;
                }

                if(doZCorrection && ch==2) {
                    IJ.log("pre-z correction radius = "+thisRadius);
                    double targetRadius = sqrt(thisRadius * pixelSize * thisRadius * pixelSize + zOffset * zOffset);
                    IJ.log("target radius is "+targetRadius);
                    double scaleFactor = targetRadius / (thisRadius * pixelSize);
                    IJ.log("scale factor is "+scaleFactor);
                    thisRadius *= scaleFactor;
                    thisRadius /= pixelSize;
                    IJ.log("radius is now "+ thisRadius);
                }
            }

            radiusPerChannel[ch - 1] = thisRadiusList;

        }
        radiusList.add(radiusPerChannel);


        double meanDiff = 0;
        int counter = 0;

        int nPeaks = 0;
        for(int i=0; i<isPeakGood.length; i++){
            if(isPeakGood[i]) nPeaks++;
        }

        double[] diffs = new double[nPeaks];

        ResultsTable rt = new ResultsTable();

        IJ.run("Tile");

        int i=0;
        int rejectionCount = 0;

        while(i<nPeaks){

            if(!isPeakGood[i]){
                diffs[i] = -1000;
                i++;
                continue;
            }

            double radCh1 = radiusList.get(0)[0][i];
            double radCh2 = radiusList.get(0)[1][i];
            double diff = radCh1 - radCh2;
            if (abs(diff) > maxDisplacementPx) {
                diffs[i] = -1000;
                i++;
                continue;
            }
            if(radCh1<0 || radCh2<0){
                diffs[i] = -1000;
                i++;
                continue;
            }
            diffs[i] = diff;

            double xC = xPositions.get(0)[i]+0.5;
            double yC = yPositions.get(0)[i]+0.5;

            EllipseRoi eRoi = new EllipseRoi(xC-radCh1, yC-radCh1, xC+radCh1, yC+radCh1, 1);
            eRoi.setStrokeColor(Color.cyan);
            eRoi.setStrokeWidth(2);
            imp.setRoi(eRoi);
            impThresholdHS.setRoi(eRoi);

            boolean reject = false;

            if(doManualCuration) {

                NonBlockingGenericDialog gdCheck = new NonBlockingGenericDialog("Check nucleus " + (i + 1) + "/" + nPeaks);
                gdCheck.addMessage("Nucleus at (" + xC + "," + yC + "):\n Radii = " + radCh1 * pixelSize + "nm (ch1) and " + radCh2 * pixelSize +
                        "nm (ch2)");//\n Difference = " + diff * pixelSize + "nm.");
                gdCheck.addCheckbox("Reject nucleus? ", false);
                gdCheck.showDialog();

                reject = gdCheck.getNextBoolean();
            }

            if(!reject) {

                meanDiff += diff;
                counter++;

                rt.incrementCounter();
                rt.addValue("x-position", xC);//*pixelSize);
                rt.addValue("y-position", yC);//*pixelSize);
                rt.addValue("radius channel 1 (nm)", radCh1 * pixelSize);
                rt.addValue("radius channel 2 (nm)", radCh2 * pixelSize);
                rt.addValue("difference", diff * pixelSize);

                rm.addRoi(eRoi);
            }
            else{
                rejectionCount++;
                eRoi.setStrokeColor(Color.RED);
                rm.addRoi(eRoi);
            }

            i++;

        }

        meanDiff/=counter;
        double std = 0;
        for(int n=0; n<nPeaks; n++){
            if(diffs[n] == -1000) continue;
            std += pow(diffs[n]-meanDiff,2);
        }
        std = sqrt(std/counter);

        IJ.log("mean difference is "+meanDiff*pixelSize+"±"+std*pixelSize);
        rt.show("Hough results ("+rejectionCount+" nuclei manually rejected)");

    }

    public static void main(String[] args) throws Exception {
        // set the plugins.dir property to make the plugin appear in the Plugins menu
        // see: https://stackoverflow.com/a/7060464/1207769
        Class<?> clazz = OffsetFinder_.class;
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
