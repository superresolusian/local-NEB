package utils;

import ij.ImageStack;
import ij.plugin.filter.Binary;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

import java.awt.*;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

import static ij.process.ImageProcessor.NO_LUT_UPDATE;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.sqrt;
import static utils.ArrayUtils.ravel;

public class ImageUtils {

    public static float[] getFpStats(FloatProcessor fp) {
        float[] pixels = (float[]) fp.getPixels();
        int nPixels = pixels.length;
        float mean = 0;
        for (int i = 0; i < nPixels; i++) mean += pixels[i] / nPixels;
        float std = 0;
        for (int i = 0; i < nPixels; i++) std += (pixels[i] - mean) * (pixels[i] - mean) / nPixels;
        std = (float) sqrt(std);

        return new float[]{mean, std};

    }

    public static FloatProcessor maxProject(ImageStack ims) {
        int w = ims.getWidth();
        int h = ims.getHeight();
        int nSlices = ims.getSize();

        float[] pixels = new float[w * h];

        for (int n = 1; n <= nSlices; n++) {
            FloatProcessor fp = ims.getProcessor(n).convertToFloatProcessor();
            float[] pixels_ = (float[]) fp.getPixels();
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    int p = x + y * w;
                    pixels[p] = max(pixels[p], pixels_[p]);
                }
            }
        }
        return new FloatProcessor(w, h, pixels);

    }

    public static FloatProcessor averageProject(ImageStack ims){
        int w = ims.getWidth();
        int h = ims.getHeight();
        int nSlices = ims.getSize();

        float[] pixels = new float[w*h];

        for(int n=1; n<=nSlices; n++){
            FloatProcessor fp = ims.getProcessor(n).convertToFloatProcessor();
            for(int y=0; y<h; y++){
                for(int x=0; x<w; x++){
                    int p = x+y*w;
                    pixels[p] += fp.getf(x, y)/nSlices;
                }
            }
        }

        FloatProcessor fp_ = new FloatProcessor(w, h, pixels);

        return fp_;
    }

    public static ImageStack groupedZProject(ImageStack ims, int groupSize){
        int nGroups = ims.getSize()/groupSize;
        int w = ims.getWidth();
        int h = ims.getHeight();
        ImageStack imsOut = new ImageStack(w, h, nGroups);

        for(int i=0; i<nGroups; i++){
            float[] pixels = new float[w*h];
            for(int n=0; n<groupSize; n++){
                FloatProcessor fp = ims.getProcessor(n+i*groupSize+1).convertToFloatProcessor();
                float[] pixels_ = (float[]) fp.getPixels();
                for(int y=0; y<h; y++){
                    for(int x=0; x<w; x++){
                        int p = x+y*w;
                        pixels[p] = max(pixels[p], pixels_[p]);
                    }
                }
            }
            imsOut.setProcessor(new FloatProcessor(w, h, pixels), i+1);
        }
        return imsOut;
    }

    public static double getThresholdForIp(ImageProcessor ip){
        ip.setAutoThreshold("Default", true, NO_LUT_UPDATE);
        double threshold = ip.getAutoThreshold();
        return threshold;
    }

    public static double getMaxOfImage(ImageProcessor ip){
        double imageMaxVal = 0;
        for(int i=0; i<ip.getPixelCount(); i++) imageMaxVal = max(imageMaxVal, ip.get(i));
        return imageMaxVal;
    }

    public static ImageProcessor getSkeletonImage(ImageProcessor ip, double blurSize, boolean invertLUT){
        Binary binary = new Binary();

        if(blurSize!=0){
            ip.blurGaussian(blurSize);
        }
        ip.setAutoThreshold("Triangle", true, NO_LUT_UPDATE);
        ip.autoThreshold();
        ImageProcessor ipMask = ip.createMask();
        if(invertLUT) ipMask.invertLut();

        binary.setup("fill holes", null);
        binary.run(ipMask);
        binary.setup("skeletonize", null);
        binary.run(ipMask);

        return ipMask;
    }

    public static ImageStack normaliseToMax(ImageStack ims){
        int w = ims.getWidth();
        int h = ims.getHeight();
        int nFrames = ims.getSize();

        ImageStack imsOut = new ImageStack(w, h, nFrames);

        for(int n=1; n<=nFrames; n++){
            FloatProcessor fp_ = ims.getProcessor(n).convertToFloatProcessor();
            float[] pixels = (float[]) fp_.getPixels();
            float maxVal = 0;
            for(int i=0; i<pixels.length; i++) maxVal = max(maxVal, pixels[i]);

            float[] newPixels = new float[pixels.length];
            for(int i=0; i<pixels.length; i++) newPixels[i] = (1000 * pixels[i]/maxVal);

            FloatProcessor fp = new FloatProcessor(w, h, newPixels);
            imsOut.setProcessor(fp, n);

        }

        return imsOut;
    }

}
