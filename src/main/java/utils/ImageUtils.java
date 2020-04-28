package utils;

import ij.ImageStack;
import ij.process.FloatProcessor;

import static java.lang.Math.max;
import static java.lang.Math.sqrt;

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
}
