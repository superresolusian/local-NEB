package utils;

import ij.io.Opener;
import ij.plugin.frame.RoiManager;

public class IOUtils {

    public static void saveRoisAsZip(RoiManager rm, String dir, String title){
        if(title.endsWith(".tif")){
            title = title.substring(0, title.length() - 4);
        }
        String savePath = dir+title+"_RoiSet.zip";
        rm.runCommand("Save", savePath);
    }

    public static void openRoisFromZip(String path){
        RoiManager rm = new RoiManager().getInstance();
        rm.close();

        new Opener().openZip(path);
    }

}
