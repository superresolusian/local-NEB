package gui;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.FolderOpener;
import ij.plugin.PlugIn;

import java.awt.*;
import java.io.File;
import java.io.IOException;

// based on NanoJ BaseDialog, but without aparapi dependencies

public abstract class BaseGUI_ implements PlugIn, DialogListener {

    public NonBlockingGenericDialog gd = null;

    public ImagePlus imp = null;

    protected boolean autoOpenImp = true;

    public String prefsHeader = null, impPath = null;
    protected Prefs prefs = new Prefs();
    protected String arg;

    //-------------------------------------
    int w, h, nChannels, nPixels, nFrames;

    //-------------------------------------

    public void run() {
        run("");
    }

    @Override
    public void run(String arg) {
        this.arg = arg;

        if (!beforeSetupDialog(arg)) return;
        if (autoOpenImp && imp == null) {
            setupImp();
            if (imp == null) return;
        }

        setupDialog();

        if (gd != null) {
            // Add listener to dialog
            gd.addDialogListener(this);
            // Show dialog
            gd.showDialog();
            if (gd.wasCanceled()) {
                return;
            }
        } else {
            if (!loadSettings()) return;
        }

        loadSettings();

        try {
            execute();
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        prefs.savePreferences();
    }

    // ~~~~~PROCESSES~~~~~
    public void getInfo(){
        w = imp.getWidth();
        h = imp.getHeight();
        nChannels = imp.getNChannels();
        nFrames = imp.getNFrames();
        nPixels = w*h;
    }

    // ~~~~~BASIC DIALOG METHODS~~~~~
    abstract public boolean beforeSetupDialog(String arg);

    abstract public void setupDialog();

    abstract public boolean loadSettings();

    abstract public void execute() throws InterruptedException, IOException;

    protected void setupImp() {
        //get open image processor
        imp = WindowManager.getCurrentImage();

        if (imp == null) {
            openImp();
        }
        if (imp == null) return;
    }

    protected void openImp() {

        impPath = IJ.getFilePath("Choose data to load...");
        if (impPath == null || impPath.equals("")) return;

        else if (impPath.endsWith(".tif")) {
            if (impPath.contains("00")) {
                //open image sequence in directory selected and show
                imp = FolderOpener.open(new File(impPath).getParent());
            } else {
                imp = IJ.openImage(impPath);
            }
        } else {
            imp = IJ.openImage(impPath);
        }
        if (imp == null) return;
        imp.show();
    }

    public String getClassName() {
        return this.getClass().getName();
    }

    // ~~~~~PREFS HANDLING~~~~~

    public int getPrefs(String key, int defaultValue) {
        if (prefsHeader == null) prefsHeader = getClassName();
        return (int) prefs.get(prefsHeader + "." + key, defaultValue);
    }

    public float getPrefs(String key, float defaultValue) {
        if (prefsHeader == null) prefsHeader = getClassName();
        return (float) prefs.get(prefsHeader + "." + key, defaultValue);
    }

    public double getPrefs(String key, double defaultValue) {
        if (prefsHeader == null) prefsHeader = getClassName();
        return (double) prefs.get(prefsHeader + "." + key, defaultValue);
    }

    public boolean getPrefs(String key, boolean defaultValue) {
        if (prefsHeader == null) prefsHeader = getClassName();
        return prefs.get(prefsHeader + "." + key, defaultValue);
    }

    public String getPrefs(String key, String defaultValue) {
        if (prefsHeader == null) prefsHeader = getClassName();
        return prefs.get(prefsHeader + "." + key, defaultValue);
    }

    public void setPrefs(String key, int value) {
        if (prefsHeader == null) prefsHeader = getClassName();
        prefs.set(prefsHeader + "." + key, value);
    }

    public void setPrefs(String key, float value) {
        if (prefsHeader == null) prefsHeader = getClassName();
        prefs.set(prefsHeader + "." + key, value);
    }

    public void setPrefs(String key, double value) {
        if (prefsHeader == null) prefsHeader = getClassName();
        prefs.set(prefsHeader + "." + key, value);
    }

    public void setPrefs(String key, boolean value) {
        if (prefsHeader == null) prefsHeader = getClassName();
        prefs.set(prefsHeader + "." + key, value);
    }

    public void setPrefs(String key, String value) {
        if (prefsHeader == null) prefsHeader = getClassName();
        prefs.set(prefsHeader + "." + key, value);
    }

    public void savePrefs() {
        prefs.savePreferences();
    }

    @Override
    public boolean dialogItemChanged(GenericDialog genericDialog, AWTEvent awtEvent) {
        return false;
    }
}
