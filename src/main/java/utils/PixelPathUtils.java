package utils;

import ij.gui.Line;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

import java.awt.*;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

import static java.lang.Math.*;
import static utils.ArrayUtils.ravel;
import static utils.ArrayUtils.unravel;

public class PixelPathUtils {

    public static short[] connectedComponents(ImageProcessor ip){
        int w = ip.getWidth();
        int h = ip.getHeight();
        int nPixels = w*h;

        short[] data = new short[nPixels];
        for(int i=0; i<nPixels; i++) data[i] = (short) ip.get(i);

        short[] components = new short[data.length];
        short id = 1;
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                int p = ravel(x, y, w);
                // foreground pixel that hasn't been labeled yet
                if (data[p] > 0 && components[p] == 0) {
                    labelComponentAt(id, x, y, components, data, w, h);
                    id += 1;
                }
            }
        }
        return components;
    }

    private static void labelComponentAt(short id, int x, int y, short[] components, short[] data, int w, int h) {
        components[ravel(x, y, w)] = id;
        LinkedList<Point> queue = new LinkedList<Point>();
        queue.add(new Point(x, y));
        while (!queue.isEmpty()) {
            Point p = queue.removeFirst();
            java.util.List<Point> neighbors = unlabeledNeighbors(p.x, p.y, components, data, w, h);
            for (Point n : neighbors)
                components[ravel(n.x, n.y, w)] = id;
            queue.addAll(neighbors);
        }
    }

    private static java.util.List<Point> unlabeledNeighbors(int x, int y, short[] components, short[] data, int w, int h) {
        assert components[ravel(x,y,w)] > 0;
        List<Point> neighbors = new LinkedList<Point>();
        for (int j = max(y - 1, 0); j < min(y + 2, h); j++) {
            for (int i = max(x - 1, 0); i < min(x + 2, w); i++) {
                int p = ravel(i, j, w);
                if (data[p] > 0 && components[p] == 0) {
                    neighbors.add(new Point(i, j));
                }
            }
        }
        return neighbors;
    }

    public static class regionProps {
        public int area;
        public Rectangle boundingBox;
        public double[] centroid;
        public int[] endPoints;
        public int[] branchPoints;
        public ArrayList<int[]> endpointsList;
        public ArrayList<int[]> branchpointsList;

        public regionProps(int area, Rectangle boundingBox, double[] centroid, int[] endPoints, int[] branchPoints){
            this.area = area;
            this.boundingBox = boundingBox;
            this.centroid = centroid;
            this.endPoints = endPoints;
            this.branchPoints = branchPoints;
        }

        public regionProps(int area, Rectangle boundingBox, double[] centroid, ArrayList<int[]> endpointsList, ArrayList<int[]> branchpointsList){
            this.area = area;
            this.boundingBox = boundingBox;
            this.centroid = centroid;
            this.endpointsList = endpointsList;
            this.branchpointsList = branchpointsList;
        }


    }

    public static regionProps getProps(int n, ImageProcessor ip){

        int w = ip.getWidth();
        int h = ip.getHeight();

        int area = 0;
        int x0 = w-1, y0 = h-1;
        int x1 = 0, y1 = 0;
        int sumX = 0, sumY = 0;
        LinkedHashMap<int[], Integer> neighbourCount = new LinkedHashMap<int[], Integer>();
        ArrayList<int[]> endPointList = new ArrayList<int[]>();
        ArrayList<int[]> branchPointList = new ArrayList<>();

        int minNeighbours = 8;

        for(int y=0; y<h; y++){
            for(int x=0; x<w; x++){
                if(ip.getf(x, y)==n){
                    area++;
                    x0 = min(x, x0);
                    y0 = min(y, y0);
                    x1 = max(x, x1);
                    y1 = max(y, y1);
                    sumX += x;
                    sumY += y;

                    int nNeighbours = 0;
                    for(int y_=max(0, y-1); y_<=min(h-1,y+1); y_++){
                        for(int x_=max(0,x-1); x_<=min(w-1, x+1); x_++){
                            if(y_==y && x_==x) continue;
                            if(ip.getf(x_,y_)>0) nNeighbours++;
                        }
                    }
                    neighbourCount.put(new int[] {x, y}, nNeighbours);
                    minNeighbours = min(minNeighbours, nNeighbours);
                }
            }
        }

        for(int[] points:neighbourCount.keySet()){
            int nNeighbours = neighbourCount.get(points);
            if(nNeighbours==minNeighbours) endPointList.add(points);
            if(nNeighbours>2) branchPointList.add(points);
        }

        Rectangle bbox = new Rectangle(x0, y0, x1-x0, y1-y0);

        return new regionProps(area, bbox, new double[]{sumX/area, sumY/area}, endPointList, branchPointList);
    }

    public static regionProps getProps(int n, int w, int h, short[] mask){
        int area = 0;
        int x0 = w-1, y0 = h-1;
        int x1 = 0, y1 = 0;
        int sumX = 0, sumY = 0;
        ArrayList<int[]> neighbourCount = new ArrayList<int[]>();
        ArrayList<Integer> endPoint = new ArrayList<Integer>();
        ArrayList<Integer> branchPoint = new ArrayList<>();

        int minNeighbours = 8;

        for(int y=0; y<h; y++){
            for(int x=0; x<w; x++){
                int p = y*w+x;
                if(mask[p]==n){
                    area++;
                    x0 = min(x, x0);
                    y0 = min(y, y0);
                    x1 = max(x, x1);
                    y1 = max(y, y1);
                    sumX += x;
                    sumY += y;

                    int nNeighbours = 0;
                    for(int y_=max(0, y-1); y_<=min(h-1,y+1); y_++){
                        for(int x_=max(0,x-1); x_<=min(w-1, x+1); x_++){
                            if(y_==y && x_==x) continue;
                            int p_ = y_*w + x_;
                            if(mask[p_]>0) nNeighbours++;
                        }
                    }
                    neighbourCount.add(new int[] {p, nNeighbours});
                    minNeighbours = min(minNeighbours, nNeighbours);
                }
            }
        }

        for(int i=0; i<area; i++){
            int[] neighbourInfo = neighbourCount.get(i);
            if(neighbourInfo[1]==minNeighbours){
                endPoint.add(neighbourInfo[0]);
            }
            if(neighbourInfo[1]>2){
                branchPoint.add(neighbourInfo[1]);
            }
        }

        int nEndpoints = endPoint.size();
        int[] endpointArray = new int[nEndpoints];
        for(int i=0; i<nEndpoints; i++) endpointArray[i] = endPoint.get(i);

        int nBranchpoints = branchPoint.size();
        int[] branchpointArray = new int[nBranchpoints];
        for(int i=0; i<nBranchpoints; i++) branchpointArray[i] = branchPoint.get(i);


        Rectangle bbox = new Rectangle(x0, y0, x1-x0, y1-y0);

        return new regionProps(area, bbox, new double[]{sumX/area, sumY/area}, endpointArray, branchpointArray);
    }

    public static LinkedHashMap<Integer, int[]> findNodeNeighbours(ShortProcessor spSegment, ArrayList<Integer> nodeList, int s, int area){
        int w = spSegment.getWidth();
        int h = spSegment.getHeight();

        // find neighbours for each node
        LinkedHashMap<Integer, int[]> nodeNeighbours = new LinkedHashMap<Integer, int[]>();

        for (int i = 0; i < area; i++) {
            int nodeIndex = nodeList.get(i);
            int x = nodeIndex % w;
            int y = nodeIndex / h;
            ArrayList<Integer> thisNeighbours = new ArrayList<Integer>();

            for (int y_ = max(0, y - 1); y_ <= min(h - 1, y + 1); y_++) {
                for (int x_ = max(0, x - 1); x_ <= min(w - 1, x + 1); x_++) {
                    int p = ravel(x_, y_, w);
                    if (spSegment.get(p) == s + 1) thisNeighbours.add(p);
                }
            }
            int nNeighbours = thisNeighbours.size();
            int[] thisNeighboursArray = new int[nNeighbours];
            for (int j = 0; j < nNeighbours; j++) thisNeighboursArray[j] = thisNeighbours.get(j);

            nodeNeighbours.put(nodeIndex, thisNeighboursArray);
        }

        return nodeNeighbours;
    }

    public static int[] getPath(regionProps props, int[] endpoint1, int[] endpoint2, ShortProcessor spSegment, int s){
        int w = spSegment.getWidth();

        // get list of pixels belonging to skeleton
        Rectangle bounds = props.boundingBox;
        ArrayList<Integer> nodeList = new ArrayList<Integer>();
        for (int y = bounds.y; y <= bounds.y + bounds.height; y++) {
            for (int x = bounds.x; x <= bounds.x + bounds.width; x++) {
                int p = ravel(x, y, w);
                if (spSegment.get(p) == s + 1) nodeList.add(p);
            }
        }

        int area = props.area;
        LinkedHashMap<Integer, int[]> nodeNeighbours = findNodeNeighbours(spSegment, nodeList, s, area);

        int p1 = ravel(endpoint1[0], endpoint1[1], w);
        int p2 = ravel(endpoint2[0], endpoint2[1], w);

        int[] path = findPath(p1, p2, nodeNeighbours);

        return path;
    }

    public static Object[] getPath(regionProps props, ShortProcessor spSegment, int s) {
        int w = spSegment.getWidth();

        // get list of pixels belonging to skeleton
        Rectangle bounds = props.boundingBox;
        ArrayList<Integer> nodeList = new ArrayList<Integer>();
        for (int y = bounds.y; y <= bounds.y + bounds.height; y++) {
            for (int x = bounds.x; x <= bounds.x + bounds.width; x++) {
                int p = ravel(x, y, w);
                if (spSegment.get(p) == s + 1) nodeList.add(p);
            }
        }


        // get endpoints of segment
        ArrayList<int[]> endpoints = props.endpointsList;
        int p1 = endpoints.get(0)[0] + endpoints.get(0)[1] * w;
        int p2 = endpoints.get(1)[0] + endpoints.get(1)[1] * w;
        int[] endpoints_ = new int[]{p1, p2};

        // find neighbours of pixels in nodes

        int area = props.area;
        LinkedHashMap<Integer, int[]> nodeNeighbours = findNodeNeighbours(spSegment, nodeList, s, area);

        int[] path = findPath(p1, p2, nodeNeighbours);

        return new Object[]{endpoints_, path};
    }

    public static int[] findPath(int start, int end, LinkedHashMap<Integer, int[]> nodeNeighbours){

        // list of visited nodes
        ArrayList<Integer> visitedNodes = new ArrayList<Integer>();

        // queue
        ArrayList<int[]> Q = new ArrayList<int[]>();

        // deal with start point
        Q.add(new int[]{start});

        // initialise path
        int[] path = new int[0];

        while (Q.size() > 0) {
            // get first path in queue then remove from queue
            path = Q.get(0);
            Q.remove(0);

            // get last node in path
            int node = path[path.length - 1];

            // check if this is the end
            if (node == end) break;

            // check we haven't already visited this node
            if (visitedNodes.contains(node)) continue;

            // loop through node neighbours
            for (int nodeNeighbour : nodeNeighbours.get(node)) {

                // add neighbour to current path as a potential route
                int[] newPath = new int[path.length + 1];
                for (int i = 0; i < path.length; i++) newPath[i] = path[i];
                newPath[path.length] = nodeNeighbour;

                // add to queue to test
                Q.add(newPath);
            }

            visitedNodes.add(node);
        }

        return path;

    }

    public static double getDistance(int[] p1, int[] p2) {
        int x1 = p1[0];
        int y1 = p1[1];
        int x2 = p2[0];
        int y2 = p2[1];

        return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
    }

    public static int[] testSkeleton(LinkedHashMap<int[], int[]> pathList, ArrayList<int[]> maximaCoordinates, int w, double linkingDistance){

        int[] endpointPair = (int[]) pathList.keySet().toArray()[0];

        int e1 = endpointPair[0];
        int e2 = endpointPair[1];

        int[] e1_ = new int[]{e1%w, e1/w};
        int[] e2_ = new int[]{e2%w, e2/w};

        boolean foundMatchesFlag = false;

        // distance from endpoint to closest hough maximum
        double minDist1 = Double.MAX_VALUE;
        double minDist2 = Double.MAX_VALUE;
        int centre1 = 0, centre2 = 0;

        for(int i=0; i<maximaCoordinates.size()-1; i++){
            int[] centre = maximaCoordinates.get(i);
            double dist1 = getDistance(centre, e1_);
            double dist2 = getDistance(centre, e2_);
            if(dist1<minDist1){
                centre1 = i;
                minDist1 = dist1;
            }
            if(dist2<minDist2){
                centre2 = i;
                minDist2 = dist2;
            }
        }

        // test if centres are close enough to centres
        if(minDist1<linkingDistance && minDist2<linkingDistance) foundMatchesFlag = true;

        if(!foundMatchesFlag){
            centre1 = -1;
            centre2 = -1;
        }

        return new int[]{centre1, centre2};
    }

    public static double getAngleFromBoundingBox(Rectangle bbox){
        int dx = bbox.width;
        int dy = bbox.height;
        return atan(dy/dx);
    }

    public static double getAngleBetweenPoints(int[] e1, int[] e2){
        // convention: work from left to right
        Line line = new Line(e1[0], e1[1], e2[0], e2[1]);
        if(e1[0]>e2[0]) line = new Line(e2[0], e2[1], e1[0], e1[1]);
        return toRadians(line.getAngle());
    }

    public static ArrayList<int[]> pruneEndpoints(regionProps skeletonprops){
        ArrayList<int[]> endPoints = skeletonprops.endpointsList;
        ArrayList<int[]> branchPoints = skeletonprops.branchpointsList;

        while(endPoints.size()>2) {
            // find endpoint with shortest distance to branchpoint
            double shortestBranchLength = Double.MAX_VALUE;
            int[] endpointOfShortestBranch = new int[0];

            for (int[] endPoint : endPoints) {
                for (int[] branchPoint : branchPoints) {
                    double dist = getDistance(endPoint, branchPoint);
                    if (dist < shortestBranchLength) {
                        shortestBranchLength = dist;
                        endpointOfShortestBranch = endPoint;
                    }

                }
            }
            endPoints.remove(endpointOfShortestBranch);
        }

        // sort by x
        ArrayList<int[]> sortedEndpoints = new ArrayList<>(2);
        if(endPoints.get(0)[0]<endPoints.get(1)[0]){
            sortedEndpoints.set(0, endPoints.get(0));
            sortedEndpoints.set(1, endPoints.get(1));
        }
        else{
            sortedEndpoints.set(0, endPoints.get(1));
            sortedEndpoints.set(1, endPoints.get(0));
        }

        return sortedEndpoints;
    }

    public static int[] pathExcludeCircle(int[] path, int[] c1, double r1, int[] c2, double r2, int w){
        ArrayList<Integer> newPath = new ArrayList<>();
        for(int p:path){
            int[] p_ = unravel(p, w);
            if(getDistance(p_, c1)<r1 || getDistance(p_, c2)<r2) continue;
            newPath.add(p);
        }
        int[] newPathArray = new int[newPath.size()];
        for(int i=0; i<newPath.size(); i++) newPathArray[i] = newPath.get(i);

        return newPathArray;
    }

}
