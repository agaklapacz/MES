public class Node {
    private double x;
    private double y;
    private double t;
    private int ownID;
    public boolean BC;

    Node(double x, double y, double t, boolean status, int ownID, Grid grid) {
        this.x = x;
        this.y = y;
        this.t = t;
        this.BC = status;
        this.ownID = ownID;
    }

    public boolean isBC() {
        return BC;
    }

    public double getX() {
        return this.x;
    }

    public double getY() {
        return y;
    }

    public double getT() {
        return t;
    }

    public void setT(double t) {
        this.t = t;
    }

}