
public class Main
{
    public static void main(String[] args) {
        Grid grid = new Grid();
        grid.gridGeneration();
        grid.showGird();
        GlobalData globalData = new GlobalData();
        Matrix matrix = new Matrix(globalData, grid);
        matrix.generateMatrices();
        matrix.nextStep();
    }
}