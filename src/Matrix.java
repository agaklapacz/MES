public class Matrix
{
    private double[][] localH;
    private double[][] localC;
    private double[][] globalH;
    private double[][] globalC;
    private double[][] H;
    private double[] P;
    private double[] t0;
    private double[] t1;
    private double[][] Hbc;
    private double[] VP;
    private double[] GP;
    private double[][] globalHbc;
    private double[][] Pcx;
    private GlobalData globalData;
    private Grid grid;
    private Node[] nodes;
    final double positiveValue = 1 / (Math.sqrt(3));
    final double negativeValue = -(1 / (Math.sqrt(3)));
    private static final double EPSILON = 1e-10;

    public Matrix(GlobalData globalData, Grid grid) {
        this.globalData = globalData;
        localH = new double[4][4];
        localC = new double[4][4];
        nodes = new Node[4];
        this.globalH = new double[globalData.getNumberOfN()][globalData.getNumberOfN()];
        this.globalC = new double[globalData.getNumberOfN()][globalData.getNumberOfN()];
        this.globalHbc = new double[globalData.getNumberOfN()][globalData.getNumberOfN()];
        this.Hbc = new double[4][4];
        this.VP = new double[4];
        this.H = new double[globalData.getNumberOfN()][globalData.getNumberOfN()];
        this.P = new double[globalData.getNumberOfN()];
        this.GP = new double[globalData.getNumberOfN()];
        this.t0 = new double[globalData.getNumberOfN()];
        this.t1 = new double[globalData.getNumberOfN()];
        this.grid = grid;
        MatrixPcx();
    }

    private void MatrixPcx() {
        Pcx = new double[4][2];
        Pcx[0][0] = negativeValue;
        Pcx[0][1] = negativeValue;
        Pcx[1][0] = positiveValue;
        Pcx[1][1] = negativeValue;
        Pcx[2][0] = positiveValue;
        Pcx[2][1] = positiveValue;
        Pcx[3][0] = negativeValue;
        Pcx[3][1] = positiveValue;
    }

    public void generateMatrices() {
        double[][] dNdKsi = new double[4][4];
        double[][] dNdEta = new double[4][4];
        double[][] N = new double[4][4];
        double[] X;
        double[] Y;
        double[] tmp;
        Element element;
        int[] ID;
        double[][] dNdX = new double[4][4];
        double[][] dNdY = new double[4][4];

        for (int i = 0; i < 4; i++) {
            double[] mk = MatrixN_Ksi(Pcx[i][1]);
            dNdKsi[i][0] = mk[0];
            dNdKsi[i][1] = mk[1];
            dNdKsi[i][2] = mk[2];
            dNdKsi[i][3] = mk[3];
        }

        for (int i = 0; i < 4; i++) {
            double[] me = MatrixN_Eta(Pcx[i][0]);
            dNdEta[i][0] = me[0];
            dNdEta[i][1] = me[1];
            dNdEta[i][2] = me[2];
            dNdEta[i][3] = me[3];
        }

        for (int i = 0; i < 4; i++) {
            double[] mx = MatrixNx(Pcx[i][0], Pcx[i][1]);
            N[i][0] = mx[0];
            N[i][1] = mx[1];
            N[i][2] = mx[2];
            N[i][3] = mx[3];
        }

        for (int i = 0; i < globalData.getNumberOfE(); i++) {
            X = new double[4];
            Y = new double[4];
            tmp = new double[4];
            element = grid.getElements()[i];
            ID = element.getID();

            for (int j = 0; j < 4; j++) {
                VP[j] = 0;
                for (int k = 0; k < 4; k++) {
                    localH[j][k] = 0;
                    localC[j][k] = 0;
                    Hbc[j][k] = 0;
                }
            }

            for (int j = 0; j < 4; j++) {
                nodes[j] = grid.getNodes()[ID[j]];
            }

            for (int j = 0; j < 4; j++) {
                X[j] = grid.getNodes()[ID[j]].getX();
                Y[j] = grid.getNodes()[ID[j]].getY();
                tmp[j] = grid.getNodes()[ID[j]].getT();
            }

            for (int integrationPoint = 0; integrationPoint < 4; integrationPoint++) {
                double dXdKsi = 0;
                double dXdEta = 0;
                double dYdKsi = 0;
                double dYdEta = 0;

                for (int j = 0; j < 4; j++) {
                    dXdKsi += dNdKsi[integrationPoint][j] * X[j];
                    dXdEta += dNdEta[integrationPoint][j] * X[j];
                    dYdKsi += dNdKsi[integrationPoint][j] * Y[j];
                    dYdEta += dNdEta[integrationPoint][j] * Y[j];
                }

                double detJ = dXdKsi * dYdEta - dXdEta * dYdKsi;

                for (int l = 0; l < 4; l++) {
                    dNdX[integrationPoint][l] = (dYdEta * 1 / detJ * dNdKsi[integrationPoint][l] + (-dYdKsi) * 1 / detJ * dNdEta[integrationPoint][l]);
                    dNdY[integrationPoint][l] = ((-dXdEta) * 1 / detJ * dNdKsi[integrationPoint][l] + dXdKsi * 1 / detJ * dNdEta[integrationPoint][l]);
                }

                for (int m = 0; m < 4; m++) {
                    for (int n = 0; n < 4; n++) {
                        localH[m][n] += (globalData.getK() * (dNdX[integrationPoint][m] * dNdX[integrationPoint][n] + dNdY[integrationPoint][m] * dNdY[integrationPoint][n]) * detJ);
                        localC[m][n] += globalData.getHeat() * globalData.getDensity() * detJ * (N[integrationPoint][m] * N[integrationPoint][n]);
                    }
                }
            }

            double w1 = 1;
            double w2 = 1;
            double detJ;

            if (nodes[0].isBC() && nodes[1].isBC()) {

                detJ = 0.5 * Math.sqrt(Math.pow(nodes[1].getX() - nodes[0].getX(), 2) + Math.pow(nodes[1].getY() - nodes[0].getY(), 2));

                double[] shapeVector1 = new double[4];
                shapeVector1 = MatrixNx(negativeValue,-1);

                double[] shapeVector2 = new double[4];
                shapeVector2 = MatrixNx(positiveValue,-1);

                for (int b = 0; b < 4; b ++){
                    VP[b] += ((shapeVector1[b] + shapeVector2[b]) * detJ * (-globalData.getAlfa()) * globalData.getAmbientT());
                    for (int c = 0; c < 4; c ++){
                        Hbc[b][c] += globalData.getAlfa() * ((w1 * shapeVector1[b] * shapeVector1[c]) + (w2 * shapeVector2[b] * shapeVector2[c])) * detJ;
                    }
                }
            }

            if (nodes[1].isBC() && nodes[2].isBC()) {
                detJ = 0.5 * Math.sqrt(Math.pow(nodes[2].getX() - nodes[1].getX(), 2) + Math.pow(nodes[2].getY() - nodes[1].getY(), 2));

                double[] shapeVector1 = new double[4];
                shapeVector1 = MatrixNx(1, negativeValue);

                double[] shapeVector2 = new double[4];
                shapeVector2 = MatrixNx(1, positiveValue);

                for (int b = 0; b < 4; b ++){
                    VP[b] += ((shapeVector1[b] + shapeVector2[b]) * detJ * (-globalData.getAlfa()) * globalData.getAmbientT());
                    for (int c = 0; c < 4; c ++){
                        Hbc[b][c] += globalData.getAlfa() * ((w1 * shapeVector1[b] * shapeVector1[c]) + (w2 * shapeVector2[b] * shapeVector2[c])) * detJ;
                    }
                }
            }

            if (nodes[2].isBC() && nodes[3].isBC()){
                detJ = 0.5 * Math.sqrt(Math.pow(nodes[2].getX() - nodes[3].getX(), 2) + Math.pow(nodes[2].getY() - nodes[3].getY(), 2));

                double[] shapeVector1 = new double[4];
                shapeVector1 = MatrixNx(negativeValue, 1);

                double[] shapeVector2 = new double[4];
                shapeVector2 = MatrixNx(positiveValue, 1);

                for (int b = 0; b < 4; b ++){
                    VP[b] += ((shapeVector1[b] + shapeVector2[b]) * detJ * (-globalData.getAlfa()) * globalData.getAmbientT());
                    for (int c = 0; c < 4; c ++){
                        Hbc[b][c] += globalData.getAlfa() * ((w1 * shapeVector1[b] * shapeVector1[c]) + (w2 * shapeVector2[b] * shapeVector2[c])) * detJ;
                    }
                }
            }

            if (nodes[3].isBC() && nodes[0].isBC()){
                detJ = 0.5 * Math.sqrt(Math.pow(nodes[3].getX() - nodes[0].getX(), 2) + Math.pow(nodes[3].getY() - nodes[0].getY(), 2));

                double[] shapeVector1 = new double[4];
                shapeVector1 = MatrixNx(-1, positiveValue);

                double[] shapeVector2 = new double[4];
                shapeVector2 = MatrixNx(-1, negativeValue);

                for (int b = 0; b < 4; b ++){
                    VP[b] += ((shapeVector1[b] + shapeVector2[b]) * detJ * (-globalData.getAlfa()) * globalData.getAmbientT());
                    for (int c = 0; c < 4; c ++){
                        Hbc[b][c] += globalData.getAlfa() * ((w1 * shapeVector1[b] * shapeVector1[c]) + (w2 * shapeVector2[b] * shapeVector2[c])) * detJ;
                    }
                }
            }
            /*System.out.println("Macierz Hbc dla elementu: " + i);
            for (int m = 0; m < 4; m++){
                for (int n = 0; n < 4; n ++){
                    System.out.print(Hbc[m][n] + " ");
                }
                System.out.println("");
            }

            System.out.println("Wektor obciążeń P dla elementu: " + i);
            for (int m = 0; m < 4; m++){
                System.out.print(VP[m] + " ");
            }
            System.out.println("");*/

            element.setLocalH(localH);
            //element.printLocalH();
            element.setLocalC(localC);
            //element.printLocalC();

            for (int k = 0; k < 4; k++) {
                GP[ID[k]] += VP[k];
                for (int l = 0; l < 4; l++) {
                    globalH[ID[k]][ID[l]] += localH[k][l];
                    globalC[ID[k]][ID[l]] += localC[k][l];
                    globalHbc[ID[k]][ID[l]] += Hbc[k][l];
                }
            }
        }
        /*for (int j = 0; j < globalData.getNumberOfN(); j++) {
            for (int k = 0; k < globalData.getNumberOfN(); k++) {
                System.out.print(globalHbc[j][k] + " ");
            }
        }
            System.out.println("");*/
        /*System.out.println("Agregacja wektora P: ");
        for (int j = 0; j < globalData.getNumberOfN(); j++){
            System.out.print(GP[j] + " ");
        }
        System.out.println("");*/

        for (int i = 0; i < globalData.getNumberOfN(); i ++){
            for (int j = 0; j < globalData.getNumberOfN(); j++){
                globalH[i][j] += globalHbc[i][j];
            }
        }
        /*System.out.println("Agregacja newH: ");
        for (int j = 0; j < globalData.getNumberOfN(); j++) {
            for (int k = 0; k < globalData.getNumberOfN(); k++) {
                System.out.print(globalH[j][k] + " ");
            }
            System.out.println("");
        }*/
    }

    void nextStep() {
        System.out.println(" ");
        for (double w = 0; w <= globalData.getSimulationTime(); w += globalData.getStepTime()){
            for (int i = 0; i < globalData.getNumberOfN(); i ++){
                for (int j = 0; j < globalData.getNumberOfN(); j++){
                    H[i][j] = globalH[i][j] + globalC[i][j]/globalData.getStepTime();
                }
            }
            for (int i = 0 ; i < globalData.getNumberOfN(); i ++){
                t0[i] = grid.getNodes()[i].getT();
            }

            for (int i = 0; i < globalData.getNumberOfN(); i ++){
                double sum=0;
                for (int j = 0; j < globalData.getNumberOfN(); j++){
                    sum += globalC[j][i] * t0[j];
                }
                P[i] = sum / globalData.getStepTime() - GP[i];
            }
            t1 = lslove(H,P);
            for (int i = 0; i < globalData.getNumberOfN(); i ++){
                grid.getNodes()[i].setT(t1[i]);
            }
            System.out.print("Simulation time:  " + w + "         Min. temperature: " + getMinTemperature() + "         Max. temperature: " + getMaxTemperature());
            System.out.println(" ");
        }
    }

    public double[] lslove(double[][] A, double[] b){
        int n = b.length;

        for (int p = 0; p < n; p++) {

            int max = p;
            for (int i = p + 1; i < n; i++) {
                if (Math.abs(A[i][p]) > Math.abs(A[max][p])) {
                    max = i;
                }
            }
            double[] temp = A[p]; A[p] = A[max]; A[max] = temp;
            double   t    = b[p]; b[p] = b[max]; b[max] = t;

            if (Math.abs(A[p][p]) <= EPSILON) {
                throw new ArithmeticException("Matrix is singular or nearly singular");
            }

            for (int i = p + 1; i < n; i++) {
                double alpha = A[i][p] / A[p][p];
                b[i] -= alpha * b[p];
                for (int j = p; j < n; j++) {
                    A[i][j] -= alpha * A[p][j];
                }
            }
        }

        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }
        return x;
    }

    double getMinTemperature(){
        double n, wynik;
        wynik = t0[0];
        for (int i = 1; i < t0.length; i ++){
            if (wynik > t0[i])
                wynik = t0[i];
        }
        return wynik;
    }

    double getMaxTemperature(){
        double n, wynik;
        wynik = t0[0];
        for (int i = 1; i < t0.length; i ++){
            if (wynik < t0[i])
                wynik = t0[i];
        }
        return wynik;
    }

    double[] MatrixN_Ksi(double eta) {
        double[] N_Ksi = new double[4];
        N_Ksi[0] = -0.25 * (1 - eta);
        N_Ksi[1] = 0.25 * (1 - eta);
        N_Ksi[2] = 0.25 * (1 + eta);
        N_Ksi[3] = -0.25 * (1 + eta);
        return (N_Ksi);
    }

    double[] MatrixN_Eta(double ksi) {
        double[] N_Eta = new double[4];
        N_Eta[0] = -0.25 * (1 - ksi);
        N_Eta[1] = -0.25 * (1 + ksi);
        N_Eta[2] = 0.25 * (1 + ksi);
        N_Eta[3] = 0.25 * (1 - ksi);
        return (N_Eta);
    }

    double[] MatrixNx(double ksi, double eta) {
        double[] Nx = new double[4];
        Nx[0] = 0.25 * (1 - ksi) * (1 - eta);
        Nx[1] = 0.25 * (1 + ksi) * (1 - eta);
        Nx[2] = 0.25 * (1 + ksi) * (1 + eta);
        Nx[3] = 0.25 * (1 - ksi) * (1 + eta);
        return (Nx);
    }

    /*public Grid getGrid() {
        return grid;
    }*/
}

