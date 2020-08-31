public class GlobalData
{
    private double initialT;  //temperatura poczatkowa
    private double alfa;   //alfa
    public static double H;  //wysokosc
    public static double W;  //szerokosc
    public static int nH;   //ilosc wezlow na wysokosc
    public static int nW;  //ilosc wezlow na szerokosc
    private int numberOfN;  //ilosc wszystkich wezłow
    private int numberOfE;  //ilosc wszystkich elementow
    private double heat;   //ciepło własciwe
    private double density; //gestosc (ro)
    private double k;  //conductivity (k) wsp.k
    private double ambientT;   //temperatura otoczenia
    private double stepTime;  //czas kroku symulacji
    private double simulationTime;  //czas symulacji

    public GlobalData()
    {
        /*alfa=300;
        H=0.1;
        W=0.1;
        nH=4;
        nW=4;
        heat=700;
        k=25;
        density=7800;
        numberOfN = nH*nW;
        numberOfE = (nH - 1)*(nW - 1);
        initialT = 100;
        ambientT = 1200;
        stepTime = 50;
        simulationTime = 500;*/

        alfa=300;
        H=0.1;
        W=0.1;
        nH=31;
        nW=31;
        heat=700;
        k=25;
        density=7800;
        numberOfN = nH*nW;
        numberOfE = (nH - 1)*(nW - 1);
        initialT = 100;
        ambientT = 1200;
        stepTime = 1;
        simulationTime = 100;
    }

    public double getAlfa() {
        return alfa;
    }

    public int getNumberOfN() {
        return numberOfN;
    }

    public int getNumberOfE() {
        return numberOfE;
    }

    public double getH() {
        return H;
    }

    public double getW() {
        return W;
    }

    public double getHeat() {
        return heat;
    }

    public double getK() {
        return k;
    }

    public double getDensity() {
        return density;
    }

    public int getnH() {
        return nH;
    }

    public int getnW() {
        return nW;
    }

    public double getInitialT(){
        return initialT;
    }

    public double getAmbientT() {return ambientT;}

    public double getStepTime() {return stepTime;}

    public double getSimulationTime() {return simulationTime;}
}
