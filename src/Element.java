public class Element
{
    public double[][] localH;
    public double[][] localC;
    public double[][] Hbc;
    private int ownID;
    private int ID[];
    private boolean BC[];

    public Element(Grid grid, int ownID)
    {
        ID = new int[4];
        BC = new boolean[4];
        this.ownID=ownID;
    }

    public int[] getID() {return ID;}

    public int getElementID() {return ownID;}

    public void setNode(int position, int newID){
        this.ID[position]=newID;
    }

    public void setLocalH(double[][] localH){
        this.localH= localH;
    }

    public void setLocalC(double[][] localC){
        this.localC= localC;
    }

    public void printLocalH(){
        System.out.println("Macierz H dla elementu " + getElementID());
        for (int j = 0; j < 4; j++){
            for (int k = 0; k < 4; k++){
                System.out.print(localH[j][k] + " ");
            }
            System.out.println(" ");
        }
    }

    public void printLocalC(){
        System.out.println("Macierz C dla elementu " + getElementID());
        for (int j = 0; j < 4; j++){
            for (int k = 0; k < 4; k++){
                System.out.print(localC[j][k] + " ");
            }
            System.out.println(" ");
        }
    }

    public void printHbc(){
        System.out.println("Macierz Hbc dla elementu " + getElementID());
        for (int j = 0; j < 4; j++){
            for (int k = 0; k < 4; k++){
                System.out.print(Hbc[j][k] + " ");
            }
            System.out.println(" ");
        }
    }

    public void show(){
        System.out.println("Element "+ownID+":");
        System.out.print(ID[3]);
        System.out.print(BC[3] ? "----" : "    ");
        System.out.print(ID[2]+"\n");
        System.out.print(BC[0] ? "|" : " ");
        System.out.print("     ");
        System.out.print(BC[2] ? "|\n" : " \n");
        System.out.print(ID[0]);
        System.out.print(BC[1] ? "----" : "    ");
        System.out.print(ID[1]+"\n\n");
    }
}

