public class Grid
{
    private Node[] nodes;
    private Element[] elements;
    private GlobalData globalData;

    public Grid()
    {
        globalData =new GlobalData();
        nodes =new Node[globalData.getNumberOfN()];
        elements = new Element[globalData.getNumberOfE()];
    }

    public void gridGeneration()
    {
        double dH = globalData.getH() / (globalData.getnH() - 1);
        double dW = globalData.getW() / (globalData.getnW() - 1);

        for (int i = 0; i < globalData.getnW(); i ++){
            for (int j = 0; j < globalData.getnH(); j ++)
            {
                int n = i * globalData.getnW() + j;
                boolean BC = false;
                if (i==0 || i==globalData.getnW()-1 || j==0 || j==globalData.getnH()-1)
                {
                    BC = true;
                }
                nodes[n] = new Node(i*dW, j*dH, globalData.getInitialT(), BC, n, this);
            }
        }

        for (int i = 0; i < globalData.getnW() - 1; i++){
            for (int j = 0; j < globalData.getnH() - 1; j ++)
            {
                int n  = i * (globalData.getnW() - 1) + j;

                elements[n] = new Element(this, n);

                elements[n].setNode(0,n+i);
                elements[n].setNode(1, n + i + globalData.getnH());
                elements[n].setNode(2, n + i + globalData.getnH() + 1);
                elements[n].setNode(3, n + i + 1);
            }
        }
    }
    public void showGird(){
        for (int i = 0; i < globalData.getNumberOfE(); i++)
            elements[i].show();
    }

    public Element[] getElements(){
        return elements;
    }

    public Node[] getNodes(){
        return nodes;
    }

}

