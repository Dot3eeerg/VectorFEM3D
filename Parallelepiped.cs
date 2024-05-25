namespace VectorFEM3D;

public class Parallelepiped
{
    public Point3D[] PointList;
    public double Sigma;

    public double hx;
    public double hy;
    public double hz;

    public void SetElement(Point3D[] pointList, double x, double y, double z, double sigma)
    {
        PointList = pointList;

        hx = x;
        hy = y;
        hz = z;

        Sigma = sigma;
    }
}