namespace VectorFEM3D;

public class Edge3D
{
    public double X0 { get; set; }
    public double Y0 { get; set; }
    public double Z0 { get; set; }
    public double X1 { get; set; }
    public double Y1 { get; set; }
    public double Z1 { get; set; }
    public double Length { get; set; }

    public Edge3D(Point3D point0, Point3D point1)
    {
        X0 = point0.X;
        Y0 = point0.Y;
        Z0 = point0.Z;
        X1 = point1.X;
        Y1 = point1.Y;
        Z1 = point1.Z;
        Length = Math.Sqrt(Math.Pow(X1 - X0, 2) + Math.Pow(Y1 - Y0, 2) + Math.Pow(Z1 - Z0, 2));
    }
}