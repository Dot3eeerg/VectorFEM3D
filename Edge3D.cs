namespace VectorFEM3D;

public class Edge3D
{
    public Point3D Point0 { get; set; }
    public Point3D Point1 { get; set; }
    public double Length { get; set; }
    public Point3D Point { get; set; }

    public Edge3D(Point3D point0, Point3D point1)
    {
        Point0 = point0;
        Point1 = point1;
        Length = Math.Sqrt(Math.Pow(Point1.X - Point0.X, 2) + Math.Pow(point1.Y - point0.Y, 2) +
                           Math.Pow(point1.Z - point0.Z, 2));

        Point = new((Point1.X + Point0.X) / 2, (Point1.Y + Point0.Y) / 2, (Point1.Z + Point0.Z) / 2);
    }

    public int GetAxis()
    {
        if (Point.X - Point0.X > 1e-12)
            return 0;
        
        if (Point.Y - Point0.Y > 1e-12)
            return 1;

        return 2;
    }
}