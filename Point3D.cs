namespace VectorFEM3D;

public class Point3D
{
    public double X { get; set; }
    public double Y { get; set; }
    public double Z { get; set; }

    public Point3D(double x, double y, double z)
    {
        X = x;
        Y = y;
        Z = z;
    }

    public static Point3D operator +(Point3D point, (double, double, double) value)
        => new(point.X + value.Item1, point.Y + value.Item2, point.Z + value.Item3);

    public static Point3D operator -(Point3D point, (double, double, double) value)
        => new(point.X - value.Item1, point.Y - value.Item2, point.Z - value.Item3);
        
    public static Point3D operator *(Point3D point, (double, double, double) value)
        => new(point.X * value.Item1, point.Y * value.Item2, point.Z * value.Item3);

    public static Point3D operator /(Point3D point, (double, double, double) value)
        => new(point.X / value.Item1, point.Y / value.Item2, point.Z / value.Item3);
}

public class Point2D
{
    public double X { get; set; }
    public double Y { get; set; }

    public Point2D(double x, double y)
    {
        X = x;
        Y = y;
    }

    public static Point2D operator +(Point2D point, (double, double) value)
        => new(point.X + value.Item1, point.Y + value.Item2);

    public static Point2D operator -(Point2D point, (double, double) value)
        => new(point.X - value.Item1, point.Y - value.Item2);
        
    public static Point2D operator *(Point2D point, (double, double) value)
        => new(point.X * value.Item1, point.Y * value.Item2);

    public static Point2D operator /(Point2D point, (double, double) value)
        => new(point.X / value.Item1, point.Y / value.Item2);
}
