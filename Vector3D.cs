namespace VectorFEM3D;

public class Vector3D
{
    private readonly double _x;
    private readonly double _y;
    private readonly double _z;

    public Vector3D(double x, double y, double z)
    {
        _x = x;
        _y = y;
        _z = z;
    }

    public static double operator *(Vector3D vector1, Vector3D vector2)
        => vector1._x * vector2._x + vector1._y * vector2._y + vector1._z * vector2._z;
}