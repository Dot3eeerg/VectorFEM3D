namespace VectorFEM3D;

public interface IBasis3D
{
    int Size { get; }
    Vector3D GetPsi(int number, Point3D point);
    Vector3D GetDPsi(int number, Point3D point);
}

public readonly record struct TriLinearVectorBasis : IBasis3D
{
    public int Size => 12;
    private readonly Vector3D _vector = new Vector3D(0, 0, 0);

    public TriLinearVectorBasis() { }
    
    public Vector3D GetPsi(int number, Point3D point)
        => number switch
        {
            0 => _vector.UpdateVector(GetXi(0, point.Y) * GetXi(0, point.Z), 0, 0),
            1 => _vector.UpdateVector(GetXi(1, point.Y) * GetXi(0, point.Z), 0, 0),
            4 => _vector.UpdateVector(0, GetXi(0, point.X) * GetXi(0, point.Z), 0),
            5 => _vector.UpdateVector(0, GetXi(1, point.X) * GetXi(0, point.Z), 0),
            8 => _vector.UpdateVector(0, 0, GetXi(0, point.X) * GetXi(0, point.Y)),
            9 => _vector.UpdateVector(0, 0, GetXi(1, point.X) * GetXi(0, point.Y)),
            10 => _vector.UpdateVector(0, 0, GetXi(0, point.X) * GetXi(1, point.Y)),
            11 => _vector.UpdateVector(0, 0, GetXi(1, point.X) * GetXi(1, point.Y)),
            2 => _vector.UpdateVector(GetXi(0, point.Y) * GetXi(1, point.Z), 0, 0),
            3 => _vector.UpdateVector(GetXi(1, point.Y) * GetXi(1, point.Z), 0, 0),
            6 => _vector.UpdateVector(0, GetXi(0, point.X) * GetXi(1, point.Z), 0),
            7 => _vector.UpdateVector(0, GetXi(1, point.X) * GetXi(1, point.Z), 0),
            _ => throw new ArgumentOutOfRangeException(nameof(number), number, "Not expected function number")
        };
    
    public Vector3D GetDPsi(int number, Point3D point)
        => number switch
        {
            0 => _vector.UpdateVector(0, -GetXi(0, point.Y), GetXi(0, point.Z)),
            1 => _vector.UpdateVector(0, -GetXi(1, point.Y), -GetXi(0, point.Z)),
            4 => _vector.UpdateVector(GetXi(0, point.X), 0, -GetXi(0, point.Z)),
            5 => _vector.UpdateVector(GetXi(1, point.X), 0, GetXi(0, point.Z)),
            8 => _vector.UpdateVector(-GetXi(0, point.X), GetXi(0, point.Y), 0),
            9 => _vector.UpdateVector(-GetXi(1, point.X), -GetXi(0, point.Y), 0),
            10 => _vector.UpdateVector(GetXi(0, point.X), GetXi(1, point.Y), 0),
            11 => _vector.UpdateVector(GetXi(1, point.X), -GetXi(1, point.Y), 0),
            2 => _vector.UpdateVector(0,  GetXi(0, point.Y), GetXi(1, point.Z)),
            3 => _vector.UpdateVector(0, GetXi(1, point.Y), -GetXi(1, point.Z)),
            6 => _vector.UpdateVector(-GetXi(0, point.X), 0, -GetXi(1, point.Z)),
            7 => _vector.UpdateVector(-GetXi(1, point.X), 0, GetXi(1, point.Z)),
            _ => throw new ArgumentOutOfRangeException(nameof(number), number, "Not expected function number")
        };

    private double GetXi(int number, double value)
        => number switch
        {
            0 => 1 - value,
            1 => value,
            _ => throw new ArgumentOutOfRangeException(nameof(number), number, "Not expected Xi member")
        };
}

public interface IBasis2D
{
    int Size { get; }
    double GetPsi(int number, Point2D point);
    double GetDPsi(int number, int dNumber, Point2D point);
}

public readonly record struct BiLinearBasis : IBasis2D
{
    public int Size => 4;

    public BiLinearBasis() { }
    
    public double GetPsi(int number, Point2D point)
        => number switch
        {
            0 => GetXi(0, point.X) * GetXi(0, point.Y),
            1 => GetXi(1, point.X) * GetXi(0, point.Y),
            2 => GetXi(0, point.X) * GetXi(1, point.Y),
            3 => GetXi(1, point.X) * GetXi(1, point.Y),
            _ => throw new ArgumentOutOfRangeException(nameof(number), number, "Not expected function number")
        };
    
    public double GetDPsi(int number, int dNumber, Point2D point)
        => dNumber switch
        {
            0 => number switch
            {
                0 => -GetXi(0, point.Y),
                1 => GetXi(0, point.Y),
                2 => -GetXi(1, point.Y),
                3 => GetXi(1, point.Y)
            },
            1 => number switch
            {
                0 => -GetXi(0, point.X),
                1 => -GetXi(1, point.X),
                2 => GetXi(0, point.X),
                3 => GetXi(1, point.X)
            },
            _ => throw new ArgumentOutOfRangeException(nameof(number), number, "Not expected function number")
        };

    private double GetXi(int number, double value)
        => number switch
        {
            0 => 1 - value,
            1 => value,
            _ => throw new ArgumentOutOfRangeException(nameof(number), number, "Not expected Xi member")
        };
}
