﻿namespace VectorFEM3D;

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
            2 => _vector.UpdateVector(0, GetXi(0, point.X) * GetXi(0, point.Z), 0),
            3 => _vector.UpdateVector(0, GetXi(1, point.X) * GetXi(0, point.Z), 0),
            4 => _vector.UpdateVector(0, 0, GetXi(0, point.X) * GetXi(0, point.Y)),
            5 => _vector.UpdateVector(0, 0, GetXi(1, point.X) * GetXi(0, point.Y)),
            6 => _vector.UpdateVector(0, 0, GetXi(0, point.X) * GetXi(1, point.Y)),
            7 => _vector.UpdateVector(0, 0, GetXi(1, point.X) * GetXi(1, point.Y)),
            8 => _vector.UpdateVector(GetXi(0, point.Y) * GetXi(1, point.Z), 0, 0),
            9 => _vector.UpdateVector(GetXi(1, point.Y) * GetXi(1, point.Z), 0, 0),
            10 => _vector.UpdateVector(0, GetXi(0, point.X) * GetXi(1, point.Z), 0),
            11 => _vector.UpdateVector(0, GetXi(1, point.X) * GetXi(1, point.Z), 0),
            _ => throw new ArgumentOutOfRangeException(nameof(number), number, "Not expected function number")
        };
    
    public Vector3D GetDPsi(int number, Point3D point)
        => number switch
        {
            0 => _vector.UpdateVector(0, -GetXi(0, point.Y), GetXi(0, point.Z)),
            1 => _vector.UpdateVector(0, -GetXi(1, point.Y), -GetXi(0, point.Z)),
            2 => _vector.UpdateVector(GetXi(0, point.X), 0, -GetXi(0, point.Z)),
            3 => _vector.UpdateVector(GetXi(1, point.X), 0, GetXi(0, point.Z)),
            4 => _vector.UpdateVector(-GetXi(0, point.X), GetXi(0, point.Y), 0),
            5 => _vector.UpdateVector(-GetXi(1, point.X), -GetXi(0, point.Y), 0),
            6 => _vector.UpdateVector(GetXi(0, point.X), GetXi(1, point.Y), 0),
            7 => _vector.UpdateVector(GetXi(1, point.X), -GetXi(1, point.Y), 0),
            8 => _vector.UpdateVector(0,  GetXi(0, point.Y), GetXi(1, point.Z)),
            9 => _vector.UpdateVector(0, GetXi(1, point.Y), -GetXi(1, point.Z)),
            10 => _vector.UpdateVector(-GetXi(0, point.X), 0, -GetXi(1, point.Z)),
            11 => _vector.UpdateVector(-GetXi(1, point.X), 0, GetXi(1, point.Z)),
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