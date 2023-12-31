﻿using System.Data;

namespace VectorFEM3D;

public class Vector3D
{
    private double _x;
    private double _y;
    private double _z;

    public Vector3D(double x, double y, double z)
    {
        _x = x;
        _y = y;
        _z = z;
    }

    public Vector3D UpdateVector(double x, double y, double z)
    {
        _x = x;
        _y = y;
        _z = z;
        return this;
    }

    public Vector3D Copy(Vector3D source)
    {
        _x = source._x;
        _y = source._y;
        _z = source._z;
        return this;
    }

    public static double operator *(Vector3D vector1, Vector3D vector2)
        => vector1._x * vector2._x + vector1._y * vector2._y + vector1._z * vector2._z;

    public static double DotProductJacob(Vector3D vector1, Vector3D vector2, double hx, double hy, double hz)
        => vector1._x * vector2._x * hy * hz / hx + vector1._y * vector2._y * hx * hz / hy +
           vector1._z * vector2._z * hx * hy / hz;
}