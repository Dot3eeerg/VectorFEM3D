using System.Runtime.CompilerServices;

namespace VectorFEM3D;

public abstract class Test
{
    protected double lambda;
    protected double sigma;

    public Test(Grid grid)
    {
        lambda = grid.Lambda;
        sigma = grid.Sigma;
    }
    
    public abstract double U(Point3D point, double t, int i);
    
    public abstract double F(Point3D point, double t, int i);

    public abstract double UValue(Point3D point, int i);

    public abstract double Theta(Point3D point, double t, ElementSide elementSide);
}

public class Test1 : Test
{
    public Test1(Grid grid) : base(grid) { }
    
    public override double U(Point3D point, double t, int i)
        => i switch
        {
            0 => UValue(point, 0),
            1 => UValue(point, 0),
            2 => UValue(point, 1),
            3 => UValue(point, 1),
            4 => UValue(point, 2),
            5 => UValue(point, 2),
            6 => UValue(point, 2),
            7 => UValue(point, 2),
            8 => UValue(point, 0),
            9 => UValue(point, 0),
            10 => UValue(point, 1),
            11 => UValue(point, 1),
        };

    public override double F(Point3D point, double t, int i)
        => i switch
        {
            0 => FValue(point, 0),
            1 => FValue(point, 0),
            2 => FValue(point, 1),
            3 => FValue(point, 1),
            4 => FValue(point, 2),
            5 => FValue(point, 2),
            6 => FValue(point, 2),
            7 => FValue(point, 2),
            8 => FValue(point, 0),
            9 => FValue(point, 0),
            10 => FValue(point, 1),
            11 => FValue(point, 1),
        };
    
    public override double UValue(Point3D point, int i)
        => i switch
        {
            0 => 1 + point.Y,
            1 => 0,
            2 => 0,
        };

    private double FValue(Point3D point, int i)
        => i switch
        {
            0 => 1 + point.Y,
            1 => 0,
            2 => 0,
        };

    public override double Theta(Point3D point, double t, ElementSide elementSide)
    {
        switch (elementSide)
        {
            case ElementSide.Left:
            case ElementSide.Right:
                return 0;
            
            case ElementSide.Bottom:
            case ElementSide.Upper:
                return 2;
            
            case ElementSide.Rear:
            case ElementSide.Front:
                return 1;
            
            default:
                throw new Exception("Can't find function type");
        }
    }
}