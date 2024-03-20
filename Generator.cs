namespace VectorFEM3D;

public class Generator
{
    public double xStart;
    public double yStart;
    public double zStart;
    public double xEnd;
    public double yEnd;
    public double zEnd;

    public double Length => xEnd - xStart;
    public double Width => yEnd - yStart;
    public double Height => zEnd - zStart;

    public Generator(double xs, double ys, double zs, double xe, double ye, double ze)
    {
        xStart = xs;
        yStart = ys;
        zStart = zs;
        xEnd = xe;
        yEnd = ye;
        zEnd = ze;
    }
}