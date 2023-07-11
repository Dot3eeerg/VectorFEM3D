namespace VectorFEM3D;

public class Grid
{
    private readonly double _xStart;
    private readonly double _xEnd;
    private readonly int _xSteps;
    private readonly double _xRaz;
    private readonly double _yStart;
    private readonly double _yEnd;
    private readonly int _ySteps;
    private readonly double _yRaz;
    private readonly double _zStart;
    private readonly double _zEnd;
    private readonly int _zSteps;
    private readonly double _zRaz;
    private readonly int[] _boundaries;
    public Point3D[] Nodes { get; private set; }
    public Edge3D[] Edges { get; private set; }
    public HashSet<int> DirichletBoundaries { get; private set; } 
    public List<(HashSet<(int, int)>, ElementSide)> NewmanBoundaries { get; private set; } 
    public int[][] Elements { get; private set; }
    public double Lambda { get; set; }
    public double Sigma { get; set; }

    public Grid(string path)
    {
        using (var sr = new StreamReader(path))
        {
            string[] data;
            data = sr.ReadLine()!.Split(" ").ToArray();
            _xStart = Convert.ToDouble(data[0]);
            _xEnd = Convert.ToDouble(data[1]);
            _xSteps = Convert.ToInt32(data[2]);
            _xRaz = Convert.ToDouble(data[3]);
            
            data = sr.ReadLine()!.Split(" ").ToArray();
            _yStart = Convert.ToDouble(data[0]);
            _yEnd = Convert.ToDouble(data[1]);
            _ySteps = Convert.ToInt32(data[2]);
            _yRaz = Convert.ToDouble(data[3]);
            
            data = sr.ReadLine()!.Split(" ").ToArray();
            _zStart = Convert.ToDouble(data[0]);
            _zEnd = Convert.ToDouble(data[1]);
            _zSteps = Convert.ToInt32(data[2]);
            _zRaz = Convert.ToDouble(data[3]);

            data = sr.ReadLine()!.Split(" ").ToArray();
            Lambda = Convert.ToDouble(data[0]);
            Sigma = Convert.ToDouble(data[1]);

            data = sr.ReadLine()!.Split(" ").ToArray();
            _boundaries = new int[6];
            _boundaries[0] = Convert.ToInt32(data[0]);
            _boundaries[1] = Convert.ToInt32(data[1]);
            _boundaries[2] = Convert.ToInt32(data[2]);
            _boundaries[3] = Convert.ToInt32(data[3]);
            _boundaries[4] = Convert.ToInt32(data[4]);
            _boundaries[5] = Convert.ToInt32(data[5]);
        }
    }

    public void BuildGrid()
    {
        Elements = new int[_xSteps * _ySteps * _zSteps].Select(_ => new int[12]).ToArray();
        Nodes = new Point3D[(_xSteps + 1) * (_ySteps + 1) * (_zSteps + 1)];

        double sumRazX = 0, sumRazY = 0, sumRazZ = 0;
        for (int i = 0; i < _xSteps; i++)
            sumRazX += Math.Pow(_xRaz, i);
        
        for (int i = 0; i < _ySteps; i++)
            sumRazY += Math.Pow(_yRaz, i);

        for (int i = 0; i < _zSteps; i++)
            sumRazZ += Math.Pow(_zRaz, i);

        int nodesInRow = _xSteps + 1;
        int nodesInSlice = nodesInRow * (_ySteps + 1);

        int zEdges = _zSteps * nodesInSlice;
        int xEdges = _xSteps;
        int yEdges = 1 + _xSteps;
        int edgesInSlice = xEdges * (1 + _ySteps) + yEdges * _ySteps;

        Edges = new Edge3D[edgesInSlice * (_zSteps + 1) + zEdges];

        double x = _xStart, y = _yStart, z = _zStart;
        double xStep = (_xEnd - _xStart) / sumRazX;
        double yStep = (_yEnd - _yStart) / sumRazY;
        double zStep = (_zEnd - _zStart) / sumRazZ;

        DirichletBoundaries = new();
        NewmanBoundaries = new();

        for (int j = 0; j < _xSteps; j++)
        {
            Nodes[j] = new(x, y, z);
            x += xStep;
            xStep *= _xRaz;
        }

        Nodes[_xSteps] = new(_xEnd, y, z);

        for (int i = 1; i <= _ySteps; i++)
        {
            y += yStep;
            yStep *= _yRaz;
            for (int j = 0; j < _xSteps + 1; j++)
            {
                Nodes[i * nodesInRow + j] = new(Nodes[j].X, y, z);
            }
        }

        for (int i = 1; i <= _zSteps; i++)
        {
            z += zStep;
            zStep *= _zRaz;
            for (int j = 0; j < _ySteps + 1; j++)
            {
                for (int k = 0; k < _xSteps + 1; k++)
                    Nodes[i * nodesInSlice + j * nodesInRow + k] = new(Nodes[k].X, Nodes[j * nodesInRow].Y, z);
            }
        }

        int index = 0;
        
        for (int j = 0; j < _zSteps + 1; j++)
        {
            int xLocal = 0;
            int yLocal = 0;
            int zLocal = 0;


            for (int i = 0; i < nodesInSlice / nodesInRow; i++)
            {
                for (int k = 0; k < nodesInRow - 1; k++)
                {
                    Edges[index++] = new Edge3D(Nodes[xLocal + nodesInSlice * j], Nodes[xLocal + nodesInSlice * j + 1]);
                    xLocal++;
                }

                xLocal++;

                if (i != nodesInSlice / nodesInRow - 1)
                {
                    for (int k = 0; k < nodesInRow; k++)
                    {
                        Edges[index++] = new Edge3D(Nodes[yLocal + nodesInSlice * j],
                            Nodes[yLocal + nodesInSlice * j + nodesInRow]);
                        yLocal++;
                    }
                }
            }

            if (j != _zSteps)
            {
                for (int k = 0; k < nodesInSlice; k++)
                {
                    Edges[index++] = new Edge3D(Nodes[zLocal + nodesInSlice * j],
                        Nodes[zLocal + nodesInSlice * j + nodesInSlice]);
                    zLocal++;
                }
            }
        }

        index = 0;
        
        for (int k = 0; k < _zSteps; k++)
        {
            for (int i = 0; i < _ySteps; i++)
            {
                for (int j = 0; j < _xSteps; j++)
                {
                    Elements[index][0] = j + (nodesInSlice + edgesInSlice) * k + (xEdges + yEdges) * i;
                    Elements[index][1] = j + (nodesInSlice + edgesInSlice) * k + (xEdges + yEdges) * (i + 1);
                    Elements[index][2] = j + (nodesInSlice + edgesInSlice) * k + (xEdges + yEdges) * i + xEdges;
                    Elements[index][3] = j + (nodesInSlice + edgesInSlice) * k + (xEdges + yEdges) * i + xEdges + 1;
                    
                    Elements[index][4] =
                        j + (nodesInSlice + edgesInSlice) * k + (xEdges + yEdges) * i + edgesInSlice - _xSteps * i;
                    Elements[index][5] =
                        j + (nodesInSlice + edgesInSlice) * k + (xEdges + yEdges) * i + edgesInSlice + 1 - _xSteps * i;
                    Elements[index][6] =
                        j + (nodesInSlice + edgesInSlice) * k + (xEdges + yEdges) * i + edgesInSlice + nodesInRow -
                        _xSteps * i;
                    Elements[index][7] =
                        j + (nodesInSlice + edgesInSlice) * k + (xEdges + yEdges) * i + edgesInSlice + 1 +
                        nodesInRow - _xSteps * i;

                    Elements[index][8] = j + (nodesInSlice + edgesInSlice) * (k + 1) + (xEdges + yEdges) * i;
                    Elements[index][9] = j + (nodesInSlice + edgesInSlice) * (k + 1) + (xEdges + yEdges) * (i + 1);
                    Elements[index][10] = j + (nodesInSlice + edgesInSlice) * (k + 1) + (xEdges + yEdges) * i + xEdges;
                    Elements[index++][11] = j + (nodesInSlice + edgesInSlice) * (k + 1) + (xEdges + yEdges) * i + xEdges + 1;
                }
            }
        }
    }

    public void AccountBoundaryConditions()
    {
        for (int ielem = 0; ielem < Elements.Length; ielem++)
        {
            if (ielem < _xSteps * _ySteps)
            {
                if (_boundaries[2] == 1) DirichletBoundary(ElementSide.Bottom, ielem);
                else if (_boundaries[2] == 2) NewmanBoundary(ElementSide.Bottom, ielem);
            }

            if (ielem >= _xSteps * _ySteps * _zSteps - _xSteps * _ySteps || _zSteps == 1)
            {
                if (_boundaries[3] == 1) DirichletBoundary(ElementSide.Upper, ielem);
                else if (_boundaries[3] == 2) NewmanBoundary(ElementSide.Upper, ielem);
            }

            if (ielem % _xSteps == 0)
            {
                if (_boundaries[0] == 1) DirichletBoundary(ElementSide.Left, ielem);
                else if (_boundaries[0] == 2) NewmanBoundary(ElementSide.Left, ielem);
            }

            if ((ielem + 1) % _xSteps == 0)
            {
                if (_boundaries[1] == 1) DirichletBoundary(ElementSide.Right, ielem);
                else if (_boundaries[1] == 2) NewmanBoundary(ElementSide.Right, ielem);
            }

            if (ielem % (_xSteps * _ySteps) < _xSteps)
            {
                if (_boundaries[5] == 1) DirichletBoundary(ElementSide.Front, ielem);
                else if (_boundaries[5] == 2) NewmanBoundary(ElementSide.Front, ielem);
            }

            if (ielem % (_xSteps * _ySteps) >= _xSteps * _ySteps - _xSteps)
            {
                if (_boundaries[4] == 1) DirichletBoundary(ElementSide.Rear, ielem);
                else if (_boundaries[4] == 2) NewmanBoundary(ElementSide.Rear, ielem);
            }
        }
    }

    private void DirichletBoundary(ElementSide elementSide, int ielem)
    {
        switch (elementSide)
        {
            case ElementSide.Bottom:
                DirichletBoundaries.Add(Elements[ielem][0]);
                DirichletBoundaries.Add(Elements[ielem][1]);
                DirichletBoundaries.Add(Elements[ielem][2]);
                DirichletBoundaries.Add(Elements[ielem][3]);
                break;
            
            case ElementSide.Upper:
                DirichletBoundaries.Add(Elements[ielem][8]);
                DirichletBoundaries.Add(Elements[ielem][9]);
                DirichletBoundaries.Add(Elements[ielem][10]);
                DirichletBoundaries.Add(Elements[ielem][11]);
                break;
            
            case ElementSide.Left:
                DirichletBoundaries.Add(Elements[ielem][2]);
                DirichletBoundaries.Add(Elements[ielem][4]);
                DirichletBoundaries.Add(Elements[ielem][5]);
                DirichletBoundaries.Add(Elements[ielem][10]);
                break;
                
            case ElementSide.Right:
                DirichletBoundaries.Add(Elements[ielem][3]);
                DirichletBoundaries.Add(Elements[ielem][5]);
                DirichletBoundaries.Add(Elements[ielem][7]);
                DirichletBoundaries.Add(Elements[ielem][11]);
                break;
            
            case ElementSide.Front:
                DirichletBoundaries.Add(Elements[ielem][0]);
                DirichletBoundaries.Add(Elements[ielem][4]);
                DirichletBoundaries.Add(Elements[ielem][5]);
                DirichletBoundaries.Add(Elements[ielem][8]);
                break;
            
            case ElementSide.Rear:
                DirichletBoundaries.Add(Elements[ielem][1]);
                DirichletBoundaries.Add(Elements[ielem][6]);
                DirichletBoundaries.Add(Elements[ielem][7]);
                DirichletBoundaries.Add(Elements[ielem][9]);
                break;
        }
    }

    private void NewmanBoundary(ElementSide elementSide, int ielem)
    {
        NewmanBoundaries.Add((new(), elementSide));

        switch (elementSide)
        {
            case ElementSide.Bottom:
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][0], 0));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][1], 1));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][2], 2));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][3], 3));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][4], 4));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][5], 5));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][6], 6));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][7], 7));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][8], 8));
                break;
            
            case ElementSide.Upper:
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][18], 18));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][19], 19));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][20], 20));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][21], 21));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][22], 22));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][23], 23));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][24], 24));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][25], 25));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][26], 26));
                break;
            
            case ElementSide.Left:
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][0], 0));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][3], 3));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][6], 6));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][9], 9));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][12], 12));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][15], 15));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][18], 18));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][21], 21));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][24], 24));
                break;
            
            case ElementSide.Right:
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][2], 2));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][5], 5));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][8], 8));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][11], 11));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][14], 14));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][17], 17));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][20], 20));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][23], 23));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][26], 26));
                break;
            
            case ElementSide.Front:
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][0], 0));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][1], 1));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][2], 2));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][9], 9));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][10], 10));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][11], 11));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][18], 18));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][19], 19));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][20], 20));
                break;
            
            case ElementSide.Rear:
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][6], 6));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][7], 7));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][8], 8));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][15], 15));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][16], 16));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][17], 17));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][24], 24));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][25], 25));
                NewmanBoundaries[^1].Item1.Add((Elements[ielem][26], 26)); 
                break;
        }
    }
}

public class TimeGrid
{
    private readonly double _tStart;
    private readonly double _tEnd;
    private readonly int _tSteps;
    private readonly double _tRaz;
    public double[] TGrid { get; set;  }

    public TimeGrid(string path)
    {
        using (var sr = new StreamReader(path))
        {
            string[] data;
            data = sr.ReadLine()!.Split(" ").ToArray();
            _tStart = Convert.ToDouble(data[0]);
            _tEnd = Convert.ToDouble(data[1]);
            _tSteps = Convert.ToInt32(data[2]);
            _tRaz = Convert.ToDouble(data[3]);
            TGrid = new double[_tSteps + 1];
        }
    }

    public double this[int index]
    {
        get => TGrid[index];
        set => TGrid[index] = value;
    }

    public void BuildTimeGrid()
    {
        double sumRaz = 0;
        for (int i = 0; i < _tSteps; i++)
            sumRaz += Math.Pow(_tRaz, i);

        double t = _tStart;
        double tStep = (_tEnd - _tStart) / sumRaz;

        for (int i = 0; i < _tSteps; i++)
        {
            TGrid[i] = t;
            t += tStep;
            tStep *= _tRaz;
        }

        TGrid[_tSteps] = _tEnd;
    }
}