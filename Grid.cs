namespace VectorFEM3D;

public class Grid
{
    private readonly List<double> _xStart = new List<double>();
    private readonly List<double> _xEnd = new List<double>();
    private readonly List<int> _xSteps = new List<int>();
    private readonly List<double> _xRaz = new List<double>();
    private readonly List<double> _yStart = new List<double>();
    private readonly List<double> _yEnd = new List<double>();
    private readonly List<int> _ySteps = new List<int>();
    private readonly List<double> _yRaz = new List<double>();
    private readonly List<double> _zStart = new List<double>();
    private readonly List<double> _zEnd = new List<double>();
    private readonly List<int> _zSteps = new List<int>();
    private readonly List<double> _zRaz = new List<double>();
    private readonly int[] _boundaries;

    private readonly double[] _xZones;
    private readonly double[] _yZones;
    private readonly double[] _zZones;
    private readonly int[][] _zones;
    private readonly double[] _sigmaValues;

    private readonly List<double> _xValues = new List<double>();
    private readonly List<double> _yValues = new List<double>();
    private readonly List<double> _zValues = new List<double>();

    public (double, double) GeneratorX;
    public (double, double) GeneratorY;
    public double GeneratorZ;
    public double SourceValue;
    public (double, double) Receiver;

    private List<int> _sumSteps;

    private List<List<Point3D>> NodeList = new List<List<Point3D>>();
    private Edge3D[][] EdgeList;
    private List<(int, HashSet<int>)> _identicalPoints = new List<(int, HashSet<int>)>();
    public Point3D[] Nodes { get; private set; }
    public Edge3D[] Edges { get; private set; }
    public HashSet<int> DirichletBoundaries { get; private set; } 
    public List<(HashSet<(int, int)>, ElementSide)> NewmanBoundaries { get; private set; } 
    public int[][] Elements { get; private set; }
    public double Mu { get; set; }
    public double Sigma { get; set; }
    public double Epsilon { get; set; }

    public Grid(string path)
    {
        using (var sr = new StreamReader(path))
        {
            string[] data;
            data = sr.ReadLine()!.Split(" ").ToArray();
            int kek = Convert.ToInt32(data[0]);
            for (int i = 0; i < kek; i++)
            {
                data = sr.ReadLine()!.Split(" ").ToArray();
                _xStart.Add(Convert.ToDouble(data[0]));
                _xEnd.Add(Convert.ToDouble(data[1]));
                _xSteps.Add(Convert.ToInt32(data[2]));
                _xRaz.Add(Convert.ToDouble(data[3]));
            }

            _xZones = sr.ReadLine()!.Split(" ").Select(x => Convert.ToDouble(x)).ToArray();
            
            data = sr.ReadLine()!.Split(" ").ToArray();
            kek = Convert.ToInt32(data[0]);
            for (int i = 0; i < kek; i++)
            {
                data = sr.ReadLine()!.Split(" ").ToArray();
                _yStart.Add(Convert.ToDouble(data[0]));
                _yEnd.Add(Convert.ToDouble(data[1]));
                _ySteps.Add(Convert.ToInt32(data[2]));
                _yRaz.Add(Convert.ToDouble(data[3]));
            }
            
            _yZones = sr.ReadLine()!.Split(" ").Select(x => Convert.ToDouble(x)).ToArray();
            
            data = sr.ReadLine()!.Split(" ").ToArray();
            kek = Convert.ToInt32(data[0]);
            for (int i = 0; i < kek; i++)
            {
                data = sr.ReadLine()!.Split(" ").ToArray();
                _zStart.Add(Convert.ToDouble(data[0]));
                _zEnd.Add(Convert.ToDouble(data[1]));
                _zSteps.Add(Convert.ToInt32(data[2]));
                _zRaz.Add(Convert.ToDouble(data[3]));
            }
            
            _zZones = sr.ReadLine()!.Split(" ").Select(x => Convert.ToDouble(x)).ToArray();

            data = sr.ReadLine()!.Split(" ").ToArray();
            _boundaries = new int[6];
            _boundaries[0] = Convert.ToInt32(data[0]);
            _boundaries[1] = Convert.ToInt32(data[1]);
            _boundaries[2] = Convert.ToInt32(data[2]);
            _boundaries[3] = Convert.ToInt32(data[3]);
            _boundaries[4] = Convert.ToInt32(data[4]);
            _boundaries[5] = Convert.ToInt32(data[5]);
            
            data = sr.ReadLine()!.Split(" ").ToArray();
            Mu = Convert.ToDouble(data[0]);
            //Sigma = Convert.ToDouble(data[1]);
            Epsilon = Convert.ToDouble(data[1]);

            kek = Convert.ToInt32(sr.ReadLine());
            _zones = new int[kek].Select(_ => new int[6]).ToArray();
            _sigmaValues = new double[kek];

            for (int i = 0; i < kek; i++)
            {
                data = sr.ReadLine()!.Split(" ").ToArray();
                _sigmaValues[i] = Convert.ToDouble(data[0]);
                _zones[i][0] = Convert.ToInt32(data[1]);
                _zones[i][1] = Convert.ToInt32(data[2]);
                _zones[i][2] = Convert.ToInt32(data[3]);
                _zones[i][3] = Convert.ToInt32(data[4]);
                _zones[i][4] = Convert.ToInt32(data[5]);
                _zones[i][5] = Convert.ToInt32(data[6]);
            }

            data = sr.ReadLine()!.Split(" ").ToArray();
            GeneratorX.Item1 = Convert.ToDouble(data[0]);
            GeneratorX.Item2 = Convert.ToDouble(data[1]);
            GeneratorY.Item1 = Convert.ToDouble(data[2]);
            GeneratorY.Item2 = Convert.ToDouble(data[3]);
            GeneratorZ = Convert.ToDouble(data[4]);
            SourceValue = Convert.ToDouble(data[5]);

            data = sr.ReadLine()!.Split(" ").ToArray();
            Receiver.Item1 = Convert.ToDouble(data[0]);
            Receiver.Item2 = Convert.ToDouble(data[1]);
        }
    }

    public void BuildGrid()
    {
        _sumSteps = new List<int>();
        int count = 1;
        int nodeCount = 1;
        _sumSteps.Add(0);
        for (int i = 0; i < _xSteps.Count; i++)
        {
            _sumSteps[0] += _xSteps[i];
        }
        count *= _sumSteps[0];
        nodeCount *= _sumSteps[0] + 1;
        int nodesInRow = _sumSteps[0] + 1;
        
        _sumSteps.Add(0);
        for (int i = 0; i < _ySteps.Count; i++)
        {
            _sumSteps[1] += _ySteps[i];
        }
        count *= _sumSteps[1];
        nodeCount *= _sumSteps[1] + 1;
        int nodesInSlice = nodesInRow * (_sumSteps[1] + 1);
        
        _sumSteps.Add(0);
        for (int i = 0; i < _zSteps.Count; i++)
        {
            _sumSteps[2] += _zSteps[i];
        }
        count *= _sumSteps[2];
        nodeCount *= _sumSteps[2] + 1;
        
        Elements = new int[count].Select(_ => new int[12]).ToArray();
        Nodes = new Point3D[nodeCount];
        //for (int i = 0; i < 8; i++)
        //{
        //    NodeList.Add(new List<Point3D>());
        //}

        List<double> sumRazX = new(), sumRazY = new(), sumRazZ = new();
        for (int i = 0; i < _xSteps.Count; i++)
        {
            sumRazX.Add(0);
            for (int j = 0; j < _xSteps[i]; j++)
            {
                sumRazX[i] += Math.Pow(_xRaz[i], j);
            }
        }
        
        for (int i = 0; i < _ySteps.Count; i++)
        {
            sumRazY.Add(0);
            for (int j = 0; j < _ySteps[i]; j++)
            {
                sumRazY[i] += Math.Pow(_yRaz[i], j);
            }
        }
        
        for (int i = 0; i < _zSteps.Count; i++)
        {
            sumRazZ.Add(0);
            for (int j = 0; j < _zSteps[i]; j++)
            {
                sumRazZ[i] += Math.Pow(_zRaz[i], j);
            }
        }
        
        int xEdges = _sumSteps[0];
        int yEdges = 1 + _sumSteps[0];
        int zEdges = _sumSteps[2] * nodesInSlice;
        int edgesInSlice = xEdges * (1 + _sumSteps[1]) + yEdges * _sumSteps[1];

        Edges = new Edge3D[edgesInSlice * (_sumSteps[2] + 1) + zEdges];
        //EdgeList = new Edge3D[8].Select(_ => new Edge3D[edgesInSlice * (_zSteps + 1) + zEdges]).ToArray();
        
        DirichletBoundaries = new();

        double x = 0, y = 0, z = 0;
        double xStep = 0, yStep = 0, zStep = 0;
        //double xStep = (_xEnd - _xStart) / sumRazX;
        //double yStep = (_yEnd - _yStart) / sumRazY;
        //double zStep = (_zEnd - _zStart) / sumRazZ;

        for (int i = 0; i < _xSteps.Count; i++)
        {
            x = _xStart[i];
            xStep = (_xEnd[i] - _xStart[i]) / sumRazX[i];
            
            for (int j = 0; j < _xSteps[i]; j++)
            {
                _xValues.Add(x);
                x += xStep;
                xStep *= _xRaz[i];
            }
            
        }
        _xValues.Add(_xEnd[^1]);
        
        for (int i = 0; i < _ySteps.Count; i++)
        {
            y = _yStart[i];
            yStep = (_yEnd[i] - _yStart[i]) / sumRazY[i];
            
            for (int j = 0; j < _ySteps[i]; j++)
            {
                _yValues.Add(y);
                y += yStep;
                yStep *= _yRaz[i];
            }
            
        }
        _yValues.Add(_yEnd[^1]);
        
        for (int i = 0; i < _zSteps.Count; i++)
        {
            z = _zStart[i];
            zStep = (_zEnd[i] - _zStart[i]) / sumRazZ[i];
            
            for (int j = 0; j < _zSteps[i]; j++)
            {
                _zValues.Add(z);
                z += zStep;
                zStep *= _zRaz[i];
            }
            
        }
        _zValues.Add(_zEnd[^1]);

        for (int i = 0; i < _zValues.Count; i++)
        {
            for (int j = 0; j < _yValues.Count; j++)
            {
                for (int k = 0; k < _xValues.Count; k++)
                {
                    Nodes[i * nodesInSlice + j * nodesInRow + k] = new Point3D(_xValues[k], _yValues[j], _zValues[i]);
                }
            }
        }

        int index = 0;
        
        for (int j = 0; j < _sumSteps[2] + 1; j++)
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

            if (j != _sumSteps[2])
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

        for (int k = 0; k < _sumSteps[2]; k++)
        {
            for (int i = 0; i < _sumSteps[1]; i++)
            {
                for (int j = 0; j < _sumSteps[0]; j++)
                {
                    // Bottom part
                    Elements[index][0] = j + (nodesInSlice + edgesInSlice) * k + (xEdges + yEdges) * i;
                    Elements[index][1] = j + (nodesInSlice + edgesInSlice) * k + (xEdges + yEdges) * (i + 1);
                    Elements[index][4] = j + (nodesInSlice + edgesInSlice) * k + (xEdges + yEdges) * i + xEdges;
                    Elements[index][5] = j + (nodesInSlice + edgesInSlice) * k + (xEdges + yEdges) * i + xEdges + 1;
                    
                    // Facet part
                    Elements[index][8] =
                        j + (nodesInSlice + edgesInSlice) * k + (xEdges + yEdges) * i + edgesInSlice - _sumSteps[0] * i;
                    Elements[index][9] =
                        j + (nodesInSlice + edgesInSlice) * k + (xEdges + yEdges) * i + edgesInSlice + 1 - _sumSteps[0] * i;
                    Elements[index][10] =
                        j + (nodesInSlice + edgesInSlice) * k + (xEdges + yEdges) * i + edgesInSlice + nodesInRow -
                        _sumSteps[0] * i;
                    Elements[index][11] =
                        j + (nodesInSlice + edgesInSlice) * k + (xEdges + yEdges) * i + edgesInSlice + 1 +
                        nodesInRow - _sumSteps[0] * i;

                    // Upper part
                    Elements[index][2] = j + (nodesInSlice + edgesInSlice) * (k + 1) + (xEdges + yEdges) * i;
                    Elements[index][3] = j + (nodesInSlice + edgesInSlice) * (k + 1) + (xEdges + yEdges) * (i + 1);
                    Elements[index][6] = j + (nodesInSlice + edgesInSlice) * (k + 1) + (xEdges + yEdges) * i + xEdges;
                    Elements[index++][7] = j + (nodesInSlice + edgesInSlice) * (k + 1) + (xEdges + yEdges) * i + xEdges + 1;
                }
            }
        }
    }

    public double GetSigma(Point3D point)
    {
        for (int i = 0; i < _zones.Length; i++)
        {
            if (point.X <= _xZones[_zones[i][1]] && point.Y <= _yZones[_zones[i][3]] && point.Z <= _zZones[_zones[i][5]] &&
                point.X >= _xZones[_zones[i][0]] && point.Y >= _yZones[_zones[i][2]] && point.Z >= _zZones[_zones[i][4]])
            {
                return _sigmaValues[i];
            }
        }

        throw new Exception("Can't find eligible zone for sigma");
    }

    public void AccountBoundaryConditions()
    {
        for (int ielem = 0; ielem < Elements.Length; ielem++)
        {
            if (ielem < _sumSteps[0] * _sumSteps[1])
            {
                if (_boundaries[2] == 1) DirichletBoundary(ElementSide.Bottom, ielem);
            }

            if (ielem >= _sumSteps[0] * _sumSteps[1] * _sumSteps[2] - _sumSteps[0] * _sumSteps[1] || _sumSteps[2] == 1)
            {
                if (_boundaries[3] == 1) DirichletBoundary(ElementSide.Upper, ielem);
            }

            if (ielem % _sumSteps[0] == 0)
            {
                if (_boundaries[0] == 1) DirichletBoundary(ElementSide.Left, ielem);
            }

            if ((ielem + 1) % _sumSteps[0] == 0)
            {
                if (_boundaries[1] == 1) DirichletBoundary(ElementSide.Right, ielem);
            }

            if (ielem % (_sumSteps[0] * _sumSteps[1]) < _sumSteps[0])
            {
                if (_boundaries[5] == 1) DirichletBoundary(ElementSide.Front, ielem);
            }

            if (ielem % (_sumSteps[0] * _sumSteps[1]) >= _sumSteps[0] * _sumSteps[1] - _sumSteps[0])
            {
                if (_boundaries[4] == 1) DirichletBoundary(ElementSide.Rear, ielem);
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
                DirichletBoundaries.Add(Elements[ielem][4]);
                DirichletBoundaries.Add(Elements[ielem][5]);
                break;
            
            case ElementSide.Upper:
                DirichletBoundaries.Add(Elements[ielem][2]);
                DirichletBoundaries.Add(Elements[ielem][3]);
                DirichletBoundaries.Add(Elements[ielem][6]);
                DirichletBoundaries.Add(Elements[ielem][7]);
                break;
            
            case ElementSide.Left:
                DirichletBoundaries.Add(Elements[ielem][4]);
                DirichletBoundaries.Add(Elements[ielem][8]);
                DirichletBoundaries.Add(Elements[ielem][10]);
                DirichletBoundaries.Add(Elements[ielem][6]);
                break;
                
            case ElementSide.Right:
                DirichletBoundaries.Add(Elements[ielem][5]);
                DirichletBoundaries.Add(Elements[ielem][9]);
                DirichletBoundaries.Add(Elements[ielem][11]);
                DirichletBoundaries.Add(Elements[ielem][7]);
                break;
            
            case ElementSide.Front:
                DirichletBoundaries.Add(Elements[ielem][0]);
                DirichletBoundaries.Add(Elements[ielem][8]);
                DirichletBoundaries.Add(Elements[ielem][9]);
                DirichletBoundaries.Add(Elements[ielem][2]);
                break;
            
            case ElementSide.Rear:
                DirichletBoundaries.Add(Elements[ielem][1]);
                DirichletBoundaries.Add(Elements[ielem][10]);
                DirichletBoundaries.Add(Elements[ielem][11]);
                DirichletBoundaries.Add(Elements[ielem][3]);
                break;
        }
    }
}

public interface ITimeGrid
{
    public double[] TGrid { get; set;  }
    
    public double this[int index]
    {
        get => TGrid[index];
        set => TGrid[index] = value;
    }
}

public class TimeGrid : ITimeGrid
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

public class GeneratedTimeGrid : ITimeGrid
{
    public double[] TGrid { get; set;  }

    public GeneratedTimeGrid(string path)
    {
        using (var sr = new StreamReader(path))
        {
            string[] data;
            data = sr.ReadToEnd()!.Split("\r\n", StringSplitOptions.RemoveEmptyEntries);
            TGrid = new double[(data.Length - 1) / 2];
            TGrid = data.Select(Convert.ToDouble).ToArray();
        }
    }

    public double this[int index]
    {
        get => TGrid[index];
        set => TGrid[index] = value;
    }

}
