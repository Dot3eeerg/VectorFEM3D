namespace VectorFEM3D;

public class FEM
{
    private SparseMatrix? _globalMatrix;
    private Vector? _globalVector;
    private Vector[]? _layers;
    private Vector? _solution;
    private Vector? _localVector;
    private Matrix? _stiffnessMatrix;
    private Matrix? _massMatrix;
    private Grid? _grid;
    private TimeGrid _timeGrid;
    private Test? _test;
    private IBasis3D? _basis;
    private Integration? _integration;
    private SLAE? _slae;
    private Scheme _scheme;

    public FEM(Grid grid, TimeGrid timeGrid)
    {
        _grid = grid;
        _timeGrid = timeGrid;
        _basis = new TriLinearVectorBasis();
        _integration = new Integration(Quadratures.SegmentGaussOrder9());
        _stiffnessMatrix = new(_basis.Size);
        _massMatrix = new(_basis.Size);
        _localVector = new(_basis.Size);
    }

    public void SetTest(Test test)
    {
        _test = test;
    }

    public void SetSolver(SLAE slae)
    {
        _slae = slae;
    }

    public void SetScheme(Scheme scheme)
    {
        _scheme = scheme;
    }

    public void Compute()
    {
        BuildPortrait();
        PrepareLayers();

        int itime = 0;
        
        switch (_scheme)
        {
            case Scheme.Three_layer_Implicit:
                itime = 2;
                break;
            
            case Scheme.Four_layer_Implicit:
                itime = 3;
                break;
        }

        for ( ; itime < _timeGrid.TGrid.Length; itime++)
        {
            AssemblySLAE(itime);
            AccountDirichletBoundaries(itime);
            
            _slae.SetSLAE(_globalVector, _globalMatrix);
            _solution = _slae.Solve();

            switch (_scheme)
            {
                case Scheme.Three_layer_Implicit:
                    Vector.Copy(_layers[1], _layers[0]);
                    Vector.Copy(_solution, _layers[1]);
                    break;
                
                case Scheme.Four_layer_Implicit:
                    Vector.Copy(_layers[1], _layers[0]);
                    Vector.Copy(_layers[2], _layers[1]);
                    Vector.Copy(_solution, _layers[2]);
                    break;
            }
            
            PrintError(itime);
        }
    }

    private void AssemblySLAE(int itime)
    {
        _globalVector.Fill(0);
        _globalMatrix.Clear();

        for (int ielem = 0; ielem < _grid.Elements.Length; ielem++)
        {
            AssemblyLocalElement(ielem, itime);
            
            _stiffnessMatrix += SchemeUsage(ielem, itime, _scheme, 0) * _massMatrix;

            for (int i = 0; i < _basis.Size; i++)
            {
                for (int j = 0; j < _basis.Size; j++)
                {
                    AddElement(_grid.Elements[ielem][i], _grid.Elements[ielem][j], _stiffnessMatrix[i, j]);
                }
            }
            
            AssemblyGlobalVector(ielem, itime);
            
            _stiffnessMatrix.Clear();
            _massMatrix.Clear();
            _localVector.Fill(0);
        }
    }

    private void AssemblyGlobalVector(int ielem, int itime)
    {
        double[] qj3 = new double[_basis.Size];
        double[] qj2 = new double[_basis.Size];
        double[] qj1 = new double[_basis.Size];
        
        switch (_scheme)
        {
            case Scheme.Three_layer_Implicit:

                for (int i = 0; i < _basis.Size; i++)
                {
                    for (int j = 0; j < _basis.Size; j++)
                    {
                        qj2[i] += _massMatrix[i, j] * _layers[1][_grid.Elements[ielem][j]];
                        qj1[i] += _massMatrix[i, j] * _layers[0][_grid.Elements[ielem][j]];
                    }
                }
                
                for (int i = 0; i < _basis.Size; i++)
                {
                    _localVector[i] += SchemeUsage(ielem, itime, _scheme, 1) * qj2[i];
                    _localVector[i] += SchemeUsage(ielem, itime, _scheme, 2) * qj1[i];
                    
                    _globalVector[_grid.Elements[ielem][i]] += _localVector[i];
                }
                
                break;
            
            case Scheme.Four_layer_Implicit:

                for (int i = 0; i < _basis.Size; i++)
                {
                    for (int j = 0; j < _basis.Size; j++)
                    {
                        qj3[i] += _massMatrix[i, j] * _layers[2][_grid.Elements[ielem][j]];
                        qj2[i] += _massMatrix[i, j] * _layers[1][_grid.Elements[ielem][j]];
                        qj1[i] += _massMatrix[i, j] * _layers[0][_grid.Elements[ielem][j]];
                    }
                }
                
                for (int i = 0; i < _basis.Size; i++)
                {
                    _localVector[i] += SchemeUsage(ielem, itime, _scheme, 1) * qj3[i];
                    _localVector[i] += SchemeUsage(ielem, itime, _scheme, 2) * qj2[i];
                    _localVector[i] += SchemeUsage(ielem, itime, _scheme, 3) * qj1[i];
                    
                    _globalVector[_grid.Elements[ielem][i]] += _localVector[i];
                }
                
                break;
        }
    }

    private void AddElement(int i, int j, double value)
    {
        if (i == j)
        {
            _globalMatrix.Di[i] += value;
            return;
        }

        if (i > j)
        {
            for (int icol = _globalMatrix.Ig[i]; icol < _globalMatrix.Ig[i + 1]; icol++)
            {
                if (_globalMatrix.Jg[icol] == j)
                {
                    _globalMatrix.Ggl[icol] += value;
                    return;
                }
            }
        }

        else
        {
            for (int icol = _globalMatrix.Ig[j]; icol < _globalMatrix.Ig[j + 1]; icol++)
            {
                if (_globalMatrix.Jg[icol] == i)
                {
                    _globalMatrix.Ggu[icol] += value;
                    return;
                }
            }
        }
    }

    private void AssemblyLocalElement(int ielem, int itime)
    {
        double hx = _grid.Edges[_grid.Elements[ielem][0]].Length;
        double hy = _grid.Edges[_grid.Elements[ielem][2]].Length;
        double hz = _grid.Edges[_grid.Elements[ielem][4]].Length;

        for (int i = 0; i < _basis.Size; i++)
        {
            for (int j = 0; j < _basis.Size; j++)
            {
                Func<Point3D, double> kek;
                Vector3D psi1 = new(0, 0, 0);
                Vector3D psi2 = new(0, 0, 0);
                Vector3D dPsi1 = new(0, 0, 0);
                Vector3D dPsi2 = new(0, 0, 0);

                int ik = i;
                int jk = j;
                kek = point =>
                {
                    psi1.Copy(_basis.GetPsi(ik, point));
                    psi2.Copy(_basis.GetPsi(jk, point));

                    return psi1 * psi2;
                };

                _massMatrix[i, j] += hx * hy * hz * _integration.Gauss3D(kek);

                kek = point =>
                {
                    dPsi1.Copy(_basis.GetDPsi(ik, point));
                    dPsi2.Copy(_basis.GetDPsi(jk, point));

                    return Vector3D.DotProductJacob(dPsi1, dPsi2, hx, hy, hz);
                };

                _stiffnessMatrix[i, j] += 1 / _grid.Mu * _integration.Gauss3D(kek);
            }

            _localVector[i] = _test.F(_grid.Edges[_grid.Elements[ielem][i]].Point, _timeGrid[itime], i, _grid.GetSigma(
                new Point3D(_grid.Edges[_grid.Elements[ielem][0]].Point.X,
                    _grid.Edges[_grid.Elements[ielem][3]].Point.Y,
                    _grid.Edges[_grid.Elements[ielem][11]].Point.Z)));
        }

        _localVector = _massMatrix * _localVector;
    }
    
    private void AccountDirichletBoundaries(int itime)
    {
        foreach (var edge in _grid.DirichletBoundaries)
        {

            _globalMatrix.Di[edge] = 1;
            _globalVector[edge] = _test.UValue(_grid.Edges[edge].Point, _timeGrid[itime], _grid.Edges[edge].GetAxis());

            for (int i = _globalMatrix.Ig[edge]; i < _globalMatrix.Ig[edge + 1]; i++)
                _globalMatrix.Ggl[i] = 0;

            for (int col = edge + 1; col < _globalMatrix.Size; col++)
                for (int j = _globalMatrix.Ig[col]; j < _globalMatrix.Ig[col + 1]; j++)
                    if (_globalMatrix.Jg[j] == edge)
                    {
                        _globalMatrix.Ggu[j] = 0;
                        break;
                    }
        }
    }

    private double SchemeUsage(int ielem, int itime, Scheme scheme, int i)
    {
       double t01 = _timeGrid[itime] - _timeGrid[itime - 1];
       double t02 = _timeGrid[itime] - _timeGrid[itime - 2];
       double t12 = _timeGrid[itime - 1] - _timeGrid[itime - 2];
                
        switch (scheme)
        {
            case Scheme.Three_layer_Implicit:
                switch (i)
                {
                    case 0:
                        return (_grid.GetSigma(new Point3D(_grid.Edges[_grid.Elements[ielem][0]].Point.X,
                                   _grid.Edges[_grid.Elements[ielem][3]].Point.Y,
                                   _grid.Edges[_grid.Elements[ielem][11]].Point.Z)) * (t01 + t02) + 2 * _grid.Epsilon) /
                               (t01 * t02);
                    
                    case 1:
                        return (_grid.GetSigma(new Point3D(_grid.Edges[_grid.Elements[ielem][0]].Point.X,
                            _grid.Edges[_grid.Elements[ielem][3]].Point.Y,
                            _grid.Edges[_grid.Elements[ielem][11]].Point.Z)) * t02 + 2 * _grid.Epsilon) / (t01 * t12);
                    
                    case 2:
                        return -(_grid.GetSigma(new Point3D(_grid.Edges[_grid.Elements[ielem][0]].Point.X,
                            _grid.Edges[_grid.Elements[ielem][3]].Point.Y,
                            _grid.Edges[_grid.Elements[ielem][11]].Point.Z)) * t01 + 2 * _grid.Epsilon) / (t02 * t12);
                    
                }

                break;
            
            case Scheme.Four_layer_Implicit:
                double t03 = _timeGrid[itime] - _timeGrid[itime - 3];
                double t13 = _timeGrid[itime - 1] - _timeGrid[itime - 3];
                double t23 = _timeGrid[itime - 2] - _timeGrid[itime - 3];

                switch (i)
                {
                    case 0:
                        return (_grid.GetSigma(new Point3D(_grid.Edges[_grid.Elements[ielem][0]].Point.X,
                                       _grid.Edges[_grid.Elements[ielem][3]].Point.Y,
                                       _grid.Edges[_grid.Elements[ielem][11]].Point.Z))
                                   * (t01 * t02 + t01 * t03 + t02 * t03) + 2 * _grid.Epsilon * (t01 + t02 + t03)) /
                               (t01 * t02 * t03);
                    
                    case 1:
                        return (_grid.GetSigma(new Point3D(_grid.Edges[_grid.Elements[ielem][0]].Point.X,
                                    _grid.Edges[_grid.Elements[ielem][3]].Point.Y,
                                    _grid.Edges[_grid.Elements[ielem][11]].Point.Z)) * t02 * t03 +
                                2 * _grid.Epsilon * (t02 + t03)) / (t01 * t12 * t13);
                    
                    case 2:
                        return -(_grid.GetSigma(new Point3D(_grid.Edges[_grid.Elements[ielem][0]].Point.X,
                                     _grid.Edges[_grid.Elements[ielem][3]].Point.Y,
                                     _grid.Edges[_grid.Elements[ielem][11]].Point.Z)) * t01 * t03 +
                                 2 * _grid.Epsilon * (t01 + t03)) / (t02 * t12 * t23);
                    
                    case 3:
                        return (_grid.GetSigma(new Point3D(_grid.Edges[_grid.Elements[ielem][0]].Point.X,
                                    _grid.Edges[_grid.Elements[ielem][3]].Point.Y,
                                    _grid.Edges[_grid.Elements[ielem][11]].Point.Z)) * t01 * t02 +
                                2 * _grid.Epsilon * (t01 + t02)) / (t03 * t13 * t23);
                }

                return 0;
            
            default:
                throw new Exception("Undefined scheme");
        }

        throw new Exception("SchemeUsage can't return a value");
    }

    private void BuildPortrait()
    {
        HashSet<int>[] list = new HashSet<int>[_grid.Edges.Length].Select(_ => new HashSet<int>()).ToArray();
        foreach (var element in _grid.Elements)
            foreach (var pos in element)
                foreach (var node in element)
                    if (pos > node)
                        list[pos].Add(node);

        list = list.Select(childlist => childlist.Order().ToHashSet()).ToArray();
        int count = list.Sum(childlist => childlist.Count);

        _globalMatrix = new(_grid.Edges.Length, count);
        _globalVector = new(_grid.Edges.Length);
        _layers = new Vector[3].Select(_ => new Vector(_grid.Edges.Length)).ToArray();

        _globalMatrix.Ig[0] = 0;

        for (int i = 0; i < list.Length; i++)
            _globalMatrix.Ig[i + 1] = _globalMatrix.Ig[i] + list[i].Count;

        int k = 0;

        foreach (var childlist in list)
            foreach (var value in childlist)
                _globalMatrix.Jg[k++] = value;
    }
    
    private void PrepareLayers()
    {
        switch (_scheme)
        {
            case Scheme.Three_layer_Implicit:
                for (int i = 0; i < _grid.Edges.Length; i++)
                {
                    _layers[0][i] = _test.UValue(_grid.Edges[i].Point, _timeGrid[0], _grid.Edges[i].GetAxis());
                    _layers[1][i] = _test.UValue(_grid.Edges[i].Point, _timeGrid[1], _grid.Edges[i].GetAxis());
                }
                
                break;
            
            case Scheme.Four_layer_Implicit:
                for (int i = 0; i < _grid.Edges.Length; i++)
                {
                    _layers[0][i] = _test.UValue(_grid.Edges[i].Point, _timeGrid[0], _grid.Edges[i].GetAxis());
                    _layers[1][i] = _test.UValue(_grid.Edges[i].Point, _timeGrid[1], _grid.Edges[i].GetAxis());
                    _layers[2][i] = _test.UValue(_grid.Edges[i].Point, _timeGrid[2], _grid.Edges[i].GetAxis());
                }
                
                break;
        }
    }

    private void PrintError(int itime)
    {
        double error = 0;
        for (int i = 0; i < _grid.Edges.Length; i++)
        {
            error += Math.Pow(
                _test.UValue(_grid.Edges[i].Point, _timeGrid[itime], _grid.Edges[i].GetAxis()) - _solution[i], 2);
        }
        Console.WriteLine($"Layer error {itime} = {Math.Sqrt(error / _grid.Edges.Length)}");
    }
}