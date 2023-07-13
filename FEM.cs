namespace VectorFEM3D;

public class FEM
{
    private SparseMatrix? _globalMatrix;
    private Vector? _globalVector;
    private Vector? _solution;
    private Vector? _localVector;
    private Matrix? _stiffnessMatrix;
    private Matrix? _massMatrix;
    private Grid? _grid;
    private Test? _test;
    private IBasis3D? _basis;
    private Integration? _integration;

    public FEM(Grid grid, TimeGrid timeGrid)
    {
        _grid = grid;
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

    private void AssemblySLAE(int itime)
    {
        _globalVector.Fill(0);
        _globalMatrix.Clear();

        for (int ielem = 0; ielem < _grid.Elements.Length; ielem++)
        {
            AssemblyLocalElement(ielem, itime);
            
            ////////////////////////////////////
            /// добавить время
            ////////////////////////////////////

            _stiffnessMatrix += _massMatrix;

            for (int i = 0; i < _basis.Size; i++)
            {
                for (int j = 0; j < _basis.Size; j++)
                {
                    AddElement(_grid.Elements[ielem][i], _grid.Elements[ielem][j], _stiffnessMatrix[i, j]);
                }
            }
            
            AssemblyGlobalVector(ielem, 0);
            
            _stiffnessMatrix.Clear();
            _massMatrix.Clear();
            _localVector.Fill(0);
        }
    }

    private void AssemblyGlobalVector(int ielem, int itime)
    {
        for (int i = 0; i < _basis.Size; i++)
        {
            _globalVector[_grid.Elements[ielem][i]] += _localVector[i];
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

                int ik = i;
                int jk = j;
                kek = point =>
                {
                    Vector3D psi1 = _basis.GetPsi(ik, point);
                    Vector3D psi2 = _basis.GetPsi(jk, point);

                    return psi1 * psi2;
                };

                _massMatrix[i, j] = hx * hy * hz * _integration.Gauss3D(kek);

                kek = point =>
                {
                    Vector3D dPsi1 = _basis.GetDPsi(ik, point);
                    Vector3D dPsi2 = _basis.GetDPsi(ik, point);

                    return Vector3D.DotProductJacob(dPsi1, dPsi2, hx, hy, hz);
                };

                _stiffnessMatrix[i, j] = _grid.Lambda * _integration.Gauss3D(kek);
            }

            _localVector[i] = _test.F(_grid.Edges[_grid.Elements[ielem][i]].Point, 0);
        }

        _localVector = _massMatrix * _localVector;

        _massMatrix = _grid.Sigma * _massMatrix;
    }

    private void BuildPortrait()
    {
        HashSet<int>[] list = new HashSet<int>[_grid.Nodes.Length].Select(_ => new HashSet<int>()).ToArray();
        foreach (var element in _grid.Elements)
        foreach (var pos in element)
        foreach (var node in element)
            if (pos > node)
                list[pos].Add(node);

        list = list.Select(childlist => childlist.Order().ToHashSet()).ToArray();
        int count = list.Sum(childlist => childlist.Count);

        _globalMatrix = new(_grid.Nodes.Length, count);
        _globalVector = new(_grid.Nodes.Length);
        //_layers = new Vector[3].Select(_ => new Vector(_grid.Nodes.Length)).ToArray();

        _globalMatrix.Ig[0] = 0;

        for (int i = 0; i < list.Length; i++)
            _globalMatrix.Ig[i + 1] = _globalMatrix.Ig[i] + list[i].Count;

        int k = 0;

        foreach (var childlist in list)
        foreach (var value in childlist)
            _globalMatrix.Jg[k++] = value;
    }
}