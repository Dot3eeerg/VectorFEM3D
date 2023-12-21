using VectorFEM3D;

Grid grid = new Grid("GridParameters");
grid.BuildGrid();
grid.AccountBoundaryConditions();

TimeGrid timeGrid = new TimeGrid("TimeGridParameters");
timeGrid.BuildTimeGrid();

FEM fem = new FEM(grid, timeGrid);

fem.SetTest(new Test2(grid));

fem.SetSolver(new BCGSTABLUSolver(1e-16, 1000));
//fem.SetSolver(new LUSolver());

fem.SetScheme(Scheme.Three_layer_Implicit);
//fem.SetScheme(Scheme.Four_layer_Implicit);

fem.Compute();