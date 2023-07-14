using VectorFEM3D;

Grid grid = new Grid("GridParameters");
grid.BuildGrid();
grid.AccountBoundaryConditions();

TimeGrid timeGrid = new TimeGrid("TimeGridParameters");
timeGrid.BuildTimeGrid();

FEM fem = new FEM(grid, timeGrid);
fem.SetTest(new Test1(grid));
fem.SetSolver(new BCGSTABLUSolver(1e-14, 1000));
//fem.SetSolver(new LUSolver());
fem.Compute();