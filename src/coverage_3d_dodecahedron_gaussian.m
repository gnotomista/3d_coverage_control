%% flags
clc
clear
close all

DOMAIN  = 'dodecahedron';
DENSITY = 'Gaussian';

PRINT_COVERAGE_COST = false;
PLOT_GENERATORS = true;
PLOT_CENTROIDS = false;
PLOT_VORONOI_CELLS = false;

%% initialization of domain
% boundary
if strcmp(DOMAIN, 'cube')
    domain = permn([-1, 1], 3);
elseif strcmp(DOMAIN, 'dodecahedron')
    [domain, ~] = platonic_solids(5, 1, false);
end
[Ad, bd] = vert2lcon(domain);

% density function
if strcmp(DENSITY, 'constant')
    PHI = @(x,y,z) 1;
elseif strcmp(DENSITY, 'Gaussian')
    PHI_C = [0.5, 0.5, 0.5];
    PHI = @(x,y,z) 1e4 * 2.^(-(x-PHI_C(1)).^2/0.5^2-(y-PHI_C(2)).^2/0.5^2-(z-PHI_C(3)).^2/0.5^2);
end

%% initialization of agents
N = 10;
DT = 0.02;
p = -0.01 + 0.02 * rand(3, N); % positions of the agents
V_MAX = 10; % max speed of the agents

%% initialization of 3D coverage control algorithm object
vc = VoronoiCoverage3D(N, domain);
vc.setDensity(PHI);
vc.setRobotsPositions(p);

%% plot initialization
initialize_plot

%% main control loop
while true
    tic
    vc.setRobotsPositions(p);
    vc.weightedCoverageControl();
    pDot = 1 / DT * (vc.G - p);
    pDot = limit_speed(pDot, V_MAX);
    p = p + pDot * DT;

    if PRINT_COVERAGE_COST, disp(['Coverage cost:', num2str(vc.coverageCost(p, PHI))]), end
    if PLOT_GENERATORS, vc.plotGenerators(), end
    if PLOT_CENTROIDS, vc.plotCentroids(), end
    if PLOT_VORONOI_CELLS, vc.plotVoronoiCells(), end

    drawnow limitrate
    pause(DT-toc)
end
