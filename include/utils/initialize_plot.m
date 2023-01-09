vc.weightedCoverageControl(); % compute centroids, Voronoi cells to initialize plot

if strcmp(DENSITY, 'constant')
    PLOT_DENSITY = false;
elseif strcmp(DENSITY, 'Gaussian')
    PLOT_DENSITY = true;
end

figure, hold on
hD = plot(alphaShape(domain(:, 1), domain(:, 2), domain(:, 3)), ...
    'FaceAlpha', 0.25, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
if PLOT_DENSITY
    step = 0.05;
    X = min(domain(:,1)) : step : max(domain(:,1));
    Y = min(domain(:,2)) : step : max(domain(:,2));
    Z = min(domain(:,3)) : step : max(domain(:,3));
    [Xm, Ym, Zm] = meshgrid(X, Y, Z);
    PHIm = PHI(Xm, Ym, Zm);
    phimv = PHIm(:);
    for l = 1e4 * [0.1 0.2 0.4 0.6 0.8 0.9]
        [~, idx] = min(abs(phimv-l));
        ll = phimv(idx);
        [ff, vv, cc] = isosurface(Xm, Ym, Zm, PHIm, ll, ll * ones(size(PHIm)));
        ff_inside = [];
        for i = 1 : size(ff, 1)
            if all(Ad * vv(ff(i,1),:)' <= bd) ...
                    && all(Ad * vv(ff(i,2),:)' <= bd) ...
                    && all(Ad * vv(ff(i,3),:)' <= bd)
                ff_inside = [ff_inside; ff(i,:)];
            end
        end
        patch('Vertices', vv, 'Faces', ff_inside, 'FaceVertexCData', cc, ...
            'FaceAlpha', 0.5, 'FaceColor', 'interp', 'EdgeColor', 'none')
    end
end
axis equal
axis(0.1 * [-1 1 -1 1 -1 1] + [min(domain(:,1)) max(domain(:,1)) min(domain(:,2)) max(domain(:,2)) min(domain(:,3)) max(domain(:,3))])
set(gca, 'Visible', 'off')
view(60, 20)
rotate3d

if PLOT_VORONOI_CELLS
    delete(hD)
end
if PLOT_GENERATORS, vc.plotGenerators(), end
if PLOT_CENTROIDS, vc.plotCentroids(), end
if PLOT_VORONOI_CELLS, vc.plotVoronoiCells(), end