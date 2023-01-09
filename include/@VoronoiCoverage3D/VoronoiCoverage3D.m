classdef VoronoiCoverage3D < handle

    properties
        p % 3xN
        N
        domain
        phi
        pMirrored
        VoronoiCells
        G
        VoronoiCellsToPlot
        % Aaug
        % baug
        colors
        hP
        hG
        hT
        numberGaussPoints = 3;
    end

    methods
        function obj = VoronoiCoverage3D(N, domain)
            obj.N = N;
            obj.domain = [domain; domain(1, :)];
            obj.VoronoiCells = {};
            obj.colors = rand(N, 3);
            obj.hP = cell(N, 1);
            obj.hG = cell(N, 1);
            obj.hT = cell(N, 1);
        end

        function setRobotsPositions(obj, p)
            obj.p = p;
        end

        function setDensity(obj,phi)
            obj.phi = phi;
        end

        function weightedCoverageControl(obj)
            P = [obj.p];
            [~, V, ~, ~] = polybnd_voronoi(P',obj.domain);
            % [~, V, obj.Aaug, obj.baug] = polybnd_voronoi(P', obj.domain);
            obj.VoronoiCellsToPlot = V;
            obj.G = zeros(3,obj.N);
            for i = 1 : obj.N
                vertex = V{i};
                obj.VoronoiCells{i} = [vertex(:,1) vertex(:,2) vertex(:,3)];
                obj.G(:,i) = obj.centroidWeighted(vertex);
            end
        end

        function H = coverageCost(obj, robotXYZ, phi)
            H = 0;
            for i = 1: size(robotXYZ, 3)
                xP = obj.VoronoiCells{i}(:,1);
                yP = obj.VoronoiCells{i}(:,2);
                zP = obj.VoronoiCells{i}(:,3);
                f = @(x,y,z) ((robotXYZ(1,i)-x).^2+(robotXYZ(2,i)-y).^2+(robotXYZ(3,i)-z).^2).* phi(x,y,z);
                trngltn = delaunay(xP, yP, zP);
                H_i = 0;
                for j = 1 : size(trngltn, 1)
                    H_i = H_i + VoronoiCoverage3D.intOfFOverT(f, obj.numberGaussPoints, obj.VoronoiCells{i}(trngltn(j,:),:));
                end
                H = H + H_i;
            end
        end

        function massWeighted = weightedCoverageControlMass(obj)
            P = [obj.p];
            [~, V, ~, ~] = polybnd_voronoi(P', obj.domain);
            massWeighted = zeros(3,obj.N);
            for i = 1 : obj.N
                vertex = V{i};
                massWeighted(:,i) = obj.massWeighted(vertex);
            end
        end

        function plotGenerators(obj)
            for i = 1 : obj.N
                if isempty(obj.hP{i})
                    obj.hP{i} = scatter3(obj.p(1, :), obj.p(2, :), obj.p(3, :), ...
                        100, 'Marker', 'o', 'MarkerFaceColor', obj.colors(i,:)*0.8, 'MarkerEdgeColor', 'none');
                else
                    obj.hP{i}.XData = obj.p(1,i);
                    obj.hP{i}.YData = obj.p(2,i);
                    obj.hP{i}.ZData = obj.p(3,i);
                end
            end
        end

        function plotCentroids(obj)
            for i = 1 : obj.N
                if isempty(obj.hG{i})
                    obj.hG{i} = scatter3(obj.G(1, i), obj.G(2, i), obj.G(3, i), ...
                        200, 'Marker', 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', obj.colors(i,:)*0.6);
                else
                    obj.hG{i}.XData = obj.G(1,i);
                    obj.hG{i}.YData = obj.G(2,i);
                    obj.hG{i}.ZData = obj.G(3,i);
                end
            end
        end

        function plotVoronoiCells(obj)
            for i = 1 : obj.N
                K = convhulln(obj.VoronoiCellsToPlot{i});
                if isempty(obj.hT{i})
                    obj.hT{i} = trisurf(K, obj.VoronoiCellsToPlot{i}(:,1), obj.VoronoiCellsToPlot{i}(:,2), obj.VoronoiCellsToPlot{i}(:,3), ...
                        'FaceColor', obj.colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                else
                    obj.hT{i}.Faces = K;
                    obj.hT{i}.Vertices = obj.VoronoiCellsToPlot{i};
                end
            end
        end
    end

    methods (Access = private)
        function G = centroidWeighted(obj, P)
            xP = P(:,1);
            yP = P(:,2);
            zP = P(:,3);
            phiA = @(x,y,z) max(eps,obj.phi(x,y,z));
            phiSx = @(x,y,z) x.*max(eps,obj.phi(x,y,z));
            phiSy = @(x,y,z) y.*max(eps,obj.phi(x,y,z));
            phiSz = @(x,y,z) z.*max(eps,obj.phi(x,y,z));
            trngltn = delaunay(xP, yP, zP);
            A = 0;
            S = 0;
            for i = 1 : size(trngltn, 1)
                A = A + VoronoiCoverage3D.intOfFOverT(phiA, obj.numberGaussPoints, P(trngltn(i,:),:));
                S = S + [VoronoiCoverage3D.intOfFOverT(phiSx, obj.numberGaussPoints, P(trngltn(i,:),:));
                    VoronoiCoverage3D.intOfFOverT(phiSy, obj.numberGaussPoints, P(trngltn(i,:),:));
                    VoronoiCoverage3D.intOfFOverT(phiSz, obj.numberGaussPoints, P(trngltn(i,:),:))];
            end
            G = S/A;
        end

        function A = massWeighted(obj, P)
            xP = P(:,1);
            yP = P(:,2);
            zP = P(:,3);
            phiA = @(x,y,z) max(0.05,obj.phi(x,y,z));
            trngltn = delaunay(xP, yP, zP);
            A = 0;
            for i = 1 : size(trngltn, 1)
                A = A + VoronoiCoverage3D.intOfFOverT(phiA, obj.numberGaussPoints, P(trngltn(i,:),:));
            end
        end
    end

    methods (Static)
        function I = intOfFOverT(f, N, T)
            x1 = T(1,1);
            x2 = T(2,1);
            x3 = T(3,1);
            x4 = T(4,1);
            y1 = T(1,2);
            y2 = T(2,2);
            y3 = T(3,2);
            y4 = T(4,2);
            z1 = T(1,3);
            z2 = T(2,3);
            z3 = T(3,3);
            z4 = T(4,3);
            xyw = VoronoiCoverage3D.TriGaussPoints(N);
            V = abs(((x1-x4)*(y2-y4)*(z3-z4))+((x2-x4)*(y3-y4)*(z1-z4))+...
                ((x3-x4)*(y1-y4)*(z2-z4))-(((x3-x4)*(y2-y4)*(z1-z4))+((x2-x4)*(y1-y4)*(z3-z4))+...
                ((x1-x4)*(y3-y4)*(z2-z4))))/6;
            NP = size(xyw(:,1), 1);
            I = 0;
            for j = 1 : NP
                x = x4*(1-xyw(j,1)-xyw(j,2)-xyw(j,3))+x3*xyw(j,1)+x2*xyw(j,2)+x1*xyw(j,3);
                y = y4*(1-xyw(j,1)-xyw(j,2)-xyw(j,3))+y3*xyw(j,1)+y2*xyw(j,2)+y1*xyw(j,3);
                z = z4*(1-xyw(j,1)-xyw(j,2)-xyw(j,3))+z3*xyw(j,1)+z2*xyw(j,2)+z1*xyw(j,3);
                I = I + f(x,y,z)*xyw(j,4);
            end
            I = V*I;
        end

        function xw = TriGaussPoints(n)
            xw = zeros(n,4);
            if n==1
                xw = [0.25 0.25 0.25 1];
            elseif n==2
                xw = [0.585410196624969 0.138196601125011 0.138196601125011  0.25
                    0.138196601125011 0.138196601125011 0.585410196624969  0.25
                    0.138196601125011 0.585410196624969 0.138196601125011 0.25
                    0.138196601125011 0.138196601125011 0.138196601125011 0.25];
            elseif n==3
                xw = [0.25 0.25 0.25 -0.8
                    0.5 0.16666666666667 0.16666666666667 0.45
                    0.16666666666667 0.16666666666667 0.5 0.45
                    0.16666666666667 0.5 0.16666666666667 0.45
                    0.16666666666667 0.16666666666667 0.16666666666667 0.45];
            elseif n==4
                xw = [0.25 0.25 0.25 -0.01315555555556
                    0.7857142857142856 0.071428571428571 0.071428571428571 0.00762222222222
                    0.071428571428571 0.7857142857142856 0.071428571428571 0.00762222222222
                    0.071428571428571 0.071428571428571 0.7857142857142856 0.00762222222222
                    0.071428571428571 0.071428571428571 0.071428571428571 0.00762222222222
                    0.399403576166799 0.399403576166799 0.100596423833201 0.02488888888889
                    0.399403576166799 0.100596423833201 0.399403576166799 0.02488888888889
                    0.100596423833201 0.399403576166799 0.399403576166799 0.02488888888889
                    0.100596423833201 0.100596423833201 0.399403576166799 0.02488888888889
                    0.399403576166799 0.100596423833201 0.100596423833201 0.02488888888889
                    0.100596423833201 0.399403576166799 0.100596423833201 0.02488888888889];
            elseif n==5
                xw = [0.25 0.25 0.25 0.030283678097089
                    0 0.33333333333333 0.33333333333333 0.006026785714286
                    0.33333333333333 0 0.33333333333333 0.006026785714286
                    0.33333333333333 0.33333333333333 0 0.006026785714286
                    0.33333333333333 0.33333333333333 0.33333333333333 0.006026785714286
                    0.727272727272727 0.090909090909091 0.090909090909091 0.011645249086029
                    0.090909090909091 0.727272727272727 0.090909090909091 0.011645249086029
                    0.090909090909091 0.090909090909091 0.727272727272727 0.011645249086029
                    0.090909090909091 0.090909090909091 0.090909090909091 0.011645249086029];
            end
        end
    end
end
