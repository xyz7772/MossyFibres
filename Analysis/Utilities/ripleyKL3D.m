function [r_vals, K_vals, L_vals] = ripleyKL3D(xyz, letplot)
    % Inputs:
    %   xyz: An N×3 matrix of point coordinates.
    %   letplot: 1 to plot figures
    %
    % Outputs:
    %   r_vals: A vector of radii.
    %   K_vals: The corresponding values of Ripley's K function.
    %   L_vals: The corresponding values of Ripley's L function.

    r_max = 50;      % Maximum search radius (μm)
    n_r = 50;         % Number of radius samples
    r_vals = linspace(0, r_max, n_r); % Radius vector

    % Study volume
    x_min = min(xyz(:,1)); x_max = max(xyz(:,1));
    y_min = min(xyz(:,2)); y_max = max(xyz(:,2));
    z_min = min(xyz(:,3)); z_max = max(xyz(:,3));
    studyVolume = (x_max - x_min) * (y_max - y_min) * (z_max - z_min);

    N = size(xyz,1);
    lambda = N / studyVolume; % Point density

    % Compute Euclidean distances (efficient pairwise)
    D = squareform(pdist(xyz));

    % Initialize K values
    K_vals = zeros(size(r_vals));

    % Compute edge correction factors (distance to boundaries)
    d_xmin = xyz(:,1) - x_min;
    d_xmax = x_max - xyz(:,1);
    d_ymin = xyz(:,2) - y_min;
    d_ymax = y_max - xyz(:,2);
    d_zmin = xyz(:,3) - z_min;
    d_zmax = z_max - xyz(:,3);

    for i = 1:length(r_vals)
        r = r_vals(i);
        count = sum(D <= r, 2) - 1; % Count neighbors within radius, exclude self

        % Compute effective volume correction
        V_full = (4/3) * pi * r^3; % Full sphere volume
        V_eff = zeros(N,1); % Effective volume for each point

        for j = 1:N
            if r > d_xmin(j) || r > d_xmax(j) || r > d_ymin(j) || r > d_ymax(j) || r > d_zmin(j) || r > d_zmax(j)
                % If the sphere is cut by any boundary, calculate the corrected volume
                V_cap_xmin = cap_volume(r, d_xmin(j));
                V_cap_xmax = cap_volume(r, d_xmax(j));
                V_cap_ymin = cap_volume(r, d_ymin(j));
                V_cap_ymax = cap_volume(r, d_ymax(j));
                V_cap_zmin = cap_volume(r, d_zmin(j));
                V_cap_zmax = cap_volume(r, d_zmax(j));

                % Wedges: where sphere intersects two planes
                V_wedge_xmin_ymin = wedge_volume(r, d_xmin(j), d_ymin(j));
                V_wedge_xmin_ymax = wedge_volume(r, d_xmin(j), d_ymax(j));
                V_wedge_xmax_ymin = wedge_volume(r, d_xmax(j), d_ymin(j));
                V_wedge_xmax_ymax = wedge_volume(r, d_xmax(j), d_ymax(j));

                V_wedge_xmin_zmin = wedge_volume(r, d_xmin(j), d_zmin(j));
                V_wedge_xmin_zmax = wedge_volume(r, d_xmin(j), d_zmax(j));
                V_wedge_xmax_zmin = wedge_volume(r, d_xmax(j), d_zmin(j));
                V_wedge_xmax_zmax = wedge_volume(r, d_xmax(j), d_zmax(j));

                % Semi-wedge: where sphere intersects three planes
                V_semiwedge = semi_wedge_volume(r, d_xmin(j), d_ymin(j), d_zmin(j));

                % Effective volume
                V_eff(j) = V_full - (V_cap_xmin + V_cap_xmax + V_cap_ymin + V_cap_ymax + V_cap_zmin + V_cap_zmax) ...
                                    + (V_wedge_xmin_ymin + V_wedge_xmin_ymax + V_wedge_xmax_ymin + V_wedge_xmax_ymax ...
                                    + V_wedge_xmin_zmin + V_wedge_xmin_zmax + V_wedge_xmax_zmin + V_wedge_xmax_zmax) ...
                                    - V_semiwedge;
            else
                % No boundary interference
                V_eff(j) = V_full;
            end
        end

        % Compute correction factor
        edge_correction = V_eff / V_full;
        K_vals(i) = sum(count ./ edge_correction) / (N * lambda);
    end

    % Compute L function
    L_vals = (K_vals / (4/3 * pi)).^(1/3);

    % Plot
    if letplot == 1
        figure(Position=[100 100 650 300]);

        subplot(1,2,1);
        plot(r_vals, K_vals, '-', 'color','k','LineWidth', 2);
        hold on;
        K_csr = (4/3)*pi*(r_vals.^3);
        plot(r_vals, K_csr, '--', 'color',[0.5,0.5,0.5], 'LineWidth', 2);
        xlabel('r (μm)'); ylabel('K(r)');
        legend('Observed', 'CSR','Location','best','box','off');
        set(gca, 'LineWidth', 1, 'FontSize', 13)

        subplot(1,2,2);
        plot(r_vals, L_vals, '-', 'color','k','LineWidth', 2);
        hold on;
        plot(r_vals, r_vals, 'r--', 'color',[0.5,0.5,0.5],'LineWidth', 2);
        xlabel('r (μm)'); ylabel('L(r)');
        legend('Observed', 'CSR','Location','best','box','off');
        set(gca, 'LineWidth', 1, 'FontSize', 13)
    end
end

function V_cap = cap_volume(r, d)
    if d > r
        V_cap = 0;
    else
        h = r - d;
        V_cap = (pi * h^2 * (3*r - h)) / 3;
    end
end

function V_wedge = wedge_volume(r, d1, d2)
    V_wedge = cap_volume(r, d1) + cap_volume(r, d2);
end

function V_semiwedge = semi_wedge_volume(r, d1, d2, d3)
    V_semiwedge = cap_volume(r, d1) + cap_volume(r, d2) + cap_volume(r, d3);
end
