function orbit_plot_sfcn(block)
    setup(block);
end

function setup(block)
    block.NumInputPorts  = 5;
    block.NumOutputPorts = 0;

    block.SetPreCompInpPortInfoToDynamic;

    % Input 1: Satellite position [x; y; z] (km)
    block.InputPort(1).Dimensions = 3;
    block.InputPort(1).DatatypeID = 0;
    block.InputPort(1).Complexity = 'Real';
    block.InputPort(1).DirectFeedthrough = true;

    % Input 2: Attitude [roll; pitch; yaw] (rad)
    block.InputPort(2).Dimensions = 3;
    block.InputPort(2).DatatypeID = 0;
    block.InputPort(2).Complexity = 'Real';
    block.InputPort(2).DirectFeedthrough = true;

    % Input 3: Sun vector [s_x; s_y; s_z] (ECI, km)
    block.InputPort(3).Dimensions = 3;
    block.InputPort(3).DatatypeID = 0;
    block.InputPort(3).Complexity = 'Real';
    block.InputPort(3).DirectFeedthrough = true;

    % Input 4: Angular velocity [wx; wy; wz] (rad/s)
    block.InputPort(4).Dimensions = 3;
    block.InputPort(4).DatatypeID = 0;
    block.InputPort(4).Complexity = 'Real';
    block.InputPort(4).DirectFeedthrough = true;

    % Input 5: Orbit Period (s)
    block.InputPort(5).Dimensions = 1;
    block.InputPort(5).DatatypeID = 0;
    block.InputPort(5).Complexity = 'Real';
    block.InputPort(5).DirectFeedthrough = true;

    block.SampleTimes = [0 0]; % Inherited sample time
    block.SimStateCompliance = 'DefaultSimState';

    block.RegBlockMethod('Start', @Start);
    block.RegBlockMethod('Outputs', @Outputs);
end

function Start(block)
    % close(findall(0, 'Type', 'figure', 'Name', 'Satellite Orbit Viewer'));
    close all

    % Parameters
    R_earth = 6378 * 0.8;   % km, scaled-down Earth
    scale = 1000;           % Arrow length (km)

    [x, y, z] = sphere(35);
    fig = figure('Name', 'Satellite Orbit Viewer', 'NumberTitle', 'off');

    % Create tiled layout
    t = tiledlayout(fig, 2, 2);
    t.TileSpacing = 'loose';
    t.Padding = 'loose';

    % Tile 1-3: 3D Orbit Plot
    ax3d = nexttile(t, 1, [2 1]);
    load('topo.mat', 'topo');
    surf(ax3d, R_earth*x, R_earth*y, R_earth*z, ...
        'CData', topo, ...
        'FaceColor', 'texturemap', ...
        'EdgeColor', [0.4 0.4 0.4], ...
        'LineStyle', '-', ...
        'LineWidth', 0.3, ...
        'FaceAlpha', 1);
    hold(ax3d, 'on');
    axis(ax3d, 'equal');
    xlim(ax3d, [-8e3 8e3]); ylim(ax3d, [-8e3 8e3]); zlim(ax3d, [-8e3 8e3]);
    xlabel(ax3d, 'X [km]'); ylabel(ax3d, 'Y [km]'); zlabel(ax3d, 'Z [km]');
    view(ax3d, 3);
    grid(ax3d, 'off');
    box(ax3d, 'on')
    title(ax3d, 'Real-Time Satellite Orbit & Attitude');

    h_x = quiver3(ax3d, 0, 0, 0, 0, 0, 0, scale, 'r', 'LineWidth', 3);
    h_y = quiver3(ax3d, 0, 0, 0, 0, 0, 0, scale, 'g', 'LineWidth', 3);
    h_z = quiver3(ax3d, 0, 0, 0, 0, 0, 0, scale, 'b', 'LineWidth', 3);
    h_trail = plot3(ax3d, NaN, NaN, NaN, 'k', 'LineWidth', 1.5);
    h_sun = quiver3(ax3d, 0, 0, 0, 0, 0, 0, 1500, 'Color', [1 0.6 0], 'LineWidth', 2);

    % Dummy handles for legend
    hx_legend = plot3(NaN, NaN, NaN, 'r-', 'LineWidth', 1);
    hy_legend = plot3(NaN, NaN, NaN, 'g-', 'LineWidth', 1);
    hz_legend = plot3(NaN, NaN, NaN, 'b-', 'LineWidth', 1);
    trail_legend = plot3(NaN, NaN, NaN, 'k-', 'LineWidth', 1.5);
    sun_legend = plot3(NaN, NaN, NaN, '-', 'LineWidth', 1.5, ...
                       'Color', [1 0.6 0]);

    legend([hx_legend, hy_legend, hz_legend, trail_legend, sun_legend], ...
           {'+X_{body}', '+Y_{body}', '+Z_{body}', 'Orbit Trail', 'Sun Vector'}, ...
           'TextColor', 'k', 'Location', 'northeast');

    % Tile 4: RPY History
    ax_rpy = nexttile(t, 2); hold(ax_rpy, 'on'); grid(ax_rpy, 'off'); box(ax_rpy, 'on');
    h_roll = plot(ax_rpy, NaN, NaN, 'r', 'LineWidth', 1, 'DisplayName', '\phi (Roll)');
    h_pitch = plot(ax_rpy, NaN, NaN, 'g', 'LineWidth', 1, 'DisplayName', '\theta (Pitch)');
    h_yaw = plot(ax_rpy, NaN, NaN, 'b', 'LineWidth', 1, 'DisplayName', '\psi (Yaw)');
    xlabel(ax_rpy, 'Orbit'); ylabel(ax_rpy, 'Angle [deg]');
    legend(ax_rpy, 'Location', 'southoutside','Orientation','horizontal');
    title(ax_rpy, 'Spacecraft Attitude');

    % Tile 6: Angular Velocity History
    ax_rates = nexttile(t, 4); hold(ax_rates, 'on'); grid(ax_rates, 'off'); box(ax_rates, 'on');
    h_wx = plot(ax_rates, NaN, NaN, 'r', 'LineWidth', 1, 'DisplayName', '\omega_x');
    h_wy = plot(ax_rates, NaN, NaN, 'g', 'LineWidth', 1, 'DisplayName', '\omega_y');
    h_wz = plot(ax_rates, NaN, NaN, 'b', 'LineWidth', 1, 'DisplayName', '\omega_z');
    xlabel(ax_rates, 'Orbit'); ylabel(ax_rates, '\omega [deg/s]');
    legend(ax_rates, 'Location', 'southoutside', 'Orientation', 'horizontal');
    title(ax_rates, 'Spacecraft Angular Rate');

    set_param(block.BlockHandle, 'UserData', struct( ...
        'fig', fig, 'ax3d', ax3d, 'h_x', h_x, 'h_y', h_y, 'h_z', h_z, ...
        'h_sun', h_sun, 'h_trail', h_trail, 'trail_data', [], 'last_update_time', -inf, ...
        'ax_rpy', ax_rpy, 'h_roll', h_roll, 'h_pitch', h_pitch, 'h_yaw', h_yaw, ...
        'time_history', [], 'rpy_history', [], ...
        'ax_rates', ax_rates, 'h_wx', h_wx, 'h_wy', h_wy, 'h_wz', h_wz, ...
        'omega_history', [] ...
    ));
end

function Outputs(block)
    t_now = block.CurrentTime;
    ud = get_param(block.BlockHandle, 'UserData');

    if t_now - ud.last_update_time < 10
        return;
    end
    ud.last_update_time = t_now;

    pos = block.InputPort(1).Data(:);
    rpy = block.InputPort(2).Data(:);
    sun_vec = block.InputPort(3).Data(:);
    omega = block.InputPort(4).Data(:);
    T = block.InputPort(5).Data(:);

    % Attitude DCM and body axes
    R = angle2dcm(rpy(3), rpy(2), rpy(1), 'ZYX');
    x_axis = R(1,:)'; y_axis = R(2,:)'; z_axis = R(3,:)';

    set(ud.h_x, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3), 'UData', x_axis(1), 'VData', x_axis(2), 'WData', x_axis(3));
    set(ud.h_y, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3), 'UData', y_axis(1), 'VData', y_axis(2), 'WData', y_axis(3));
    set(ud.h_z, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3), 'UData', z_axis(1), 'VData', z_axis(2), 'WData', z_axis(3));

    sun_dir = sun_vec / norm(sun_vec);
    set(ud.h_sun, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3), 'UData', sun_dir(1), 'VData', sun_dir(2), 'WData', sun_dir(3));

    ud.trail_data = [ud.trail_data; pos'];
    set(ud.h_trail, 'XData', ud.trail_data(:,1), 'YData', ud.trail_data(:,2), 'ZData', ud.trail_data(:,3));

    % RPY & Omega history
    ud.time_history = [ud.time_history; t_now];
    ud.rpy_history = [ud.rpy_history; rpy(:)' * 180/pi];
    ud.omega_history = [ud.omega_history; omega(:)'];

    set(ud.h_roll, 'XData', ud.time_history/T, 'YData', ud.rpy_history(:,1));
    set(ud.h_pitch, 'XData', ud.time_history/T, 'YData', ud.rpy_history(:,2));
    set(ud.h_yaw,  'XData', ud.time_history/T, 'YData', ud.rpy_history(:,3));

    set(ud.h_wx, 'XData', ud.time_history/T, 'YData', ud.omega_history(:,1));
    set(ud.h_wy, 'XData', ud.time_history/T, 'YData', ud.omega_history(:,2));
    set(ud.h_wz, 'XData', ud.time_history/T, 'YData', ud.omega_history(:,3));

    set_param(block.BlockHandle, 'UserData', ud);
    drawnow limitrate;
end
