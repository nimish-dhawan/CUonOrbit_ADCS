function cubesat_attitude_sfcn(block)
    setup(block);
end

function setup(block)
    block.NumInputPorts  = 2;
    block.NumOutputPorts = 0;

    block.SetPreCompInpPortInfoToDynamic;

    block.InputPort(1).Dimensions = 3;  % RPY [rad]
    block.InputPort(2).Dimensions = 3;  % Sun vector [km]

    block.SampleTimes = [0 0];
    block.SimStateCompliance = 'DefaultSimState';

    block.RegBlockMethod('Start', @attitudeViewerStart);
    block.RegBlockMethod('Outputs', @attitudeViewerUpdate);
end

function attitudeViewerStart(block)
    % close all;

    % Parameters
    arrow_len = 10;             % in km
    cubesat_origin = [0; 0; 0]; % fixed position at origin

    att_fig = figure('Name', 'CubeSat Attitude', 'NumberTitle', 'off');
    hold on; axis equal;
    xlim([-20 20]); ylim([-20 20]); zlim([-20 20]);
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    view(3); grid on;
    title('3U CubeSat Attitude in ECI Frame');

    % Body axes
    h_att_x = quiver3(0,0,0,0,0,0,arrow_len,'r','LineWidth',2);
    h_att_y = quiver3(0,0,0,0,0,0,arrow_len,'g','LineWidth',2);
    h_att_z = quiver3(0,0,0,0,0,0,arrow_len,'b','LineWidth',2);

    % Sun vector
    h_sun_vec = quiver3(0,0,0,0,0,0,arrow_len, ...
        'Color',[1 0.6 0],'LineWidth',2);

    % 3U CubeSat body patch
    Lx = 30; Ly = 10; Lz = 10;
    [xv, yv, zv] = ndgrid([-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5]);
    cube_verts = [ ...
        xv(:) * Lx, ...
        yv(:) * Ly, ...
        zv(:) * Lz ...
    ];
    cube_faces = [1 3 7 5; 2 4 8 6; 1 2 6 5;
                  3 4 8 7; 1 2 4 3; 5 6 8 7];
    h_cubesat = patch('Vertices', cube_verts, 'Faces', cube_faces, ...
        'FaceColor', 'cyan', 'EdgeColor', 'k', ...
        'FaceAlpha', 0.3, 'LineWidth', 1);

    % Dummy legend handles
    l1 = plot3(NaN,NaN,NaN,'r-','LineWidth',1);
    l2 = plot3(NaN,NaN,NaN,'g-','LineWidth',1);
    l3 = plot3(NaN,NaN,NaN,'b-','LineWidth',1);
    l4 = plot3(NaN,NaN,NaN,'c-','LineWidth',1.5);
    l5 = plot3(NaN,NaN,NaN,'-','LineWidth',1.5,'Color',[1 0.6 0]);

    legend([l1 l2 l3 l5 l4], ...
        {'+X_{body}', '+Y_{body}', '+Z_{body}', 'Sun Vector', '3U CubeSat'}, ...
        'TextColor','k','Location','northeast');

    set_param(block.BlockHandle, 'UserData', struct( ...
        'attFig', att_fig, ...
        'qX', h_att_x, ...
        'qY', h_att_y, ...
        'qZ', h_att_z, ...
        'qSun', h_sun_vec, ...
        'cubePatch', h_cubesat, ...
        'cubeVerts', cube_verts, ...
        'cubePos', cubesat_origin, ...
        't_last', -inf ...
    ));
end

function attitudeViewerUpdate(block)
    ud = get_param(block.BlockHandle, 'UserData');
    t_now = block.CurrentTime;

    % Only update every 5s
    if t_now - ud.t_last < 5
        return;
    end
    ud.t_last = t_now;

    rpy = block.InputPort(1).Data(:);
    sun_vec = block.InputPort(2).Data(:);

    R_bi = angle2dcm(rpy(3), rpy(2), rpy(1), 'ZYX');
    xb = R_bi(1,:)'; yb = R_bi(2,:)'; zb = R_bi(3,:)';

    % Update triads
    p = ud.cubePos;
    set(ud.qX, 'XData', p(1), 'YData', p(2), 'ZData', p(3), ...
               'UData', xb(1), 'VData', xb(2), 'WData', xb(3));
    set(ud.qY, 'XData', p(1), 'YData', p(2), 'ZData', p(3), ...
               'UData', yb(1), 'VData', yb(2), 'WData', yb(3));
    set(ud.qZ, 'XData', p(1), 'YData', p(2), 'ZData', p(3), ...
               'UData', zb(1), 'VData', zb(2), 'WData', zb(3));

    % Sun vector
    sun_unit = sun_vec / norm(sun_vec);
    set(ud.qSun, 'XData', p(1), 'YData', p(2), 'ZData', p(3), ...
                 'UData', sun_unit(1), 'VData', sun_unit(2), 'WData', sun_unit(3));

    % Cube rotation
    rotated_verts = (R_bi * ud.cubeVerts')';
    set(ud.cubePatch, 'Vertices', rotated_verts + p');

    set_param(block.BlockHandle, 'UserData', ud);
    drawnow limitrate
end
