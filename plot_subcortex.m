function plot_subcortex(cifti_data, output_plot_path, color_map, min_scale, max_scale, structure_label)

    if isstruct(cifti_data)
        cifti_data = cifti_data.cdata;
    end
    
    if isempty(color_map)
        % Define the number of points for each segment
        n = 128; % Half of 256

        % -------------------------------
        % Define Red to Grey Transition
        % -------------------------------

        % Red channel: 1 to 0.5
        R_red_to_grey = linspace(1, 0.5, n);

        % Green channel: 0 to 0.5
        G_red_to_grey = linspace(0, 0.5, n);

        % Blue channel: 0 to 0.5
        B_red_to_grey = linspace(0, 0.5, n);

        % -------------------------------
        % Define Grey to Blue Transition
        % -------------------------------

        % Red channel: 0.5 to 0
        R_grey_to_blue = linspace(0.5, 0, n);

        % Green channel: 0.5 to 0
        G_grey_to_blue = linspace(0.5, 0, n);

        % Blue channel: 0.5 to 1
        B_grey_to_blue = linspace(0.5, 1, n);

        % -------------------------------
        % Combine the Segments
        % -------------------------------

        % Concatenate the segments for each RGB channel
        R = [R_red_to_grey, R_grey_to_blue];
        G = [G_red_to_grey, G_grey_to_blue];
        B = [B_red_to_grey, B_grey_to_blue];

        % Combine into a single colormap matrix
        cmap = [R', G', B'];
        color_map = flipud(cmap);
    end

    if strcmp(structure_label,'hippocampus')
        structure_model_idx = [14,15];
        view_pos = [-150.6371   14.3183];
    elseif strcmp(structure_label,'striatum')
        structure_model_idx = [3,4,8,9,18,19];
        view_pos = [-165.3716   28.1047];
    elseif strcmp(structure_label,'thalamus')
        structure_model_idx = [20,21];
        view_pos = [-33.4297   33.1598];
    elseif strcmp(structure_label,'amygdala')
        structure_model_idx = [5,6];
        view_pos = [-150.6371   14.3183];
    elseif strcmp(structure_label,'cerebellum')
        structure_model_idx = [10,11];
        view_pos = [-23.0028   13.5360];
    elseif strcmp(structure_label,'brainstem')
        structure_model_idx = 7;
        view_pos = [-150.6371   14.3183];
    else
        error('No such subcortical area!!');
    end
    
    labels_subcortical = cifti_read('atlas/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii');
    
    cifti_data_select = [];
    label_select = [];
    vox_position = [];
    
    models = labels_subcortical.diminfo{1}.models;
    sform = labels_subcortical.diminfo{1, 1}.vol.sform;
    
    for n_model = structure_model_idx
        cifti_data_select = [cifti_data_select;
            cifti_data(models{1,n_model}.start:(models{1,n_model}.start+models{1,n_model}.count-1))];
        label_select = [label_select;
            labels_subcortical.cdata(models{1,n_model}.start:(models{1,n_model}.start+models{1,n_model}.count-1))];
        MNI_coords_temp = sform  * [models{1,n_model}.voxlist;ones(1,models{1,n_model}.count)];
        vox_position = [vox_position, MNI_coords_temp(1:3,:)];
    end
    
    color = zeros(size(vox_position,2), 3);
    for n_vox = 1:size(vox_position,2)
        if cifti_data_select(n_vox) <= max_scale && cifti_data_select(n_vox) > min_scale
            val = cifti_data_select(n_vox);
        elseif cifti_data_select(n_vox) <= min_scale
            val = min_scale + eps;
        elseif cifti_data_select(n_vox) > max_scale
            val = max_scale;
        end

        color(n_vox, :) = color_map(ceil((val-min_scale)/(max_scale-min_scale) * 256), :);
    end
    
    scatter_size = 50;
    
    figure;
    scatter3(vox_position(1,:)', vox_position(2,:)', vox_position(3,:)', ...
        scatter_size, color, 'filled');

    % Adjusting figure to get better visualization
    axis equal; axis vis3d;
    % Adjusting figure size.
    set(gcf, 'Position', [0, 0, 800, 800]);
    xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis');
    % Hidding the axis.
    set(gca, 'xticklabel', []); set(gca, 'yticklabel', []); set(gca, 'zticklabel', []);
    set(gca, 'xtick', [], 'ytick', [], 'ztick', [], 'xcolor', 'w', 'ycolor', 'w', 'zcolor', 'w', 'Visible', 0);
    % Setting background transparent.
    set(gca, 'color', 'none'); set(gcf, 'color', 'none');
    
    ax = gca;
    ax.View = view_pos;
    
    hold on; % Enable adding more graphics to the current plot
    
    % Retrieve Camera Properties
    cam_pos = ax.CameraPosition;
    
    if abs(cam_pos(2) - min(ylim)) < abs(cam_pos(2) - max(ylim))
        arrow_origin = [0,min(ylim),min(zlim)-10];
        if strcmp(structure_label,'cerebellum')
            arrow_origin = arrow_origin - [0,20,10];
        end
        posterior_dx = -2;
        posterior_dy = -3;
        right_dz = -2;
    else
        arrow_origin = [0,max(ylim),min(zlim)-15];
        posterior_dx = -2;
        posterior_dy = 2;
        right_dz = -4;
    end
    
    % Determine the range of the data to position the arrows appropriately
    x_limits = xlim;
    y_limits = ylim;
    z_limits = zlim;

    % Define arrow length as a fraction of the axis limits
    arrow_length = min([diff(x_limits), diff(y_limits), diff(z_limits)]) * 0.3;

%     % Define starting points for arrows (e.g., from one corner)
%     arrow_origin = [x_limits(1), y_limits(1), z_limits(1)];
    
    % Right Arrow (+X)
    quiver3(arrow_origin(1), arrow_origin(2), arrow_origin(3), arrow_length, 0, 0, ...
        'LineWidth', 4, 'Color', [0.5,0.5,0.5], 'MaxHeadSize', 0.5, 'AutoScale','off');
    text(arrow_origin(1) + arrow_length, arrow_origin(2), arrow_origin(3) + right_dz, 'Right', 'FontSize', 24, 'Color', [0.5,0.5,0.5]);

    % Posterior Arrow (+Y)
    quiver3(arrow_origin(1), arrow_origin(2), arrow_origin(3), 0, -arrow_length, 0, ...
        'LineWidth', 4, 'Color', [0.5,0.5,0.5], 'MaxHeadSize', 0.5, 'AutoScale','off');
    text(arrow_origin(1) + posterior_dx, arrow_origin(2) - arrow_length + posterior_dy, arrow_origin(3), 'Posterior', 'FontSize', 24, 'Color', [0.5,0.5,0.5]);

    % Dorsal Arrow (+Z)
    quiver3(arrow_origin(1), arrow_origin(2), arrow_origin(3), 0, 0, arrow_length, ...
        'LineWidth', 4, 'Color', [0.5,0.5,0.5], 'MaxHeadSize', 0.5, 'AutoScale','off');
    text(arrow_origin(1), arrow_origin(2), arrow_origin(3) + arrow_length + 3, 'Dorsal', 'FontSize', 24, 'Color', [0.5,0.5,0.5]);

    
%     figure('Name', 'Custom Red-Grey-Blue Colormap', 'NumberTitle', 'off');
%     imagesc(1:256);
%     colormap(color_map);
%     colorbar;
%     title('Custom Red-Grey-Blue Colormap with 256 Points');
    
    if ~isempty(output_plot_path)
        saveas(gcf,output_plot_path)
    end

end

