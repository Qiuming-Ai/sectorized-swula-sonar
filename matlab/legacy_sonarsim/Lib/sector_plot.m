function sector_plot(angles, ranges, img, dynRange, type)
    % 扇面图绘制（平面波几何：x= r*sinθ, z= r*cosθ）
    % angles : 1×Na (rad)
    % ranges : Nr×1 (m)
    % img    : Nr×Na (复值或实值强度)
    % dynRange : (可选) 显示动态范围，默认 60 dB
    
    if nargin < 4 || isempty(dynRange), dynRange = 60;type = "sector"; end
    if nargin < 5 || isempty(type), type = "sector"; end
    
    angles = angles(:).';                % 1×Na
    ranges = ranges(:);                  % Nr×1
    [Nr, Na] = size(img);
    

    I = abs(img);
    I = I ./ (max(I(:)) + eps);
    IdB = 20*log10(I + eps);
    IdB = max(IdB, -dynRange);           % 限动态范围
    

    
    % 极坐标→笛卡尔 (x 横向, z 距离/深度)
    [TH, RR] = meshgrid(angles, ranges); % Nr×Na
    X = RR .* sin(TH);
    Z = RR .* cos(TH);
    
    % 绘图
    % figure; 
    if type == "sector"
        surf(X, Z, IdB, 'EdgeColor','none'); 
        view(2); 
        % axis equal tight;
        set(gca, 'YDir','reverse');          % 让近处在上、远处在下（声呐常用）
        colormap("jet"); c = colorbar; 
        c.Label.String = 'Amplitude (dB)';
        caxis([-dynRange 0]);
        xlabel('Lateral X (m)'); ylabel('Range Z (m)');
        title('DAS Sector Image');
        set(gcf,'Color','w');                % 黑背景可改为 'k'
    elseif type == "rect"
        imagesc(angles/pi*180,ranges, IdB);
        colormap("hot");
        c = colorbar;
        c.Label.String = 'Amplitude (dB)';
        caxis([-dynRange 0]);
        % xlabel('Angle X (°)'); ylabel('Range Z (m)');
        % title('DAS Image');
    else
        close
    end
    
    
    
    end
    