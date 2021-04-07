clear all; format long; close all
load('3dIsingData_overleaf.mat');
T = T(chi < 60); D = D(chi < 60); xi2 = xi2(chi < 60); m = m(chi < 60); chi = chi(chi < 60);
idx = [848, 849];
T(idx) = []; D(idx) = []; xi2(idx) = []; m(idx) = []; chi(idx) = [];

m = abs(m .* xi2.^(beta / nu));
T = T(m < 0.5); D = D(m < 0.5); chi = chi(m < 0.5); m = m(m < 0.5);

ft = @(a, b, c, d, e, f, g, h, x)(a + (c + g .* x + h .* x.^2) .* abs(x).^(abs(nu / beta + abs(d) ./ (abs(e) + abs(x).^abs(b)).^abs(f))));

f = {};
for d = unique(D)
    for ksi = unique(chi)
        t = T(D == d & chi == ksi);
        [mm, ord] = sort(m(D == d & chi == ksi), 'ascend');
        t = t(ord);
        t = t(mm > 1e-2); mm = mm(mm > 1e-2);
        if numel(t) > 10 && max(mm) > 0.4 && min(mm) < 0.3
            f{end + 1} = fit(mm(:), t(:), ft, 'StartPoint', [4.55, 7, -0.2, -0.6, 4.8, 27, 0.5, -1]);
            % plot(mm,t,'.')
            % hold on
            % xx = linspace(min(mm),max(mm),101);
            % plot(xx,f{end}(xx),'-');
            % hold off
        end
    end
end

numel(f)

crossT = []; crossm = [];
for i = 1:numel(f)
    for j = i + 1:numel(f)
        [x, fval] = fzero(@(x)(f{i}(x) - f{j}(x)), 0.2);
        if ~isnan(fval)
            crossT = [crossT; f{i}(x)]; crossm = [crossm; x];
        end
    end
end

plot(crossT, crossm, '.')

mean(crossT)
var(crossT)

mean(crossm)
var(crossm)
