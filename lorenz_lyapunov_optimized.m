% =========================================================
% اختيار الشروط الابتدائية لتحقيق أقصى فوضى في نظام لورنز
% مع حساب أسية ليابونوف وعرضها في نافذة مخصصة
% =========================================================

clear; clc; close all;

% -----------------------------
% معاملات نظام لورنز
% -----------------------------
sigma = 10;
rho   = 28;
beta  = 8/3;

% -----------------------------
% إعدادات المحاكاة
% -----------------------------
tspan = [0 50];  % تقليل الفترة الزمنية للسرعة
dt    = 0.01;

% عدد الحالات الابتدائية المراد فحصها
Ncases = 100;

% حدود الشروط الابتدائية
x_min = -20; x_max = 20;
y_min = -20; y_max = 20;
z_min = -20; z_max = 20;

% تجاوز الشرط الأولي
transient = 10;  % تجاهل أول 10 وحدات زمنية

% -----------------------------
% توليد شروط ابتدائية عشوائية
% -----------------------------
rng(1);  % لضمان نفس النتائج في كل تشغيل
ICs = [rand(Ncases,1)*(x_max-x_min)+x_min, ...
       rand(Ncases,1)*(y_max-y_min)+y_min, ...
       rand(Ncases,1)*(z_max-z_min)+z_min];

% فترة زمنية موحدة لجميع الحالات
t_eval = linspace(tspan(1), tspan(2), 5000);

% -----------------------------
% حساب أسية ليابونوف لكل حالة
% -----------------------------
lambda_values = zeros(Ncases,1);
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);

fprintf('جاري حساب أسية ليابونوف...\n');

for i = 1:Ncases
    IC = ICs(i,:);
    
    try
        % حل نظام لورنز للشرط الأساسي
        [t1, X] = ode45(@(t, x) lorenz_ode(t, x, sigma, rho, beta), ...
                         tspan, IC, options);
        
        % إعادة استيفاء النقاط لجعلها موحدة
        X_interp = interp1(t1, X, t_eval, 'spline');
        
        % شرط ابتدائي مضطرب قليلاً
        perturbation = 1e-8;
        IC_perturbed = IC + perturbation * randn(1,3);
        
        % حل نظام لورنز للشرط المضطرب
        [t2, X_perturbed] = ode45(@(t, x) lorenz_ode(t, x, sigma, rho, beta), ...
                                   tspan, IC_perturbed, options);
        
        % إعادة استيفاء النقاط
        X_perturbed_interp = interp1(t2, X_perturbed, t_eval, 'spline');
        
        % حساب الفارق بين المسارين
        % تجنب فترة التحول الأولية
        idx_start = find(t_eval >= transient, 1);
        idx_start = max(idx_start, 1);
        
        % حساب المسافة بين المسارين
        divergence = sqrt(sum((X_interp(idx_start:end,:) - X_perturbed_interp(idx_start:end,:)).^2, 2));
        
        % تجنب صفر واللانهاية
        divergence(divergence < 1e-15) = 1e-15;
        divergence(isinf(divergence)) = [];
        
        if ~isempty(divergence)
            % أسية ليابونوف = متوسط (log(divergence) - log(perturbation)) / t
            t_active = t_eval(idx_start:idx_start+length(divergence)-1) - t_eval(idx_start);
            t_active = max(t_active, 1e-10);
            
            log_divergence = log(divergence / perturbation) ./ t_active;
            lambda_values(i) = mean(log_divergence);
        else
            lambda_values(i) = 0;
        end
        
    catch
        lambda_values(i) = 0;  % إذا حدث خطأ، استخدم صفر
    end
    
    % طباعة التقدم
    if mod(i, 25) == 0
        fprintf('تم معالجة %d من %d حالة\n', i, Ncases);
    end
end

% تصفية القيم غير المنطقية
lambda_values(isinf(lambda_values) | isnan(lambda_values) | lambda_values < -10) = 0;

% تحويل القيم السالبة إلى صفر إذا لزم الأمر (الفوضى تتطلب قيم موجبة)
lambda_values = max(lambda_values, 0);

% -----------------------------
% اختيار الحالة ذات أعلى قيمة Lyapunov
% -----------------------------
[max_lambda, best_idx] = max(lambda_values);
best_IC = ICs(best_idx,:);

fprintf('\n========================================\n');
fprintf('نتائج التحليل:\n');
fprintf('========================================\n');
fprintf('أكبر قيمة لأسية ليابونوف = %.6f\n', max_lambda);
fprintf('الشروط الابتدائية المختارة:\n');
fprintf('  x0 = %.4f\n  y0 = %.4f\n  z0 = %.4f\n', best_IC(1), best_IC(2), best_IC(3));
fprintf('رقم الحالة = %d من %d\n', best_idx, Ncases);
fprintf('========================================\n\n');

% -----------------------------
% حل نظام لورنز للحالة المختارة (نطاق زمني أطول)
% -----------------------------
tspan_long = [0 100];
[t_best, X_best] = ode45(@(t, x) lorenz_ode(t, x, sigma, rho, beta), ...
                  tspan_long, best_IC, options);

% =========================================================
% الرسوم البيانية
% =========================================================

% الرسم 1: مسار لورنز ثلاثي الأبعاد
fig1 = figure('Name', 'مسار لورنز 3D', 'NumberTitle', 'off', 'Color','w', 'Position', [100 100 900 700]);
plot3(X_best(:,1), X_best(:,2), X_best(:,3), 'b', 'LineWidth', 1.5);
grid on; 
xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('z', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf(['مسار لورنز الأمثل\nIC = [%.3f, %.3f, %.3f], λ = %.6f'], ...
    best_IC(1), best_IC(2), best_IC(3), max_lambda), ...
    'FontSize', 12, 'FontWeight', 'bold');
view(45, 25);
axis tight;
rotate3d on;

% الرسم 2: توزيع أسية ليابونوف
fig2 = figure('Name', 'توزيع Lyapunov', 'NumberTitle', 'off', 'Color','w', 'Position', [1050 100 900 700]);

% رسم شريطي للقيم
subplot(2,1,1);
[sorted_lambda, ~] = sort(lambda_values, 'descend');
num_top = min(50, length(sorted_lambda));
bar(sorted_lambda(1:num_top), 0.8, 'FaceColor', [0.2 0.6 0.9], 'EdgeColor', 'b');
hold on;
plot(1, max_lambda, 'r*', 'MarkerSize', 15, 'LineWidth', 2.5);
xlabel('الترتيب (أكبر القيم أولاً)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('λ (Lyapunov Exponent)', 'FontSize', 11, 'FontWeight', 'bold');
title('أفضل 50 قيمة لأسية ليابونوف', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
ylim([0, max(sorted_lambda)*1.1]);

% رسم خطي لجميع القيم
subplot(2,1,2);
plot(1:Ncases, lambda_values, 'b-o', 'LineWidth', 1.2, 'MarkerSize', 4);
hold on;
plot(best_idx, max_lambda, 'r*', 'MarkerSize', 15, 'LineWidth', 2.5);
xlabel('رقم الحالة الابتدائية', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('λ', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('أسية ليابونوف لكل حالة ابتدائية (الأكبر = %.6f)', max_lambda), ...
    'FontSize', 12, 'FontWeight', 'bold');
grid on;
legend('كل الحالات', 'القيمة القصوى', 'Location', 'best');

% الرسم 3: سحابة الشروط الابتدائية
fig3 = figure('Name', 'توزيع الشروط الابتدائية', 'NumberTitle', 'off', 'Color','w', 'Position', [100 850 900 600]);
scatter3(ICs(:,1), ICs(:,2), ICs(:,3), 50, lambda_values, 'filled', 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.3);
hold on;
plot3(best_IC(1), best_IC(2), best_IC(3), 'r*', 'MarkerSize', 20, 'LineWidth', 2.5);
colormap(jet);
cbar = colorbar;
set(cbar, 'Label', 'λ (Lyapunov Exponent)', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('x0', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y0', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('z0', 'FontSize', 12, 'FontWeight', 'bold');
title('توزيع الشروط الابتدائية ملونة حسب Lyapunov', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
view(45, 30);
rotate3d on;

% الرسم 4: تطور المتغيرات مع الزمن
fig4 = figure('Name', 'تطور المتغيرات الزمني', 'NumberTitle', 'off', 'Color','w', 'Position', [1050 850 900 600]);

subplot(2,2,1);
plot(t_best, X_best(:,1), 'b', 'LineWidth', 1.5);
grid on;
xlabel('الزمن', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('x(t)', 'FontSize', 11, 'FontWeight', 'bold');
title('تطور المتغير x', 'FontSize', 11, 'FontWeight', 'bold');

subplot(2,2,2);
plot(t_best, X_best(:,2), 'r', 'LineWidth', 1.5);
grid on;
xlabel('الزمن', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('y(t)', 'FontSize', 11, 'FontWeight', 'bold');
title('��طور المتغير y', 'FontSize', 11, 'FontWeight', 'bold');

subplot(2,2,3);
plot(t_best, X_best(:,3), 'g', 'LineWidth', 1.5);
grid on;
xlabel('الزمن', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('z(t)', 'FontSize', 11, 'FontWeight', 'bold');
title('تطور المتغير z', 'FontSize', 11, 'FontWeight', 'bold');

subplot(2,2,4);
% حساب norm للمسار
norm_X = sqrt(X_best(:,1).^2 + X_best(:,2).^2 + X_best(:,3).^2);
plot(t_best, norm_X, 'purple', 'LineWidth', 1.5);
grid on;
xlabel('الزمن', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('||X(t)||', 'FontSize', 11, 'FontWeight', 'bold');
title('حجم المتجه (Norm)', 'FontSize', 11, 'FontWeight', 'bold');

% =========================================================
% رسم إضافي: مقارنة عدة حالات
% =========================================================
fig5 = figure('Name', 'مقارنة أفضل 5 حالات', 'NumberTitle', 'off', 'Color','w', 'Position', [500 500 1000 700]);

% اختيار أفضل 5 حالات
[sorted_lambda_vals, sorted_indices] = sort(lambda_values, 'descend');
top_5_idx = sorted_indices(1:min(5, length(sorted_indices)));

colors = {'b', 'r', 'g', 'm', 'c'};

for j = 1:length(top_5_idx)
    idx = top_5_idx(j);
    [~, X_j] = ode45(@(t, x) lorenz_ode(t, x, sigma, rho, beta), ...
                      tspan_long, ICs(idx,:), options);
    
    plot3(X_j(:,1), X_j(:,2), X_j(:,3), 'Color', colors{j}, 'LineWidth', 1.5, ...
          'DisplayName', sprintf('λ=%.6f', lambda_values(idx)));
    hold on;
end

grid on;
xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('z', 'FontSize', 12, 'FontWeight', 'bold');
title('مقارنة أفضل 5 حالات ابتدائية', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');
view(45, 25);
rotate3d on;

fprintf('اكتمل التحليل بنجاح!\n');

% =========================================================
% دالة نظام لورنز
% =========================================================
function dx = lorenz_ode(~, x, sigma, rho, beta)
    dx = zeros(3,1);
    dx(1) = sigma * (x(2) - x(1));
    dx(2) = x(1) * (rho - x(3)) - x(2);
    dx(3) = x(1) * x(2) - beta * x(3);
end
