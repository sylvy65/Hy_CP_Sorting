clear all
close all
clc

n_track = 1; 
step_time = 10; 
step_avg = 200; 
peak_lap = 0.2; 
threshold_norm = 0.2; 
threshold_burst_cp = 10; 

file_names = dir('*.csv');

for q = 1:length(file_names)
  
input = readtable(file_names(q).name);

time = table2array(input(:,1));
signal = table2array(input(:,n_track+1));

time_in_seconds = time / 1e6; 
track_length(q) = max(time_in_seconds) - min(time_in_seconds);

disp(['Length of ', file_names(q).name, '(s): ', num2str(track_length(q))]);

signal_norm = signal/max(signal);
signal_mean = movmean(signal_norm,step_avg);
derivative_mean = zeros(1,length(time));
second_der_mean = derivative_mean;

peak = 0; 
j = 1; 
i = 2*step_time; 

for k = i:(length(time)-i)     
    if (i >= (length(time)-peak_lap*10^4))
        break;
    end

   derivative_mean(i) = (signal_mean(i+step_time/2)-signal_mean(i-step_time/2))/(time(i+step_time/2)-time(i-step_time/2));
   second_der_mean(i) = (derivative_mean(i+step_time/2)-derivative_mean(i-step_time/2))/(time(i+step_time/2)-time(i-step_time/2));    
   if (signal_norm(i) >= threshold_norm && derivative_mean(i) > 0 && peak == 0)
        
     peak = 1;
     signal_peak(j,1,q) = time(i);
     signal_peak(j,2,q) = signal(i);
     j = j+1;
     i = i+1;
        
    elseif (second_der_mean(i-step_time) <= 0 && peak == 1)                
        peak = 1;
        j = j-1;
        signal_peak(j,1,q) = time(i);
        signal_peak(j,2,q) = signal(i);
        j = j+1;
        i = i+1;
        
    elseif (second_der_mean(i-step_time) > 0 && peak == 1)
        i = i + peak_lap*10^4;
        peak = 0;
        i = i+1;
             
        else 
        peak = 0;
        i = i+1;       
    end   
end

figure('Name', file_names(q).name);
plot(time/1E6, signal, 'k', 'LineWidth', 1.2);
hold on
plot(signal_peak(:, 1, q)/1E6, signal_peak(:, 2, q), 'o', ...
   'MarkerFaceColor', 'none', ...   
    'MarkerEdgeColor', 'm', ...      
      'MarkerEdgeColor', 'm', ...
      'LineWidth', 1.0, ...
    'MarkerSize', 5);

xlabel('Relative Time (s)');
ylabel('Voltage (uV)');
legend('Original Signal','Extrapolated Peak');
end
 
figure('Name', 'Peaks Temporal Distribution (Raster Plot)');
integral_time_peak = zeros(length(signal_peak),length(file_names));
delta_time_peak = integral_time_peak;
signal_reset = zeros(length(signal_peak),length(file_names));
test_number = zeros(length(signal_peak),length(file_names));

for q = 1:length(file_names) 
signal_reset(:,q) = signal_peak(:,1,q) - signal_peak(1,1,q); 
integral_time_peak = signal_reset/1E6;  
test_number(:,q) = q; 
first_tract(q) = signal_peak(1,1,q)/1E6;
final_track_length(q) = track_length(q) - first_tract(q);
disp(['Final Length of ', file_names(q).name, '(s): ', num2str(final_track_length(q))]);

plot(integral_time_peak(:, q), test_number(:, q), '|', ...
        'MarkerSize', 20, ...
        'LineWidth', 0.5, ...
        'Color', 'k'); 
    hold on
end

ylim([0.6 length(file_names)+0.4])
xlim([0 max(max(integral_time_peak))])
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Hy experiments', 'FontSize', 12, 'FontWeight', 'bold');
title('Peaks Temporal Distribution (Raster Plot)');
if isstruct(file_names)
    yticklabels({file_names.name});
else
    yticklabels(file_names);
end

min_final_track_length = min(final_track_length); 
valid_integral_time_peak = zeros(length(integral_time_peak) - 1, size(integral_time_peak, 2)); 

for s = 1:size(integral_time_peak, 2) 
    i_cp = 1; 
    i_burst = 1; 
  
    for z = 1:(length(integral_time_peak) - 1)
                i_valid = integral_time_peak(z, s);
       
        if i_valid > 0 && i_valid <= min_final_track_length
            
            valid_integral_time_peak(z, s) = i_valid;
                    
                  delta_time_peak(z, s) = (integral_time_peak(z+1,s) - integral_time_peak(z,s)); % (s)
                
                      if delta_time_peak(z, s) > 0 && delta_time_peak(z, s) > threshold_burst_cp
                   
                    IcBI_distribution(i_burst, s) = delta_time_peak(z, s);
                    i_burst = i_burst + 1; 
                    i_cp = 1; 
                     IcBI_distribution(IcBI_distribution ==0) = NaN;
               else
                   if  delta_time_peak(z, s) > 0 && delta_time_peak(z, s) < threshold_burst_cp
                  
                  cp_distribution(i_cp, i_burst, s) = delta_time_peak(z, s);                   i_cp = i_cp + 1; 
                    
                   end
                                end
                     
            burst_counter(s) = i_burst; 
            cp_counter_burst(i_burst, s) = i_cp; 
            
        end            
              
            end
    max_burst(s) = max(burst_counter(s));

end

cp_distribution_file = cell(1, size(cp_distribution, 3));  
for q = 1:size(cp_distribution, 3)
   cp_per_file = cp_distribution(:, :, q);   
   cp_per_file = cp_per_file(:);      
   cp_per_file = cp_per_file(cp_per_file ~= 0);  
   cp_per_file = cp_per_file(~isnan(cp_per_file)); 
   cp_distribution_file{q} = cp_per_file; 
end
cp_distribution_file = cell(1, size(cp_distribution, 3));  

max_cp = 0;
for q = 1:size(cp_distribution, 3)
    cp_per_file = cp_distribution(:, :, q); 
    cp_per_file = cp_per_file(:);  
    cp_per_file = cp_per_file(cp_per_file ~= 0); 
    cp_per_file = cp_per_file(~isnan(cp_per_file));  
    if length(cp_per_file) > max_cp
        max_cp = length(cp_per_file); 
    end
end
cp_counter_burst(cp_counter_burst == 0) = NaN;
disp('(nCPBurst) NUMBER OF CPs PER BURST (REFER TO: cp_counter_burst):');
disp(cp_counter_burst);
disp('(nBurst) NUMBER OF BURSTs PER FILE (REFER TO: burst_counter):');
disp(burst_counter);
cp_distribution_file = NaN(max_cp, size(cp_distribution, 3));  


for q = 1:size(cp_distribution, 3)   
    cp_per_file = cp_distribution(:, :, q);         
    cp_per_file = cp_per_file(:);         
    cp_per_file = cp_per_file(cp_per_file ~= 0);  
    cp_per_file = cp_per_file(~isnan(cp_per_file));  
    cp_distribution_file(1:length(cp_per_file), q) = cp_per_file;
end
disp('(CPI) DURATION (seconds) OF A SINGLE CP PER FILE (REFER TO: cp_distribution_file:');
disp(cp_distribution_file);

figure('Name','Delta Time vs Relative Time');
hold on
colors = lines(length(file_names)); 
for q = 1:length(file_names)
    plot(integral_time_peak(:, q), delta_time_peak(:, q), 'o', 'Color', colors(q, :));
end
legend({file_names.name}, 'Interpreter', 'none');
xlabel('Relative Time (s)')
ylabel('Delta Time (s)')
ylim([0 100])
xlim([0 max(max(integral_time_peak) + 10)])
title('Delta Time vs Relative Time');
grid on
colors = lines(length(file_names)); 
[n_cp, n_burst, n_files] = size(cp_distribution); 
sum_cp_per_burst = zeros(n_burst, n_files); 

for q = 1:n_files
     sum_n_cp_q = sum(cp_counter_burst, 'omitnan');  
       for y = 1:n_burst
               sum_cp_per_burst(y,q) = sum(cp_distribution(:,y,q), 'omitnan');
               sum_cp_per_burst(sum_cp_per_burst == 0) = NaN;
                     end
     sum_interval_cp_per_q(q) = sum(sum_cp_per_burst(:, q), 'omitnan');
end


disp('(nCP) TOTAL NUMBER OF CPs PER FILE (REFER TO: sum_n_cp_q):' );
disp(sum_n_cp_q);
disp('(BurstTime) SUM OF THE DURATION (seconds) OF ALL CPs IN EACH BURST PER FILE (REFER TO: sum_cp_per_burst):');
disp(sum_cp_per_burst); 
disp('IcBIs PER FILE (seconds) (REFER TO: IcBI_distribution)');
disp(IcBI_distribution);   
disp('(CTime) SUM OF THE DURATION OF ALL CPs PER FILE (Contraction Time in seconds) (REFER TO: sum_interval_cp_per_q):');
disp(sum_interval_cp_per_q);


[n_files] = size(IcBI_distribution);
sum_IcBI_per_q = zeros(n_files); 
sum_IcBI_per_q= sum(IcBI_distribution, 'omitnan'); 
disp(' (ETime) SUM OF IcBIs PER FILE (Elongation Time in seconds) (REFER TO: sum_IcBI_per_q)');
disp(sum_IcBI_per_q);   
Hy_Activity_Index = zeros(1, length(sum_interval_cp_per_q));  % Inizializza il vettore

for q = 1:length(sum_interval_cp_per_q)
        Hy_Activity_Index(q) = sum_interval_cp_per_q(q) / (sum_IcBI_per_q(q));
end

disp('(Hy Activity Index) HYDRA ACTIVITY INDEX PER FILE (REFER TO: Hy_Activity_Index)');
disp(Hy_Activity_Index);   

function make_boxplot(data, file_names, y_label, title_text)
    if isempty(data)
        warning(['Dataset vuoto per: ', title_text]);
        return
    end
    [~, n_files] = size(data);
    figure('Name', title_text);
    b = boxplot(data, ...
    'Labels', {file_names.name}, ...
    'Colors', 'k', ...
    'Symbol', 'k+');
    ylabel(y_label);
    xlabel('File');
    title(title_text);
    hold on

    means_data = mean(data, 1, 'omitnan');
    stds_data  = std(data, 0, 1, 'omitnan');
    means_data = means_data(:)';
    stds_data  = stds_data(:)';

    max_val = max(data(:), [], 'omitnan');
    ymax = max([max_val, means_data + stds_data]);
    ymin = min([0, min(means_data - stds_data)]);
    ylim([floor(ymin*1.1) ceil(ymax*1.1)]);

    colors = lines(n_files);
    h = findobj(gca, 'Tag', 'Box');
    for j = 1:length(h)
        patch(get(h(j), 'XData'), get(h(j), 'YData'), ...
            colors(length(h)-j+1, :), ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', colors(length(h)-j+1, :));
    end

    for i = 1:n_files
        x_jitter = (rand(size(data(:, i))) - 0.5) * 0.2;
        scatter(i + x_jitter, data(:, i), ...
            25, 'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', colors(i,:), ...
            'MarkerFaceAlpha', 0.5, 'LineWidth', 0.5);
    end

    errorbar(1:n_files, means_data, stds_data, ...
        'LineStyle', 'none', ...
        'Color', 'r', ...
        'LineWidth', 1.2, ...
        'CapSize', 12);
    plot(1:n_files, means_data, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
end

make_boxplot(cp_distribution_file, file_names, 'CPI (s)', 'CP Duration');
make_boxplot(IcBI_distribution, file_names, 'IcBI (s)', 'Interval between bursts (IcBI)');
make_boxplot(sum_cp_per_burst, file_names, 'Burst Time (s)', 'Sum of CPs per burst in each file');
make_boxplot(cp_counter_burst, file_names, 'nCPs per burst', 'Number of CPs per burst in each file');
durate_tracce = NaN(1, length(file_names)); 
for q = 1:length(file_names)
  
    input = readtable(file_names(q).name, 'VariableNamingRule', 'preserve');
    time = table2array(input(:,1));   
    signal = table2array(input(:,2));  
    recorded_times = time / 1E6;      

    max_valid_integral_time_peak(q) = max(valid_integral_time_peak(:,q));
    max_time_1(q) = first_tract(q);  
    min_time_2(q) = max_valid_integral_time_peak(q) + first_tract(q); 

    valid_indices = (recorded_times >= max_time_1(q)) & (recorded_times <= min_time_2(q));
    time_filtered = recorded_times(valid_indices);
    signal_filtered = signal(valid_indices);

     final_target_time = max_time_1(q) + min_final_track_length;
    if ~isempty(time_filtered) && time_filtered(end) < final_target_time
        extra_indices = (recorded_times > time_filtered(end)) & (recorded_times <= final_target_time);
        time_extra = recorded_times(extra_indices);
        signal_extra = signal(extra_indices);

        time_filtered = [time_filtered; time_extra];
        signal_filtered = [signal_filtered; signal_extra];
       
    end

    time_filtered_aligned = time_filtered - max_time_1(q);
    if ~isempty(time_filtered_aligned)
        durate_tracce(q) = time_filtered_aligned(end);  
    end

    output_table = table(time_filtered_aligned * 1E6, signal_filtered, ...
                         'VariableNames', input.Properties.VariableNames);
    new_file_name = [file_names(q).name(1:end-4), ' - Analyzed_Trace'];
    writetable(output_table, new_file_name);
    disp(['SAVED: ', new_file_name]);
    figure('Name', ['Overlapped traces - ', file_names(q).name]);
    plot(recorded_times, signal, 'Color', [0 0 0], 'LineWidth', 3); 
hold on;

if q <= 5
    plot(time_filtered, signal_filtered, 'Color', [0.00, 0.45, 0.74], 'LineWidth', 0.8); 
else
    plot(time_filtered, signal_filtered, 'Color', [0.85, 0.33, 0.10], 'LineWidth', 0.8); 
end

xlabel('Time (s)');
ylabel('Amplitude');
title(['Traces overlapping - ', file_names(q).name], 'Interpreter', 'none');
legend({'Original Trace', 'Analyzed Trace'}, 'Location', 'best');
grid on

end
disp(['Minimum duration set for all tracks (s): ', num2str(min_final_track_length)]);
disp('Check that all analyzed traces are of the same length (s):');
disp(durate_tracce);

max_time_2 = min_final_track_length; 
n_files = length(file_names);  
if ~exist('valid_integral_time_peak', 'var')
    error('The valid_integral_time_peak matrix is ​​not present in the file.mat');
end
if ~exist('file_names', 'var')
    error('The file_names variable is not present in the file .mat');
end
figure('Name','Peaks alignment to the maximum time value of the shorter trace');
hold on;
set(gcf, 'Color', 'w', 'Position', [100 100 800 600]);

colors = lines(n_files); 
if n_files >= 1
    num_blue = min(5, n_files);
    colors(1:num_blue, :) = repmat([0 0 1], num_blue, 1);
end
if n_files > 5
    start_index = n_files - 4;
    colors(start_index:n_files, :) = repmat([1 0 0], 5, 1);
end
for q = 1:n_files
   
    current_spikes = valid_integral_time_peak(:,q);
    current_spikes = current_spikes(~isnan(current_spikes)); 
    current_spikes = current_spikes(current_spikes >= 0);    
    valid_spikes = current_spikes(current_spikes <= max_time_2);
    
    if ~isempty(valid_spikes)
       plot(valid_spikes, q * ones(size(valid_spikes)), '|', ...
            'MarkerSize', 20, 'LineWidth', 1.0, 'Color', colors(q,:));
    end  
end


xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Hy experiments', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Inset: Raster Plot of analyzed traces', n_files), 'FontSize', 14, 'FontWeight', 'bold');
xlim([0 max_time_2 * 1.02]);  
ylim([0.5 n_files + 0.5]);
yticks(1:n_files);

if isstruct(file_names)
    yticklabels({file_names.name});
else
    yticklabels(file_names);
end

set(gca, 'FontSize', 10, 'LineWidth', 1.5, 'GridAlpha', 0.3);

line([max_time_2 max_time_2], ylim, 'Color', 'r', 'LineStyle', ':', 'LineWidth', 2);
text(max_time_2, n_files+0.2, sprintf('%.3f s', max_time_2), ...
    'Color', 'r', 'FontSize', 10, 'HorizontalAlignment', 'right');
figure('Name','Control traces vs analized traces');
hold on;
num_tests = size(integral_time_peak, 2);

for q = 1:num_tests
    
    full_peaks = integral_time_peak(:, q);
    valid_peaks = valid_integral_time_peak(:, q);
       
    full_peaks = full_peaks(~isnan(full_peaks));
    valid_peaks = valid_peaks(~isnan(valid_peaks));
    
    plot(full_peaks, q * ones(size(full_peaks)), '|', ...
         'MarkerSize', 20, 'LineWidth', 1, 'Color', [0 0 0]);
  
    if q <= 5
        plot(valid_peaks, q * ones(size(valid_peaks)), '|', ...
             'MarkerSize', 10, 'LineWidth', 1.5, 'Color', [0.00, 0.45, 0.74]); % blu
    else
        plot(valid_peaks, q * ones(size(valid_peaks)), '|', ...
             'MarkerSize', 10, 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.10]); % rosso
    end
end

ylim([0.6 n_files + 0.4]);
xlim([0 max(max(integral_time_peak,[],'omitnan')) + 1]);
yticks(1:n_files);
line([max_time_2 max_time_2], ylim, 'Color', 'r', 'LineStyle', ':', 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Hy experiments', 'FontSize', 12, 'FontWeight', 'bold');
title('Raster Plot: Control Traces vs Analyzed Traces');
if isstruct(file_names)
    yticklabels({file_names.name});
else
    yticklabels(file_names);
end

