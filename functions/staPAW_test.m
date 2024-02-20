function [test_lls] = staPAW_test(Data_test, mmhat, num_subsamp, len_data)

%% staPAW_test
% input the Data structure for testing, the fitted parameters, and the
% number and length of subsamples
% the output is a vector of testLLs from subsamples, normalized in length

%% extract test Data
[xxf, yyf, alltrials, time] = data2xy(Data_test);  %Data
alldis = extractfield(Data_test, 'dis');  % Data

% creating new xx,yy for fitting if needed
yyf = [yyf; alldis];
% xxf = [xxf; allopto];
maskf = alltrials;  %((alltrials)==0) = false;%

%% sampling tests
n_samp = num_subsamp;  % repeats
l_samp = len_data;  % length
test_lls = zeros(1, n_samp);
for nn = 1:n_samp
    randomInteger = randi([1, length(yyf)-l_samp]);%randi([1,min(wind)]);%
    wind_i = randomInteger:randomInteger+l_samp;
    yy = yyf(:, wind_i);
    xx = xxf(:, wind_i);
    mask = maskf(wind_i);
    [logp_test,gams_test,xis_test, xisum_test] = runFB_GLMHMM_xi(mmhat,xx,yy,mask);
    test_lls(nn) = logp_test/length(wind_i);
end

end