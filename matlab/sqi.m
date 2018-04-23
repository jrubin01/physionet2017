function result = sqi(directory, fname)

	data_dir = [pwd filesep directory filesep];

	%% Add this directory to the MATLAB path.
	addpath(pwd)

    result = randn(1)

