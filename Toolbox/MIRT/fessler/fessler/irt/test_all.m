% test_all.m
%
% Attempt to run all test routines
% to verify completeness of software distribution.
%
% I recommend running these tests to verify the completeness
% of your installation!
%
% Jeff Fessler

%double6 double
%disp 'Note: forcing double precision because sparse(double)*single fails :-('

prompt draw % do all plots without pausing/prompting
%printm 'Note: hit "r" at the prompt to disable subsequent prompts'

% uncomment the suite(s) you want to test

test_all_mex	% see if all the mex files work (will fail on Windows)
test_all_util
test_all_reg
test_all_nufft
test_all_systems
test_all_emission
test_all_transmission
test_all_wls
test_all_ct
test_all_example
