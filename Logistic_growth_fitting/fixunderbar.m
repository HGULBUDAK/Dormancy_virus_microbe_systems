function outstring = fixunderbar(instring)
% outstring = fixunderbar(instring)
% replace "_" by "\_"
% in a string
%
% for use with matgraph
%
% uses sed
% and is cool
outstring = strrep(instring,'_','\_');