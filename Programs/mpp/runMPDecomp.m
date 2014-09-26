% Run MP decomposition on a ready local.ctl file
% Vinay Shirhatti, 03 Sep 2014
%==========================================================================

function runMPDecomp(filename)

sourcefolder = '~/Documents/MP/source/'; % Vinay
commandline = [sourcefolder 'gabord ' filename];
unix(commandline)

end