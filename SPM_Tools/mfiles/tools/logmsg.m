function logmsg(file,msg)
% Save comment in log file

fid = fopen(file,'at');
if fid == -1
    error('Log file cannot be edited.');
end
t = clock;
wall_clock = sprintf('%4d/%02d/%02d - %02d:%02d:%02d',...
		t(1),t(2),t(3),t(4),t(5),floor(t(6)));
string = [wall_clock '  : ' strrep(msg,'\','\\') '\n'] ;
fprintf(fid, string);
fclose(fid);
fprintf(string);
