function n = days_distance(filestr,formatstr)
% function n = days_distance(filestr,formatstr)
% distance in days between date in filestr and nu
% filestr: date string, formatted using formatstr
% formatstr: something like 'dd-mmm-yyyy HH:MM:SS'
% output: n, a float
% Nikos, February 2014
start_n = datenum(filestr,formatstr);
finish_n = datenum(datestr(now,formatstr),formatstr);
n = finish_n - start_n;