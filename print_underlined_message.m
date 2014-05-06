function print_underlined_message(str1,str2)
% print a message
screenwidth = 100;
s = sprintf(' %s _____________ %s\n',str2,str1);
underscoreslength = screenwidth-length(s);
underscoresstr = '';
for ii = 1:underscoreslength
    underscoresstr = strcat(underscoresstr,'_');
end
fprintf(underscoresstr);
fprintf(s);