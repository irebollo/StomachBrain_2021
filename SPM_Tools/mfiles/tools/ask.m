function ok = ask(confirm, string)

switch confirm
    case 'no', ok = 1;
    case 'yes',
        resp = input([string ' (y/n)? '],'s');
        if strcmp(resp,'y'), ok = 1;
        else ok = 0;
        end
    otherwise
        error('ask: confirm should be either "yes" or "no".')
end
