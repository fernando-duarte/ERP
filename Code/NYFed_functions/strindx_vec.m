function indx = strindx(str_vec,search_str,exact);

N = size(search_str,1);
indx = [];

for n = 1:N,
    
    test = strfind(str_vec,deblank(search_str(n,:)));
    indx_n = [];

    if isempty(test),
        indx_n = [];
        return;
    else,
        for i = 1:length(test),
            if ~isempty(test{i}),
                if nargin == 3 && ~isempty(strmatch(str_vec(i,:),deblank(search_str(n,:)),exact)),
                    indx_n = [indx_n i];
                elseif nargin == 2,
                    indx_n = [indx_n i];
                end;            
            end;
        end;
    end;    
    
    indx{n} = indx_n;
    
end; % for n = 1:N





