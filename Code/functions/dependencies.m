%%% Getting and Identifying dependencies functions

file_name = 'strindx_vec.m';

[fList, pList] = matlab.codetools.requiredFilesAndProducts(file_name);
g = getcallinfo(file_name);
f = g(1).calls.fcnCalls.names;
cellfun(@(x) which(x), unique(f))