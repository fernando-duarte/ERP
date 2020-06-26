function [ output_args, info_args ] = haver_pull(haver_series,x_min,x_max,freq)
%% This function pulls data from haver and formats it into a matrix
% the matlab dates will be the leftmost column, and the rest are the series
% if you specify present, it will pull to present data, otherwise the
% format for the dates is 'mm/dd/yyyy' and the freq matches the information
% from the haver.fetch command frequency
% example
% my_series =
% haver_pull({'diff%(GDP@USECON)','LR@USECON'},'01/01/1981','present','Q');

%% Versioning %%
% Latest Version: Dec 2013
% By: Michael Kubiske

if strcmp(x_max,'present')
    x_max=date;
end
fprintf('Connecting to Haver...\n');
for i=1:length(haver_series)
    idx1=strfind(haver_series{i},'(');
    idx2=strfind(haver_series{i},'@');
    idx3=strfind(haver_series{i},')');
    if isempty(idx3)
        idx3=length(haver_series{i});
    else
        idx3=length(haver_series{i})-1;
    end
    if isempty(idx1)
        idx1=1;
    else
        idx1=idx1+1;
    end
    series{i}=haver_series{i}(idx1:idx2-1);
    database{i}=haver_series{i}(idx2+1:idx3);
    if idx1==1
        transform{i}='no';
    else
        transform{i}=haver_series{i}(1:idx1-2);
    end
end
for i=1:length(series)
    error=1;
    while ~isempty(error)
        try hav_con=haver(cat(2, 'K:\DLX\DATA\',database{i},'.dat'));
            error=[];
        catch error
            fprintf('Database connection problem, trying again.');
            fprintf('    %s\n',error.message);
        end
    end
    fprintf('Pulling Series %s with %s transformation.\n',series{i},transform{i});
    error=1;
    while ~isempty(error)
        try eval(cat(2, strrep(transform{i},'%','_'),series{i},'=fetch(hav_con,char(series{i}),x_min,x_max,freq);'));
            error=[];
        catch error
            fprintf('Database connection problem, trying again.\n');
            fprintf('    %s\n',error.message);
        end
    end
    error=1;
    while ~isempty(error)
        try info_args{i} = info(hav_con,char(series{i}));
            error=[];
        catch error
            fprintf('Database connection problem, trying again.\n');
            fprintf('    %s\n',error.message);
        end
    end
    close(hav_con);
end
fprintf('Connection Closed.\n');
for i=1:length(transform)
    temp='none';
    if ~strcmp(transform{i},'no')
        eval(cat(2, strrep(transform{i},'%','_'),series{i},'=haver_transform(',strrep(transform{i},'%','_'),series{i},',transform{i});'));
    else
        eval(cat(2, strrep(transform{i},'%','_'),series{i},'=haver_transform(',strrep(transform{i},'%','_'),series{i},',temp);'));
    end
    eval(cat(2, 'startdates(i)=',strrep(transform{i},'%','_'),series{i},'(1,1);'));
    eval(cat(2, 'enddates(i)=',strrep(transform{i},'%','_'),series{i},'(end,1);'));
    eval(cat(2, strrep(transform{i},'%','_'),series{i},'(:,1)=datenum(year(',strrep(transform{i},'%','_'),series{i},'(:,1)),month(',strrep(transform{i},'%','_'),series{i},'(:,1)),day(',strrep(transform{i},'%','_'),series{i},'(:,1)));'));
end

if length(series)~=1
    for i=1:length(series)-1
        if i==1
            eval(cat(2, 'final_dates=union(',strrep(transform{i},'%','_'),series{i},'(:,1),',strrep(transform{i},'%','_'),series{i+1},'(:,1));'));
        else
            eval(cat(2, 'final_dates=union(final_dates,',strrep(transform{i+1},'%','_'),series{i+1},'(:,1));'));
        end
    end
else
    eval(cat(2, 'final_dates=',strrep(transform{i},'%','_'),series{i},'(:,1);'));
end
output_args(:,1)=final_dates;
warn_toggle=0;
if length(haver_series)~=1
    for j=1:length(haver_series)
        for i=1:length(final_dates)
            temp=[];
            eval(cat(2, 'temp=find(final_dates(i)==',strrep(transform{j},'%','_'),series{j},'(:,1));'));
            if isempty(temp)
                if warn_toggle==0
                    fprintf('A series does not have data for all dates, the matrix will have NaNs.\n');
                    warn_toggle=1;
                end
                output_args(i,j+1)=NaN(1,1);
            else
                eval(cat(2, 'output_args(i,j+1)=',strrep(transform{j},'%','_'),series{j},'(temp,2);'));
            end
        end
    end
else
    eval(cat(2, 'output_args=',strrep(transform{1},'%','_'),series{1},';'));
end
% final end of program end
end

