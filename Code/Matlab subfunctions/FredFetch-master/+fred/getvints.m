function [vints] = getvints(series)

  vints = fred.dispatch_(0, 0, @fred.getvints_, series);

end
