function Rsq = RSquaredValue(ydata, yest) % finds the R squared value of your data
  % �Copyright 2023 University of Florida Research Foundation, Inc. All Commercial Rights Reserved.
    SSres = sum((ydata - yest).^2);
    SStot = sum((ydata - mean(ydata)).^2);
    Rsq = 1 - SSres / SStot;
end
