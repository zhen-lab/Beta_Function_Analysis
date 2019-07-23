headPos = 35;

curvdatafilteredCutD = (curvdatafiltered>0) .* curvdatafiltered;
curvdatafilteredCutV = (curvdatafiltered<0) .* curvdatafiltered;

[corr2(dorsal_smd(headPos:end, :), curvdatafilteredCutD(headPos+1:end, :));  corr2(ventral_smd(headPos:end, :), -curvdatafilteredCutV(headPos+1:end, :))]