function spec_lowres = get_lowresspec(wvl,spec,FWHM)
% spec_lowres = get_lowresspec(wvl,spec,FHWM)

DL = abs(wvl(2)-wvl(1));

instr_wvl = -3*FWHM:DL:3*FWHM;
instr_function = FWHM/pi./(instr_wvl.^2+FWHM^2);

spec_lowres = conv(spec,instr_function,'same');