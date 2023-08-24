function [params_pH] = pH_interpolation(pH_for_calc, pHs, pKas, pKas_err, pA, pA_err, pT, pT_err, pGS, pGS_err,...
    kGStES1, kGStES1_err, kES1tGS, kES1tGS_err, kGStES2, kGStES2_err, kES2tGS, kES2tGS_err,...
    kES2tES1, kES2tES1_err, kES1tES2, kES1tES2_err, slopes_GStES2, slopes_minor)

%make calculations with fractions (not percents)

pT_pH      = pT./100;
pT_pH_err  = pT_err./100;

pA_pH      = (10.^(pH_for_calc - pKas))./(1+10.^(pH_for_calc - pKas));
pA_pH_err  = pA_pH.*sqrt((pKas_err.*log(10)).^2+((10.^(pH_for_calc - pKas).*(pKas_err.*log(10)))./(1+10.^(pH_for_calc - pKas))).^2);

pA_pH(pHs==pH_for_calc,:)     = pA(pHs==pH_for_calc,:)./100;
pA_pH_err(pHs==pH_for_calc,:) = pA_err(pHs==pH_for_calc,:)./100;

pGS_pH     = 1 - pA_pH - pT_pH;
pGS_pH_err = sqrt((pA_pH_err).^2+(pT_pH_err).^2);

pGS_pH(pHs==pH_for_calc,:)     = pGS(pHs==pH_for_calc,:)./100;
pGS_pH_err(pHs==pH_for_calc,:) = pGS_err(pHs==pH_for_calc,:)./100;

kGStES1_pH               = kGStES1;
kGStES1_pH_err           = kGStES1_err;
kES1tGS_pH               = kES1tGS;
kES1tGS_pH_err           = kES1tGS_err;


kGStES2_for_pH_ln        = log(kGStES2);
kGStES2_for_pH_ln_err    = kGStES2_err./kGStES2;
inter_GStES2_for_pH      = kGStES2_for_pH_ln - slopes_GStES2(1,3).*pHs;
inter_GStES2_for_pH_err  = sqrt((kGStES2_for_pH_ln_err).^2+(pHs).^2.*(slopes_GStES2(2,3)).^2);
kGStES2_pH_ln            = slopes_GStES2(1,3).*pH_for_calc + inter_GStES2_for_pH;
kGStES2_pH_ln_err        = sqrt((inter_GStES2_for_pH_err).^2+(pH_for_calc).^2.*(slopes_GStES2(2,3)).^2);
kGStES2_pH               = exp(kGStES2_pH_ln);
kGStES2_pH_err           = kGStES2_pH.*kGStES2_pH_ln_err;

kGStES2_pH(pHs==pH_for_calc,:)     = kGStES2(pHs==pH_for_calc,:);
kGStES2_pH_err(pHs==pH_for_calc,:) = kGStES2_err(pHs==pH_for_calc,:);

kES2tGS_pH               = kES2tGS;
kES2tGS_pH_err           = kES2tGS_err;


kES2tES1_for_pH_ln       = log(kES2tES1);
kES2tES1_for_pH_ln_err   = kES2tES1_err./kES2tES1;
inter_ES2tES1_for_pH     = kES2tES1_for_pH_ln - slopes_minor(1,1).*pHs;
inter_ES2tES1_for_pH_err = sqrt((kES2tES1_for_pH_ln_err).^2+(pHs).^2.*(slopes_minor(2,1)).^2);
kES2tES1_pH_ln           = slopes_minor(1,1).*pH_for_calc + inter_ES2tES1_for_pH;
kES2tES1_pH_ln_err       = sqrt((inter_ES2tES1_for_pH_err).^2+(pH_for_calc).^2.*(slopes_minor(2,1)).^2);
kES2tES1_pH              = exp(kES2tES1_pH_ln);
kES2tES1_pH_err          = kES2tES1_pH.*kES2tES1_pH_ln_err;

kES2tES1_pH(pHs==pH_for_calc,:)     = kES2tES1(pHs==pH_for_calc,:);
kES2tES1_pH_err(pHs==pH_for_calc,:) = kES2tES1_err(pHs==pH_for_calc,:);

kES1tES2_for_pH_ln       = log(kES1tES2);
kES1tES2_for_pH_ln_err   = kES1tES2_err./kES1tES2;
inter_ES1tES2_for_pH     = kES1tES2_for_pH_ln - slopes_minor(1,2).*pHs;
inter_ES1tES2_for_pH_err = sqrt((kES1tES2_for_pH_ln_err).^2+(pHs).^2.*(slopes_minor(2,2)).^2);
kES1tES2_pH_ln           = slopes_minor(1,2).*pH_for_calc + inter_ES1tES2_for_pH;
kES1tES2_pH_ln_err       = sqrt((inter_ES1tES2_for_pH_err).^2+(pH_for_calc).^2.*(slopes_minor(2,2)).^2);
kES1tES2_pH              = exp(kES1tES2_pH_ln);
kES1tES2_pH_err          = kES1tES2_pH.*kES1tES2_pH_ln_err;

kES1tES2_pH(pHs==pH_for_calc,:)     = kES1tES2(pHs==pH_for_calc,:);
kES1tES2_pH_err(pHs==pH_for_calc,:) = kES1tES2_err(pHs==pH_for_calc,:);


kES1_pH     = kGStES1_pH + kES1tGS_pH;
kES1_pH_err = sqrt(kGStES1_pH_err.^2 + kES1tGS_pH_err.^2);

kES2_pH     = kGStES2_pH + kES2tGS_pH;
kES2_pH_err = sqrt(kGStES2_pH_err.^2 + kES2tGS_pH_err.^2);

kMinor_pH     = kES1tES2_pH + kES2tES1_pH;
kMinor_pH_err = sqrt(kES1tES2_pH_err.^2 + kES2tES1_pH_err.^2);

%return parameters (populations are fractions, not percents)
params_pH = [pGS_pH,pT_pH,pA_pH,...
    kES1_pH,kES2_pH,kMinor_pH,...
    kGStES1_pH,kES1tGS_pH,kGStES2_pH,kES2tGS_pH,kES1tES2_pH,kES2tES1_pH,...
    pGS_pH_err,pT_pH_err,pA_pH_err,...
    kES1_pH_err,kES2_pH_err,kMinor_pH_err,...
    kGStES1_pH_err,kES1tGS_pH_err,kGStES2_pH_err,kES2tGS_pH_err,kES1tES2_pH_err,kES2tES1_pH_err];
