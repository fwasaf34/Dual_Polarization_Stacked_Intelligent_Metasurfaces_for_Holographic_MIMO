clc
clear
close all



for seed =1:100
    rng(seed)
    c=3*10^8;
    f=28*10^9;
    lambda=c/f;
    zTx=0;
    zRx=250;
    noisedBm=-110;
    noise=db2pow(noisedBm)*10^-3;
    PdBm=20;
    P=db2pow(PdBm)*10^-3;
    alp=0.2;
    rx=32;
    tx=256;
    GTilde=zeros(rx, tx);
    PL=20*log10(4*pi/lambda)+10*3.5*log10(zRx)+9;
    rho_sq=1/db2pow(PL);
    GTilde = sqrt((1-alp)*rho_sq/2) * (randn(rx,tx) + 1j*randn(rx,tx));
    G=GTilde;
    SE(seed) = calculate_spectral_efficiency(G, P, noise);
    EE(seed) = calculate_spectral_efficiency(G, P, noise)/P;
end

mean(SE)
mean(EE)



